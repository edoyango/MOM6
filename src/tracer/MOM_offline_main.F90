! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> The routines here implement the offline tracer algorithm used in MOM6. These are called from step_offline
!! Some routines called here can be found in the MOM_offline_aux module.
module MOM_offline_main

use MOM_ALE,                  only : ALE_CS, ALE_regrid, ALE_offline_inputs
use MOM_ALE,                  only : pre_ALE_adjustments, ALE_update_regrid_weights
use MOM_ALE,                  only : ALE_remap_tracers
use MOM_checksums,            only : hchksum, uvchksum
use MOM_coms,                 only : reproducing_sum
use MOM_cpu_clock,            only : cpu_clock_id, cpu_clock_begin, cpu_clock_end
use MOM_cpu_clock,            only : CLOCK_COMPONENT, CLOCK_SUBCOMPONENT
use MOM_cpu_clock,            only : CLOCK_MODULE_DRIVER, CLOCK_MODULE, CLOCK_ROUTINE
use MOM_diabatic_aux,         only : diabatic_aux_CS, set_pen_shortwave
use MOM_diabatic_driver,      only : diabatic_CS, extract_diabatic_member
use MOM_diabatic_aux,         only : tridiagTS
use MOM_diag_mediator,        only : diag_ctrl, post_data, register_diag_field
use MOM_domains,              only : pass_var, pass_vector
use MOM_error_handler,        only : MOM_error, MOM_mesg, FATAL, WARNING
use MOM_error_handler,        only : callTree_enter, callTree_leave
use MOM_file_parser,          only : read_param, get_param, log_version, param_file_type
use MOM_forcing_type,         only : forcing
use MOM_grid,                 only : ocean_grid_type
use MOM_interface_heights,    only : calc_derived_thermo, thickness_to_dz
use MOM_io,                   only : MOM_read_data, MOM_read_vector, CENTER
use MOM_offline_aux,          only : update_offline_from_arrays, update_offline_from_files
use MOM_offline_aux,          only : next_modulo_time, offline_add_diurnal_sw
use MOM_offline_aux,          only : update_h_horizontal_flux, update_h_vertical_flux, limit_mass_flux_3d
use MOM_offline_aux,          only : distribute_residual_uh_barotropic, distribute_residual_vh_barotropic
use MOM_offline_aux,          only : distribute_residual_uh_upwards, distribute_residual_vh_upwards
use MOM_opacity,              only : opacity_CS, optics_type
use MOM_open_boundary,        only : ocean_OBC_type
use MOM_time_manager,         only : time_type, real_to_time
use MOM_tracer_advect,        only : tracer_advect_CS, advect_tracer
use MOM_tracer_diabatic,      only : applyTracerBoundaryFluxesInOut
use MOM_tracer_flow_control,  only : tracer_flow_control_CS, call_tracer_column_fns, call_tracer_stocks
use MOM_tracer_registry,      only : tracer_registry_type, MOM_tracer_chksum, MOM_tracer_chkinv
use MOM_unit_scaling,         only : unit_scale_type
use MOM_variables,            only : thermo_var_ptrs
use MOM_verticalGrid,         only : verticalGrid_type, get_thickness_units

implicit none ; private

#include "MOM_memory.h"

!> The control structure for the offline transport module
type, public :: offline_transport_CS ; private

  ! Pointers to relevant fields from the main MOM control structure
  type(ALE_CS),                  pointer :: ALE_CSp         => NULL()
          !< A pointer to the ALE control structure
  type(diabatic_CS),             pointer :: diabatic_CSp    => NULL()
          !< A pointer to the diabatic control structure
  type(diag_ctrl),               pointer :: diag            => NULL()
          !< Structure that regulates diagnostic output
  type(ocean_OBC_type),          pointer :: OBC             => NULL()
          !< A pointer to the open boundary condition control structure
  type(tracer_advect_CS),        pointer :: tracer_adv_CSp  => NULL()
          !< A pointer to the tracer advection control structure
  type(opacity_CS),              pointer :: opacity_CSp     => NULL()
          !< A pointer to the opacity control structure
  type(tracer_flow_control_CS),  pointer :: tracer_flow_CSp => NULL()
          !< A pointer to control structure that orchestrates the calling of tracer packages
  type(tracer_registry_type),    pointer :: tracer_Reg      => NULL()
          !< A pointer to the tracer registry
  type(thermo_var_ptrs),         pointer :: tv              => NULL()
          !< A structure pointing to various thermodynamic variables
  type(optics_type),             pointer :: optics          => NULL()
          !< Pointer to the optical properties type
  type(diabatic_aux_CS),         pointer :: diabatic_aux_CSp => NULL()
          !< Pointer to the diabatic_aux control structure

  !> Variables related to reading in fields from online run
  integer :: start_index  !< Timelevel to start
  integer :: iter_no      !< Timelevel to start
  integer :: numtime      !< How many timelevels in the input fields
  type(time_type) :: accumulated_time !< Length of time accumulated in the current offline interval
  type(time_type) :: vertical_time !< The next value of accumulate_time at which to apply vertical processes
  ! Index of each of the variables to be read in with separate indices for each variable if they
  ! are set off from each other in time
  integer :: ridx_sum = -1 !< Read index offset of the summed variables
  integer :: ridx_snap = -1 !< Read index offset of the snapshot variables
  integer :: nk_input     !< Number of input levels in the input fields
  character(len=200) :: offlinedir  !< Directory where offline fields are stored
  character(len=200) :: & ! Names of input files
    surf_file,  &         !< Contains surface fields (2d arrays)
    snap_file,  &         !< Snapshotted fields (layer thicknesses)
    sum_file,   &         !< Fields which are accumulated over time
    mean_file             !< Fields averaged over time
  character(len=20)  :: redistribute_method !< 'barotropic' if evenly distributing extra flow
                                            !! throughout entire watercolumn, 'upwards',
                                            !! if trying to do it just in the layers above
                                            !! 'both' if both methods are used
  character(len=20) :: mld_var_name !< Name of the mixed layer depth variable to use
  logical :: fields_are_offset !< True if the time-averaged fields and snapshot fields are
                               !! offset by one time level
  logical :: x_before_y        !< Which horizontal direction is advected first
  logical :: print_adv_offline !< Prints out some updates each advection sub interation
  logical :: skip_diffusion    !< Skips horizontal diffusion of tracers
  logical :: read_sw           !< Read in averaged values for shortwave radiation
  logical :: read_mld          !< Check to see whether mixed layer depths should be read in
  real    :: Hmix_fixed        !< A fixed mixed layer depth to use when read_mld is false [Z ~> m]
  logical :: diurnal_sw        !< Adds a synthetic diurnal cycle on shortwave radiation
  logical :: debug             !< If true, write verbose debugging messages
  logical :: redistribute_barotropic !< Redistributes column-summed residual transports throughout
                                     !! a column weighted by thickness
  logical :: redistribute_upwards    !< Redistributes remaining fluxes only in layers above
                                     !! the current one based as the max allowable transport
                                     !! in that cell
  logical :: read_all_ts_uvh  !< If true, then all timelevels of temperature, salinity, mass transports, and
                              !! Layer thicknesses are read during initialization
  !! Variables controlling some of the numerical considerations of offline transport
  integer :: num_off_iter   !< Number of advection iterations per offline step
  integer :: num_vert_iter  !< Number of vertical iterations per offline step
  integer :: off_ale_mod    !< Sets how frequently the ALE step is done during the advection
  real :: dt_offline        !< Timestep used for offline tracers [T ~> s]
  real :: dt_offline_vertical !< Timestep used for calls to tracer vertical physics [T ~> s]
  real :: evap_CFL_limit    !< Limit on the fraction of the water that can be fluxed out of the top
                            !! layer in a timestep [nondim].  This is Copied from diabatic_CS controlling
                            !! how tracers follow freshwater fluxes
  real :: minimum_forcing_depth !< The smallest depth over which fluxes can be applied [H ~> m or kg m-2].
                            !! This is copied from diabatic_CS controlling how tracers follow freshwater fluxes

  real :: Kd_max        !< Runtime parameter specifying the maximum value of vertical diffusivity
                        !! [H Z T-1 ~> m2 s-1 or kg m-1 s-1]
  real :: min_residual  !< The minimum amount of total mass flux before exiting the main advection
                        !! routine [H L2 ~> m3 or kg]
  !>@{ Diagnostic manager IDs for some fields that may be of interest when doing offline transport
  integer :: &
    id_uhr = -1, id_vhr = -1, &
    ! Unused: id_ear = -1, id_ebr = -1, &
    id_hr = -1,  &
    id_hdiff = -1, &
    id_uhr_redist = -1, &
    id_vhr_redist = -1, &
    id_uhr_end = -1, &
    id_vhr_end = -1, &
    id_eta_pre_distribute  = -1, &
    id_eta_post_distribute = -1, &
    id_h_redist = -1, &
    id_eta_diff_end = -1

  ! Diagnostic IDs for the regridded/remapped input fields
  integer :: &
    id_uhtr_regrid = -1, &
    id_vhtr_regrid = -1, &
    id_temp_regrid = -1, &
    id_salt_regrid = -1, &
    id_h_regrid = -1
  !>@}

  ! IDs for timings of various offline components
  integer :: id_clock_read_fields = -1   !< A CPU time clock
  integer :: id_clock_offline_diabatic = -1  !< A CPU time clock
  integer :: id_clock_offline_adv  = -1  !< A CPU time clock
  integer :: id_clock_redistribute = -1  !< A CPU time clock

  !> Zonal transport that may need to be stored between calls to step_MOM [H L2 ~> m3 or kg]
  real, allocatable, dimension(:,:,:) :: uhtr
  !> Meridional transport that may need to be stored between calls to step_MOM [H L2 ~> m3 or kg]
  real, allocatable, dimension(:,:,:) :: vhtr

  ! Fields at T-point
  real, allocatable, dimension(:,:,:) :: eatr
                   !< Amount of fluid entrained from the layer above within
                   !! one time step [H ~> m or kg m-2]
  real, allocatable, dimension(:,:,:) :: ebtr
                   !< Amount of fluid entrained from the layer below within
                   !! one time step [H ~> m or kg m-2]
  ! Fields at T-points on interfaces
  real, allocatable, dimension(:,:,:) :: Kd     !< Vertical diffusivity [H Z T-1 ~> m2 s-1 or kg m-1 s-1]
  real, allocatable, dimension(:,:,:) :: h_end  !< Thicknesses at the end of offline timestep [H ~> m or kg m-2]

  real, allocatable, dimension(:,:) :: mld        !< Mixed layer depths at thickness points [Z ~> m]

  ! Allocatable arrays to read in entire fields during initialization
  real, allocatable, dimension(:,:,:,:) :: uhtr_all !< Entire field of zonal transport [H L2 ~> m3 or kg]
  real, allocatable, dimension(:,:,:,:) :: vhtr_all !< Entire field of meridional transport [H L2 ~> m3 or kg]
  real, allocatable, dimension(:,:,:,:) :: hend_all !< Entire field of layer thicknesses [H ~> m or kg m-2]
  real, allocatable, dimension(:,:,:,:) :: temp_all !< Entire field of temperatures [C ~> degC]
  real, allocatable, dimension(:,:,:,:) :: salt_all !< Entire field of salinities [S ~> ppt]

end type offline_transport_CS

public offline_advection_ale
public offline_redistribute_residual
public offline_diabatic_ale
public offline_fw_fluxes_into_ocean
public offline_fw_fluxes_out_ocean
public offline_advection_layer
public register_diags_offline_transport
public update_offline_fields
public insert_offline_main
public extract_offline_main
public post_offline_convergence_diags
public offline_transport_init
public offline_transport_end


  interface
module subroutine offline_advection_ale(fluxes, Time_start, time_interval, G, GV, US, CS, id_clock_ale, &
                                 h_pre, uhtr, vhtr, converged)
  type(forcing),           intent(inout) :: fluxes        !< pointers to forcing fields
  type(time_type),         intent(in)    :: Time_start    !< starting time of a segment, as a time type
  real,                    intent(in)    :: time_interval !< time interval covered by this call [T ~> s]
  type(ocean_grid_type),   intent(inout) :: G             !< Ocean grid structure
  type(verticalGrid_type), intent(in)    :: GV            !< Vertical grid structure
  type(unit_scale_type),   intent(in)    :: US            !< A dimensional unit scaling type
  type(offline_transport_CS), pointer    :: CS            !< control structure for offline module
  integer,                 intent(in)    :: id_clock_ALE  !< Clock for ALE routines
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(inout) :: h_pre         !< layer thicknesses before advection
                                                          !! [H ~> m or kg m-2]
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                           intent(inout) :: uhtr          !< Zonal mass transport [H L2 ~> m3 or kg]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                           intent(inout) :: vhtr          !< Meridional mass transport [H L2 ~> m3 or kg]
  logical,                 intent(  out) :: converged     !< True if the iterations have converged

  ! Local variables


  ! Variables used to keep track of layer thicknesses at various points in the code
                                               ! in the same units as thicknesses [H ~> m or kg m-2]
                          ! top layer in a timestep [nondim]
end subroutine offline_advection_ale
module subroutine offline_redistribute_residual(CS, G, GV, US, h_pre, uhtr, vhtr, converged)
  type(offline_transport_CS), pointer       :: CS    !< control structure from initialize_MOM
  type(ocean_grid_type),      intent(inout) :: G     !< Ocean grid structure
  type(verticalGrid_type),    intent(in)    :: GV    !< Vertical grid structure
  type(unit_scale_type),      intent(in)    :: US    !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                              intent(inout) :: h_pre !< layer thicknesses before advection [H ~> m or kg m-2]
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                              intent(inout) :: uhtr  !< Zonal mass transport [H L2 ~> m3 or kg]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                              intent(inout) :: vhtr  !< Meridional mass transport [H L2 ~> m3 or kg]
  logical,                    intent(in   ) :: converged !< True if the iterations have converged

  ! Variables used to keep track of layer thicknesses at various points in the code

  ! Used to calculate the eta diagnostics


end subroutine offline_redistribute_residual
real module function remaining_transport_sum(G, GV, US, uhtr, vhtr, h_new)
  type(ocean_grid_type),      intent(in)    :: G     !< Ocean grid structure
  type(verticalGrid_type),    intent(in)    :: GV    !< Vertical grid structure
  type(unit_scale_type),      intent(in)    :: US    !< A dimensional unit scaling type
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                              intent(in   ) :: uhtr  !< Zonal mass transport [H L2 ~> m3 or kg]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                              intent(in   ) :: vhtr  !< Meridional mass transport [H L2 ~> m3 or kg]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                              intent(in   ) :: h_new !< Layer thicknesses [H ~> m or kg m-2]

  ! Local variables
                     !! transports through the faces of a column [R Z L2 ~> kg].
                     !! of a tracer cell [H L2 ~> m3 or kg]

end function remaining_transport_sum
module subroutine offline_diabatic_ale(fluxes, Time_start, Time_end, G, GV, US, CS, h_pre, tv, eatr, ebtr)

  type(forcing),           intent(inout) :: fluxes     !< pointers to forcing fields
  type(time_type),         intent(in)    :: Time_start !< starting time of a segment, as a time type
  type(time_type),         intent(in)    :: Time_end   !< ending time of a segment, as a time type
  type(ocean_grid_type),   intent(in)    :: G          !< Ocean grid structure
  type(verticalGrid_type), intent(in)    :: GV         !< Vertical grid structure
  type(unit_scale_type),   intent(in)    :: US         !< A dimensional unit scaling type
  type(offline_transport_CS), pointer    :: CS         !< control structure from initialize_MOM
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(inout) :: h_pre      !< layer thicknesses before advection [H ~> m or kg m-2]
  type(thermo_var_ptrs),   intent(in   ) :: tv         !< A structure pointing to various thermodynamic variables
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(inout) :: eatr       !< Entrainment from layer above [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(inout) :: ebtr       !< Entrainment from layer below [H ~> m or kg m-2]

  ! Local variables
end subroutine offline_diabatic_ale
module subroutine offline_fw_fluxes_into_ocean(G, GV, CS, fluxes, h, in_flux_optional)
  type(offline_transport_CS), intent(inout) :: CS !< Offline control structure
  type(ocean_grid_type),      intent(in)    :: G  !< Grid structure
  type(verticalGrid_type),    intent(in)    :: GV !< ocean vertical grid structure
  type(forcing),              intent(inout) :: fluxes !< Surface fluxes container
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                              intent(inout) :: h  !< Layer thickness [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJ_(G)), &
                    optional, intent(in)    :: in_flux_optional !< The total time-integrated amount
                                                  !! of tracer that leaves with freshwater
                                                  !! [CU H ~> Conc m or Conc kg m-2]


end subroutine offline_fw_fluxes_into_ocean
module subroutine offline_fw_fluxes_out_ocean(G, GV, CS, fluxes, h, out_flux_optional)
  type(offline_transport_CS), intent(inout) :: CS !< Offline control structure
  type(ocean_grid_type),      intent(in)    :: G  !< Grid structure
  type(verticalGrid_type),    intent(in)    :: GV !< ocean vertical grid structure
  type(forcing),              intent(inout) :: fluxes !< Surface fluxes container
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                              intent(inout) :: h  !< Layer thickness [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJ_(G)), &
                    optional, intent(in)    :: out_flux_optional !< The total time-integrated amount
                                                  !! of tracer that leaves with freshwater
                                                  !! [CU H ~> Conc m or Conc kg m-2]


end subroutine offline_fw_fluxes_out_ocean
module subroutine offline_advection_layer(fluxes, Time_start, time_interval, G, GV, US, CS, h_pre, eatr, ebtr, uhtr, vhtr)
  type(forcing),              intent(inout) :: fluxes        !< pointers to forcing fields
  type(time_type),            intent(in)    :: Time_start    !< starting time of a segment, as a time type
  real,                       intent(in)    :: time_interval !< Offline transport time interval [T ~> s]
  type(ocean_grid_type),      intent(inout) :: G             !< Ocean grid structure
  type(verticalGrid_type),    intent(in)    :: GV            !< Vertical grid structure
  type(unit_scale_type),      intent(in)    :: US            !< A dimensional unit scaling type
  type(offline_transport_CS), pointer       :: CS            !< Control structure for offline module
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                              intent(inout) :: h_pre !< layer thicknesses before advection [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                              intent(inout) :: eatr !< Entrainment from layer above [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                              intent(inout) :: ebtr !< Entrainment from layer below [H ~> m or kg m-2]
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                              intent(inout) :: uhtr  !< Zonal mass transport [H L2 ~> m3 or kg]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                              intent(inout) :: vhtr  !< Meridional mass transport [H L2 ~> m3 or kg]

  ! Local variables


                         ! mass fluxes through the faces of a column or within a column [R Z L2 ~> kg]
                         ! used to keep track of how close to convergence we are.

  ! Variables used to keep track of layer thicknesses at various points in the code
  ! Work arrays for temperature and salinity

end subroutine offline_advection_layer
module subroutine update_offline_fields(CS, G, GV, US, h, fluxes, do_ale)
  type(offline_transport_CS), pointer       :: CS !< Control structure for offline module
  type(ocean_grid_type),      intent(inout) :: G  !< Ocean grid structure
  type(verticalGrid_type),    intent(in)    :: GV !< Vertical grid structure
  type(unit_scale_type),      intent(in)    :: US !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                              intent(inout) :: h !< The regridded layer thicknesses [H ~> m or kg m-2]
  type(forcing),              intent(inout) :: fluxes !< Pointers to forcing fields
  logical,                    intent(in   ) :: do_ale !< True if using ALE
  ! Local variables
end subroutine update_offline_fields
module subroutine register_diags_offline_transport(Time, diag, CS, GV, US)

  type(offline_transport_CS), pointer :: CS   !< Control structure for offline module
  type(verticalGrid_type), intent(in) :: GV   !< Vertical grid structure
  type(unit_scale_type),   intent(in) :: US   !< A dimensional unit scaling type
  type(time_type),         intent(in) :: Time !< current model time
  type(diag_ctrl),         intent(in) :: diag !< Structure that regulates diagnostic output

  ! U-cell fields
end subroutine register_diags_offline_transport
module subroutine post_offline_convergence_diags(G, GV, CS, h_off, h_end, uhtr, vhtr)
  type(ocean_grid_type),      intent(in)    :: G      !< Ocean grid structure
  type(verticalGrid_type),    intent(in)    :: GV     !< Vertical grid structure
  type(offline_transport_CS), intent(in   ) :: CS     !< Offline control structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                              intent(inout) :: h_off  !< Thicknesses at end of offline step [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                              intent(inout) :: h_end  !< Stored thicknesses [H ~> m or kg m-2]
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                              intent(inout) :: uhtr   !< Remaining zonal mass transport [H L2 ~> m3 or kg]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                              intent(inout) :: vhtr   !< Remaining meridional mass transport [H L2 ~> m3 or kg]


end subroutine post_offline_convergence_diags
module subroutine extract_offline_main(CS, uhtr, vhtr, eatr, ebtr, h_end, accumulated_time, vertical_time, &
                                dt_offline, dt_offline_vertical, skip_diffusion)
  type(offline_transport_CS), target, intent(in   ) :: CS !< Offline control structure
  ! Returned optional arguments
  real, dimension(:,:,:), optional, pointer       :: uhtr !< Remaining zonal mass transport [H L2 ~> m3 or kg]
  real, dimension(:,:,:), optional, pointer       :: vhtr !< Remaining meridional mass transport [H L2 ~> m3 or kg]
  real, dimension(:,:,:), optional, pointer       :: eatr !< Amount of fluid entrained from the layer above within
                                                          !! one time step [H ~> m or kg m-2]
  real, dimension(:,:,:), optional, pointer       :: ebtr !< Amount of fluid entrained from the layer below within
                                                          !! one time step [H ~> m or kg m-2]
  real, dimension(:,:,:), optional, pointer       :: h_end !< Thicknesses at the end of offline timestep
                                                          !! [H ~> m or kg m-2]
  type(time_type),        optional, pointer       :: accumulated_time !< Length of time accumulated in the
                                                          !! current offline interval
  type(time_type),        optional, pointer       :: vertical_time !< The next value of accumulate_time at which to
                                                          !! vertical processes
  real,                   optional, intent(  out) :: dt_offline !< Timestep used for offline tracers [T ~> s]
  real,                   optional, intent(  out) :: dt_offline_vertical !< Timestep used for calls to tracer
                                                          !! vertical physics [T ~> s]
  logical,                optional, intent(  out) :: skip_diffusion !< Skips horizontal diffusion of tracers

  ! Pointers to 3d members
end subroutine extract_offline_main
module subroutine insert_offline_main(CS, ALE_CSp, diabatic_CSp, diag, OBC, tracer_adv_CSp, &
                               tracer_flow_CSp, tracer_Reg, tv, x_before_y, debug)
  type(offline_transport_CS), intent(inout) :: CS  !< Offline control structure
  ! Inserted optional arguments
  type(ALE_CS), &
            target, optional, intent(in   ) :: ALE_CSp  !< A pointer to the ALE control structure
  type(diabatic_CS), &
            target, optional, intent(in   ) :: diabatic_CSp !< A pointer to the diabatic control structure
  type(diag_ctrl), &
            target, optional, intent(in   ) :: diag     !< A pointer to the structure that regulates diagnostic output
  type(ocean_OBC_type), &
            target, optional, intent(in   ) :: OBC      !< A pointer to the open boundary condition control structure
  type(tracer_advect_CS), &
            target, optional, intent(in   ) :: tracer_adv_CSp !< A pointer to the tracer advection control structure
  type(tracer_flow_control_CS), &
            target, optional, intent(in   ) :: tracer_flow_CSp !< A pointer to the tracer flow control control structure
  type(tracer_registry_type), &
            target, optional, intent(in   ) :: tracer_Reg !< A pointer to the tracer registry
  type(thermo_var_ptrs), &
            target, optional, intent(in   ) :: tv       !< A structure pointing to various thermodynamic variables
  logical,          optional, intent(in   ) :: x_before_y !< Indicates which horizontal direction is advected first
  logical,          optional, intent(in   ) :: debug    !< If true, write verbose debugging messages


end subroutine insert_offline_main
module subroutine offline_transport_init(param_file, CS, diabatic_CSp, G, GV, US)

  type(param_file_type),           intent(in) :: param_file !< A structure to parse for run-time parameters
  type(offline_transport_CS),      pointer    :: CS !< Offline control structure
  type(diabatic_CS),               intent(in) :: diabatic_CSp !< The diabatic control structure
  type(ocean_grid_type),   target, intent(in) :: G  !< ocean grid structure
  type(verticalGrid_type), target, intent(in) :: GV !< ocean vertical grid structure
  type(unit_scale_type),   target, intent(in) :: US !< A dimensional unit scaling type

  ! This include declares and sets the variable "version".

end subroutine offline_transport_init
module subroutine read_all_input(CS, G, GV, US)
  type(offline_transport_CS), intent(inout) :: CS    !< Control structure for offline module
  type(ocean_grid_type),      intent(in)    :: G     !< Ocean grid structure
  type(verticalGrid_type),    intent(in)    :: GV    !< Vertical grid structure
  type(unit_scale_type),      intent(in)    :: US    !< A dimensional unit scaling type


end subroutine read_all_input
module subroutine offline_transport_end(CS)
  type(offline_transport_CS), pointer :: CS !< Control structure for offline module

  ! Explicitly allocate all allocatable arrays
end subroutine offline_transport_end
  end interface

end module MOM_offline_main
