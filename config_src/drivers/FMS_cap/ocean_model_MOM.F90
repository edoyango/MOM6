! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> Top-level module for the MOM6 ocean model in coupled mode.
module ocean_model_mod

! This is the top level module for the MOM6 ocean model.  It contains routines
! for initialization, termination and update of ocean model state.  This
! particular version wraps all of the calls for MOM6 in the calls that had
! been used for MOM4.
!
! This code is a stop-gap wrapper of the MOM6 code to enable it to be called
! in the same way as MOM4.

use MOM, only : initialize_MOM, step_MOM, MOM_control_struct, MOM_end
use MOM, only : extract_surface_state, allocate_surface_state, finish_MOM_initialization
use MOM, only : get_MOM_state_elements, MOM_state_is_synchronized
use MOM, only : get_ocean_stocks, step_offline
use MOM, only : save_MOM_restart
use MOM_coms,      only : field_chksum
use MOM_constants, only : CELSIUS_KELVIN_OFFSET, hlf
use MOM_coupler_types, only : coupler_1d_bc_type, coupler_2d_bc_type
use MOM_coupler_types, only : coupler_type_spawn, coupler_type_write_chksums
use MOM_coupler_types, only : coupler_type_initialized, coupler_type_copy_data
use MOM_coupler_types, only : coupler_type_set_diags, coupler_type_send_data
use MOM_diag_mediator, only : diag_ctrl, enable_averages, disable_averaging
use MOM_diag_mediator, only : diag_mediator_close_registration, diag_mediator_end
use MOM_domains, only : MOM_domain_type, domain2d, clone_MOM_domain, get_domain_extent
use MOM_domains, only : pass_var, pass_vector, AGRID, BGRID_NE, CGRID_NE, TO_ALL, Omit_Corners
use MOM_error_handler, only : MOM_error, MOM_mesg, FATAL, WARNING, is_root_pe
use MOM_error_handler, only : callTree_enter, callTree_leave
use MOM_EOS, only : gsw_sp_from_sr, gsw_pt_from_ct
use MOM_file_parser, only : get_param, log_version, close_param_file, param_file_type
use MOM_forcing_type, only : forcing, mech_forcing, allocate_forcing_type
use MOM_forcing_type, only : fluxes_accumulate, get_net_mass_forcing
use MOM_forcing_type, only : forcing_diagnostics, mech_forcing_diags
use MOM_get_input, only : Get_MOM_Input, directories
use MOM_grid, only : ocean_grid_type
use MOM_io, only : write_version_number, stdout_if_root
use MOM_marine_ice, only : iceberg_forces, iceberg_fluxes, marine_ice_init, marine_ice_CS
use MOM_string_functions, only : uppercase
use MOM_surface_forcing_gfdl, only : surface_forcing_init, convert_IOB_to_fluxes
use MOM_surface_forcing_gfdl, only : convert_IOB_to_forces, ice_ocn_bnd_type_chksum
use MOM_surface_forcing_gfdl, only : ice_ocean_boundary_type, surface_forcing_CS
use MOM_surface_forcing_gfdl, only : forcing_save_restart
use MOM_time_manager, only : time_type, operator(>), operator(+), operator(-)
use MOM_time_manager, only : operator(*), operator(/), operator(/=)
use MOM_time_manager, only : operator(<=), operator(>=), operator(<)
use MOM_time_manager, only : real_to_time, time_to_real
use MOM_tracer_flow_control, only : call_tracer_register, tracer_flow_control_init
use MOM_tracer_flow_control, only : call_tracer_flux_init
use MOM_unit_scaling, only : unit_scale_type
use MOM_variables, only : surface
use MOM_verticalGrid, only : verticalGrid_type
use MOM_ice_shelf, only : initialize_ice_shelf, shelf_calc_flux, ice_shelf_CS
use MOM_ice_shelf, only : initialize_ice_shelf_fluxes, initialize_ice_shelf_forces
use MOM_ice_shelf, only : add_shelf_forces, ice_shelf_end, ice_shelf_save_restart
use MOM_ice_shelf, only : ice_sheet_calving_to_ocean_sfc, adjust_ice_sheet_frazil
use MOM_wave_interface, only: wave_parameters_CS, MOM_wave_interface_init
use MOM_wave_interface, only: Update_Surface_Waves
use iso_fortran_env, only : int64

#include <MOM_memory.h>

#ifdef _USE_GENERIC_TRACER
use MOM_generic_tracer, only : MOM_generic_tracer_fluxes_accumulate
#endif

implicit none ; private

public ocean_model_init, ocean_model_end, update_ocean_model
public ocean_model_save_restart, Ocean_stock_pe
public ice_ocean_boundary_type
public ocean_model_init_sfc, ocean_model_flux_init
public ocean_model_restart
public ice_ocn_bnd_type_chksum
public ocean_public_type_chksum
public ocean_model_data_get
public get_ocean_grid
public ocean_model_get_UV_surf

!> This interface extracts a named scalar field or array from the ocean surface or public type
interface ocean_model_data_get
  module procedure ocean_model_data1D_get
  module procedure ocean_model_data2D_get
end interface


!> This type is used for communication with other components via the FMS coupler.
!! The element names and types can be changed only with great deliberation, hence
!! the persistence of things like the cutesy element name "avg_kount".
type, public ::  ocean_public_type
  type(domain2d) :: Domain    !< The domain for the surface fields.
  logical :: is_ocean_pe      !< .true. on processors that run the ocean model.
  character(len=32) :: instance_name = '' !< A name that can be used to identify
                                 !! this instance of an ocean model, for example
                                 !! in ensembles when writing messages.
  integer, pointer, dimension(:) :: pelist => NULL()   !< The list of ocean PEs.
  logical, pointer, dimension(:,:) :: maskmap =>NULL() !< A pointer to an array
                    !! indicating which logical processors are actually used for
                    !! the ocean code. The other logical processors would be all
                    !! land points and are not assigned to actual processors.
                    !! This need not be assigned if all logical processors are used.

  integer :: stagger = -999   !< The staggering relative to the tracer points
                    !! points of the two velocity components. Valid entries
                    !! include AGRID, BGRID_NE, CGRID_NE, BGRID_SW, and CGRID_SW,
                    !! corresponding to the community-standard Arakawa notation.
                    !! (These are named integers taken from the MOM_domains module.)
                    !! Following MOM5, stagger is BGRID_NE by default when the
                    !! ocean is initialized, but here it is set to -999 so that
                    !! a global max across ocean and non-ocean processors can be
                    !! used to determine its value.
  real, pointer, dimension(:,:)  :: &
    t_surf => NULL(), & !< SST on t-cell [degrees Kelvin]
    s_surf => NULL(), & !< SSS on t-cell [ppt]
    u_surf => NULL(), & !< i-velocity at the locations indicated by stagger [m s-1].
    v_surf => NULL(), & !< j-velocity at the locations indicated by stagger [m s-1].
    sea_lev => NULL(), & !< Sea level in m after correction for surface pressure,
                        !! i.e. dzt(1) + eta_t + patm/rho0/grav [m]
    frazil =>NULL(), &  !< Accumulated heating [J m-2] from frazil
                        !! formation in the ocean.
    melt_potential => NULL(), & !< Instantaneous heat used to melt sea ice [J m-2].
    OBLD => NULL(),   & !< Ocean boundary layer depth [m].
    area => NULL(),   & !< cell area of the ocean surface [m2].
    calving => NULL(), &!< The mass per unit area of the ice shelf to convert to
                        !! bergs [kg m-2].
    calving_hflx => NULL() !< Calving heat flux [W m-2].
  type(coupler_2d_bc_type) :: fields    !< A structure that may contain named
                                        !! arrays of tracer-related surface fields.
  integer                  :: avg_kount !< A count of contributions to running
                                        !! sums, used externally by the FMS coupler
                                        !! for accumulating averages of this type.
  integer, dimension(2)    :: axes = 0  !< Axis numbers that are available
                                        !! for I/O using this surface data.
end type ocean_public_type


!> The ocean_state_type contains all information about the state of the ocean,
!! with a format that is private so it can be readily changed without disrupting
!! other coupled components.
type, public :: ocean_state_type ; private
  ! This type is private, and can therefore vary between different ocean models.
  logical :: is_ocean_PE = .false.  !< True if this is an ocean PE.
  type(time_type) :: Time     !< The ocean model's time and master clock.
  type(time_type) :: Time_dyn !< The ocean model's time for the dynamics.  Time and Time_dyn
                              !! should be the same after a full time step.
  integer :: Restart_control  !< An integer that is bit-tested to determine whether
                              !! incremental restart files are saved and whether they
                              !! have a time stamped name.  +1 (bit 0) for generic
                              !! files and +2 (bit 1) for time-stamped files.  A
                              !! restart file is saved at the end of a run segment
                              !! unless Restart_control is negative.

  integer :: nstep = 0        !< The number of calls to update_ocean that update the dynamics.
  integer :: nstep_thermo = 0 !< The number of calls to update_ocean that update the thermodynamics.
  logical :: use_ice_shelf    !< If true, the ice shelf model is enabled.
  logical :: use_waves        !< If true use wave coupling.

  logical :: icebergs_alter_ocean !< If true, the icebergs can change ocean the
                              !! ocean dynamics and forcing fluxes.
  real :: press_to_z          !< A conversion factor between pressure and ocean depth,
                              !! usually 1/(rho_0*g) [Z T2 R-1 L-2 ~> m Pa-1].
  logical :: calve_ice_shelf_bergs = .false. !< If true, bergs are initialized according to
                              !! ice shelf flux through the ice front
  real :: C_p                 !< The heat capacity of seawater [J degC-1 kg-1].
  logical :: offline_tracer_mode = .false. !< If false, use the model in prognostic mode
                              !! with the barotropic and baroclinic dynamics, thermodynamics,
                              !! etc. stepped forward integrated in time.
                              !! If true, all of the above are bypassed with all
                              !! fields necessary to integrate only the tracer advection
                              !! and diffusion equation read in from files stored from
                              !! a previous integration of the prognostic model.

  logical :: single_step_call !< If true, advance the state of MOM with a single
                              !! step including both dynamics and thermodynamics.
                              !! If false, the two phases are advanced with
                              !! separate calls. The default is true.
  ! The following 3 variables are only used here if single_step_call is false.
  real    :: dt               !< (baroclinic) dynamics time step [T ~> s]
  real    :: dt_therm         !< thermodynamics time step [T ~> s]
  logical :: thermo_spans_coupling !< If true, thermodynamic and tracer time
                              !! steps can span multiple coupled time steps.
  logical :: diabatic_first   !< If true, apply diabatic and thermodynamic
                              !! processes before time stepping the dynamics.

  type(directories) :: dirs   !< A structure containing several relevant directory paths.
  type(mech_forcing)          :: forces  !< A structure with the driving mechanical surface forces
  type(forcing)               :: fluxes  !< A structure containing pointers to
                                                    !! the thermodynamic ocean forcing fields.
  type(forcing)               :: flux_tmp !< A secondary structure containing pointers to the
                              !! ocean forcing fields for when multiple coupled
                              !! timesteps are taken per thermodynamic step.
  type(surface)               :: sfc_state   !< A structure containing pointers to
                              !! the ocean surface state fields.
  type(ocean_grid_type), pointer :: &
    grid => NULL()            !< A pointer to a grid structure containing metrics
                              !! and related information.
  type(verticalGrid_type), pointer :: &
    GV => NULL()              !< A pointer to a structure containing information
                              !! about the vertical grid.
  type(unit_scale_type), pointer :: &
    US => NULL()              !< A pointer to a structure containing dimensional
                              !! unit scaling factors.
  type(MOM_control_struct) :: MOM_CSp
                              !< MOM control structure
  type(ice_shelf_CS), pointer :: &
    Ice_shelf_CSp => NULL()   !< A pointer to the control structure for the
                              !! ice shelf model that couples with MOM6.  This
                              !! is null if there is no ice shelf.
  type(marine_ice_CS), pointer :: &
    marine_ice_CSp => NULL()  !< A pointer to the control structure for the
                              !! marine ice effects module.
  type(wave_parameters_cs), pointer :: &
    Waves => NULL()           !< A pointer to the surface wave control structure
  type(surface_forcing_CS), pointer :: &
    forcing_CSp => NULL()     !< A pointer to the MOM forcing control structure
  type(diag_ctrl), pointer :: &
    diag => NULL()            !< A pointer to the diagnostic regulatory structure
end type ocean_state_type


  interface
module subroutine ocean_model_init(Ocean_sfc, OS, Time_init, Time_in, wind_stagger, gas_fields_ocn, calve_ice_shelf_bergs)
  type(ocean_public_type), target, &
                       intent(inout) :: Ocean_sfc !< A structure containing various publicly
                                !! visible ocean surface properties after initialization,
                                !! the data in this type is intent out.
  type(ocean_state_type), pointer    :: OS        !< A structure whose internal
                                !! contents are private to ocean_model_mod that may be used to
                                !! contain all information about the ocean's interior state.
  type(time_type),     intent(in)    :: Time_init !< The start time for the coupled model's calendar
  type(time_type),     intent(in)    :: Time_in   !< The time at which to initialize the ocean model.
  integer, optional,   intent(in)    :: wind_stagger !< If present, the staggering of the winds that are
                                                     !! being provided in calls to update_ocean_model
  type(coupler_1d_bc_type), &
             optional, intent(in)    :: gas_fields_ocn !< If present, this type describes the
                                              !! ocean and surface-ice fields that will participate
                                              !! in the calculation of additional gas or other
                                              !! tracer fluxes, and can be used to spawn related
                                              !! internal variables in the ice model.
  logical, optional,   intent(in)    :: calve_ice_shelf_bergs !< If true, track ice shelf flux through a
                                              !! static ice shelf, so that it can be converted into icebergs
  ! Local variables
                      !! The actual depth over which melt potential is computed will
                      !! min(HFrz, OBLD), where OBLD is the boundary layer depth.
                      !! If HFrz <= 0 (default), melt potential will not be computed.

  ! This include declares and sets the variable "version".
                                ! surface velocities returned to the coupler.

end subroutine ocean_model_init
module subroutine update_ocean_model(Ice_ocean_boundary, OS, Ocean_sfc, time_start_update, &
                              Ocean_coupling_time_step, update_dyn, update_thermo, &
                              Ocn_fluxes_used, start_cycle, end_cycle, cycle_length)
  type(ice_ocean_boundary_type), &
                     intent(in)    :: Ice_ocean_boundary !< A structure containing the various
                                              !! forcing fields coming from the ice and atmosphere.
  type(ocean_state_type), &
                     pointer       :: OS      !< A pointer to a private structure containing the
                                              !! internal ocean state.
  type(ocean_public_type), &
                     intent(inout) :: Ocean_sfc !< A structure containing all the publicly visible
                                              !! ocean surface fields after a coupling time step.
                                              !! The data in this type is intent out.
  type(time_type),   intent(in)    :: time_start_update  !< The time at the beginning of the update step.
  type(time_type),   intent(in)    :: Ocean_coupling_time_step !< The amount of time over which to
                                              !! advance the ocean.
  logical, optional, intent(in)    :: update_dyn !< If present and false, do not do updates
                                              !! due to the ocean dynamics.
  logical, optional, intent(in)    :: update_thermo !< If present and false, do not do updates
                                              !! due to the ocean thermodynamics or remapping.
  logical, optional, intent(in)    :: Ocn_fluxes_used !< If present, this indicates whether the
                                              !! cumulative thermodynamic fluxes from the ocean,
                                              !! like frazil, have been used and should be reset.
  logical, optional, intent(in)    :: start_cycle !< This indicates whether this call is to be
                                              !! treated as the first call to step_MOM in a
                                              !! time-stepping cycle; missing is like true.
  logical, optional, intent(in)    :: end_cycle   !< This indicates whether this call is to be
                                              !! treated as the last call to step_MOM in a
                                              !! time-stepping cycle; missing is like true.
  real,    optional, intent(in)    :: cycle_length !< The duration of a coupled time stepping cycle [s].

  ! Local variables
                            ! start of this call to allow step_MOM to temporarily change the time
                            ! as seen by internal modules.
                            ! this call to allow step_MOM to temporarily change the time as seen by
                            ! internal modules.

end subroutine update_ocean_model
module subroutine ocean_model_restart(OS, timestamp)
  type(ocean_state_type),     pointer    :: OS !< A pointer to the structure containing the
                                               !! internal ocean state being saved to a restart file
  character(len=*), optional, intent(in) :: timestamp !< An optional timestamp string that should be
                                               !! prepended to the file name. (Currently this is unused.)

end subroutine ocean_model_restart
module subroutine ocean_model_end(Ocean_sfc, Ocean_state, Time)
  type(ocean_public_type), intent(inout) :: Ocean_sfc   !< An ocean_public_type structure that is
                                                        !! to be deallocated upon termination.
  type(ocean_state_type),  pointer       :: Ocean_state !< A pointer to the structure containing
                                                        !! the internal ocean state to be deallocated
                                                        !! upon termination.
  type(time_type),         intent(in)    :: Time        !< The model time, used for writing restarts.

end subroutine ocean_model_end
module subroutine ocean_model_save_restart(OS, Time, directory, filename_suffix)
  type(ocean_state_type),     pointer    :: OS  !< A pointer to the structure containing the
                                                !! internal ocean state (in).
  type(time_type),            intent(in) :: Time !< The model time at this call, needed for writing files.
  character(len=*), optional, intent(in) :: directory  !<  An optional directory into which to
                                                !! write these restart files.
  character(len=*), optional, intent(in) :: filename_suffix !< An optional suffix (e.g., a time-stamp)
                                                !! to append to the restart file names.
! Note: This is a new routine - it will need to exist for the new incremental
!   checkpointing.  It will also be called by ocean_model_end, giving the same
!   restart behavior as now in FMS.

end subroutine ocean_model_save_restart
module subroutine initialize_ocean_public_type(input_domain, Ocean_sfc, diag, gas_fields_ocn)
  type(MOM_domain_type),   intent(in)    :: input_domain !< The ocean model domain description
  type(ocean_public_type), intent(inout) :: Ocean_sfc !< A structure containing various publicly
                                              !! visible ocean surface properties after
                                              !! initialization, whose elements are allocated here.
  type(diag_ctrl),         intent(in)    :: diag  !< A structure that regulates diagnostic output
  type(coupler_1d_bc_type), &
                 optional, intent(in)    :: gas_fields_ocn !< If present, this type describes the
                                              !! ocean and surface-ice fields that will participate
                                              !! in the calculation of additional gas or other
                                              !! tracer fluxes.

  ! ice-ocean-boundary fields are always allocated using absolute indices
  ! and have no halos.

end subroutine initialize_ocean_public_type
module subroutine convert_state_to_ocean_type(sfc_state, Ocean_sfc, G, US, patm, press_to_z)
  type(surface),         intent(inout) :: sfc_state !< A structure containing fields that
                                               !! describe the surface state of the ocean.
  type(ocean_public_type), &
                 target, intent(inout) :: Ocean_sfc !< A structure containing various publicly
                                               !! visible ocean surface fields, whose elements
                                               !! have their data set here.
  type(ocean_grid_type), intent(inout) :: G    !< The ocean's grid structure
  type(unit_scale_type), intent(in)    :: US   !< A dimensional unit scaling type
  real,        optional, intent(in)    :: patm(:,:)  !< The pressure at the ocean surface [R L2 T-2 ~> Pa]
  real,        optional, intent(in)    :: press_to_z !< A conversion factor between pressure and ocean
                                               !! depth, usually 1/(rho_0*g) [Z T2 R-1 L-2 ~> m Pa-1]
  ! Local variables

end subroutine convert_state_to_ocean_type
module subroutine convert_shelf_state_to_ocean_type(Ocean_sfc, CS, US)
  type(ocean_public_type), &
               target, intent(inout) :: Ocean_sfc !< A structure containing various publicly
                                                  !! visible ocean surface fields, whose elements
                                                  !! have their data set here.
  type(ice_shelf_CS),      pointer :: CS        !< A pointer to the ice shelf control structure
  type(unit_scale_type), intent(in)    :: US   !< A dimensional unit scaling type

end subroutine convert_shelf_state_to_ocean_type
module subroutine ocean_model_init_sfc(OS, Ocean_sfc)
  type(ocean_state_type),  pointer       :: OS  !< The structure with the complete ocean state
  type(ocean_public_type), intent(inout) :: Ocean_sfc !< A structure containing various publicly
                                !! visible ocean surface properties after initialization, whose
                                !! elements have their data set here.

end subroutine ocean_model_init_sfc
module subroutine ocean_model_flux_init(OS, verbosity)
  type(ocean_state_type), optional, pointer :: OS  !< An optional pointer to the ocean state,
                                             !! used to figure out if this is an ocean PE that
                                             !! has already been initialized.
  integer, optional, intent(in) :: verbosity !< A 0-9 integer indicating a level of verbosity.


end subroutine ocean_model_flux_init
module subroutine Ocean_stock_pe(OS, index, value, time_index)
  use stock_constants_mod, only : ISTOCK_WATER, ISTOCK_HEAT,ISTOCK_SALT
  type(ocean_state_type), pointer     :: OS         !< A structure containing the internal ocean state.
                                                    !! The data in OS is intent in.
  integer,                intent(in)  :: index      !< The stock index for the quantity of interest.
  real,                   intent(out) :: value      !< Sum returned for the conservation quantity of interest [various]
  integer,      optional, intent(in)  :: time_index !< An unused optional argument, present only for
                                                    !! interfacial compatibility with other models.
! Arguments: OS - A structure containing the internal ocean state.
!  (in)      index - Index of conservation quantity of interest.
!  (in)      value -  Sum returned for the conservation quantity of interest.
!  (in,opt)  time_index - Index for time level to use if this is necessary.


end subroutine Ocean_stock_pe
module subroutine ocean_model_data2D_get(OS, Ocean, name, array2D, isc, jsc)
  use MOM_constants, only : CELSIUS_KELVIN_OFFSET
  type(ocean_state_type),     pointer    :: OS    !< A pointer to the structure containing the
                                                  !! internal ocean state (intent in).
  type(ocean_public_type),    intent(in) :: Ocean !< A structure containing various publicly
                                                  !! visible ocean surface fields.
  character(len=*)          , intent(in) :: name  !< The name of the field to extract
  integer                   , intent(in) :: isc   !< The starting i-index of array2D
  integer                   , intent(in) :: jsc   !< The starting j-index of array2D
  real, dimension(isc:,jsc:), intent(out):: array2D !< The values of the named field, it must
                                                  !! cover only the computational domain [various]


end subroutine ocean_model_data2D_get
module subroutine ocean_model_data1D_get(OS, Ocean, name, value)
  type(ocean_state_type),     pointer    :: OS    !< A pointer to the structure containing the
                                                  !! internal ocean state (intent in).
  type(ocean_public_type),    intent(in) :: Ocean !< A structure containing various publicly
                                                  !! visible ocean surface fields.
  character(len=*),           intent(in) :: name  !< The name of the field to extract
  real,                       intent(out):: value !< The value of the named field [various]

end subroutine ocean_model_data1D_get
module subroutine ocean_public_type_chksum(id, timestep, ocn)

  character(len=*),        intent(in) :: id  !< An identifying string for this call
  integer,                 intent(in) :: timestep !< The number of elapsed timesteps
  type(ocean_public_type), intent(in) :: ocn !< A structure containing various publicly
                                             !! visible ocean surface fields.
  ! Local variables

end subroutine ocean_public_type_chksum
module subroutine get_ocean_grid(OS, Gridp)
  ! Obtain the ocean grid.
  type(ocean_state_type) :: OS              !< A structure containing the
                                            !! internal ocean state
  type(ocean_grid_type) , pointer :: Gridp  !< The ocean's grid structure

end subroutine get_ocean_grid
module subroutine ocean_model_get_UV_surf(OS, Ocean, name, array2D, isc, jsc)

  type(ocean_state_type),     pointer    :: OS    !< A pointer to the structure containing the
                                                  !! internal ocean state (intent in).
  type(ocean_public_type),    intent(in) :: Ocean !< A structure containing various publicly
                                                  !! visible ocean surface fields.
  character(len=*)          , intent(in) :: name  !< The name of the current (ua or va) to extract
  integer                   , intent(in) :: isc   !< The starting i-index of array2D
  integer                   , intent(in) :: jsc   !< The starting j-index of array2D
  real, dimension(isc:,jsc:), intent(out):: array2D !< The values of the named field, it must
                                                  !! cover only the computational domain [L T-1 ~> m s-1]

                                                  !! describe the surface state of the ocean.


end subroutine ocean_model_get_UV_surf
  end interface

end module ocean_model_mod
