! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> Implements the thermodynamic aspects of ocean / ice-shelf interactions,
!!  along with a crude placeholder for a later implementation of full
!!  ice shelf dynamics, all using the MOM framework and coding style.
module MOM_ice_shelf

use MOM_array_transform,      only : rotate_array
use MOM_constants, only : hlf
use MOM_cpu_clock, only : cpu_clock_id, cpu_clock_begin, cpu_clock_end
use MOM_cpu_clock, only : CLOCK_COMPONENT, CLOCK_ROUTINE
use MOM_coms,                 only : num_PEs, reproducing_sum
use MOM_data_override,       only : data_override
use MOM_diag_mediator, only    : MOM_diag_ctrl=>diag_ctrl
use MOM_IS_diag_mediator, only : post_data=>post_IS_data, post_scalar_data=>post_IS_data_0d
use MOM_IS_diag_mediator, only : register_diag_field=>register_MOM_IS_diag_field, safe_alloc_ptr
use MOM_IS_diag_mediator, only : register_scalar_field=>register_MOM_IS_scalar_field
use MOM_IS_diag_mediator, only : set_IS_axes_info, diag_ctrl, time_type
use MOM_IS_diag_mediator, only : MOM_IS_diag_mediator_init, MOM_IS_diag_mediator_end
use MOM_IS_diag_mediator, only : set_IS_diag_mediator_grid
use MOM_IS_diag_mediator, only : enable_averages, disable_averaging
use MOM_IS_diag_mediator, only : MOM_IS_diag_mediator_infrastructure_init
use MOM_IS_diag_mediator, only : MOM_IS_diag_mediator_close_registration
use MOM_domains, only : MOM_domains_init, pass_var, pass_vector, clone_MOM_domain
use MOM_domains, only : TO_ALL, CGRID_NE, BGRID_NE, CORNER
use MOM_dyn_horgrid, only : dyn_horgrid_type, create_dyn_horgrid, destroy_dyn_horgrid
use MOM_error_handler, only : MOM_error, MOM_mesg, FATAL, WARNING, is_root_pe
use MOM_error_handler, only : callTree_showQuery
use MOM_error_handler, only : callTree_enter, callTree_leave, callTree_waypoint
use MOM_file_parser, only : read_param, get_param, log_param, log_version, param_file_type
use MOM_grid, only : MOM_grid_init, ocean_grid_type
use MOM_grid_initialize, only : initialize_masks, set_grid_metrics
use MOM_hor_index,             only : hor_index_type, hor_index_init
use MOM_hor_index,             only : rotate_hor_index
use MOM_fixed_initialization, only : MOM_initialize_topography
use MOM_fixed_initialization, only : MOM_initialize_rotation
use user_initialization, only : user_initialize_topography
use MOM_io, only : field_exists, file_exists, MOM_read_data, write_version_number
use MOM_io, only : slasher, fieldtype, vardesc, var_desc
use MOM_io, only : close_file, SINGLE_FILE, MULTIPLE
use MOM_restart, only : register_restart_field, save_restart
use MOM_restart, only : restart_init, restore_state, MOM_restart_CS, register_restart_pair
use MOM_time_manager, only : time_type, time_to_real, real_to_time, operator(>), operator(-)
use MOM_transcribe_grid, only : copy_dyngrid_to_MOM_grid, copy_MOM_grid_to_dyngrid
use MOM_transcribe_grid, only : rotate_dyngrid
use MOM_unit_scaling, only : unit_scale_type, unit_scaling_init, fix_restart_unit_scaling
use MOM_variables, only : surface, allocate_surface_state, deallocate_surface_state
use MOM_variables, only : rotate_surface_state
use MOM_forcing_type, only : forcing, allocate_forcing_type, deallocate_forcing_type, MOM_forcing_chksum
use MOM_forcing_type, only : mech_forcing, allocate_mech_forcing, deallocate_mech_forcing, MOM_mech_forcing_chksum
use MOM_forcing_type, only : copy_common_forcing_fields, rotate_forcing, rotate_mech_forcing
use MOM_get_input, only : directories, Get_MOM_input
use MOM_EOS, only : calculate_density, calculate_density_derivs, calculate_TFreeze, EOS_domain
use MOM_EOS, only : EOS_type, EOS_init
use MOM_ice_shelf_dynamics, only : ice_shelf_dyn_CS, update_ice_shelf, write_ice_shelf_energy
use MOM_ice_shelf_dynamics, only : register_ice_shelf_dyn_restarts, initialize_ice_shelf_dyn
use MOM_ice_shelf_dynamics, only : ice_shelf_min_thickness_calve, change_in_draft
use MOM_ice_shelf_dynamics, only : ice_time_step_CFL, ice_shelf_dyn_end, IS_dynamics_post_data
use MOM_ice_shelf_dynamics, only : volume_above_floatation, masked_var_grounded
use MOM_ice_shelf_initialize, only : initialize_ice_thickness
!MJH use MOM_ice_shelf_initialize, only : initialize_ice_shelf_boundary
use MOM_ice_shelf_state, only : ice_shelf_state, ice_shelf_state_end, ice_shelf_state_init
use user_shelf_init, only : USER_initialize_shelf_mass, USER_update_shelf_mass
use user_shelf_init, only : user_ice_shelf_CS
use MOM_spatial_means, only : global_area_integral
use MOM_checksums, only : hchksum, qchksum, chksum, uchksum, vchksum, uvchksum
use MOM_interpolate, only : init_external_field, time_interp_external, time_interp_external_init
use MOM_interpolate, only : external_field

implicit none ; private

#include <MOM_memory.h>
#ifdef SYMMETRIC_MEMORY_
#  define GRID_SYM_ .true.
#else
#  define GRID_SYM_ .false.
#endif

public shelf_calc_flux, initialize_ice_shelf, ice_shelf_end, ice_shelf_query
public ice_shelf_save_restart, solo_step_ice_shelf, add_shelf_forces
public initialize_ice_shelf_fluxes, initialize_ice_shelf_forces
public ice_sheet_calving_to_ocean_sfc
public adjust_ice_sheet_frazil

! A note on unit descriptions in comments: MOM6 uses units that can be rescaled for dimensional
! consistency testing. These are noted in comments with units like Z, H, L, and T, along with
! their mks counterparts with notation like "a velocity [Z T-1 ~> m s-1]".  If the units
! vary with the Boussinesq approximation, the Boussinesq variant is given first.

!> Control structure that contains ice shelf parameters and diagnostics handles
type, public :: ice_shelf_CS ; private
  ! Parameters
  type(MOM_restart_CS), pointer :: restart_CSp => NULL() !< A pointer to the restart control
                                                  !! structure for the ice shelves
  type(ocean_grid_type), pointer :: Grid_in => NULL() !< un-rotated input grid metric
  logical :: rotate_index = .false.   !< True if index map is rotated
  integer :: turns                    !< The number of quarter turns for rotation testing.
  type(ocean_grid_type), pointer :: Grid => NULL() !< Grid for the ice-shelf model
  type(unit_scale_type), pointer :: &
    US => NULL()       !< A structure containing various unit conversion factors
  type(ocean_grid_type), pointer :: ocn_grid => NULL() !< A pointer to the ocean model grid
                                          !! The rest is private
  real ::   flux_factor = 1.0             !< A factor that can be used to turn off ice shelf
                                          !! melting (flux_factor = 0) [nondim].
  character(len=128) :: restart_output_dir = ' ' !< The directory in which to write restart files
  type(ice_shelf_state), pointer :: ISS => NULL() !< A structure with elements that describe
                                          !! the ice-shelf state
  type(ice_shelf_dyn_CS), pointer :: dCS => NULL() !< The control structure for the ice-shelf dynamics.

  real, pointer, dimension(:,:) :: &
    utide   => NULL()  !< An unresolved tidal velocity [L T-1 ~> m s-1]

  real :: ustar_bg     !< A minimum value for ustar under ice shelves [Z T-1 ~> m s-1].
  real :: ustar_max    !< A maximum value for ustar under ice shelves, or a negative value to
                       !! have no limit [Z T-1 ~> m s-1].
  real :: cdrag        !< drag coefficient under ice shelves [nondim].
  real :: g_Earth      !< The gravitational acceleration [L2 Z-1 T-2 ~> m s-2]
  real :: Cp           !< The heat capacity of sea water [Q C-1 ~> J kg-1 degC-1].
  real :: Rho_ocn      !< A reference ocean density [R ~> kg m-3].
  real :: Cp_ice       !< The heat capacity of fresh ice [Q C-1 ~> J kg-1 degC-1].
  real :: gamma_t      !< The (fixed) turbulent exchange velocity in the
                       !< 2-equation formulation [Z T-1 ~> m s-1].
  real :: Salin_ice    !< The salinity of shelf ice [S ~> ppt].
  real :: Temp_ice     !< The core temperature of shelf ice [C ~> degC].
  real :: kv_ice       !< The viscosity of ice [L4 Z-2 T-1 ~> m2 s-1].
  real :: density_ice  !< A typical density of ice [R ~> kg m-3].
  real :: kv_molec     !< The molecular kinematic viscosity of sea water [Z2 T-1 ~> m2 s-1].
  real :: kd_molec_salt!< The molecular diffusivity of salt [Z2 T-1 ~> m2 s-1].
  real :: kd_molec_temp!< The molecular diffusivity of heat [Z2 T-1 ~> m2 s-1].
  real :: Lat_fusion   !< The latent heat of fusion [Q ~> J kg-1].
  real :: Gamma_T_3EQ  !< Nondimensional heat-transfer coefficient, used in the 3Eq. formulation [nondim]
  real :: Gamma_S_3EQ  !< Nondimensional salt-transfer coefficient, used in the 3Eq. formulation [nondim]
                       !< This number should be specified by the user.
  real :: col_mass_melt_threshold !< An ocean column mass below the iceshelf below which melting
                       !! does not occur [R Z ~> kg m-2]
  logical :: mass_from_file !< Read the ice shelf mass from a file every dt
  logical :: ustar_shelf_from_vel !< If true, use the surface velocities, and not the previous
                       !! values of the stresses to set ustar.

  !!!! PHYSICAL AND NUMERICAL PARAMETERS FOR ICE DYNAMICS !!!!!!

  real :: time_step    !< this is the shortest timestep that the ice shelf sees [T ~> s], and
                       !! is equal to the forcing timestep (it is passed in when the shelf
                       !! is initialized - so need to reorganize MOM driver.
                       !! it will be the prognostic timestep ... maybe.

  logical :: solo_ice_sheet !< whether the ice model is running without being
                            !! coupled to the ocean
  logical :: GL_regularize  !< whether to regularize the floatation condition
                            !! at the grounding line a la Goldberg Holland Schoof 2009
  logical :: GL_couple      !< whether to let the floatation condition be
                            !!determined by ocean column thickness means update_OD_ffrac
                            !! will be called (note: GL_regularize and GL_couple
                            !! should be exclusive)
  logical :: calve_to_mask  !< If true, calve any ice that passes outside of a masked area
  logical :: calve_ice_shelf_bergs=.false. !< If true, flux through a static ice front is converted
                                           !! to point bergs
  real :: min_thickness_simple_calve !< min. ice shelf thickness criteria for calving [Z ~> m].
  real :: T0                !< temperature at ocean surface in the restoring region [C ~> degC]
  real :: S0                !< Salinity at ocean surface in the restoring region [S ~> ppt].
  real :: input_flux        !< The vertically integrated inward ice thickness flux per
                            !! unit face length at an upstream boundary [Z L T-1 ~> m2 s-1]
  real :: input_thickness   !< Ice thickness at an upstream open boundary [Z ~> m].

  type(time_type) :: Time                !< The component's time.
  type(EOS_type) :: eqn_of_state         !< Type that indicates the equation of state to use.
  logical :: active_shelf_dynamics       !< True if the ice shelf mass changes as a result
                                         !! the dynamic ice-shelf model.
  logical :: shelf_mass_is_dynamic       !< True if ice shelf mass changes over time. If true, ice
                                         !! shelf dynamics will be initialized
  logical :: data_override_shelf_fluxes  !< True if the ice shelf surface mass fluxes can be
                                         !! written using the data_override feature (only for MOSAIC grids)
  logical :: override_shelf_movement     !< If true, user code specifies the shelf movement
                                         !! instead of using the dynamic ice-shelf mode.
  logical :: isthermo                    !< True if the ice shelf can exchange heat and
                                         !! mass with the underlying ocean.
  logical :: threeeq                     !< If true, the 3 equation consistency equations are
                                         !! used to calculate the flux at the ocean-ice
                                         !! interface.
  logical :: insulator                   !< If true, ice shelf is a perfect insulator
  logical :: const_gamma                 !< If true, gamma_T is specified by the user.
  logical :: constant_sea_level          !< if true, apply an evaporative, heat and salt
                                         !! fluxes. It will avoid large increase in sea level.
  logical :: constant_sea_level_misomip  !< If true, constant_sea_level fluxes are applied only over
                                         !! the surface sponge cells from the ISOMIP/MISOMIP configuration
  logical :: smb_diag                    !< If true, calculate diagnostics related to surface mass balance
  logical :: bmb_diag                    !< If true, calculate diagnostics related to basal mass balance
  real    :: min_ocean_mass_float        !< The minimum ocean mass per unit area before the ice
                                         !! shelf is considered to float when constant_sea_level
                                         !! is used [R Z ~> kg m-2]
  real    :: cutoff_depth                !< Depth above which melt is set to zero (>= 0) [Z ~> m].
  logical :: find_salt_root              !< If true, if true find Sbdry using a quadratic eq.
  real    :: TFr_0_0                     !< The freezing point at 0 pressure and 0 salinity [C ~> degC]
  real    :: dTFr_dS                     !< Partial derivative of freezing temperature with
                                         !! salinity [C S-1 ~> degC ppt-1]
  real    :: dTFr_dp                     !< Partial derivative of freezing temperature with
                                         !! pressure [C T2 R-1 L-2 ~> degC Pa-1]
  real    :: Zeta_N                      !< The stability constant xi_N = 0.052 from Holland & Jenkins '99
                                         !! divided by the von Karman constant VK [nondim]. Was 1/8.
  real :: Vk                             !< Von Karman's constant [nondim]
  real :: Rc                             !< critical flux Richardson number [nondim]
  logical :: ustar_from_vel_bugfix       !< If true, fixes ustar from ocean velocity bug
  logical :: buoy_flux_itt_bugfix        !< If true, fixes buoyancy iteration bug
  logical :: salt_flux_itt_bugfix        !< If true, fixes salt iteration bug
  real :: buoy_flux_tol                  !< Fractional buoyancy iteration tolerance for convergence [nondim]

  !>@{ Diagnostic handles
  integer :: id_melt = -1, id_exch_vel_s = -1, id_exch_vel_t = -1, &
             id_tfreeze = -1, id_tfl_shelf = -1, &
             id_thermal_driving = -1, id_haline_driving = -1, &
             id_u_ml = -1, id_v_ml = -1, id_sbdry = -1, &
             id_h_shelf = -1, id_dhdt_shelf = -1, id_h_mask = -1, id_frazil = -1, &
             id_surf_elev = -1, id_bathym = -1, &
             id_area_shelf_h = -1, &
             id_ustar_shelf = -1, id_shelf_mass = -1, id_mass_flux = -1, &
             id_shelf_sfc_mass_flux = -1, &
             id_vaf = -1, id_g_adott = -1, id_f_adott = -1, id_adott = -1, &
             id_bdott_melt = -1, id_bdott_accum = -1, id_bdott = -1, &
             id_dvafdt = -1, id_g_adot = -1, id_f_adot = -1, id_adot = -1, &
             id_bdot_melt = -1, id_bdot_accum = -1, id_bdot = -1, &
             id_t_area = -1, id_g_area = -1, id_f_area = -1, &
             id_Ant_vaf = -1, id_Ant_g_adott = -1, id_Ant_f_adott = -1, id_Ant_adott = -1, &
             id_Ant_bdott_melt = -1, id_Ant_bdott_accum = -1, id_Ant_bdott = -1, &
             id_Ant_dvafdt = -1, id_Ant_g_adot = -1, id_Ant_f_adot = -1, id_Ant_adot = -1, &
             id_Ant_bdot_melt = -1, id_Ant_bdot_accum = -1, id_Ant_bdot = -1, &
             id_Ant_t_area = -1, id_Ant_g_area = -1, id_Ant_f_area = -1, &
             id_Gr_vaf = -1, id_Gr_g_adott = -1, id_Gr_f_adott = -1, id_Gr_adott = -1, &
             id_Gr_bdott_melt = -1, id_Gr_bdott_accum = -1, id_Gr_bdott = -1, &
             id_Gr_dvafdt = -1, id_Gr_g_adot = -1, id_Gr_f_adot = -1, id_Gr_adot = -1, &
             id_Gr_bdot_melt = -1, id_Gr_bdot_accum = -1, id_Gr_bdot = -1, &
             id_Gr_t_area = -1, id_Gr_g_area = -1, id_Gr_f_area = -1
  !>@}

  type(external_field) :: mass_handle
    !< Handle for reading the time interpolated ice shelf mass from a file
  type(external_field) :: area_handle
    !< Handle for reading the time interpolated ice shelf area from a file

  type(diag_ctrl), pointer :: diag => NULL() !< A structure that is used to control diagnostic output.
  type(user_ice_shelf_CS), pointer :: user_CS => NULL() !< A pointer to the control structure for
                                  !! user-supplied modifications to the ice shelf code.

  logical :: debug                !< If true, write verbose checksums for debugging purposes
                                  !! and use reproducible sums
end type ice_shelf_CS

!>@{ CPU time clock IDs
integer :: id_clock_shelf=-1 !< CPU Clock for the ice shelf code
integer :: id_clock_pass=-1  !< CPU Clock for ice shelf group pass calls
!>@}


  interface
module subroutine shelf_calc_flux(sfc_state_in, fluxes_in, Time, time_step_in, CS)
  type(surface), target,  intent(inout) :: sfc_state_in !< A structure containing fields that
                                                 !! describe the surface state of the ocean.  The
                                                 !! intent is only inout to allow for halo updates.
  type(forcing),  target, intent(inout) :: fluxes_in !< structure containing pointers to any
                                                 !! possible thermodynamic or mass-flux forcing fields.
  type(time_type),        intent(in)    :: Time  !< Start time of the fluxes.
  real,                   intent(in)    :: time_step_in !< Length of time over which these fluxes
                                                 !! will be applied [T ~> s].
  type(ice_shelf_CS),     pointer       :: CS    !< A pointer to the control structure returned
                                                 !! by a previous call to initialize_ice_shelf.

  ! Local variables
                                                 !! various unit conversion factors
                                                 !! the ice-shelf state


end subroutine shelf_calc_flux
module subroutine adjust_ice_sheet_frazil(sfc_state_in, fluxes_in, CS)
  type(surface), target,  intent(inout) :: sfc_state_in !< A structure containing fields that
                                                 !! describe the surface state of the ocean.  The
                                                 !! intent is only inout to allow for halo updates.
  type(forcing),  target, intent(in)    :: fluxes_in !< structure containing pointers to any
                                                 !! possible thermodynamic or mass-flux forcing fields.
  type(ice_shelf_CS),     pointer       :: CS    !< A pointer to the control structure returned
                                                 !! by a previous call to initialize_ice_shelf.
  ! Local variables
                                                 !! the ice-shelf state

end subroutine adjust_ice_sheet_frazil
module function integrate_over_ice_sheet_area(G, ISS, var, unscale, hemisphere) result(var_out)
  type(ocean_grid_type), intent(in) :: G  !< The grid structure used by the ice shelf.
  type(ice_shelf_state), intent(in) :: ISS  !< A structure with elements that describe the ice-shelf state
  real, dimension(SZI_(G),SZJ_(G)), intent(in)  :: var !< Ice variable to integrate in arbitrary units [A ~> a]
  real, intent(in) :: unscale !< Dimensional scaling for variable to integrate [a A-1 ~> 1]
  integer, optional, intent(in) :: hemisphere !< 0 for Antarctica only, 1 for Greenland only. Otherwise, all ice sheets
  real :: var_out !< Variable integrated over the area of the ice sheet in arbitrary scaled units [A L2 ~> a m2]

  ! Local variables
                                                !! in arbitrary units [A L2 ~> a m2]

end function integrate_over_ice_sheet_area
module subroutine ice_sheet_calving_to_ocean_sfc(CS,US,calving,calving_hflx)
  type(ice_shelf_CS),      pointer :: CS        !< A pointer to the ice shelf control structure
  type(unit_scale_type), intent(in)    :: US   !< A dimensional unit scaling type
  real, dimension(:,:), intent(inout) :: calving      !< The mass flux per unit area of the ice shelf
                                                      !! to convert to bergs [R Z T-1 ~> kg m-2 s-1].
  real, dimension(:,:), intent(inout) :: calving_hflx !< Calving heat flux [Q R Z T-1 ~> W m-2].
  ! Local variables
                                                  !! the ice-shelf state

end subroutine ice_sheet_calving_to_ocean_sfc
module subroutine change_thickness_using_melt(ISS, G, US, time_step, fluxes, density_ice, debug)
  type(ocean_grid_type), intent(inout) :: G  !< The ocean's grid structure.
  type(ice_shelf_state), intent(inout) :: ISS !< A structure with elements that describe
                                              !! the ice-shelf state
  type(unit_scale_type), intent(in)    :: US   !< A dimensional unit scaling type
  real,                  intent(in)    :: time_step !< The time step for this update [T ~> s].
  type(forcing),         intent(inout) :: fluxes !< structure containing pointers to any possible
                                                 !! thermodynamic or mass-flux forcing fields.
  real,                  intent(in)    :: density_ice !< The density of ice-shelf ice [R ~> kg m-3].
  logical,     optional, intent(in)    :: debug !< If present and true, write chksums

  ! locals

end subroutine change_thickness_using_melt
module subroutine add_shelf_forces(Ocn_grid, US, CS, forces_in, do_shelf_area, external_call)
  type(ocean_grid_type), intent(in)    :: Ocn_grid !< The ocean's grid structure.
  type(unit_scale_type), intent(in)    :: US   !< A dimensional unit scaling type
  type(ice_shelf_CS),    pointer       :: CS   !< This module's control structure.
  type(mech_forcing), target, intent(inout) :: forces_in !< A structure with the
                                               !! driving mechanical forces
  logical, optional,     intent(in)    :: do_shelf_area !< If true find the shelf-covered areas.
  logical, optional,     intent(in)    :: external_call !< If true the incoming forcing type
                                               !! is using the unrotated input grid and may need
                                               !! to be rotated.
                                          ! the ice-shelf state


end subroutine add_shelf_forces
module subroutine add_shelf_pressure(Ocn_grid, US, CS, fluxes)
  type(ocean_grid_type), intent(in) :: Ocn_grid  !< The ocean's grid structure.
  type(unit_scale_type), intent(in)    :: US     !< A dimensional unit scaling type
  type(ice_shelf_CS),    intent(in)    :: CS     !< This module's control structure.
  type(forcing),         intent(inout) :: fluxes  !< A structure of surface fluxes that may be updated.


end subroutine add_shelf_pressure
module subroutine add_shelf_flux(G, US, CS, sfc_state, fluxes, time_step)
  type(ocean_grid_type), intent(inout) :: G    !< The ocean's grid structure.
  type(unit_scale_type), intent(in)    :: US   !< A dimensional unit scaling type
  type(ice_shelf_CS),    pointer       :: CS   !< This module's control structure.
  type(surface),         intent(inout) :: sfc_state !< Surface ocean state
  type(forcing),         intent(inout) :: fluxes  !< A structure of surface fluxes that may be used/updated.
  real,                  intent(in)    :: time_step !< Time step over which fluxes are applied [T ~> s]
  ! local variables
                          !! balancing the net melt flux occurs, 0 to 1 [nondim]
                          !! at at previous time (Time-dt) [R Z ~> kg m-2]
                          !! the two timesteps at (Time) and (Time-dt) [R Z ~> kg m-2].
                          !! at at previous time (Time-dt)
                          !! at at previous time (Time-dt)
                          !! at at previous time (Time-dt)
                          !! since previous time (Time-dt)
                                          !! the ice-shelf state

end subroutine add_shelf_flux
module subroutine initialize_ice_shelf(param_file, ocn_grid, Time, CS, diag, Time_init, directory, forces_in, &
                                fluxes_in, sfc_state_in, solo_ice_sheet_in, calve_ice_shelf_bergs)
  type(param_file_type),        intent(in)    :: param_file !< A structure to parse for run-time parameters
  type(ocean_grid_type),        pointer       :: ocn_grid   !< The calling ocean model's horizontal grid structure
  type(time_type),              intent(inout) :: Time !< The clock that that will indicate the model time
  type(ice_shelf_CS),           pointer       :: CS   !< A pointer to the ice shelf control structure
  type(MOM_diag_ctrl),          pointer       :: diag !< This is a pointer to the MOM diag CS
                                                      !! which will be discarded
  type(time_type),              intent(in)    :: Time_init !< The time at initialization.
  character(len=*),             intent(in)    :: directory  !< The directory where the energy file goes.

  type(mech_forcing), optional, target, intent(inout) :: forces_in !< A structure with the driving mechanical forces
  type(forcing),      optional, target, intent(inout) :: fluxes_in !< A structure containing pointers to any
                                                           !!  possible thermodynamic or mass-flux forcing fields.
  type(surface), target, optional, intent(inout) :: sfc_state_in !< A structure containing fields that
                                                !! describe the surface state of the ocean.  The
                                                !! intent is only inout to allow for halo updates.
  logical,            optional, intent(in)    :: solo_ice_sheet_in !< If present, this indicates whether
                                                   !! a solo ice-sheet driver.
  logical, optional :: calve_ice_shelf_bergs !< If true, will add point iceberg calving variables to the ice
                                             !! shelf restart

                                                 ! various unit conversion factors
                                          !! the ice-shelf state
                              ! [T kg R-1 Z-1 m-2 s-1 ~> nondim]
                        ! to be floating when CONST_SEA_LEVEL = True [Z ~> m].
  !This include declares and sets the variable "version".
                                   ! does not occur [Z ~> m]


end subroutine initialize_ice_shelf
module subroutine initialize_ice_shelf_fluxes(CS, ocn_grid, US, fluxes_in)
  type(ice_shelf_CS),           pointer       :: CS   !< A pointer to the ice shelf control structure
  type(ocean_grid_type),        pointer       :: ocn_grid   !< The calling ocean model's horizontal grid structure
  type(unit_scale_type),        intent(in)    :: US  !< A dimensional unit scaling type
  type(forcing),        target, intent(inout) :: fluxes_in !< A structure containing pointers to any
                                                           !!  possible thermodynamic or mass-flux forcing fields.

  ! Local variables

end subroutine initialize_ice_shelf_fluxes
module subroutine initialize_ice_shelf_forces(CS, ocn_grid, US, forces_in)
  type(ice_shelf_CS),           pointer       :: CS   !< A pointer to the ice shelf control structure
  type(ocean_grid_type),        pointer       :: ocn_grid   !< The calling ocean model's horizontal grid structure
  type(unit_scale_type),        intent(in)    :: US   !< A dimensional unit scaling type
  type(mech_forcing),   target, intent(inout) :: forces_in !< A structure with the driving mechanical forces

  ! Local variables

end subroutine initialize_ice_shelf_forces
module subroutine initialize_shelf_mass(G, param_file, CS, ISS, new_sim)

  type(ocean_grid_type), intent(in) :: G   !< The ocean's grid structure.
  type(param_file_type), intent(in) :: param_file !< A structure to parse for run-time parameters
  type(ice_shelf_CS),    pointer    :: CS !< A pointer to the ice shelf control structure
  type(ice_shelf_state), intent(inout) :: ISS !< The ice shelf state type that is being updated
  logical,     optional, intent(in) :: new_sim !< If present and false, this run is being restarted

end subroutine initialize_shelf_mass
module subroutine change_thickness_using_precip(CS, ISS, G, US, fluxes, time_step, Time)
  type(ice_shelf_CS),    intent(in)    :: CS  !< A pointer to the ice shelf control structure
  type(ocean_grid_type), intent(inout) :: G  !< The ocean's grid structure.
  type(ice_shelf_state), intent(inout) :: ISS !< A structure with elements that describe
                                              !! the ice-shelf state
  type(forcing),         intent(in)    :: fluxes  !< A structure of surface fluxes that
                                                  !! includes surface mass flux
  type(time_type),       intent(in)    :: Time !< The current model time
  type(unit_scale_type), intent(in)    :: US   !< A dimensional unit scaling type
  real,                  intent(in)    :: time_step !< The time step for this update [T ~> s].

  ! locals

end subroutine change_thickness_using_precip
module subroutine update_shelf_mass(G, US, CS, ISS, Time)
  type(ocean_grid_type), intent(inout) :: G   !< The ocean's grid structure.
  type(unit_scale_type), intent(in)    :: US  !< A dimensional unit scaling type
  type(ice_shelf_CS),    intent(in)    :: CS  !< A pointer to the ice shelf control structure
  type(ice_shelf_state), intent(inout) :: ISS !< The ice shelf state type that is being updated
  type(time_type),       intent(in)    :: Time !< The current model time

  ! local variables

end subroutine update_shelf_mass
module subroutine ice_shelf_query(CS, G, frac_shelf_h, mass_shelf, data_override_shelf_fluxes)
  type(ice_shelf_CS),         pointer    :: CS !< ice shelf control structure
  type(ocean_grid_type), intent(in)      :: G  !< A pointer to an ocean grid control structure.
  real, optional, dimension(SZI_(G),SZJ_(G)), intent(out)  :: frac_shelf_h !< Ice shelf area fraction [nondim].
  real, optional, dimension(SZI_(G),SZJ_(G)), intent(out)  :: mass_shelf !< Ice shelf mass [R Z ~> kg m-2]
  logical, optional                      :: data_override_shelf_fluxes !< If true, shelf fluxes can be written using
                                               !! the data_override capability (only for MOSAIC grids)


end subroutine ice_shelf_query
module subroutine ice_shelf_save_restart(CS, Time, directory, time_stamped, filename_suffix)
  type(ice_shelf_CS),         pointer    :: CS !< ice shelf control structure
  type(time_type),            intent(in) :: Time !< model time at this call
  character(len=*), optional, intent(in) :: directory !< An optional directory into which to write
                                               !! these restart files.
  logical,          optional, intent(in) :: time_stamped !< f true, the restart file names include
                                               !! a unique time stamp.  The default is false.
  character(len=*), optional, intent(in) :: filename_suffix !< An optional suffix (e.g., a
                                               !! time-stamp) to append to the restart file names.
  ! local variables

end subroutine ice_shelf_save_restart
module subroutine ice_shelf_end(CS)
  type(ice_shelf_CS), pointer   :: CS !< A pointer to the ice shelf control structure

end subroutine ice_shelf_end
module subroutine solo_step_ice_shelf(CS, time_interval, nsteps, Time, min_time_step_in, fluxes_in)
  type(ice_shelf_CS), pointer    :: CS      !< A pointer to the ice shelf control structure
  type(time_type), intent(in)    :: time_interval !< The time interval for this update [s].
  integer,         intent(inout) :: nsteps  !< The running number of ice shelf steps.
  type(time_type), intent(inout) :: Time    !< The current model time
  real,  optional, intent(in)    :: min_time_step_in !< The minimum permitted time step [T ~> s].
  type(forcing),      optional, target, intent(inout) :: fluxes_in !< A structure containing pointers to any
                                                         !!  possible thermodynamic or mass-flux forcing fields.
                                                 ! various unit conversion factors
                                          !! the ice-shelf state
                            ! coupled ice-ocean dynamics.
                               !for all ice sheets, Antarctica only, or Greenland only [Z L2 ~> m3]

end subroutine solo_step_ice_shelf
module subroutine process_and_post_scalar_data(CS, vaf0, vaf0_A, vaf0_G, Itime_step, dh_adott, dh_bdott)
  type(ice_shelf_CS), pointer    :: CS      !< A pointer to the ice shelf control structure
  real :: vaf0   !< The previous volumes above floatation for all ice sheets [Z L2 ~> m3]
  real :: vaf0_A !< The previous volumes above floatation for the Antarctic ice sheet [Z L2 ~> m3]
  real :: vaf0_G !< The previous volumes above floatation for the Greenland ice sheet [Z L2 ~> m3]
  real :: Itime_step !< Inverse of the time step [T-1 ~> s-1]
  real, dimension(SZI_(CS%grid),SZJ_(CS%grid)) :: dh_adott !< Surface (plus basal if solo shelf mode)
                               !! melt/accumulation over a time step  [Z ~> m]
  real, dimension(SZI_(CS%grid),SZJ_(CS%grid)) :: dh_bdott !< Surface (plus basal if solo shelf mode)
                               !! melt/accumulation over a time step  [Z ~> m]

  ! Local variables

end subroutine process_and_post_scalar_data
  end interface

end module MOM_ice_shelf
