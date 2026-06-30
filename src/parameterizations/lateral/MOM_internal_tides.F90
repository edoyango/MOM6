! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> Subroutines that use the ray-tracing equations to propagate the internal tide energy density.
!!
!! \author Benjamin Mater & Robert Hallberg, 2015
module MOM_internal_tides

use MOM_checksums,     only : hchksum
use MOM_debugging,     only : is_NaN
use MOM_diag_mediator, only : post_data, query_averaging_enabled, diag_axis_init
use MOM_diag_mediator, only : disable_averaging, enable_averages
use MOM_diag_mediator, only : register_diag_field, diag_ctrl, safe_alloc_ptr
use MOM_diag_mediator, only : axes_grp, define_axes_group
use MOM_domains, only       : AGRID, To_South, To_West, To_All, CGRID_NE
use MOM_domains, only       : create_group_pass, do_group_pass, pass_var, pass_vector
use MOM_domains, only       : group_pass_type, start_group_pass, complete_group_pass
use MOM_error_handler, only : MOM_error, FATAL, WARNING, MOM_mesg, is_root_pe
use MOM_file_parser, only   : read_param, get_param, log_param, log_version, param_file_type
use MOM_forcing_type,only   : forcing
use MOM_grid, only          : ocean_grid_type
use MOM_int_tide_input, only: int_tide_input_CS, get_input_TKE, get_barotropic_tidal_vel
use MOM_io, only            : slasher, MOM_read_data, file_exists, axis_info
use MOM_io, only            : set_axis_info, get_axis_info, stdout
use MOM_restart, only       : register_restart_field, MOM_restart_CS, restart_init, save_restart
use MOM_restart, only       : lock_check, restart_registry_lock
use MOM_spatial_means, only : global_area_integral
use MOM_string_functions, only: extract_real, uppercase
use MOM_time_manager, only  : time_type, time_type_to_real, operator(+), operator(/), operator(-)
use MOM_unit_scaling, only  : unit_scale_type
use MOM_variables, only     : surface, thermo_var_ptrs, vertvisc_type
use MOM_verticalGrid, only  : verticalGrid_type
use MOM_wave_speed, only    : wave_speeds, wave_speed_CS, wave_speed_init
use mpp_domains_mod, only : NORTH_FACE => NORTH, EAST_FACE => EAST

implicit none ; private

#include <MOM_memory.h>

public propagate_int_tide, register_int_tide_restarts
public internal_tides_init, internal_tides_end
public get_lowmode_loss, get_lowmode_diffusivity

!> This control structure has parameters for the MOM_internal_tides module
type, public :: int_tide_CS ; private
  logical :: initialized = .false.   !< True if this control structure has been initialized.
  logical :: do_int_tides    !< If true, use the internal tide code.
  integer :: nFreq = 0       !< The number of internal tide frequency bands
  integer :: nMode = 1       !< The number of internal tide vertical modes
  integer :: nAngle = 24     !< The number of internal tide angular orientations
  integer :: energized_angle = -1 !< If positive, only this angular band is energized for debugging purposes
  real    :: dt_itides       !< The timestep for internal tides ray-tracing [T ~> s]
  real    :: uniform_test_cg !< Uniform group velocity of internal tide
                             !! for testing internal tides [L T-1 ~> m s-1]
  logical :: corner_adv      !< If true, use a corner advection rather than PPM.
  logical :: upwind_1st      !< If true, use a first-order upwind scheme.
  logical :: simple_2nd      !< If true, use a simple second order (arithmetic mean) interpolation
                             !! of the edge values instead of the higher order interpolation.
  logical :: vol_CFL         !< If true, use the ratio of the open face lengths to the tracer cell
                             !! areas when estimating CFL numbers.  Without aggress_adjust,
                             !! the default is false; it is always true with aggress_adjust.
  logical :: use_PPMang      !< If true, use PPM for advection of energy in angular space.
  logical :: update_Kd       !< If true, the scheme will modify the diffusivities seen by the dynamics
  logical :: apply_refraction  !< If false, skip refraction (for debugging)
  logical :: apply_propagation !< If False, do not propagate energy (for debugging)
  logical :: turn_critical_lat !< If True, rays change direction at critical latitude instead
                               !! of being trapped
  logical :: reflect_critical_lat !< If True, rays reflect at the critical latitude instead
                               !! of turning parallel to it
  logical :: debug             !< If true, use debugging prints
  logical :: init_forcing_only !< if True, add TKE forcing only at first step (for debugging)
  logical :: force_posit_En    !< if True, remove subroundoff negative values (needs enhancement)
  logical :: add_tke_forcing = .true. !< Whether to add forcing, used by init_forcing_only

  real, allocatable, dimension(:,:) :: fraction_tidal_input
                        !< how the energy from one tidal component is distributed
                        !! over the various vertical modes, 2d in frequency and mode [nondim]
  real, allocatable, dimension(:,:) :: refl_angle
                        !< local coastline/ridge/shelf angles read from file [rad]
                        ! (could be in G control structure)
  real :: nullangle = -999.9 !< placeholder value in cells with no reflection [rad]
  real, allocatable, dimension(:,:) :: refl_pref
                        !< partial reflection coeff for each "coast cell" [nondim]
                        ! (could be in G control structure)
  logical, allocatable, dimension(:,:) :: refl_pref_logical
                        !< true if reflecting cell with partial reflection
                        ! (could be in G control structure)
  logical, allocatable, dimension(:,:) :: refl_dbl
                        !< identifies reflection cells where double reflection
                        !! is possible (i.e. ridge cells)
                        ! (could be in G control structure)
  real, allocatable, dimension(:,:) :: trans
                        !< partial transmission coeff for each "coast cell" [nondim]
  real, allocatable, dimension(:,:) :: residual
                        !< residual of reflection and transmission coeff for each "coast cell" [nondim]
  real, allocatable, dimension(:,:,:,:) :: cp
                        !< horizontal phase speed [L T-1 ~> m s-1]
  real, allocatable, dimension(:,:,:,:,:) :: TKE_leak_loss
                        !< energy lost due to misc background processes [H Z2 T-3 ~> m3 s-3 or W m-2]
  real, allocatable, dimension(:,:,:,:,:) :: TKE_quad_loss
                        !< energy lost due to quadratic bottom drag [H Z2 T-3 ~> m3 s-3 or W m-2]
  real, allocatable, dimension(:,:,:,:,:) :: TKE_Froude_loss
                        !< energy lost due to wave breaking [H Z2 T-3 ~> m3 s-3 or W m-2]
  real, allocatable, dimension(:,:) :: TKE_itidal_loss_fixed
                        !< Fixed part of the energy lost due to small-scale drag [H Z2 L-2 ~> kg m-2] here.
                        !! This will be multiplied by N and the squared near-bottom velocity (and by
                        !! the near-bottom density in non-Boussinesq mode) to get the energy losses
                        !! in [R Z4 H-1 L-2 ~> kg m-2 or m]
  real, allocatable, dimension(:,:,:,:,:) :: TKE_itidal_loss
                        !< energy lost due to small-scale wave drag [H Z2 T-3 ~> m3 s-3 or W m-2]
  real, allocatable, dimension(:,:,:,:,:) :: TKE_residual_loss
                        !< internal tide energy loss due to the residual at slopes [H Z2 T-3 ~> m3 s-3 or W m-2]
  real, allocatable, dimension(:,:,:,:,:) :: TKE_slope_loss
                        !< internal tide energy loss due to the residual at slopes [H Z2 T-3 ~> m3 s-3 or W m-2]
  real, allocatable, dimension(:,:) :: TKE_input_glo_dt
                        !< The integrated energy input to the internal waves [H Z2 L2 T-2 ~> m5 s-2 or J]
  real, allocatable, dimension(:,:) :: TKE_leak_loss_glo_dt
                        !< Integrated energy lost due to misc background processes [H Z2 L2 T-2 ~> m5 s-2 or J]
  real, allocatable, dimension(:,:) :: TKE_quad_loss_glo_dt
                        !< Integrated energy lost due to quadratic bottom drag [H Z2 L2 T-2 ~> m5 s-2 or J]
  real, allocatable, dimension(:,:) :: TKE_Froude_loss_glo_dt
                        !< Integrated energy lost due to wave breaking [H Z2 L2 T-2 ~> m5 s-2 or J]
  real, allocatable, dimension(:,:) :: TKE_itidal_loss_glo_dt
                        !< energy lost due to small-scale wave drag [H Z2 T-2 ~> m3 s-2 or J m-2]
  real, allocatable, dimension(:,:) :: TKE_residual_loss_glo_dt
                        !< internal tide energy loss due to the residual at slopes [H Z2 L2 T-2 ~> m5 s-2 or J]
  real, allocatable, dimension(:,:) :: error_mode
                        !< internal tide energy budget error for each mode [H Z2 L2 T-2 ~> m5 s-2 or J]
  real, allocatable, dimension(:,:) :: tot_leak_loss !< Energy loss rates due to misc background processes,
                        !! summed over angle, frequency and mode [H Z2 T-3 ~> m3 s-3 or W m-2]
  real, allocatable, dimension(:,:) :: tot_quad_loss !< Energy loss rates due to quadratic bottom drag,
                        !! summed over angle, frequency and mode [H Z2 T-3 ~> m3 s-3 or W m-2]
  real, allocatable, dimension(:,:) :: tot_itidal_loss !< Energy loss rates due to small-scale drag,
                        !! summed over angle, frequency and mode [H Z2 T-3 ~> m3 s-3 or W m-2]
  real, allocatable, dimension(:,:) :: tot_Froude_loss !< Energy loss rates due to wave breaking,
                        !! summed over angle, frequency and mode [H Z2 T-3 ~> m3 s-3 or W m-2]
  real, allocatable, dimension(:,:) :: tot_residual_loss !< Energy loss rates due to residual on slopes,
                        !! summed over angle, frequency and mode [H Z2 T-3 ~> m3 s-3 or W m-2]
  real, allocatable, dimension(:,:) :: tot_allprocesses_loss !< Energy loss rates due to all processes,
                        !! summed over angle, frequency and mode [H Z2 T-3 ~> m3 s-3 or W m-2]
  real, allocatable, dimension(:,:,:,:) :: w_struct !< Vertical structure of vertical velocity (normalized)
                        !! for each frequency and each mode [nondim]
  real, allocatable, dimension(:,:,:,:) :: u_struct !< Vertical structure of horizontal velocity (normalized and
                        !! divided by layer thicknesses) for each frequency and each mode [Z-1 ~> m-1]
  real, allocatable, dimension(:,:,:) :: u_struct_max !< Maximum of u_struct,
                        !! for each mode [Z-1 ~> m-1]
  real, allocatable, dimension(:,:,:) :: u_struct_bot !< Bottom value of u_struct,
                        !! for each mode [Z-1 ~> m-1]
  real, allocatable, dimension(:,:,:) :: int_w2 !< Vertical integral of w_struct squared,
                        !! for each mode [H ~> m or kg m-2]
  real, allocatable, dimension(:,:,:) :: int_U2 !< Vertical integral of u_struct squared,
                        !! for each mode [H Z-2 ~> m-1 or kg m-4]
  real, allocatable, dimension(:,:,:) :: int_N2w2 !< Depth-integrated Brunt Vaissalla freqency times
                        !! vertical profile squared, for each mode [H T-2 ~> m s-2 or kg m-2 s-2]
  real :: q_itides      !< fraction of local dissipation [nondim]
  real :: mixing_effic  !< mixing efficiency [nondim]
  real :: En_sum        !< global sum of energy for use in debugging, in MKS units [m5 s-2 or J]
  real :: En_underflow  !< A minuscule amount of energy [H Z2 T-2 ~> m3 s-2 or J m-2]
  integer :: En_restart_power !< A power factor of 2 by which to multiply the energy in restart [nondim]
  type(time_type), pointer :: Time => NULL() !< A pointer to the model's clock.
  type(group_pass_type) :: pass_En !< Pass 5d array Energy as a group of 3d arrays
  character(len=200) :: inputdir !< directory to look for coastline angle file
  integer :: itides_adv_limiter !< The type of limiter to use for the energy advection scheme
  real, allocatable, dimension(:,:,:,:) :: decay_rate_2d !< rate at which internal tide energy is
                                                         !! lost to the interior ocean internal wave field
                                                         !! as a function of longitude, latitude, frequency
                                                         !! and vertical mode [T-1 ~> s-1].
  real :: cdrag         !< The bottom drag coefficient [nondim].
  real :: drag_min_depth !< The minimum total ocean thickness that will be used in the denominator
                        !! of the quadratic drag terms for internal tides when
                        !! INTERNAL_TIDE_QUAD_DRAG is true [H ~> m or kg m-2]
  real :: gamma_osborn  !< Mixing efficiency from Osborn 1980 [nondim]
  real :: Kd_min        !< The minimum diapycnal diffusivity. [L2 T-1 ~> m2 s-1]
  real :: max_TKE_to_Kd !< Maximum allowed value for TKE_to_kd [H Z2 T-3 ~> m3 s-3 or W m-2]
  real :: min_thick_layer_Kd !< minimum layer thickness allowed to use with TKE_to_kd [H ~> m or kg m-2]
  logical :: apply_background_drag
                        !< If true, apply a drag due to background processes as a sink.
  logical :: apply_bottom_drag
                        !< If true, apply a quadratic bottom drag as a sink.
  logical :: apply_wave_drag
                        !< If true, apply scattering due to small-scale roughness as a sink.
  logical :: apply_Froude_drag
                        !< If true, apply wave breaking as a sink.
  real :: En_check_tol  !< An energy density tolerance for flagging points with small negative
                        !! internal tide energy [H Z2 T-2 ~> m3 s-2 or J m-2]
  logical :: apply_residual_drag
                        !< If true, apply sink from residual term of reflection/transmission.
  logical :: use_2d_decay_rate
                        !< If true, use a spatially varying decay rate for each harmonic.
  real, allocatable :: En(:,:,:,:,:)
                        !< The internal wave energy density as a function of (i,j,angle,frequency,mode)
                        !! integrated within an angular and frequency band [H Z2 T-2 ~> m3 s-2 or J m-2]
  real, allocatable :: En_ini_glo(:,:)
                        !< The internal wave energy density as a function of (frequency,mode) spatially
                        !! integrated within an angular and frequency band [H Z2 L2 T-2 ~> m5 s-2 or J]
                        !! only at the start of the routine (for diags)
  real, allocatable :: En_end_glo(:,:)
                        !< The internal wave energy density as a function of (frequency,mode) spatially
                        !! integrated within an angular and frequency band [H Z2 L2 T-2 ~> m5 s-2 or J]
                        !! only at the end of the routine (for diags)
  real, allocatable :: En_restart_mode1(:,:,:,:)
                        !< The internal wave energy density as a function of (i,j,angle,freq)
                        !! for mode 1 [H Z2 T-2 ~> m3 s-2 or J m-2]
  real, allocatable :: En_restart_mode2(:,:,:,:)
                        !< The internal wave energy density as a function of (i,j,angle,freq)
                        !! for mode 2 [H Z2 T-2 ~> m3 s-2 or J m-2]
  real, allocatable :: En_restart_mode3(:,:,:,:)
                        !< The internal wave energy density as a function of (i,j,angle,freq)
                        !! for mode 3 [H Z2 T-2 ~> m3 s-2 or J m-2]
  real, allocatable :: En_restart_mode4(:,:,:,:)
                        !< The internal wave energy density as a function of (i,j,angle,freq)
                        !! for mode 4 [H Z2 T-2 ~> m3 s-2 or J m-2]
  real, allocatable :: En_restart_mode5(:,:,:,:)
                        !< The internal wave energy density as a function of (i,j,angle,freq)
                        !! for mode 5 [H Z2 T-2 ~> m3 s-2 or J m-2]

  real, allocatable, dimension(:) :: frequency  !< The frequency of each band [T-1 ~> s-1].
  real :: Int_tide_decay_scale  !< vertical decay scale for St Laurent profile [Z ~> m]
  real :: Int_tide_decay_scale_slope  !< vertical decay scale for St Laurent profile on slopes [Z ~> m]

  type(wave_speed_CS) :: wave_speed  !< Wave speed control structure
  type(diag_ctrl), pointer :: diag => NULL() !< A structure that is used to regulate the
                        !! timing of diagnostic output.

  !>@{ Diag handles
  ! Diag handles relevant to all modes, frequencies, and angles
  integer :: id_cg1      = -1                 ! diagnostic handle for mode-1 speed
  integer, allocatable, dimension(:) :: id_cn ! diagnostic handle for all mode speeds
  integer :: id_tot_En = -1
  integer :: id_refl_pref = -1, id_refl_ang = -1, id_land_mask = -1
  integer :: id_trans = -1, id_residual = -1
  integer :: id_dx_Cv = -1, id_dy_Cu = -1
  ! Diag handles considering: sums over all modes, frequencies, and angles
  integer :: id_tot_leak_loss = -1, id_tot_quad_loss = -1, id_tot_itidal_loss = -1
  integer :: id_tot_Froude_loss = -1, id_tot_residual_loss = -1, id_tot_allprocesses_loss = -1
  ! Diag handles considering: all modes & frequencies; summed over angles
  integer, allocatable, dimension(:,:) :: &
             id_En_mode, &
             id_itidal_loss_mode, &
             id_leak_loss_mode, &
             id_quad_loss_mode, &
             id_Froude_loss_mode, &
             id_residual_loss_mode, &
             id_allprocesses_loss_mode, &
             id_itide_drag, &
             id_Ub_mode, &
             id_cp_mode
  ! Diag handles considering: all modes, frequencies, and angles
  integer, allocatable, dimension(:,:) :: &
             id_En_ang_mode, &
             id_itidal_loss_ang_mode
  integer, allocatable, dimension(:) :: &
             id_TKE_itidal_input, &
             id_Ustruct_mode, &
             id_Wstruct_mode, &
             id_int_w2_mode, &
             id_int_U2_mode, &
             id_int_N2w2_mode
  !>@}

end type int_tide_CS

!> A structure with the active energy loop bounds.
type :: loop_bounds_type ; private
  !>@{ The active loop bounds
  integer :: ish, ieh, jsh, jeh
  !>@}
end type loop_bounds_type

!>@{ Enumeration values for numerical schemes
integer, parameter :: LIMITER_ADV_MINMOD = 1
integer, parameter :: LIMITER_ADV_POSITIVE = 2
character*(20), parameter :: LIMITER_ADV_MINMOD_STRING = "MINMOD"
character*(20), parameter :: LIMITER_ADV_POSITIVE_STRING = "POSITIVE"
!>@}


  interface
module subroutine propagate_int_tide(h, tv, Nb, Rho_bot, dt, G, GV, US, inttide_input_CSp, CS)
  type(ocean_grid_type),            intent(inout) :: G  !< The ocean's grid structure.
  type(verticalGrid_type),          intent(in)    :: GV !< The ocean's vertical grid structure.
  type(unit_scale_type),            intent(in)    :: US !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                                    intent(in)    :: h  !< Layer thicknesses [H ~> m or kg m-2]
  type(thermo_var_ptrs),            intent(in)    :: tv !< Pointer to thermodynamic variables
                                                        !! (needed for wave structure).
  real, dimension(SZI_(G),SZJ_(G)), intent(inout) :: Nb !< Near-bottom buoyancy frequency [T-1 ~> s-1].
                                                        !! In some cases the input values are used, but in
                                                        !! others this is set along with the wave speeds.
  real, dimension(SZI_(G),SZJ_(G)), intent(in)    :: Rho_bot !< Near-bottom density or the Boussinesq
                                                        !! reference density [R ~> kg m-3].
  real,                             intent(in)    :: dt !< Length of time over which to advance
                                                        !! the internal tides [T ~> s].
  type(int_tide_input_CS),          intent(in)    :: inttide_input_CSp !< Internal tide input control structure
  type(int_tide_CS),                intent(inout) :: CS !< Internal tide control structure

  ! Local variables

end subroutine propagate_int_tide
module subroutine sum_En(G, GV, US, CS, En, label)
  type(ocean_grid_type),  intent(in) :: G  !< The ocean's grid structure.
  type(verticalGrid_type),intent(in) :: GV !< The ocean's vertical grid structure.
  type(unit_scale_type),  intent(in) :: US !< A dimensional unit scaling type
  type(int_tide_CS),      intent(inout) :: CS !< Internal tide control structure
  real, dimension(G%isd:G%ied,G%jsd:G%jed,CS%NAngle), &
                          intent(in) :: En !< The energy density of the internal tides [H Z2 T-2 ~> m3 s-2 or J m-2].
  character(len=*),       intent(in) :: label !< A label to use in error messages
  ! Local variables
  ! real :: En_sum_diff  ! Change in energy from the expected value [m5 s-2 or J]
  ! real :: En_sum_pdiff ! Percentage change in energy from the expected value [nondim]
  ! character(len=160) :: mesg  ! The text of an error message
  ! real :: days          ! The time in days for use in output messages [days]

end subroutine sum_En
module subroutine itidal_lowmode_loss(G, GV, US, CS, Nb, Rho_bot, Ub, En, TKE_loss_fixed, TKE_loss, dt, halo_size)
  type(ocean_grid_type),     intent(in)    :: G  !< The ocean's grid structure.
  type(verticalGrid_type),   intent(in)    :: GV !< The ocean's vertical grid structure.
  type(unit_scale_type),     intent(in)    :: US !< A dimensional unit scaling type
  type(int_tide_CS),         intent(in)    :: CS !< Internal tide control structure
  real, dimension(G%isd:G%ied,G%jsd:G%jed), &
                             intent(in)    :: Nb !< Near-bottom stratification [T-1 ~> s-1].
  real, dimension(G%isd:G%ied,G%jsd:G%jed), &
                             intent(in)    :: Rho_bot !< Near-bottom density [R ~> kg m-3].
  real, dimension(G%isd:G%ied,G%jsd:G%jed,CS%nFreq,CS%nMode), &
                             intent(inout) :: Ub !< RMS (over one period) near-bottom horizontal
                                                 !! mode velocity [L T-1 ~> m s-1].
  real, dimension(G%isd:G%ied,G%jsd:G%jed), &
                             intent(in) :: TKE_loss_fixed !< Fixed part of energy loss [R Z4 H-1 L-2 ~> kg m-2 or m]
                                                 !! (rho*kappa*h^2) or (kappa*h^2).
  real, dimension(G%isd:G%ied,G%jsd:G%jed,CS%NAngle,CS%nFreq,CS%nMode), &
                             intent(inout) :: En !< Energy density of the internal waves [H Z2 T-2 ~> m3 s-2 or J m-2].
  real, dimension(G%isd:G%ied,G%jsd:G%jed,CS%NAngle,CS%nFreq,CS%nMode), &
                             intent(out)   :: TKE_loss    !< Energy loss rate [H Z2 T-3 ~> m3 s-3 or W m-2]
                                                 !! (q*rho*kappa*h^2*N*U^2).
  real,                      intent(in)    :: dt !< Time increment [T ~> s].
  integer, optional,         intent(in)    :: halo_size !< The halo size over which to do the calculations
  ! Local variables
                             ! and point summed over angles [H Z2 T-2 ~> m3 s-2 or J m-2]
                             ! and point summed over angles [H Z2 T-3 ~> m3 s-3 or W m-2]
                             ! assumed to stay in propagating mode for now - BDM) [nondim]
                             ! units [H Z2 s2 T-2 kg-1 ~> m3 kg-1 or 1]

end subroutine itidal_lowmode_loss
module subroutine get_lowmode_loss(i,j,G,CS,mechanism,TKE_loss_sum)
  integer,               intent(in)  :: i   !< The i-index of the value to be reported.
  integer,               intent(in)  :: j   !< The j-index of the value to be reported.
  type(ocean_grid_type), intent(in)  :: G   !< The ocean's grid structure
  type(int_tide_CS),     intent(in)  :: CS  !< Internal tide control structure
  character(len=*),      intent(in)  :: mechanism    !< The named mechanism of loss to return
  real,                  intent(out) :: TKE_loss_sum !< Total energy loss rate due to specified
                                                     !! mechanism [H Z2 T-3 ~> m3 s-3 or W m-2].

end subroutine get_lowmode_loss
module subroutine get_lowmode_diffusivity(G, GV, h, tv, US, h_bot, k_bot, j, N2_lay, N2_int, TKE_to_Kd, Kd_max, CS, &
                                   Kd_leak, Kd_quad, Kd_itidal, Kd_Froude, Kd_slope, &
                                   Kd_lay, Kd_int, profile_leak, profile_quad, profile_itidal, &
                                   profile_Froude, profile_slope)

  type(ocean_grid_type),               intent(in)  :: G       !< The ocean's grid structure
  type(verticalGrid_type),             intent(in)  :: GV      !< The ocean's vertical grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                                       intent(in)  :: h       !< Layer thicknesses [H ~> m or kg m-2]
  type(thermo_var_ptrs),               intent(in)  :: tv      !< Structure containing pointers to any available
  type(unit_scale_type),               intent(in)  :: US      !< A dimensional unit scaling type
  real, dimension(SZI_(G)),            intent(in)  :: h_bot   !< Bottom boundary layer thickness [H ~> m or kg m-2]
  integer, dimension(SZI_(G)),         intent(in)  :: k_bot   !< Bottom boundary layer top layer index
  integer,                             intent(in)  :: j       !< The j-index to work on
  real, dimension(SZI_(G),SZK_(GV)),   intent(in)  :: N2_lay  !< The squared buoyancy frequency of the
                                                              !! layers [T-2 ~> s-2].
  real, dimension(SZI_(G),SZK_(GV)+1), intent(in)  :: N2_int  !< The squared buoyancy frequency of the
                                                              !! interfaces [T-2 ~> s-2].
  real, dimension(SZI_(G),SZK_(GV)),   intent(in)  :: TKE_to_Kd !< The conversion rate between the TKE
                                                              !! dissipated within a layer and the
                                                              !! diapycnal diffusivity within that layer,
                                                              !! usually (~Rho_0 / (G_Earth * dRho_lay))
                                                              !! [T2 Z-1 ~> s2 m-1]
  real,                                 intent(in) :: Kd_max  !< The maximum increment for diapycnal
                                                              !! diffusivity due to TKE-based processes
                                                              !! [H Z T-1 ~> m2 s-1 or kg m-1 s-1].
                                                              !! Set this to a negative value to have no limit.
                                                              !! [H Z T-1 ~> m2 s-1 or kg m-1 s-1].
  type(int_tide_cs),                    intent(in)    :: CS   !< The control structure for this module

  real, dimension(SZI_(G),SZK_(GV)+1),  intent(out) :: Kd_leak        !< Diffusivity due to background drag
                                                                      !! [H Z T-1 ~> m2 s-1 or kg m-1 s-1].
  real, dimension(SZI_(G),SZK_(GV)+1),  intent(out) :: Kd_quad        !< Diffusivity due to bottom drag
                                                                      !! [H Z T-1 ~> m2 s-1 or kg m-1 s-1].
  real, dimension(SZI_(G),SZK_(GV)+1),  intent(out) :: Kd_itidal      !< Diffusivity due to wave drag
                                                                      !! [H Z T-1 ~> m2 s-1 or kg m-1 s-1].
  real, dimension(SZI_(G),SZK_(GV)+1),  intent(out) :: Kd_Froude      !< Diffusivity due to high Froude breaking
                                                                      !! [H Z T-1 ~> m2 s-1 or kg m-1 s-1].
  real, dimension(SZI_(G),SZK_(GV)+1),  intent(out) :: Kd_slope       !< Diffusivity due to critical slopes
                                                                      !! [H Z T-1 ~> m2 s-1 or kg m-1 s-1].
  real, dimension(SZI_(G),SZK_(GV)),    intent(inout) :: Kd_lay       !< The diapycnal diffusivity in layers
                                                                      !! [H Z T-1 ~> m2 s-1 or kg m-1 s-1].
  real, dimension(SZI_(G),SZK_(GV)+1),  intent(inout) :: Kd_int       !< The diapycnal diffusivity at interfaces
                                                                      !! [H Z T-1 ~> m2 s-1 or kg m-1 s-1].
  real, dimension(SZI_(G), SZK_(GV)),   intent(out) :: profile_leak   !< Normalized profile for background drag
                                                                      !! [H-1 ~> m-1 or m2 kg-1]
  real, dimension(SZI_(G), SZK_(GV)),   intent(out) :: profile_quad   !< Normalized profile for  bottom drag
                                                                      !! [H-1 ~> m-1 or m2 kg-1]
  real, dimension(SZI_(G), SZK_(GV)),   intent(out) :: profile_itidal !< Normalized profile for wave drag
                                                                      !! [H-1 ~> m-1 or m2 kg-1]
  real, dimension(SZI_(G), SZK_(GV)),   intent(out) :: profile_Froude !< Normalized profile for Froude drag
                                                                      !! [H-1 ~> m-1 or m2 kg-1]
  real, dimension(SZI_(G), SZK_(GV)),   intent(out) :: profile_slope  !< Normalized profile for critical slopes
                                                                      !! [H-1 ~> m-1 or m2 kg-1]

  ! local variables

  ! vertical profiles have units Z-1 for conversion to Kd to be dim correct (see eq 2 of St Laurent GRL 2002)
                                                       ! [H-1 ~> m-1 or m2 kg-1]
                                                       ! [H-1 ~> m-1 or m2 kg-1]




end subroutine get_lowmode_diffusivity
module subroutine refract(En, cn, freq, dt, G, US, NAngle, use_PPMang)
  type(ocean_grid_type), intent(in)    :: G    !< The ocean's grid structure.
  integer,               intent(in)    :: NAngle !< The number of wave orientations in the
                                               !! discretized wave energy spectrum.
  real, dimension(G%isd:G%ied,G%jsd:G%jed,NAngle), &
                         intent(inout) :: En   !< The internal gravity wave energy density as a
                                               !! function of space and angular resolution,
                                               !! [H Z2 T-2 ~> m3 s-2 or J m-2].
  real, dimension(G%isd:G%ied,G%jsd:G%jed),        &
                         intent(in)    :: cn   !< Baroclinic mode speed [L T-1 ~> m s-1].
  real,                  intent(in)    :: freq !< Wave frequency [T-1 ~> s-1].
  real,                  intent(in)    :: dt   !< Time step [T ~> s].
  type(unit_scale_type), intent(in)    :: US   !< A dimensional unit scaling type
  logical,               intent(in)    :: use_PPMang !< If true, use PPM for advection rather
                                               !! than upwind.
  ! Local variables
                          ! within a timestep [H Z2 T-2 ~> m3 s-2 or J m-2]

end subroutine refract
module subroutine PPM_angular_advect(En2d, CFL_ang, Flux_En, NAngle, dt, halo_ang)
  integer,                   intent(in)    :: NAngle  !< The number of wave orientations in the
                                                      !! discretized wave energy spectrum [nondim]
  real,                      intent(in)    :: dt      !< Time increment [T ~> s].
  integer,                   intent(in)    :: halo_ang !< The halo size in angular space
  real, dimension(1-halo_ang:NAngle+halo_ang),   &
                             intent(in)    :: En2d    !< The internal gravity wave energy density as a
                                                      !! function of angular resolution [H Z2 T-2 ~> m3 s-2 or J m-2].
  real, dimension(1-halo_ang:NAngle+halo_ang),   &
                             intent(in)    :: CFL_ang !< The CFL number of the energy advection across angles [nondim]
  real, dimension(0:NAngle), intent(out)   :: Flux_En !< The time integrated internal wave energy flux
                                                      !! across angles  [H Z2 T-2 ~> m3 s-2 or J m-2].
  ! Local variables
                       ! orientation [H Z2 T-2 rad-1 ~> m3 s-2 rad-1 or J m-2 rad-1]

end subroutine PPM_angular_advect
module subroutine propagate(En, cn, freq, dt, G, GV, US, CS, NAngle, test, halo_size, residual_loss)
  type(ocean_grid_type), intent(inout) :: G    !< The ocean's grid structure.
  type(verticalGrid_type), intent(in)  :: GV   !< The ocean's vertical grid structure.
  integer,               intent(in)    :: NAngle !< The number of wave orientations in the
                                               !! discretized wave energy spectrum.
  real, dimension(G%isd:G%ied,G%jsd:G%jed,NAngle), &
                         intent(inout) :: En   !< The internal gravity wave energy density as a
                                               !! function of space and angular resolution,
                                               !! [H Z2 T-2 ~> m3 s-2 or J m-2].
  real, dimension(G%isd:G%ied,G%jsd:G%jed),        &
                         intent(in)    :: cn   !< Baroclinic mode speed [L T-1 ~> m s-1].
  real,                  intent(in)    :: freq !< Wave frequency [T-1 ~> s-1].
  real,                  intent(in)    :: dt   !< Time step [T ~> s].
  type(unit_scale_type), intent(in)    :: US   !< A dimensional unit scaling type
  real, dimension(G%isd:G%ied,G%jsd:G%jed,2), intent(in) :: test !< test rotation vector
  type(int_tide_CS),     intent(inout)    :: CS   !< Internal tide control structure
  integer, intent(in) :: halo_size  !< halo size for correct rotation
  real, dimension(G%isd:G%ied,G%jsd:G%jed,NAngle), &
                         intent(inout) :: residual_loss !< internal tide energy loss due
                                                        !! to the residual at slopes [H Z2 T-3 ~> m3 s-3 or W m-2].
  ! Local variables

end subroutine propagate
module subroutine propagate_x(En, speed_x, Cgx_av, dCgx, dt, G, US, Nangle, CS, LB, residual_loss, freq2)
  type(ocean_grid_type),   intent(in)    :: G  !< The ocean's grid structure.
  integer,                 intent(in)    :: NAngle !< The number of wave orientations in the
                                               !! discretized wave energy spectrum.
  real, dimension(G%isd:G%ied,G%jsd:G%jed,Nangle),   &
                           intent(inout) :: En !< The energy density integrated over an angular
                                               !! band [H Z2 T-2 ~> m3 s-2 or J m-2].
  real, dimension(G%IsdB:G%IedB,G%jsd:G%jed),        &
                           intent(in)    :: speed_x !< The magnitude of the group velocity at the
                                               !! Cu points [L T-1 ~> m s-1].
  real, dimension(Nangle), intent(in)    :: Cgx_av !< The average x-projection in each angular band [nondim]
  real, dimension(Nangle), intent(in)    :: dCgx !< The difference in x-projections between the
                                               !! edges of each angular band [nondim].
  real,                    intent(in)    :: dt !< Time increment [T ~> s].
  type(unit_scale_type),   intent(in)    :: US !< A dimensional unit scaling type
  type(int_tide_CS),       intent(in)    :: CS !< Internal tide control structure
  type(loop_bounds_type),  intent(in)    :: LB !< A structure with the active energy loop bounds.
  real, dimension(G%isd:G%ied,G%jsd:G%jed,Nangle),   &
                           intent(inout) :: residual_loss !< internal tide energy loss due
                                                          !! to the residual at slopes [H Z2 T-3 ~> m3 s-3 or W m-2].
  real, intent(in) :: freq2 !< The square of internal tides frequency [T-2 ~> s-2].

  ! Local variables

end subroutine propagate_x
module subroutine propagate_y(En, speed_y, Cgy_av, dCgy, dt, G, US, Nangle, CS, LB, residual_loss, freq2)
  type(ocean_grid_type),   intent(in)    :: G  !< The ocean's grid structure.
  integer,                 intent(in)    :: NAngle !< The number of wave orientations in the
                                               !! discretized wave energy spectrum.
  real, dimension(G%isd:G%ied,G%jsd:G%jed,Nangle), &
                           intent(inout) :: En !< The energy density integrated over an angular
                                               !! band [H Z2 T-2 ~> m3 s-2 or J m-2].
  real, dimension(G%isd:G%ied,G%JsdB:G%JedB),      &
                           intent(in)    :: speed_y !< The magnitude of the group velocity at the
                                               !! Cv points [L T-1 ~> m s-1].
  real, dimension(Nangle), intent(in)    :: Cgy_av !< The average y-projection in each angular band [nondim]
  real, dimension(Nangle), intent(in)    :: dCgy !< The difference in y-projections between the
                                               !! edges of each angular band [nondim]
  real,                    intent(in)    :: dt !< Time increment [T ~> s].
  type(unit_scale_type),   intent(in)    :: US !< A dimensional unit scaling type
  type(int_tide_CS),       intent(in)    :: CS !< Internal tide control structure
  type(loop_bounds_type),  intent(in)    :: LB !< A structure with the active energy loop bounds.
  real, dimension(G%isd:G%ied,G%jsd:G%jed,Nangle),   &
                           intent(inout) :: residual_loss !< internal tide energy loss due
                                                          !! to the residual at slopes [H Z2 T-3 ~> m3 s-3 or W m-2].
  real, intent(in) :: freq2 !< The square of internal tides frequency [T-2 ~> s-2].

  ! Local variables

end subroutine propagate_y
module subroutine zonal_flux_En(u, h, hL, hR, uh, dt, G, US, j, ish, ieh, vol_CFL)
  type(ocean_grid_type),     intent(in)    :: G  !< The ocean's grid structure.
  real, dimension(SZIB_(G)), intent(in)    :: u  !< The zonal velocity [L T-1 ~> m s-1].
  real, dimension(SZI_(G)),  intent(in)    :: h  !< Energy density used to calculate the fluxes
                                                 !! [H Z2 T-2 ~> m3 s-2 or J m-2].
  real, dimension(SZI_(G)),  intent(in)    :: hL !< Left- Energy densities in the reconstruction
                                                 !! [H Z2 T-2 ~> m3 s-2 or J m-2].
  real, dimension(SZI_(G)),  intent(in)    :: hR !< Right- Energy densities in the reconstruction
                                                 !! [H Z2 T-2 ~> m3 s-2 or J m-2].
  real, dimension(SZIB_(G)), intent(out) :: uh !< The zonal energy transport [H Z2 L2 T-3 ~> m5 s-3 or J s-1].
  real,                      intent(in)    :: dt !< Time increment [T ~> s].
  type(unit_scale_type),     intent(in)    :: US !< A dimensional unit scaling type
  integer,                   intent(in)    :: j  !< The j-index to work on.
  integer,                   intent(in)    :: ish !< The start i-index range to work on.
  integer,                   intent(in)    :: ieh !< The end i-index range to work on.
  logical,                   intent(in)    :: vol_CFL !< If true, rescale the ratio of face areas to
                                                 !! the cell areas when estimating the CFL number.
  ! Local variables

end subroutine zonal_flux_En
module subroutine merid_flux_En(v, h, hL, hR, vh, dt, G, US, J, ish, ieh, vol_CFL)
  type(ocean_grid_type),            intent(in)    :: G  !< The ocean's grid structure.
  real, dimension(SZI_(G)),         intent(in)    :: v  !< The meridional velocity [L T-1 ~> m s-1].
  real, dimension(SZI_(G),SZJ_(G)), intent(in)    :: h  !< Energy density used to calculate the
                                                        !! fluxes [H Z2 T-2 ~> m3 s-2 or J m-2].
  real, dimension(SZI_(G),SZJ_(G)), intent(in)    :: hL !< Left- Energy densities in the
                                                        !! reconstruction [H Z2 T-2 ~> m3 s-2 or J m-2].
  real, dimension(SZI_(G),SZJ_(G)), intent(in)    :: hR !< Right- Energy densities in the
                                                        !! reconstruction [H Z2 T-2 ~> m3 s-2 or J m-2].
  real, dimension(SZI_(G)),         intent(out) :: vh !< The meridional energy transport
                                                      !! [H Z2 L2 T-3 ~> m5 s-3 or J s-1].
  real,                             intent(in)    :: dt !< Time increment [T ~> s].
  type(unit_scale_type),            intent(in)    :: US !< A dimensional unit scaling type
  integer,                          intent(in)    :: J  !< The j-index to work on.
  integer,                          intent(in)    :: ish !< The start i-index range to work on.
  integer,                          intent(in)    :: ieh !< The end i-index range to work on.
  logical,                          intent(in)    :: vol_CFL !< If true, rescale the ratio of face
                                                        !! areas to the cell areas when estimating
                                                        !! the CFL number.
  ! Local variables

end subroutine merid_flux_En
module subroutine reflect(En, NAngle, CS, G, LB)
  type(ocean_grid_type),  intent(in)    :: G  !< The ocean's grid structure
  integer,                intent(in)    :: NAngle !< The number of wave orientations in the
                                              !! discretized wave energy spectrum.
  real, dimension(G%isd:G%ied,G%jsd:G%jed,NAngle), &
                          intent(inout) :: En !< The internal gravity wave energy density as a
                                              !! function of space and angular resolution
                                              !! [H Z2 T-2 ~> m3 s-2 or J m-2].
  type(int_tide_CS),      intent(in)    :: CS !< Internal tide control structure
  type(loop_bounds_type), intent(in)    :: LB !< A structure with the active energy loop bounds.

  ! Local variables
                                           ! angle of boundary wrt equator [rad]
                                           ! fraction of wave energy reflected
                                           ! values should collocate with angle_c [nondim]
                                           ! tags of cells with double reflection

                                  ! (values exclude halos)
                                  ! leaving out outdated halo points (march in)

end subroutine reflect
module subroutine turning_latitude(En, NAngle, freq2, CS, G, LB)
  type(ocean_grid_type),  intent(in)    :: G  !< The ocean's grid structure
  integer,                intent(in)    :: NAngle !< The number of wave orientations in the
                                              !! discretized wave energy spectrum.
  real, dimension(G%isd:G%ied,G%jsd:G%jed,NAngle), &
                          intent(inout) :: En !< The internal gravity wave energy density as a
                                              !! function of space and angular resolution
                                              !! [H Z2 T-2 ~> m3 s-2 or J m-2].
  type(int_tide_CS),      intent(in)    :: CS !< Internal tide control structure
  type(loop_bounds_type), intent(in)    :: LB !< A structure with the active energy loop bounds.
  real, intent(in)                      :: freq2 !< The square of the internal tide frequency [T-2 ~> s-2]

  ! Local variables
                                           ! angle of boundary wrt equator [rad]


                                  ! (values exclude halos)
                                  ! leaving out outdated halo points (march in)

end subroutine turning_latitude
module subroutine teleport(En, NAngle, CS, G, LB)
  type(ocean_grid_type),  intent(in)    :: G  !< The ocean's grid structure.
  integer,                intent(in)    :: NAngle !< The number of wave orientations in the
                                              !! discretized wave energy spectrum.
  real, dimension(G%isd:G%ied,G%jsd:G%jed,NAngle), &
                          intent(inout) :: En !< The internal gravity wave energy density as a
                                              !! function of space and angular resolution
                                              !! [H Z2 T-2 ~> m3 s-2 or J m-2].
  type(int_tide_CS),      intent(in)    :: CS !< Internal tide control structure
  type(loop_bounds_type), intent(in)    :: LB !< A structure with the active energy loop bounds.
  ! Local variables
                                              ! angle of boundary wrt equator [rad]
                                              ! fraction of wave energy reflected
                                              ! values should collocate with angle_c [nondim]
                                              ! flag for partial reflection
                                              ! tags of cells with double reflection
                                    ! leaving out outdated halo points (march in)

end subroutine teleport
module subroutine correct_halo_rotation(En, test, G, NAngle, halo)
  type(ocean_grid_type),      intent(in)    :: G    !< The ocean's grid structure
  real, dimension(:,:,:,:,:), intent(inout) :: En   !< The internal gravity wave energy density as a
                                       !! function of space, angular orientation, frequency,
                                       !! and vertical mode [H Z2 T-2 ~> m3 s-2 or J m-2].
  real, dimension(SZI_(G),SZJ_(G),2), &
                              intent(in)    :: test !< An x-unit vector that has been passed through
                                       !! the halo updates, to enable the rotation of the
                                       !! wave energies in the halo region to be corrected [nondim].
  integer,                    intent(in)    :: NAngle !< The number of wave orientations in the
                                                      !! discretized wave energy spectrum.
  integer,                    intent(in)    :: halo   !< The halo size over which to do the calculations
  ! Local variables
                                              ! in a frequency band and mode [H Z2 T-2 ~> m3 s-2 or J m-2].
end subroutine correct_halo_rotation
module subroutine correct_halo_rotation_2d(En, test, G, NAngle, halo)
  type(ocean_grid_type),      intent(in)    :: G    !< The ocean's grid structure
  real, dimension(:,:,:), intent(inout) :: En   !< The internal gravity wave energy density as a
                                       !! function of space, angular orientation, frequency,
                                       !! and vertical mode [H Z2 T-2 ~> m3 s-2 or J m-2].
  real, dimension(SZI_(G),SZJ_(G),2), &
                              intent(in)    :: test !< An x-unit vector that has been passed through
                                       !! the halo updates, to enable the rotation of the
                                       !! wave energies in the halo region to be corrected [nondim].
  integer,                    intent(in)    :: NAngle !< The number of wave orientations in the
                                                      !! discretized wave energy spectrum.
  integer,                    intent(in)    :: halo   !< The halo size over which to do the calculations
  ! Local variables
                                              ! in a frequency band and mode [H Z2 T-2 ~> m3 s-2 or J m-2].
end subroutine correct_halo_rotation_2d
module subroutine PPM_reconstruction_x(h_in, h_l, h_r, G, LB, simple_2nd, adv_limiter)
  type(ocean_grid_type),            intent(in)  :: G    !< The ocean's grid structure.
  real, dimension(SZI_(G),SZJ_(G)), intent(in)  :: h_in !< Energy density in a sector (2D)
                                                        !! [H Z2 T-2 ~> m3 s-2 or J m-2]
  real, dimension(SZI_(G),SZJ_(G)), intent(out) :: h_l  !< Left edge value of reconstruction (2D)
                                                        !! [H Z2 T-2 ~> m3 s-2 or J m-2]
  real, dimension(SZI_(G),SZJ_(G)), intent(out) :: h_r  !< Right edge value of reconstruction (2D)
                                                        !! [H Z2 T-2 ~> m3 s-2 or J m-2]
  type(loop_bounds_type),           intent(in)  :: LB   !< A structure with the active loop bounds.
  logical,                          intent(in)  :: simple_2nd !< If true, use the arithmetic mean
                                                        !! energy densities as default edge values
                                                        !! for a simple 2nd order scheme.
  integer,                          intent(in)  :: adv_limiter !< The type of limiter used

  ! Local variables
                                           ! [H Z2 T-2 ~> m3 s-2 or J m-2]
                   ! relative to the center point [H Z2 T-2 ~> m3 s-2 or J m-2]

end subroutine PPM_reconstruction_x
module subroutine PPM_reconstruction_y(h_in, h_l, h_r, G, LB, simple_2nd, adv_limiter)
  type(ocean_grid_type),            intent(in)  :: G    !< The ocean's grid structure.
  real, dimension(SZI_(G),SZJ_(G)), intent(in)  :: h_in !< Energy density in a sector (2D)
                                                        !! [H Z2 T-2 ~> m3 s-2 or J m-2]
  real, dimension(SZI_(G),SZJ_(G)), intent(out) :: h_l  !< Left edge value of reconstruction (2D)
                                                        !! [H Z2 T-2 ~> m3 s-2 or J m-2]
  real, dimension(SZI_(G),SZJ_(G)), intent(out) :: h_r  !< Right edge value of reconstruction (2D)
                                                        !! [H Z2 T-2 ~> m3 s-2 or J m-2]
  type(loop_bounds_type),           intent(in)  :: LB   !< A structure with the active loop bounds.
  logical,                          intent(in)  :: simple_2nd !< If true, use the arithmetic mean
                                                        !! energy densities as default edge values
                                                        !! for a simple 2nd order scheme.
  integer,                          intent(in)  :: adv_limiter !< The type of limiter used

  ! Local variables
                                           ! [H Z2 T-2 ~> m3 s-2 or J m-2]
                   ! relative to the center point [H Z2 T-2 ~> m3 s-2 or J m-2]

end subroutine PPM_reconstruction_y
module subroutine PPM_limit_pos(h_in, h_L, h_R, h_min, G, iis, iie, jis, jie)
  type(ocean_grid_type),            intent(in)     :: G     !< The ocean's grid structure.
  real, dimension(SZI_(G),SZJ_(G)), intent(in)     :: h_in  !< Energy density in each sector (2D)
                                                            !! [H Z2 T-2 ~> m3 s-2 or J m-2]
  real, dimension(SZI_(G),SZJ_(G)), intent(inout)  :: h_L   !< Left edge value of reconstruction
                                                            !!  [H Z2 T-2 ~> m3 s-2 or J m-2]
  real, dimension(SZI_(G),SZJ_(G)), intent(inout)  :: h_R   !< Right edge value of reconstruction
                                                            !! [H Z2 T-2 ~> m3 s-2 or J m-2]
  real,                             intent(in)     :: h_min !< The minimum value that can be
                                                            !! obtained by a concave parabolic fit
                                                            !! [H Z2 T-2 ~> m3 s-2 or J m-2]
  integer,                          intent(in)     :: iis   !< Start i-index for computations
  integer,                          intent(in)     :: iie   !< End i-index for computations
  integer,                          intent(in)     :: jis   !< Start j-index for computations
  integer,                          intent(in)     :: jie   !< End j-index for computations
  ! Local variables

end subroutine PPM_limit_pos
module subroutine minmod_limiter(h_in, h_L, h_R, G, iis, iie, jis, jie)
  type(ocean_grid_type),            intent(in)     :: G     !< The ocean's grid structure.
  real, dimension(SZI_(G),SZJ_(G)), intent(in)     :: h_in  !< Energy density in each sector (2D)
                                                            !! [H Z2 T-2 ~> m3 s-2 or J m-2]
  real, dimension(SZI_(G),SZJ_(G)), intent(inout)  :: h_L   !< Left edge value of reconstruction
                                                            !!  [H Z2 T-2 ~> m3 s-2 or J m-2]
  real, dimension(SZI_(G),SZJ_(G)), intent(inout)  :: h_R   !< Right edge value of reconstruction
                                                            !! [H Z2 T-2 ~> m3 s-2 or J m-2]
  integer,                          intent(in)     :: iis   !< Start i-index for computations
  integer,                          intent(in)     :: iie   !< End i-index for computations
  integer,                          intent(in)     :: jis   !< Start j-index for computations
  integer,                          intent(in)     :: jie   !< End j-index for computations
  ! Local variables

end subroutine minmod_limiter
module subroutine register_int_tide_restarts(G, GV, US, param_file, CS, restart_CS)
  type(ocean_grid_type), intent(in) :: G          !< The ocean's grid structure
  type(verticalGrid_type),intent(in):: GV         !< The ocean's vertical grid structure.
  type(unit_scale_type), intent(in) :: US         !< A dimensional unit scaling type
  type(param_file_type), intent(in) :: param_file !< A structure to parse for run-time parameters
  type(int_tide_CS),     pointer    :: CS         !< Internal tide control structure
  type(MOM_restart_CS),  pointer    :: restart_CS !< MOM restart control structure

  ! This subroutine is used to allocate and register any fields in this module
  ! that should be written to or read from the restart file.


end subroutine register_int_tide_restarts
module subroutine internal_tides_init(Time, G, GV, US, param_file, diag, CS)
  type(time_type), target,   intent(in)    :: Time !< The current model time.
  type(ocean_grid_type),     intent(in)    :: G    !< The ocean's grid structure.
  type(verticalGrid_type),   intent(in)    :: GV   !< The ocean's vertical grid structure.
  type(unit_scale_type),     intent(in)    :: US   !< A dimensional unit scaling type
  type(param_file_type),     intent(in)    :: param_file !< A structure to parse for run-time
                                                   !! parameters.
  type(diag_ctrl), target,   intent(in)    :: diag !< A structure that is used to regulate
                                                   !! diagnostic output.
  type(int_tide_CS),         pointer :: CS         !< Internal tide control structure

  ! Local variables
                                                  ! of cells with double-reflecting ridges [nondim]
                                                 ! lost to the interior ocean internal wave field [T-1 ~> s-1].
                 ! mode speeds are not calculated but simply assigned a speed of 0 [L T-1 ~> m s-1].
                                ! nominal ocean depth, or a negative value for no limit [nondim]
                                ! to mks [T2 kg H-1 Z-2 s-2 ~> kg m-3 or 1]
                                ! to mks [T3 kg H-1 Z-2 s-3 ~> kg m-3 or 1]
                                ! units [H Z2 s2 T-2 kg-1 ~> m3 kg-1 or 1]
  ! This include declares and sets the variable "version".


end subroutine internal_tides_init
module subroutine internal_tides_end(CS)
  type(int_tide_CS), intent(inout) :: CS  !<  Internal tide control structure

end subroutine internal_tides_end
  end interface

end module MOM_internal_tides
