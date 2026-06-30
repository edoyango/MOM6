! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> Interface for surface waves
module MOM_wave_interface

use MOM_data_override, only : data_override_init, data_override
use MOM_diag_mediator, only : post_data, register_diag_field, safe_alloc_alloc
use MOM_diag_mediator, only : diag_ctrl
use MOM_domains,       only : pass_var, pass_vector, AGRID
use MOM_domains,       only : To_South, To_West, To_All
use MOM_error_handler, only : MOM_error, FATAL, WARNING
use MOM_file_parser,   only : get_param, log_version, param_file_type
use MOM_forcing_type,  only : mech_forcing
use MOM_grid,          only : ocean_grid_type
use MOM_hor_index,     only : hor_index_type
use MOM_io,            only : file_exists, get_var_sizes, read_variable
use MOM_io,            only : vardesc, var_desc
use MOM_safe_alloc,    only : safe_alloc_ptr
use MOM_spatial_means, only : global_area_mean
use MOM_time_manager,  only : time_type, operator(+), operator(/)
use MOM_unit_scaling,  only : unit_scale_type
use MOM_variables,     only : thermo_var_ptrs, surface
use MOM_verticalgrid,  only : verticalGrid_type
use MOM_restart,       only : register_restart_pair, MOM_restart_CS

implicit none ; private

#include <MOM_memory.h>

public MOM_wave_interface_init ! Public interface to fully initialize the wave routines.
public query_wave_properties ! Public interface to obtain information from the waves control structure.
public Update_Surface_Waves ! Public interface to update wave information at the
                            ! coupler/driver level.
public Update_Stokes_Drift ! Public interface to update the Stokes drift profiles
                           ! called in step_mom.
public get_Langmuir_Number ! Public interface to compute Langmuir number called from
                           ! ePBL or KPP routines.
public Stokes_PGF ! Public interface to compute Stokes-shear induced pressure gradient force anomaly
public StokesMixing ! NOT READY - Public interface to add down-Stokes gradient
                    ! momentum mixing (e.g. the approach of Harcourt 2013/2015)
public CoriolisStokes ! NOT READY - Public interface to add Coriolis-Stokes acceleration
                      ! of the mean currents, needed for comparison with LES.  It is
                      ! presently advised against implementing in non-1d settings without
                      ! serious consideration of the full 3d wave-averaged Navier-Stokes
                      ! CL2 effects.
public Waves_end ! public interface to deallocate and free wave related memory.
public get_wave_method ! public interface to obtain the wave method string
public waves_register_restarts ! public interface to register wave restart fields

! A note on unit descriptions in comments: MOM6 uses units that can be rescaled for dimensional
! consistency testing. These are noted in comments with units like Z, H, L, and T, along with
! their mks counterparts with notation like "a velocity [Z T-1 ~> m s-1]".  If the units
! vary with the Boussinesq approximation, the Boussinesq variant is given first.

!> Container for all surface wave related parameters
type, public :: wave_parameters_CS ; private

  ! Main surface wave options and publicly visible variables
  logical, public :: UseWaves = .false.     !< Flag to enable surface gravity wave feature
  logical, public :: Stokes_VF = .false.    !< True if Stokes vortex force is used
  logical, public :: Passive_Stokes_VF = .false. !< Computes Stokes VF, but doesn't affect dynamics
  logical, public :: Stokes_PGF = .false.   !< True if Stokes shear pressure Gradient force is used
  logical, public :: robust_Stokes_PGF = .false.  !< If true, use expressions to calculate the
                                            !! Stokes-induced pressure gradient anomalies that are
                                            !! more accurate in the limit of thin layers.
  logical, public :: Passive_Stokes_PGF = .false. !< Keeps Stokes_PGF on, but doesn't affect dynamics
  logical, public :: Stokes_DDT = .false.   !< Developmental:
                                            !! True if Stokes d/dt is used
  logical, public :: Passive_Stokes_DDT = .false.   !< Keeps Stokes_DDT on, but doesn't affect dynamics
  logical :: Homogenize_Surfbands !< True to homogenize surface band Stokes drift in the horizontal

  real, allocatable, dimension(:,:,:), public :: &
    Us_x               !< 3d zonal Stokes drift profile [L T-1 ~> m s-1]
                       !! Horizontal -> U points
                       !! Vertical -> Mid-points
  real, allocatable, dimension(:,:,:), public :: &
    Us_y               !< 3d meridional Stokes drift profile [L T-1 ~> m s-1]
                       !! Horizontal -> V points
                       !! Vertical -> Mid-points
  real, allocatable, dimension(:,:,:), public :: &
    ddt_Us_x           !< 3d time tendency of zonal Stokes drift profile [L T-2 ~> m s-2]
                       !! Horizontal -> U points
                       !! Vertical -> Mid-points
  real, allocatable, dimension(:,:,:), public :: &
    ddt_Us_y           !< 3d time tendency of meridional Stokes drift profile [L T-2 ~> m s-2]
                       !! Horizontal -> V points
                       !! Vertical -> Mid-points
  real, allocatable, dimension(:,:,:), public :: &
    Us_x_from_ddt      !< Check of 3d zonal Stokes drift profile [L T-1 ~> m s-1]
                       !! Horizontal -> U points
                       !! Vertical -> Mid-points
  real, allocatable, dimension(:,:,:), public :: &
    Us_y_from_ddt      !< Check of 3d meridional Stokes drift profile [L T-1 ~> m s-1]
                       !! Horizontal -> V points
                       !! Vertical -> Mid-points
  real, allocatable, dimension(:,:,:), public :: &
    Us_x_prev          !< 3d zonal Stokes drift profile, previous dynamics call [L T-1 ~> m s-1]
                       !! Horizontal -> U points
                       !! Vertical -> Mid-points
  real, allocatable, dimension(:,:,:), public :: &
    Us_y_prev          !< 3d meridional Stokes drift profile, previous dynamics call [L T-1 ~> m s-1]
                       !! Horizontal -> V points
                       !! Vertical -> Mid-points
  real, allocatable, dimension(:,:,:), public :: &
    KvS                !< Viscosity for Stokes Drift shear [H Z T-1 ~> m2 s-1 or Pa s]
  real, allocatable, dimension(:), public :: &
    WaveNum_Cen        !< Wavenumber bands for read/coupled [Z-1 ~> m-1]
  real, allocatable, dimension(:,:,:), public :: &
    UStk_Hb            !< Surface Stokes Drift spectrum (zonal) [L T-1 ~> m s-1]
                       !! Horizontal -> H-points
                       !! 3rd dimension -> Freq/Wavenumber
  real, allocatable, dimension(:,:,:), public :: &
    VStk_Hb            !< Surface Stokes Drift spectrum (meridional) [L T-1 ~> m s-1]
                       !! Horizontal -> H-points
                       !! 3rd dimension -> Freq/Wavenumber
  real, allocatable, dimension(:,:), public :: &
    Omega_w2x          !< wind direction ccw from model x- axis   [nondim radians]
  integer, public :: NumBands = 0   !< Number of wavenumber/frequency partitions
                                    !! Must match the number of bands provided
                                    !! via either coupling or file.

  ! The remainder of this control structure is private
  integer :: WaveMethod = -99 !< Options for including wave information
                              !! Valid (tested) choices are:
                              !!   0 - Test Profile
                              !!   1 - Surface Stokes Drift Bands
                              !!   2 - DHH85
                              !!   3 - LF17
                              !! -99 - No waves computed, but empirical Langmuir number used.
  logical :: LagrangianMixing !< This feature is in development and not ready
                              !! True if Stokes drift is present and mixing
                              !! should be applied to Lagrangian current
                              !! (mean current + Stokes drift).
                              !! See Reichl et al., 2016 KPP-LT approach
  logical :: StokesMixing     !< This feature is in development and not ready.
                              !! True if vertical mixing of momentum
                              !! should be applied directly to Stokes current
                              !! (with separate mixing parameter for Eulerian
                              !! mixing contribution).
                              !! See Harcourt 2013, 2015 Second-Moment approach
  logical :: CoriolisStokes   !< This feature is in development and not ready.
                              ! True if Coriolis-Stokes acceleration should be applied.
  real :: Stokes_min_thick_avg !< A layer thickness below which the cell-center Stokes drift is
                              !! used instead of the cell average [Z ~> m].  This is only used if
                              !! WAVE_INTERFACE_ANSWER_DATE < 20230101.
  integer :: answer_date      !< The vintage of the order of arithmetic and expressions in the
                              !! surface wave calculations.  Values below 20230101 recover the
                              !! answers from the end of 2022, while higher values use updated
                              !! and more robust forms of the same expressions.

  ! Options if WaveMethod is Surface Stokes Drift Bands (1)
  integer :: PartitionMode  !< Method for partition mode (meant to check input)
                            !! 0 - wavenumbers
                            !! 1 - frequencies
  integer :: DataSource !< Integer that specifies where the model Looks for data
                        !! Valid choices are:
                        !! 1 - FMS DataOverride Routine
                        !! 2 - Reserved For Coupler
                        !! 3 - User input (fixed values, useful for 1d testing)

  ! Options if using FMS DataOverride Routine
  character(len=40)  :: SurfBandFileName !< Filename if using DataOverride
  real :: land_speed    !< A large Stokes velocity that can be used to indicate land values in
                        !! a data override file [L T-1 ~> m s-1].  Stokes drift components larger
                        !! than this are set to zero in data override calls for the Stokes drift.
  logical :: DataOver_initialized !< Flag for DataOverride Initialization

  ! Options for computing Langmuir number
  real :: LA_FracHBL         !< Fraction of OSBL for averaging Langmuir number [nondim]
  real :: LA_HBL_min         !< Minimum boundary layer depth for averaging Langmuir number [Z ~> m]
  logical :: LA_Misalignment = .false. !< Flag to use misalignment in Langmuir number
  logical :: LA_misalign_bug = .false. !< Flag to use code with a sign error when calculating the
                       !! misalignment between the shear and waves in the Langmuir number calculation.
  real :: g_Earth      !< The gravitational acceleration, equivalent to GV%g_Earth but with
                       !! different dimensional rescaling appropriate for deep-water gravity
                       !! waves [Z T-2 ~> m s-2]
  real :: I_g_Earth    !< The inverse of the gravitational acceleration, with dimensional rescaling
                       !! appropriate for deep-water gravity waves [T2 Z-1 ~> s2 m-1]
  ! Surface Wave Dependent 1d/2d/3d vars
  real, allocatable, dimension(:) :: &
    Freq_Cen           !< Central frequency for wave bands, including a factor of 2*pi [T-1 ~> s-1]
  real, allocatable, dimension(:) :: &
    PrescribedSurfStkX !< Surface Stokes drift if prescribed [L T-1 ~> m s-1]
  real, allocatable, dimension(:) :: &
    PrescribedSurfStkY !< Surface Stokes drift if prescribed [L T-1 ~> m s-1]
  real, allocatable, dimension(:,:) :: &
    La_Turb            !< Aligned Turbulent Langmuir number [nondim]
                       !! Horizontal -> H points
  real, allocatable, dimension(:,:) :: &
    US0_x              !< Surface Stokes Drift (zonal) [L T-1 ~> m s-1]
                       !! Horizontal -> U points
  real, allocatable, dimension(:,:) :: &
    US0_y              !< Surface Stokes Drift (meridional) [L T-1 ~> m s-1]
                       !! Horizontal -> V points
  real, allocatable, dimension(:,:,:) :: &
    STKx0              !< Stokes Drift spectrum (zonal) [L T-1 ~> m s-1]
                       !! Horizontal -> U points
                       !! 3rd dimension -> Freq/Wavenumber
  real, allocatable, dimension(:,:,:) :: &
    STKy0              !< Stokes Drift spectrum (meridional) [L T-1 ~> m s-1]
                       !! Horizontal -> V points
                       !! 3rd dimension -> Freq/Wavenumber

  real :: La_min       !< An arbitrary lower-bound on the Langmuir number [nondim].
                       !! Langmuir number is sqrt(u_star/u_stokes).  When both are small
                       !! but u_star is orders of magnitude smaller, the Langmuir number could
                       !! have unintended consequences.  Since both are small it can be safely
                       !! capped to avoid such consequences.
  real :: La_Stk_backgnd !< A small background Stokes velocity used in the denominator of
                       !! some expressions for the Langmuir number [L T-1 ~> m s-1]

  ! Parameters used in estimating the wind speed or wave properties from the friction velocity
  real :: VonKar = -1.0 !< The von Karman coefficient as used in the MOM_wave_interface module [nondim]
  real :: rho_air  !< A typical density of air at sea level, as used in wave calculations [R ~> kg m-3]
  real :: nu_air   !< The viscosity of air, as used in wave calculations [Z2 T-1 ~> m2 s-1]
  real :: rho_ocn  !< A typical surface density of seawater, as used in wave calculations in
                   !! comparison with the density of air [R ~> kg m-3].  The default is RHO_0.
  real :: SWH_from_u10sq !< A factor for converting the square of the 10 m wind speed to the
                   !! significant wave height [Z T2 L-2 ~> s2 m-1]
  real :: Charnock_min !< The minimum value of the Charnock coefficient, which relates the square of
                   !! the air friction velocity divided by the gravitational acceleration to the
                   !! wave roughness length [nondim]
  real :: Charnock_slope_U10 !< The partial derivative of the Charnock coefficient with the 10 m wind
                   !! speed [T L-1 ~> s m-1].   Note that in eq. 13 of the Edson et al. 2013 describing
                   !! the COARE 3.5 bulk flux algorithm, this slope is given as 0.017.  However, 0.0017
                   !! reproduces the curve in their figure 6, so that is the default value used in MOM6.
  real :: Charnock_intercept !< The intercept of the fit for the Charnock coefficient in the limit of
                   !! no wind [nondim].  Note that this can be negative because CHARNOCK_MIN will keep
                   !! the final value for the Charnock coefficient from being from being negative.

  ! Options used with the test profile
  real    :: TP_STKX0     !< Test profile x-stokes drift amplitude [L T-1 ~> m s-1]
  real    :: TP_STKY0     !< Test profile y-stokes drift amplitude [L T-1 ~> m s-1]
  real    :: TP_WVL       !< Test profile wavelength [Z ~> m]

  ! Options for use with the Donelan et al., 1985 (DHH85) spectrum
  logical :: WaveAgePeakFreq !< Flag to use wave age to determine the peak frequency with DHH85
  logical :: StaticWaves  !< Flag to disable updating DHH85 Stokes drift
  logical :: DHH85_is_set !< The if the wave properties have been set when WaveMethod = DHH85.
  real    :: WaveAge      !< The fixed wave age used with the DHH85 spectrum [nondim]
  real    :: WaveWind     !< Wind speed for the DHH85 spectrum [L T-1 ~> m s-1]
  real    :: omega_min    !< Minimum wave frequency with the DHH85 spectrum [T-1 ~> s-1]
  real    :: omega_max    !< Maximum wave frequency with the DHH85 spectrum [T-1 ~> s-1]

  type(time_type), pointer :: Time !< A pointer to the ocean model's clock.
  type(diag_ctrl), pointer :: diag !< A structure that is used to regulate the
                                   !! timing of diagnostic output.

  !>@{ Diagnostic handles
  integer, public :: id_PFu_Stokes = -1 , id_PFv_Stokes = -1
  integer, public :: id_3dstokes_x_from_ddt = -1 , id_3dstokes_y_from_ddt = -1
  integer :: id_P_deltaStokes_L = -1, id_P_deltaStokes_i = -1
  integer :: id_surfacestokes_x = -1 , id_surfacestokes_y = -1
  integer :: id_3dstokes_x = -1 , id_3dstokes_y = -1
  integer :: id_ddt_3dstokes_x = -1 , id_ddt_3dstokes_y = -1
  integer :: id_La_turb = -1
  !>@}

end type wave_parameters_CS

! Switches needed in import_stokes_drift
!>@{ Enumeration values for the wave method
integer, parameter :: TESTPROF = 0, SURFBANDS = 1, DHH85 = 2, LF17 = 3, EFACTOR = 4, NULL_WaveMethod = -99
!>@}
!>@{ Enumeration values for the wave data source
integer, parameter :: DATAOVR = 1, COUPLER = 2, INPUT = 3
!>@}

! Strings for the wave method
character*(5), parameter  :: NULL_STRING      = "EMPTY"         !< null wave method string
character*(12), parameter :: TESTPROF_STRING  = "TEST_PROFILE"  !< test profile string
character*(13), parameter :: SURFBANDS_STRING = "SURFACE_BANDS" !< surface bands string
character*(5), parameter  :: DHH85_STRING     = "DHH85"         !< DHH85 wave method string
character*(4), parameter  :: LF17_STRING      = "LF17"          !< LF17 wave method string
character*(7), parameter  :: EFACTOR_STRING   = "EFACTOR"       !< EFACTOR (based on vr12-ma) wave method string


  interface
module subroutine MOM_wave_interface_init(time, G, GV, US, param_file, CS, diag)
  type(time_type), target, intent(in)    :: Time       !< Model time
  type(ocean_grid_type),   intent(inout) :: G          !< Grid structure
  type(verticalGrid_type), intent(in)    :: GV         !< Vertical grid structure
  type(unit_scale_type),   intent(in)    :: US         !< A dimensional unit scaling type
  type(param_file_type),   intent(in)    :: param_file !< Input parameter structure
  type(wave_parameters_CS), pointer      :: CS         !< Wave parameter control structure
  type(diag_ctrl), target, intent(inout) :: diag       !< Diagnostic Pointer

  ! Local variables
  ! This include declares and sets the variable "version".

  ! Dummy Check
end subroutine MOM_wave_interface_init
module subroutine set_LF17_wave_params(param_file, mdl, GV, US, CS)
  type(param_file_type),   intent(in)    :: param_file !< Input parameter structure
  character(len=*),        intent(in)    :: mdl        !< A module name to use in the get_param calls
  type(verticalGrid_type), intent(in)    :: GV         !< Vertical grid structure
  type(unit_scale_type),   intent(in)    :: US         !< A dimensional unit scaling type
  type(wave_parameters_CS), pointer      :: CS         !< Wave parameter control structure

  ! A separate routine is used to set these parameters because there are multiple ways that the
  ! underlying parameterizations are enabled.

end subroutine set_LF17_wave_params
module subroutine query_wave_properties(CS, NumBands, WaveNumbers, US)
  type(wave_parameters_CS),        pointer     :: CS   !< Wave parameter Control structure
  integer,               optional, intent(out) :: NumBands    !< If present, this returns the number of
                                                       !!< wavenumber partitions in the wave discretization
  real, dimension(:),    optional, intent(out) :: Wavenumbers !< If present this returns the characteristic
                                                       !! wavenumbers of the wave discretization [m-1] or [Z-1 ~> m-1]
  type(unit_scale_type), optional, intent(in)  :: US   !< A dimensional unit scaling type that is used to undo
                                                       !! the dimensional scaling of the output variables, if present

end subroutine query_wave_properties
module subroutine Update_Surface_Waves(G, GV, US, Time_present, dt, CS, forces)
  type(wave_parameters_CS), pointer    :: CS  !< Wave parameter Control structure
  type(ocean_grid_type), intent(inout) :: G   !< Grid structure
  type(verticalGrid_type), intent(in)  :: GV  !< Vertical grid structure
  type(unit_scale_type),   intent(in)  :: US  !< A dimensional unit scaling type
  type(time_type),         intent(in)  :: Time_present !< Model Time
  type(time_type),         intent(in)  :: dt  !< Time increment as a time-type
  type(mech_forcing),      intent(in), optional  :: forces !< MOM_forcing_type
  ! Local variables

end subroutine Update_Surface_Waves
module subroutine Update_Stokes_Drift(G, GV, US, CS, dz, ustar, dt, dynamics_step)
  type(wave_parameters_CS), pointer       :: CS    !< Wave parameter Control structure
  type(ocean_grid_type),    intent(inout) :: G     !< Grid structure
  type(verticalGrid_type),  intent(in)    :: GV    !< Vertical grid structure
  type(unit_scale_type),    intent(in)    :: US    !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                            intent(in)    :: dz    !< Thickness in height units [Z ~> m]
  real, dimension(SZI_(G),SZJ_(G)), &
                            intent(in)    :: ustar !< Wind friction velocity [Z T-1 ~> m s-1].
  real, intent(in)                        :: dt    !< Time-step for computing Stokes-tendency [T ~> s]
  logical, intent(in)                     :: dynamics_step !< True if this call is on a dynamics step

  ! Local Variables

end subroutine Update_Stokes_Drift
real module function one_minus_exp_x(x)
  real, intent(in) :: x !< The argument of the function ((1 - exp(-x))/x) [nondim]
end function one_minus_exp_x
real module function one_minus_exp(x)
  real, intent(in) :: x !< The argument of the function ((1 - exp(-x))/x) [nondim]
end function one_minus_exp
module subroutine Surface_Bands_by_data_override(Time, G, GV, US, CS)
  type(time_type),          intent(in) :: Time       !< Time to get Stokes drift bands
  type(wave_parameters_CS), pointer    :: CS         !< Wave structure
  type(ocean_grid_type), intent(inout) :: G          !< Grid structure
  type(verticalGrid_type),  intent(in) :: GV         !< Vertical grid structure
  type(unit_scale_type),    intent(in) :: US         !< A dimensional unit scaling type

  ! Local variables

end subroutine Surface_Bands_by_data_override
module subroutine get_Langmuir_Number( LA, G, GV, US, HBL, ustar, i, j, dz, Waves, &
                                U_H, V_H, Override_MA )
  type(ocean_grid_type),     intent(in)  :: G     !< Ocean grid structure
  type(verticalGrid_type),   intent(in)  :: GV    !< Ocean vertical grid structure
  real,                      intent(out) :: LA    !< Langmuir number [nondim]
  type(unit_scale_type),     intent(in)  :: US    !< A dimensional unit scaling type
  real,                      intent(in)  :: HBL   !< (Positive) thickness of boundary layer [Z ~> m]
  real,                      intent(in)  :: ustar !< Friction velocity [Z T-1 ~> m s-1]
  integer,                   intent(in)  :: i     !< Meridional index of h-point
  integer,                   intent(in)  :: j     !< Zonal index of h-point
  real, dimension(SZK_(GV)), intent(in)  :: dz    !< Grid layer thickness [Z ~> m]
  type(Wave_parameters_CS),  pointer     :: Waves !< Surface wave control structure.
  real, dimension(SZK_(GV)), &
                   optional, intent(in)  :: U_H   !< Zonal velocity at H point [L T-1 ~> m s-1] or [m s-1]
  real, dimension(SZK_(GV)), &
                   optional, intent(in)  :: V_H   !< Meridional velocity at H point [L T-1 ~> m s-1] or [m s-1]
  logical,         optional, intent(in)  :: Override_MA !< Override to use misalignment in LA
                                                  !! calculation. This can be used if diagnostic
                                                  !! LA outputs are desired that are different than
                                                  !! those used by the dynamical model.


!Local Variables

  ! Compute averaging depth for Stokes drift (negative)
end subroutine get_Langmuir_Number
module function get_wave_method(CS)
  character(:), allocatable :: get_wave_method
  type(wave_parameters_CS), pointer :: CS !< Control structure

end function get_wave_method
module subroutine get_StokesSL_LiFoxKemper(ustar, hbl, GV, US, CS, UStokes_SL, LA)
  real, intent(in)  :: ustar !< water-side surface friction velocity [Z T-1 ~> m s-1].
  real, intent(in)  :: hbl   !< boundary layer depth [Z ~> m].
  type(verticalGrid_type), intent(in) :: GV !< Ocean vertical grid structure
  type(unit_scale_type),   intent(in) :: US !< A dimensional unit scaling type
  type(wave_parameters_CS), pointer   :: CS  !< Wave parameter Control structure
  real, intent(out) :: UStokes_SL !< Surface layer averaged Stokes drift [L T-1 ~> m s-1]
  real, intent(out) :: LA    !< Langmuir number [nondim]
  ! Local variables
  ! parameters
                   ! Pierson-Moskowitz spectrum (Webb, 2011) [nondim]
                   ! boundary layer depth [nondim]

end subroutine get_StokesSL_LiFoxKemper
module subroutine Get_SL_Average_Prof( GV, AvgDepth, dz, Profile, Average )
  type(verticalGrid_type),  &
       intent(in)   :: GV       !< Ocean vertical grid structure
  real, intent(in)  :: AvgDepth !< Depth to average over (negative) [Z ~> m]
  real, dimension(SZK_(GV)), &
       intent(in)   :: dz       !< Grid thickness [Z ~> m]
  real, dimension(SZK_(GV)), &
       intent(in)   :: Profile  !< Profile of quantity to be averaged in arbitrary units [A]
                                !! (used here for Stokes drift)
  real, intent(out) :: Average  !< Output quantity averaged over depth AvgDepth [A]
                                !! (used here for Stokes drift)
  !Local variables

  ! Initializing sum
end subroutine Get_SL_Average_Prof
module subroutine Get_SL_Average_Band( GV, AvgDepth, NB, WaveNumbers, SurfStokes, Average )
  type(verticalGrid_type),  &
       intent(in)     :: GV          !< Ocean vertical grid
  real, intent(in)    :: AvgDepth    !< Depth to average over [Z ~> m].
  integer, intent(in) :: NB          !< Number of bands used
  real, dimension(NB), &
       intent(in)     :: WaveNumbers !< Wavenumber corresponding to each band [Z-1 ~> m-1]
  real, dimension(NB), &
       intent(in)     :: SurfStokes  !< Surface Stokes drift for each band [L T-1 ~> m s-1]
  real, intent(out)   :: Average     !< Output average Stokes drift over depth AvgDepth [L T-1 ~> m s-1]

  ! Local variables

  ! Loop over bands
end subroutine Get_SL_Average_Band
module subroutine DHH85_mid(GV, US, CS, zpt, UStokes)
  type(verticalGrid_type), intent(in)  :: GV  !< Ocean vertical grid
  type(unit_scale_type),   intent(in)  :: US  !< A dimensional unit scaling type
  type(wave_parameters_CS), pointer    :: CS  !< Wave parameter Control structure
  real, intent(in)  :: zpt   !< Depth to get Stokes drift [Z ~> m].
  real, intent(out) :: UStokes !< Stokes drift [L T-1 ~> m s-1]
  !

end subroutine DHH85_mid
module subroutine StokesMixing(G, GV, dt, h, dz, u, v, Waves )
  type(ocean_grid_type), &
       intent(in)    :: G     !< Ocean grid
  type(verticalGrid_type), &
       intent(in)    :: GV    !< Ocean vertical grid
  real, intent(in)   :: dt    !< Time step of MOM6 [T ~> s] for explicit solver
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
       intent(in)    :: h     !< Layer thicknesses [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
       intent(in)    :: dz    !< Vertical distance between interfaces around a layer [Z ~> m]
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
       intent(inout) :: u     !< Velocity i-component [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
       intent(inout) :: v     !< Velocity j-component [L T-1 ~> m s-1]
  type(Wave_parameters_CS), &
       pointer       :: Waves !< Surface wave related control structure.
  ! Local variables

! This is a template to think about down-Stokes mixing.
! This is not ready for use...

end subroutine StokesMixing
module subroutine CoriolisStokes(G, GV, dt, h, u, v, Waves)
  type(ocean_grid_type), &
       intent(in)    :: G     !< Ocean grid
  type(verticalGrid_type), &
       intent(in)   :: GV     !< Ocean vertical grid
  real, intent(in)  :: dt     !< Time step of MOM6 [T ~> s]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
       intent(in)    :: h     !< Layer thicknesses [H ~> m or kg m-2]
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
       intent(inout) :: u     !< Velocity i-component [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
       intent(inout) :: v     !< Velocity j-component [L T-1 ~> m s-1]
  type(Wave_parameters_CS), &
       pointer       :: Waves !< Surface wave related control structure.

  ! Local variables

end subroutine CoriolisStokes
module subroutine Stokes_PGF(G, GV, US, dz, u, v, PFu_Stokes, PFv_Stokes, CS )
  type(ocean_grid_type), &
       intent(in)    :: G     !< Ocean grid
  type(verticalGrid_type), &
       intent(in)    :: GV    !< Ocean vertical grid
  type(unit_scale_type), &
       intent(in)    :: US    !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)),&
       intent(in)    :: dz      !< Layer thicknesses in height units [Z ~> m]
  real, dimension(SZIB_(G),SZJ_(G),SZK_(G)), &
       intent(in) :: u          !< Lagrangian Velocity i-component [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJB_(G),SZK_(G)), &
       intent(in) :: v          !< Lagrangian Velocity j-component [L T-1 ~> m s-1]
  real, dimension(SZIB_(G),SZJ_(G),SZK_(G)), &
       intent(out) :: PFu_Stokes !< PGF Stokes-shear i-component [L T-2 ~> m s-2]
  real, dimension(SZI_(G),SZJB_(G),SZK_(G)), &
       intent(out) :: PFv_Stokes !< PGF Stokes-shear j-component [L T-2 ~> m s-2]
  type(Wave_parameters_CS), &
       pointer       :: CS !< Surface wave related control structure.

  ! Local variables
                                                              ! layer averaged [L2 T-2 ~> m2 s-2]
                                                              ! at interfaces [L2 T-2 ~> m2 s-2]
                                   ! (left/right of point) [L2 T-2 ~> m2 s-2]
                                         ! (left/right of point) [Z L2 T-2 ~> m3 s-2]
                                         ! (left/right of point) [L2 T-2 ~> m2 s-2]
                                   ! (left/right of point) [L2 T-2 ~> m2 s-2]
                                              ! contribution to Stokes pressure anomalies [nondim].


  !---------------------------------------------------------------
  ! Compute the Stokes contribution to the pressure gradient force
  !---------------------------------------------------------------
  ! Notes on the algorithm/code:
  ! This code requires computing velocities at bounding h points
  ! of the u/v points to get the pressure-gradient. In this
  ! implementation there are several redundant calculations as the
  ! left/right points are computed at each cell while integrating
  ! in the vertical, requiring about twice the calculations.  The
  ! velocities at the tracer points could be precomputed and
  ! stored, but this would require more memory and cycling through
  ! large 3d arrays while computing the pressures. This could be
  ! explored as a way to speed up this code.
  !---------------------------------------------------------------

end subroutine Stokes_PGF
module subroutine ust_2_u10_coare3p5(USTair, U10, GV, US, CS)
  real, intent(in)                    :: USTair !< Wind friction velocity [Z T-1 ~> m s-1]
  real, intent(out)                   :: U10    !< 10-m neutral wind speed [L T-1 ~> m s-1]
  type(verticalGrid_type), intent(in) :: GV     !< vertical grid type
  type(unit_scale_type),   intent(in) :: US     !< A dimensional unit scaling type
  type(wave_parameters_CS), pointer   :: CS     !< Wave parameter Control structure

  ! Local variables
                ! roughness length [nondim]

  ! Uses empirical formula for z0 to convert ustar_air to u10 based on the
  !  COARE 3.5 paper (Edson et al., 2013)
  ! alpha=m*U10+b
  ! Note in Edson et al. 2013, eq. 13 m is given as 0.017.  However,
  ! m=0.0017 reproduces the curve in their figure 6.

end subroutine ust_2_u10_coare3p5
module subroutine Waves_end(CS)
  type(wave_parameters_CS), pointer :: CS !< Control structure

end subroutine Waves_end
module subroutine waves_register_restarts(CS, HI, GV, US, param_file, restart_CSp)
  type(wave_parameters_CS), pointer       :: CS           !< Wave parameter Control structure
  type(hor_index_type),     intent(inout) :: HI           !< Grid structure
  type(verticalGrid_type),  intent(in)    :: GV           !< Vertical grid structure
  type(unit_scale_type),    intent(in)    :: US           !< A dimensional unit scaling type
  type(param_file_type),    intent(in)    :: param_file   !< Input parameter structure
  type(MOM_restart_CS),     pointer       :: restart_CSp  !< Restart structure, data intent(inout)
  ! Local variables

end subroutine waves_register_restarts
  end interface

end module MOM_wave_interface
