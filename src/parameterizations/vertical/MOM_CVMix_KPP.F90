! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> Provides the K-Profile Parameterization (KPP) of Large et al., 1994, via CVMix.
module MOM_CVMix_KPP

use MOM_coms,           only : max_across_PEs
use MOM_debugging,      only : hchksum, is_NaN
use MOM_diag_mediator,  only : time_type, diag_ctrl, safe_alloc_ptr, post_data
use MOM_diag_mediator,  only : query_averaging_enabled, register_diag_field
use MOM_error_handler,  only : MOM_error, MOM_mesg, FATAL, WARNING, is_root_PE
use MOM_EOS,            only : EOS_type, calculate_density
use MOM_file_parser,    only : get_param, log_param, log_version, param_file_type
use MOM_file_parser,    only : openParameterBlock, closeParameterBlock
use MOM_grid,           only : ocean_grid_type, isPointInCell
use MOM_interface_heights, only : thickness_to_dz
use MOM_restart,        only : MOM_restart_CS, register_restart_field
use MOM_unit_scaling,   only : unit_scale_type
use MOM_variables,      only : thermo_var_ptrs
use MOM_verticalGrid,   only : verticalGrid_type
use MOM_wave_interface, only : wave_parameters_CS, Get_Langmuir_Number, get_wave_method
use MOM_domains,        only : pass_var
use MOM_cpu_clock,      only : cpu_clock_id, cpu_clock_begin, cpu_clock_end
use MOM_cpu_clock,      only : CLOCK_MODULE, CLOCK_ROUTINE
use MOM_tracer_types,   only : tracer_type

use CVMix_kpp, only : CVMix_init_kpp, CVMix_put_kpp, CVMix_get_kpp_real
use CVMix_kpp, only : CVMix_coeffs_kpp
use CVMix_kpp, only : CVMix_kpp_compute_OBL_depth
use CVMix_kpp, only : CVMix_kpp_compute_turbulent_scales
use CVMix_kpp, only : CVMix_kpp_compute_bulk_Richardson
use CVMix_kpp, only : CVMix_kpp_compute_unresolved_shear
use CVMix_kpp, only : CVMix_kpp_params_type
use CVMix_kpp, only : CVMix_kpp_compute_kOBL_depth
use CVMix_kpp, only : CVMix_kpp_compute_StokesXi

implicit none ; private

#include "MOM_memory.h"

public :: register_KPP_restarts
public :: KPP_init
public :: KPP_compute_BLD
public :: KPP_calculate
public :: KPP_end
public :: KPP_NonLocalTransport_temp
public :: KPP_NonLocalTransport_saln
public :: KPP_NonLocalTransport
public :: KPP_get_BLD

! Enumerated constants
integer, private, parameter :: NLT_SHAPE_CVMix     = 0 !< Use the CVMix profile
integer, private, parameter :: NLT_SHAPE_LINEAR    = 1 !< Linear, \f$ G(\sigma) = 1-\sigma \f$
integer, private, parameter :: NLT_SHAPE_PARABOLIC = 2 !< Parabolic, \f$ G(\sigma) = (1-\sigma)^2 \f$
integer, private, parameter :: NLT_SHAPE_CUBIC     = 3 !< Cubic, \f$ G(\sigma) = 1 + (2\sigma-3) \sigma^2\f$
integer, private, parameter :: NLT_SHAPE_CUBIC_LMD = 4 !< Original shape,
                                                       !!    \f$ G(\sigma) = \frac{27}{4} \sigma (1-\sigma)^2 \f$

integer, private, parameter :: SW_METHOD_ALL_SW = 0 !< Use all shortwave radiation
integer, private, parameter :: SW_METHOD_MXL_SW = 1 !< Use shortwave radiation absorbed in mixing layer
integer, private, parameter :: SW_METHOD_LV1_SW = 2 !< Use shortwave radiation absorbed in layer 1
integer, private, parameter :: LT_K_CONSTANT = 1,        & !< Constant enhance K through column
                               LT_K_SCALED = 2,          & !< Enhance K scales with G(sigma)
                               LT_K_MODE_CONSTANT = 1,   & !< Prescribed enhancement for K
                               LT_K_MODE_VR12 = 2,       & !< Enhancement for K based on
                                                           !! Van Roekel et al., 2012
                               LT_K_MODE_RW16 = 3,       & !< Enhancement for K based on
                                                           !! Reichl et al., 2016
                               LT_VT2_MODE_CONSTANT = 1, & !< Prescribed enhancement for Vt2
                               LT_VT2_MODE_VR12 = 2,     & !< Enhancement for Vt2 based on
                                                           !! Van Roekel et al., 2012
                               LT_VT2_MODE_RW16 = 3,     & !< Enhancement for Vt2 based on
                                                           !! Reichl et al., 2016
                               LT_VT2_MODE_LF17 = 4        !< Enhancement for Vt2 based on
                                                           !! Li and Fox-Kemper, 2017

!> Control structure for containing KPP parameters/data
type, public :: KPP_CS ; private

  ! Parameters
  real    :: Ri_crit                   !< Critical bulk Richardson number (defines OBL depth) [nondim]
  real    :: vonKarman                 !< von Karman constant (dimensionless) [nondim]
  real    :: cs                        !< Parameter for computing velocity scale function (dimensionless) [nondim]
  real    :: cs2                       !< Parameter for multiplying by non-local term [nondim]
                                       !   This is active for NLT_SHAPE_CUBIC_LMD only
  logical :: enhance_diffusion         !< If True, add enhanced diffusivity at base of boundary layer.
  character(len=32) :: interpType      !< Type of interpolation to compute bulk Richardson number
  character(len=32) :: interpType2     !< Type of interpolation to compute diff and visc at OBL_depth
  logical :: StokesMOST                !< If True, use Stokes similarity package
  logical :: computeEkman              !< If True, compute Ekman depth limit for OBLdepth
  logical :: computeMoninObukhov       !< If True, compute Monin-Obukhov limit for OBLdepth
  logical :: passiveMode               !< If True, makes KPP passive meaning it does NOT alter the diffusivity
  real    :: deepOBLoffset             !< If non-zero, is a distance from the bottom that the OBL can not
                                       !! penetrate through [Z ~> m]
  real    :: minOBLdepth               !< If non-zero, is a minimum depth for the OBL [Z ~> m]
  real    :: surf_layer_ext            !< Fraction of OBL depth considered in the surface layer [nondim]
  real    :: minVtsqr                  !< Min for the squared unresolved velocity used in Rib CVMix
                                       !! calculation [L2 T-2 ~> m2 s-2]
  logical :: fixedOBLdepth             !< If True, will fix the OBL depth at fixedOBLdepth_value
  real    :: fixedOBLdepth_value       !< value for the fixed OBL depth when fixedOBLdepth==True [Z ~> m]
  logical :: debug                     !< If True, calculate checksums and write debugging information
  character(len=30) :: MatchTechnique  !< Method used in CVMix for setting diffusivity and NLT profile functions
  integer :: NLT_shape                 !< MOM6 over-ride of CVMix NLT shape function
  logical :: applyNonLocalTrans        !< If True, apply non-local transport to all tracers
  integer :: n_smooth                  !< Number of times smoothing operator is applied on OBLdepth.
  logical :: deepen_only               !< If true, apply OBLdepth smoothing at a cell only if the OBLdepth gets deeper.
  logical :: KPPzeroDiffusivity        !< If True, will set diffusivity and viscosity from KPP to zero
                                       !! for testing purposes.
  logical :: KPPisAdditive             !< If True, will add KPP diffusivity to initial diffusivity.
                                       !! If False, will replace initial diffusivity wherever KPP diffusivity
                                       !! is non-zero.
  real    :: min_thickness             !< A minimum thickness used to avoid division by small numbers
                                       !! in the vicinity of vanished layers [Z ~> m]
  integer :: SW_METHOD                 !< Sets method for using shortwave radiation in surface buoyancy flux
  logical :: LT_K_Enhancement          !< Flags if enhancing mixing coefficients due to LT
  integer :: LT_K_Shape                !< Integer for constant or shape function enhancement
  integer :: LT_K_Method               !< Integer for mixing coefficients LT method
  real    :: KPP_CVt2                  !< Parameter for Stokes MOST convection entrainment [nondim]
  real    :: KPP_K_ENH_FAC             !< Factor to multiply by K if Method is CONSTANT [nondim]
  logical :: LT_Vt2_Enhancement        !< Flags if enhancing Vt2 due to LT
  integer :: LT_VT2_METHOD             !< Integer for Vt2 LT method
  real    :: KPP_VT2_ENH_FAC           !< Factor to multiply by VT2 if Method is CONSTANT [nondim]
  real    :: MLD_guess_min             !< The minimum estimate of the mixed layer depth used to
                                       !! calculate the Langmuir number for Langmuir turbulence
                                       !! enhancement with KPP [Z ~> m]
  logical :: STOKES_MIXING             !< Flag if model is mixing down Stokes gradient
                                       !! This is relevant for which current to use in RiB
  integer :: answer_date               !< The vintage of the order of arithmetic in the CVMix KPP
                                       !! calculations.  Values below 20240501 recover the answers
                                       !! from early in 2024, while higher values use expressions
                                       !! that have been refactored for rotational symmetry.

  !> CVMix parameters
  type(CVMix_kpp_params_type), pointer :: KPP_params => NULL()

  type(diag_ctrl), pointer :: diag => NULL() !< Pointer to diagnostics control structure
  !>@{ Diagnostic handles
  integer :: id_OBLdepth = -1, id_BulkRi   = -1
  integer :: id_N        = -1, id_N2       = -1
  integer :: id_Ws       = -1, id_Vt2      = -1
  integer :: id_BulkUz2  = -1, id_BulkDrho = -1
  integer :: id_uStar    = -1, id_buoyFlux = -1
  integer :: id_sigma    = -1, id_Kv_KPP   = -1
  integer :: id_Kt_KPP   = -1, id_Ks_KPP   = -1
  integer :: id_Tsurf    = -1, id_Ssurf    = -1
  integer :: id_Usurf    = -1, id_Vsurf    = -1
  integer :: id_Kd_in    = -1
  integer :: id_NLTt     = -1
  integer :: id_NLTs     = -1
  integer :: id_EnhK     = -1, id_EnhVt2   = -1
  integer :: id_EnhW     = -1
  integer :: id_La_SL    = -1
  integer :: id_OBLdepth_original = -1
  integer :: id_StokesXI = -1
  integer :: id_Lam2     = -1
  !>@}

  ! Diagnostics arrays
  real, pointer,     dimension(:,:)   :: OBLdepth  !< Depth (positive) of ocean boundary layer (OBL) [Z ~> m]
  real, allocatable, dimension(:,:)   :: OBLdepth_original  !< Depth (positive) of OBL without smoothing [Z ~> m]
  real, allocatable, dimension(:,:)   :: StokesParXI !< Stokes similarity parameter [nondim]
  real, allocatable, dimension(:,:)   :: Lam2      !< La^(-2) = Ustk0/u* [nondim]
  real, allocatable, dimension(:,:)   :: kOBL      !< Level (+fraction) of OBL extent [nondim]
  real, allocatable, dimension(:,:)   :: OBLdepthprev !< previous Depth (positive) of OBL [Z ~> m]
  real, allocatable, dimension(:,:)   :: La_SL     !< Langmuir number used in KPP [nondim]
  real, allocatable, dimension(:,:,:) :: dRho      !< Bulk difference in density [R ~> kg m-3]
  real, allocatable, dimension(:,:,:) :: Uz2       !< Square of bulk difference in resolved velocity [L2 T-2 ~> m2 s-2]
  real, allocatable, dimension(:,:,:) :: BulkRi    !< Bulk Richardson number for each layer [nondim]
  real, allocatable, dimension(:,:,:) :: sigma     !< Sigma coordinate (dimensionless) [nondim]
  real, allocatable, dimension(:,:,:) :: Ws        !< Turbulent velocity scale for scalars [Z T-1 ~> m s-1]
  real, allocatable, dimension(:,:,:) :: N         !< Brunt-Vaisala frequency [T-1 ~> s-1]
  real, allocatable, dimension(:,:,:) :: N2        !< Squared Brunt-Vaisala frequency [T-2 ~> s-2]
  real, allocatable, dimension(:,:,:) :: Vt2       !< Unresolved squared turbulence velocity for
                                                   !! bulk Ri [Z2 T-2 ~> m2 s-2]
  real, allocatable, dimension(:,:,:) :: Kt_KPP    !< Temp diffusivity from KPP [Z2 T-1 ~> m2 s-1]
  real, allocatable, dimension(:,:,:) :: Ks_KPP    !< Scalar diffusivity from KPP [Z2 T-1 ~> m2 s-1]
  real, allocatable, dimension(:,:,:) :: Kv_KPP    !< Viscosity due to KPP [Z2 T-1 ~> m2 s-1]
  real, allocatable, dimension(:,:)   :: Tsurf     !< Temperature of surface layer [C ~> degC]
  real, allocatable, dimension(:,:)   :: Ssurf     !< Salinity of surface layer [S ~> ppt]
  real, allocatable, dimension(:,:)   :: Usurf     !< i-velocity of surface layer [L T-1 ~> m s-1]
  real, allocatable, dimension(:,:)   :: Vsurf     !< j-velocity of surface layer [L T-1 ~> m s-1]
  real, allocatable, dimension(:,:,:) :: EnhK      !< Enhancement for mixing coefficient [nondim]
  real, allocatable, dimension(:,:,:) :: EnhVt2    !< Enhancement for Vt2 [nondim]

end type KPP_CS

!>@{ CPU time clocks
integer :: id_clock_KPP_calc, id_clock_KPP_compute_BLD, id_clock_KPP_smoothing
!>@}

#define __DO_SAFETY_CHECKS__


  interface
module subroutine register_KPP_restarts(G, param_file, restart_CSp, CS)
  type(ocean_grid_type), intent(in)    :: G           !< The ocean's grid structure
  type(param_file_type), intent(in)    :: param_file  !< A structure to parse for run-time parameters
  type(MOM_restart_CS),  pointer       :: restart_CSp  !< MOM restart control structure
  type(KPP_CS),         pointer        :: CS           !< module control structure


end subroutine register_KPP_restarts
logical module function KPP_init(paramFile, G, GV, US, diag, Time, CS, passive)

  ! Arguments
  type(param_file_type),   intent(in)    :: paramFile !< File parser
  type(ocean_grid_type),   intent(in)    :: G         !< Ocean grid
  type(verticalGrid_type), intent(in)    :: GV        !< Vertical grid structure
  type(unit_scale_type),   intent(in)    :: US        !< A dimensional unit scaling type
  type(diag_ctrl), target, intent(in)    :: diag      !< Diagnostics
  type(time_type),         intent(in)    :: Time      !< Model time
  type(KPP_CS),            pointer       :: CS        !< Control structure
  logical,       optional, intent(out)   :: passive   !< Copy of %passiveMode

  ! Local variables
                                       !! passed to CVMix, e.g., LWF16
                                       !! False => compute G'(1) as in LMD94
  ! Read parameters
end function KPP_init
module subroutine KPP_calculate(CS, G, GV, US, h, tv, uStar, buoyFlux, Kt, Ks, Kv, &
                         nonLocalTransHeat, nonLocalTransScalar, Waves, lamult)

  ! Arguments
  type(KPP_CS),                                pointer       :: CS    !< Control structure
  type(ocean_grid_type),                       intent(in)    :: G     !< Ocean grid
  type(verticalGrid_type),                     intent(in)    :: GV    !< Ocean vertical grid
  type(unit_scale_type),                       intent(in)    :: US    !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),   intent(in)    :: h     !< Layer thicknesses [H ~> m or kg m-2]
  type(thermo_var_ptrs),                       intent(in)    :: tv    !< Thermodynamics structure.
  real, dimension(SZI_(G),SZJ_(G)),            intent(in)    :: uStar !< Surface friction velocity [Z T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)+1), intent(in)    :: buoyFlux !< Surface buoyancy flux [L2 T-3 ~> m2 s-3]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)+1), intent(inout) :: Kt  !< (in)  Vertical diffusivity of heat w/o KPP
                                                                    !! (out) Vertical diffusivity including KPP
                                                                    !!       [H Z T-1 ~> m2 s-1 or kg m-1 s-1]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)+1), intent(inout) :: Ks  !< (in)  Vertical diffusivity of salt w/o KPP
                                                                    !! (out) Vertical diffusivity including KPP
                                                                    !!       [H Z T-1 ~> m2 s-1 or kg m-1 s-1]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)+1), intent(inout) :: Kv  !< (in)  Vertical viscosity w/o KPP
                                                                    !! (out) Vertical viscosity including KPP
                                                                    !!       [H Z T-1 ~> m2 s-1 or Pa s]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)+1), intent(inout) :: nonLocalTransHeat   !< Temp non-local transport [nondim]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)+1), intent(inout) :: nonLocalTransScalar !< scalar non-local trans. [nondim]
  type(wave_parameters_CS),                    pointer       :: Waves   !< Wave CS for Langmuir turbulence
  real, dimension(SZI_(G),SZJ_(G)),  optional, intent(in)    :: lamult  !< Langmuir enhancement multiplier [nondim]

  ! Local variables


  ! For Langmuir Calculations

end subroutine KPP_calculate
module subroutine KPP_compute_BLD(CS, G, GV, US, h, Temp, Salt, u, v, tv, uStar, buoyFlux, Waves, lamult)

  ! Arguments
  type(KPP_CS),                               pointer       :: CS    !< Control structure
  type(ocean_grid_type),                      intent(inout) :: G     !< Ocean grid
  type(verticalGrid_type),                    intent(in)    :: GV    !< Ocean vertical grid
  type(unit_scale_type),                      intent(in)    :: US    !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),  intent(in)    :: h     !< Layer thicknesses [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),  intent(in)    :: Temp  !< potential/cons temp [C ~> degC]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),  intent(in)    :: Salt  !< Salinity [S ~> ppt]
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), intent(in)    :: u     !< Velocity i-component [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), intent(in)    :: v     !< Velocity j-component [L T-1 ~> m s-1]
  type(thermo_var_ptrs),                      intent(in)    :: tv    !< Thermodynamics structure.
  real, dimension(SZI_(G),SZJ_(G)),           intent(in)    :: uStar !< Surface friction velocity [Z T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)+1), intent(in)   :: buoyFlux !< Surface buoyancy flux [L2 T-3 ~> m2 s-3]
  type(wave_parameters_CS),                   pointer       :: Waves !< Wave CS for Langmuir turbulence
  real, dimension(SZI_(G),SZJ_(G)), optional, intent(in)    :: lamult !< Langmuir enhancement factor [nondim]

  ! Local variables
  ! Variables for passing to CVMix routines, often in MKS units

  ! Variables for EOS calculations

                        ! rescaling [H T-2 R-1 ~> m4 kg-1 s-2 or m s-2]

  ! For Langmuir Calculations


                                                                 ! [L T-1 ~> m s-1]

end subroutine KPP_compute_BLD
module subroutine KPP_smooth_BLD(CS, G, GV, US, dz)
  ! Arguments
  type(KPP_CS),                           pointer       :: CS   !< Control structure
  type(ocean_grid_type),                  intent(inout) :: G    !< Ocean grid
  type(verticalGrid_type),                intent(in)    :: GV   !< Ocean vertical grid
  type(unit_scale_type),                  intent(in)    :: US   !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in) :: dz   !< Layer thicknesses [Z ~> m]

  ! local variables
                                                        ! for the minimum layer thickness [Z ~> m]
                                                        ! (negative in the ocean)
                                                        ! (negative in the ocean)

end subroutine KPP_smooth_BLD
module subroutine KPP_get_BLD(CS, BLD, G, US, m_to_BLD_units)
  type(KPP_CS),                     pointer     :: CS  !< Control structure for
                                                       !! this module
  type(ocean_grid_type),            intent(in)  :: G   !< Grid structure
  type(unit_scale_type),            intent(in)  :: US  !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G)), intent(inout) :: BLD !< Boundary layer depth [Z ~> m] or other units
  real,                   optional, intent(in)  :: m_to_BLD_units !< A conversion factor from meters
                                                       !! to the desired units for BLD [various]
  ! Local variables

end subroutine KPP_get_BLD
module subroutine KPP_NonLocalTransport(CS, G, GV, h, nonLocalTrans, surfFlux, &
                                 dt, diag, tr_ptr, scalar, flux_scale)
  type(KPP_CS),                               intent(in)    :: CS            !< Control structure
  type(ocean_grid_type),                      intent(in)    :: G             !< Ocean grid
  type(verticalGrid_type),                    intent(in)    :: GV            !< Ocean vertical grid
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),  intent(in)    :: h             !< Layer/level thickness [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)+1), intent(in)   :: nonLocalTrans !< Non-local transport [nondim]
  real, dimension(SZI_(G),SZJ_(G)),           intent(in)    :: surfFlux      !< Surface flux of scalar
                                                                        !! [conc H T-1 ~> conc m s-1 or conc kg m-2 s-1]
  real,                                       intent(in)    :: dt            !< Time-step [T ~> s]
  type(diag_ctrl), target,                    intent(in)    :: diag          !< Diagnostics
  type(tracer_type), pointer,                 intent(in)    :: tr_ptr        !< tracer_type has diagnostic ids on it
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),  intent(inout) :: scalar        !< Scalar (scalar units [conc])
  real, optional,                             intent(in)    :: flux_scale    !< Scale factor to get surfFlux
                                                                             !! into proper units [various]

                                                   ! in [conc H T-1 ~> conc m s-1 or conc kg m-2 s-1] or other units

  ! term used to scale
end subroutine KPP_NonLocalTransport
module subroutine KPP_NonLocalTransport_temp(CS, G, GV, h, nonLocalTrans, surfFlux, dt, tr_ptr, scalar, C_p)
  type(KPP_CS),                               intent(in)    :: CS     !< Control structure
  type(ocean_grid_type),                      intent(in)    :: G      !< Ocean grid
  type(verticalGrid_type),                    intent(in)    :: GV     !< Ocean vertical grid
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),  intent(in)    :: h      !< Layer/level thickness [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)+1), intent(in)   :: nonLocalTrans !< Non-local transport [nondim]
  real, dimension(SZI_(G),SZJ_(G)),           intent(in)    :: surfFlux  !< Surface flux of temperature
                                                                      !! [C H T-1 ~> degC m s-1 or degC kg m-2 s-1]
  real,                                       intent(in)    :: dt     !< Time-step [T ~> s]
  type(tracer_type), pointer,                 intent(in)    :: tr_ptr !< tracer_type has diagnostic ids on it
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),  intent(inout) :: scalar !< temperature [C ~> degC]
  real,                                       intent(in)    :: C_p    !< Seawater specific heat capacity
                                                                      !! [Q C-1 ~> J kg-1 degC-1]

end subroutine KPP_NonLocalTransport_temp
module subroutine KPP_NonLocalTransport_saln(CS, G, GV, h, nonLocalTrans, surfFlux, dt, tr_ptr, scalar)
  type(KPP_CS),                               intent(in)    :: CS            !< Control structure
  type(ocean_grid_type),                      intent(in)    :: G             !< Ocean grid
  type(verticalGrid_type),                    intent(in)    :: GV            !< Ocean vertical grid
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),  intent(in)    :: h             !< Layer/level thickness [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)+1), intent(in)   :: nonLocalTrans !< Non-local transport [nondim]
  real, dimension(SZI_(G),SZJ_(G)),           intent(in)    :: surfFlux      !< Surface flux of salt
                                                                             !! [S H T-1 ~> ppt m s-1 or ppt kg m-2 s-1]
  real,                                       intent(in)    :: dt            !< Time-step [T ~> s]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),  intent(inout) :: scalar        !< Salinity [S ~> ppt]
  type(tracer_type), pointer,                 intent(in)    :: tr_ptr        !< tracer_type has diagnostic ids on it

end subroutine KPP_NonLocalTransport_saln
module subroutine Compute_StokesDrift(i ,j, ztop, zbot, uS_i, vS_i, uS_k, vS_k, uSbar, vSbar, waves)

  type(wave_parameters_CS), pointer  :: waves  !< Wave CS for Langmuir turbulence
  real,                intent(in)    :: ztop   !< cell top
  real,                intent(in)    :: zbot   !< cell bottom
  real,                intent(inout) :: uS_i   !< Stokes u velocity at zbot interface
  real,                intent(inout) :: vS_i   !< Stokes v velocity at zbot interface
  real,                intent(inout) :: uS_k   !< Stokes u velocity at zk center
  real,                intent(inout) :: vS_k   !< Stokes v at zk =0.5(ztop+zbot)
  real,                intent(inout) :: uSbar  !< mean Stokes u (ztop to zbot)
  real,                intent(inout) :: vSbar  !< mean Stokes v (ztop to zbot)
  integer,             intent(in)    :: i      !< Meridional index of H-point
  integer,             intent(in)    :: j      !< Zonal index of H-point

  ! local variables

end subroutine Compute_StokesDrift
module subroutine KPP_end(CS)
  type(KPP_CS), pointer :: CS !< Control structure

end subroutine KPP_end
  end interface

end module MOM_CVMix_KPP
