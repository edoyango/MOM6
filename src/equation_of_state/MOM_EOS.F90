! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> Provides subroutines for quantities specific to the equation of state
module MOM_EOS

use MOM_EOS_base_type, only : EOS_base
use MOM_EOS_linear, only : linear_EOS, avg_spec_vol_linear
use MOM_EOS_linear, only : int_density_dz_linear, int_spec_vol_dp_linear
use MOM_EOS_Wright, only : buggy_Wright_EOS, avg_spec_vol_buggy_Wright
use MOM_EOS_Wright, only : int_density_dz_wright, int_spec_vol_dp_wright
use MOM_EOS_Wright_full, only : Wright_full_EOS, avg_spec_vol_Wright_full
use MOM_EOS_Wright_full, only : int_density_dz_wright_full, int_spec_vol_dp_wright_full
use MOM_EOS_Wright_red,  only : Wright_red_EOS, avg_spec_vol_Wright_red
use MOM_EOS_Wright_red,  only : int_density_dz_wright_red, int_spec_vol_dp_wright_red
use MOM_EOS_Jackett06, only : Jackett06_EOS
use MOM_EOS_UNESCO, only : UNESCO_EOS
use MOM_EOS_Roquet_rho, only : Roquet_rho_EOS
use MOM_EOS_Roquet_SpV, only : Roquet_SpV_EOS
use MOM_EOS_TEOS10, only : TEOS10_EOS
use MOM_EOS_TEOS10, only : gsw_sp_from_sr, gsw_pt_from_ct, gsw_sr_from_sp, gsw_ct_from_pt
use MOM_temperature_convert, only : poTemp_to_consTemp, consTemp_to_poTemp
use MOM_TFreeze,    only : calculate_TFreeze_linear, calculate_TFreeze_Millero
use MOM_TFreeze,    only : calculate_TFreeze_teos10, calculate_TFreeze_TEOS_poly
use MOM_error_handler, only : MOM_error, FATAL, WARNING, MOM_mesg
use MOM_file_parser, only : get_param, log_version, param_file_type
use MOM_hor_index,   only : hor_index_type
use MOM_io,          only : stdout, stderr
use MOM_string_functions, only : uppercase
use MOM_unit_scaling, only : unit_scale_type

implicit none ; private

public EOS_domain
public EOS_init
public EOS_manual_init
public EOS_quadrature
! public EOS_use_linear
public EOS_fit_range
public EOS_unit_tests
public analytic_int_density_dz
public analytic_int_specific_vol_dp
public average_specific_vol
public calculate_compress
public calculate_density_elem
public calculate_density
public calculate_density_derivs
public calculate_density_second_derivs
public calculate_spec_vol
public calculate_specific_vol_derivs
public calculate_TFreeze
public convert_temp_salt_for_TEOS10
public cons_temp_to_pot_temp
public pot_temp_to_cons_temp
public abs_saln_to_prac_saln
public prac_saln_to_abs_saln
public gsw_sp_from_sr
public gsw_sr_from_sp
public gsw_pt_from_ct
public query_compressible
public get_EOS_name

! A note on unit descriptions in comments: MOM6 uses units that can be rescaled for dimensional
! consistency testing. These are noted in comments with units like Z, H, L, and T, along with
! their mks counterparts with notation like "a velocity [Z T-1 ~> m s-1]".  If the units
! vary with the Boussinesq approximation, the Boussinesq variant is given first.

!> Calculates density of sea water from T, S and P
interface calculate_density
  module procedure calculate_density_scalar
  module procedure calculate_density_1d
  module procedure calculate_stanley_density_scalar
  module procedure calculate_stanley_density_1d
end interface calculate_density

!> Calculates specific volume of sea water from T, S and P
interface calculate_spec_vol
  module procedure calc_spec_vol_scalar
  module procedure calc_spec_vol_1d
end interface calculate_spec_vol

!> Calculate the derivatives of density with temperature and salinity from T, S, and P
interface calculate_density_derivs
  module procedure calculate_density_derivs_scalar, calculate_density_derivs_array
  module procedure calculate_density_derivs_1d
end interface calculate_density_derivs

!> Calculate the derivatives of specific volume with temperature and salinity from T, S, and P
interface calculate_specific_vol_derivs
  module procedure calc_spec_vol_derivs_1d
end interface calculate_specific_vol_derivs

!> Calculates the second derivatives of density with various combinations of temperature,
!! salinity, and pressure from T, S and P
interface calculate_density_second_derivs
  module procedure calculate_density_second_derivs_scalar, calculate_density_second_derivs_1d
end interface calculate_density_second_derivs

!> Calculates the freezing point of sea water from T, S and P
interface calculate_TFreeze
  module procedure calculate_TFreeze_scalar, calculate_TFreeze_1d, calculate_TFreeze_array
end interface calculate_TFreeze

!> Calculates the compressibility of water from T, S, and P
interface calculate_compress
  module procedure calculate_compress_scalar, calculate_compress_1d
end interface calculate_compress

!> A control structure for the equation of state
type, public :: EOS_type ; private
  integer :: form_of_EOS = 0 !< The equation of state to use.
  integer :: form_of_TFreeze = 0 !< The expression for the potential temperature
                             !! of the freezing point.
  logical :: EOS_quadrature  !< If true, always use the generic (quadrature)
                             !! code for the integrals of density.
  logical :: Compressible = .true. !< If true, in situ density is a function of pressure.
! The following parameters are used with the linear equation of state only.
  real :: Rho_T0_S0 !< The density at T=0, S=0 [kg m-3]
  real :: dRho_dT   !< The partial derivative of density with temperature [kg m-3 degC-1]
  real :: dRho_dS   !< The partial derivative of density with salinity [kg m-3 ppt-1]
  real :: dRho_dp   !< The partial derivative of density with pressure [s2 m-2]
! The following parameters are use with the linear expression for the freezing
! point only.
  real :: TFr_S0_P0 !< The freezing potential temperature at S=0, P=0 [degC]
  real :: dTFr_dS   !< The derivative of freezing point with salinity [degC ppt-1]
  real :: dTFr_dp   !< The derivative of freezing point with pressure [degC Pa-1]
! The following are logicals pertaining to definitions of the thermodynamic state variables
  logical :: use_conT_absS =.false. !< True if the model internal temperature is the conservative temperature and
                           !! the salinity is absolute salinity.  These could be separated into two flags,
                           !! but right now it is controlled by one input parameter and there is no known
                           !! need to have one True and one False.
  logical :: TFreeze_S_is_pracS =.true. !< True if the freezing point expression is formulated from practical salinity
  logical :: TFreeze_T_is_potT = .true. !< True if the freezing point expression yields a potential temperature

  logical :: use_Wright_2nd_deriv_bug = .false.  !< If true, use a separate subroutine that
                           !! retains a buggy version of the calculations of the second
                           !! derivative of density with temperature and with temperature and
                           !! pressure.  This bug is corrected in the default version.

! Unit conversion factors (normally used for dimensional testing but could also allow for
! change of units of arguments to functions)
  real :: m_to_Z = 1.      !< A constant that translates distances in meters to the units of depth [Z m-1 ~> 1]
  real :: kg_m3_to_R = 1.  !< A constant that translates kilograms per meter cubed to the
                           !! units of density [R m3 kg-1 ~> 1]
  real :: R_to_kg_m3 = 1.  !< A constant that translates the units of density to
                           !! kilograms per meter cubed [kg m-3 R-1 ~> 1]
  real :: RL2_T2_to_Pa = 1.!< Convert pressures from R L2 T-2 to Pa [Pa T2 R-1 L-2 ~> 1]
  real :: L_T_to_m_s = 1.  !< Convert lateral velocities from L T-1 to m s-1 [m T s-1 L-1 ~> 1]
  real :: degC_to_C = 1.   !< A constant that translates degrees Celsius to the units of temperature [C degC-1 ~> 1]
  real :: C_to_degC = 1.   !< A constant that translates the units of temperature to degrees Celsius [degC C-1 ~> 1]
  real :: ppt_to_S = 1.    !< A constant that translates parts per thousand to the units of salinity [S ppt-1 ~> 1]
  real :: S_to_ppt = 1.    !< A constant that translates the units of salinity to parts per thousand [ppt S-1 ~> 1]

  !> The instance of the actual equation of state
  class(EOS_base), allocatable :: type

end type EOS_type

! The named integers that might be stored in eqn_of_state_type%form_of_EOS.
integer, parameter, public :: EOS_LINEAR = 1 !< A named integer specifying an equation of state
integer, parameter, public :: EOS_UNESCO = 2 !< A named integer specifying an equation of state
integer, parameter, public :: EOS_WRIGHT = 3 !< A named integer specifying an equation of state
integer, parameter, public :: EOS_WRIGHT_FULL = 4 !< A named integer specifying an equation of state
integer, parameter, public :: EOS_WRIGHT_REDUCED = 5 !< A named integer specifying an equation of state
integer, parameter, public :: EOS_TEOS10 = 6 !< A named integer specifying an equation of state
integer, parameter, public :: EOS_ROQUET_RHO = 7 !< A named integer specifying an equation of state
integer, parameter, public :: EOS_ROQUET_SPV = 8 !< A named integer specifying an equation of state
integer, parameter, public :: EOS_JACKETT06 = 9 !< A named integer specifying an equation of state
!> A list of all the available EOS
integer, dimension(9), public :: list_of_EOS = (/ EOS_LINEAR, EOS_UNESCO, &
            EOS_WRIGHT, EOS_WRIGHT_FULL, EOS_WRIGHT_REDUCED, &
            EOS_TEOS10, EOS_ROQUET_RHO, EOS_ROQUET_SPV, EOS_JACKETT06 /)

character*(12), parameter :: EOS_LINEAR_STRING = "LINEAR" !< A string for specifying the equation of state
character*(12), parameter :: EOS_UNESCO_STRING = "UNESCO" !< A string for specifying the equation of state
character*(12), parameter :: EOS_JACKETT_STRING = "JACKETT_MCD" !< A string for specifying the equation of state
character*(12), parameter :: EOS_WRIGHT_STRING = "WRIGHT" !< A string for specifying the equation of state
character*(16), parameter :: EOS_WRIGHT_RED_STRING = "WRIGHT_REDUCED" !< A string for specifying the equation of state
character*(12), parameter :: EOS_WRIGHT_FULL_STRING = "WRIGHT_FULL" !< A string for specifying the equation of state
character*(12), parameter :: EOS_TEOS10_STRING = "TEOS10" !< A string for specifying the equation of state
character*(12), parameter :: EOS_NEMO_STRING   = "NEMO"   !< A string for specifying the equation of state
character*(12), parameter :: EOS_ROQUET_RHO_STRING = "ROQUET_RHO"   !< A string for specifying the equation of state
character*(12), parameter :: EOS_ROQUET_SPV_STRING = "ROQUET_SPV"   !< A string for specifying the equation of state
character*(12), parameter :: EOS_JACKETT06_STRING = "JACKETT_06" !< A string for specifying the equation of state
character*(12), parameter :: EOS_DEFAULT = EOS_WRIGHT_FULL_STRING !< The default equation of state

integer, parameter :: TFREEZE_LINEAR = 1  !< A named integer specifying a freezing point expression
integer, parameter :: TFREEZE_MILLERO = 2 !< A named integer specifying a freezing point expression
integer, parameter :: TFREEZE_TEOS10 = 3  !< A named integer specifying a freezing point expression
integer, parameter :: TFREEZE_TEOSPOLY = 4 !< A named integer specifying a freezing point expression
character*(10), parameter :: TFREEZE_LINEAR_STRING = "LINEAR" !< A string for specifying the freezing point expression
character*(10), parameter :: TFREEZE_MILLERO_STRING = "MILLERO_78" !< A string for specifying the
                                                              !! freezing point expression
character*(10), parameter :: TFREEZE_TEOSPOLY_STRING = "TEOS_POLY" !< A string for specifying the
                                                              !! freezing point expression
character*(10), parameter :: TFREEZE_TEOS10_STRING = "TEOS10" !< A string for specifying the freezing point expression


  interface
real elemental module function calculate_density_elem(EOS, T, S, pressure, rho_ref, scale)
  type(EOS_type), intent(in)  :: EOS      !< Equation of state structure
  real,           intent(in)  :: T        !< Potential temperature referenced to the surface [C ~> degC]
  real,           intent(in)  :: S        !< Salinity [S ~> ppt]
  real,           intent(in)  :: pressure !< Pressure [R L2 T-2 ~> Pa]
  real, optional, intent(in)  :: rho_ref  !< A reference density [R ~> kg m-3]
  real, optional, intent(in)  :: scale    !< A multiplicative factor by which to scale output density in
                                          !! combination with scaling stored in EOS [various]

end function calculate_density_elem
module subroutine calculate_density_scalar(T, S, pressure, rho, EOS, rho_ref, scale)
  real,           intent(in)  :: T        !< Potential temperature referenced to the surface [C ~> degC]
  real,           intent(in)  :: S        !< Salinity [S ~> ppt]
  real,           intent(in)  :: pressure !< Pressure [R L2 T-2 ~> Pa]
  real,           intent(out) :: rho      !< Density (in-situ if pressure is local) [R ~> kg m-3]
  type(EOS_type), intent(in)  :: EOS      !< Equation of state structure
  real, optional, intent(in)  :: rho_ref  !< A reference density [R ~> kg m-3]
  real, optional, intent(in)  :: scale    !< A multiplicative factor by which to scale output density in
                                          !! combination with scaling stored in EOS [various]


end subroutine calculate_density_scalar
module subroutine calculate_stanley_density_scalar(T, S, pressure, Tvar, TScov, Svar, rho, EOS, rho_ref, scale)
  real,           intent(in)  :: T        !< Potential temperature referenced to the surface [C ~> degC]
  real,           intent(in)  :: S        !< Salinity [S ~> ppt]
  real,           intent(in)  :: Tvar     !< Variance of potential temperature referenced to the surface [C2 ~> degC2]
  real,           intent(in)  :: TScov    !< Covariance of potential temperature and salinity [C S ~> degC ppt]
  real,           intent(in)  :: Svar     !< Variance of salinity [S2 ~> ppt2]
  real,           intent(in)  :: pressure !< Pressure [R L2 T-2 ~> Pa]
  real,           intent(out) :: rho      !< Density (in-situ if pressure is local) [R ~> kg m-3]
  type(EOS_type), intent(in)  :: EOS      !< Equation of state structure
  real, optional, intent(in)  :: rho_ref  !< A reference density [R ~> kg m-3].
  real, optional, intent(in)  :: scale    !< A multiplicative factor by which to scale output density in
                                          !! combination with scaling stored in EOS [various]
  ! Local variables

end subroutine calculate_stanley_density_scalar
module subroutine calculate_density_1d(T, S, pressure, rho, EOS, dom, rho_ref, scale)
  real, dimension(:),    intent(in)    :: T        !< Potential temperature referenced to the surface [C ~> degC]
  real, dimension(:),    intent(in)    :: S        !< Salinity [S ~> ppt]
  real, dimension(:),    intent(in)    :: pressure !< Pressure [R L2 T-2 ~> Pa]
  real, dimension(:),    intent(inout) :: rho      !< Density (in-situ if pressure is local) [R ~> kg m-3]
  type(EOS_type),        intent(in)    :: EOS      !< Equation of state structure
  integer, dimension(2), optional, intent(in) :: dom   !< The domain of indices to work on, taking
                                                       !! into account that arrays start at 1.
  real,                  optional, intent(in) :: rho_ref !< A reference density [R ~> kg m-3]
  real,                  optional, intent(in) :: scale !< A multiplicative factor by which to scale density
                                                   !! in combination with scaling stored in EOS [various]
  ! Local variables

end subroutine calculate_density_1d
module subroutine calculate_stanley_density_1d(T, S, pressure, Tvar, TScov, Svar, rho, EOS, dom, rho_ref, scale)
  real, dimension(:),    intent(in)    :: T        !< Potential temperature referenced to the surface [C ~> degC]
  real, dimension(:),    intent(in)    :: S        !< Salinity [S ~> ppt]
  real, dimension(:),    intent(in)    :: pressure !< Pressure [R L2 T-2 ~> Pa]
  real, dimension(:),    intent(in)    :: Tvar     !< Variance of potential temperature [C2 ~> degC2]
  real, dimension(:),    intent(in)    :: TScov    !< Covariance of potential temperature and salinity [C S ~> degC ppt]
  real, dimension(:),    intent(in)    :: Svar     !< Variance of salinity [S2 ~> ppt2]
  real, dimension(:),    intent(inout) :: rho      !< Density (in-situ if pressure is local) [R ~> kg m-3]
  type(EOS_type),        intent(in)    :: EOS      !< Equation of state structure
  integer, dimension(2), optional, intent(in) :: dom   !< The domain of indices to work on, taking
                                                       !! into account that arrays start at 1.
  real,                  optional, intent(in) :: rho_ref !< A reference density [R ~> kg m-3]
  real,                  optional, intent(in) :: scale !< A multiplicative factor by which to scale density
                                                   !! in combination with scaling stored in EOS [various]
  ! Local variables

end subroutine calculate_stanley_density_1d
module subroutine calculate_spec_vol_array(T, S, pressure, specvol, start, npts, EOS, spv_ref, scale)
  real, dimension(:), intent(in)    :: T        !< potential temperature relative to the surface [degC]
  real, dimension(:), intent(in)    :: S        !< salinity [ppt]
  real, dimension(:), intent(in)    :: pressure !< pressure [Pa]
  real, dimension(:), intent(inout) :: specvol  !< in situ specific volume [kg m-3]
  integer,            intent(in)    :: start    !< the starting point in the arrays.
  integer,            intent(in)    :: npts     !< the number of values to calculate.
  type(EOS_type),     intent(in)    :: EOS      !< Equation of state structure
  real,     optional, intent(in)    :: spv_ref  !< A reference specific volume [m3 kg-1]
  real,     optional, intent(in)    :: scale    !< A multiplicative factor by which to scale specific
                                                !! volume in combination with scaling stored in EOS [various]


end subroutine calculate_spec_vol_array
module subroutine calc_spec_vol_scalar(T, S, pressure, specvol, EOS, spv_ref, scale)
  real,           intent(in)  :: T        !< Potential temperature referenced to the surface [C ~> degC]
  real,           intent(in)  :: S        !< Salinity [S ~> ppt]
  real,           intent(in)  :: pressure !< Pressure [R L2 T-2 ~> Pa]
  real,           intent(out) :: specvol  !< In situ or potential specific volume [R-1 ~> m3 kg-1]
                                          !! or other units determined by the scale argument
  type(EOS_type), intent(in)  :: EOS      !< Equation of state structure
  real, optional, intent(in)  :: spv_ref  !< A reference specific volume [R-1 ~> m3 kg-1]
  real, optional, intent(in)  :: scale    !< A multiplicative factor by which to scale specific
                                          !! volume in combination with scaling stored in EOS [various]


end subroutine calc_spec_vol_scalar
module subroutine calc_spec_vol_1d(T, S, pressure, specvol, EOS, dom, spv_ref, scale)
  real, dimension(:),    intent(in)    :: T        !< Potential temperature referenced to the surface [C ~> degC]
  real, dimension(:),    intent(in)    :: S        !< Salinity [S ~> ppt]
  real, dimension(:),    intent(in)    :: pressure !< Pressure [R L2 T-2 ~> Pa]
  real, dimension(:),    intent(inout) :: specvol  !< In situ specific volume [R-1 ~> m3 kg-1]
  type(EOS_type),        intent(in)    :: EOS      !< Equation of state structure
  integer, dimension(2), optional, intent(in) :: dom   !< The domain of indices to work on, taking
                                                       !! into account that arrays start at 1.
  real,                  optional, intent(in) :: spv_ref !< A reference specific volume [R-1 ~> m3 kg-1]
  real,                  optional, intent(in) :: scale !< A multiplicative factor by which to scale
                                                       !! output specific volume in combination with
                                                       !! scaling stored in EOS [various]
  ! Local variables

end subroutine calc_spec_vol_1d
module subroutine calculate_TFreeze_scalar(S, pressure, T_fr, EOS, pres_scale, scale_from_EOS)
  real,           intent(in)  :: S    !< Salinity, [ppt] or [S ~> ppt] depending on scale_from_EOS
  real,           intent(in)  :: pressure !< Pressure, in [Pa] or [R L2 T-2 ~> Pa] depending on
                                      !! pres_scale or scale_from_EOS
  real,           intent(out) :: T_fr !< Freezing point potential temperature referenced to the
                                      !! surface [degC] or [C ~> degC] depending on scale_from_EOS
  type(EOS_type), intent(in)  :: EOS  !< Equation of state structure
  real, optional, intent(in)  :: pres_scale  !< A multiplicative factor to convert pressure
                                      !! into Pa [Pa T2 R-1 L-2 ~> 1].
  logical, optional, intent(in)  :: scale_from_EOS !< If present true use the dimensional scaling
                                      !! factors stored in EOS.  Omission is the same .false.

  ! Local variables

end subroutine calculate_TFreeze_scalar
module subroutine calculate_TFreeze_array(S, pressure, T_fr, start, npts, EOS, pres_scale)
  real, dimension(:), intent(in)    :: S        !< Salinity [ppt]
  real, dimension(:), intent(in)    :: pressure !< Pressure, in [Pa] or [R L2 T-2 ~> Pa] depending on pres_scale
  real, dimension(:), intent(inout) :: T_fr     !< Freezing point, either potential temperature referenced to the
                                                !! surface or conservative temperature depending on settings [degC]
  integer,            intent(in)    :: start    !< Starting index within the array
  integer,            intent(in)    :: npts     !< The number of values to calculate
  type(EOS_type),     intent(in)    :: EOS      !< Equation of state structure
  real,     optional, intent(in)    :: pres_scale !< A multiplicative factor to convert pressure
                                                !! into Pa [Pa T2 R-1 L-2 ~> 1].

  ! Local variables

end subroutine calculate_TFreeze_array
module subroutine calculate_TFreeze_1d(S, pressure, T_fr, EOS, dom)
  real, dimension(:), intent(in)    :: S        !< Salinity [S ~> ppt]
  real, dimension(:), intent(in)    :: pressure !< Pressure [R L2 T-2 ~> Pa]
  real, dimension(:), intent(inout) :: T_fr     !< Freezing point, either potential temperature referenced to the
                                                !! surface or conservative temperature depending on settings
                                                !! [C ~> degC]
  type(EOS_type),     intent(in)    :: EOS      !< Equation of state structure
  integer, dimension(2), optional, intent(in) :: dom   !< The domain of indices to work on, taking
                                                       !! into account that arrays start at 1.

  ! Local variables

end subroutine calculate_TFreeze_1d
module subroutine calculate_density_derivs_array(T, S, pressure, drho_dT, drho_dS, start, npts, EOS, scale)
  real, dimension(:), intent(in)    :: T        !< Potential temperature referenced to the surface [degC]
  real, dimension(:), intent(in)    :: S        !< Salinity [ppt]
  real, dimension(:), intent(in)    :: pressure !< Pressure [Pa]
  real, dimension(:), intent(inout) :: drho_dT  !< The partial derivative of density with potential
                                                !! temperature [kg m-3 degC-1] or other units determined
                                                !! by the optional scale argument
  real, dimension(:), intent(inout) :: drho_dS  !< The partial derivative of density with salinity,
                                                !! in [kg m-3 ppt-1] or other units determined
                                                !! by the optional scale argument
  integer,            intent(in)    :: start    !< Starting index within the array
  integer,            intent(in)    :: npts     !< The number of values to calculate
  type(EOS_type),     intent(in)    :: EOS      !< Equation of state structure
  real,     optional, intent(in)    :: scale    !< A multiplicative factor by which to scale density
                                                !! in combination with scaling stored in EOS [various]

  ! Local variables

end subroutine calculate_density_derivs_array
module subroutine calculate_density_derivs_1d(T, S, pressure, drho_dT, drho_dS, EOS, dom, scale)
  real, dimension(:),    intent(in)    :: T        !< Potential temperature referenced to the surface [C ~> degC]
  real, dimension(:),    intent(in)    :: S        !< Salinity [S ~> ppt]
  real, dimension(:),    intent(in)    :: pressure !< Pressure [R L2 T-2 ~> Pa]
  real, dimension(:),    intent(inout) :: drho_dT  !< The partial derivative of density with potential
                                                   !! temperature [R C-1 ~> kg m-3 degC-1]
  real, dimension(:),    intent(inout) :: drho_dS  !< The partial derivative of density with salinity
                                                   !! [R S-1 ~> kg m-3 ppt-1]
  type(EOS_type),        intent(in)    :: EOS      !< Equation of state structure
  integer, dimension(2), optional, intent(in) :: dom   !< The domain of indices to work on, taking
                                                       !! into account that arrays start at 1.
  real,                  optional, intent(in) :: scale !< A multiplicative factor by which to scale density
                                                       !! in combination with scaling stored in EOS [various]
  ! Local variables

end subroutine calculate_density_derivs_1d
module subroutine calculate_density_derivs_scalar(T, S, pressure, drho_dT, drho_dS, EOS, scale)
  real,           intent(in)  :: T !< Potential temperature referenced to the surface [C ~> degC]
  real,           intent(in)  :: S !< Salinity [S ~> ppt]
  real,           intent(in)  :: pressure !< Pressure [R L2 T-2 ~> Pa]
  real,           intent(out) :: drho_dT !< The partial derivative of density with potential
                                         !! temperature [R C-1 ~> kg m-3 degC-1] or other
                                         !! units determined by the optional scale argument
  real,           intent(out) :: drho_dS !< The partial derivative of density with salinity,
                                         !! in [R S-1 ~> kg m-3 ppt-1] or other units
                                         !! determined by the optional scale argument
  type(EOS_type), intent(in)  :: EOS     !< Equation of state structure
  real, optional, intent(in)  :: scale   !< A multiplicative factor by which to scale density
                                         !! in combination with scaling stored in EOS [various]
  ! Local variables

end subroutine calculate_density_derivs_scalar
module subroutine calculate_density_second_derivs_1d(T, S, pressure, drho_dS_dS, drho_dS_dT, drho_dT_dT, &
                                              drho_dS_dP, drho_dT_dP, EOS, dom, scale)
  real, dimension(:), intent(in)  :: T !< Potential temperature referenced to the surface [C ~> degC]
  real, dimension(:), intent(in)  :: S !< Salinity [S ~> ppt]
  real, dimension(:), intent(in)  :: pressure   !< Pressure [R L2 T-2 ~> Pa]
  real, dimension(:), intent(inout) :: drho_dS_dS !< Partial derivative of beta with respect to S
                                                  !! [R S-2 ~> kg m-3 ppt-2]
  real, dimension(:), intent(inout) :: drho_dS_dT !< Partial derivative of beta with respect to T
                                                  !! [R S-1 C-1 ~> kg m-3 ppt-1 degC-1]
  real, dimension(:), intent(inout) :: drho_dT_dT !< Partial derivative of alpha with respect to T
                                                  !! [R C-2 ~> kg m-3 degC-2]
  real, dimension(:), intent(inout) :: drho_dS_dP !< Partial derivative of beta with respect to pressure
                                                  !! [T2 S-1 L-2 ~> kg m-3 ppt-1 Pa-1]
  real, dimension(:), intent(inout) :: drho_dT_dP !< Partial derivative of alpha with respect to pressure
                                                  !! [T2 C-1 L-2 ~> kg m-3 degC-1 Pa-1]
  type(EOS_type),     intent(in)    :: EOS        !< Equation of state structure
  integer, dimension(2), optional, intent(in) :: dom   !< The domain of indices to work on, taking
                                                  !! into account that arrays start at 1.
  real,     optional, intent(in)    :: scale      !< A multiplicative factor by which to scale density
                                                  !! in combination with scaling stored in EOS [various]
  ! Local variables

end subroutine calculate_density_second_derivs_1d
module subroutine calculate_density_second_derivs_scalar(T, S, pressure, drho_dS_dS, drho_dS_dT, drho_dT_dT, &
                                                  drho_dS_dP, drho_dT_dP, EOS, scale)
  real, intent(in)  :: T !< Potential temperature referenced to the surface [C ~> degC]
  real, intent(in)  :: S !< Salinity [S ~> ppt]
  real, intent(in)  :: pressure   !< Pressure [R L2 T-2 ~> Pa]
  real, intent(out) :: drho_dS_dS !< Partial derivative of beta with respect to S
                                  !! [R S-2 ~> kg m-3 ppt-2]
  real, intent(out) :: drho_dS_dT !< Partial derivative of beta with respect to T
                                  !! [R S-1 C-1 ~> kg m-3 ppt-1 degC-1]
  real, intent(out) :: drho_dT_dT !< Partial derivative of alpha with respect to T
                                  !! [R C-2 ~> kg m-3 degC-2]
  real, intent(out) :: drho_dS_dP !< Partial derivative of beta with respect to pressure
                                  !! [T2 S-1 L-2 ~> kg m-3 ppt-1 Pa-1]
  real, intent(out) :: drho_dT_dP !< Partial derivative of alpha with respect to pressure
                                  !! [T2 C-1 L-2 ~> kg m-3 degC-1 Pa-1]
  type(EOS_type), intent(in) :: EOS !< Equation of state structure
  real, optional, intent(in) :: scale !< A multiplicative factor by which to scale density
                                  !! in combination with scaling stored in EOS [various]
  ! Local variables

end subroutine calculate_density_second_derivs_scalar
module subroutine calculate_spec_vol_derivs_array(T, S, pressure, dSV_dT, dSV_dS, start, npts, EOS)
  real, dimension(:), intent(in)  :: T !< Potential temperature referenced to the surface [degC]
  real, dimension(:), intent(in)  :: S !< Salinity [ppt]
  real, dimension(:), intent(in)  :: pressure !< Pressure [Pa]
  real, dimension(:), intent(inout) :: dSV_dT !< The partial derivative of specific volume with potential
                                              !! temperature [m3 kg-1 degC-1]
  real, dimension(:), intent(inout) :: dSV_dS !< The partial derivative of specific volume with salinity
                                              !! [m3 kg-1 ppt-1]
  integer,            intent(in)  :: start  !< Starting index within the array
  integer,            intent(in)  :: npts   !< The number of values to calculate
  type(EOS_type),     intent(in)  :: EOS    !< Equation of state structure

end subroutine calculate_spec_vol_derivs_array
module subroutine calc_spec_vol_derivs_1d(T, S, pressure, dSV_dT, dSV_dS, EOS, dom, scale)
  real, dimension(:), intent(in)    :: T        !< Potential temperature referenced to the surface [C ~> degC]
  real, dimension(:), intent(in)    :: S        !< Salinity [S ~> ppt]
  real, dimension(:), intent(in)    :: pressure !< Pressure [R L2 T-2 ~> Pa]
  real, dimension(:), intent(inout) :: dSV_dT   !< The partial derivative of specific volume with potential
                                                !! temperature [R-1 C-1 ~> m3 kg-1 degC-1]
  real, dimension(:), intent(inout) :: dSV_dS   !< The partial derivative of specific volume with salinity
                                                !! [R-1 S-1 ~> m3 kg-1 ppt-1]
  type(EOS_type),     intent(in)    :: EOS      !< Equation of state structure
  integer, dimension(2), optional, intent(in) :: dom   !< The domain of indices to work on, taking
                                                       !! into account that arrays start at 1.
  real,                  optional, intent(in) :: scale !< A multiplicative factor by which to scale specific
                                                !! volume in combination with scaling stored in EOS [various]

  ! Local variables

end subroutine calc_spec_vol_derivs_1d
module subroutine calculate_compress_1d(T, S, pressure, rho, drho_dp, EOS, dom)
  real, dimension(:), intent(in)    :: T        !< Potential temperature referenced to the surface [C ~> degC]
  real, dimension(:), intent(in)    :: S        !< Salinity [S ~> ppt]
  real, dimension(:), intent(in)    :: pressure !< Pressure [R L2 T-2 ~> Pa]
  real, dimension(:), intent(inout) :: rho      !< In situ density [R ~> kg m-3]
  real, dimension(:), intent(inout) :: drho_dp  !< The partial derivative of density with pressure
                                                !! (also the inverse of the square of sound speed)
                                                !! [T2 L-2 ~> s2 m-2]
  type(EOS_type),     intent(in)  :: EOS      !< Equation of state structure
  integer, dimension(2), optional, intent(in) :: dom   !< The domain of indices to work on, taking
                                                       !! into account that arrays start at 1.

  ! Local variables

end subroutine calculate_compress_1d
module subroutine calculate_compress_scalar(T, S, pressure, rho, drho_dp, EOS)
  real, intent(in)        :: T        !< Potential temperature referenced to the surface [C ~> degC]
  real, intent(in)        :: S        !< Salinity [S ~> ppt]
  real, intent(in)        :: pressure !< Pressure [R L2 T-2 ~> Pa]
  real, intent(out)       :: rho      !< In situ density [R ~> kg m-3]
  real, intent(out)       :: drho_dp  !< The partial derivative of density with pressure (also the
                                      !! inverse of the square of sound speed) [T2 L-2 ~> s2 m-2]
  type(EOS_type), intent(in) :: EOS   !< Equation of state structure

  ! Local variables
  ! These arrays use the same units as their counterparts in calculate_compress_1d.
                              ! inverse of the square of sound speed) in a 1d array [T2 L-2 ~> s2 m-2]

end subroutine calculate_compress_scalar
module subroutine average_specific_vol(T, S, p_t, dp, SpV_avg, EOS, dom, scale)
  real, dimension(:),    intent(in)    :: T   !< Potential temperature referenced to the surface [C ~> degC]
  real, dimension(:),    intent(in)    :: S   !< Salinity [S ~> ppt]
  real, dimension(:),    intent(in)    :: p_t !< Pressure at the top of the layer [R L2 T-2 ~> Pa]
  real, dimension(:),    intent(in)    :: dp  !< Pressure change in the layer [R L2 T-2 ~> Pa]
  real, dimension(:),    intent(inout) :: SpV_avg !< The vertical average specific volume
                                              !! in the layer [R-1 ~> m3 kg-1]
  type(EOS_type),        intent(in)    :: EOS !< Equation of state structure
  integer, dimension(2), optional, intent(in) :: dom   !< The domain of indices to work on, taking
                                                       !! into account that arrays start at 1.
  real,                  optional, intent(in) :: scale !< A multiplicative factor by which to scale
                                                       !! output specific volume in combination with
                                                       !! scaling stored in EOS [various]

  ! Local variables

end subroutine average_specific_vol
module subroutine EoS_fit_range(EOS, T_min, T_max, S_min, S_max, p_min, p_max)
  type(EOS_type), intent(in) :: EOS   !< Equation of state structure
  real, optional, intent(out) :: T_min !< The minimum temperature over which this EoS is fitted [degC]
  real, optional, intent(out) :: T_max !< The maximum temperature over which this EoS is fitted [degC]
  real, optional, intent(out) :: S_min !< The minimum salinity over which this EoS is fitted [ppt]
  real, optional, intent(out) :: S_max !< The maximum salinity over which this EoS is fitted [ppt]
  real, optional, intent(out) :: p_min !< The minimum pressure over which this EoS is fitted [Pa]
  real, optional, intent(out) :: p_max !< The maximum pressure over which this EoS is fitted [Pa]

end subroutine EoS_fit_range
module function EOS_domain(HI, halo) result(EOSdom)
  type(hor_index_type), intent(in)  :: HI    !< The horizontal index structure
  integer,    optional, intent(in)  :: halo  !< The halo size to work on; missing is equivalent to 0.
  integer, dimension(2) :: EOSdom   !< The index domain that the EOS will work on, taking into account
                                    !! that the arrays inside the EOS routines will start at 1.

  ! Local variables

end function EOS_domain
module subroutine analytic_int_specific_vol_dp(T, S, p_t, p_b, alpha_ref, HI, EOS, &
                               dza, intp_dza, intx_dza, inty_dza, halo_size, &
                               bathyP, P_surf, dP_tiny, MassWghtInterp)
  type(hor_index_type), intent(in)  :: HI  !< The horizontal index structure
  real, dimension(HI%isd:HI%ied,HI%jsd:HI%jed), &
                        intent(in)  :: T   !< Potential temperature referenced to the surface [C ~> degC]
  real, dimension(HI%isd:HI%ied,HI%jsd:HI%jed), &
                        intent(in)  :: S   !< Salinity [S ~> ppt]
  real, dimension(HI%isd:HI%ied,HI%jsd:HI%jed), &
                        intent(in)  :: p_t !< Pressure at the top of the layer [R L2 T-2 ~> Pa]
  real, dimension(HI%isd:HI%ied,HI%jsd:HI%jed), &
                        intent(in)  :: p_b !< Pressure at the bottom of the layer [R L2 T-2 ~> Pa]
  real,                 intent(in)  :: alpha_ref !< A mean specific volume that is subtracted out
                            !! to reduce the magnitude of each of the integrals [R-1 ~> m3 kg-1]
                            !! The calculation is mathematically identical with different values of
                            !! alpha_ref, but this reduces the effects of roundoff.
  type(EOS_type),       intent(in)  :: EOS !< Equation of state structure
  real, dimension(HI%isd:HI%ied,HI%jsd:HI%jed), &
                        intent(inout) :: dza !< The change in the geopotential anomaly across
                            !! the layer [L2 T-2 ~> m2 s-2]
  real, dimension(HI%isd:HI%ied,HI%jsd:HI%jed), &
              optional, intent(inout) :: intp_dza !< The integral in pressure through the layer of the
                            !! geopotential anomaly relative to the anomaly at the bottom of the
                            !! layer [R L4 T-4 ~> Pa m2 s-2]
  real, dimension(HI%IsdB:HI%IedB,HI%jsd:HI%jed), &
              optional, intent(inout) :: intx_dza !< The integral in x of the difference between the
                            !! geopotential anomaly at the top and bottom of the layer divided by
                            !! the x grid spacing [L2 T-2 ~> m2 s-2]
  real, dimension(HI%isd:HI%ied,HI%JsdB:HI%JedB), &
              optional, intent(inout) :: inty_dza !< The integral in y of the difference between the
                            !! geopotential anomaly at the top and bottom of the layer divided by
                            !! the y grid spacing [L2 T-2 ~> m2 s-2]
  integer,    optional, intent(in)  :: halo_size !< The width of halo points on which to calculate dza.
  real, dimension(HI%isd:HI%ied,HI%jsd:HI%jed), &
              optional, intent(in)  :: bathyP  !< The pressure at the bathymetry [R L2 T-2 ~> Pa]
  real, dimension(HI%isd:HI%ied,HI%jsd:HI%jed), &
              optional, intent(in)  :: P_surf !< The pressure at the ocean surface [R L2 T-2 ~> Pa]
  real,       optional, intent(in)  :: dP_tiny !< A miniscule pressure change with
                            !! the same units as p_t [R L2 T-2 ~> Pa]
  integer,    optional, intent(in)  :: MassWghtInterp !< A flag indicating whether and how to use
                            !! mass weighting to interpolate T/S in integrals

  ! Local variables

  ! We should never reach this point with quadrature. EOS_quadrature indicates that numerical
  ! integration be used instead of analytic. This is a safety check.
end subroutine analytic_int_specific_vol_dp
module subroutine analytic_int_density_dz(T, S, z_t, z_b, rho_ref, rho_0, G_e, HI, EOS, dpa, &
                          intz_dpa, intx_dpa, inty_dpa, bathyT, SSH, dz_neglect, MassWghtInterp, Z_0p)
  type(hor_index_type), intent(in)  :: HI !< Ocean horizontal index structure
  real, dimension(HI%isd:HI%ied,HI%jsd:HI%jed), &
                        intent(in)  :: T   !< Potential temperature referenced to the surface [C ~> degC]
  real, dimension(HI%isd:HI%ied,HI%jsd:HI%jed), &
                        intent(in)  :: S   !< Salinity [S ~> ppt]
  real, dimension(HI%isd:HI%ied,HI%jsd:HI%jed), &
                        intent(in)  :: z_t !< Height at the top of the layer in depth units [Z ~> m]
  real, dimension(HI%isd:HI%ied,HI%jsd:HI%jed), &
                        intent(in)  :: z_b !< Height at the bottom of the layer [Z ~> m]
  real,                 intent(in)  :: rho_ref !< A mean density [R ~> kg m-3], that is
                                           !! subtracted out to reduce the magnitude of each of the
                                           !! integrals.
  real,                 intent(in)  :: rho_0 !< A density [R ~> kg m-3], that is used
                                           !! to calculate the pressure (as p~=-z*rho_0*G_e)
                                           !! used in the equation of state.
  real,                 intent(in)  :: G_e !< The Earth's gravitational acceleration
                                           !! [L2 Z-1 T-2 ~> m s-2]
  type(EOS_type),       intent(in)  :: EOS !< Equation of state structure
  real, dimension(HI%isd:HI%ied,HI%jsd:HI%jed), &
                      intent(inout) :: dpa !< The change in the pressure anomaly
                                           !! across the layer [R L2 T-2 ~> Pa]
  real, dimension(HI%isd:HI%ied,HI%jsd:HI%jed), &
            optional, intent(inout) :: intz_dpa !< The integral through the thickness of the
                                           !! layer of the pressure anomaly relative to the
                                           !! anomaly at the top of the layer [R L2 Z T-2 ~> Pa m]
  real, dimension(HI%IsdB:HI%IedB,HI%jsd:HI%jed), &
            optional, intent(inout) :: intx_dpa !< The integral in x of the difference between
                                          !! the pressure anomaly at the top and bottom of the
                                          !! layer divided by the x grid spacing [R L2 T-2 ~> Pa]
  real, dimension(HI%isd:HI%ied,HI%JsdB:HI%JedB), &
            optional, intent(inout) :: inty_dpa !< The integral in y of the difference between
                                          !! the pressure anomaly at the top and bottom of the
                                          !! layer divided by the y grid spacing [R L2 T-2 ~> Pa]
  real, dimension(HI%isd:HI%ied,HI%jsd:HI%jed), &
              optional, intent(in)  :: bathyT !< The depth of the bathymetry [Z ~> m]
  real, dimension(HI%isd:HI%ied,HI%jsd:HI%jed), &
              optional, intent(in)  :: SSH !< The sea surface height [Z ~> m]
  real,       optional, intent(in)  :: dz_neglect !< A miniscule thickness change [Z ~> m]
  integer,    optional, intent(in)  :: MassWghtInterp !< A flag indicating whether and how to use
                                          !! mass weighting to interpolate T/S in integrals
  real, dimension(HI%isd:HI%ied,HI%jsd:HI%jed), &
              optional, intent(in)  :: Z_0p !< The height at which the pressure is 0 [Z ~> m]

  ! Local variables
                     ! desired units [R m3 kg-1 ~> 1]

  ! We should never reach this point with quadrature. EOS_quadrature indicates that numerical
  ! integration be used instead of analytic. This is a safety check.
end subroutine analytic_int_density_dz
logical module function query_compressible(EOS)
  type(EOS_type), intent(in) :: EOS !< Equation of state structure

end function query_compressible
module function get_EOS_name(id) result (eos_name)
  integer,        optional, intent(in) :: id !< Enumerated ID
  character(:), allocatable :: eos_name !< The name of the EOS

end function get_EOS_name
module subroutine EOS_init(param_file, EOS, US, use_conT_absS)
  type(param_file_type), intent(in) :: param_file !< Parameter file structure
  type(EOS_type), intent(inout)     :: EOS !< Equation of state structure
  type(unit_scale_type), intent(in) :: US  !< A dimensional unit scaling type
  logical, intent(in), optional     :: use_conT_absS !< True if the model is formulated for
                                                     !! conservative temp and absolute salinity
  optional :: US
  ! Local variables

  ! Read all relevant parameters and write them to the model log.
end subroutine EOS_init
module subroutine EOS_manual_init(EOS, form_of_EOS, form_of_TFreeze, EOS_quadrature, Compressible, &
                           Rho_T0_S0, drho_dT, dRho_dS, dRho_dp, TFr_S0_P0, dTFr_dS, dTFr_dp, &
                           use_Wright_2nd_deriv_bug)
  type(EOS_type),    intent(inout) :: EOS !< Equation of state structure
  integer, optional, intent(in) :: form_of_EOS !< A coded integer indicating the equation of state to use.
  integer, optional, intent(in) :: form_of_TFreeze !< A coded integer indicating the expression for
                                       !! the potential temperature of the freezing point.
  logical, optional, intent(in) :: EOS_quadrature !< If true, always use the generic (quadrature)
                                       !! code for the integrals of density.
  logical, optional, intent(in) :: Compressible  !< If true, in situ density is a function of pressure.
  real   , optional, intent(in) :: Rho_T0_S0 !< Density at T=0 degC and S=0 ppt [kg m-3]
  real   , optional, intent(in) :: drho_dT   !< Partial derivative of density with temperature
                                             !! in [kg m-3 degC-1]
  real   , optional, intent(in) :: dRho_dS   !< Partial derivative of density with salinity
                                             !! in [kg m-3 ppt-1]
  real   , optional, intent(in) :: dRho_dp   !< Partial derivative of density with pressure
                                             !! in [s2 m-2]
  real   , optional, intent(in) :: TFr_S0_P0 !< The freezing potential temperature at S=0, P=0 [degC]
  real   , optional, intent(in) :: dTFr_dS   !< The derivative of freezing point with salinity
                                             !! in [degC ppt-1]
  real   , optional, intent(in) :: dTFr_dp   !< The derivative of freezing point with pressure
                                             !! in [degC Pa-1]
  logical, optional, intent(in) :: use_Wright_2nd_deriv_bug !< Allow the Wright 2nd deriv bug

end subroutine EOS_manual_init
module subroutine convert_temp_salt_for_TEOS10(T, S, HI, kd, mask_z, EOS)
  integer,               intent(in)    :: kd  !< The number of layers to work on
  type(hor_index_type),  intent(in)    :: HI       !< The horizontal index structure
  real, dimension(HI%isd:HI%ied,HI%jsd:HI%jed,kd), &
                         intent(inout) :: T   !< Potential temperature referenced to the surface [C ~> degC]
  real, dimension(HI%isd:HI%ied,HI%jsd:HI%jed,kd), &
                         intent(inout) :: S   !< Salinity [S ~> ppt]
  real, dimension(HI%isd:HI%ied,HI%jsd:HI%jed,kd), &
                         intent(in)    :: mask_z !< 3d mask regulating which points to convert [nondim]
  type(EOS_type),        intent(in)    :: EOS !< Equation of state structure

                                    ! practical salinity to reference salinity [PSU ppt-1]

end subroutine convert_temp_salt_for_TEOS10
module subroutine cons_temp_to_pot_temp(T, S, poTemp, EOS, dom, scale)
  real, dimension(:), intent(in)    :: T        !< Conservative temperature [C ~> degC]
  real, dimension(:), intent(in)    :: S        !< Absolute salinity [S ~> ppt]
  real, dimension(:), intent(inout) :: poTemp   !< The potential temperature with a reference pressure
                                                !! of 0 Pa, [C ~> degC]
  type(EOS_type),     intent(in)    :: EOS      !< Equation of state structure
  integer, dimension(2), optional, intent(in) :: dom  !< The domain of indices to work on, taking
                                                !! into account that arrays start at 1.
  real,     optional, intent(in)    :: scale    !< A multiplicative factor by which to scale the output
                                                !! potential temperature in place of with scaling stored
                                                !! in EOS.  A value of 1.0 returns temperatures in [degC],
                                                !! while the default is equivalent to EOS%degC_to_C.

  ! Local variables

end subroutine cons_temp_to_pot_temp
module subroutine pot_temp_to_cons_temp(T, S, consTemp, EOS, dom, scale)
  real, dimension(:), intent(in)    :: T        !< Potential temperature [C ~> degC]
  real, dimension(:), intent(in)    :: S        !< Absolute salinity [S ~> ppt]
  real, dimension(:), intent(inout) :: consTemp !< The conservative temperature [C ~> degC]
  type(EOS_type),     intent(in)    :: EOS      !< Equation of state structure
  integer, dimension(2), optional, intent(in) :: dom  !< The domain of indices to work on, taking
                                                !! into account that arrays start at 1.
  real,     optional, intent(in)    :: scale    !< A multiplicative factor by which to scale the output
                                                !! potential temperature in place of with scaling stored
                                                !! in EOS.  A value of 1.0 returns temperatures in [degC],
                                                !! while the default is equivalent to EOS%degC_to_C.

  ! Local variables

end subroutine pot_temp_to_cons_temp
module subroutine abs_saln_to_prac_saln(S, prSaln, EOS, dom, scale)
  real, dimension(:), intent(in)    :: S        !< Absolute salinity [S ~> ppt]
  real, dimension(:), intent(inout) :: prSaln   !< Practical salinity [S ~> PSU]
  type(EOS_type),     intent(in)    :: EOS      !< Equation of state structure
  integer, dimension(2), optional, intent(in) :: dom  !< The domain of indices to work on, taking
                                                !! into account that arrays start at 1.
  real,     optional, intent(in)    :: scale    !< A multiplicative factor by which to scale the output
                                                !! practical salinities in place of with scaling stored
                                                !! in EOS.  A value of 1.0 returns salinities in [PSU],
                                                !! while the default is equivalent to EOS%ppt_to_S.

  ! Local variables
                                    ! reference salinity to practical salinity [PSU ppt-1]

end subroutine abs_saln_to_prac_saln
module subroutine prac_saln_to_abs_saln(S, absSaln, EOS, dom, scale)
  real, dimension(:), intent(in)    :: S        !< Practical salinity [S ~> PSU]
  real, dimension(:), intent(inout) :: absSaln  !< Absolute salinity [S ~> ppt]
  type(EOS_type),     intent(in)    :: EOS      !< Equation of state structure
  integer, dimension(2), optional, intent(in) :: dom  !< The domain of indices to work on, taking
                                                !! into account that arrays start at 1.
  real,     optional, intent(in)    :: scale    !< A multiplicative factor by which to scale the output
                                                !! absolute salnities in place of with scaling stored
                                                !! in EOS.  A value of 1.0 returns salinities in [ppt],
                                                !! while the default is equivalent to EOS%ppt_to_S.

  ! Local variables
                                    ! practical salinity to reference salinity [PSU ppt-1]

end subroutine prac_saln_to_abs_saln
logical module function EOS_quadrature(EOS)
  type(EOS_type), intent(in) :: EOS   !< Equation of state structure

end function EOS_quadrature
logical module function EOS_unit_tests(verbose)
  logical, intent(in) :: verbose !< If true, write results to stdout
  ! Local variables

end function EOS_unit_tests
logical module function test_TS_conversion_consistency(T_cons, S_abs, T_pot, S_prac, EOS, verbose) &
                                      result(inconsistent)
  real,              intent(in) :: T_cons    !< Conservative temperature [degC]
  real,              intent(in) :: S_abs     !< Absolute salinity [g kg-1]
  real,              intent(in) :: T_pot     !< Potential temperature [degC]
  real,              intent(in) :: S_prac    !< Practical salinity [PSU]
  type(EOS_type),    intent(in) :: EOS      !< Equation of state structure
  logical,           intent(in) :: verbose  !< If true, write results to stdout

  ! Local variables

end function test_TS_conversion_consistency
logical module function test_TFr_consistency(S_test, p_test, EOS, verbose, EOS_name, TFr_check) &
                                      result(inconsistent)
  real,              intent(in) :: S_test   !< Salinity or absolute salinity [S ~> ppt]
  real,              intent(in) :: p_test   !< Pressure [R L2 T-2 ~> Pa]
  type(EOS_type),    intent(in) :: EOS      !< Equation of state structure
  logical,           intent(in) :: verbose  !< If true, write results to stdout
  character(len=*),  intent(in) :: EOS_name !< A name used in error messages to describe the EoS
  real,    optional, intent(in) :: TFr_check  !< A check value for the Freezing point [C ~> degC]

  ! Local variables
  ! real :: tol        ! The nondimensional tolerance from roundoff [nondim]

end function test_TFr_consistency
module subroutine write_check_msg(var_name, val, val_chk, val_tol, test_OK)
  character(len=*), intent(in) :: var_name !< The name of the variable being tested.
  real,             intent(in) :: val      !< The value being checked [various]
  real,             intent(in) :: val_chk  !< The value being checked [various]
  real,             intent(in) :: val_tol  !< The value being checked [various]
  logical,          intent(in) :: test_OK  !< True if the values are within their tolerance


end subroutine write_check_msg
logical module function test_EOS_consistency(T_test, S_test, p_test, EOS, verbose, &
                                      EOS_name, rho_check, spv_check, skip_2nd, avg_Sv_check) result(inconsistent)
  real,             intent(in) :: T_test   !< Potential temperature or conservative temperature [C ~> degC]
  real,             intent(in) :: S_test   !< Salinity or absolute salinity [S ~> ppt]
  real,             intent(in) :: p_test   !< Pressure [R L2 T-2 ~> Pa]
  type(EOS_type),   intent(in) :: EOS      !< Equation of state structure
  logical,          intent(in) :: verbose  !< If true, write results to stdout
  character(len=*), intent(in) :: EOS_name !< A name used in error messages to describe the EoS
  real,   optional, intent(in) :: rho_check  !< A check value for the density [R ~> kg m-3]
  real,   optional, intent(in) :: spv_check  !< A check value for the specific volume [R-1 ~> m3 kg-1]
  logical, optional, intent(in) :: skip_2nd  !< If present and true, do not check the 2nd derivatives.
  logical, optional, intent(in) :: avg_Sv_check !< If present and true, compare analytical and numerical
                                             !! quadrature estimates of the layer-averaged specific volume.

  ! Local variables
                                       ! perturbed points [R ~> kg m-3]
                                       ! perturbed points [R-1 ~> m3 kg-1]
                    ! temperature [R C-1 ~> kg m-3 degC-1]
                    ! in [R S-1 ~> kg m-3 ppt-1]
                    ! inverse of the square of sound speed) [T2 L-2 ~> s2 m-2]
                    ! temperature [R-1 C-1 ~> m3 kg-1 degC-1]
                    ! [R-1 S-1 ~> m3 kg-1 ppt-1]
                     ! [T2 S-1 L-2 ~> kg m-3 ppt-1 Pa-1]
                     ! [T2 C-1 L-2 ~> kg m-3 degC-1 Pa-1]

                        ! with potential temperature [R C-1 ~> kg m-3 degC-1]
                        ! with salinity [R S-1 ~> kg m-3 ppt-1]
                        ! with pressure (also the inverse of the square of sound speed) [T2 L-2 ~> s2 m-2]
                        ! specific volume with potential temperature [R-1 C-1 ~> m3 kg-1 degC-1]
                        ! specific volume with salinity [R-1 S-1 ~> m3 kg-1 ppt-1]
                            ! density with respect to salinity [R S-2 ~> kg m-3 ppt-2]
                            ! with respect to temperature and salinity [R S-1 C-1 ~> kg m-3 ppt-1 degC-1]
                            ! density with respect to temperature [R C-2 ~> kg m-3 degC-2]
                            ! with respect to salinity and pressure [T2 S-1 L-2 ~> kg m-3 ppt-1 Pa-1]
                            ! with respect to temperature and pressure [T2 C-1 L-2 ~> kg m-3 degC-1 Pa-1]
                     ! denominator in the finite difference derivative expression [nondim]
                     ! denominator in the finite difference second derivative expression [nondim]
                     ! averaged specific volume

end function test_EOS_consistency
  end interface

end module MOM_EOS
