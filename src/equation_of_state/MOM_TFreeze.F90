! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> Freezing point expressions
module MOM_TFreeze

!********+*********+*********+*********+*********+*********+*********+**
!*  The subroutines in this file determine the potential temperature   *
!* or conservative temperature at which sea-water freezes.             *
!********+*********+*********+*********+*********+*********+*********+**
use gsw_mod_toolbox, only : gsw_ct_freezing_exact

implicit none ; private

public calculate_TFreeze_linear, calculate_TFreeze_Millero, calculate_TFreeze_teos10
public calculate_TFreeze_TEOS_poly

!> Compute the freezing point potential temperature [degC] from salinity [ppt] and
!! pressure [Pa] using a simple linear expression, with coefficients passed in as arguments.
interface calculate_TFreeze_linear
  module procedure calculate_TFreeze_linear_scalar, calculate_TFreeze_linear_array
end interface calculate_TFreeze_linear

!> Compute the freezing point potential temperature [degC] from salinity [PSU] and
!! pressure [Pa] using the expression from Millero (1978) (and in appendix A of Gill 1982),
!! but with the of the pressure dependence changed from 7.53e-8 to 7.75e-8 to make this an
!! expression for potential temperature (not in situ temperature), using a
!! value that is correct at the freezing point at 35 PSU and 5e6 Pa (500 dbar).
interface calculate_TFreeze_Millero
  module procedure calculate_TFreeze_Millero_scalar, calculate_TFreeze_Millero_array
end interface calculate_TFreeze_Millero

!> Compute the freezing point conservative temperature [degC] from absolute salinity [g kg-1]
!! and pressure [Pa] using the TEOS10 package.
interface calculate_TFreeze_teos10
  module procedure calculate_TFreeze_teos10_scalar, calculate_TFreeze_teos10_array
end interface calculate_TFreeze_teos10

!> Compute the freezing point conservative temperature [degC] from absolute salinity [g kg-1] and
!! pressure [Pa] using a rescaled and refactored version of the expressions from the TEOS10 package.
interface calculate_TFreeze_TEOS_poly
  module procedure calculate_TFreeze_TEOS_poly_scalar, calculate_TFreeze_TEOS_poly_array
end interface calculate_TFreeze_TEOS_poly


  interface
module subroutine calculate_TFreeze_linear_scalar(S, pres, T_Fr, TFr_S0_P0, &
                                           dTFr_dS, dTFr_dp)
  real,  intent(in)  :: S         !< salinity [ppt].
  real,  intent(in)  :: pres      !< pressure [Pa].
  real,  intent(out) :: T_Fr      !< Freezing point potential temperature [degC].
  real,  intent(in)  :: TFr_S0_P0 !< The freezing point at S=0, p=0 [degC].
  real,  intent(in)  :: dTFr_dS   !< The derivative of freezing point with salinity,
                                  !! [degC ppt-1].
  real,  intent(in)  :: dTFr_dp   !< The derivative of freezing point with pressure,
                                  !! [degC Pa-1].

end subroutine calculate_TFreeze_linear_scalar
module subroutine calculate_TFreeze_linear_array(S, pres, T_Fr, start, npts, &
                                          TFr_S0_P0, dTFr_dS, dTFr_dp)
  real,  dimension(:), intent(in)  :: S         !< salinity [ppt].
  real,  dimension(:), intent(in)  :: pres      !< pressure [Pa].
  real,  dimension(:), intent(out) :: T_Fr      !< Freezing point potential temperature [degC].
  integer,             intent(in)  :: start     !< the starting point in the arrays.
  integer,             intent(in)  :: npts      !< the number of values to calculate.
  real,                intent(in)  :: TFr_S0_P0 !< The freezing point at S=0, p=0, [degC].
  real,                intent(in)  :: dTFr_dS   !< The derivative of freezing point with salinity,
                                                !! [degC ppt-1].
  real,                intent(in)  :: dTFr_dp   !< The derivative of freezing point with pressure,
                                                !! [degC Pa-1].

end subroutine calculate_TFreeze_linear_array
module subroutine calculate_TFreeze_Millero_scalar(S, pres, T_Fr)
  real,    intent(in)  :: S    !< Salinity [PSU]
  real,    intent(in)  :: pres !< Pressure [Pa]
  real,    intent(out) :: T_Fr !< Freezing point potential temperature [degC]

  ! Local variables

end subroutine calculate_TFreeze_Millero_scalar
module subroutine calculate_TFreeze_Millero_array(S, pres, T_Fr, start, npts)
  real,  dimension(:), intent(in)  :: S     !< Salinity [PSU].
  real,  dimension(:), intent(in)  :: pres  !< Pressure [Pa].
  real,  dimension(:), intent(out) :: T_Fr  !< Freezing point potential temperature [degC].
  integer,             intent(in)  :: start !< The starting point in the arrays.
  integer,             intent(in)  :: npts  !< The number of values to calculate.

  ! Local variables

end subroutine calculate_TFreeze_Millero_array
module subroutine calculate_TFreeze_TEOS_poly_scalar(S, pres, T_Fr)
  real,    intent(in)  :: S    !< Absolute salinity [g kg-1].
  real,    intent(in)  :: pres !< Pressure [Pa].
  real,    intent(out) :: T_Fr !< Freezing point conservative temperature [degC].

  ! Local variables

end subroutine calculate_TFreeze_TEOS_poly_scalar
module subroutine calculate_TFreeze_TEOS_poly_array(S, pres, T_Fr, start, npts)
  real, dimension(:), intent(in)  :: S     !< absolute salinity [g kg-1].
  real, dimension(:), intent(in)  :: pres  !< Pressure [Pa].
  real, dimension(:), intent(out) :: T_Fr  !< Freezing point conservative temperature [degC].
  integer,            intent(in)  :: start !< The starting point in the arrays
  integer,            intent(in)  :: npts  !< The number of values to calculate

  ! Local variables
  ! The coefficients here use the notation TFab for contributions proportional to S**a/2 * P**b.

end subroutine calculate_TFreeze_TEOS_poly_array
module subroutine calculate_TFreeze_teos10_scalar(S, pres, T_Fr)
  real,    intent(in)  :: S    !< Absolute salinity [g kg-1].
  real,    intent(in)  :: pres !< Pressure [Pa].
  real,    intent(out) :: T_Fr !< Freezing point conservative temperature [degC].

  ! Local variables

end subroutine calculate_TFreeze_teos10_scalar
module subroutine calculate_TFreeze_teos10_array(S, pres, T_Fr, start, npts)
  real, dimension(:), intent(in)  :: S     !< absolute salinity [g kg-1].
  real, dimension(:), intent(in)  :: pres  !< pressure [Pa].
  real, dimension(:), intent(out) :: T_Fr  !< Freezing point conservative temperature [degC].
  integer,            intent(in)  :: start !< the starting point in the arrays.
  integer,            intent(in)  :: npts  !< the number of values to calculate.

  ! Local variables
  ! Assume sea-water contains no dissolved air.

end subroutine calculate_TFreeze_teos10_array
  end interface

end module MOM_TFreeze
