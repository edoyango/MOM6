! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> The equation of state using the TEOS10 expressions
module MOM_EOS_TEOS10

use gsw_mod_toolbox, only : gsw_sp_from_sr, gsw_pt_from_ct, gsw_sr_from_sp, gsw_ct_from_pt
use gsw_mod_toolbox, only : gsw_rho, gsw_specvol
use gsw_mod_toolbox, only : gsw_rho_first_derivatives, gsw_specvol_first_derivatives
use gsw_mod_toolbox, only : gsw_rho_second_derivatives
use MOM_EOS_base_type, only : EOS_base

implicit none ; private

public gsw_sp_from_sr, gsw_pt_from_ct, gsw_sr_from_sp, gsw_ct_from_pt
public TEOS10_EOS

real, parameter :: Pa2db  = 1.e-4  !< The conversion factor from Pa to dbar [dbar Pa-1]

!> The EOS_base implementation of the TEOS10 equation of state
type, extends (EOS_base) :: TEOS10_EOS

contains
  !> Implementation of the in-situ density as an elemental function [kg m-3]
  procedure :: density_elem => density_elem_TEOS10
  !> Implementation of the in-situ density anomaly as an elemental function [kg m-3]
  procedure :: density_anomaly_elem => density_anomaly_elem_TEOS10
  !> Implementation of the in-situ specific volume as an elemental function [m3 kg-1]
  procedure :: spec_vol_elem => spec_vol_elem_TEOS10
  !> Implementation of the in-situ specific volume anomaly as an elemental function [m3 kg-1]
  procedure :: spec_vol_anomaly_elem => spec_vol_anomaly_elem_TEOS10
  !> Implementation of the calculation of derivatives of density
  procedure :: calculate_density_derivs_elem => calculate_density_derivs_elem_TEOS10
  !> Implementation of the calculation of second derivatives of density
  procedure :: calculate_density_second_derivs_elem => calculate_density_second_derivs_elem_TEOS10
  !> Implementation of the calculation of derivatives of specific volume
  procedure :: calculate_specvol_derivs_elem => calculate_specvol_derivs_elem_TEOS10
  !> Implementation of the calculation of compressibility
  procedure :: calculate_compress_elem => calculate_compress_elem_TEOS10
  !> Implementation of the range query function
  procedure :: EOS_fit_range => EOS_fit_range_TEOS10

end type TEOS10_EOS


  interface
real elemental module function density_elem_TEOS10(this, T, S, pressure)
  class(TEOS10_EOS), intent(in) :: this     !< This EOS
  real,              intent(in) :: T        !< Conservative temperature [degC].
  real,              intent(in) :: S        !< Absolute salinity [g kg-1].
  real,              intent(in) :: pressure !< pressure [Pa].

end function density_elem_TEOS10
real elemental module function density_anomaly_elem_TEOS10(this, T, S, pressure, rho_ref)
  class(TEOS10_EOS), intent(in) :: this     !< This EOS
  real,              intent(in) :: T        !< Conservative temperature [degC].
  real,              intent(in) :: S        !< Absolute salinity [g kg-1].
  real,              intent(in) :: pressure !< pressure [Pa].
  real,              intent(in) :: rho_ref  !< A reference density [kg m-3].

end function density_anomaly_elem_TEOS10
real elemental module function spec_vol_elem_TEOS10(this, T, S, pressure)
  class(TEOS10_EOS), intent(in) :: this     !< This EOS
  real,              intent(in) :: T        !< Conservative temperature [degC].
  real,              intent(in) :: S        !< Absolute salinity [g kg-1].
  real,              intent(in) :: pressure !< pressure [Pa].

end function spec_vol_elem_TEOS10
real elemental module function spec_vol_anomaly_elem_TEOS10(this, T, S, pressure, spv_ref)
  class(TEOS10_EOS), intent(in) :: this     !< This EOS
  real,              intent(in) :: T        !< Conservative temperature [degC].
  real,              intent(in) :: S        !< Absolute salinity [g kg-1].
  real,              intent(in) :: pressure !< pressure [Pa].
  real,              intent(in) :: spv_ref  !< A reference specific volume [m3 kg-1].

end function spec_vol_anomaly_elem_TEOS10
elemental module subroutine calculate_density_derivs_elem_TEOS10(this, T, S, pressure, drho_dT, drho_dS)
  class(TEOS10_EOS), intent(in)  :: this     !< This EOS
  real,              intent(in)  :: T        !< Conservative temperature [degC]
  real,              intent(in)  :: S        !< Absolute salinity [g kg-1] = [ppt]
  real,              intent(in)  :: pressure !< Pressure [Pa]
  real,              intent(out) :: drho_dT  !< The partial derivative of density with conservative
                                             !! temperature [kg m-3 degC-1]
  real,              intent(out) :: drho_dS  !< The partial derivative of density with salinity,
                                             !! in [kg m-3 ppt-1]
  ! Local variables

  ! Conversions
end subroutine calculate_density_derivs_elem_TEOS10
elemental module subroutine calculate_density_second_derivs_elem_TEOS10(this, T, S, pressure, &
                       drho_dS_dS, drho_dS_dT, drho_dT_dT, drho_dS_dP, drho_dT_dP)
  class(TEOS10_EOS), intent(in)    :: this !< This EOS
  real,              intent(in)    :: T        !< Conservative temperature [degC]
  real,              intent(in)    :: S        !< Absolute salinity [g kg-1] = [ppt]
  real,              intent(in)    :: pressure !< Pressure [Pa]
  real,              intent(inout) :: drho_ds_ds !< Partial derivative of beta with respect
                                                 !! to S [kg m-3 ppt-2]
  real,              intent(inout) :: drho_ds_dt !< Partial derivative of beta with respect
                                                 !! to T [kg m-3 ppt-1 degC-1]
  real,              intent(inout) :: drho_dt_dt !< Partial derivative of alpha with respect
                                                 !! to T [kg m-3 degC-2]
  real,              intent(inout) :: drho_ds_dp !< Partial derivative of beta with respect
                                                 !! to pressure [kg m-3 ppt-1 Pa-1] = [s2 m-2 ppt-1]
  real,              intent(inout) :: drho_dt_dp !< Partial derivative of alpha with respect
                                                 !! to pressure [kg m-3 degC-1 Pa-1] = [s2 m-2 degC-1]
  ! Local variables

  ! Conversions
end subroutine calculate_density_second_derivs_elem_TEOS10
elemental module subroutine calculate_specvol_derivs_elem_TEOS10(this, T, S, pressure, dSV_dT, dSV_dS)
  class(TEOS10_EOS),  intent(in)    :: this     !< This EOS
  real,               intent(in)    :: T        !< Conservative temperature [degC]
  real,               intent(in)    :: S        !< Absolute salinity [g kg-1] = [ppt]
  real,               intent(in)    :: pressure !< Pressure [Pa]
  real,               intent(inout) :: dSV_dT   !< The partial derivative of specific volume with
                                                !! conservative temperature [m3 kg-1 degC-1]
  real,               intent(inout) :: dSV_dS   !< The partial derivative of specific volume with
                                                !! absolute salinity [m3 kg-1 ppt-1]
  ! Local variables

  ! Conversions
end subroutine calculate_specvol_derivs_elem_TEOS10
elemental module subroutine calculate_compress_elem_TEOS10(this, T, S, pressure, rho, drho_dp)
  class(TEOS10_EOS),  intent(in)  :: this     !< This EOS
  real,               intent(in)  :: T        !< Conservative temperature [degC]
  real,               intent(in)  :: S        !< Absolute salinity [g kg-1]
  real,               intent(in)  :: pressure !< Pressure [Pa]
  real,               intent(out) :: rho      !< In situ density [kg m-3]
  real,               intent(out) :: drho_dp  !< The partial derivative of density with pressure
                                              !! (also the inverse of the square of sound speed)
                                              !! [s2 m-2]

  ! Local variables

  ! Conversions
end subroutine calculate_compress_elem_TEOS10
module subroutine EoS_fit_range_teos10(this, T_min, T_max, S_min, S_max, p_min, p_max)
  class(TEOS10_EOS),  intent(in)  :: this     !< This EOS
  real, optional, intent(out) :: T_min !< The minimum conservative temperature over which this EoS is fitted [degC]
  real, optional, intent(out) :: T_max !< The maximum conservative temperature over which this EoS is fitted [degC]
  real, optional, intent(out) :: S_min !< The minimum absolute salinity over which this EoS is fitted [g kg-1]
  real, optional, intent(out) :: S_max !< The maximum absolute salinity over which this EoS is fitted [g kg-1]
  real, optional, intent(out) :: p_min !< The minimum pressure over which this EoS is fitted [Pa]
  real, optional, intent(out) :: p_max !< The maximum pressure over which this EoS is fitted [Pa]

end subroutine EoS_fit_range_teos10
  end interface

end module MOM_EOS_TEOS10
