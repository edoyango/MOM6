! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> The equation of state using the Jackett et al 2006 expressions that are often used in Hycom
module MOM_EOS_Jackett06

use MOM_EOS_base_type, only : EOS_base

implicit none ; private

public Jackett06_EOS

!>@{ Parameters in the Jackett et al. equation of state, which is a fit to the Fiestel (2003)
!    equation of state for the range: -2 < theta < 40 [degC], 0 < S < 42 [PSU], 0 < p < 1e8 [Pa].
!    The notation here is for terms in the numerator of the expression for density of
!    RNabc for terms proportional to S**a * T**b * P**c, and terms in the denominator as RDabc.
!    For terms proportional to S**1.5, 6 is used in this notation.

! --- coefficients for 25-term rational function sigloc().
real, parameter :: &
  RN000 =  9.9984085444849347d+02, & ! Density numerator constant coefficient [kg m-3]
  RN001 =  1.1798263740430364d-06, & ! Density numerator P       coefficient [kg m-3 Pa-1]
  RN002 = -2.5862187075154352d-16, & ! Density numerator P^2     coefficient [kg m-3 Pa-2]
  RN010 =  7.3471625860981584d+00, & ! Density numerator T       coefficient [kg m-3 degC-1]
  RN020 = -5.3211231792841769d-02, & ! Density numerator T^2     coefficient [kg m-3 degC-2]
  RN021 =  9.8920219266399117d-12, & ! Density numerator T^2 P   coefficient [kg m-3 degC-2 Pa-1]
  RN022 = -3.2921414007960662d-20, & ! Density numerator T^2 P^2 coefficient [kg m-3 degC-2 Pa-2]
  RN030 =  3.6492439109814549d-04, & ! Density numerator T^3     coefficient [kg m-3 degC-3]
  RN100 =  2.5880571023991390d+00, & ! Density numerator S       coefficient [kg m-3 PSU-1]
  RN101 =  4.6996642771754730d-10, & ! Density numerator S P     coefficient [kg m-3 PSU-1 Pa-1]
  RN110 = -6.7168282786692355d-03, & ! Density numerator S T     coefficient [kg m-3 degC-1 PSU-1]
  RN200 =  1.9203202055760151d-03, & ! Density numerator S^2      coefficient [kg m-3]

  RD001 =  6.7103246285651894d-10, & ! Density denominator P       coefficient [Pa-1]
  RD010 =  7.2815210113327091d-03, & ! Density denominator T       coefficient [degC-1]
  RD013 = -9.1534417604289062d-30, & ! Density denominator T P^3   coefficient [degC-1 Pa-3]
  RD020 = -4.4787265461983921d-05, & ! Density denominator T^2     coefficient [degC-2]
  RD030 =  3.3851002965802430d-07, & ! Density denominator T^3     coefficient [degC-3]
  RD032 = -2.4461698007024582d-25, & ! Density denominator T^3 P^2 coefficient [degC-3  Pa-2]
  RD040 =  1.3651202389758572d-10, & ! Density denominator T^4     coefficient [degC-4]
  RD100 =  1.7632126669040377d-03, & ! Density denominator S       coefficient [PSU-1]
  RD110 = -8.8066583251206474d-06, & ! Density denominator S T     coefficient [degC-1 PSU-1]
  RD130 = -1.8832689434804897d-10, & ! Density denominator S T^3   coefficient [degC-3 PSU-1]
  RD600 =  5.7463776745432097d-06, & ! Density denominator S^1.5   coefficient [PSU-1.5]
  RD620 =  1.4716275472242334d-09    ! Density denominator S^1.5 T^2 coefficient [degC-2 PSU-1.5]
!>@}

!> The EOS_base implementation of the Jackett et al, 2006, equation of state
type, extends (EOS_base) :: Jackett06_EOS

contains
  !> Implementation of the in-situ density as an elemental function [kg m-3]
  procedure :: density_elem => density_elem_Jackett06
  !> Implementation of the in-situ density anomaly as an elemental function [kg m-3]
  procedure :: density_anomaly_elem => density_anomaly_elem_Jackett06
  !> Implementation of the in-situ specific volume as an elemental function [m3 kg-1]
  procedure :: spec_vol_elem => spec_vol_elem_Jackett06
  !> Implementation of the in-situ specific volume anomaly as an elemental function [m3 kg-1]
  procedure :: spec_vol_anomaly_elem => spec_vol_anomaly_elem_Jackett06
  !> Implementation of the calculation of derivatives of density
  procedure :: calculate_density_derivs_elem => calculate_density_derivs_elem_Jackett06
  !> Implementation of the calculation of second derivatives of density
  procedure :: calculate_density_second_derivs_elem => calculate_density_second_derivs_elem_Jackett06
  !> Implementation of the calculation of derivatives of specific volume
  procedure :: calculate_specvol_derivs_elem => calculate_specvol_derivs_elem_Jackett06
  !> Implementation of the calculation of compressibility
  procedure :: calculate_compress_elem => calculate_compress_elem_Jackett06
  !> Implementation of the range query function
  procedure :: EOS_fit_range => EOS_fit_range_Jackett06

end type Jackett06_EOS


  interface
real elemental module function density_elem_Jackett06(this, T, S, pressure)
  class(Jackett06_EOS), intent(in) :: this !< This EOS
  real, intent(in) :: T        !< Potential temperature relative to the surface [degC].
  real, intent(in) :: S        !< Salinity [PSU].
  real, intent(in) :: pressure !< Pressure [Pa].

  ! Local variables
                  ! for density [kg m-3]
                  ! for density [nondim]

end function density_elem_Jackett06
real elemental module function density_anomaly_elem_Jackett06(this, T, S, pressure, rho_ref)
  class(Jackett06_EOS), intent(in) :: this !< This EOS
  real, intent(in) :: T        !< Potential temperature relative to the surface [degC].
  real, intent(in) :: S        !< Salinity [PSU].
  real, intent(in) :: pressure !< Pressure [Pa].
  real, intent(in) :: rho_ref  !< A reference density [kg m-3].

  ! Local variables
                  ! for density [kg m-3]
                  ! for density [nondim]

end function density_anomaly_elem_Jackett06
real elemental module function spec_vol_elem_Jackett06(this, T, S, pressure)
  class(Jackett06_EOS), intent(in) :: this !< This EOS
  real,           intent(in) :: T        !< potential temperature relative to the surface [degC].
  real,           intent(in) :: S        !< salinity [PSU].
  real,           intent(in) :: pressure !< pressure [Pa].

  ! Local variables
                  ! for density (not specific volume) [kg m-3]
                  ! for density (not specific volume) [nondim]

end function spec_vol_elem_Jackett06
real elemental module function spec_vol_anomaly_elem_Jackett06(this, T, S, pressure, spv_ref)
  class(Jackett06_EOS), intent(in) :: this !< This EOS
  real,           intent(in) :: T        !< potential temperature relative to the surface [degC].
  real,           intent(in) :: S        !< salinity [PSU].
  real,           intent(in) :: pressure !< pressure [Pa].
  real,           intent(in) :: spv_ref  !< A reference specific volume [m3 kg-1].

  ! Local variables
                  ! for density (not specific volume) [kg m-3]
                  ! for density (not specific volume) [nondim]

end function spec_vol_anomaly_elem_Jackett06
elemental module subroutine calculate_density_derivs_elem_Jackett06(this, T, S, pressure, drho_dT, drho_dS)
  class(Jackett06_EOS), intent(in) :: this    !< This EOS
  real,                 intent(in)  :: T        !< Potential temperature relative to the surface [degC]
  real,                 intent(in)  :: S        !< Salinity [PSU]
  real,                 intent(in)  :: pressure !< Pressure [Pa]
  real,                 intent(out) :: drho_dT  !< The partial derivative of density with potential
                                                !! temperature [kg m-3 degC-1]
  real,                 intent(out) :: drho_dS  !< The partial derivative of density with salinity,
                                                !! in [kg m-3 PSU-1]
  ! Local variables
                  ! for density [nondim]

end subroutine calculate_density_derivs_elem_Jackett06
elemental module subroutine calculate_density_second_derivs_elem_Jackett06(this, T, S, pressure, &
                       drho_ds_ds, drho_ds_dt, drho_dt_dt, drho_ds_dp, drho_dt_dp)
  class(Jackett06_EOS), intent(in)  :: this     !< This EOS
  real,               intent(in)    :: T !< Potential temperature referenced to 0 dbar [degC]
  real,               intent(in)    :: S !< Salinity [PSU]
  real,               intent(in)    :: pressure !< Pressure [Pa]
  real,               intent(inout) :: drho_ds_ds !< Partial derivative of beta with respect
                                                  !! to S [kg m-3 PSU-2]
  real,               intent(inout) :: drho_ds_dt !< Partial derivative of beta with respect
                                                  !! to T [kg m-3 PSU-1 degC-1]
  real,               intent(inout) :: drho_dt_dt !< Partial derivative of alpha with respect
                                                  !! to T [kg m-3 degC-2]
  real,               intent(inout) :: drho_ds_dp !< Partial derivative of beta with respect
                                                  !! to pressure [kg m-3 PSU-1 Pa-1] = [s2 m-2 PSU-1]
  real,               intent(inout) :: drho_dt_dp !< Partial derivative of alpha with respect
                                                  !! to pressure [kg m-3 degC-1 Pa-1] = [s2 m-2 degC-1]

  ! Local variables
                      ! for density [nondim]
                      ! salinity [kg m-3 degC-1 PSU-1]
                      ! pressure [kg m-3 degC-1 dbar-1]
                      ! pressure [kg m-3 PSU-1 dbar-1]
                      ! for density [nondim]
                      ! for density [nondim]

end subroutine calculate_density_second_derivs_elem_Jackett06
elemental module subroutine calculate_specvol_derivs_elem_Jackett06(this, T, S, pressure, dSV_dT, dSV_dS)
  class(Jackett06_EOS), intent(in)  :: this !< This EOS
  real,               intent(in)    :: T        !< Potential temperature [degC]
  real,               intent(in)    :: S        !< Salinity [PSU]
  real,               intent(in)    :: pressure !< Pressure [Pa]
  real,               intent(inout) :: dSV_dT   !< The partial derivative of specific volume with
                                                !! potential temperature [m3 kg-1 degC-1]
  real,               intent(inout) :: dSV_dS   !< The partial derivative of specific volume with
                                                !! salinity [m3 kg-1 PSU-1]

  ! Local variables
                  ! for density [nondim]

end subroutine calculate_specvol_derivs_elem_Jackett06
elemental module subroutine calculate_compress_elem_Jackett06(this, T, S, pressure, rho, drho_dp)
  class(Jackett06_EOS), intent(in) :: this    !< This EOS
  real,               intent(in)  :: T        !< Potential temperature relative to the surface [degC]
  real,               intent(in)  :: S        !< Salinity [PSU]
  real,               intent(in)  :: pressure !< Pressure [Pa]
  real,               intent(out) :: rho      !< In situ density [kg m-3]
  real,               intent(out) :: drho_dp  !< The partial derivative of density with pressure
                                              !! (also the inverse of the square of sound speed)
                                              !! [s2 m-2]
  ! Local variables

end subroutine calculate_compress_elem_Jackett06
module subroutine EoS_fit_range_Jackett06(this, T_min, T_max, S_min, S_max, p_min, p_max)
  class(Jackett06_EOS), intent(in) :: this !< This EOS
  real, optional, intent(out) :: T_min !< The minimum potential temperature over which this EoS is fitted [degC]
  real, optional, intent(out) :: T_max !< The maximum potential temperature over which this EoS is fitted [degC]
  real, optional, intent(out) :: S_min !< The minimum practical salinity over which this EoS is fitted [PSU]
  real, optional, intent(out) :: S_max !< The maximum practical salinity over which this EoS is fitted [PSU]
  real, optional, intent(out) :: p_min !< The minimum pressure over which this EoS is fitted [Pa]
  real, optional, intent(out) :: p_max !< The maximum pressure over which this EoS is fitted [Pa]

  ! Note that the actual fit range is given for the surface range of temperatures and salinities,
  ! but Jackett et al. use a more limited range of properties at higher pressures.
end subroutine EoS_fit_range_Jackett06
  end interface

end module MOM_EOS_Jackett06
