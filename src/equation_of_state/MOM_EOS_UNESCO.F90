! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> The equation of state using the Jackett and McDougall fits to the UNESCO EOS
module MOM_EOS_UNESCO

use MOM_EOS_base_type, only : EOS_base

implicit none ; private

public UNESCO_EOS

!>@{ Parameters in the UNESCO equation of state, as published in appendix A3 of Gill, 1982.
! The following constants are used to calculate rho0, the density of seawater at 1 atmosphere pressure.
! The notation is Rab for the contribution to rho0 from S^a*T^b, with 6 used for the 1.5 power.
real, parameter :: R00 = 999.842594   ! A coefficient in the fit for rho0 [kg m-3]
real, parameter :: R01 = 6.793952e-2  ! A coefficient in the fit for rho0 [kg m-3 degC-1]
real, parameter :: R02 = -9.095290e-3 ! A coefficient in the fit for rho0 [kg m-3 degC-2]
real, parameter :: R03 = 1.001685e-4  ! A coefficient in the fit for rho0 [kg m-3 degC-3]
real, parameter :: R04 = -1.120083e-6 ! A coefficient in the fit for rho0 [kg m-3 degC-4]
real, parameter :: R05 = 6.536332e-9  ! A coefficient in the fit for rho0 [kg m-3 degC-5]
real, parameter :: R10 = 0.824493     ! A coefficient in the fit for rho0 [kg m-3 PSU-1]
real, parameter :: R11 = -4.0899e-3   ! A coefficient in the fit for rho0 [kg m-3 degC-1 PSU-1]
real, parameter :: R12 = 7.6438e-5    ! A coefficient in the fit for rho0 [kg m-3 degC-2 PSU-1]
real, parameter :: R13 = -8.2467e-7   ! A coefficient in the fit for rho0 [kg m-3 degC-3 PSU-1]
real, parameter :: R14 = 5.3875e-9    ! A coefficient in the fit for rho0 [kg m-3 degC-4 PSU-1]
real, parameter :: R60 = -5.72466e-3  ! A coefficient in the fit for rho0 [kg m-3 PSU-1.5]
real, parameter :: R61 = 1.0227e-4    ! A coefficient in the fit for rho0 [kg m-3 degC-1 PSU-1.5]
real, parameter :: R62 = -1.6546e-6   ! A coefficient in the fit for rho0 [kg m-3 degC-2 PSU-1.5]
real, parameter :: R20 = 4.8314e-4    ! A coefficient in the fit for rho0 [kg m-3 PSU-2]

! The following constants are used to calculate the secant bulk modulus.
! The notation here is Sabc for terms proportional to S^a*T^b*P^c, with 6 used for the 1.5 power.
!   Note that these values differ from those in Appendix 3 of Gill (1982) because the expressions
! from Jackett and MacDougall (1995) use potential temperature, rather than in situ temperature.
real, parameter :: S000 = 1.965933e4   ! A coefficient in the secant bulk modulus fit [bar]
real, parameter :: S010 = 1.444304e2   ! A coefficient in the secant bulk modulus fit [bar degC-1]
real, parameter :: S020 = -1.706103    ! A coefficient in the secant bulk modulus fit [bar degC-2]
real, parameter :: S030 = 9.648704e-3  ! A coefficient in the secant bulk modulus fit [bar degC-3]
real, parameter :: S040 = -4.190253e-5 ! A coefficient in the secant bulk modulus fit [bar degC-4]
real, parameter :: S100 = 52.84855     ! A coefficient in the secant bulk modulus fit [bar PSU-1]
real, parameter :: S110 = -3.101089e-1 ! A coefficient in the secant bulk modulus fit [bar degC-1 PSU-1]
real, parameter :: S120 = 6.283263e-3  ! A coefficient in the secant bulk modulus fit [bar degC-2 PSU-1]
real, parameter :: S130 = -5.084188e-5 ! A coefficient in the secant bulk modulus fit [bar degC-3 PSU-1]
real, parameter :: S600 = 3.886640e-1  ! A coefficient in the secant bulk modulus fit [bar PSU-1.5]
real, parameter :: S610 = 9.085835e-3  ! A coefficient in the secant bulk modulus fit [bar degC-1 PSU-1.5]
real, parameter :: S620 = -4.619924e-4 ! A coefficient in the secant bulk modulus fit [bar degC-2 PSU-1.5]

real, parameter :: S001 = 3.186519     ! A coefficient in the secant bulk modulus fit [nondim]
real, parameter :: S011 = 2.212276e-2  ! A coefficient in the secant bulk modulus fit [degC-1]
real, parameter :: S021 = -2.984642e-4 ! A coefficient in the secant bulk modulus fit [degC-2]
real, parameter :: S031 = 1.956415e-6  ! A coefficient in the secant bulk modulus fit [degC-3]
real, parameter :: S101 = 6.704388e-3  ! A coefficient in the secant bulk modulus fit [PSU-1]
real, parameter :: S111 = -1.847318e-4 ! A coefficient in the secant bulk modulus fit [degC-1 PSU-1]
real, parameter :: S121 = 2.059331e-7  ! A coefficient in the secant bulk modulus fit [degC-2 PSU-1]
real, parameter :: S601 = 1.480266e-4  ! A coefficient in the secant bulk modulus fit [PSU-1.5]

real, parameter :: S002 = 2.102898e-4  ! A coefficient in the secant bulk modulus fit [bar-1]
real, parameter :: S012 = -1.202016e-5 ! A coefficient in the secant bulk modulus fit [bar-1 degC-1]
real, parameter :: S022 = 1.394680e-7  ! A coefficient in the secant bulk modulus fit [bar-1 degC-2]
real, parameter :: S102 = -2.040237e-6 ! A coefficient in the secant bulk modulus fit [bar-1 PSU-1]
real, parameter :: S112 = 6.128773e-8  ! A coefficient in the secant bulk modulus fit [bar-1 degC-1 PSU-1]
real, parameter :: S122 = 6.207323e-10 ! A coefficient in the secant bulk modulus fit [bar-1 degC-2 PSU-1]
!>@}

!> The EOS_base implementation of the UNESCO equation of state
type, extends (EOS_base) :: UNESCO_EOS

contains
  !> Implementation of the in-situ density as an elemental function [kg m-3]
  procedure :: density_elem => density_elem_UNESCO
  !> Implementation of the in-situ density anomaly as an elemental function [kg m-3]
  procedure :: density_anomaly_elem => density_anomaly_elem_UNESCO
  !> Implementation of the in-situ specific volume as an elemental function [m3 kg-1]
  procedure :: spec_vol_elem => spec_vol_elem_UNESCO
  !> Implementation of the in-situ specific volume anomaly as an elemental function [m3 kg-1]
  procedure :: spec_vol_anomaly_elem => spec_vol_anomaly_elem_UNESCO
  !> Implementation of the calculation of derivatives of density
  procedure :: calculate_density_derivs_elem => calculate_density_derivs_elem_UNESCO
  !> Implementation of the calculation of second derivatives of density
  procedure :: calculate_density_second_derivs_elem => calculate_density_second_derivs_elem_UNESCO
  !> Implementation of the calculation of derivatives of specific volume
  procedure :: calculate_specvol_derivs_elem => calculate_specvol_derivs_elem_UNESCO
  !> Implementation of the calculation of compressibility
  procedure :: calculate_compress_elem => calculate_compress_elem_UNESCO
  !> Implementation of the range query function
  procedure :: EOS_fit_range => EOS_fit_range_UNESCO

end type UNESCO_EOS


  interface
real elemental module function density_elem_UNESCO(this, T, S, pressure)
  class(UNESCO_EOS), intent(in) :: this     !< This EOS
  real,              intent(in) :: T        !< Potential temperature relative to the surface [degC]
  real,              intent(in) :: S        !< Salinity [PSU]
  real,              intent(in) :: pressure !< Pressure [Pa]

  ! Local variables

end function density_elem_UNESCO
real elemental module function density_anomaly_elem_UNESCO(this, T, S, pressure, rho_ref)
  class(UNESCO_EOS), intent(in) :: this     !< This EOS
  real,              intent(in) :: T        !< Potential temperature relative to the surface [degC]
  real,              intent(in) :: S        !< Salinity [PSU]
  real,              intent(in) :: pressure !< Pressure [Pa]
  real,              intent(in) :: rho_ref  !< A reference density [kg m-3]

  ! Local variables

end function density_anomaly_elem_UNESCO
real elemental module function spec_vol_elem_UNESCO(this, T, S, pressure)
  class(UNESCO_EOS), intent(in) :: this    !< This EOS
  real,           intent(in) :: T        !< Potential temperature relative to the surface [degC]
  real,           intent(in) :: S        !< Salinity [PSU]
  real,           intent(in) :: pressure !< Pressure [Pa]

  ! Local variables

end function spec_vol_elem_UNESCO
real elemental module function spec_vol_anomaly_elem_UNESCO(this, T, S, pressure, spv_ref)
  class(UNESCO_EOS), intent(in) :: this    !< This EOS
  real,           intent(in) :: T        !< Potential temperature relative to the surface [degC]
  real,           intent(in) :: S        !< Salinity [PSU]
  real,           intent(in) :: pressure !< Pressure [Pa]
  real,           intent(in) :: spv_ref  !< A reference specific volume [m3 kg-1]

  ! Local variables

end function spec_vol_anomaly_elem_UNESCO
elemental module subroutine calculate_density_derivs_elem_UNESCO(this, T, S, pressure, drho_dT, drho_dS)
  class(UNESCO_EOS), intent(in)  :: this     !< This EOS
  real,              intent(in)  :: T        !< Potential temperature relative to the surface [degC]
  real,              intent(in)  :: S        !< Salinity [PSU]
  real,              intent(in)  :: pressure !< Pressure [Pa]
  real,              intent(out) :: drho_dT  !< The partial derivative of density with potential
                                             !! temperature [kg m-3 degC-1]
  real,              intent(out) :: drho_dS  !< The partial derivative of density with salinity,
                                             !! in [kg m-3 PSU-1]
  ! Local variables

end subroutine calculate_density_derivs_elem_UNESCO
elemental module subroutine calculate_density_second_derivs_elem_UNESCO(this, T, S, pressure, &
                            drho_ds_ds, drho_ds_dt, drho_dt_dt, drho_ds_dp, drho_dt_dp)
  class(UNESCO_EOS), intent(in)    :: this !< This EOS
  real,              intent(in)    :: T !< Potential temperature referenced to 0 dbar [degC]
  real,              intent(in)    :: S !< Salinity [PSU]
  real,              intent(in)    :: pressure !< Pressure [Pa]
  real,              intent(inout) :: drho_ds_ds !< Partial derivative of beta with respect
                                                 !! to S [kg m-3 PSU-2]
  real,              intent(inout) :: drho_ds_dt !< Partial derivative of beta with respect
                                                 !! to T [kg m-3 PSU-1 degC-1]
  real,              intent(inout) :: drho_dt_dt !< Partial derivative of alpha with respect
                                                 !! to T [kg m-3 degC-2]
  real,              intent(inout) :: drho_ds_dp !< Partial derivative of beta with respect
                                                 !! to pressure [kg m-3 PSU-1 Pa-1] = [s2 m-2 PSU-1]
  real,              intent(inout) :: drho_dt_dp !< Partial derivative of alpha with respect
                                                 !! to pressure [kg m-3 degC-1 Pa-1] = [s2 m-2 degC-1]

  ! Local variables

end subroutine calculate_density_second_derivs_elem_UNESCO
elemental module subroutine calculate_specvol_derivs_elem_UNESCO(this, T, S, pressure, dSV_dT, dSV_dS)
  class(UNESCO_EOS), intent(in)    :: this     !< This EOS
  real,              intent(in)    :: T        !< Potential temperature [degC]
  real,              intent(in)    :: S        !< Salinity [PSU]
  real,              intent(in)    :: pressure !< Pressure [Pa]
  real,              intent(inout) :: dSV_dT   !< The partial derivative of specific volume with
                                               !! potential temperature [m3 kg-1 degC-1]
  real,              intent(inout) :: dSV_dS   !< The partial derivative of specific volume with
                                               !! salinity [m3 kg-1 PSU-1]
  ! Local variables

end subroutine calculate_specvol_derivs_elem_UNESCO
elemental module subroutine calculate_compress_elem_UNESCO(this, T, S, pressure, rho, drho_dp)
  class(UNESCO_EOS), intent(in)  :: this     !< This EOS
  real,              intent(in)  :: T        !< Potential temperature relative to the surface [degC]
  real,              intent(in)  :: S        !< Salinity [PSU]
  real,              intent(in)  :: pressure !< Pressure [Pa]
  real,              intent(out) :: rho      !< In situ density [kg m-3]
  real,              intent(out) :: drho_dp  !< The partial derivative of density with pressure
                                             !! (also the inverse of the square of sound speed)
                                             !! [s2 m-2]
  ! Local variables

end subroutine calculate_compress_elem_UNESCO
module subroutine EoS_fit_range_UNESCO(this, T_min, T_max, S_min, S_max, p_min, p_max)
  class(UNESCO_EOS), intent(in) :: this !< This EOS
  real, optional, intent(out) :: T_min !< The minimum potential temperature over which this EoS is fitted [degC]
  real, optional, intent(out) :: T_max !< The maximum potential temperature over which this EoS is fitted [degC]
  real, optional, intent(out) :: S_min !< The minimum practical salinity over which this EoS is fitted [PSU]
  real, optional, intent(out) :: S_max !< The maximum practical salinity over which this EoS is fitted [PSU]
  real, optional, intent(out) :: p_min !< The minimum pressure over which this EoS is fitted [Pa]
  real, optional, intent(out) :: p_max !< The maximum pressure over which this EoS is fitted [Pa]

end subroutine EoS_fit_range_UNESCO
  end interface

end module MOM_EOS_UNESCO
