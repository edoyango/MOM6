! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> The equation of state using the expressions of Roquet et al. (2015) that are used in NEMO
module MOM_EOS_Roquet_rho

use MOM_EOS_base_type, only : EOS_base

implicit none ; private

public Roquet_rho_EOS

real, parameter :: Pa2kb  = 1.e-8 !< Conversion factor between Pa and kbar [kbar Pa-1]
!>@{ Parameters in the Roquet_rho (Roquet density) equation of state
real, parameter :: rdeltaS = 32.          ! An offset to salinity before taking its square root [g kg-1]
real, parameter :: r1_S0 = 0.875/35.16504 ! The inverse of a plausible range of oceanic salinities [kg g-1]
real, parameter :: I_Ts = 0.025           ! The inverse of a plausible range of oceanic temperatures [degC-1]

! The following are the coefficients of the fit to the reference density profile (rho00p) as a function of
! pressure (P), with a contribution R0c * P**(c+1).  The nomenclature follows Roquet.
real, parameter :: R00 = 4.6494977072e+01*Pa2kb     ! rho00p P coef.    [kg m-3 Pa-1]
real, parameter :: R01 = -5.2099962525*Pa2kb**2     ! rho00p P**2 coef. [kg m-3 Pa-2]
real, parameter :: R02 = 2.2601900708e-01*Pa2kb**3  ! rho00p P**3 coef. [kg m-3 Pa-3]
real, parameter :: R03 = 6.4326772569e-02*Pa2kb**4  ! rho00p P**4 coef. [kg m-3 Pa-4]
real, parameter :: R04 = 1.5616995503e-02*Pa2kb**5  ! rho00p P**5 coef. [kg m-3 Pa-5]
real, parameter :: R05 = -1.7243708991e-03*Pa2kb**6 ! rho00p P**6 coef. [kg m-3 Pa-6]

! The following are coefficients of contributions to density as a function of the square root
! of normalized salinity with an offset (zs), temperature (T) and pressure (P), with a contribution
! EOSabc * zs**a * T**b * P**c.  The numbers here are copied directly from Roquet et al. (2015), but
! the expressions here do not use the same nondimensionalization for pressure or temperature as they do.
real, parameter :: EOS000 = 8.0189615746e+02                  ! A constant density contribution [kg m-3]
real, parameter :: EOS100 = 8.6672408165e+02                  ! EoS zs coef.                [kg m-3]
real, parameter :: EOS200 = -1.7864682637e+03                 ! EoS zs**2 coef.             [kg m-3]
real, parameter :: EOS300 = 2.0375295546e+03                  ! EoS zs**3 coef.             [kg m-3]
real, parameter :: EOS400 = -1.2849161071e+03                 ! EoS zs**4 coef.             [kg m-3]
real, parameter :: EOS500 = 4.3227585684e+02                  ! EoS zs**5 coef.             [kg m-3]
real, parameter :: EOS600 = -6.0579916612e+01                 ! EoS zs**6 coef.             [kg m-3]
real, parameter :: EOS010 = 2.6010145068e+01*I_Ts             ! EoS T coef.          [kg m-3 degC-1]
real, parameter :: EOS110 = -6.5281885265e+01*I_Ts            ! EoS zs * T coef.     [kg m-3 degC-1]
real, parameter :: EOS210 = 8.1770425108e+01*I_Ts             ! EoS zs**2 * T coef.  [kg m-3 degC-1]
real, parameter :: EOS310 = -5.6888046321e+01*I_Ts            ! EoS zs**3 * T coef.  [kg m-3 degC-1]
real, parameter :: EOS410 = 1.7681814114e+01*I_Ts             ! EoS zs**2 * T coef.  [kg m-3 degC-1]
real, parameter :: EOS510 = -1.9193502195*I_Ts                ! EoS zs**5 * T coef.  [kg m-3 degC-1]
real, parameter :: EOS020 = -3.7074170417e+01*I_Ts**2         ! EoS T**2 coef.       [kg m-3 degC-2]
real, parameter :: EOS120 = 6.1548258127e+01*I_Ts**2          ! EoS zs * T**2 coef.  [kg m-3 degC-2]
real, parameter :: EOS220 = -6.0362551501e+01*I_Ts**2         ! EoS zs**2 * T**2 coef. [kg m-3 degC-2]
real, parameter :: EOS320 = 2.9130021253e+01*I_Ts**2          ! EoS zs**3 * T**2 coef. [kg m-3 degC-2]
real, parameter :: EOS420 = -5.4723692739*I_Ts**2             ! EoS zs**4 * T**2 coef. [kg m-3 degC-2]
real, parameter :: EOS030 = 2.1661789529e+01*I_Ts**3          ! EoS T**3 coef.       [kg m-3 degC-3]
real, parameter :: EOS130 = -3.3449108469e+01*I_Ts**3         ! EoS zs * T**3 coef.  [kg m-3 degC-3]
real, parameter :: EOS230 = 1.9717078466e+01*I_Ts**3          ! EoS zs**2 * T**3 coef. [kg m-3 degC-3]
real, parameter :: EOS330 = -3.1742946532*I_Ts**3             ! EoS zs**3 * T**3 coef. [kg m-3 degC-3]
real, parameter :: EOS040 = -8.3627885467*I_Ts**4             ! EoS T**4 coef.       [kg m-3 degC-4]
real, parameter :: EOS140 = 1.1311538584e+01*I_Ts**4          ! EoS zs * T**4 coef.  [kg m-3 degC-4]
real, parameter :: EOS240 = -5.3563304045*I_Ts**4             ! EoS zs**2 * T**4 coef. [kg m-3 degC-4]
real, parameter :: EOS050 = 5.4048723791e-01*I_Ts**5          ! EoS T**5 coef.       [kg m-3 degC-5]
real, parameter :: EOS150 = 4.8169980163e-01*I_Ts**5          ! EoS zs * T**5 coef.  [kg m-3 degC-5]
real, parameter :: EOS060 = -1.9083568888e-01*I_Ts**6         ! EoS T**6             [kg m-3 degC-6]
real, parameter :: EOS001 = 1.9681925209e+01*Pa2kb            ! EoS P coef.            [kg m-3 Pa-1]
real, parameter :: EOS101 = -4.2549998214e+01*Pa2kb           ! EoS zs * P coef.       [kg m-3 Pa-1]
real, parameter :: EOS201 = 5.0774768218e+01*Pa2kb            ! EoS zs**2 * P coef.    [kg m-3 Pa-1]
real, parameter :: EOS301 = -3.0938076334e+01*Pa2kb           ! EoS zs**3 * P coef.    [kg m-3 Pa-1]
real, parameter :: EOS401 = 6.6051753097*Pa2kb                ! EoS zs**4 * P coef.    [kg m-3 Pa-1]
real, parameter :: EOS011 = -1.3336301113e+01*(I_Ts*Pa2kb)    ! EoS T * P coef. [kg m-3 degC-1 Pa-1]
real, parameter :: EOS111 = -4.4870114575*(I_Ts*Pa2kb)        ! EoS zs * T * P coef. [kg m-3 degC-1 Pa-1]
real, parameter :: EOS211 = 5.0042598061*(I_Ts*Pa2kb)         ! EoS zs**2 * T * P coef. [kg m-3 degC-1 Pa-1]
real, parameter :: EOS311 = -6.5399043664e-01*(I_Ts*Pa2kb)    ! EoS zs**3 * T * P coef. [kg m-3 degC-1 Pa-1]
real, parameter :: EOS021 = 6.7080479603*(I_Ts**2*Pa2kb)      ! EoS T**2 * P coef. [kg m-3 degC-2 Pa-1]
real, parameter :: EOS121 = 3.5063081279*(I_Ts**2*Pa2kb)      ! EoS zs * T**2 * P coef. [kg m-3 degC-2 Pa-1]
real, parameter :: EOS221 = -1.8795372996*(I_Ts**2*Pa2kb)     ! EoS zs**2 * T**2 * P coef. [kg m-3 degC-2 Pa-1]
real, parameter :: EOS031 = -2.4649669534*(I_Ts**3*Pa2kb)     ! EoS T**3 * P coef. [kg m-3 degC-3 Pa-1]
real, parameter :: EOS131 = -5.5077101279e-01*(I_Ts**3*Pa2kb) ! EoS zs * T**3 * P coef. [kg m-3 degC-3 Pa-1]
real, parameter :: EOS041 = 5.5927935970e-01*(I_Ts**4*Pa2kb)  ! EoS T**4 * P coef. [kg m-3 degC-4 Pa-1]
real, parameter :: EOS002 = 2.0660924175*Pa2kb**2             ! EoS P**2 coef.         [kg m-3 Pa-2]
real, parameter :: EOS102 = -4.9527603989*Pa2kb**2            ! EoS zs * P**2 coef.    [kg m-3 Pa-2]
real, parameter :: EOS202 = 2.5019633244*Pa2kb**2             ! EoS zs**2 * P**2 coef. [kg m-3 Pa-2]
real, parameter :: EOS012 = 2.0564311499*(I_Ts*Pa2kb**2)      ! EoS T * P**2 coef. [kg m-3 degC-1 Pa-2]
real, parameter :: EOS112 = -2.1311365518e-01*(I_Ts*Pa2kb**2) ! EoS zs * T * P**2 coef. [kg m-3 degC-1 Pa-2]
real, parameter :: EOS022 = -1.2419983026*(I_Ts**2*Pa2kb**2)  ! EoS T**2 * P**2 coef. [kg m-3 degC-2 Pa-2]
real, parameter :: EOS003 = -2.3342758797e-02*Pa2kb**3        ! EoS P**3 coef.         [kg m-3 Pa-3]
real, parameter :: EOS103 = -1.8507636718e-02*Pa2kb**3        ! EoS zs * P**3 coef.    [kg m-3 Pa-3]
real, parameter :: EOS013 = 3.7969820455e-01*(I_Ts*Pa2kb**3)  ! EoS T * P**3 coef. [kg m-3 degC-1 Pa-3]

real, parameter :: ALP000 =    EOS010   ! Constant in the drho_dT fit                [kg m-3 degC-1]
real, parameter :: ALP100 =    EOS110   ! drho_dT fit zs coef.                       [kg m-3 degC-1]
real, parameter :: ALP200 =    EOS210   ! drho_dT fit zs**2 coef.                    [kg m-3 degC-1]
real, parameter :: ALP300 =    EOS310   ! drho_dT fit zs**3 coef.                    [kg m-3 degC-1]
real, parameter :: ALP400 =    EOS410   ! drho_dT fit zs**4 coef.                    [kg m-3 degC-1]
real, parameter :: ALP500 =    EOS510   ! drho_dT fit zs**5 coef.                    [kg m-3 degC-1]
real, parameter :: ALP010 = 2.*EOS020   ! drho_dT fit T coef.                        [kg m-3 degC-2]
real, parameter :: ALP110 = 2.*EOS120   ! drho_dT fit zs * T coef.                   [kg m-3 degC-2]
real, parameter :: ALP210 = 2.*EOS220   ! drho_dT fit zs**2 * T coef.                [kg m-3 degC-2]
real, parameter :: ALP310 = 2.*EOS320   ! drho_dT fit zs**3 * T coef.                [kg m-3 degC-2]
real, parameter :: ALP410 = 2.*EOS420   ! drho_dT fit zs**4 * T coef.                [kg m-3 degC-2]
real, parameter :: ALP020 = 3.*EOS030   ! drho_dT fit T**2 coef.                     [kg m-3 degC-3]
real, parameter :: ALP120 = 3.*EOS130   ! drho_dT fit zs * T**2 coef.                [kg m-3 degC-3]
real, parameter :: ALP220 = 3.*EOS230   ! drho_dT fit zs**2 * T**2 coef.             [kg m-3 degC-3]
real, parameter :: ALP320 = 3.*EOS330   ! drho_dT fit zs**3 * T**2 coef.             [kg m-3 degC-3]
real, parameter :: ALP030 = 4.*EOS040   ! drho_dT fit T**3 coef.                     [kg m-3 degC-4]
real, parameter :: ALP130 = 4.*EOS140   ! drho_dT fit zs * T**3 coef.                [kg m-3 degC-4]
real, parameter :: ALP230 = 4.*EOS240   ! drho_dT fit zs**2 * T**3 coef.             [kg m-3 degC-4]
real, parameter :: ALP040 = 5.*EOS050   ! drho_dT fit T**4 coef.                     [kg m-3 degC-5]
real, parameter :: ALP140 = 5.*EOS150   ! drho_dT fit zs* * T**4 coef.               [kg m-3 degC-5]
real, parameter :: ALP050 = 6.*EOS060   ! drho_dT fit T**5 coef.                     [kg m-3 degC-6]
real, parameter :: ALP001 =    EOS011   ! drho_dT fit P coef.                   [kg m-3 degC-1 Pa-1]
real, parameter :: ALP101 =    EOS111   ! drho_dT fit zs * P coef.              [kg m-3 degC-1 Pa-1]
real, parameter :: ALP201 =    EOS211   ! drho_dT fit zs**2 * P coef.           [kg m-3 degC-1 Pa-1]
real, parameter :: ALP301 =    EOS311   ! drho_dT fit zs**3 * P coef.           [kg m-3 degC-1 Pa-1]
real, parameter :: ALP011 = 2.*EOS021   ! drho_dT fit T * P coef.               [kg m-3 degC-2 Pa-1]
real, parameter :: ALP111 = 2.*EOS121   ! drho_dT fit zs * T * P coef.          [kg m-3 degC-2 Pa-1]
real, parameter :: ALP211 = 2.*EOS221   ! drho_dT fit zs**2 * T * P coef.       [kg m-3 degC-2 Pa-1]
real, parameter :: ALP021 = 3.*EOS031   ! drho_dT fit T**2 * P coef.            [kg m-3 degC-3 Pa-1]
real, parameter :: ALP121 = 3.*EOS131   ! drho_dT fit zs * T**2 * P coef.       [kg m-3 degC-3 Pa-1]
real, parameter :: ALP031 = 4.*EOS041   ! drho_dT fit T**3 * P coef.            [kg m-3 degC-4 Pa-1]
real, parameter :: ALP002 =    EOS012   ! drho_dT fit P**2 coef.                [kg m-3 degC-1 Pa-2]
real, parameter :: ALP102 =    EOS112   ! drho_dT fit zs * P**2 coef.           [kg m-3 degC-1 Pa-2]
real, parameter :: ALP012 = 2.*EOS022   ! drho_dT fit T * P**2 coef.            [kg m-3 degC-2 Pa-2]
real, parameter :: ALP003 =    EOS013   ! drho_dT fit P**3 coef.                [kg m-3 degC-1 Pa-3]

real, parameter :: BET000 = 0.5*EOS100*r1_S0  ! Constant in the drho_dS fit           [kg m-3 ppt-1]
real, parameter :: BET100 =     EOS200*r1_S0  ! drho_dS fit zs coef.                  [kg m-3 ppt-1]
real, parameter :: BET200 = 1.5*EOS300*r1_S0  ! drho_dS fit zs**2 coef.               [kg m-3 ppt-1]
real, parameter :: BET300 = 2.0*EOS400*r1_S0  ! drho_dS fit zs**3 coef.               [kg m-3 ppt-1]
real, parameter :: BET400 = 2.5*EOS500*r1_S0  ! drho_dS fit zs**4 coef.               [kg m-3 ppt-1]
real, parameter :: BET500 = 3.0*EOS600*r1_S0  ! drho_dS fit zs**5 coef.               [kg m-3 ppt-1]
real, parameter :: BET010 = 0.5*EOS110*r1_S0  ! drho_dS fit T coef.            [kg m-3 ppt-1 degC-1]
real, parameter :: BET110 =     EOS210*r1_S0  ! drho_dS fit zs * T coef.       [kg m-3 ppt-1 degC-1]
real, parameter :: BET210 = 1.5*EOS310*r1_S0  ! drho_dS fit zs**2 * T coef.    [kg m-3 ppt-1 degC-1]
real, parameter :: BET310 = 2.0*EOS410*r1_S0  ! drho_dS fit zs**3 * T coef.    [kg m-3 ppt-1 degC-1]
real, parameter :: BET410 = 2.5*EOS510*r1_S0  ! drho_dS fit zs**4 * T coef.    [kg m-3 ppt-1 degC-1]
real, parameter :: BET020 = 0.5*EOS120*r1_S0  ! drho_dS fit T**2 coef.         [kg m-3 ppt-1 degC-2]
real, parameter :: BET120 =     EOS220*r1_S0  ! drho_dS fit zs * T**2 coef.    [kg m-3 ppt-1 degC-2]
real, parameter :: BET220 = 1.5*EOS320*r1_S0  ! drho_dS fit zs**2 * T**2 coef. [kg m-3 ppt-1 degC-2]
real, parameter :: BET320 = 2.0*EOS420*r1_S0  ! drho_dS fit zs**3 * T**2 coef. [kg m-3 ppt-1 degC-2]
real, parameter :: BET030 = 0.5*EOS130*r1_S0  ! drho_dS fit T**3 coef.         [kg m-3 ppt-1 degC-3]
real, parameter :: BET130 =     EOS230*r1_S0  ! drho_dS fit zs * T**3 coef.    [kg m-3 ppt-1 degC-3]
real, parameter :: BET230 = 1.5*EOS330*r1_S0  ! drho_dS fit zs**2 * T**3 coef. [kg m-3 ppt-1 degC-3]
real, parameter :: BET040 = 0.5*EOS140*r1_S0  ! drho_dS fit T**4 coef.         [kg m-3 ppt-1 degC-4]
real, parameter :: BET140 =     EOS240*r1_S0  ! drho_dS fit zs * T**4 coef.    [kg m-3 ppt-1 degC-4]
real, parameter :: BET050 = 0.5*EOS150*r1_S0  ! drho_dS fit T**5 coef.         [kg m-3 ppt-1 degC-5]
real, parameter :: BET001 = 0.5*EOS101*r1_S0  ! drho_dS fit P coef.              [kg m-3 ppt-1 Pa-1]
real, parameter :: BET101 =     EOS201*r1_S0  ! drho_dS fit zs * P coef.         [kg m-3 ppt-1 Pa-1]
real, parameter :: BET201 = 1.5*EOS301*r1_S0  ! drho_dS fit zs**2 * P coef.      [kg m-3 ppt-1 Pa-1]
real, parameter :: BET301 = 2.0*EOS401*r1_S0  ! drho_dS fit zs**3 * P coef.      [kg m-3 ppt-1 Pa-1]
real, parameter :: BET011 = 0.5*EOS111*r1_S0  ! drho_dS fit T * P coef.   [kg m-3 ppt-1 degC-1 Pa-1]
real, parameter :: BET111 =     EOS211*r1_S0  ! drho_dS fit zs * T * P coef. [kg m-3 ppt-1 degC-1 Pa-1]
real, parameter :: BET211 = 1.5*EOS311*r1_S0  ! drho_dS fit zs**2 * T * P coef. [kg m-3 ppt-1 degC-1 Pa-1]
real, parameter :: BET021 = 0.5*EOS121*r1_S0  ! drho_dS fit T**2 * P coef. [kg m-3 ppt-1 degC-2 Pa-1]
real, parameter :: BET121 =     EOS221*r1_S0  ! drho_dS fit zs * T**2 * P coef. [kg m-3 ppt-1 degC-2 Pa-1]
real, parameter :: BET031 = 0.5*EOS131*r1_S0  ! drho_dS fit T**3 * P coef. [kg m-3 ppt-1 degC-3 Pa-1]
real, parameter :: BET002 = 0.5*EOS102*r1_S0  ! drho_dS fit P**2 coef.           [kg m-3 ppt-1 Pa-2]
real, parameter :: BET102 =     EOS202*r1_S0  ! drho_dS fit zs * P**2 coef.      [kg m-3 ppt-1 Pa-2]
real, parameter :: BET012 = 0.5*EOS112*r1_S0  ! drho_dS fit T * P**2 coef. [kg m-3 ppt-1 degC-1 Pa-2]
real, parameter :: BET003 = 0.5*EOS103*r1_S0  ! drho_dS fit P**3 coef.           [kg m-3 ppt-1 Pa-3]
!>@}

!> The EOS_base implementation of the Roquet et al., 2015, equation of state
type, extends (EOS_base) :: Roquet_rho_EOS

contains
  !> Implementation of the in-situ density as an elemental function [kg m-3]
  procedure :: density_elem => density_elem_Roquet_rho
  !> Implementation of the in-situ density anomaly as an elemental function [kg m-3]
  procedure :: density_anomaly_elem => density_anomaly_elem_Roquet_rho
  !> Implementation of the in-situ specific volume as an elemental function [m3 kg-1]
  procedure :: spec_vol_elem => spec_vol_elem_Roquet_rho
  !> Implementation of the in-situ specific volume anomaly as an elemental function [m3 kg-1]
  procedure :: spec_vol_anomaly_elem => spec_vol_anomaly_elem_Roquet_rho
  !> Implementation of the calculation of derivatives of density
  procedure :: calculate_density_derivs_elem => calculate_density_derivs_elem_Roquet_rho
  !> Implementation of the calculation of second derivatives of density
  procedure :: calculate_density_second_derivs_elem => calculate_density_second_derivs_elem_Roquet_rho
  !> Implementation of the calculation of derivatives of specific volume
  procedure :: calculate_specvol_derivs_elem => calculate_specvol_derivs_elem_Roquet_rho
  !> Implementation of the calculation of compressibility
  procedure :: calculate_compress_elem => calculate_compress_elem_Roquet_rho
  !> Implementation of the range query function
  procedure :: EOS_fit_range => EOS_fit_range_Roquet_rho

  !> Local implementation of generic calculate_density_array for efficiency
  procedure :: calculate_density_array => calculate_density_array_Roquet_rho
  !> Local implementation of generic calculate_spec_vol_array for efficiency
  procedure :: calculate_spec_vol_array => calculate_spec_vol_array_Roquet_rho

end type Roquet_rho_EOS


  interface
real elemental module function density_elem_Roquet_rho(this, T, S, pressure)
  class(Roquet_rho_EOS), intent(in) :: this     !< This EOS
  real,                  intent(in) :: T        !< Conservative temperature [degC]
  real,                  intent(in) :: S        !< Absolute salinity [g kg-1]
  real,                  intent(in) :: pressure !< Pressure [Pa]

  ! Local variables
                 ! by an assumed salinity range [nondim]
                 ! density at the reference temperature and salinity [kg m-3]
                 ! surface pressure [kg m-3]

  ! The following algorithm was published by Roquet et al. (2015), intended for use with NEMO.

  ! Conversions to the units used here.
end function density_elem_Roquet_rho
real elemental module function density_anomaly_elem_Roquet_rho(this, T, S, pressure, rho_ref)
  class(Roquet_rho_EOS), intent(in) :: this     !< This EOS
  real,                  intent(in) :: T        !< Conservative temperature [degC]
  real,                  intent(in) :: S        !< Absolute salinity [g kg-1]
  real,                  intent(in) :: pressure !< Pressure [Pa]
  real,                  intent(in) :: rho_ref  !< A reference density [kg m-3]

  ! Local variables
                 ! by an assumed salinity range [nondim]
                 ! density at the reference temperature and salinity [kg m-3]
                 ! surface pressure [kg m-3]

  ! The following algorithm was published by Roquet et al. (2015), intended for use with NEMO.

  ! Conversions to the units used here.
end function density_anomaly_elem_Roquet_rho
real elemental module function spec_vol_elem_Roquet_rho(this, T, S, pressure)
  class(Roquet_rho_EOS), intent(in) :: this     !< This EOS
  real,                  intent(in) :: T        !< Conservative temperature [degC]
  real,                  intent(in) :: S        !< Absolute salinity [g kg-1]
  real,                  intent(in) :: pressure !< Pressure [Pa]

end function spec_vol_elem_Roquet_rho
real elemental module function spec_vol_anomaly_elem_Roquet_rho(this, T, S, pressure, spv_ref)
  class(Roquet_rho_EOS), intent(in) :: this     !< This EOS
  real,                  intent(in) :: T        !< Conservative temperature [degC]
  real,                  intent(in) :: S        !< Absolute salinity [g kg-1]
  real,                  intent(in) :: pressure !< Pressure [Pa]
  real,                  intent(in) :: spv_ref  !< A reference specific volume [m3 kg-1]

end function spec_vol_anomaly_elem_Roquet_rho
elemental module subroutine calculate_density_derivs_elem_Roquet_rho(this, T, S, pressure, drho_dT, drho_dS)
  class(Roquet_rho_EOS), intent(in)  :: this     !< This EOS
  real,                  intent(in)  :: T        !< Conservative temperature [degC]
  real,                  intent(in)  :: S        !< Absolute salinity [g kg-1]
  real,                  intent(in)  :: pressure !< Pressure [Pa]
  real,                  intent(out) :: drho_dT  !< The partial derivative of density with potential
                                                 !! temperature [kg m-3 degC-1]
  real,                  intent(out) :: drho_dS  !< The partial derivative of density with salinity,
                                                 !! in [kg m-3 ppt-1]

  ! Local variables
                  ! by an assumed salinity range [nondim]
                  ! from temperature anomalies at the surface pressure
                  ! proportional to pressure
                  ! proportional to pressure**2
                  ! proportional to pressure**3
                  ! salinity [kg m-3 ppt-1] from temperature anomalies at the surface pressure
                  ! salinity [kg m-3 ppt-1 Pa-1] proportional to pressure
                  ! salinity [kg m-3 ppt-1 Pa-2] proportional to pressure**2
                  ! salinity [kg m-3 ppt-1 Pa-3] proportional to pressure**3

  ! Conversions to the units used here.
end subroutine calculate_density_derivs_elem_Roquet_rho
elemental module subroutine calculate_density_second_derivs_elem_Roquet_rho(this, T, S, pressure, &
                       drho_ds_ds, drho_ds_dt, drho_dt_dt, drho_ds_dp, drho_dt_dp)
  class(Roquet_rho_EOS), intent(in) :: this !< This EOS
  real,               intent(in)    :: T !< Conservative temperature [degC]
  real,               intent(in)    :: S !< Absolute salinity [g kg-1]
  real,               intent(in)    :: pressure !< Pressure [Pa]
  real,               intent(inout) :: drho_ds_ds !< Partial derivative of beta with respect
                                                  !! to S [kg m-3 ppt-2]
  real,               intent(inout) :: drho_ds_dt !< Partial derivative of beta with respect
                                                  !! to T [kg m-3 ppt-1 degC-1]
  real,               intent(inout) :: drho_dt_dt !< Partial derivative of alpha with respect
                                                  !! to T [kg m-3 degC-2]
  real,               intent(inout) :: drho_ds_dp !< Partial derivative of beta with respect
                                                  !! to pressure [kg m-3 ppt-1 Pa-1] = [s2 m-2 ppt-1]
  real,               intent(inout) :: drho_dt_dp !< Partial derivative of alpha with respect
                                                  !! to pressure [kg m-3 degC-1 Pa-1] = [s2 m-2 degC-1]

  ! Local variables
                 ! by an assumed salinity range [nondim]

  ! Conversions to the units used here.
end subroutine calculate_density_second_derivs_elem_Roquet_rho
elemental module subroutine calculate_specvol_derivs_elem_Roquet_rho(this, T, S, pressure, dSV_dT, dSV_dS)
  class(Roquet_rho_EOS), intent(in)    :: this     !< This EOS
  real,                  intent(in)    :: T        !< Conservative temperature [degC]
  real,                  intent(in)    :: S        !< Absolute salinity [g kg-1]
  real,                  intent(in)    :: pressure !< Pressure [Pa]
  real,                  intent(inout) :: dSV_dT   !< The partial derivative of specific volume with
                                                   !! potential temperature [m3 kg-1 degC-1]
  real,                  intent(inout) :: dSV_dS   !< The partial derivative of specific volume with
                                                   !! salinity [m3 kg-1 ppt-1]
  ! Local variables

end subroutine calculate_specvol_derivs_elem_Roquet_rho
elemental module subroutine calculate_compress_elem_Roquet_rho(this, T, S, pressure, rho, drho_dp)
  class(Roquet_rho_EOS), intent(in)  :: this !< This EOS
  real,                  intent(in)  :: T        !< Conservative temperature [degC]
  real,                  intent(in)  :: S        !< Absolute salinity [g kg-1]
  real,                  intent(in)  :: pressure !< Pressure [Pa]
  real,                  intent(out) :: rho      !< In situ density [kg m-3]
  real,                  intent(out) :: drho_dp  !< The partial derivative of density with pressure
                                                 !! (also the inverse of the square of sound speed)
                                                 !! [s2 m-2]
  ! Local variables
                 ! by an assumed salinity range [nondim]
                 ! density profile [kg m-3]
                 ! surface pressure [kg m-3]

  ! The following algorithm was published by Roquet et al. (2015), intended for use with NEMO.
  ! Conversions to the units used here.
end subroutine calculate_compress_elem_Roquet_rho
module subroutine EoS_fit_range_Roquet_rho(this, T_min, T_max, S_min, S_max, p_min, p_max)
  class(Roquet_rho_EOS), intent(in) :: this !< This EOS
  real, optional, intent(out) :: T_min !< The minimum conservative temperature over which this EoS is fitted [degC]
  real, optional, intent(out) :: T_max !< The maximum conservative temperature over which this EoS is fitted [degC]
  real, optional, intent(out) :: S_min !< The minimum absolute salinity over which this EoS is fitted [g kg-1]
  real, optional, intent(out) :: S_max !< The maximum absolute salinity over which this EoS is fitted [g kg-1]
  real, optional, intent(out) :: p_min !< The minimum pressure over which this EoS is fitted [Pa]
  real, optional, intent(out) :: p_max !< The maximum pressure over which this EoS is fitted [Pa]

end subroutine EoS_fit_range_Roquet_rho
module subroutine calculate_density_array_Roquet_rho(this, T, S, pressure, rho, start, npts, rho_ref)
  class(Roquet_rho_EOS),  intent(in) :: this  !< This EOS
  real, dimension(:), intent(in)  :: T        !< Potential temperature relative to the surface [degC]
  real, dimension(:), intent(in)  :: S        !< Salinity [PSU]
  real, dimension(:), intent(in)  :: pressure !< Pressure [Pa]
  real, dimension(:), intent(out) :: rho      !< In situ density [kg m-3]
  integer,            intent(in)  :: start    !< The starting index for calculations
  integer,            intent(in)  :: npts     !< The number of values to calculate
  real,     optional, intent(in)  :: rho_ref  !< A reference density [kg m-3]

  ! Local variables

end subroutine calculate_density_array_Roquet_rho
module subroutine calculate_spec_vol_array_Roquet_rho(this, T, S, pressure, specvol, start, npts, spv_ref)
  class(Roquet_rho_EOS),  intent(in) :: this  !< This EOS
  real, dimension(:), intent(in)  :: T        !< Potential temperature relative to the surface [degC]
  real, dimension(:), intent(in)  :: S        !< Salinity [PSU]
  real, dimension(:), intent(in)  :: pressure !< Pressure [Pa]
  real, dimension(:), intent(out) :: specvol  !< In situ specific volume [m3 kg-1]
  integer,            intent(in)  :: start    !< The starting index for calculations
  integer,            intent(in)  :: npts     !< The number of values to calculate
  real,     optional, intent(in)  :: spv_ref  !< A reference specific volume [m3 kg-1]

  ! Local variables

end subroutine calculate_spec_vol_array_Roquet_rho
  end interface

end module MOM_EOS_Roquet_rho
