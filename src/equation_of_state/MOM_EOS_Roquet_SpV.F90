! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> The equation of state for specific volume (SpV) using the expressions of Roquet et al. 2015
module MOM_EOS_Roquet_Spv

use MOM_EOS_base_type, only : EOS_base

implicit none ; private

public Roquet_SpV_EOS

real, parameter :: Pa2kb  = 1.e-8 !< Conversion factor between Pa and kbar [kbar Pa-1]
!>@{ Parameters in the Roquet specific volume polynomial equation of state
real, parameter :: rdeltaS = 24.          ! An offset to salinity before taking its square root [g kg-1]
real, parameter :: r1_S0 = 0.875/35.16504 ! The inverse of a plausible range of oceanic salinities [kg g-1]
real, parameter :: I_Ts = 0.025           ! The inverse of a plausible range of oceanic temperatures [degC-1]
! The following are the coefficients of the fit to the reference density profile (rho00p) as a function of
! pressure (P), with a contribution R0c * P**(c+1).  The nomenclature follows Roquet.
real, parameter :: V00 = -4.4015007269e-05*Pa2kb    ! SpV00p P coef.    [m3 kg-1 Pa-1]
real, parameter :: V01 = 6.9232335784e-06*Pa2kb**2  ! SpV00p P**2 coef. [m3 kg-1 Pa-2]
real, parameter :: V02 = -7.5004675975e-07*Pa2kb**3 ! SpV00p P**3 coef. [m3 kg-1 Pa-3]
real, parameter :: V03 = 1.7009109288e-08*Pa2kb**4  ! SpV00p P**4 coef. [m3 kg-1 Pa-4]
real, parameter :: V04 = -1.6884162004e-08*Pa2kb**5 ! SpV00p P**5 coef. [m3 kg-1 Pa-5]
real, parameter :: V05 = 1.9613503930e-09*Pa2kb**6  ! SpV00p P**6 coef. [m3 kg-1 Pa-6]

! The following terms are contributions to specific volume (SpV) as a function of the square root of
! normalized absolute salinity with an offset (zs), temperature (T) and pressure (P), with a contribution
! SPVabc * zs**a * T**b * P**c.  The numbers here are copied directly from Roquet et al. (2015), but
! the expressions here do not use the same nondimensionalization for pressure or temperature as they do.
real, parameter :: SPV000 = 1.0772899069e-03                  ! Constant SpV contribution  [m3 kg-1]
real, parameter :: SPV100 = -3.1263658781e-04                 ! SpV zs coef.               [m3 kg-1]
real, parameter :: SPV200 = 6.7615860683e-04                  ! SpV zs**2 coef.            [m3 kg-1]
real, parameter :: SPV300 = -8.6127884515e-04                 ! SpV zs**3 coef.            [m3 kg-1]
real, parameter :: SPV400 = 5.9010812596e-04                  ! SpV zs**4 coef.            [m3 kg-1]
real, parameter :: SPV500 = -2.1503943538e-04                 ! SpV zs**5 coef.            [m3 kg-1]
real, parameter :: SPV600 = 3.2678954455e-05                  ! SpV zs**6 coef.            [m3 kg-1]
real, parameter :: SPV010 = -1.4949652640e-05*I_Ts            ! SpV T coef.         [m3 kg-1 degC-1]
real, parameter :: SPV110 = 3.1866349188e-05*I_Ts             ! SpV zs * T coef.    [m3 kg-1 degC-1]
real, parameter :: SPV210 = -3.8070687610e-05*I_Ts            ! SpV zs**2 * T coef. [m3 kg-1 degC-1]
real, parameter :: SPV310 = 2.9818473563e-05*I_Ts             ! SpV zs**3 * T coef. [m3 kg-1 degC-1]
real, parameter :: SPV410 = -1.0011321965e-05*I_Ts            ! SpV zs**4 * T coef. [m3 kg-1 degC-1]
real, parameter :: SPV510 = 1.0751931163e-06*I_Ts             ! SpV zs**5 * T coef. [m3 kg-1 degC-1]
real, parameter :: SPV020 = 2.7546851539e-05*I_Ts**2          ! SpV T**2 coef.      [m3 kg-1 degC-2]
real, parameter :: SPV120 = -3.6597334199e-05*I_Ts**2         ! SpV zs * T**2 coef. [m3 kg-1 degC-2]
real, parameter :: SPV220 = 3.4489154625e-05*I_Ts**2          ! SpV zs**2 * T**2 coef. [m3 kg-1 degC-2]
real, parameter :: SPV320 = -1.7663254122e-05*I_Ts**2         ! SpV zs**3 * T**2 coef. [m3 kg-1 degC-2]
real, parameter :: SPV420 = 3.5965131935e-06*I_Ts**2          ! SpV zs**4 * T**2 coef. [m3 kg-1 degC-2]
real, parameter :: SPV030 = -1.6506828994e-05*I_Ts**3         ! SpV T**3 coef.      [m3 kg-1 degC-3]
real, parameter :: SPV130 = 2.4412359055e-05*I_Ts**3          ! SpV zs * T**3 coef. [m3 kg-1 degC-3]
real, parameter :: SPV230 = -1.4606740723e-05*I_Ts**3         ! SpV zs**2 * T**3 coef. [m3 kg-1 degC-3]
real, parameter :: SPV330 = 2.3293406656e-06*I_Ts**3          ! SpV zs**3 * T**3 coef. [m3 kg-1 degC-3]
real, parameter :: SPV040 = 6.7896174634e-06*I_Ts**4          ! SpV T**4 coef.      [m3 kg-1 degC-4]
real, parameter :: SPV140 = -8.7951832993e-06*I_Ts**4         ! SpV zs * T**4 coef. [m3 kg-1 degC-4]
real, parameter :: SPV240 = 4.4249040774e-06*I_Ts**4          ! SpV zs**2 * T**4 coef. [m3 kg-1 degC-4]
real, parameter :: SPV050 = -7.2535743349e-07*I_Ts**5         ! SpV T**5 coef.      [m3 kg-1 degC-5]
real, parameter :: SPV150 = -3.4680559205e-07*I_Ts**5         ! SpV zs * T**5 coef. [m3 kg-1 degC-5]
real, parameter :: SPV060 = 1.9041365570e-07*I_Ts**6          ! SpV T**6 coef.      [m3 kg-1 degC-6]
real, parameter :: SPV001 = -1.6889436589e-05*Pa2kb           ! SpV P coef.           [m3 kg-1 Pa-1]
real, parameter :: SPV101 = 2.1106556158e-05*Pa2kb            ! SpV zs * P coef.      [m3 kg-1 Pa-1]
real, parameter :: SPV201 = -2.1322804368e-05*Pa2kb           ! SpV zs**2 * P coef.   [m3 kg-1 Pa-1]
real, parameter :: SPV301 = 1.7347655458e-05*Pa2kb            ! SpV zs**3 * P coef.   [m3 kg-1 Pa-1]
real, parameter :: SPV401 = -4.3209400767e-06*Pa2kb           ! SpV zs**4 * P coef.   [m3 kg-1 Pa-1]
real, parameter :: SPV011 = 1.5355844621e-05*(I_Ts*Pa2kb)     ! SpV T * P coef. [m3 kg-1 degC-1 Pa-1]
real, parameter :: SPV111 = 2.0914122241e-06*(I_Ts*Pa2kb)     ! SpV zs * T * P coef. [m3 kg-1 degC-1 Pa-1]
real, parameter :: SPV211 = -5.7751479725e-06*(I_Ts*Pa2kb)    ! SpV zs**2 * T * P coef. [m3 kg-1 degC-1 Pa-1]
real, parameter :: SPV311 = 1.0767234341e-06*(I_Ts*Pa2kb)     ! SpV zs**3 * T * P coef. [m3 kg-1 degC-1 Pa-1]
real, parameter :: SPV021 = -9.6659393016e-06*(I_Ts**2*Pa2kb) ! SpV T**2 * P coef. [m3 kg-1 degC-2 Pa-1]
real, parameter :: SPV121 = -7.0686982208e-07*(I_Ts**2*Pa2kb) ! SpV zs * T**2 * P coef. [m3 kg-1 degC-2 Pa-1]
real, parameter :: SPV221 = 1.4488066593e-06*(I_Ts**2*Pa2kb)  ! SpV zs**2 * T**2 * P coef. [m3 kg-1 degC-2 Pa-1]
real, parameter :: SPV031 = 3.1134283336e-06*(I_Ts**3*Pa2kb)  ! SpV T**3 * P coef. [m3 kg-1 degC-3 Pa-1]
real, parameter :: SPV131 = 7.9562529879e-08*(I_Ts**3*Pa2kb)  ! SpV zs * T**3 * P coef. [m3 kg-1 degC-3 Pa-1]
real, parameter :: SPV041 = -5.6590253863e-07*(I_Ts**4*Pa2kb) ! SpV T**4 * P coef. [m3 kg-1 degC-4 Pa-1]
real, parameter :: SPV002 = 1.0500241168e-06*Pa2kb**2         ! SpV P**2 coef.        [m3 kg-1 Pa-2]
real, parameter :: SPV102 = 1.9600661704e-06*Pa2kb**2         ! SpV zs * P**2 coef.   [m3 kg-1 Pa-2]
real, parameter :: SPV202 = -2.1666693382e-06*Pa2kb**2        ! SpV zs**2 * P**2 coef. [m3 kg-1 Pa-2]
real, parameter :: SPV012 = -3.8541359685e-06*(I_Ts*Pa2kb**2) ! SpV T * P**2 coef. [m3 kg-1 degC-1 Pa-2]
real, parameter :: SPV112 = 1.0157632247e-06*(I_Ts*Pa2kb**2)  ! SpV zs * T * P**2 coef. [m3 kg-1 degC-1 Pa-2]
real, parameter :: SPV022 = 1.7178343158e-06*(I_Ts**2*Pa2kb**2) ! SpV T**2 * P**2 coef. [m3 kg-1 degC-2 Pa-2]
real, parameter :: SPV003 = -4.1503454190e-07*Pa2kb**3        ! SpV P**3 coef.        [m3 kg-1 Pa-3]
real, parameter :: SPV103 = 3.5627020989e-07*Pa2kb**3         ! SpV zs * P**3 coef.   [m3 kg-1 Pa-3]
real, parameter :: SPV013 = -1.1293871415e-07*(I_Ts*Pa2kb**3) ! SpV T * P**3 coef. [m3 kg-1 degC-1 Pa-3]

real, parameter :: ALP000 =    SPV010   ! Constant in the dSpV_dT fit               [m3 kg-1 degC-1]
real, parameter :: ALP100 =    SPV110   ! dSpV_dT fit zs coef.                      [m3 kg-1 degC-1]
real, parameter :: ALP200 =    SPV210   ! dSpV_dT fit zs**2 coef.                   [m3 kg-1 degC-1]
real, parameter :: ALP300 =    SPV310   ! dSpV_dT fit zs**3 coef.                   [m3 kg-1 degC-1]
real, parameter :: ALP400 =    SPV410   ! dSpV_dT fit zs**4 coef.                   [m3 kg-1 degC-1]
real, parameter :: ALP500 =    SPV510   ! dSpV_dT fit zs**5 coef.                   [m3 kg-1 degC-1]
real, parameter :: ALP010 = 2.*SPV020   ! dSpV_dT fit T coef.                       [m3 kg-1 degC-2]
real, parameter :: ALP110 = 2.*SPV120   ! dSpV_dT fit zs * T coef.                  [m3 kg-1 degC-2]
real, parameter :: ALP210 = 2.*SPV220   ! dSpV_dT fit zs**2 * T coef.               [m3 kg-1 degC-2]
real, parameter :: ALP310 = 2.*SPV320   ! dSpV_dT fit zs**3 * T coef.               [m3 kg-1 degC-2]
real, parameter :: ALP410 = 2.*SPV420   ! dSpV_dT fit zs**4 * T coef.               [m3 kg-1 degC-2]
real, parameter :: ALP020 = 3.*SPV030   ! dSpV_dT fit T**2 coef.                    [m3 kg-1 degC-3]
real, parameter :: ALP120 = 3.*SPV130   ! dSpV_dT fit zs * T**2 coef.               [m3 kg-1 degC-3]
real, parameter :: ALP220 = 3.*SPV230   ! dSpV_dT fit zs**2 * T**2 coef.            [m3 kg-1 degC-3]
real, parameter :: ALP320 = 3.*SPV330   ! dSpV_dT fit zs**3 * T**2 coef.            [m3 kg-1 degC-3]
real, parameter :: ALP030 = 4.*SPV040   ! dSpV_dT fit T**3 coef.                    [m3 kg-1 degC-4]
real, parameter :: ALP130 = 4.*SPV140   ! dSpV_dT fit zs * T**3 coef.               [m3 kg-1 degC-4]
real, parameter :: ALP230 = 4.*SPV240   ! dSpV_dT fit zs**2 * T**3 coef.            [m3 kg-1 degC-4]
real, parameter :: ALP040 = 5.*SPV050   ! dSpV_dT fit T**4 coef.                    [m3 kg-1 degC-5]
real, parameter :: ALP140 = 5.*SPV150   ! dSpV_dT fit zs* * T**4 coef.              [m3 kg-1 degC-5]
real, parameter :: ALP050 = 6.*SPV060   ! dSpV_dT fit T**5 coef.                    [m3 kg-1 degC-6]
real, parameter :: ALP001 =    SPV011   ! dSpV_dT fit P coef.                  [m3 kg-1 degC-1 Pa-1]
real, parameter :: ALP101 =    SPV111   ! dSpV_dT fit zs * P coef.             [m3 kg-1 degC-1 Pa-1]
real, parameter :: ALP201 =    SPV211   ! dSpV_dT fit zs**2 * P coef.          [m3 kg-1 degC-1 Pa-1]
real, parameter :: ALP301 =    SPV311   ! dSpV_dT fit zs**3 * P coef.          [m3 kg-1 degC-1 Pa-1]
real, parameter :: ALP011 = 2.*SPV021   ! dSpV_dT fit T * P coef.              [m3 kg-1 degC-2 Pa-1]
real, parameter :: ALP111 = 2.*SPV121   ! dSpV_dT fit zs * T * P coef.         [m3 kg-1 degC-2 Pa-1]
real, parameter :: ALP211 = 2.*SPV221   ! dSpV_dT fit zs**2 * T * P coef.      [m3 kg-1 degC-2 Pa-1]
real, parameter :: ALP021 = 3.*SPV031   ! dSpV_dT fit T**2 * P coef.           [m3 kg-1 degC-3 Pa-1]
real, parameter :: ALP121 = 3.*SPV131   ! dSpV_dT fit zs * T**2 * P coef.      [m3 kg-1 degC-3 Pa-1]
real, parameter :: ALP031 = 4.*SPV041   ! dSpV_dT fit T**3 * P coef.           [m3 kg-1 degC-4 Pa-1]
real, parameter :: ALP002 =    SPV012   ! dSpV_dT fit P**2 coef.               [m3 kg-1 degC-1 Pa-2]
real, parameter :: ALP102 =    SPV112   ! dSpV_dT fit zs * P**2 coef.          [m3 kg-1 degC-1 Pa-2]
real, parameter :: ALP012 = 2.*SPV022   ! dSpV_dT fit T * P**2 coef.           [m3 kg-1 degC-2 Pa-2]
real, parameter :: ALP003 =    SPV013   ! dSpV_dT fit P**3 coef.               [m3 kg-1 degC-1 Pa-3]

real, parameter :: BET000 = 0.5*SPV100*r1_S0  ! Constant in the dSpV_dS fit          [m3 kg-1 ppt-1]
real, parameter :: BET100 =     SPV200*r1_S0  ! dSpV_dS fit zs coef.                 [m3 kg-1 ppt-1]
real, parameter :: BET200 = 1.5*SPV300*r1_S0  ! dSpV_dS fit zs**2 coef.              [m3 kg-1 ppt-1]
real, parameter :: BET300 = 2.0*SPV400*r1_S0  ! dSpV_dS fit zs**3 coef.              [m3 kg-1 ppt-1]
real, parameter :: BET400 = 2.5*SPV500*r1_S0  ! dSpV_dS fit zs**4 coef.              [m3 kg-1 ppt-1]
real, parameter :: BET500 = 3.0*SPV600*r1_S0  ! dSpV_dS fit zs**5 coef.              [m3 kg-1 ppt-1]
real, parameter :: BET010 = 0.5*SPV110*r1_S0  ! dSpV_dS fit T coef.           [m3 kg-1 ppt-1 degC-1]
real, parameter :: BET110 =     SPV210*r1_S0  ! dSpV_dS fit zs * T coef.      [m3 kg-1 ppt-1 degC-1]
real, parameter :: BET210 = 1.5*SPV310*r1_S0  ! dSpV_dS fit zs**2 * T coef.   [m3 kg-1 ppt-1 degC-1]
real, parameter :: BET310 = 2.0*SPV410*r1_S0  ! dSpV_dS fit zs**3 * T coef.   [m3 kg-1 ppt-1 degC-1]
real, parameter :: BET410 = 2.5*SPV510*r1_S0  ! dSpV_dS fit zs**4 * T coef.   [m3 kg-1 ppt-1 degC-1]
real, parameter :: BET020 = 0.5*SPV120*r1_S0  ! dSpV_dS fit T**2 coef.        [m3 kg-1 ppt-1 degC-2]
real, parameter :: BET120 =     SPV220*r1_S0  ! dSpV_dS fit zs * T**2 coef.   [m3 kg-1 ppt-1 degC-2]
real, parameter :: BET220 = 1.5*SPV320*r1_S0  ! dSpV_dS fit zs**2 * T**2 coef. [m3 kg-1 ppt-1 degC-2]
real, parameter :: BET320 = 2.0*SPV420*r1_S0  ! dSpV_dS fit zs**3 * T**2 coef. [m3 kg-1 ppt-1 degC-2]
real, parameter :: BET030 = 0.5*SPV130*r1_S0  ! dSpV_dS fit T**3 coef.        [m3 kg-1 ppt-1 degC-3]
real, parameter :: BET130 =     SPV230*r1_S0  ! dSpV_dS fit zs * T**3 coef.   [m3 kg-1 ppt-1 degC-3]
real, parameter :: BET230 = 1.5*SPV330*r1_S0  ! dSpV_dS fit zs**2 * T**3 coef. [m3 kg-1 ppt-1 degC-3]
real, parameter :: BET040 = 0.5*SPV140*r1_S0  ! dSpV_dS fit T**4 coef.        [m3 kg-1 ppt-1 degC-4]
real, parameter :: BET140 =     SPV240*r1_S0  ! dSpV_dS fit zs * T**4 coef.   [m3 kg-1 ppt-1 degC-4]
real, parameter :: BET050 = 0.5*SPV150*r1_S0  ! dSpV_dS fit T**5 coef.        [m3 kg-1 ppt-1 degC-5]
real, parameter :: BET001 = 0.5*SPV101*r1_S0  ! dSpV_dS fit P coef.             [m3 kg-1 ppt-1 Pa-1]
real, parameter :: BET101 =     SPV201*r1_S0  ! dSpV_dS fit zs * P coef.        [m3 kg-1 ppt-1 Pa-1]
real, parameter :: BET201 = 1.5*SPV301*r1_S0  ! dSpV_dS fit zs**2 * P coef.     [m3 kg-1 ppt-1 Pa-1]
real, parameter :: BET301 = 2.0*SPV401*r1_S0  ! dSpV_dS fit zs**3 * P coef.     [m3 kg-1 ppt-1 Pa-1]
real, parameter :: BET011 = 0.5*SPV111*r1_S0  ! dSpV_dS fit T * P coef.  [m3 kg-1 ppt-1 degC-1 Pa-1]
real, parameter :: BET111 =     SPV211*r1_S0  ! dSpV_dS fit zs * T * P coef. [m3 kg-1 ppt-1 degC-1 Pa-1]
real, parameter :: BET211 = 1.5*SPV311*r1_S0  ! dSpV_dS fit zs**2 * T * P coef. [m3 kg-1 ppt-1 degC-1 Pa-1]
real, parameter :: BET021 = 0.5*SPV121*r1_S0  ! dSpV_dS fit T**2 * P coef. [m3 kg-1 ppt-1 degC-2 Pa-1]
real, parameter :: BET121 =     SPV221*r1_S0  ! dSpV_dS fit zs * T**2 * P coef. [m3 kg-1 ppt-1 degC-2 Pa-1]
real, parameter :: BET031 = 0.5*SPV131*r1_S0  ! dSpV_dS fit T**3 * P coef. [m3 kg-1 ppt-1 degC-3 Pa-1]
real, parameter :: BET002 = 0.5*SPV102*r1_S0  ! dSpV_dS fit P**2 coef.          [m3 kg-1 ppt-1 Pa-2]
real, parameter :: BET102 =     SPV202*r1_S0  ! dSpV_dS fit zs * P**2 coef.     [m3 kg-1 ppt-1 Pa-2]
real, parameter :: BET012 = 0.5*SPV112*r1_S0  ! dSpV_dS fit T * P**2 coef. [m3 kg-1 ppt-1 degC-1 Pa-2]
real, parameter :: BET003 = 0.5*SPV103*r1_S0  ! dSpV_dS fit P**3 coef.          [m3 kg-1 ppt-1 Pa-3]
!>@}

!> The EOS_base implementation of the Roquet et al., 2015, equation of state
type, extends (EOS_base) :: Roquet_SpV_EOS

contains
  !> Implementation of the in-situ density as an elemental function [kg m-3]
  procedure :: density_elem => density_elem_Roquet_SpV
  !> Implementation of the in-situ density anomaly as an elemental function [kg m-3]
  procedure :: density_anomaly_elem => density_anomaly_elem_Roquet_SpV
  !> Implementation of the in-situ specific volume as an elemental function [m3 kg-1]
  procedure :: spec_vol_elem => spec_vol_elem_Roquet_SpV
  !> Implementation of the in-situ specific volume anomaly as an elemental function [m3 kg-1]
  procedure :: spec_vol_anomaly_elem => spec_vol_anomaly_elem_Roquet_SpV
  !> Implementation of the calculation of derivatives of density
  procedure :: calculate_density_derivs_elem => calculate_density_derivs_elem_Roquet_SpV
  !> Implementation of the calculation of second derivatives of density
  procedure :: calculate_density_second_derivs_elem => calculate_density_second_derivs_elem_Roquet_SpV
  !> Implementation of the calculation of derivatives of specific volume
  procedure :: calculate_specvol_derivs_elem => calculate_specvol_derivs_elem_Roquet_SpV
  !> Implementation of the calculation of compressibility
  procedure :: calculate_compress_elem => calculate_compress_elem_Roquet_SpV
  !> Implementation of the range query function
  procedure :: EOS_fit_range => EOS_fit_range_Roquet_SpV

  !> Local implementation of generic calculate_density_array for efficiency
  procedure :: calculate_density_array => calculate_density_array_Roquet_SpV
  !> Local implementation of generic calculate_spec_vol_array for efficiency
  procedure :: calculate_spec_vol_array => calculate_spec_vol_array_Roquet_SpV

end type Roquet_SpV_EOS


  interface
real elemental module function spec_vol_elem_Roquet_SpV(this, T, S, pressure)
  class(Roquet_SpV_EOS), intent(in) :: this     !< This EOS
  real,                  intent(in) :: T        !< Conservative temperature [degC]
  real,                  intent(in) :: S        !< Absolute salinity [g kg-1]
  real,                  intent(in) :: pressure !< pressure [Pa]

  ! Local variables
                 ! by an assumed salinity range [nondim]
                 ! specific volume at the reference temperature and salinity [m3 kg-1]
                 ! the surface pressure [m3 kg-1]
                 ! proportional to pressure [m3 kg-1 Pa-1]
                 ! proportional to pressure**2 [m3 kg-1 Pa-2]
                 ! proportional to pressure**3 [m3 kg-1 Pa-3]

  ! The following algorithm was published by Roquet et al. (2015), intended for use in non-Boussinesq ocean models.

  ! Conversions to the units used here.
end function spec_vol_elem_Roquet_SpV
real elemental module function spec_vol_anomaly_elem_Roquet_SpV(this, T, S, pressure, spv_ref)
  class(Roquet_SpV_EOS), intent(in) :: this     !< This EOS
  real,                  intent(in) :: T        !< Conservative temperature [degC]
  real,                  intent(in) :: S        !< Absolute salinity [g kg-1]
  real,                  intent(in) :: pressure !< pressure [Pa]
  real,                  intent(in) :: spv_ref  !< A reference specific volume [m3 kg-1]

  ! Local variables
                 ! by an assumed salinity range [nondim]
                 ! specific volume at the reference temperature and salinity [m3 kg-1]
                 ! the surface pressure [m3 kg-1]
                 ! proportional to pressure [m3 kg-1 Pa-1]
                 ! proportional to pressure**2 [m3 kg-1 Pa-2]
                 ! proportional to pressure**3 [m3 kg-1 Pa-3]

  ! The following algorithm was published by Roquet et al. (2015), intended for use in non-Boussinesq ocean models.

  ! Conversions to the units used here.
end function spec_vol_anomaly_elem_Roquet_SpV
real elemental module function density_elem_Roquet_SpV(this, T, S, pressure)
  class(Roquet_SpV_EOS), intent(in) :: this     !< This EOS
  real,                  intent(in) :: T        !< Conservative temperature [degC]
  real,                  intent(in) :: S        !< Absolute salinity [g kg-1]
  real,                  intent(in) :: pressure !< Pressure [Pa]

  ! Local variables

end function density_elem_Roquet_SpV
real elemental module function density_anomaly_elem_Roquet_SpV(this, T, S, pressure, rho_ref)
  class(Roquet_SpV_EOS), intent(in) :: this     !< This EOS
  real,                  intent(in) :: T        !< Conservative temperature [degC]
  real,                  intent(in) :: S        !< Absolute salinity [g kg-1]
  real,                  intent(in) :: pressure !< Pressure [Pa]
  real,                  intent(in) :: rho_ref  !< A reference density [kg m-3]

  ! Local variables

end function density_anomaly_elem_Roquet_SpV
elemental module subroutine calculate_specvol_derivs_elem_Roquet_SpV(this, T, S, pressure, dSV_dT, dSV_dS)
  class(Roquet_SpV_EOS), intent(in)    :: this     !< This EOS
  real,                  intent(in)    :: T        !< Conservative temperature [degC]
  real,                  intent(in)    :: S        !< Absolute salinity [g kg-1]
  real,                  intent(in)    :: pressure !< Pressure [Pa]
  real,                  intent(inout) :: dSV_dT   !< The partial derivative of specific volume with
                                                   !! conservative temperature [m3 kg-1 degC-1]
  real,                  intent(inout) :: dSV_dS   !< The partial derivative of specific volume with
                                                   !! absolute salinity [m3 kg-1 ppt-1]

                  ! by an assumed salinity range [nondim]
                  ! from temperature anomalies at the surface pressure [m3 kg-1 degC-1]
                  ! that is proportional to pressure [m3 kg-1 degC-1 Pa-1]
                  ! that is proportional to pressure**2 [m3 kg-1 degC-1 Pa-2]
                  ! that is proportional to pressure**3 [m3 kg-1 degC-1 Pa-3]
                  ! salinity [m3 kg-1 ppt-1] from temperature anomalies at the surface pressure
                  ! salinity [m3 kg-1 ppt-1 Pa-1] proportional to pressure
                  ! salinity [m3 kg-1 ppt-1 Pa-2] proportional to pressure**2
                  ! salinity [m3 kg-1 ppt-1 Pa-3] proportional to pressure**3

  ! Conversions to the units used here.
end subroutine calculate_specvol_derivs_elem_Roquet_SpV
elemental module subroutine calculate_density_derivs_elem_Roquet_SpV(this, T, S, pressure, drho_dT, drho_dS)
  class(Roquet_SpV_EOS), intent(in)  :: this     !< This EOS
  real,                  intent(in)  :: T        !< Conservative temperature [degC]
  real,                  intent(in)  :: S        !< Absolute salinity [g kg-1]
  real,                  intent(in)  :: pressure !< pressure [Pa]
  real,                  intent(out) :: drho_dT  !< The partial derivative of density with
                                                  !! conservative temperature [kg m-3 degC-1]
  real,                  intent(out) :: drho_dS  !< The partial derivative of density with
                                                 !! absolute salinity [kg m-3 ppt-1]

  ! Local variables
                                       ! conservative temperature [m3 kg-1 degC-1]
                                       ! absolute salinity [m3 kg-1 ppt-1]

end subroutine calculate_density_derivs_elem_Roquet_SpV
elemental module subroutine calculate_compress_elem_Roquet_SpV(this, T, S, pressure, rho, drho_dp)
  class(Roquet_SpV_EOS), intent(in)  :: this     !< This EOS
  real,                  intent(in)  :: T        !< Conservative temperature [degC]
  real,                  intent(in)  :: S        !< Absolute salinity [g kg-1]
  real,                  intent(in)  :: pressure !< pressure [Pa]
  real,                  intent(out) :: rho      !< In situ density [kg m-3]
  real,                  intent(out) :: drho_dp  !< The partial derivative of density with pressure
                                                 !! (also the inverse of the square of sound speed)
                                                 !! [s2 m-2]

  ! Local variables
                 ! by an assumed salinity range [nondim]
                 ! pressure [m3 kg-1 Pa-1]
                 ! pressure [m3 kg-1 Pa-1]
                 ! specific volume at the reference temperature and salinity [m3 kg-1]
                 ! the surface pressure [m3 kg-1]
                 ! proportional to pressure [m3 kg-1 Pa-1]
                 ! proportional to pressure**2 [m3 kg-1 Pa-2]
                 ! proportional to pressure**3 [m3 kg-1 Pa-3]

  ! The following algorithm was published by Roquet et al. (2015), intended for use
  ! with NEMO, but it is not necessarily the algorithm used in NEMO ocean model.

  ! Conversions to the units used here.
end subroutine calculate_compress_elem_Roquet_SpV
elemental module subroutine calc_spec_vol_second_derivs_elem_Roquet_SpV(T, S, P, &
                           dSV_ds_ds, dSV_ds_dt, dSV_dt_dt, dSV_ds_dp, dSV_dt_dp)
  real, intent(in)    :: T          !< Conservative temperature [degC]
  real, intent(in)    :: S          !< Absolute salinity [g kg-1]
  real, intent(in)    :: P          !< Pressure [Pa]
  real, intent(inout) :: dSV_ds_ds  !< Second derivative of specific volume with respect
                                    !! to salinity [m3 kg-1 ppt-2]
  real, intent(inout) :: dSV_ds_dt  !< Second derivative of specific volume with respect
                                    !! to salinity and temperature [m3 kg-1 ppt-1 degC-1]
  real, intent(inout) :: dSV_dt_dt  !< Second derivative of specific volume with respect
                                    !! to temperature [m3 kg-1 degC-2]
  real, intent(inout) :: dSV_ds_dp  !< Second derivative of specific volume with respect to pressure
                                    !! and salinity [m3 kg-1 ppt-1 Pa-1]
  real, intent(inout) :: dSV_dt_dp  !< Second derivative of specific volume with respect to pressure
                                    !! and temperature [m3 kg-1 degC-1 Pa-1]
  ! Local variables
                  ! by an assumed salinity range [nondim]

  ! Conversions to the units used here.
end subroutine calc_spec_vol_second_derivs_elem_Roquet_SpV
elemental module subroutine calculate_density_second_derivs_elem_Roquet_SpV(this, T, S, pressure, &
                               drho_ds_ds, drho_ds_dt, drho_dt_dt, drho_ds_dp, drho_dt_dp)
  class(Roquet_SpV_EOS), intent(in)    :: this       !< This EOS
  real,                  intent(in)    :: T          !< Conservative temperature [degC]
  real,                  intent(in)    :: S          !< Absolute salinity [g kg-1]
  real,                  intent(in)    :: pressure   !< Pressure [Pa]
  real,                  intent(inout) :: drho_ds_ds !< Second derivative of density with respect
                                                     !! to salinity [kg m-3 ppt-2]
  real,                  intent(inout) :: drho_ds_dt !< Second derivative of density with respect
                                                     !! to salinity and temperature [kg m-3 ppt-1 degC-1]
  real,                  intent(inout) :: drho_dt_dt !< Second derivative of density with respect
                                                     !! to temperature [kg m-3 degC-2]
  real,                  intent(inout) :: drho_ds_dp !< Second derivative of density with respect to pressure
                                                     !! and salinity [kg m-3 ppt-1 Pa-1] = [s2 m-2 ppt-1]
  real,                  intent(inout) :: drho_dt_dp !< Second derivative of density with respect to pressure
                                                     !! and temperature [kg m-3 degC-1 Pa-1] = [s2 m-2 degC-1]

  ! Local variables
                    ! (also the inverse of the square of sound speed)  [s2 m-2]
                    ! conservative temperature [m3 kg-1 degC-1]
                    ! absolute salinity [m3 kg-1 ppt-1]
                    ! to salinity [m3 kg-1 ppt-2]
                    ! to salinity and temperature [m3 kg-1 ppt-1 degC-1]
                    ! to temperature [m3 kg-1 degC-2]
                    ! and salinity [m3 kg-1 ppt-1 Pa-1]
                    ! and temperature [m3 kg-1 degC-1 Pa-1]

end subroutine calculate_density_second_derivs_elem_Roquet_SpV
module subroutine EoS_fit_range_Roquet_SpV(this, T_min, T_max, S_min, S_max, p_min, p_max)
  class(Roquet_SpV_EOS), intent(in)    :: this       !< This EOS
  real, optional, intent(out) :: T_min !< The minimum conservative temperature over which this EoS is fitted [degC]
  real, optional, intent(out) :: T_max !< The maximum conservative temperature over which this EoS is fitted [degC]
  real, optional, intent(out) :: S_min !< The minimum absolute salinity over which this EoS is fitted [g kg-1]
  real, optional, intent(out) :: S_max !< The maximum absolute salinity over which this EoS is fitted [g kg-1]
  real, optional, intent(out) :: p_min !< The minimum pressure over which this EoS is fitted [Pa]
  real, optional, intent(out) :: p_max !< The maximum pressure over which this EoS is fitted [Pa]

end subroutine EoS_fit_range_Roquet_SpV
module subroutine calculate_density_array_Roquet_SpV(this, T, S, pressure, rho, start, npts, rho_ref)
  class(Roquet_SpV_EOS),  intent(in) :: this  !< This EOS
  real, dimension(:), intent(in)  :: T        !< Potential temperature relative to the surface [degC]
  real, dimension(:), intent(in)  :: S        !< Salinity [PSU]
  real, dimension(:), intent(in)  :: pressure !< Pressure [Pa]
  real, dimension(:), intent(out) :: rho      !< In situ density [kg m-3]
  integer,            intent(in)  :: start    !< The starting index for calculations
  integer,            intent(in)  :: npts     !< The number of values to calculate
  real,     optional, intent(in)  :: rho_ref  !< A reference density [kg m-3]

  ! Local variables

end subroutine calculate_density_array_Roquet_SpV
module subroutine calculate_spec_vol_array_Roquet_SpV(this, T, S, pressure, specvol, start, npts, spv_ref)
  class(Roquet_SpV_EOS),  intent(in) :: this  !< This EOS
  real, dimension(:), intent(in)  :: T        !< Potential temperature relative to the surface [degC]
  real, dimension(:), intent(in)  :: S        !< Salinity [PSU]
  real, dimension(:), intent(in)  :: pressure !< Pressure [Pa]
  real, dimension(:), intent(out) :: specvol  !< In situ specific volume [m3 kg-1]
  integer,            intent(in)  :: start    !< The starting index for calculations
  integer,            intent(in)  :: npts     !< The number of values to calculate
  real,     optional, intent(in)  :: spv_ref  !< A reference specific volume [m3 kg-1]

  ! Local variables

end subroutine calculate_spec_vol_array_Roquet_SpV
  end interface

end module MOM_EOS_Roquet_Spv
