! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> Functions to convert between conservative and potential temperature
module MOM_temperature_convert

implicit none ; private

public poTemp_to_consTemp, consTemp_to_poTemp

!>@{ Parameters in the temperature conversion code
real, parameter :: Sprac_Sref = (35.0/35.16504) ! The TEOS 10 conversion factor to go from
                                    ! reference salinity to practical salinity [nondim]
real, parameter :: I_S0 = 0.025*Sprac_Sref ! The inverse of a plausible range of oceanic salinities [kg g-1]
real, parameter :: I_Ts = 0.025          ! The inverse of a plausible range of oceanic temperatures [degC-1]
real, parameter :: I_cp0 = 1.0/3991.86795711963 ! The inverse of the "specific heat" for use
                      !  with Conservative Temperature, as defined with TEOS10 [degC kg J-1]

! The following are coefficients of contributions to conservative temperature as a function of the square root
! of normalized absolute salinity with an offset (zS) and potential temperature (T) with a contribution
! Hab * zS**a * T**b.  The numbers here are copied directly from the corresponding gsw module, but
! the expressions here do not use the same nondimensionalization for pressure or temperature as they do.

real, parameter :: H00 = 61.01362420681071*I_cp0             ! Tp to Tc fit constant [degC]
real, parameter :: H01 = 168776.46138048015*(I_cp0*I_Ts)     ! Tp to Tc fit T coef. [nondim]
real, parameter :: H02 = -2735.2785605119625*(I_cp0*I_Ts**2) ! Tp to Tc fit T**2 coef. [degC-1]
real, parameter :: H03 = 2574.2164453821433*(I_cp0*I_Ts**3)  ! Tp to Tc fit T**3 coef. [degC-2]
real, parameter :: H04 = -1536.6644434977543*(I_cp0*I_Ts**4) ! Tp to Tc fit T**4 coef. [degC-3]
real, parameter :: H05 = 545.7340497931629*(I_cp0*I_Ts**5)   ! Tp to Tc fit T**5 coef. [degC-4]
real, parameter :: H06 = -50.91091728474331*(I_cp0*I_Ts**6)  ! Tp to Tc fit T**6 coef. [degC-5]
real, parameter :: H07 = -18.30489878927802*(I_cp0*I_Ts**7)  ! Tp to Tc fit T**7 coef. [degC-6]
real, parameter :: H20 = 268.5520265845071*I_cp0             ! Tp to Tc fit zS**2 coef. [degC]
real, parameter :: H21 = -12019.028203559312*(I_cp0*I_Ts)    ! Tp to Tc fit zS**2 * T coef. [nondim]
real, parameter :: H22 = 3734.858026725145*(I_cp0*I_Ts**2)   ! Tp to Tc fit zS**2 * T**2 coef. [degC-1]
real, parameter :: H23 = -2046.7671145057618*(I_cp0*I_Ts**3) ! Tp to Tc fit zS**2 * T**3 coef. [degC-2]
real, parameter :: H24 = 465.28655623826234*(I_cp0*I_Ts**4)  ! Tp to Tc fit zS**2 * T**4 coef. [degC-3]
real, parameter :: H25 = -0.6370820302376359*(I_cp0*I_Ts**5) ! Tp to Tc fit zS**2 * T**5 coef. [degC-4]
real, parameter :: H26 = -10.650848542359153*(I_cp0*I_Ts**6) ! Tp to Tc fit zS**2 * T**6 coef. [degC-5]
real, parameter :: H30 = 937.2099110620707*I_cp0             ! Tp to Tc fit zS**3 coef. [degC]
real, parameter :: H31 = 588.1802812170108*(I_cp0*I_Ts)      ! Tp to Tc fit zS** 3* T coef. [nondim]
real, parameter :: H32 = 248.39476522971285*(I_cp0*I_Ts**2)  ! Tp to Tc fit zS**3 * T**2 coef. [degC-1]
real, parameter :: H33 = -3.871557904936333*(I_cp0*I_Ts**3)  ! Tp to Tc fit zS**3 * T**3 coef. [degC-2]
real, parameter :: H34 = -2.6268019854268356*(I_cp0*I_Ts**4) ! Tp to Tc fit zS**3 * T**4 coef. [degC-3]
real, parameter :: H40 = -1687.914374187449*I_cp0            ! Tp to Tc fit zS**4 coef. [degC]
real, parameter :: H41 = 936.3206544460336*(I_cp0*I_Ts)      ! Tp to Tc fit zS**4 * T coef. [nondim]
real, parameter :: H42 = -942.7827304544439*(I_cp0*I_Ts**2)  ! Tp to Tc fit zS**4 * T**2 coef. [degC-1]
real, parameter :: H43 = 369.4389437509002*(I_cp0*I_Ts**3)   ! Tp to Tc fit zS**4 * T**3 coef. [degC-2]
real, parameter :: H44 = -33.83664947895248*(I_cp0*I_Ts**4)  ! Tp to Tc fit zS**4 * T**4 coef. [degC-3]
real, parameter :: H45 = -9.987880382780322*(I_cp0*I_Ts**5)  ! Tp to Tc fit zS**4 * T**5 coef. [degC-4]
real, parameter :: H50 = 246.9598888781377*I_cp0             ! Tp to Tc fit zS**5 coef. [degC]
real, parameter :: H60 = 123.59576582457964*I_cp0            ! Tp to Tc fit zS**6 coef. [degC]
real, parameter :: H70 = -48.5891069025409*I_cp0             ! Tp to Tc fit zS**7 coef. [degC]

!>@}


  interface
elemental real module function poTemp_to_consTemp(T, Sa) result(Tc)
  real, intent(in) :: T  !< Potential temperature [degC]
  real, intent(in) :: Sa !< Absolute salinity [g kg-1]

  ! Local variables

end function poTemp_to_consTemp
elemental real module function dTc_dTp(T, Sa)
  real, intent(in) :: T  !< Potential temperature [degC]
  real, intent(in) :: Sa !< Absolute salinity [g kg-1]

  ! Local variables

end function dTc_dTp
elemental real module function consTemp_to_poTemp(Tc, Sa) result(Tp)
  real, intent(in) :: Tc !< Conservative temperature [degC]
  real, intent(in) :: Sa !< Absolute salinity [g kg-1]

  ! The following are coefficients in the nominator (TPNxx) or denominator (TPDxx) of a simple rational
  ! expression that approximately converts conservative temperature to potential temperature.

  ! Estimate the potential temperature and its derivative from an approximate rational function fit.
end function consTemp_to_poTemp
  end interface

end module MOM_temperature_convert
