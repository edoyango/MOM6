submodule (MOM_TFreeze) MOM_TFreeze_s
  implicit none
contains
module procedure calculate_TFreeze_linear_scalar
  T_Fr = (TFr_S0_P0 + dTFr_dS*S) + dTFr_dp*pres

end procedure calculate_TFreeze_linear_scalar
module procedure calculate_TFreeze_linear_array
  integer :: j
  do j=start,start+npts-1
    T_Fr(j) = (TFr_S0_P0 + dTFr_dS*S(j)) + dTFr_dp*pres(j)
  enddo

end procedure calculate_TFreeze_linear_array
module procedure calculate_TFreeze_Millero_scalar
  real, parameter :: cS1 = -0.0575      ! A term in the freezing point fit [degC PSU-1]
  real, parameter :: cS3_2 = 1.710523e-3 ! A term in the freezing point fit [degC PSU-3/2]
  real, parameter :: cS2 = -2.154996e-4 ! A term in the freezing point fit [degC PSU-2]
  real, parameter :: dTFr_dp = -7.75e-8 ! Derivative of freezing point with pressure [degC Pa-1]
  T_Fr = S*(cS1 + (cS3_2 * sqrt(max(S, 0.0)) + cS2 * S)) + dTFr_dp*pres

end procedure calculate_TFreeze_Millero_scalar
module procedure calculate_TFreeze_Millero_array
  real, parameter :: cS1 = -0.0575      ! A term in the freezing point fit [degC PSU-1]
  real, parameter :: cS3_2 = 1.710523e-3 ! A term in the freezing point fit [degC PSU-3/2]
  real, parameter :: cS2 = -2.154996e-4 ! A term in the freezing point fit [degC PSU-2]
  real, parameter :: dTFr_dp = -7.75e-8 ! Derivative of freezing point with pressure [degC Pa-1]
  integer :: j
  do j=start,start+npts-1
    T_Fr(j) = S(j)*(cS1 + (cS3_2 * sqrt(max(S(j), 0.0)) + cS2 * S(j))) + &
              dTFr_dp*pres(j)
  enddo

end procedure calculate_TFreeze_Millero_array
module procedure calculate_TFreeze_TEOS_poly_scalar
  real, dimension(1) :: S0    ! Salinity at a point [g kg-1]
  real, dimension(1) :: pres0 ! Pressure at a point [Pa]
  real, dimension(1) :: tfr0  ! The freezing temperature [degC]
  S0(1) = S
  pres0(1) = pres

  call calculate_TFreeze_TEOS_poly_array(S0, pres0, tfr0, 1, 1)
  T_Fr = tfr0(1)

end procedure calculate_TFreeze_TEOS_poly_scalar
module procedure calculate_TFreeze_TEOS_poly_array
  real :: Sa    ! Absolute salinity [g kg-1] = [ppt]
  real :: rS    ! Square root of salinity [ppt1/2]
  real, parameter :: TF00 =  0.017947064327968736  ! Freezing point coefficient [degC]
  real, parameter :: TF20 = -6.076099099929818e-2  ! Freezing point coefficient [degC ppt-1]
  real, parameter :: TF30 =  4.883198653547851e-3  ! Freezing point coefficient [degC ppt-3/2]
  real, parameter :: TF40 = -1.188081601230542e-3  ! Freezing point coefficient [degC ppt-2]
  real, parameter :: TF50 =  1.334658511480257e-4  ! Freezing point coefficient [degC ppt-5/2]
  real, parameter :: TF60 = -8.722761043208607e-6  ! Freezing point coefficient [degC ppt-3]
  real, parameter :: TF70 =  2.082038908808201e-7  ! Freezing point coefficient [degC ppt-7/2]
  real, parameter :: TF01 = -7.389420998107497e-8  ! Freezing point coefficient [degC Pa-1]
  real, parameter :: TF21 = -9.891538123307282e-11 ! Freezing point coefficient [degC ppt-1 Pa-1]
  real, parameter :: TF31 = -8.987150128406496e-13 ! Freezing point coefficient [degC ppt-3/2 Pa-1]
  real, parameter :: TF41 =  1.054318231187074e-12 ! Freezing point coefficient [degC ppt-2 Pa-1]
  real, parameter :: TF51 =  3.850133554097069e-14 ! Freezing point coefficient [degC ppt-5/2 Pa-1]
  real, parameter :: TF61 = -2.079022768390933e-14 ! Freezing point coefficient [degC ppt-3 Pa-1]
  real, parameter :: TF71 =  1.242891021876471e-15 ! Freezing point coefficient [degC ppt-7/2 Pa-1]
  real, parameter :: TF02 = -2.110913185058476e-16 ! Freezing point coefficient [degC Pa-2]
  real, parameter :: TF22 =  3.831132432071728e-19 ! Freezing point coefficient [degC ppt-1 Pa-2]
  real, parameter :: TF32 =  1.065556599652796e-19 ! Freezing point coefficient [degC ppt-3/2 Pa-2]
  real, parameter :: TF42 = -2.078616693017569e-20 ! Freezing point coefficient [degC ppt-2 Pa-2]
  real, parameter :: TF52 =  1.596435439942262e-21 ! Freezing point coefficient [degC ppt-5/2 Pa-2]
  real, parameter :: TF03 =  2.295491578006229e-25 ! Freezing point coefficient [degC Pa-3]
  real, parameter :: TF23 = -7.997496801694032e-27 ! Freezing point coefficient [degC ppt-1 Pa-3]
  real, parameter :: TF33 =  8.756340772729538e-28 ! Freezing point coefficient [degC ppt-3/2 Pa-3]
  real, parameter :: TF43 =  1.338002171109174e-29 ! Freezing point coefficient [degC ppt-2 Pa-3]
  integer :: j
  do j=start,start+npts-1
    rS = sqrt(max(S(j), 0.0))
    T_Fr(j) =       (TF00 + S(j)*(TF20 + rS*(TF30 + rS*(TF40 + rS*(TF50 + rS*(TF60 + rS*TF70)))))) &
        + pres(j)*( (TF01 + S(j)*(TF21 + rS*(TF31 + rS*(TF41 + rS*(TF51 + rS*(TF61 + rS*TF71)))))) &
         + pres(j)*((TF02 + S(j)*(TF22 + rS*(TF32 + rS*(TF42 + rS* TF52)))) &
          + pres(j)*(TF03 + S(j)*(TF23 + rS*(TF33 + rS* TF43))) ) )
  enddo

end procedure calculate_TFreeze_TEOS_poly_array
module procedure calculate_TFreeze_teos10_scalar
  real, dimension(1) :: S0    ! Salinity at a point [g kg-1]
  real, dimension(1) :: pres0 ! Pressure at a point [Pa]
  real, dimension(1) :: tfr0  ! The freezing temperature [degC]
  S0(1) = S
  pres0(1) = pres

  call calculate_TFreeze_teos10_array(S0, pres0, tfr0, 1, 1)
  T_Fr = tfr0(1)

end procedure calculate_TFreeze_teos10_scalar
module procedure calculate_TFreeze_teos10_array
  real, parameter :: Pa2db  = 1.e-4  ! The conversion factor from Pa to dbar [dbar Pa-1]
  real :: zp    ! Pressures in [dbar]
  integer :: j
  real, parameter :: saturation_fraction = 0.0 ! Air saturation fraction in seawater [nondim]
  do j=start,start+npts-1
    !Conversions
    zp = pres(j)* Pa2db         !Convert pressure from Pascal to decibar

    if (S(j) < -1.0e-10) cycle !Can we assume safely that this is a missing value?
    T_Fr(j) = gsw_ct_freezing_exact(S(j), zp, saturation_fraction)
  enddo

end procedure calculate_TFreeze_teos10_array
end submodule MOM_TFreeze_s
