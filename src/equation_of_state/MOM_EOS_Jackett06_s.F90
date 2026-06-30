submodule (MOM_EOS_Jackett06) MOM_EOS_Jackett06_s
  implicit none
contains
module procedure density_elem_Jackett06
  real :: num_STP ! State dependent part of the numerator of the rational expresion
  real :: den     ! Denominator of the rational expresion for density [nondim]
  real :: den_STP ! State dependent part of the denominator of the rational expresion
  real :: I_den   ! The inverse of the denominator of the rational expresion for density [nondim]
  real :: T2      ! Temperature squared [degC2]
  real :: S1_2    ! Limited square root of salinity [PSU1/2]
  S1_2 = sqrt(max(0.0,s))
  T2 = T*T

  num_STP = (T*(RN010 + T*(RN020 + T*RN030)) + &
             S*(RN100 + (T*RN110 + S*RN200)) ) + &
            pressure*(RN001 + ((T2*RN021 + S*RN101) + pressure*(RN002 + T2*RN022)))
  den = 1.0 + ((T*(RD010 + T*(RD020 + T*(RD030 + T* RD040))) + &
                S*(RD100 + (T*(RD110 + T2*RD130) + S1_2*(RD600 + T2*RD620))) ) + &
               pressure*(RD001 + pressure*T*(T2*RD032 + pressure*RD013)) )
  I_den = 1.0 / den

  density_elem_Jackett06 = (RN000 + num_STP)*I_den

end procedure density_elem_Jackett06
module procedure density_anomaly_elem_Jackett06
  real :: num_STP ! State dependent part of the numerator of the rational expresion
  real :: den     ! Denominator of the rational expresion for density [nondim]
  real :: den_STP ! State dependent part of the denominator of the rational expresion
  real :: I_den   ! The inverse of the denominator of the rational expresion for density [nondim]
  real :: T2      ! Temperature squared [degC2]
  real :: S1_2    ! Limited square root of salinity [PSU1/2]
  real :: rho0    ! The surface density of fresh water at 0 degC, perhaps less the refernce density [kg m-3]
  S1_2 = sqrt(max(0.0,s))
  T2 = T*T

  num_STP = (T*(RN010 + T*(RN020 + T*RN030)) + &
             S*(RN100 + (T*RN110 + S*RN200)) ) + &
            pressure*(RN001 + ((T2*RN021 + S*RN101) + pressure*(RN002 + T2*RN022)))
  den = 1.0 + ((T*(RD010 + T*(RD020 + T*(RD030 + T* RD040))) + &
                S*(RD100 + (T*(RD110 + T2*RD130) + S1_2*(RD600 + T2*RD620))) ) + &
               pressure*(RD001 + pressure*T*(T2*RD032 + pressure*RD013)) )
  I_den = 1.0 / den

  rho0 = RN000 - rho_ref*den

  density_anomaly_elem_Jackett06 = (rho0 + num_STP)*I_den

end procedure density_anomaly_elem_Jackett06
module procedure spec_vol_elem_Jackett06
  real :: num_STP ! State dependent part of the numerator of the rational expresion
  real :: den_STP ! State dependent part of the denominator of the rational expresion
  real :: I_num   ! The inverse of the numerator of the rational expresion for density [nondim]
  real :: T2      ! Temperature squared [degC2]
  real :: S1_2    ! Limited square root of salinity [PSU1/2]
  S1_2 = sqrt(max(0.0,s))
  T2 = T*T

  num_STP = (T*(RN010 + T*(RN020 + T*RN030)) + &
             S*(RN100 + (T*RN110 + S*RN200)) ) + &
            pressure*(RN001 + ((T2*RN021 + S*RN101) + pressure*(RN002 + T2*RN022)))
  den_STP = (T*(RD010 + T*(RD020 + T*(RD030 + T* RD040))) + &
             S*(RD100 + (T*(RD110 + T2*RD130) + S1_2*(RD600 + T2*RD620))) ) + &
            pressure*(RD001 + pressure*T*(T2*RD032 + pressure*RD013))
  I_num = 1.0 / (RN000 + num_STP)

  spec_vol_elem_Jackett06 = (1.0 + den_STP) * I_num

end procedure spec_vol_elem_Jackett06
module procedure spec_vol_anomaly_elem_Jackett06
  real :: num_STP ! State dependent part of the numerator of the rational expresion
  real :: den_STP ! State dependent part of the denominator of the rational expresion
  real :: I_num   ! The inverse of the numerator of the rational expresion for density [nondim]
  real :: T2      ! Temperature squared [degC2]
  real :: S1_2    ! Limited square root of salinity [PSU1/2]
  S1_2 = sqrt(max(0.0,s))
  T2 = T*T

  num_STP = (T*(RN010 + T*(RN020 + T*RN030)) + &
             S*(RN100 + (T*RN110 + S*RN200)) ) + &
            pressure*(RN001 + ((T2*RN021 + S*RN101) + pressure*(RN002 + T2*RN022)))
  den_STP = (T*(RD010 + T*(RD020 + T*(RD030 + T* RD040))) + &
             S*(RD100 + (T*(RD110 + T2*RD130) + S1_2*(RD600 + T2*RD620))) ) + &
            pressure*(RD001 + pressure*T*(T2*RD032 + pressure*RD013))
  I_num = 1.0 / (RN000 + num_STP)

  ! This form is slightly more complicated, but it cancels the leading terms better.
  spec_vol_anomaly_elem_Jackett06 = ((1.0 - spv_ref*RN000) + (den_STP - spv_ref*num_STP)) * I_num

end procedure spec_vol_anomaly_elem_Jackett06
module procedure calculate_density_derivs_elem_Jackett06
  real :: num     ! Numerator of the rational expresion for density [kg m-3]
  real :: den     ! Denominator of the rational expresion for density [nondim]
  real :: I_denom2 ! The inverse of the square of the denominator of the rational expression
  real :: dnum_dT ! The derivative of num with potential temperature [kg m-3 degC-1]
  real :: dnum_dS ! The derivative of num with salinity [kg m-3 PSU-1]
  real :: dden_dT ! The derivative of den with potential temperature [degC-1]
  real :: dden_dS ! The derivative of den with salinity PSU-1]
  real :: T2      ! Temperature squared [degC2]
  real :: S1_2    ! Limited square root of salinity [PSU1/2]
  S1_2 = sqrt(max(0.0,s))
  T2 = T*T

  num = RN000 + ((T*(RN010 + T*(RN020 + T*RN030)) + &
                  S*(RN100 + (T*RN110 + S*RN200)) ) + &
                 pressure*(RN001 + ((T2*RN021 + S*RN101) + pressure*(RN002 + T2*RN022))) )
  den = 1.0 + ((T*(RD010 + T*(RD020 + T*(RD030 + T* RD040))) + &
                S*(RD100 + (T*(RD110 + T2*RD130) + S1_2*(RD600 + T2*RD620))) ) + &
               pressure*(RD001 + pressure*T*(T2*RD032 + pressure*RD013)) )

  dnum_dT = ((RN010 + T*(2.*RN020 + T*(3.*RN030))) + S*RN110) + &
            pressure*T*(2.*RN021 + pressure*(2.*RN022))
  dnum_dS = (RN100 + (T*RN110 + S*(2.*RN200))) + pressure*RN101
  dden_dT = ((RD010 + T*((2.*RD020) + T*((3.*RD030) + T*(4.*RD040)))) + &
             S*((RD110 + T2*(3.*RD130)) + S1_2*T*(2.*RD620)) ) + &
            pressure**2*(T2*3.*RD032 + pressure*RD013)
  dden_dS = RD100 + (T*(RD110 + T2*RD130) + S1_2*(1.5*RD600 + T2*(1.5*RD620)))
  I_denom2 = 1.0 / den**2

  ! rho = num / den
  drho_dT = (dnum_dT * den - num * dden_dT) * I_denom2
  drho_dS = (dnum_dS * den - num * dden_dS) * I_denom2

end procedure calculate_density_derivs_elem_Jackett06
module procedure calculate_density_second_derivs_elem_Jackett06
  real :: num         ! Numerator of the rational expresion for density [kg m-3]
  real :: den         ! Denominator of the rational expresion for density [nondim]
  real :: I_num2      ! The inverse of the square of the numerator of the rational expression
  real :: dnum_dT     ! The derivative of num with potential temperature [kg m-3 degC-1]
  real :: dnum_dS     ! The derivative of num with salinity [kg m-3 PSU-1]
  real :: dden_dT     ! The derivative of den with potential temperature [degC-1]
  real :: dden_dS     ! The derivative of den with salinity PSU-1]
  real :: dnum_dp     ! The derivative of num with pressure [kg m-3 dbar-1]
  real :: dden_dp     ! The derivative of det with pressure [dbar-1]
  real :: d2num_dT2   ! The second derivative of num with potential temperature [kg m-3 degC-2]
  real :: d2num_dT_dS ! The second derivative of num with potential temperature and
  real :: d2num_dS2   ! The second derivative of num with salinity [kg m-3 PSU-2]
  real :: d2num_dT_dp ! The second derivative of num with potential temperature and
  real :: d2num_dS_dp ! The second derivative of num with salinity and
  real :: d2den_dT2   ! The second derivative of den with potential temperature [degC-2]
  real :: d2den_dT_dS ! The second derivative of den with potential temperature and salinity [degC-1 PSU-1]
  real :: d2den_dS2   ! The second derivative of den with salinity [PSU-2]
  real :: d2den_dT_dp ! The second derivative of den with potential temperature and pressure [degC-1 dbar-1]
  real :: d2den_dS_dp ! The second derivative of den with salinity and pressure [PSU-1 dbar-1]
  real :: T2          ! Temperature squared [degC2]
  real :: S1_2        ! Limited square root of salinity [PSU1/2]
  real :: I_s12       ! The inverse of the square root of salinity [PSU-1/2]
  real :: I_denom2    ! The inverse of the square of the denominator of the rational expression
  real :: I_denom3    ! The inverse of the cube of the denominator of the rational expression
  S1_2 = sqrt(max(0.0,s))
  T2 = T*T

  num = RN000 + ((T*(RN010 + T*(RN020 + T*RN030)) + &
                  S*(RN100 + (T*RN110 + S*RN200)) ) + &
                 pressure*(RN001 + ((T2*RN021 + S*RN101) + pressure*(RN002 + T2*RN022))) )
  den = 1.0 + ((T*(RD010 + T*(RD020 + T*(RD030 + T* RD040))) + &
                S*(RD100 + (T*(RD110 + T2*RD130) + S1_2*(RD600 + T2*RD620))) ) + &
               pressure*(RD001 + pressure*T*(T2*RD032 + pressure*RD013)) )
  ! rho = num*I_den

  dnum_dT = ((RN010 + T*(2.*RN020 + T*(3.*RN030))) + S*RN110) + &
            pressure*T*(2.*RN021 + pressure*(2.*RN022))
  dnum_dS = (RN100 + (T*RN110 + S*(2.*RN200))) + pressure*RN101
  dnum_dp = RN001 + ((T2*RN021 + S*RN101) + pressure*(2.*RN002 + T2*(2.*RN022)))
  d2num_dT2 = 2.*RN020 + T*(6.*RN030) + pressure*(2.*RN021 + pressure*(2.*RN022))
  d2num_dT_dS = RN110
  d2num_dS2 = 2.*RN200
  d2num_dT_dp = T*(2.*RN021 + pressure*(4.*RN022))
  d2num_dS_dp = RN101

  dden_dT = ((RD010 + T*((2.*RD020) + T*((3.*RD030) + T*(4.*RD040)))) + &
             S*((RD110 + T2*(3.*RD130)) + S1_2*T*(2.*RD620)) ) + &
            pressure**2*(T2*3.*RD032 + pressure*RD013)
  dden_dS = RD100 + (T*(RD110 + T2*RD130) + S1_2*(1.5*RD600 + T2*(1.5*RD620)))
  dden_dp = RD001 + pressure*T*(T2*(2.*RD032) + pressure*(3.*RD013))

  d2den_dT2 = (((2.*RD020) + T*((6.*RD030) + T*(12.*RD040))) + &
               S*(T*(6.*RD130) + S1_2*(2.*RD620)) ) + pressure**2*(T*(6.*RD032))
  d2den_dT_dS = (RD110 + T2*3.*RD130) + (T*S1_2)*(3.0*RD620)
  d2den_dT_dp = pressure*(T2*(6.*RD032) + pressure*(3.*RD013))
  d2den_dS_dp = 0.0

  ! The Jackett et al. 2006 equation of state is a fit to density, but it chooses a form that
  ! exhibits a singularity in the second derivatives with salinity for fresh water.  To avoid
  ! this, the square root of salinity can be treated with a floor such that the contribution from
  ! the S**1.5 terms to both the surface density and the secant bulk modulus are lost to roundoff.
  ! This salinity is given by (~1e-16/RD600)**(2/3) ~= 7e-8 PSU, or S1_2 ~= 2.6e-4
  I_S12 = 1.0 / (max(S1_2, 1.0e-4))
  d2den_dS2 = (0.75*RD600 + T2*(0.75*RD620)) * I_S12

  I_denom3 = 1.0 / den**3

  ! In deriving the following, it is useful to note that:
  !   drho_dp = (dnum_dp * den - num * dden_dp) / den**2
  !   drho_dT = (dnum_dT * den - num * dden_dT) / den**2
  !   drho_dS = (dnum_dS * den - num * dden_dS) / den**2
  drho_dS_dS = (den*(den*d2num_dS2 - 2.*dnum_dS*dden_dS) + num*(2.*dden_dS**2 - den*d2den_dS2)) * I_denom3
  drho_dS_dt = (den*(den*d2num_dT_dS - (dnum_dT*dden_dS + dnum_dS*dden_dT)) + &
                   num*(2.*dden_dT*dden_dS - den*d2den_dT_dS)) * I_denom3
  drho_dT_dT = (den*(den*d2num_dT2 - 2.*dnum_dT*dden_dT) + num*(2.*dden_dT**2 - den*d2den_dT2)) * I_denom3

  drho_dS_dp = (den*(den*d2num_dS_dp - (dnum_dp*dden_dS + dnum_dS*dden_dp)) + &
                   num*(2.*dden_dS*dden_dp - den*d2den_dS_dp)) * I_denom3
  drho_dT_dp = (den*(den*d2num_dT_dp - (dnum_dp*dden_dT + dnum_dT*dden_dp)) + &
                   num*(2.*dden_dT*dden_dp - den*d2den_dT_dp)) * I_denom3

end procedure calculate_density_second_derivs_elem_Jackett06
module procedure calculate_specvol_derivs_elem_Jackett06
  real :: num     ! Numerator of the rational expresion for density (not specific volume) [kg m-3]
  real :: den     ! Denominator of the rational expresion for density (not specific volume) [nondim]
  real :: I_num2  ! The inverse of the square of the numerator of the rational expression
  real :: dnum_dT ! The derivative of num with potential temperature [kg m-3 degC-1]
  real :: dnum_dS ! The derivative of num with salinity [kg m-3 PSU-1]
  real :: dden_dT ! The derivative of den with potential temperature [degC-1]
  real :: dden_dS ! The derivative of den with salinity PSU-1]
  real :: T2      ! Temperature squared [degC2]
  real :: S1_2    ! Limited square root of salinity [PSU1/2]
  S1_2 = sqrt(max(0.0,s))
  T2 = T*T

  num = RN000 + ((T*(RN010 + T*(RN020 + T*RN030)) + &
                  S*(RN100 + (T*RN110 + S*RN200)) ) + &
                 pressure*(RN001 + ((T2*RN021 + S*RN101) + pressure*(RN002 + T2*RN022))) )
  den = 1.0 + ((T*(RD010 + T*(RD020 + T*(RD030 + T* RD040))) + &
                S*(RD100 + (T*(RD110 + T2*RD130) + S1_2*(RD600 + T2*RD620))) ) + &
               pressure*(RD001 + pressure*T*(T2*RD032 + pressure*RD013)) )

  dnum_dT = ((RN010 + T*(2.*RN020 + T*(3.*RN030))) + S*RN110) + &
            pressure*T*(2.*RN021 + pressure*(2.*RN022))
  dnum_dS = (RN100 + (T*RN110 + S*(2.*RN200))) + pressure*RN101
  dden_dT = ((RD010 + T*((2.*RD020) + T*((3.*RD030) + T*(4.*RD040)))) + &
             S*((RD110 + T2*(3.*RD130)) + S1_2*T*(2.*RD620)) ) + &
            pressure**2*(T2*3.*RD032 + pressure*RD013)
  dden_dS = RD100 + (T*(RD110 + T2*RD130) + S1_2*(1.5*RD600 + T2*(1.5*RD620)))
  I_num2 = 1.0 / num**2

  ! SV = den / num
  dSV_dT = (num * dden_dT - dnum_dT * den) * I_num2
  dSV_dS = (num * dden_dS - dnum_dS * den) * I_num2

end procedure calculate_specvol_derivs_elem_Jackett06
module procedure calculate_compress_elem_Jackett06
  real :: num     ! Numerator of the rational expresion for density [kg m-3]
  real :: den     ! Denominator of the rational expresion for density [nondim]
  real :: I_den   ! The inverse of the denominator of the rational expression for density [nondim]
  real :: dnum_dp ! The derivative of num with pressure [kg m-3 dbar-1]
  real :: dden_dp ! The derivative of den with pressure [dbar-1]
  real :: T2      ! Temperature squared [degC2]
  real :: S1_2    ! Limited square root of salinity [PSU1/2]
  integer :: j
  S1_2 = sqrt(max(0.0,s))
  T2 = T*T

  num = RN000 + ((T*(RN010 + T*(RN020 + T*RN030)) + &
                  S*(RN100 + (T*RN110 + S*RN200)) ) + &
                 pressure*(RN001 + ((T2*RN021 + S*RN101) + pressure*(RN002 + T2*RN022))) )
  den = 1.0 + ((T*(RD010 + T*(RD020 + T*(RD030 + T* RD040))) + &
                S*(RD100 + (T*(RD110 + T2*RD130) + S1_2*(RD600 + T2*RD620))) ) + &
               pressure*(RD001 + pressure*T*(T2*RD032 + pressure*RD013)) )
  dnum_dp = RN001 + ((T2*RN021 + S*RN101) + pressure*(2.*RN002 + T2*(2.*RN022)))
  dden_dp = RD001 + pressure*T*(T2*(2.*RD032) + pressure*(3.*RD013))

  I_den  = 1.0 / den
  rho = num * I_den
  drho_dp = (dnum_dp * den - num * dden_dp) * I_den**2

end procedure calculate_compress_elem_Jackett06
module procedure EoS_fit_range_Jackett06
  if (present(T_min)) T_min = -4.5
  if (present(T_max)) T_max = 40.0
  if (present(S_min)) S_min =  0.0
  if (present(S_max)) S_max = 42.0
  if (present(p_min)) p_min = 0.0
  if (present(p_max)) p_max = 8.5e7

end procedure EoS_fit_range_Jackett06
end submodule MOM_EOS_Jackett06_s
