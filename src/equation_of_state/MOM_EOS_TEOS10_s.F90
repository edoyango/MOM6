submodule (MOM_EOS_TEOS10) MOM_EOS_TEOS10_s
  implicit none
contains
module procedure density_elem_TEOS10
  density_elem_TEOS10 = gsw_rho(S, T, pressure * Pa2db)

end procedure density_elem_TEOS10
module procedure density_anomaly_elem_TEOS10
  density_anomaly_elem_TEOS10 = gsw_rho(S, T, pressure * Pa2db)
  density_anomaly_elem_TEOS10 = density_anomaly_elem_TEOS10 - rho_ref

end procedure density_anomaly_elem_TEOS10
module procedure spec_vol_elem_TEOS10
  spec_vol_elem_TEOS10 = gsw_specvol(S, T, pressure * Pa2db)

end procedure spec_vol_elem_TEOS10
module procedure spec_vol_anomaly_elem_TEOS10
  spec_vol_anomaly_elem_TEOS10 = gsw_specvol(S, T, pressure * Pa2db) - spv_ref

end procedure spec_vol_anomaly_elem_TEOS10
module procedure calculate_density_derivs_elem_TEOS10
  real :: zs  ! Absolute salinity [g kg-1]
  real :: zt  ! Conservative temperature [degC]
  real :: zp  ! Pressure converted to decibars [dbar]
  zs = S
  zt = T
  zp = pressure * Pa2db      ! Convert pressure from Pascal to decibar
  ! The following conversions are unnecessary because the arguments are already the right variables.
  ! zs = gsw_sr_from_sp(S)   ! Uncomment to convert practical salinity to absolute salinity
  ! zt = gsw_ct_from_pt(S,T) ! Uncomment to convert potential temp to conservative temp

  call gsw_rho_first_derivatives(zs, zt, zp, drho_dsa=drho_dS, drho_dct=drho_dT)

end procedure calculate_density_derivs_elem_TEOS10
module procedure calculate_density_second_derivs_elem_TEOS10
  real :: zs  ! Absolute salinity [g kg-1]
  real :: zt  ! Conservative temperature [degC]
  real :: zp  ! Pressure converted to decibars [dbar]
  zs = S
  zt = T
  zp = pressure * Pa2db      ! Convert pressure from Pascal to decibar
  ! The following conversions are unnecessary because the arguments are already the right variables.
  ! zs = gsw_sr_from_sp(S)   ! Uncomment to convert practical salinity to absolute salinity
  ! zt = gsw_ct_from_pt(S,T) ! Uncomment to convert potential temp to conservative temp

  call gsw_rho_second_derivatives(zs, zt, zp, rho_sa_sa=drho_dS_dS, rho_sa_ct=drho_dS_dT, &
                                  rho_ct_ct=drho_dT_dT, rho_sa_p=drho_dS_dP, rho_ct_p=drho_dT_dP)

end procedure calculate_density_second_derivs_elem_TEOS10
module procedure calculate_specvol_derivs_elem_TEOS10
  real :: zs  ! Absolute salinity [g kg-1]
  real :: zt  ! Conservative temperature [degC]
  real :: zp  ! Pressure converted to decibars [dbar]
  zs = S
  zt = T
  zp = pressure * Pa2db      ! Convert pressure from Pascal to decibar
  ! The following conversions are unnecessary because the arguments are already the right variables.
  ! zs = gsw_sr_from_sp(S)   ! Uncomment to convert practical salinity to absolute salinity
  ! zt = gsw_ct_from_pt(S,T) ! Uncomment to convert potential temp to conservative temp

  call gsw_specvol_first_derivatives(zs, zt, zp, v_sa=dSV_dS, v_ct=dSV_dT)

end procedure calculate_specvol_derivs_elem_TEOS10
module procedure calculate_compress_elem_TEOS10
  real :: zs  ! Absolute salinity [g kg-1]
  real :: zt  ! Conservative temperature [degC]
  real :: zp  ! Pressure converted to decibars [dbar]
  zs = S
  zt = T
  zp = pressure * Pa2db      ! Convert pressure from Pascal to decibar
  ! The following conversions are unnecessary because the arguments are already the right variables.
  ! zs = gsw_sr_from_sp(S)   ! Uncomment to convert practical salinity to absolute salinity
  ! zt = gsw_ct_from_pt(S,T) ! Uncomment to convert potential temp to conservative temp

  rho = gsw_rho(zs, zt, zp)
  call gsw_rho_first_derivatives(zs, zt, zp, drho_dp=drho_dp)

end procedure calculate_compress_elem_TEOS10
module procedure EoS_fit_range_teos10
  if (present(T_min)) T_min = -6.0
  if (present(T_max)) T_max = 40.0
  if (present(S_min)) S_min =  0.0
  if (present(S_max)) S_max = 42.0
  if (present(p_min)) p_min = 0.0
  if (present(p_max)) p_max = 1.0e8

end procedure EoS_fit_range_teos10
end submodule MOM_EOS_TEOS10_s
