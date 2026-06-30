submodule (MOM_EOS_Roquet_Spv) MOM_EOS_Roquet_Spv_s
  implicit none
contains
module procedure spec_vol_elem_Roquet_SpV
  real :: zp     ! Pressure [Pa]
  real :: zt     ! Conservative temperature [degC]
  real :: zs     ! The square root of absolute salinity with an offset normalized
  real :: SV_00p ! A pressure-dependent but temperature and salinity independent contribution to
  real :: SV_TS  ! Specific volume without a pressure-dependent contribution [m3 kg-1]
  real :: SV_TS0 ! A contribution to specific volume from temperature and salinity anomalies at
  real :: SV_TS1 ! A temperature and salinity dependent specific volume contribution that is
  real :: SV_TS2 ! A temperature and salinity dependent specific volume contribution that is
  real :: SV_TS3 ! A temperature and salinity dependent specific volume contribution that is
  real :: SV_0S0 ! Salinity dependent specific volume at the surface pressure and zero temperature [m3 kg-1]
  zt = T
  zs = SQRT( ABS( S + rdeltaS ) * r1_S0 )  ! square root of normalized salinity plus an offset [nondim]
  zp = pressure

  ! The next two lines should be used if it is necessary to convert potential temperature and
  ! practical salinity to conservative temperature and absolute salinity.
  ! zt = gsw_ct_from_pt(S,T) ! Convert potential temp to conservative temp [degC]
  ! zs = SQRT( ABS( gsw_sr_from_sp(S) + rdeltaS ) * r1_S0 ) ! Convert S from practical to absolute salinity.

  SV_TS3 = SPV003 + (zs*SPV103 + zt*SPV013)
  SV_TS2 = SPV002 + (zs*(SPV102 +  zs*SPV202) &
                   + zt*(SPV012 + (zs*SPV112 + zt*SPV022)) )
  SV_TS1 = SPV001 + (zs*(SPV101 +  zs*(SPV201 +  zs*(SPV301 +  zs*SPV401))) &
                   + zt*(SPV011 + (zs*(SPV111 +  zs*(SPV211 +  zs*SPV311)) &
                                 + zt*(SPV021 + (zs*(SPV121 +  zs*SPV221) &
                                               + zt*(SPV031 + (zs*SPV131 + zt*SPV041)) )) )) )
  SV_TS0 = zt*(SPV010 &
             + (zs*(SPV110 +  zs*(SPV210 +  zs*(SPV310 +  zs*(SPV410 +  zs*SPV510)))) &
              + zt*(SPV020 + (zs*(SPV120 +  zs*(SPV220 +  zs*(SPV320 +  zs*SPV420))) &
                            + zt*(SPV030 + (zs*(SPV130 +  zs*(SPV230 +  zs*SPV330)) &
                                          + zt*(SPV040 + (zs*(SPV140 +  zs*SPV240) &
                                                        + zt*(SPV050 + (zs*SPV150 + zt*SPV060)) )) )) )) ) )

  SV_0S0 = SPV000 + zs*(SPV100 + zs*(SPV200 + zs*(SPV300 + zs*(SPV400 + zs*(SPV500 + zs*SPV600)))))

  SV_00p = zp*(V00 + zp*(V01 + zp*(V02 + zp*(V03 + zp*(V04 + zp*V05)))))

  SV_TS  = (SV_TS0 + SV_0S0) + zp*(SV_TS1 + zp*(SV_TS2 +  zp*SV_TS3))
  spec_vol_elem_Roquet_SpV = SV_TS + SV_00p  ! In situ specific volume [m3 kg-1]

end procedure spec_vol_elem_Roquet_SpV
module procedure spec_vol_anomaly_elem_Roquet_SpV
  real :: zp     ! Pressure [Pa]
  real :: zt     ! Conservative temperature [degC]
  real :: zs     ! The square root of absolute salinity with an offset normalized
  real :: SV_00p ! A pressure-dependent but temperature and salinity independent contribution to
  real :: SV_TS  ! Specific volume without a pressure-dependent contribution [m3 kg-1]
  real :: SV_TS0 ! A contribution to specific volume from temperature and salinity anomalies at
  real :: SV_TS1 ! A temperature and salinity dependent specific volume contribution that is
  real :: SV_TS2 ! A temperature and salinity dependent specific volume contribution that is
  real :: SV_TS3 ! A temperature and salinity dependent specific volume contribution that is
  real :: SV_0S0 ! Salinity dependent specific volume at the surface pressure and zero temperature [m3 kg-1]
  zt = T
  zs = SQRT( ABS( S + rdeltaS ) * r1_S0 )  ! square root of normalized salinity plus an offset [nondim]
  zp = pressure

  ! The next two lines should be used if it is necessary to convert potential temperature and
  ! practical salinity to conservative temperature and absolute salinity.
  ! zt = gsw_ct_from_pt(S,T) ! Convert potential temp to conservative temp [degC]
  ! zs = SQRT( ABS( gsw_sr_from_sp(S) + rdeltaS ) * r1_S0 ) ! Convert S from practical to absolute salinity.

  SV_TS3 = SPV003 + (zs*SPV103 + zt*SPV013)
  SV_TS2 = SPV002 + (zs*(SPV102 +  zs*SPV202) &
                   + zt*(SPV012 + (zs*SPV112 + zt*SPV022)) )
  SV_TS1 = SPV001 + (zs*(SPV101 +  zs*(SPV201 +  zs*(SPV301 +  zs*SPV401))) &
                   + zt*(SPV011 + (zs*(SPV111 +  zs*(SPV211 +  zs*SPV311)) &
                                 + zt*(SPV021 + (zs*(SPV121 +  zs*SPV221) &
                                               + zt*(SPV031 + (zs*SPV131 + zt*SPV041)) )) )) )
  SV_TS0 = zt*(SPV010 &
             + (zs*(SPV110 +  zs*(SPV210 +  zs*(SPV310 +  zs*(SPV410 +  zs*SPV510)))) &
              + zt*(SPV020 + (zs*(SPV120 +  zs*(SPV220 +  zs*(SPV320 +  zs*SPV420))) &
                            + zt*(SPV030 + (zs*(SPV130 +  zs*(SPV230 +  zs*SPV330)) &
                                          + zt*(SPV040 + (zs*(SPV140 +  zs*SPV240) &
                                                        + zt*(SPV050 + (zs*SPV150 + zt*SPV060)) )) )) )) ) )

  SV_0S0 = SPV000 + zs*(SPV100 + zs*(SPV200 + zs*(SPV300 + zs*(SPV400 + zs*(SPV500 + zs*SPV600)))))

  SV_00p = zp*(V00 + zp*(V01 + zp*(V02 + zp*(V03 + zp*(V04 + zp*V05)))))

  SV_0S0 = SV_0S0 - spv_ref

  SV_TS  = (SV_TS0 + SV_0S0) + zp*(SV_TS1 + zp*(SV_TS2 +  zp*SV_TS3))
  spec_vol_anomaly_elem_Roquet_SpV = SV_TS + SV_00p  ! In situ specific volume [m3 kg-1]

end procedure spec_vol_anomaly_elem_Roquet_SpV
module procedure density_elem_Roquet_SpV
  real :: spv ! The specific volume [m3 kg-1]
  spv = spec_vol_elem_Roquet_SpV(this, T, S, pressure)
  density_elem_Roquet_SpV = 1.0 / spv  ! In situ density [kg m-3]

end procedure density_elem_Roquet_SpV
module procedure density_anomaly_elem_Roquet_SpV
  real :: spv ! The specific volume [m3 kg-1]
  spv = spec_vol_anomaly_elem_Roquet_SpV(this, T, S, pressure, spv_ref=1.0/rho_ref)
  density_anomaly_elem_Roquet_SpV = -rho_ref**2*spv / (rho_ref*spv + 1.0)  ! In situ density [kg m-3]

end procedure density_anomaly_elem_Roquet_SpV
module procedure calculate_specvol_derivs_elem_Roquet_SpV
  real :: zp      ! Pressure [Pa]
  real :: zt      ! Conservative temperature [degC]
  real :: zs      ! The square root of absolute salinity with an offset normalized
  real :: dSVdzt0 ! A contribution to the partial derivative of specific volume with temperature
  real :: dSVdzt1 ! A contribution to the partial derivative of specific volume with temperature
  real :: dSVdzt2 ! A contribution to the partial derivative of specific volume with temperature
  real :: dSVdzt3 ! A contribution to the partial derivative of specific volume with temperature
  real :: dSVdzs0 ! A contribution to the partial derivative of specific volume with
  real :: dSVdzs1 ! A contribution to the partial derivative of specific volume with
  real :: dSVdzs2 ! A contribution to the partial derivative of specific volume with
  real :: dSVdzs3 ! A contribution to the partial derivative of specific volume with
  zt = T
  zs = SQRT( ABS( S + rdeltaS ) * r1_S0 )  ! square root of normalized salinity plus an offset [nondim]
  zp = pressure

  ! The next two lines should be used if it is necessary to convert potential temperature and
  ! practical salinity to conservative temperature and absolute salinity.
  ! zt = gsw_ct_from_pt(S,T) ! Convert potential temp to conservative temp [degC]
  ! zs = SQRT( ABS( gsw_sr_from_sp(S) + rdeltaS ) * r1_S0 ) ! Convert S from practical to absolute salinity.

  ! Find the partial derivative of specific volume with temperature
  dSVdzt3 = ALP003
  dSVdzt2 = ALP002 + (zs*ALP102 + zt*ALP012)
  dSVdzt1 = ALP001 + (zs*(ALP101 + zs*(ALP201 + zs*ALP301)) &
                    + zt*(ALP011 + (zs*(ALP111 + zs*ALP211) &
                                  + zt*(ALP021 + (zs*ALP121 + zt*ALP031)) )) )
  dSVdzt0 = ALP000 + (zs*(ALP100 +  zs*(ALP200 +  zs*(ALP300 + zs*(ALP400 + zs*ALP500)))) &
                    + zt*(ALP010 + (zs*(ALP110 +  zs*(ALP210 + zs*(ALP310 + zs*ALP410))) &
                                  + zt*(ALP020 + (zs*(ALP120 + zs*(ALP220 + zs*ALP320)) &
                                                + zt*(ALP030 + (zt*(ALP040 + (zs*ALP140 + zt*ALP050)) &
                                                              + zs*(ALP130 + zs*ALP230) )) )) )) )

  dSV_dT = dSVdzt0 + zp*(dSVdzt1 + zp*(dSVdzt2 + zp*dSVdzt3))

  ! Find the partial derivative of specific volume with salinity
  dSVdzs3 = BET003
  dSVdzs2 = BET002 + (zs*BET102 + zt*BET012)
  dSVdzs1 = BET001 + (zs*(BET101 + zs*(BET201 + zs*BET301)) &
                    + zt*(BET011 + (zs*(BET111 + zs*BET211) &
                                  + zt*(BET021 + (zs*BET121 + zt*BET031)) )) )
  dSVdzs0 = BET000 + (zs*(BET100 + zs*(BET200 + zs*(BET300 + zs*(BET400 + zs*BET500)))) &
                    + zt*(BET010 + (zs*(BET110 + zs*(BET210 + zs*(BET310 + zs*BET410))) &
                                  + zt*(BET020 + (zs*(BET120 + zs*(BET220 + zs*BET320)) &
                                                + zt*(BET030 + (zt*(BET040 + (zs*BET140 + zt*BET050)) &
                                                              + zs*(BET130 + zs*BET230) )) )) )) )

  ! The division by zs here is because zs = sqrt(S + S0), so dSV_dS = dzs_dS * dSV_dzs = (0.5 / zs) * dSV_dzs
  dSV_dS = (dSVdzs0 + zp*(dSVdzs1 + zp*(dSVdzs2 + zp * dSVdzs3))) / zs

end procedure calculate_specvol_derivs_elem_Roquet_SpV
module procedure calculate_density_derivs_elem_Roquet_SpV
  real :: dSV_dT   ! The partial derivative of specific volume with
  real :: dSV_dS   ! The partial derivative of specific volume with
  real :: specvol  ! The specific volume [m3 kg-1]
  real :: rho  ! The in situ density [kg m-3]
  call this%calculate_specvol_derivs_elem(T, S, pressure, dSV_dT, dSV_dS)

  specvol = this%spec_vol_elem(T, S, pressure)
  rho = 1.0 / specvol
  drho_dT = -dSv_dT * rho**2
  drho_dS = -dSv_dS * rho**2

end procedure calculate_density_derivs_elem_Roquet_SpV
module procedure calculate_compress_elem_Roquet_SpV
  real :: zp     ! Pressure [Pa]
  real :: zt     ! Conservative temperature [degC]
  real :: zs     ! The square root of absolute salinity with an offset normalized
  real :: dSV_00p_dp ! Derivative of the pressure-dependent reference specific volume profile with
  real :: dSV_TS_dp  ! Derivative of the specific volume anomaly from the reference profile with
  real :: SV_00p ! A pressure-dependent but temperature and salinity independent contribution to
  real :: SV_TS  ! specific volume without a pressure-dependent contribution [m3 kg-1]
  real :: SV_TS0 ! A contribution to specific volume from temperature and salinity anomalies at
  real :: SV_TS1 ! A temperature and salinity dependent specific volume contribution that is
  real :: SV_TS2 ! A temperature and salinity dependent specific volume contribution that is
  real :: SV_TS3 ! A temperature and salinity dependent specific volume contribution that is
  real :: SV_0S0 ! Salinity dependent specific volume at the surface pressure and zero temperature [m3 kg-1]
  real :: dSpecVol_dp ! The partial derivative of specific volume with pressure [m3 kg-1 Pa-1]
  zt = T
  zs = SQRT( ABS( S + rdeltaS ) * r1_S0 )  ! square root of normalized salinity plus an offset [nondim]
  zp = pressure

  ! The next two lines should be used if it is necessary to convert potential temperature and
  ! practical salinity to conservative temperature and absolute salinity.
  ! zt = gsw_ct_from_pt(S,T) ! Convert potential temp to conservative temp [degC]
  ! zs = SQRT( ABS( gsw_sr_from_sp(S) + rdeltaS ) * r1_S0 ) ! Convert S from practical to absolute salinity.

  SV_TS3 = SPV003 + (zs*SPV103 + zt*SPV013)
  SV_TS2 = SPV002 + (zs*(SPV102 +  zs*SPV202) &
                   + zt*(SPV012 + (zs*SPV112 + zt*SPV022)) )
  SV_TS1 = SPV001 + (zs*(SPV101 +  zs*(SPV201 +  zs*(SPV301 +  zs*SPV401))) &
                   + zt*(SPV011 + (zs*(SPV111 +  zs*(SPV211 +  zs*SPV311)) &
                                 + zt*(SPV021 + (zs*(SPV121 +  zs*SPV221) &
                                               + zt*(SPV031 + (zs*SPV131 + zt*SPV041)) )) )) )

  SV_TS0 = zt*(SPV010 &
             + (zs*(SPV110 +  zs*(SPV210 +  zs*(SPV310 +  zs*(SPV410 +  zs*SPV510)))) &
              + zt*(SPV020 + (zs*(SPV120 +  zs*(SPV220 +  zs*(SPV320 +  zs*SPV420))) &
                            + zt*(SPV030 + (zs*(SPV130 +  zs*(SPV230 +  zs*SPV330)) &
                                          + zt*(SPV040 + (zs*(SPV140 +  zs*SPV240) &
                                                        + zt*(SPV050 + (zs*SPV150 + zt*SPV060)) )) )) )) ) )

  SV_0S0 = SPV000 + zs*(SPV100 + zs*(SPV200 + zs*(SPV300 + zs*(SPV400 + zs*(SPV500 + zs*SPV600)))))

  SV_00p = zp*(V00 + zp*(V01 + zp*(V02 + zp*(V03 + zp*(V04 + zp*V05)))))

  SV_TS  = (SV_TS0 + SV_0S0) + zp*(SV_TS1 + zp*(SV_TS2 +  zp*SV_TS3))
  ! specvol = SV_TS + SV_00p ! In situ specific volume [m3 kg-1]
  rho = 1.0 / (SV_TS + SV_00p) ! In situ density [kg m-3]

  dSV_00p_dp = V00 + zp*(2.*V01 + zp*(3.*V02 + zp*(4.*V03 + zp*(5.*V04 + zp*(6.*V05)))))
  dSV_TS_dp  = SV_TS1 + zp*(2.*SV_TS2 + zp*(3.*SV_TS3))
  dSpecVol_dp = dSV_TS_dp + dSV_00p_dp  !  [m3 kg-1 Pa-1]
  drho_dp = -dSpecVol_dp * rho**2 ! Compressibility [s2 m-2]

end procedure calculate_compress_elem_Roquet_SpV
module procedure calc_spec_vol_second_derivs_elem_Roquet_SpV
  real :: zp      ! Pressure [Pa]
  real :: zt      ! Conservative temperature [degC]
  real :: zs      ! The square root of absolute salinity with an offset normalized
  real :: I_s     ! The inverse of zs [nondim]
  real :: d2SV_p0 ! A contribution to one of the second derivatives that is independent of pressure [various]
  real :: d2SV_p1 ! A contribution to one of the second derivatives that is proportional to pressure [various]
  real :: d2SV_p2 ! A contribution to one of the second derivatives that is proportional to pressure**2 [various]
  real :: d2SV_p3 ! A contribution to one of the second derivatives that is proportional to pressure**3 [various]
  zt = T
  zs = SQRT( ABS( S + rdeltaS ) * r1_S0 )  ! square root of normalized salinity plus an offset [nondim]
  zp = P

  ! The next two lines should be used if it is necessary to convert potential temperature and
  ! practical salinity to conservative temperature and absolute salinity.
  ! zt = gsw_ct_from_pt(S,T) ! Convert potential temp to conservative temp [degC]
  ! zs = SQRT( ABS( gsw_sr_from_sp(S) + rdeltaS ) * r1_S0 )  ! Convert S from practical to absolute salinity.

  I_s = 1.0 / zs

  ! Find dSV_ds_ds
  d2SV_p3 = -SPV103*I_s**2
  d2SV_p2 = -(SPV102 + zt*SPV112)*I_s**2
  d2SV_p1 = (3.*SPV301 + (zt*(3.*SPV311) + zs*(8.*SPV401))) &
            - ( SPV101 + zt*(SPV111 + zt*(SPV121 + zt*SPV131)) )*I_s**2
  d2SV_p0 = (3.*SPV300 + (zs*(8.*SPV400 + zs*(15.*SPV500 + zs*(24.*SPV600))) &
                        + zt*(3.*SPV310 + (zs*(8.*SPV410 + zs*(15.*SPV510)) &
                                         + zt*(3.*SPV320 + (zs*(8.*SPV420) + zt*(3.*SPV330))) )) )) &
            - (SPV100 + zt*(SPV110 + zt*(SPV120 + zt*(SPV130 + zt*(SPV140 + zt*SPV150)))) )*I_s**2
  dSV_dS_dS = (0.5*r1_S0)**2 * ((d2SV_p0 + zp*(d2SV_p1 + zp*(d2SV_p2 + zp*d2SV_p3))) * I_s)

  ! Find dSV_ds_dt
  d2SV_p2 = SPV112
  d2SV_p1 = SPV111 + (zs*(2.*SPV211 +  zs*(3.*SPV311)) &
                    + zt*(2.*SPV121 + (zs*(4.*SPV221) + zt*(3.*SPV131))) )
  d2SV_p0 = SPV110 + (zs*(2.*SPV210 +  zs*(3.*SPV310 +  zs*(4.*SPV410 +  zs*(5.*SPV510)))) &
                    + zt*(2.*SPV120 + (zs*(4.*SPV220 +  zs*(6.*SPV320 +  zs*(8.*SPV420))) &
                                     + zt*(3.*SPV130 + (zs*(6.*SPV230 +  zs*(9.*SPV330)) &
                                                      + zt*(4.*SPV140 + (zs*(8.*SPV240) &
                                                                       + zt*(5.*SPV150))) )) )) )
  dSV_ds_dt = (0.5*r1_S0) * ((d2SV_p0 + zp*(d2SV_p1 + zp*d2SV_p2)) * I_s)

  ! Find dSV_dt_dt
  d2SV_p2 = 2.*SPV022
  d2SV_p1 = 2.*SPV021 + (zs*(2.*SPV121 +  zs*(2.*SPV221)) &
                       + zt*(6.*SPV031 + (zs*(6.*SPV131) + zt*(12.*SPV041))) )
  d2SV_p0 = 2.*SPV020 + (zs*(2.*SPV120 +  zs*( 2.*SPV220 +  zs*( 2.*SPV320 + zs * (2.*SPV420)))) &
                       + zt*(6.*SPV030 + (zs*( 6.*SPV130 +  zs*( 6.*SPV230 + zs * (6.*SPV330))) &
                                        + zt*(12.*SPV040 + (zs*(12.*SPV140 + zs *(12.*SPV240)) &
                                                          + zt*(20.*SPV050 + (zs*(20.*SPV150) &
                                                                            + zt*(30.*SPV060) )) )) )) )
  dSV_dt_dt = d2SV_p0 + zp*(d2SV_p1 + zp*d2SV_p2)

  ! Find dSV_ds_dp
  d2SV_p2 = 3.*SPV103
  d2SV_p1 = 2.*SPV102 + (zs*(4.*SPV202) + zt*(2.*SPV112))
  d2SV_p0 = SPV101 + (zs*(2.*SPV201 + zs*(3.*SPV301 +  zs*(4.*SPV401))) &
                    + zt*(SPV111 +   (zs*(2.*SPV211 +  zs*(3.*SPV311)) &
                                    + zt*(   SPV121 + (zs*(2.*SPV221) + zt*SPV131)) )) )
  dSV_ds_dp =  ((d2SV_p0 + zp*(d2SV_p1 + zp*d2SV_p2)) * I_s) * (0.5*r1_S0)

  ! Find dSV_dt_dp
  d2SV_p2 = 3.*SPV013
  d2SV_p1 = 2.*SPV012 + (zs*(2.*SPV112) + zt*(4.*SPV022))
  d2SV_p0 = SPV011 + (zs*(SPV111     + zs*(   SPV211 +  zs*    SPV311)) &
                    + zt*(2.*SPV021 + (zs*(2.*SPV121 +  zs*(2.*SPV221)) &
                                     + zt*(3.*SPV031 + (zs*(3.*SPV131) + zt*(4.*SPV041))) )) )
  dSV_dt_dp =  d2SV_p0 + zp*(d2SV_p1 + zp*d2SV_p2)

end procedure calc_spec_vol_second_derivs_elem_Roquet_SpV
module procedure calculate_density_second_derivs_elem_Roquet_SpV
  real :: rho       ! The in situ density [kg m-3]
  real :: drho_dp   ! The partial derivative of density with pressure
  real :: dSV_dT    ! The partial derivative of specific volume with
  real :: dSV_dS    ! The partial derivative of specific volume with
  real :: dSV_ds_ds ! Second derivative of specific volume with respect
  real :: dSV_ds_dt ! Second derivative of specific volume with respect
  real :: dSV_dt_dt ! Second derivative of specific volume with respect
  real :: dSV_ds_dp ! Second derivative of specific volume with respect to pressure
  real :: dSV_dt_dp ! Second derivative of specific volume with respect to pressure
  call calc_spec_vol_second_derivs_elem_Roquet_SpV(T, S, pressure, &
                 dSV_ds_ds, dSV_ds_dt, dSV_dt_dt, dSV_ds_dp, dSV_dt_dp)
  call this%calculate_specvol_derivs_elem(T, S, pressure, dSV_dT, dSV_dS)
  call this%calculate_compress_elem(T, S, pressure, rho, drho_dp)

  ! Find drho_ds_ds
  drho_dS_dS = rho**2 * (2.0*rho*dSV_dS**2 - dSV_dS_dS)

  ! Find drho_ds_dt
  drho_ds_dt = rho**2 * (2.0*rho*(dSV_dT*dSV_dS) - dSV_dS_dT)

  ! Find drho_dt_dt
  drho_dT_dT = rho**2 * (2.0*rho*dSV_dT**2 - dSV_dT_dT)

  ! Find drho_ds_dp
  drho_ds_dp =  -rho * (2.0*dSV_dS * drho_dp + rho * dSV_dS_dp)

  ! Find drho_dt_dp
  drho_dt_dp =  -rho * (2.0*dSV_dT * drho_dp + rho * dSV_dT_dp)

end procedure calculate_density_second_derivs_elem_Roquet_SpV
module procedure EoS_fit_range_Roquet_SpV
  if (present(T_min)) T_min = -6.0
  if (present(T_max)) T_max = 40.0
  if (present(S_min)) S_min =  0.0
  if (present(S_max)) S_max = 42.0
  if (present(p_min)) p_min = 0.0
  if (present(p_max)) p_max = 1.0e8

end procedure EoS_fit_range_Roquet_SpV
module procedure calculate_density_array_Roquet_SpV
  integer :: j
  if (present(rho_ref)) then
    do j = start, start+npts-1
      rho(j) = density_anomaly_elem_Roquet_SpV(this, T(j), S(j), pressure(j), rho_ref)
    enddo
  else
    do j = start, start+npts-1
      rho(j) = density_elem_Roquet_SpV(this, T(j), S(j), pressure(j))
    enddo
  endif

end procedure calculate_density_array_Roquet_SpV
module procedure calculate_spec_vol_array_Roquet_SpV
  integer :: j
  if (present(spv_ref)) then
    do j = start, start+npts-1
      specvol(j) = spec_vol_anomaly_elem_Roquet_SpV(this, T(j), S(j), pressure(j), spv_ref)
    enddo
  else
    do j = start, start+npts-1
      specvol(j) = spec_vol_elem_Roquet_SpV(this, T(j), S(j), pressure(j) )
    enddo
  endif

end procedure calculate_spec_vol_array_Roquet_SpV
end submodule MOM_EOS_Roquet_Spv_s
