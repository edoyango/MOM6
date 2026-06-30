submodule (MOM_EOS_Roquet_rho) MOM_EOS_Roquet_rho_s
  implicit none
contains
module procedure density_elem_Roquet_rho
  real :: zp     ! Pressure [Pa]
  real :: zt     ! Conservative temperature [degC]
  real :: zs     ! The square root of absolute salinity with an offset normalized
  real :: rho00p ! A pressure-dependent but temperature and salinity independent contribution to
  real :: rhoTS  ! Density without a pressure-dependent contribution [kg m-3]
  real :: rhoTS0 ! A contribution to density from temperature and salinity anomalies at the
  real :: rhoTS1 ! A density contribution proportional to pressure [kg m-3 Pa-1]
  real :: rhoTS2 ! A density contribution proportional to pressure**2 [kg m-3 Pa-2]
  real :: rhoTS3 ! A density contribution proportional to pressure**3 [kg m-3 Pa-3]
  real :: rho0S0 ! Salinity dependent density at the surface pressure and zero temperature [kg m-3]
  zt = T
  zs = SQRT( ABS( S + rdeltaS ) * r1_S0 )  ! square root of normalized salinity plus an offset [nondim]
  zp = pressure

  ! The next two lines should be used if it is necessary to convert potential temperature and
  ! practical salinity to conservative temperature and absolute salinity.
  ! zt = gsw_ct_from_pt(S,T) ! Convert potential temp to conservative temp [degC]
  ! zs = SQRT( ABS( gsw_sr_from_sp(S) + rdeltaS ) * r1_S0 ) ! Convert S from practical to absolute salinity.

  rhoTS3 = EOS003 + (zs*EOS103 + zt*EOS013)
  rhoTS2 = EOS002 + (zs*(EOS102 +  zs*EOS202) &
                   + zt*(EOS012 + (zs*EOS112 + zt*EOS022)) )
  rhoTS1 = EOS001 + (zs*(EOS101 +  zs*(EOS201 +  zs*(EOS301 +  zs*EOS401))) &
                   + zt*(EOS011 + (zs*(EOS111 +  zs*(EOS211 +  zs*EOS311)) &
                                 + zt*(EOS021 + (zs*(EOS121 +  zs*EOS221) &
                                               + zt*(EOS031 + (zs*EOS131 + zt*EOS041)) )) )) )
  rhoTS0 = zt*(EOS010 &
             + (zs*(EOS110 +  zs*(EOS210 +  zs*(EOS310 +  zs*(EOS410 +  zs*EOS510)))) &
              + zt*(EOS020 + (zs*(EOS120 +  zs*(EOS220 +  zs*(EOS320 +  zs*EOS420))) &
                            + zt*(EOS030 + (zs*(EOS130 +  zs*(EOS230 +  zs*EOS330)) &
                                          + zt*(EOS040 + (zs*(EOS140 +  zs*EOS240) &
                                                        + zt*(EOS050 + (zs*EOS150 + zt*EOS060)) )) )) )) ) )

  rho0S0 = EOS000 + zs*(EOS100 + zs*(EOS200 + zs*(EOS300 + zs*(EOS400 + zs*(EOS500 + zs*EOS600)))))

  rho00p = zp*(R00 + zp*(R01 + zp*(R02 + zp*(R03 + zp*(R04 + zp*R05)))))

  rhoTS  = (rhoTS0 + rho0S0) + zp*(rhoTS1 + zp*(rhoTS2 +  zp*rhoTS3))
  density_elem_Roquet_rho = rhoTS + rho00p  ! In situ density [kg m-3]

end procedure density_elem_Roquet_rho
module procedure density_anomaly_elem_Roquet_rho
  real :: zp     ! Pressure [Pa]
  real :: zt     ! Conservative temperature [degC]
  real :: zs     ! The square root of absolute salinity with an offset normalized
  real :: rho00p ! A pressure-dependent but temperature and salinity independent contribution to
  real :: rhoTS  ! Density without a pressure-dependent contribution [kg m-3]
  real :: rhoTS0 ! A contribution to density from temperature and salinity anomalies at the
  real :: rhoTS1 ! A density contribution proportional to pressure [kg m-3 Pa-1]
  real :: rhoTS2 ! A density contribution proportional to pressure**2 [kg m-3 Pa-2]
  real :: rhoTS3 ! A density contribution proportional to pressure**3 [kg m-3 Pa-3]
  real :: rho0S0 ! Salinity dependent density at the surface pressure and zero temperature [kg m-3]
  zt = T
  zs = SQRT( ABS( S + rdeltaS ) * r1_S0 )  ! square root of normalized salinity plus an offset [nondim]
  zp = pressure

  ! The next two lines should be used if it is necessary to convert potential temperature and
  ! practical salinity to conservative temperature and absolute salinity.
  ! zt = gsw_ct_from_pt(S,T) ! Convert potential temp to conservative temp [degC]
  ! zs = SQRT( ABS( gsw_sr_from_sp(S) + rdeltaS ) * r1_S0 ) ! Convert S from practical to absolute salinity.

  rhoTS3 = EOS003 + (zs*EOS103 + zt*EOS013)
  rhoTS2 = EOS002 + (zs*(EOS102 +  zs*EOS202) &
                   + zt*(EOS012 + (zs*EOS112 + zt*EOS022)) )
  rhoTS1 = EOS001 + (zs*(EOS101 +  zs*(EOS201 +  zs*(EOS301 +  zs*EOS401))) &
                   + zt*(EOS011 + (zs*(EOS111 +  zs*(EOS211 +  zs*EOS311)) &
                                 + zt*(EOS021 + (zs*(EOS121 +  zs*EOS221) &
                                               + zt*(EOS031 + (zs*EOS131 + zt*EOS041)) )) )) )
  rhoTS0 = zt*(EOS010 &
             + (zs*(EOS110 +  zs*(EOS210 +  zs*(EOS310 +  zs*(EOS410 +  zs*EOS510)))) &
              + zt*(EOS020 + (zs*(EOS120 +  zs*(EOS220 +  zs*(EOS320 +  zs*EOS420))) &
                            + zt*(EOS030 + (zs*(EOS130 +  zs*(EOS230 +  zs*EOS330)) &
                                          + zt*(EOS040 + (zs*(EOS140 +  zs*EOS240) &
                                                        + zt*(EOS050 + (zs*EOS150 + zt*EOS060)) )) )) )) ) )

  rho0S0 = EOS000 + zs*(EOS100 + zs*(EOS200 + zs*(EOS300 + zs*(EOS400 + zs*(EOS500 + zs*EOS600)))))

  rho00p = zp*(R00 + zp*(R01 + zp*(R02 + zp*(R03 + zp*(R04 + zp*R05)))))

  rho0S0 = rho0S0 - rho_ref

  rhoTS  = (rhoTS0 + rho0S0) + zp*(rhoTS1 + zp*(rhoTS2 +  zp*rhoTS3))
  density_anomaly_elem_Roquet_rho = rhoTS + rho00p  ! In situ density [kg m-3]

end procedure density_anomaly_elem_Roquet_rho
module procedure spec_vol_elem_Roquet_rho
  spec_vol_elem_Roquet_rho = 1. / density_elem_Roquet_rho(this, T, S, pressure)

end procedure spec_vol_elem_Roquet_rho
module procedure spec_vol_anomaly_elem_Roquet_rho
  spec_vol_anomaly_elem_Roquet_rho = 1. / density_elem_Roquet_rho(this, T, S, pressure)
  spec_vol_anomaly_elem_Roquet_rho = spec_vol_anomaly_elem_Roquet_rho - spv_ref

end procedure spec_vol_anomaly_elem_Roquet_rho
module procedure calculate_density_derivs_elem_Roquet_rho
  real :: zp      ! Pressure [Pa]
  real :: zt      ! Conservative temperature [degC]
  real :: zs      ! The square root of absolute salinity with an offset normalized
  real :: dRdzt0  ! A contribution to the partial derivative of density with temperature [kg m-3 degC-1]
  real :: dRdzt1  ! A contribution to the partial derivative of density with temperature [kg m-3 degC-1 Pa-1]
  real :: dRdzt2  ! A contribution to the partial derivative of density with temperature [kg m-3 degC-1 Pa-2]
  real :: dRdzt3  ! A contribution to the partial derivative of density with temperature [kg m-3 degC-1 Pa-3]
  real :: dRdzs0  ! A contribution to the partial derivative of density with
  real :: dRdzs1  ! A contribution to the partial derivative of density with
  real :: dRdzs2  ! A contribution to the partial derivative of density with
  real :: dRdzs3  ! A contribution to the partial derivative of density with
  zt = T
  zs = SQRT( ABS( S + rdeltaS ) * r1_S0 )  ! square root of normalized salinity plus an offset [nondim]
  zp = pressure

  ! The next two lines should be used if it is necessary to convert potential temperature and
  ! practical salinity to conservative temperature and absolute salinity.
  ! zt = gsw_ct_from_pt(S,T) ! Convert potential temp to conservative temp [degC]
  ! zs = SQRT( ABS( gsw_sr_from_sp(S) + rdeltaS ) * r1_S0 ) ! Convert S from practical to absolute salinity.

  ! Find the partial derivative of density with temperature
  dRdzt3 = ALP003
  dRdzt2 = ALP002 + (zs*ALP102 + zt*ALP012)
  dRdzt1 = ALP001 + (zs*(ALP101 + zs*(ALP201 + zs*ALP301)) &
                   + zt*(ALP011 + (zs*(ALP111 + zs*ALP211) &
                                 + zt*(ALP021 + (zs*ALP121 + zt*ALP031)) )) )
  dRdzt0 = ALP000 + (zs*(ALP100 +  zs*(ALP200 +  zs*(ALP300 + zs*(ALP400 + zs*ALP500)))) &
                   + zt*(ALP010 + (zs*(ALP110 +  zs*(ALP210 + zs*(ALP310 + zs*ALP410))) &
                                 + zt*(ALP020 + (zs*(ALP120 + zs*(ALP220 + zs*ALP320)) &
                                               + zt*(ALP030 + (zt*(ALP040 + (zs*ALP140 + zt*ALP050)) &
                                                             + zs*(ALP130 + zs*ALP230) )) )) )) )

  drho_dT = dRdzt0 + zp*(dRdzt1 + zp*(dRdzt2 + zp*dRdzt3))

  ! Find the partial derivative of density with salinity
  dRdzs3 = BET003
  dRdzs2 = BET002 + (zs*BET102 + zt*BET012)
  dRdzs1 = BET001 + (zs*(BET101 + zs*(BET201 + zs*BET301)) &
                   + zt*(BET011 + (zs*(BET111 + zs*BET211) &
                                 + zt*(BET021 + (zs*BET121 + zt*BET031)) )) )
  dRdzs0 = BET000 + (zs*(BET100 + zs*(BET200 + zs*(BET300 + zs*(BET400 + zs*BET500)))) &
                   + zt*(BET010 + (zs*(BET110 + zs*(BET210 + zs*(BET310 + zs*BET410))) &
                                 + zt*(BET020 + (zs*(BET120 + zs*(BET220 + zs*BET320)) &
                                               + zt*(BET030 + (zt*(BET040 + (zs*BET140 + zt*BET050)) &
                                                             + zs*(BET130 + zs*BET230) )) )) )) )

  ! The division by zs here is because zs = sqrt(S + S0), so drho_dS = dzs_dS * drho_dzs = (0.5 / zs) * drho_dzs
  drho_dS = (dRdzs0 + zp*(dRdzs1 + zp*(dRdzs2 + zp * dRdzs3))) / zs

end procedure calculate_density_derivs_elem_Roquet_rho
module procedure calculate_density_second_derivs_elem_Roquet_rho
  real :: zp     ! Pressure [Pa]
  real :: zt     ! Conservative temperature [degC]
  real :: zs     ! The square root of absolute salinity with an offset normalized
  real :: I_s    ! The inverse of zs [nondim]
  real :: d2R_p0 ! A contribution to one of the second derivatives that is independent of pressure [various]
  real :: d2R_p1 ! A contribution to one of the second derivatives that is proportional to pressure [various]
  real :: d2R_p2 ! A contribution to one of the second derivatives that is proportional to pressure**2 [various]
  real :: d2R_p3 ! A contribution to one of the second derivatives that is proportional to pressure**3 [various]
  zt = T
  zs = SQRT( ABS( S + rdeltaS ) * r1_S0 )  ! square root of normalized salinity plus an offset [nondim]
  zp = pressure

  ! The next two lines should be used if it is necessary to convert potential temperature and
  ! practical salinity to conservative temperature and absolute salinity.
  ! zt = gsw_ct_from_pt(S,T) ! Convert potential temp to conservative temp [degC]
  ! zs = SQRT( ABS( gsw_sr_from_sp(S) + rdeltaS ) * r1_S0 )  ! Convert S from practical to absolute salinity.

  I_s = 1.0 / zs

  ! Find drho_ds_ds
  d2R_p3 = -EOS103*I_s**2
  d2R_p2 = -(EOS102 + zt*EOS112)*I_s**2
  d2R_p1 = (3.*EOS301 + (zt*(3.*EOS311) + zs*(8.*EOS401))) &
           - ( EOS101 + zt*(EOS111 + zt*(EOS121 + zt*EOS131)) )*I_s**2
  d2R_p0 = (3.*EOS300 + (zs*(8.*EOS400 + zs*(15.*EOS500 + zs*(24.*EOS600))) &
                       + zt*(3.*EOS310 + (zs*(8.*EOS410 + zs*(15.*EOS510)) &
                                        + zt*(3.*EOS320 + (zs*(8.*EOS420) + zt*(3.*EOS330))) )) )) &
           - (EOS100 + zt*(EOS110 + zt*(EOS120 + zt*(EOS130 + zt*(EOS140 + zt*EOS150)))) )*I_s**2
  drho_dS_dS = (0.5*r1_S0)**2 * ((d2R_p0 + zp*(d2R_p1 + zp*(d2R_p2 + zp*d2R_p3))) * I_s)

  ! Find drho_ds_dt
  d2R_p2 = EOS112
  d2R_p1 = EOS111 + (zs*(2.*EOS211 +  zs*(3.*EOS311)) &
                   + zt*(2.*EOS121 + (zs*(4.*EOS221) + zt*(3.*EOS131))) )
  d2R_p0 = EOS110 + (zs*(2.*EOS210 +  zs*(3.*EOS310 +  zs*(4.*EOS410 +  zs*(5.*EOS510)))) &
                   + zt*(2.*EOS120 + (zs*(4.*EOS220 +  zs*(6.*EOS320 +  zs*(8.*EOS420))) &
                                    + zt*(3.*EOS130 + (zs*(6.*EOS230 +  zs*(9.*EOS330)) &
                                                     + zt*(4.*EOS140 + (zs*(8.*EOS240) &
                                                                      + zt*(5.*EOS150))) )) )) )
  drho_ds_dt = (0.5*r1_S0) * ((d2R_p0 + zp*(d2R_p1 + zp*d2R_p2)) * I_s)

  ! Find drho_dt_dt
  d2R_p2 = 2.*EOS022
  d2R_p1 = 2.*EOS021 + (zs*(2.*EOS121 +  zs*(2.*EOS221)) &
                      + zt*(6.*EOS031 + (zs*(6.*EOS131) + zt*(12.*EOS041))) )
  d2R_p0 = 2.*EOS020 + (zs*(2.*EOS120 +  zs*( 2.*EOS220 +  zs*( 2.*EOS320 + zs * (2.*EOS420)))) &
                      + zt*(6.*EOS030 + (zs*( 6.*EOS130 +  zs*( 6.*EOS230 + zs * (6.*EOS330))) &
                                       + zt*(12.*EOS040 + (zs*(12.*EOS140 + zs *(12.*EOS240)) &
                                                         + zt*(20.*EOS050 + (zs*(20.*EOS150) &
                                                                           + zt*(30.*EOS060) )) )) )) )
  drho_dt_dt = (d2R_p0 + zp*(d2R_p1 + zp*d2R_p2))

  ! Find drho_ds_dp
  d2R_p2 = 3.*EOS103
  d2R_p1 = 2.*EOS102 + (zs*(4.*EOS202) + zt*(2.*EOS112))
  d2R_p0 = EOS101 + (zs*(2.*EOS201 + zs*(3.*EOS301 +  zs*(4.*EOS401))) &
                   + zt*(EOS111 +   (zs*(2.*EOS211 +  zs*(3.*EOS311)) &
                                   + zt*(   EOS121 + (zs*(2.*EOS221) + zt*EOS131)) )) )
  drho_ds_dp =  ((d2R_p0 + zp*(d2R_p1 + zp*d2R_p2)) * I_s) * (0.5*r1_S0)

  ! Find drho_dt_dp
  d2R_p2 = 3.*EOS013
  d2R_p1 = 2.*EOS012 + (zs*(2.*EOS112) + zt*(4.*EOS022))
  d2R_p0 = EOS011 + (zs*(EOS111     + zs*(   EOS211 +  zs*    EOS311)) &
                   + zt*(2.*EOS021 + (zs*(2.*EOS121 +  zs*(2.*EOS221)) &
                                    + zt*(3.*EOS031 + (zs*(3.*EOS131) + zt*(4.*EOS041))) )) )
  drho_dt_dp =  (d2R_p0 + zp*(d2R_p1 + zp*d2R_p2))

end procedure calculate_density_second_derivs_elem_Roquet_rho
module procedure calculate_specvol_derivs_elem_Roquet_rho
  real :: rho     ! In situ density [kg m-3]
  real :: dRho_dT ! Derivative of density with temperature [kg m-3 degC-1]
  real :: dRho_dS ! Derivative of density with salinity [kg m-3 ppt-1]
  call this%calculate_density_derivs_elem(T, S, pressure, drho_dT, drho_dS)
  rho = this%density_elem(T, S, pressure)
  dSV_dT = -dRho_DT/(rho**2)
  dSV_dS = -dRho_DS/(rho**2)

end procedure calculate_specvol_derivs_elem_Roquet_rho
module procedure calculate_compress_elem_Roquet_rho
  real :: zp     ! Pressure [Pa]
  real :: zt     ! Conservative temperature [degC]
  real :: zs     ! The square root of absolute salinity with an offset normalized
  real :: drho00p_dp ! Derivative of the pressure-dependent reference density profile with pressure [kg m-3 Pa-1]
  real :: drhoTS_dp  ! Derivative of the density anomaly from the reference profile with pressure [kg m-3 Pa-1]
  real :: rho00p ! The pressure-dependent (but temperature and salinity independent) reference
  real :: rhoTS  ! Density anomaly from the reference profile [kg m-3]
  real :: rhoTS0 ! A contribution to density from temperature and salinity anomalies at the
  real :: rhoTS1 ! A density contribution proportional to pressure [kg m-3 Pa-1]
  real :: rhoTS2 ! A density contribution proportional to pressure**2 [kg m-3 Pa-2]
  real :: rhoTS3 ! A density contribution proportional to pressure**3 [kg m-3 Pa-3]
  real :: rho0S0 ! Salinity dependent density at the surface pressure and zero temperature [kg m-3]
  zt = T
  zs = SQRT( ABS( S + rdeltaS ) * r1_S0 )  ! square root of normalized salinity plus an offset [nondim]
  zp = pressure

  ! The next two lines should be used if it is necessary to convert potential temperature and
  ! practical salinity to conservative temperature and absolute salinity.
  ! zt = gsw_ct_from_pt(S,T) ! Convert potential temp to conservative temp [degC]
  ! zs = SQRT( ABS( gsw_sr_from_sp(S) + rdeltaS ) * r1_S0 ) ! Convert S from practical to absolute salinity.

  rhoTS3 = EOS003 + (zs*EOS103 + zt*EOS013)
  rhoTS2 = EOS002 + (zs*(EOS102 +  zs*EOS202) &
                   + zt*(EOS012 + (zs*EOS112 + zt*EOS022)) )
  rhoTS1 = EOS001 + (zs*(EOS101 +  zs*(EOS201 +  zs*(EOS301 +  zs*EOS401))) &
                   + zt*(EOS011 + (zs*(EOS111 +  zs*(EOS211 +  zs*EOS311)) &
                                 + zt*(EOS021 + (zs*(EOS121 +  zs*EOS221) &
                                               + zt*(EOS031 + (zs*EOS131 + zt*EOS041)) )) )) )

  rhoTS0 = zt*(EOS010 &
             + (zs*(EOS110 +  zs*(EOS210 +  zs*(EOS310 +  zs*(EOS410 +  zs*EOS510)))) &
              + zt*(EOS020 + (zs*(EOS120 +  zs*(EOS220 +  zs*(EOS320 +  zs*EOS420))) &
                            + zt*(EOS030 + (zs*(EOS130 +  zs*(EOS230 +  zs*EOS330)) &
                                          + zt*(EOS040 + (zs*(EOS140 +  zs*EOS240) &
                                                        + zt*(EOS050 + (zs*EOS150 + zt*EOS060)) )) )) )) ) )

  rho0S0 = EOS000 + zs*(EOS100 + zs*(EOS200 + zs*(EOS300 + zs*(EOS400 + zs*(EOS500 + zs*EOS600)))))

  rho00p = zp*(R00 + zp*(R01 + zp*(R02 + zp*(R03 + zp*(R04 + zp*R05)))))

  rhoTS  = (rhoTS0 + rho0S0) + zp*(rhoTS1 + zp*(rhoTS2 +  zp*rhoTS3))
  rho = rhoTS + rho00p ! In situ density [kg m-3]

  drho00p_dp = R00 + zp*(2.*R01 + zp*(3.*R02 + zp*(4.*R03 + zp*(5.*R04 + zp*(6.*R05)))))
  drhoTS_dp  = rhoTS1 + zp*(2.*rhoTS2 + zp*(3.*rhoTS3))
  drho_dp = drhoTS_dp + drho00p_dp ! Compressibility [s2 m-2]

end procedure calculate_compress_elem_Roquet_rho
module procedure EoS_fit_range_Roquet_rho
  if (present(T_min)) T_min = -6.0
  if (present(T_max)) T_max = 40.0
  if (present(S_min)) S_min =  0.0
  if (present(S_max)) S_max = 42.0
  if (present(p_min)) p_min = 0.0
  if (present(p_max)) p_max = 1.0e8

end procedure EoS_fit_range_Roquet_rho
module procedure calculate_density_array_Roquet_rho
  integer :: j
  if (present(rho_ref)) then
    do j = start, start+npts-1
      rho(j) = density_anomaly_elem_Roquet_rho(this, T(j), S(j), pressure(j), rho_ref)
    enddo
  else
    do j = start, start+npts-1
      rho(j) = density_elem_Roquet_rho(this, T(j), S(j), pressure(j))
    enddo
  endif

end procedure calculate_density_array_Roquet_rho
module procedure calculate_spec_vol_array_Roquet_rho
  integer :: j
  if (present(spv_ref)) then
    do j = start, start+npts-1
      specvol(j) = spec_vol_anomaly_elem_Roquet_rho(this, T(j), S(j), pressure(j), spv_ref)
    enddo
  else
    do j = start, start+npts-1
      specvol(j) = spec_vol_elem_Roquet_rho(this, T(j), S(j), pressure(j) )
    enddo
  endif

end procedure calculate_spec_vol_array_Roquet_rho
end submodule MOM_EOS_Roquet_rho_s
