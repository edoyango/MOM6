submodule (MOM_EOS_UNESCO) MOM_EOS_UNESCO_s
  implicit none
contains
module procedure density_elem_UNESCO
  real :: t1   ! A copy of the temperature at a point [degC]
  real :: s1   ! A copy of the salinity at a point [PSU]
  real :: p1   ! Pressure converted to bars [bar]
  real :: s12  ! The square root of salinity [PSU1/2]
  real :: rho0 ! Density at 1 bar pressure [kg m-3]
  real :: sig0 ! The anomaly of rho0 from R00 [kg m-3]
  real :: ks   ! The secant bulk modulus [bar]
  p1 = pressure*1.0e-5 ; t1 = T
  s1 = max(S, 0.0) ; s12 = sqrt(s1)

  ! Compute rho(s,theta,p=0) - (same as rho(s,t_insitu,p=0) ).
  sig0 = ( t1*(R01 + t1*(R02 + t1*(R03 + t1*(R04 + t1*R05)))) + &
           s1*((R10 + t1*(R11 + t1*(R12 + t1*(R13 + t1*R14)))) + &
               (s12*(R60 + t1*(R61 + t1*R62)) + s1*R20)) )
  rho0 = R00 + sig0

  ! Compute rho(s,theta,p), first calculating the secant bulk modulus.
  ks = (S000 + ( t1*(S010 + t1*(S020 + t1*(S030 + t1*S040))) + &
                 s1*((S100 + t1*(S110 + t1*(S120 + t1*S130))) + s12*(S600 + t1*(S610 + t1*S620))) )) + &
       p1*( (S001 + ( t1*(S011 + t1*(S021 + t1*S031)) + &
                      s1*((S101 + t1*(S111 + t1*S121)) + s12*S601) )) + &
            p1*(S002 + ( t1*(S012 + t1*S022) + s1*(S102 + t1*(S112 + t1*S122)) )) )

  density_elem_UNESCO = rho0*ks / (ks - p1)

end procedure density_elem_UNESCO
module procedure density_anomaly_elem_UNESCO
  real :: t1   ! A copy of the temperature at a point [degC]
  real :: s1   ! A copy of the salinity at a point [PSU]
  real :: p1   ! Pressure converted to bars [bar]
  real :: s12  ! The square root of salinity [PSU1/2]
  real :: rho0 ! Density at 1 bar pressure [kg m-3]
  real :: sig0 ! The anomaly of rho0 from R00 [kg m-3]
  real :: ks   ! The secant bulk modulus [bar]
  p1 = pressure*1.0e-5 ; t1 = T
  s1 = max(S, 0.0) ; s12 = sqrt(s1)

  ! Compute rho(s,theta,p=0) - (same as rho(s,t_insitu,p=0) ).
  sig0 = ( t1*(R01 + t1*(R02 + t1*(R03 + t1*(R04 + t1*R05)))) + &
           s1*((R10 + t1*(R11 + t1*(R12 + t1*(R13 + t1*R14)))) + &
               (s12*(R60 + t1*(R61 + t1*R62)) + s1*R20)) )
  rho0 = R00 + sig0

  ! Compute rho(s,theta,p), first calculating the secant bulk modulus.
  ks = (S000 + ( t1*(S010 + t1*(S020 + t1*(S030 + t1*S040))) + &
                 s1*((S100 + t1*(S110 + t1*(S120 + t1*S130))) + s12*(S600 + t1*(S610 + t1*S620))) )) + &
       p1*( (S001 + ( t1*(S011 + t1*(S021 + t1*S031)) + &
                      s1*((S101 + t1*(S111 + t1*S121)) + s12*S601) )) + &
            p1*(S002 + ( t1*(S012 + t1*S022) + s1*(S102 + t1*(S112 + t1*S122)) )) )

  density_anomaly_elem_UNESCO = ((R00 - rho_ref)*ks + (sig0*ks + p1*rho_ref)) / (ks - p1)

end procedure density_anomaly_elem_UNESCO
module procedure spec_vol_elem_UNESCO
  real :: t1   ! A copy of the temperature at a point [degC]
  real :: s1   ! A copy of the salinity at a point [PSU]
  real :: p1   ! Pressure converted to bars [bar]
  real :: s12  ! The square root of salinity [PSU1/2]l553
  real :: rho0 ! Density at 1 bar pressure [kg m-3]
  real :: ks   ! The secant bulk modulus [bar]
  p1 = pressure*1.0e-5 ; t1 = T
  s1 = max(S, 0.0) ; s12 = sqrt(s1)

  ! Compute rho(s,theta,p=0), which is the same as rho(s,t_insitu,p=0).
  rho0 = R00 + ( t1*(R01 + t1*(R02 + t1*(R03 + t1*(R04 + t1*R05)))) + &
                 s1*((R10 + t1*(R11 + t1*(R12 + t1*(R13 + t1*R14)))) + &
                     (s12*(R60 + t1*(R61 + t1*R62)) + s1*R20)) )

  ! Compute rho(s,theta,p), first calculating the secant bulk modulus.
  ks = (S000 + ( t1*(S010 + t1*(S020 + t1*(S030 + t1*S040))) + &
                 s1*((S100 + t1*(S110 + t1*(S120 + t1*S130))) + s12*(S600 + t1*(S610 + t1*S620))) )) + &
       p1*( (S001 + ( t1*(S011 + t1*(S021 + t1*S031)) + &
                      s1*((S101 + t1*(S111 + t1*S121)) + s12*S601) )) + &
            p1*(S002 + ( t1*(S012 + t1*S022) + s1*(S102 + t1*(S112 + t1*S122)) )) )

  spec_vol_elem_UNESCO = (ks - p1) / (rho0*ks)

end procedure spec_vol_elem_UNESCO
module procedure spec_vol_anomaly_elem_UNESCO
  real :: t1   ! A copy of the temperature at a point [degC]
  real :: s1   ! A copy of the salinity at a point [PSU]
  real :: p1   ! Pressure converted to bars [bar]
  real :: s12  ! The square root of salinity [PSU1/2]
  real :: rho0 ! Density at 1 bar pressure [kg m-3]
  real :: ks   ! The secant bulk modulus [bar]
  p1 = pressure*1.0e-5 ; t1 = T
  s1 = max(S, 0.0) ; s12 = sqrt(s1)

  ! Compute rho(s,theta,p=0), which is the same as rho(s,t_insitu,p=0).
  rho0 = R00 + ( t1*(R01 + t1*(R02 + t1*(R03 + t1*(R04 + t1*R05)))) + &
                 s1*((R10 + t1*(R11 + t1*(R12 + t1*(R13 + t1*R14)))) + &
                     (s12*(R60 + t1*(R61 + t1*R62)) + s1*R20)) )

  ! Compute rho(s,theta,p), first calculating the secant bulk modulus.
  ks = (S000 + ( t1*(S010 + t1*(S020 + t1*(S030 + t1*S040))) + &
                 s1*((S100 + t1*(S110 + t1*(S120 + t1*S130))) + s12*(S600 + t1*(S610 + t1*S620))) )) + &
       p1*( (S001 + ( t1*(S011 + t1*(S021 + t1*S031)) + &
                      s1*((S101 + t1*(S111 + t1*S121)) + s12*S601) )) + &
            p1*(S002 + ( t1*(S012 + t1*S022) + s1*(S102 + t1*(S112 + t1*S122)) )) )

  spec_vol_anomaly_elem_UNESCO = (ks*(1.0 - (rho0*spv_ref)) - p1) / (rho0*ks)

end procedure spec_vol_anomaly_elem_UNESCO
module procedure calculate_density_derivs_elem_UNESCO
  real :: t1       ! A copy of the temperature at a point [degC]
  real :: s1       ! A copy of the salinity at a point [PSU]
  real :: p1       ! Pressure converted to bars [bar]
  real :: s12      ! The square root of salinity [PSU1/2]
  real :: rho0     ! Density at 1 bar pressure [kg m-3]
  real :: ks       ! The secant bulk modulus [bar]
  real :: drho0_dT ! Derivative of rho0 with T [kg m-3 degC-1]
  real :: drho0_dS ! Derivative of rho0 with S [kg m-3 PSU-1]
  real :: dks_dT   ! Derivative of ks with T [bar degC-1]
  real :: dks_dS   ! Derivative of ks with S [bar psu-1]
  real :: I_denom  ! 1.0 / (ks - p1) [bar-1]
  p1 = pressure*1.0e-5 ; t1 = T
  s1 = max(S, 0.0) ; s12 = sqrt(s1)

  ! Compute rho(s,theta,p=0) and its derivatives with temperature and salinity
  rho0 = R00 + ( t1*(R01 + t1*(R02 + t1*(R03 + t1*(R04 + t1*R05)))) + &
                 s1*((R10 + t1*(R11 + t1*(R12 + t1*(R13 + t1*R14)))) + &
                     (s12*(R60 + t1*(R61 + t1*R62)) + s1*R20)) )
  drho0_dT = R01 + ( t1*(2.0*R02 + t1*(3.0*R03 + t1*(4.0*R04 + t1*(5.0*R05)))) + &
                     s1*(R11 + (t1*(2.0*R12 + t1*(3.0*R13 + t1*(4.0*R14))) + &
                                s12*(R61 + t1*(2.0*R62)) )) )
  drho0_dS = R10 + ( t1*(R11 + t1*(R12 + t1*(R13 + t1*R14))) + &
                     (1.5*(s12*(R60 + t1*(R61 + t1*R62))) + s1*(2.0*R20)) )

  ! Compute the secant bulk modulus and its derivatives with temperature and salinity
  ks = ( S000 + (t1*(S010 + t1*(S020 + t1*(S030 + t1*S040))) + &
                 s1*((S100 + t1*(S110 + t1*(S120 + t1*S130))) + s12*(S600 + t1*(S610 + t1*S620)))) ) + &
       p1*( (S001 + ( t1*(S011 + t1*(S021 + t1*S031)) + &
                      s1*((S101 + t1*(S111 + t1*S121)) + s12*S601) )) + &
            p1*(S002 + ( t1*(S012 + t1*S022) + s1*(S102 + t1*(S112 + t1*S122)) )) )
  dks_dT = ( S010 + (t1*(2.0*S020 + t1*(3.0*S030 + t1*(4.0*S040))) + &
                     s1*((S110 + t1*(2.0*S120 + t1*(3.0*S130))) + s12*(S610 + t1*(2.0*S620)))) ) + &
           p1*(((S011 + t1*(2.0*S021 + t1*(3.0*S031))) + s1*(S111 + t1*(2.0*S121)) ) + &
               p1*(S012 + t1*(2.0*S022) + s1*(S112 + t1*(2.0*S122))) )
  dks_dS = ( S100 + (t1*(S110 + t1*(S120 + t1*S130)) + 1.5*(s12*(S600 + t1*(S610 + t1*S620)))) ) + &
           p1*((S101 + t1*(S111 + t1*S121) + s12*(1.5*S601)) + &
               p1*(S102 + t1*(S112 + t1*S122)) )

  I_denom = 1.0 / (ks - p1)
  drho_dT = (ks*drho0_dT - dks_dT*((rho0*p1)*I_denom)) * I_denom
  drho_dS = (ks*drho0_dS - dks_dS*((rho0*p1)*I_denom)) * I_denom

end procedure calculate_density_derivs_elem_UNESCO
module procedure calculate_density_second_derivs_elem_UNESCO
  real :: t1          ! A copy of the temperature at a point [degC]
  real :: s1          ! A copy of the salinity at a point [PSU]
  real :: p1          ! Pressure converted to bars [bar]
  real :: s12         ! The square root of salinity [PSU1/2]
  real :: I_s12       ! The inverse of the square root of salinity [PSU-1/2]
  real :: rho0        ! Density at 1 bar pressure [kg m-3]
  real :: drho0_dT    ! Derivative of rho0 with T [kg m-3 degC-1]
  real :: drho0_dS    ! Derivative of rho0 with S [kg m-3 PSU-1]
  real :: d2rho0_dS2  ! Second derivative of rho0 with salinity [kg m-3 PSU-1]
  real :: d2rho0_dSdT ! Second derivative of rho0 with temperature and salinity [kg m-3 degC-1 PSU-1]
  real :: d2rho0_dT2  ! Second derivative of rho0 with temperature [kg m-3 degC-2]
  real :: ks          ! The secant bulk modulus [bar]
  real :: ks_0        ! The secant bulk modulus at zero pressure [bar]
  real :: ks_1        ! The linear pressure dependence of the secant bulk modulus at zero pressure [nondim]
  real :: ks_2        ! The quadratic pressure dependence of the secant bulk modulus at zero pressure [bar-1]
  real :: dks_dp      ! The derivative of the secant bulk modulus with pressure [nondim]
  real :: dks_dT      ! Derivative of the secant bulk modulus with temperature [bar degC-1]
  real :: dks_dS      ! Derivative of the secant bulk modulus with salinity [bar psu-1]
  real :: d2ks_dT2    ! Second derivative of the secant bulk modulus with temperature [bar degC-2]
  real :: d2ks_dSdT   ! Second derivative of the secant bulk modulus with salinity and temperature [bar psu-1 degC-1]
  real :: d2ks_dS2    ! Second derivative of the secant bulk modulus with salinity [bar psu-2]
  real :: d2ks_dSdp   ! Second derivative of the secant bulk modulus with salinity and pressure [psu-1]
  real :: d2ks_dTdp   ! Second derivative of the secant bulk modulus with temperature and pressure [degC-1]
  real :: I_denom     ! The inverse of the denominator of the expression for density [bar-1]
  p1 = pressure*1.0e-5 ; t1 = T
  s1 = max(S, 0.0) ; s12 = sqrt(s1)
  ! The UNESCO equation of state is a fit to density, but it chooses a form that exhibits a
  ! singularity in the second derivatives with salinity for fresh water.  To avoid this, the
  ! square root of salinity can be treated with a floor such that the contribution from the
  ! S**1.5 terms to both the surface density and the secant bulk modulus are lost to roundoff.
  ! This salinity is given by (~1e-16*S000/S600)**(2/3) ~= 3e-8 PSU, or S12 ~= 1.7e-4
  I_s12 = 1.0 / (max(s12, 1.0e-4))

  ! Calculate the density at sea level pressure and its derivatives
  rho0 = R00 + ( t1*(R01 + t1*(R02 + t1*(R03 + t1*(R04 + t1*R05)))) + &
                 s1*((R10 + t1*(R11 + t1*(R12 + t1*(R13 + t1*R14)))) + &
                     (s12*(R60 + t1*(R61 + t1*R62)) + s1*R20)) )
  drho0_dT = R01 + ( t1*(2.0*R02 + t1*(3.0*R03 + t1*(4.0*R04 + t1*(5.0*R05)))) + &
                     s1*(R11 + ( t1*(2.0*R12 + t1*(3.0*R13 + t1*(4.0*R14))) + &
                                 s12*(R61 + t1*(2.0*R62)) ) ) )
  drho0_dS = R10 + ( t1*(R11 + t1*(R12 + t1*(R13 + t1*R14))) + &
                     (1.5*(s12*(R60 + t1*(R61 + t1*R62))) + s1*(2.0*R20)) )
  d2rho0_dS2 = 0.75*(R60 + t1*(R61 + t1*R62))*I_s12 + 2.0*R20
  d2rho0_dSdT = R11 + ( t1*(2.0*R12 + t1*(3.0*R13 + t1*(4.0*R14))) + s12*(1.5*R61 + t1*(3.0*R62)) )
  d2rho0_dT2 = 2.0*R02 + ( t1*(6.0*R03 + t1*(12.0*R04 + t1*(20.0*R05))) + &
                           s1*((2.0*R12 + t1*(6.0*R13 + t1*(12.0*R14))) + s12*(2.0*R62)) )

  !  Calculate the secant bulk modulus and its derivatives
  ks_0 = S000 + ( t1*( S010 + t1*(S020 + t1*(S030 + t1*S040))) + &
                  s1*((S100 + t1*(S110 + t1*(S120 + t1*S130))) + s12*(S600 + t1*(S610 + t1*S620))) )
  ks_1 = S001 + ( t1*( S011 + t1*(S021 + t1*S031)) + &
                  s1*((S101 + t1*(S111 + t1*S121)) + s12*S601) )
  ks_2 = S002 + ( t1*( S012 + t1*S022) + s1*(S102 + t1*(S112 + t1*S122)) )

  ks = ks_0 + p1*(ks_1 + p1*ks_2)
  dks_dp = ks_1 + 2.0*p1*ks_2
  dks_dT = (S010 + ( t1*(2.0*S020 + t1*(3.0*S030 + t1*(4.0*S040))) + &
                     s1*((S110 + t1*(2.0*S120 + t1*(3.0*S130))) + s12*(S610 + t1*(2.0*S620))) )) + &
           p1*((S011 + t1*(2.0*S021 + t1*(3.0*S031)) + s1*(S111 + t1*(2.0*S121))) + &
               p1*(S012 + t1*(2.0*S022) + s1*(S112 + t1*(2.0*S122))))
  dks_dS = (S100 + ( t1*(S110 + t1*(S120 + t1*S130)) + 1.5*(s12*(S600 + t1*(S610 + t1*S620))) )) + &
           p1*((S101 + t1*(S111 + t1*S121) + s12*(1.5*S601)) + &
               p1*(S102 + t1*(S112 + t1*S122)))
  d2ks_dS2 = 0.75*((S600 + t1*(S610 + t1*S620)) + p1*S601)*I_s12
  d2ks_dSdT = (S110 + ( t1*(2.0*S120 + t1*(3.0*S130)) + s12*(1.5*S610 + t1*(3.0*S620)) )) + &
              p1*((S111 + t1*(2.0*S121)) +  p1*(S112 + t1*(2.0*S122)))
  d2ks_dT2 = 2.0*(S020 + ( t1*(3.0*S030 + t1*(6.0*S040)) + s1*((S120 + t1*(3.0*S130)) + s12*S620) )) + &
             2.0*p1*((S021 + (t1*(3.0*S031) + s1*S121)) + p1*(S022 + s1*S122))

  d2ks_dSdp = (S101 + (t1*(S111 + t1*S121) + s12*(1.5*S601))) + &
              2.0*p1*(S102 + t1*(S112 + t1*S122))
  d2ks_dTdp = (S011 + (t1*(2.0*S021 + t1*(3.0*S031)) + s1*(S111 + t1*(2.0*S121)))) + &
              2.0*p1*(S012 + t1*(2.0*S022) + s1*(S112 + t1*(2.0*S122)))
  I_denom = 1.0 / (ks - p1)

  ! Expressions for density and its first derivatives are copied here for reference:
  !   rho = rho0*ks * I_denom
  !   drho_dT = I_denom*(ks*drho0_dT - p1*rho0*I_denom*dks_dT)
  !   drho_dS = I_denom*(ks*drho0_dS - p1*rho0*I_denom*dks_dS)
  !   drho_dp = 1.0e-5 * (rho0 * I_denom**2) * (ks - dks_dp*p1)

  ! Finally calculate the second derivatives
  drho_dS_dS = I_denom * ( ks*d2rho0_dS2 - (p1*I_denom) * &
                    (2.0*drho0_dS*dks_dS + rho0*(d2ks_dS2 - 2.0*dks_dS**2*I_denom)) )
  drho_dS_dT = I_denom * (ks * d2rho0_dSdT - (p1*I_denom) * &
                      ((drho0_dT*dks_dS + drho0_dS*dks_dT) + &
                       rho0*(d2ks_dSdT - 2.0*(dks_dS*dks_dT)*I_denom)) )
  drho_dT_dT = I_denom * ( ks*d2rho0_dT2 - (p1*I_denom) * &
                    (2.0*drho0_dT*dks_dT + rho0*(d2ks_dT2 - 2.0*dks_dT**2*I_denom)) )

  ! The factor of 1.0e-5 is because pressure here is in bars, not Pa.
  drho_dS_dp = (1.0e-5 * I_denom**2) * ( (ks*drho0_dS - rho0*dks_dS) - &
                    p1*( (dks_dp*drho0_dS + rho0*d2ks_dSdp) - &
                         2.0*(rho0*dks_dS) * ((dks_dp - 1.0)*I_denom) ) )
  drho_dT_dp = (1.0e-5 * I_denom**2) * ( (ks*drho0_dT - rho0*dks_dT) - &
                    p1*( (dks_dp*drho0_dT + rho0*d2ks_dTdp) - &
                         2.0*(rho0*dks_dT) * ((dks_dp - 1.0)*I_denom) ) )

end procedure calculate_density_second_derivs_elem_UNESCO
module procedure calculate_specvol_derivs_elem_UNESCO
  real :: t1       ! A copy of the temperature at a point [degC]
  real :: s1       ! A copy of the salinity at a point [PSU]
  real :: p1       ! Pressure converted to bars [bar]
  real :: s12      ! The square root of salinity [PSU1/2]
  real :: rho0     ! Density at 1 bar pressure [kg m-3]
  real :: ks       ! The secant bulk modulus [bar]
  real :: drho0_dT ! Derivative of rho0 with T [kg m-3 degC-1]
  real :: drho0_dS ! Derivative of rho0 with S [kg m-3 PSU-1]
  real :: dks_dT   ! Derivative of ks with T [bar degC-1]
  real :: dks_dS   ! Derivative of ks with S [bar psu-1]
  real :: I_denom2 ! 1.0 / (rho0*ks)**2 [m6 kg-2 bar-2]
  p1 = pressure*1.0e-5 ; t1 = T
  s1 = max(S, 0.0) ; s12 = sqrt(s1)

  ! Compute rho(s,theta,p=0) and its derivatives with temperature and salinity
  rho0 = R00 + ( t1*(R01 + t1*(R02 + t1*(R03 + t1*(R04 + t1*R05)))) + &
                 s1*((R10 + t1*(R11 + t1*(R12 + t1*(R13 + t1*R14)))) + &
                     (s12*(R60 + t1*(R61 + t1*R62)) + s1*R20)) )
  drho0_dT = R01 + ( t1*(2.0*R02 + t1*(3.0*R03 + t1*(4.0*R04 + t1*(5.0*R05)))) + &
                     s1*(R11 + (t1*(2.0*R12 + t1*(3.0*R13 + t1*(4.0*R14))) + &
                                s12*(R61 + t1*(2.0*R62)) )) )
  drho0_dS = R10 + ( t1*(R11 + t1*(R12 + t1*(R13 + t1*R14))) + &
                     (1.5*(s12*(R60 + t1*(R61 + t1*R62))) + s1*(2.0*R20)) )

  ! Compute the secant bulk modulus and its derivatives with temperature and salinity
  ks = ( S000 + (t1*(S010 + t1*(S020 + t1*(S030 + t1*S040))) + &
                 s1*((S100 + t1*(S110 + t1*(S120 + t1*S130))) + s12*(S600 + t1*(S610 + t1*S620)))) ) + &
       p1*( (S001 + ( t1*(S011 + t1*(S021 + t1*S031)) + &
                      s1*((S101 + t1*(S111 + t1*S121)) + s12*S601) )) + &
            p1*(S002 + ( t1*(S012 + t1*S022) + s1*(S102 + t1*(S112 + t1*S122)) )) )
  dks_dT = ( S010 + (t1*(2.0*S020 + t1*(3.0*S030 + t1*(4.0*S040))) + &
                     s1*((S110 + t1*(2.0*S120 + t1*(3.0*S130))) + s12*(S610 + t1*(2.0*S620)))) ) + &
           p1*(((S011 + t1*(2.0*S021 + t1*(3.0*S031))) + s1*(S111 + t1*(2.0*S121)) ) + &
               p1*(S012 + t1*(2.0*S022) + s1*(S112 + t1*(2.0*S122))) )
  dks_dS = ( S100 + (t1*(S110 + t1*(S120 + t1*S130)) + 1.5*(s12*(S600 + t1*(S610 + t1*S620)))) ) + &
           p1*((S101 + t1*(S111 + t1*S121) + s12*(1.5*S601)) + &
               p1*(S102 + t1*(S112 + t1*S122)) )

  ! specvol = (ks - p1) / (rho0*ks) = 1/rho0 - p1/(rho0*ks)
  I_denom2 = 1.0 / (rho0*ks)**2
  dSV_dT = ((p1*rho0)*dks_dT + ((p1 - ks)*ks)*drho0_dT) * I_denom2
  dSV_dS = ((p1*rho0)*dks_dS + ((p1 - ks)*ks)*drho0_dS) * I_denom2

end procedure calculate_specvol_derivs_elem_UNESCO
module procedure calculate_compress_elem_UNESCO
  real :: t1      ! A copy of the temperature at a point [degC]
  real :: s1      ! A copy of the salinity at a point [PSU]
  real :: p1      ! Pressure converted to bars [bar]
  real :: s12     ! The square root of salinity [PSU1/2]
  real :: rho0    ! Density at 1 bar pressure [kg m-3]
  real :: ks      ! The secant bulk modulus [bar]
  real :: ks_0    ! The secant bulk modulus at zero pressure [bar]
  real :: ks_1    ! The linear pressure dependence of the secant bulk modulus at zero pressure [nondim]
  real :: ks_2    ! The quadratic pressure dependence of the secant bulk modulus at zero pressure [bar-1]
  real :: dks_dp  ! The derivative of the secant bulk modulus with pressure [nondim]
  real :: I_denom  ! 1.0 / (ks - p1) [bar-1]
  p1 = pressure*1.0e-5 ; t1 = T
  s1 = max(S, 0.0) ; s12 = sqrt(s1)

  ! Compute rho(s,theta,p=0), which is the same as rho(s,t_insitu,p=0).

  rho0 = R00 + ( t1*(R01 + t1*(R02 + t1*(R03 + t1*(R04 + t1*R05)))) + &
                 s1*((R10 + t1*(R11 + t1*(R12 + t1*(R13 + t1*R14)))) + &
                     (s12*(R60 + t1*(R61 + t1*R62)) + s1*R20)) )

  ! Calculate the secant bulk modulus and its derivative with pressure.
  ks_0 = S000 + ( t1*( S010 + t1*(S020 + t1*(S030 + t1*S040))) + &
                  s1*((S100 + t1*(S110 + t1*(S120 + t1*S130))) + s12*(S600 + t1*(S610 + t1*S620))) )
  ks_1 = S001 + ( t1*( S011 + t1*(S021 + t1*S031)) + &
                  s1*((S101 + t1*(S111 + t1*S121)) + s12*S601) )
  ks_2 = S002 + ( t1*( S012 + t1*S022) + s1*(S102 + t1*(S112 + t1*S122)) )

  ks = ks_0 + p1*(ks_1 + p1*ks_2)
  dks_dp = ks_1 + 2.0*p1*ks_2
  I_denom = 1.0 / (ks - p1)

  ! Compute the in situ density, rho(s,theta,p), and its derivative with pressure.
  rho = rho0*ks * I_denom
  ! The factor of 1.0e-5 is because pressure here is in bars, not Pa.
  drho_dp = 1.0e-5 * ((rho0 * (ks - p1*dks_dp)) * I_denom**2)

end procedure calculate_compress_elem_UNESCO
module procedure EoS_fit_range_UNESCO
  if (present(T_min)) T_min = -2.5
  if (present(T_max)) T_max = 40.0
  if (present(S_min)) S_min =  0.0
  if (present(S_max)) S_max = 42.0
  if (present(p_min)) p_min = 0.0
  if (present(p_max)) p_max = 1.0e8

end procedure EoS_fit_range_UNESCO
end submodule MOM_EOS_UNESCO_s
