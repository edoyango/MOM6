submodule (MOM_EOS_Wright_red) MOM_EOS_Wright_red_s
  implicit none
contains
module procedure density_elem_Wright_red
  real :: al0     ! The specific volume at 0 lambda in the Wright EOS [m3 kg-1]
  real :: p0      ! The pressure offset in the Wright EOS [Pa]
  real :: lambda  ! The sound speed squared at 0 alpha in the Wright EOS [m2 s-2]
  al0 = a0 + (a1*T + a2*S)
  p0 = b0 + ( b4*S + T * (b1 + (T*(b2 + b3*T) + b5*S)) )
  lambda = c0 + ( c4*S + T * (c1 + (T*(c2 + c3*T) + c5*S)) )
  density_elem_Wright_red = (pressure + p0) / (lambda + al0*(pressure + p0))

end procedure density_elem_Wright_red
module procedure density_anomaly_elem_Wright_red
  real :: al0     ! The specific volume at 0 lambda in the Wright EOS [m3 kg-1]
  real :: al_TS   ! The contributions of temperature and salinity to al0 [m3 kg-1]
  real :: p_TSp   ! A combination of the pressure and the temperature and salinity contributions to p0 [Pa]
  real :: lam_TS  ! The contributions of temperature and salinity to lambda [m2 s-2]
  real :: pa_000  ! A corrected offset to the pressure, including contributions from rho_ref [Pa]
  pa_000 = b0*(1.0 - a0*rho_ref) - rho_ref*c0
  al_TS = a1*T + a2*S
  al0 = a0 + al_TS
  p_TSp = pressure + (b4*S + T * (b1 + (T*(b2 + b3*T) + b5*S)))
  lam_TS = c4*S + T * (c1 + (T*(c2 + c3*T) + c5*S))

  ! The following two expressions are mathematically equivalent.
  ! rho = (b0 + p0_TSp) / ((c0 + lam_TS) + al0*(b0 + p0_TSp)) - rho_ref
  density_anomaly_elem_Wright_red = &
           (pa_000 + (p_TSp - rho_ref*(p_TSp*al0 + (b0*al_TS + lam_TS)))) / &
           ( (c0 + lam_TS) + al0*(b0 + p_TSp) )

end procedure density_anomaly_elem_Wright_red
module procedure spec_vol_elem_Wright_red
  real :: al0     ! The specific volume at 0 lambda in the Wright EOS [m3 kg-1]
  real :: p0      ! The pressure offset in the Wright EOS [Pa]
  real :: lambda  ! The sound speed squared at 0 alpha in the Wright EOS [m2 s-2], perhaps with
  real :: al_TS   ! The contributions of temperature and salinity to al0 [m3 kg-1]
  real :: p_TSp   ! A combination of the pressure and the temperature and salinity contributions to p0 [Pa]
  real :: lam_000 ! A corrected offset to lambda, including contributions from spv_ref [m2 s-2]
  al0 = a0 + (a1*T + a2*S)
  p0 = b0 + ( b4*S + T * (b1 + (T*(b2 + b3*T) + b5*S)) )
  lambda = c0 + ( c4*S + T * (c1 + (T*(c2 + c3*T) + c5*S)) )
  spec_vol_elem_Wright_red = al0 + lambda / (pressure + p0)

end procedure spec_vol_elem_Wright_red
module procedure spec_vol_anomaly_elem_Wright_red
  real :: al0     ! The specific volume at 0 lambda in the Wright EOS [m3 kg-1]
  real :: p0      ! The pressure offset in the Wright EOS [Pa]
  real :: lambda  ! The sound speed squared at 0 alpha in the Wright EOS [m2 s-2], perhaps with
  real :: al_TS   ! The contributions of temperature and salinity to al0 [m3 kg-1]
  real :: p_TSp   ! A combination of the pressure and the temperature and salinity contributions to p0 [Pa]
  real :: lam_000 ! A corrected offset to lambda, including contributions from spv_ref [m2 s-2]
    lam_000 = c0 + (a0 - spv_ref)*b0
    al_TS = a1*T + a2*S
    p_TSp = pressure + (b4*S + T * (b1 + (T*(b2 + b3*T) + b5*S)))
    lambda = lam_000 + ( c4*S + T * (c1 + (T*(c2 + c3*T) + c5*S)) )
    ! This is equivalent to the expression below minus spv_ref, but less sensitive to roundoff.
    spec_vol_anomaly_elem_Wright_red = al_TS + (lambda + (a0 - spv_ref)*p_TSp) / (b0 + p_TSp)

end procedure spec_vol_anomaly_elem_Wright_red
module procedure calculate_density_derivs_elem_Wright_red
  real :: al0     ! The specific volume at 0 lambda in the Wright EOS [m3 kg-1]
  real :: p0      ! The pressure offset in the Wright EOS [Pa]
  real :: lambda  ! The sound speed squared at 0 alpha in the Wright EOS [m2 s-2]
  real :: I_denom2 ! The inverse of the square of the denominator of density in the Wright EOS [s4 m-4]
  al0 = a0 + (a1*T + a2*S)
  p0 = b0 + ( b4*S + T * (b1 + (T*(b2 + b3*T) + b5*S)) )
  lambda = c0 + ( c4*S + T * (c1 + (T*(c2 + c3*T) + c5*S)) )

  I_denom2 = 1.0 / (lambda + al0*(pressure + p0))**2
  drho_dT = I_denom2 * (lambda * (b1 + (T*(2.0*b2 + 3.0*b3*T) + b5*S)) - &
     (pressure+p0) * ( (pressure+p0)*a1 + (c1 + (T*(c2*2.0 + c3*3.0*T) + c5*S)) ))
  drho_dS = I_denom2 * (lambda * (b4 + b5*T) - &
     (pressure+p0) * ( (pressure+p0)*a2 + (c4 + c5*T) ))

end procedure calculate_density_derivs_elem_Wright_red
module procedure calculate_density_second_derivs_elem_Wright_red
  real :: al0     ! The specific volume at 0 lambda in the Wright EOS [m3 kg-1]
  real :: lambda  ! The sound speed squared at 0 alpha in the Wright EOS [m2 s-2]
  real :: p_p0    ! A local work variable combining the pressure and pressure
  real :: dp0_dT  ! The partial derivative of p0 with temperature [Pa degC-1]
  real :: dp0_dS  ! The partial derivative of p0 with salinity [Pa PSU-1]
  real :: dlam_dT ! The partial derivative of lambda with temperature [m2 s-2 degC-1]
  real :: dlam_dS ! The partial derivative of lambda with salinity [m2 s-2 degC-1]
  real :: dRdT_num  ! The numerator in the expression for drho_dT [Pa m2 s-2 degC-1] = [kg m s-4 degC-1]
  real :: dRdS_num  ! The numerator in the expression for drho_ds [Pa m2 s-2 PSU-1] = [kg m s-4 PSU-1]
  real :: ddenom_dT ! The derivative of the denominator of density in the Wright EOS with temperature [m2 s-2 deg-1]
  real :: ddenom_dS ! The derivative of the denominator of density in the Wright EOS with salinity [m2 s-2 PSU-1]
  real :: I_denom   ! The inverse of the denominator of density in the Wright EOS [s2 m-2]
  real :: I_denom2  ! The inverse of the square of the denominator of density in the Wright EOS [s4 m-4]
  real :: I_denom3  ! The inverse of the cube of the denominator of density in the Wright EOS [s6 m-6]
  al0 = a0 + (a1*T + a2*S)
  p_p0 = pressure + ( b0 + (b4*S + T*(b1 + (b5*S + T*(b2 + b3*T)))) ) ! P + p0
  lambda = c0 + ( c4*S + T * (c1 + (T*(c2 + c3*T) + c5*S)) )
  dp0_dT = b1 + (b5*S + T*(2.*b2 + 3.*b3*T))
  dp0_dS = b4 + b5*T
  dlam_dT = c1 + (c5*S + T*(2.*c2 + 3.*c3*T))
  dlam_dS = c4 + c5*T
  I_denom = 1.0 / (lambda + al0*p_p0)
  I_denom2 = I_denom*I_denom
  I_denom3 = I_denom*I_denom2

  ddenom_dS = (dlam_dS + a2*p_p0) + al0*dp0_dS
  ddenom_dT = (dlam_dT + a1*p_p0) + al0*dp0_dT
  dRdS_num = dp0_dS*lambda - p_p0*(dlam_dS + a2*p_p0)
  dRdT_num = dp0_dT*lambda - p_p0*(dlam_dT + a1*p_p0)

  ! In deriving the following, it is useful to note that:
  !   rho = p_p0 / (lambda + al0*p_p0)
  !   drho_dp = lambda * I_denom2
  !   drho_dT = (dp0_dT*lambda - p_p0*(dlam_dT + a1*p_p0)) * I_denom2 = dRdT_num * I_denom2
  !   drho_dS = (dp0_dS*lambda - p_p0*(dlam_dS + a2*p_p0)) * I_denom2 = dRdS_num * I_denom2
  drho_ds_ds = -2.*(p_p0*(a2*dp0_dS)) * I_denom2 - 2.*(dRdS_num*ddenom_dS) * I_denom3
  drho_ds_dt = ((b5*lambda - p_p0*(c5 + 2.*a2*dp0_dT)) + (dp0_dS*dlam_dT - dp0_dT*dlam_dS))*I_denom2 - &
                  2.*(ddenom_dT*dRdS_num) * I_denom3
  drho_dt_dt = 2.*((b2 + 3.*b3*T)*lambda - p_p0*((c2 + 3.*c3*T) + a1*dp0_dT))*I_denom2 - &
                  2.*(dRdT_num * ddenom_dT) * I_denom3

  ! The following is a rearranged form that is equivalent to
  ! drho_ds_dp = dlam_dS * I_denom2 - 2.0 * lambda * (dlam_dS + a2*p_p0 + al0*dp0_ds) * Idenom3
  drho_ds_dp = (-dlam_dS - 2.*a2*p_p0) * I_denom2 - (2.*al0*dRdS_num) * I_denom3
  drho_dt_dp = (-dlam_dT - 2.*a1*p_p0) * I_denom2 - (2.*al0*dRdT_num) * I_denom3

end procedure calculate_density_second_derivs_elem_Wright_red
module procedure calculate_specvol_derivs_elem_Wright_red
  real :: p0      ! The pressure offset in the Wright EOS [Pa]
  real :: lambda  ! The sound speed squared at 0 alpha in the Wright EOS [m2 s-2]
  real :: I_denom ! The inverse of the denominator of specific volume in the Wright EOS [Pa-1]
  p0 = b0 + ( b4*S + T * (b1 + (T*(b2 + b3*T) + b5*S)) )
  lambda = c0 + ( c4*S + T * (c1 + (T*(c2 + c3*T) + c5*S)) )

  ! SV = al0 + lambda / (pressure + p0)

  I_denom = 1.0 / (pressure + p0)
  dSV_dT = a1 + I_denom * ((c1 + (T*(2.0*c2 + 3.0*c3*T) + c5*S)) - &
                              (I_denom * lambda) * (b1 + (T*(2.0*b2 + 3.0*b3*T) + b5*S)))
  dSV_dS = a2 + I_denom * ((c4 + c5*T) - &
                              (I_denom * lambda) * (b4 + b5*T))

end procedure calculate_specvol_derivs_elem_Wright_red
module procedure calculate_compress_elem_Wright_red
  real :: al0     ! The specific volume at 0 lambda in the Wright EOS [m3 kg-1]
  real :: p0      ! The pressure offset in the Wright EOS [Pa]
  real :: lambda  ! The sound speed squared at 0 alpha in the Wright EOS [m2 s-2]
  real :: I_denom ! The inverse of the denominator of density in the Wright EOS [s2 m-2]
  al0 = a0 + (a1*T + a2*S)
  p0 = b0 + ( b4*S + T * (b1 + (T*(b2 + b3*T) + b5*S)) )
  lambda = c0 + ( c4*S + T * (c1 + (T*(c2 + c3*T) + c5*S)) )

  I_denom = 1.0 / (lambda + al0*(pressure + p0))
  rho = (pressure + p0) * I_denom
  drho_dp = lambda * I_denom**2

end procedure calculate_compress_elem_Wright_red
module procedure avg_spec_vol_Wright_red
  real :: al0        ! A term in the Wright EOS [m3 kg-1]
  real :: p0         ! A term in the Wright EOS [Pa]
  real :: lambda     ! A term in the Wright EOS [m2 s-2]
  real :: eps2       ! The square of a nondimensional ratio [nondim]
  real :: I_pterm    ! The inverse of p0 plus p_ave [Pa-1].
  real, parameter :: C1_3 = 1.0/3.0, C1_7 = 1.0/7.0, C1_9 = 1.0/9.0 ! Rational constants [nondim]
  integer :: j
  do j=start,start+npts-1
    al0 = a0 + (a1*T(j) + a2*S(j))
    p0 = b0 + ( b4*S(j) + T(j) * (b1 + (T(j)*(b2 + b3*T(j)) + b5*S(j))) )
    lambda = c0 + ( c4*S(j) + T(j) * (c1 + (T(j)*(c2 + c3*T(j)) + c5*S(j))) )

    I_pterm = 1.0 / (p0 + (p_t(j) + 0.5*dp(j)))
    eps2 = (0.5 * dp(j) * I_pterm)**2
    SpV_avg(j) = al0 + (lambda * I_pterm) * &
                         (1.0 + eps2*(C1_3 + eps2*(0.2 + eps2*(C1_7 + eps2*C1_9))))
  enddo
end procedure avg_spec_vol_Wright_red
module procedure EoS_fit_range_Wright_red
  if (present(T_min)) T_min = -2.0
  if (present(T_max)) T_max = 30.0
  if (present(S_min)) S_min = 28.0
  if (present(S_max)) S_max = 38.0
  if (present(p_min)) p_min = 0.0
  if (present(p_max)) p_max = 5.0e7

end procedure EoS_fit_range_Wright_red
module procedure int_density_dz_wright_red
  real, dimension(HI%isd:HI%ied,HI%jsd:HI%jed) :: al0_2d ! A term in the Wright EOS [m3 kg-1]
  real, dimension(HI%isd:HI%ied,HI%jsd:HI%jed) :: p0_2d  ! A term in the Wright EOS [Pa]
  real, dimension(HI%isd:HI%ied,HI%jsd:HI%jed) :: lambda_2d ! A term in the Wright EOS [m2 s-2]
  real, dimension(HI%isd:HI%ied,HI%jsd:HI%jed) :: z0pres ! The height at which the pressure is zero [Z ~> m]
  real :: al0        ! A term in the Wright EOS [m3 kg-1]
  real :: p0         ! A term in the Wright EOS [Pa]
  real :: lambda     ! A term in the Wright EOS [m2 s-2]
  real :: rho_anom   ! The density anomaly from rho_ref [kg m-3].
  real :: eps, eps2  ! A nondimensional ratio and its square [nondim]
  real :: rem        ! [kg m-1 s-2]
  real :: GxRho      ! The gravitational acceleration times density and unit conversion factors [Pa Z-1 ~> kg m-2 s-2]
  real :: g_Earth    ! The gravitational acceleration [m2 Z-1 s-2 ~> m s-2]
  real :: I_Rho      ! The inverse of the Boussinesq density [m3 kg-1]
  real :: rho_ref_mks ! The reference density in MKS units [kg m-3]
  real :: p_ave      ! The layer averaged pressure [Pa]
  real :: I_al0      ! The inverse of al0 [kg m-3]
  real :: I_Lzz      ! The inverse of the denominator [Pa-1]
  real :: dz         ! The layer thickness [Z ~> m].
  real :: hWght      ! A pressure-thickness below topography [Z ~> m].
  real :: hL, hR     ! Pressure-thicknesses of the columns to the left and right [Z ~> m].
  real :: iDenom     ! The inverse of the denominator in the weights [Z-2 ~> m-2].
  real :: hWt_LL, hWt_LR ! hWt_LA is the weighted influence of A on the left column [nondim].
  real :: hWt_RL, hWt_RR ! hWt_RA is the weighted influence of A on the right column [nondim].
  real :: wt_L, wt_R ! The linear weights of the left and right columns [nondim].
  real :: wtT_L, wtT_R ! The weights for tracers from the left and right columns [nondim].
  real :: intz(5)    ! The gravitational acceleration times the integrals of density
  real :: Pa_to_RL2_T2 ! A conversion factor of pressures from Pa to the output units indicated by
  real :: a1s        ! Partly rescaled version of a1 [m3 kg-1 C-1 ~> m3 kg-1 degC-1]
  real :: a2s        ! Partly rescaled version of a2 [m3 kg-1 S-1 ~> m3 kg-1 PSU-1]
  real :: b1s        ! Partly rescaled version of b1 [Pa C-1 ~> Pa degC-1]
  real :: b2s        ! Partly rescaled version of b2 [Pa C-2 ~> Pa degC-2]
  real :: b3s        ! Partly rescaled version of b3 [Pa C-3 ~> Pa degC-3]
  real :: b4s        ! Partly rescaled version of b4 [Pa S-1 ~> Pa PSU-1]
  real :: b5s        ! Partly rescaled version of b5 [Pa C-1 S-1 ~> Pa degC-1 PSU-1]
  real :: c1s        ! Partly rescaled version of c1 [m2 s-2 C-1 ~> m2 s-2 degC-1]
  real :: c2s        ! Partly rescaled version of c2 [m2 s-2 C-2 ~> m2 s-2 degC-2]
  real :: c3s        ! Partly rescaled version of c3 [m2 s-2 C-3 ~> m2 s-2 degC-3]
  real :: c4s        ! Partly rescaled version of c4 [m2 s-2 S-1 ~> m2 s-2 PSU-1]
  real :: c5s        ! Partly rescaled version of c5 [m2 s-2 C-1 S-1 ~> m2 s-2 degC-1 PSU-1]
  logical :: do_massWeight ! Indicates whether to do mass weighting.
  logical :: top_massWeight ! Indicates whether to do mass weighting the sea surface
  real, parameter :: C1_3 = 1.0/3.0, C1_7 = 1.0/7.0    ! Rational constants [nondim]
  real, parameter :: C1_9 = 1.0/9.0, C1_90 = 1.0/90.0  ! Rational constants [nondim]
  integer :: is, ie, js, je, Isq, Ieq, Jsq, Jeq, i, j, m
  Isq = HI%IscB ; Ieq = HI%IecB
  Jsq = HI%JscB ; Jeq = HI%JecB
  is = HI%isc ; ie = HI%iec
  js = HI%jsc ; je = HI%jec

  if (present(pres_scale)) then
    GxRho = pres_scale * G_e * rho_0 ; g_Earth = pres_scale * G_e
    Pa_to_RL2_T2 = 1.0 / pres_scale
  else
    GxRho = G_e * rho_0 ; g_Earth = G_e
    Pa_to_RL2_T2 = 1.0
  endif
  if (present(rho_scale)) then
    g_Earth = g_Earth * rho_scale
    rho_ref_mks = rho_ref / rho_scale ; I_Rho = rho_scale / rho_0
  else
    rho_ref_mks = rho_ref ; I_Rho = 1.0 / rho_0
  endif
  if (present(Z_0p)) then
    do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
      z0pres(i,j) = Z_0p(i,j)
    enddo ; enddo
  else
    z0pres(:,:) = 0.0
  endif

  a1s = a1 ; a2s = a2
  b1s = b1 ; b2s = b2 ; b3s = b3 ; b4s = b4 ; b5s = b5
  c1s = c1 ; c2s = c2 ; c3s = c3 ; c4s = c4 ; c5s = c5

  if (present(temp_scale)) then ; if (temp_scale /= 1.0) then
    a1s = a1s * temp_scale
    b1s = b1s * temp_scale    ; b2s = b2s * temp_scale**2
    b3s = b3s * temp_scale**3 ; b5s = b5s * temp_scale
    c1s = c1s * temp_scale    ; c2s = c2s * temp_scale**2
    c3s = c3s * temp_scale**3 ; c5s = c5s * temp_scale
  endif ; endif

  if (present(saln_scale)) then ; if (saln_scale /= 1.0) then
    a2s = a2s * saln_scale
    b4s = b4s * saln_scale ; b5s = b5s * saln_scale
    c4s = c4s * saln_scale ; c5s = c5s * saln_scale
  endif ; endif

  do_massWeight = .false. ; top_massWeight = .false.
  if (present(MassWghtInterp)) then
    do_massWeight = BTEST(MassWghtInterp, 0) ! True for odd values
    top_massWeight = BTEST(MassWghtInterp, 1) ! True if the 2 bit is set
  endif

  do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
    al0_2d(i,j) = a0 + (a1s*T(i,j) + a2s*S(i,j))
    p0_2d(i,j) = b0 + ( b4s*S(i,j) + T(i,j) * (b1s + (T(i,j)*(b2s + b3s*T(i,j)) + b5s*S(i,j))) )
    lambda_2d(i,j) = c0 + ( c4s*S(i,j) + T(i,j) * (c1s + (T(i,j)*(c2s + c3s*T(i,j)) + c5s*S(i,j))) )

    al0 = al0_2d(i,j) ; p0 = p0_2d(i,j) ; lambda = lambda_2d(i,j)

    dz = z_t(i,j) - z_b(i,j)
    p_ave = -GxRho*(0.5*(z_t(i,j)+z_b(i,j)) - z0pres(i,j))

    I_al0 = 1.0 / al0
    I_Lzz = 1.0 / ((p0 + p_ave) + lambda * I_al0)
    eps = 0.5*(GxRho*dz)*I_Lzz ; eps2 = eps*eps

!     rho(j) = (pressure(j) + p0) / (lambda + al0*(pressure(j) + p0))

    rho_anom = (p0 + p_ave)*(I_Lzz*I_al0) - rho_ref_mks
    rem = (I_Rho * (lambda * I_al0**2)) * (eps2 * (C1_3 + eps2*(0.2 + eps2*(C1_7 + C1_9*eps2))))
    dpa(i,j) = Pa_to_RL2_T2 * ((g_Earth*rho_anom)*dz - 2.0*eps*rem)
    if (present(intz_dpa)) &
      intz_dpa(i,j) = Pa_to_RL2_T2 * (0.5*(g_Earth*rho_anom)*dz**2 - dz*((1.0+eps)*rem))
  enddo ; enddo

  if (present(intx_dpa)) then ; do j=js,je ; do I=Isq,Ieq
    ! hWght is the distance measure by which the cell is violation of
    ! hydrostatic consistency. For large hWght we bias the interpolation of
    ! T & S along the top and bottom integrals, akin to thickness weighting.
    hWght = 0.0
    if (do_massWeight) &
      hWght = max(0., -bathyT(i,j)-z_t(i+1,j), -bathyT(i+1,j)-z_t(i,j))
    if (top_massWeight) &
      hWght = max(hWght, z_b(i+1,j)-SSH(i,j), z_b(i,j)-SSH(i+1,j))
    if (hWght > 0.) then
      hL = (z_t(i,j) - z_b(i,j)) + dz_neglect
      hR = (z_t(i+1,j) - z_b(i+1,j)) + dz_neglect
      hWght = hWght * ( (hL-hR)/(hL+hR) )**2
      iDenom = 1.0 / ( hWght*(hR + hL) + hL*hR )
      hWt_LL = (hWght*hL + hR*hL) * iDenom ; hWt_LR = (hWght*hR) * iDenom
      hWt_RR = (hWght*hR + hR*hL) * iDenom ; hWt_RL = (hWght*hL) * iDenom
    else
      hWt_LL = 1.0 ; hWt_LR = 0.0 ; hWt_RR = 1.0 ; hWt_RL = 0.0
    endif

    intz(1) = dpa(i,j) ; intz(5) = dpa(i+1,j)
    do m=2,4
      wt_L = 0.25*real(5-m) ; wt_R = 1.0-wt_L
      wtT_L = (wt_L*hWt_LL) + (wt_R*hWt_RL) ; wtT_R = (wt_L*hWt_LR) + (wt_R*hWt_RR)

      al0 = (wtT_L*al0_2d(i,j)) + (wtT_R*al0_2d(i+1,j))
      p0 = (wtT_L*p0_2d(i,j)) + (wtT_R*p0_2d(i+1,j))
      lambda = (wtT_L*lambda_2d(i,j)) + (wtT_R*lambda_2d(i+1,j))

      dz = (wt_L*(z_t(i,j) - z_b(i,j))) + (wt_R*(z_t(i+1,j) - z_b(i+1,j)))
      p_ave = -GxRho*((wt_L * (0.5*(z_t(i,j)+z_b(i,j)) - z0pres(i,j))) + &
                      (wt_R * (0.5*(z_t(i+1,j)+z_b(i+1,j)) - z0pres(i+1,j))))

      I_al0 = 1.0 / al0
      I_Lzz = 1.0 / ((p0 + p_ave) + lambda * I_al0)
      eps = 0.5*(GxRho*dz)*I_Lzz ; eps2 = eps*eps

      intz(m) = Pa_to_RL2_T2 * ( (g_Earth*dz) * ((p0 + p_ave)*(I_Lzz*I_al0) - rho_ref_mks) - 2.0*eps * &
                  (I_Rho * (lambda * I_al0**2)) * (eps2 * (C1_3 + eps2*(0.2 + eps2*(C1_7 + C1_9*eps2)))) )
    enddo
    ! Use Boole's rule to integrate the values.
    intx_dpa(i,j) = C1_90*(7.0*(intz(1)+intz(5)) + 32.0*(intz(2)+intz(4)) + 12.0*intz(3))
  enddo ; enddo ; endif

  if (present(inty_dpa)) then ; do J=Jsq,Jeq ; do i=is,ie
    ! hWght is the distance measure by which the cell is violation of
    ! hydrostatic consistency. For large hWght we bias the interpolation of
    ! T & S along the top and bottom integrals, akin to thickness weighting.
    hWght = 0.0
    if (do_massWeight) &
      hWght = max(0., -bathyT(i,j)-z_t(i,j+1), -bathyT(i,j+1)-z_t(i,j))
    if (top_massWeight) &
      hWght = max(hWght, z_b(i,j+1)-SSH(i,j), z_b(i,j)-SSH(i,j+1))
    if (hWght > 0.) then
      hL = (z_t(i,j) - z_b(i,j)) + dz_neglect
      hR = (z_t(i,j+1) - z_b(i,j+1)) + dz_neglect
      hWght = hWght * ( (hL-hR)/(hL+hR) )**2
      iDenom = 1.0 / ( hWght*(hR + hL) + hL*hR )
      hWt_LL = (hWght*hL + hR*hL) * iDenom ; hWt_LR = (hWght*hR) * iDenom
      hWt_RR = (hWght*hR + hR*hL) * iDenom ; hWt_RL = (hWght*hL) * iDenom
    else
      hWt_LL = 1.0 ; hWt_LR = 0.0 ; hWt_RR = 1.0 ; hWt_RL = 0.0
    endif

    intz(1) = dpa(i,j) ; intz(5) = dpa(i,j+1)
    do m=2,4
      wt_L = 0.25*real(5-m) ; wt_R = 1.0-wt_L
      wtT_L = (wt_L*hWt_LL) + (wt_R*hWt_RL) ; wtT_R = (wt_L*hWt_LR) + (wt_R*hWt_RR)

      al0 = (wtT_L*al0_2d(i,j)) + (wtT_R*al0_2d(i,j+1))
      p0 = (wtT_L*p0_2d(i,j)) + (wtT_R*p0_2d(i,j+1))
      lambda = (wtT_L*lambda_2d(i,j)) + (wtT_R*lambda_2d(i,j+1))

      dz = (wt_L*(z_t(i,j) - z_b(i,j))) + (wt_R*(z_t(i,j+1) - z_b(i,j+1)))
      p_ave = -GxRho*((wt_L*(0.5*(z_t(i,j)+z_b(i,j))-z0pres(i,j))) + &
                      (wt_R*(0.5*(z_t(i,j+1)+z_b(i,j+1))-z0pres(i,j+1))))

      I_al0 = 1.0 / al0
      I_Lzz = 1.0 / ((p0 + p_ave) + lambda * I_al0)
      eps = 0.5*(GxRho*dz)*I_Lzz ; eps2 = eps*eps

      intz(m) = Pa_to_RL2_T2 * ( (g_Earth*dz) * ((p0 + p_ave)*(I_Lzz*I_al0) - rho_ref_mks) - 2.0*eps * &
                  (I_Rho * (lambda * I_al0**2)) * (eps2 * (C1_3 + eps2*(0.2 + eps2*(C1_7 + C1_9*eps2)))) )
    enddo
    ! Use Boole's rule to integrate the values.
    inty_dpa(i,j) = C1_90*(7.0*(intz(1)+intz(5)) + 32.0*(intz(2)+intz(4)) + 12.0*intz(3))
  enddo ; enddo ; endif

end procedure int_density_dz_wright_red
module procedure int_spec_vol_dp_wright_red
  real, dimension(HI%isd:HI%ied,HI%jsd:HI%jed) :: al0_2d ! A term in the Wright EOS [R-1 ~> m3 kg-1]
  real, dimension(HI%isd:HI%ied,HI%jsd:HI%jed) :: p0_2d  ! A term in the Wright EOS [R L2 T-2 ~> Pa]
  real, dimension(HI%isd:HI%ied,HI%jsd:HI%jed) :: lambda_2d ! A term in the Wright EOS [L2 T-2 ~> m2 s-2]
  real :: al0        ! A term in the Wright EOS [R-1 ~> m3 kg-1]
  real :: p0         ! A term in the Wright EOS [R L2 T-2 ~> Pa]
  real :: lambda     ! A term in the Wright EOS [L2 T-2 ~> m2 s-2]
  real :: al0_scale  ! Scaling factor to convert al0 from MKS units [kg m-3 R-1 ~> 1]
  real :: p0_scale   ! Scaling factor to convert p0 from MKS units [R L2 T-2 Pa-1 ~> 1]
  real :: lam_scale  ! Scaling factor to convert lambda from MKS units [L2 s2 T-2 m-2 ~> 1]
  real :: p_ave      ! The layer average pressure [R L2 T-2 ~> Pa]
  real :: rem        ! [L2 T-2 ~> m2 s-2]
  real :: eps, eps2  ! A nondimensional ratio and its square [nondim]
  real :: alpha_anom ! The depth averaged specific volume anomaly [R-1 ~> m3 kg-1].
  real :: dp         ! The pressure change through a layer [R L2 T-2 ~> Pa].
  real :: hWght      ! A pressure-thickness below topography [R L2 T-2 ~> Pa].
  real :: hL, hR     ! Pressure-thicknesses of the columns to the left and right [R L2 T-2 ~> Pa].
  real :: iDenom     ! The inverse of the denominator in the weights [T4 R-2 L-4 ~> Pa-2].
  real :: I_pterm    ! The inverse of p0 plus p_ave [T2 R-1 L-2 ~> Pa-1].
  real :: hWt_LL, hWt_LR ! hWt_LA is the weighted influence of A on the left column [nondim].
  real :: hWt_RL, hWt_RR ! hWt_RA is the weighted influence of A on the right column [nondim].
  real :: wt_L, wt_R ! The linear weights of the left and right columns [nondim].
  real :: wtT_L, wtT_R ! The weights for tracers from the left and right columns [nondim].
  real :: intp(5)    ! The integrals of specific volume with pressure at the
  real :: a1s        ! Partly rescaled version of a1 [m3 kg-1 C-1 ~> m3 kg-1 degC-1]
  real :: a2s        ! Partly rescaled version of a2 [m3 kg-1 S-1 ~> m3 kg-1 PSU-1]
  real :: b1s        ! Partly rescaled version of b1 [Pa C-1 ~> Pa degC-1]
  real :: b2s        ! Partly rescaled version of b2 [Pa C-2 ~> Pa degC-2]
  real :: b3s        ! Partly rescaled version of b3 [Pa C-3 ~> Pa degC-3]
  real :: b4s        ! Partly rescaled version of b4 [Pa S-1 ~> Pa PSU-1]
  real :: b5s        ! Partly rescaled version of b5 [Pa C-1 S-1 ~> Pa degC-1 PSU-1]
  real :: c1s        ! Partly rescaled version of c1 [m2 s-2 C-1 ~> m2 s-2 degC-1]
  real :: c2s        ! Partly rescaled version of c2 [m2 s-2 C-2 ~> m2 s-2 degC-2]
  real :: c3s        ! Partly rescaled version of c3 [m2 s-2 C-3 ~> m2 s-2 degC-3]
  real :: c4s        ! Partly rescaled version of c4 [m2 s-2 S-1 ~> m2 s-2 PSU-1]
  real :: c5s        ! Partly rescaled version of c5 [m2 s-2 C-1 S-1 ~> m2 s-2 degC-1 PSU-1]
  logical :: do_massWeight ! Indicates whether to do mass weighting.
  logical :: top_massWeight ! Indicates whether to do mass weighting the sea surface
  logical :: massWeight_bug ! If true, use an incorrect expression to determine where to apply mass weighting
  real, parameter :: C1_3 = 1.0/3.0, C1_7 = 1.0/7.0    ! Rational constants [nondim]
  real, parameter :: C1_9 = 1.0/9.0, C1_90 = 1.0/90.0  ! Rational constants [nondim]
  integer :: Isq, Ieq, Jsq, Jeq, ish, ieh, jsh, jeh, i, j, m, halo
  Isq = HI%IscB ; Ieq = HI%IecB ; Jsq = HI%JscB ; Jeq = HI%JecB
  halo = 0 ; if (present(halo_size)) halo = MAX(halo_size,0)
  ish = HI%isc-halo ; ieh = HI%iec+halo ; jsh = HI%jsc-halo ; jeh = HI%jec+halo
  if (present(intx_dza)) then ; ish = MIN(Isq,ish) ; ieh = MAX(Ieq+1,ieh) ; endif
  if (present(inty_dza)) then ; jsh = MIN(Jsq,jsh) ; jeh = MAX(Jeq+1,jeh) ; endif


  al0_scale = 1.0 ; if (present(SV_scale)) al0_scale = SV_scale
  p0_scale = 1.0
  if (present(pres_scale)) then ; if (pres_scale /= 1.0) then
    p0_scale = 1.0 / pres_scale
  endif ; endif
  lam_scale = al0_scale * p0_scale

  a1s = a1 ; a2s = a2
  b1s = b1 ; b2s = b2 ; b3s = b3 ; b4s = b4 ; b5s = b5
  c1s = c1 ; c2s = c2 ; c3s = c3 ; c4s = c4 ; c5s = c5

  if (present(temp_scale)) then ; if (temp_scale /= 1.0) then
    a1s = a1s * temp_scale
    b1s = b1s * temp_scale    ; b2s = b2s * temp_scale**2
    b3s = b3s * temp_scale**3 ; b5s = b5s * temp_scale
    c1s = c1s * temp_scale    ; c2s = c2s * temp_scale**2
    c3s = c3s * temp_scale**3 ; c5s = c5s * temp_scale
  endif ; endif

  if (present(saln_scale)) then ; if (saln_scale /= 1.0) then
    a2s = a2s * saln_scale
    b4s = b4s * saln_scale ; b5s = b5s * saln_scale
    c4s = c4s * saln_scale ; c5s = c5s * saln_scale
  endif ; endif

  do_massWeight = .false. ; massWeight_bug = .false. ; top_massWeight = .false.
  if (present(MassWghtInterp)) then
    do_massWeight = BTEST(MassWghtInterp, 0) ! True for odd values
    top_massWeight = BTEST(MassWghtInterp, 1) ! True if the 2 bit is set
    massWeight_bug = BTEST(MassWghtInterp, 3) ! True if the 8 bit is set
  endif

  !  alpha(j) = (lambda + al0*(pressure(j) + p0)) / (pressure(j) + p0)
  do j=jsh,jeh ; do i=ish,ieh
    al0_2d(i,j) = al0_scale * ( a0 + (a1s*T(i,j) + a2s*S(i,j)) )
    p0_2d(i,j) = p0_scale * ( b0 + ( b4s*S(i,j) + T(i,j) * (b1s + (T(i,j)*(b2s + b3s*T(i,j)) + b5s*S(i,j))) ) )
    lambda_2d(i,j) = lam_scale * ( c0 + ( c4s*S(i,j) + T(i,j) * (c1s + (T(i,j)*(c2s + c3s*T(i,j)) + c5s*S(i,j))) ) )

    al0 = al0_2d(i,j) ; p0 = p0_2d(i,j) ; lambda = lambda_2d(i,j)
    dp = p_b(i,j) - p_t(i,j)
    p_ave = 0.5*(p_t(i,j)+p_b(i,j))
    I_pterm = 1.0 / (p0 + p_ave)

    eps = 0.5 * dp * I_pterm ; eps2 = eps*eps
    alpha_anom = (al0 - spv_ref) + lambda * I_pterm
    rem = (lambda * eps2) * (C1_3 + eps2*(0.2 + eps2*(C1_7 + C1_9*eps2)))
    dza(i,j) = alpha_anom*dp + 2.0*eps*rem
    if (present(intp_dza)) &
      intp_dza(i,j) = 0.5*alpha_anom*dp**2 - dp*((1.0-eps)*rem)
  enddo ; enddo

  if (present(intx_dza)) then ; do j=HI%jsc,HI%jec ; do I=Isq,Ieq
    ! hWght is the distance measure by which the cell is violation of
    ! hydrostatic consistency. For large hWght we bias the interpolation of
    ! T & S along the top and bottom integrals, akin to thickness weighting.
    hWght = 0.0
    if (do_massWeight .and. massWeight_bug) then
      hWght = max(0., bathyP(i,j)-p_t(i+1,j), bathyP(i+1,j)-p_t(i,j))
    elseif (do_massWeight) then
      hWght = max(0., p_t(i+1,j)-bathyP(i,j), p_t(i,j)-bathyP(i+1,j))
    endif
    if (top_massWeight) &
      hWght = max(hWght, P_surf(i,j)-p_b(i+1,j), P_surf(i+1,j)-p_b(i,j))
    if (hWght > 0.) then
      hL = (p_b(i,j) - p_t(i,j)) + dP_neglect
      hR = (p_b(i+1,j) - p_t(i+1,j)) + dP_neglect
      hWght = hWght * ( (hL-hR)/(hL+hR) )**2
      iDenom = 1.0 / ( hWght*(hR + hL) + hL*hR )
      hWt_LL = (hWght*hL + hR*hL) * iDenom ; hWt_LR = (hWght*hR) * iDenom
      hWt_RR = (hWght*hR + hR*hL) * iDenom ; hWt_RL = (hWght*hL) * iDenom
    else
      hWt_LL = 1.0 ; hWt_LR = 0.0 ; hWt_RR = 1.0 ; hWt_RL = 0.0
    endif

    intp(1) = dza(i,j) ; intp(5) = dza(i+1,j)
    do m=2,4
      wt_L = 0.25*real(5-m) ; wt_R = 1.0-wt_L
      wtT_L = (wt_L*hWt_LL) + (wt_R*hWt_RL) ; wtT_R = (wt_L*hWt_LR) + (wt_R*hWt_RR)

      ! T, S, and p are interpolated in the horizontal.  The p interpolation
      ! is linear, but for T and S it may be thickness weighted.
      al0 = (wtT_L*al0_2d(i,j)) + (wtT_R*al0_2d(i+1,j))
      p0 = (wtT_L*p0_2d(i,j)) + (wtT_R*p0_2d(i+1,j))
      lambda = (wtT_L*lambda_2d(i,j)) + (wtT_R*lambda_2d(i+1,j))

      dp = (wt_L*(p_b(i,j) - p_t(i,j))) + (wt_R*(p_b(i+1,j) - p_t(i+1,j)))
      p_ave = 0.5*((wt_L*(p_t(i,j)+p_b(i,j))) + (wt_R*(p_t(i+1,j)+p_b(i+1,j))))
      I_pterm = 1.0 / (p0 + p_ave)

      eps = 0.5 * dp * I_pterm ; eps2 = eps*eps
      intp(m) = ((al0 - spv_ref) + lambda * I_pterm)*dp + 2.0*eps* &
               lambda * eps2 * (C1_3 + eps2*(0.2 + eps2*(C1_7 + C1_9*eps2)))
    enddo
    ! Use Boole's rule to integrate the values.
    intx_dza(i,j) = C1_90*(7.0*(intp(1)+intp(5)) + 32.0*(intp(2)+intp(4)) + &
                           12.0*intp(3))
  enddo ; enddo ; endif

  if (present(inty_dza)) then ; do J=Jsq,Jeq ; do i=HI%isc,HI%iec
    ! hWght is the distance measure by which the cell is violation of
    ! hydrostatic consistency. For large hWght we bias the interpolation of
    ! T & S along the top and bottom integrals, akin to thickness weighting.
    hWght = 0.0
    if (do_massWeight .and. massWeight_bug) then
      hWght = max(0., bathyP(i,j)-p_t(i,j+1), bathyP(i,j+1)-p_t(i,j))
    elseif (do_massWeight) then
      hWght = max(0., p_t(i,j+1)-bathyP(i,j), p_t(i,j)-bathyP(i,j+1))
    endif
    if (top_massWeight) &
      hWght = max(hWght, P_surf(i,j)-p_b(i,j+1), P_surf(i,j+1)-p_b(i,j))
    if (hWght > 0.) then
      hL = (p_b(i,j) - p_t(i,j)) + dP_neglect
      hR = (p_b(i,j+1) - p_t(i,j+1)) + dP_neglect
      hWght = hWght * ( (hL-hR)/(hL+hR) )**2
      iDenom = 1.0 / ( hWght*(hR + hL) + hL*hR )
      hWt_LL = (hWght*hL + hR*hL) * iDenom ; hWt_LR = (hWght*hR) * iDenom
      hWt_RR = (hWght*hR + hR*hL) * iDenom ; hWt_RL = (hWght*hL) * iDenom
    else
      hWt_LL = 1.0 ; hWt_LR = 0.0 ; hWt_RR = 1.0 ; hWt_RL = 0.0
    endif

    intp(1) = dza(i,j) ; intp(5) = dza(i,j+1)
    do m=2,4
      wt_L = 0.25*real(5-m) ; wt_R = 1.0-wt_L
      wtT_L = (wt_L*hWt_LL) + (wt_R*hWt_RL) ; wtT_R = (wt_L*hWt_LR) + (wt_R*hWt_RR)

      ! T, S, and p are interpolated in the horizontal.  The p interpolation
      ! is linear, but for T and S it may be thickness weighted.
      al0 = (wt_L*al0_2d(i,j)) + (wt_R*al0_2d(i,j+1))
      p0 = (wt_L*p0_2d(i,j)) + (wt_R*p0_2d(i,j+1))
      lambda = (wt_L*lambda_2d(i,j)) + (wt_R*lambda_2d(i,j+1))

      dp = (wt_L*(p_b(i,j) - p_t(i,j))) + (wt_R*(p_b(i,j+1) - p_t(i,j+1)))
      p_ave = 0.5*((wt_L*(p_t(i,j)+p_b(i,j))) + (wt_R*(p_t(i,j+1)+p_b(i,j+1))))
      I_pterm = 1.0 / (p0 + p_ave)

      eps = 0.5 * dp * I_pterm ; eps2 = eps*eps
      intp(m) = ((al0 - spv_ref) + lambda * I_pterm)*dp + 2.0*eps* &
               lambda * eps2 * (C1_3 + eps2*(0.2 + eps2*(C1_7 + C1_9*eps2)))
    enddo
    ! Use Boole's rule to integrate the values.
    inty_dza(i,j) = C1_90*(7.0*(intp(1)+intp(5)) + 32.0*(intp(2)+intp(4)) + &
                           12.0*intp(3))
  enddo ; enddo ; endif
end procedure int_spec_vol_dp_wright_red
module procedure calculate_density_array_Wright_red
  integer :: j
  if (present(rho_ref)) then
    do j = start, start+npts-1
      rho(j) = density_anomaly_elem_Wright_red(this, T(j), S(j), pressure(j), rho_ref)
    enddo
  else
    do j = start, start+npts-1
      rho(j) = density_elem_Wright_red(this, T(j), S(j), pressure(j))
    enddo
  endif

end procedure calculate_density_array_Wright_red
module procedure calculate_spec_vol_array_Wright_red
  integer :: j
  if (present(spv_ref)) then
    do j = start, start+npts-1
      specvol(j) = spec_vol_anomaly_elem_Wright_red(this, T(j), S(j), pressure(j), spv_ref)
    enddo
  else
    do j = start, start+npts-1
      specvol(j) = spec_vol_elem_Wright_red(this, T(j), S(j), pressure(j) )
    enddo
  endif

end procedure calculate_spec_vol_array_Wright_red
end submodule MOM_EOS_Wright_red_s
