submodule (MOM_EOS_linear) MOM_EOS_linear_s
  implicit none
contains
module procedure density_elem_linear
  density_elem_linear = this%Rho_T0_S0 + this%dRho_dT*T + this%dRho_dS*S + this%dRho_dp*pressure

end procedure density_elem_linear
module procedure density_anomaly_elem_linear
  density_anomaly_elem_linear = &
      (this%Rho_T0_S0 - rho_ref) + ((this%dRho_dT*T + this%dRho_dS*S) + this%dRho_dp*pressure)

end procedure density_anomaly_elem_linear
module procedure spec_vol_elem_linear
  spec_vol_elem_linear = &
      1.0 / ( this%Rho_T0_S0 + ((this%dRho_dT*T + this%dRho_dS*S) + this%dRho_dp*pressure) )

end procedure spec_vol_elem_linear
module procedure spec_vol_anomaly_elem_linear
  spec_vol_anomaly_elem_linear = &
      ((1.0 - this%Rho_T0_S0*spv_ref) - &
        spv_ref*((this%dRho_dT*T + this%dRho_dS*S) + this%dRho_dp*pressure)) / &
      ( this%Rho_T0_S0 + ((this%dRho_dT*T + this%dRho_dS*S) + this%dRho_dp*pressure) )

end procedure spec_vol_anomaly_elem_linear
module procedure calculate_density_derivs_elem_linear
  drho_dT = this%dRho_dT
  drho_dS = this%dRho_dS

end procedure calculate_density_derivs_elem_linear
module procedure calculate_density_second_derivs_elem_linear
  drho_dS_dS = 0.
  drho_dS_dT = 0.
  drho_dT_dT = 0.
  drho_dS_dP = 0.
  drho_dT_dP = 0.

end procedure calculate_density_second_derivs_elem_linear
module procedure calculate_specvol_derivs_elem_linear
  real :: I_rho2  ! The inverse of density squared [m6 kg-2]
  I_rho2 = 1.0 / (this%Rho_T0_S0 + ((this%dRho_dT*T + this%dRho_dS*S) + this%dRho_dp*pressure))**2
  dSV_dT = -this%dRho_dT * I_rho2
  dSV_dS = -this%dRho_dS * I_rho2

end procedure calculate_specvol_derivs_elem_linear
module procedure calculate_compress_elem_linear
  rho = this%Rho_T0_S0 + this%dRho_dT*T + this%dRho_dS*S + this%dRho_dp*pressure
  drho_dp = this%dRho_dp

end procedure calculate_compress_elem_linear
module procedure avg_spec_vol_linear
  real :: eps2        ! The square of a nondimensional ratio [nondim]
  real :: alpha_p_ave ! The specific volume at pressure mid-point [R-1 ~> m3 kg-1]
  real, parameter :: C1_3 = 1.0/3.0, C1_7 = 1.0/7.0, C1_9 = 1.0/9.0 ! Rational constants [nondim]
  integer :: j
  do j=start,start+npts-1
    alpha_p_ave = &
      1.0 / (Rho_T0_S0 + ((dRho_dT*T(j) + dRho_dS*S(j)) + dRho_dp*(p_t(j) + 0.5 * dp(j))))
    eps2 = (0.5 * (dRho_dp * dp(j)) * alpha_p_ave)**2
    SpV_avg(j) = alpha_p_ave * (1.0 + eps2 * (C1_3 + eps2 * (0.2 + eps2 * (C1_7 + C1_9 * eps2))))
  enddo
end procedure avg_spec_vol_linear
module procedure EoS_fit_range_linear
  if (present(T_min)) T_min = -273.0
  if (present(T_max)) T_max = 100.0
  if (present(S_min)) S_min = 0.0
  if (present(S_max)) S_max = 1000.0
  if (present(p_min)) p_min = 0.0
  if (present(p_max)) p_max = 1.0e9

end procedure EoS_fit_range_linear
module procedure set_params_linear
  if (present(Rho_T0_S0)) this%Rho_T0_S0 = Rho_T0_S0
  if (present(dRho_dT)) this%dRho_dT = dRho_dT
  if (present(dRho_dS)) this%dRho_dS = dRho_dS
  if (present(dRho_dp)) this%dRho_dp = dRho_dp

end procedure set_params_linear
module procedure int_density_dz_linear
  real, dimension(HI%isd:HI%ied,HI%jsd:HI%jed) :: z0pres ! The height at which the pressure is zero [Z ~> m]
  real :: rho_anom      ! The density anomaly from rho_ref [R ~> kg m-3].
  real :: raL, raR      ! rho_anom to the left and right [R ~> kg m-3].
  real :: dz, dzL, dzR  ! Layer thicknesses [Z ~> m].
  real :: GxRho      ! The gravitational acceleration times mean ocean density [R L2 Z-1 T-2 ~> Pa m-1]
  real :: p_ave      ! The layer averaged pressure [R L2 T-2 ~> Pa]
  real :: hWght      ! A pressure-thickness below topography [Z ~> m].
  real :: hL, hR     ! Pressure-thicknesses of the columns to the left and right [Z ~> m].
  real :: iDenom     ! The inverse of the denominator in the weights [Z-2 ~> m-2].
  real :: hWt_LL, hWt_LR ! hWt_LA is the weighted influence of A on the left column [nondim].
  real :: hWt_RL, hWt_RR ! hWt_RA is the weighted influence of A on the right column [nondim].
  real :: wt_L, wt_R ! The linear weights of the left and right columns [nondim].
  real :: wtT_L, wtT_R ! The weights for tracers from the left and right columns [nondim].
  real :: intz(5)    ! The integrals of density with height at the
  logical :: do_massWeight ! Indicates whether to do mass weighting.
  logical :: top_massWeight ! Indicates whether to do mass weighting the sea surface
  real, parameter :: C1_6 = 1.0/6.0, C1_90 = 1.0/90.0  ! Rational constants [nondim].
  integer :: is, ie, js, je, Isq, Ieq, Jsq, Jeq, i, j, m
  Isq = HI%IscB ; Ieq = HI%IecB
  Jsq = HI%JscB ; Jeq = HI%JecB
  is = HI%isc ; ie = HI%iec
  js = HI%jsc ; je = HI%jec

  GxRho = G_e * rho_0

  if (present(Z_0p)) then
    do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
      z0pres(i,j) = Z_0p(i,j)
    enddo ; enddo
  else
    z0pres(:,:) = 0.0
  endif

  do_massWeight = .false. ; top_massWeight = .false.
  if (present(MassWghtInterp)) then
    do_massWeight = BTEST(MassWghtInterp, 0) ! True for odd values
    top_massWeight = BTEST(MassWghtInterp, 1) ! True if the 2 bit is set
  endif

  do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
    dz = z_t(i,j) - z_b(i,j)
    p_ave = -GxRho * (0.5 * (z_t(i,j) + z_b(i,j)) - z0pres(i,j))
    rho_anom = (Rho_T0_S0 - rho_ref) + dRho_dT * T(i,j) + dRho_dS * S(i,j) + dRho_dp * p_ave
    dpa(i,j) = G_e * rho_anom * dz
    if (present(intz_dpa)) &
      intz_dpa(i,j) = 0.5 * G_e * (rho_anom - C1_6 * dRho_dp * (GxRho * dz)) * dz**2
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

    if (hWght <= 0.0) then
      dzL = z_t(i,j) - z_b(i,j) ; dzR = z_t(i+1,j) - z_b(i+1,j)

      p_ave = -GxRho * (0.5 * (z_t(i,j) + z_b(i,j)) - z0pres(i,j))
      raL = (Rho_T0_S0 - rho_ref) + ((dRho_dT*T(i,j) + dRho_dS*S(i,j)) + dRho_dp*p_ave)

      p_ave = -GxRho * (0.5 * (z_t(i+1,j) + z_b(i+1,j)) - z0pres(i+1,j))
      raR = (Rho_T0_S0 - rho_ref) + ((dRho_dT*T(i+1,j) + dRho_dS*S(i+1,j)) + dRho_dp*p_ave)

      intx_dpa(i,j) = G_e*C1_6 * ((dzL*(2.0*raL + raR)) + (dzR*(2.0*raR + raL)))
    else
      hL = (z_t(i,j) - z_b(i,j)) + dz_neglect
      hR = (z_t(i+1,j) - z_b(i+1,j)) + dz_neglect
      hWght = hWght * ( (hL-hR)/(hL+hR) )**2
      iDenom = 1.0 / ( hWght*(hR + hL) + hL*hR )
      hWt_LL = (hWght*hL + hR*hL) * iDenom ; hWt_LR = (hWght*hR) * iDenom
      hWt_RR = (hWght*hR + hR*hL) * iDenom ; hWt_RL = (hWght*hL) * iDenom

      intz(1) = dpa(i,j) ; intz(5) = dpa(i+1,j)
      do m=2,4
        wt_L = 0.25*real(5-m) ; wt_R = 1.0-wt_L
        wtT_L = (wt_L*hWt_LL) + (wt_R*hWt_RL) ; wtT_R = (wt_L*hWt_LR) + (wt_R*hWt_RR)

        dz = (wt_L*(z_t(i,j) - z_b(i,j))) + (wt_R*(z_t(i+1,j) - z_b(i+1,j)))
        p_ave = -GxRho * ((wt_L * (0.5 * (z_t(i,j) + z_b(i,j)) - z0pres(i,j))) + &
                          (wt_R * (0.5 * (z_t(i+1,j) + z_b(i+1,j)) - z0pres(i+1,j))))
        rho_anom = (Rho_T0_S0 - rho_ref) + &
                   ((dRho_dT * ((wtT_L*T(i,j)) + (wtT_R*T(i+1,j))) + &
                     dRho_dS * ((wtT_L*S(i,j)) + (wtT_R*S(i+1,j)))) + dRho_dp * p_ave)
        intz(m) = G_e*rho_anom*dz
      enddo
      ! Use Boole's rule to integrate the values.
      intx_dpa(i,j) = C1_90*(7.0*(intz(1)+intz(5)) + 32.0*(intz(2)+intz(4)) + &
                             12.0*intz(3))
    endif
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

    if (hWght <= 0.0) then
      dzL = z_t(i,j) - z_b(i,j) ; dzR = z_t(i,j+1) - z_b(i,j+1)

      p_ave = -GxRho * (0.5 * (z_t(i,j) + z_b(i,j)) - z0pres(i,j))
      raL = (Rho_T0_S0 - rho_ref) + ((dRho_dT*T(i,j) + dRho_dS*S(i,j)) + dRho_dp*p_ave)

      p_ave = -GxRho * (0.5 * (z_t(i,j+1) + z_b(i,j+1)) - z0pres(i,j+1))
      raR = (Rho_T0_S0 - rho_ref) + ((dRho_dT*T(i,j+1) + dRho_dS*S(i,j+1)) + dRho_dp*p_ave)

      inty_dpa(i,j) = G_e*C1_6 * ((dzL*(2.0*raL + raR)) + (dzR*(2.0*raR + raL)))
    else
      hL = (z_t(i,j) - z_b(i,j)) + dz_neglect
      hR = (z_t(i,j+1) - z_b(i,j+1)) + dz_neglect
      hWght = hWght * ( (hL-hR)/(hL+hR) )**2
      iDenom = 1.0 / ( hWght*(hR + hL) + hL*hR )
      hWt_LL = (hWght*hL + hR*hL) * iDenom ; hWt_LR = (hWght*hR) * iDenom
      hWt_RR = (hWght*hR + hR*hL) * iDenom ; hWt_RL = (hWght*hL) * iDenom

      intz(1) = dpa(i,j) ; intz(5) = dpa(i,j+1)
      do m=2,4
        wt_L = 0.25*real(5-m) ; wt_R = 1.0-wt_L
        wtT_L = (wt_L*hWt_LL) + (wt_R*hWt_RL) ; wtT_R = (wt_L*hWt_LR) + (wt_R*hWt_RR)

        dz = (wt_L*(z_t(i,j) - z_b(i,j))) + (wt_R*(z_t(i,j+1) - z_b(i,j+1)))
        p_ave = -GxRho * ((wt_L * (0.5 * (z_t(i,j) + z_b(i,j)) - z0pres(i,j))) + &
                          (wt_R * (0.5 * (z_t(i,j+1) + z_b(i,j+1)) - z0pres(i,j+1))))
        rho_anom = (Rho_T0_S0 - rho_ref) + &
                   ((dRho_dT * ((wtT_L*T(i,j)) + (wtT_R*T(i,j+1))) + &
                     dRho_dS * ((wtT_L*S(i,j)) + (wtT_R*S(i,j+1)))) + dRho_dp * p_ave)
        intz(m) = G_e*rho_anom*dz
      enddo
      ! Use Boole's rule to integrate the values.
      inty_dpa(i,j) = C1_90*(7.0*(intz(1)+intz(5)) + 32.0*(intz(2)+intz(4)) + &
                             12.0*intz(3))
    endif

  enddo ; enddo ; endif
end procedure int_density_dz_linear
module procedure int_spec_vol_dp_linear
  real :: dRho          ! The density anomaly due to T, S and p [R ~> kg m-3]
  real :: lambda        ! The sound speed squared [L2 T-2 ~> m2 s-2]
  real :: eps, eps2     ! A nondimensional ratio and its square [nondim]
  real :: rem           ! [L2 T-2 ~> m2 s-2]
  real :: p_ave         ! The layer averaged pressure [R L2 T-2 ~> Pa]
  real :: alpha_p_ave   ! The specific volume at p_ave [R-1 ~> m3 kg-1]
  real :: alpha_anom    ! The specific volume anomaly from 1/rho_ref [R-1 ~> m3 kg-1]
  real :: aaL, aaR      ! The specific volume anomaly to the left and right [R-1 ~> m3 kg-1]
  real :: dp, dpL, dpR  ! Layer pressure thicknesses [R L2 T-2 ~> Pa]
  real :: hWght      ! A pressure-thickness below topography [R L2 T-2 ~> Pa]
  real :: hL, hR     ! Pressure-thicknesses of the columns to the left and right [R L2 T-2 ~> Pa]
  real :: iDenom     ! The inverse of the denominator in the weights [T4 R-2 L-4 ~> Pa-2]
  real :: hWt_LL, hWt_LR ! hWt_LA is the weighted influence of A on the left column [nondim].
  real :: hWt_RL, hWt_RR ! hWt_RA is the weighted influence of A on the right column [nondim].
  real :: wt_L, wt_R ! The linear weights of the left and right columns [nondim].
  real :: wtT_L, wtT_R ! The weights for tracers from the left and right columns [nondim].
  real :: intp(5)    ! The integrals of specific volume with pressure at the
  logical :: do_massWeight ! Indicates whether to do mass weighting.
  logical :: top_massWeight ! Indicates whether to do mass weighting the sea surface
  logical :: massWeight_bug ! If true, use an incorrect expression to determine where to apply mass weighting
  real, parameter :: C1_3 = 1.0/3.0, C1_7 = 1.0/7.0, C1_9 = 1.0/9.0  ! Rational constants [nondim]
  real, parameter :: C1_6 = 1.0/6.0, C1_90 = 1.0/90.0  ! Rational constants [nondim].
  integer :: Isq, Ieq, Jsq, Jeq, ish, ieh, jsh, jeh, i, j, m, halo
  Isq = HI%IscB ; Ieq = HI%IecB ; Jsq = HI%JscB ; Jeq = HI%JecB
  halo = 0 ; if (present(halo_size)) halo = MAX(halo_size,0)
  ish = HI%isc-halo ; ieh = HI%iec+halo ; jsh = HI%jsc-halo ; jeh = HI%jec+halo
  if (present(intx_dza)) then ; ish = MIN(Isq,ish) ; ieh = MAX(Ieq+1,ieh) ; endif
  if (present(inty_dza)) then ; jsh = MIN(Jsq,jsh) ; jeh = MAX(Jeq+1,jeh) ; endif

  do_massWeight = .false. ; massWeight_bug = .false. ; top_massWeight = .false.
  if (present(MassWghtInterp)) then
    do_massWeight = BTEST(MassWghtInterp, 0) ! True for odd values
    top_massWeight = BTEST(MassWghtInterp, 1) ! True if the 2 bit is set
    massWeight_bug = BTEST(MassWghtInterp, 3) ! True if the 8 bit is set
  endif

  lambda = 0.0 ; if (dRho_dp/=0.0) lambda = 1.0 / dRho_dp
  do j=jsh,jeh ; do i=ish,ieh
    dp = p_b(i,j) - p_t(i,j)
    p_ave = 0.5 * (p_t(i,j) + p_b(i,j))

    drho = (dRho_dT * T(i,j) + dRho_dS * S(i,j)) + dRho_dp * p_ave
    alpha_p_ave = 1.0 / (Rho_T0_S0 + drho)

    ! A realistic upbound of eps is ~0.02, using dRho_dp ~ (1500 m/s)**(-2), alpha_p_ave ~ 1/(1030 kg/m3)
    ! and dp ~ 1e8 Pa [~dz=10000m]. And if we use dp ~ 1e6 [~dz=100m], eps ~ 2e-4.
    ! Analytically dza = 1/dRho_dp * ln[(1+eps)/(1-eps)] - alpha_ref * dp, and the expression here gives the first
    ! five terms from its Taylor series with a truncation error of O(eps**11), which is beyond double floating
    ! point precision.
    eps = 0.5 * (dRho_dp * dp) * alpha_p_ave ; eps2 = eps * eps
    ! alpha_anom = 1.0/(Rho_T0_S0 + dRho) - alpha_ref
    alpha_anom = ((1.0 - Rho_T0_S0 * alpha_ref) - drho * alpha_ref) / (Rho_T0_S0 + drho)
    ! The following expression would be more efficient but I suspect it changes answer.
    ! alpha_anom = ((1.0 - Rho_T0_S0 * alpha_ref) - drho * alpha_ref) * alpha_p_ave
    rem = (lambda * eps2) * (C1_3 + eps2 * (0.2 + eps2 * (C1_7 + C1_9 * eps2)))
    dza(i,j) = alpha_anom * dp + 2.0 * eps * rem
    if (present(intp_dza)) &
      intp_dza(i,j) = 0.5 * alpha_anom * dp**2 - dp * ((1.0 - eps) * rem)
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

    if (hWght <= 0.0) then
      dpL = p_b(i,j) - p_t(i,j) ; dpR = p_b(i+1,j) - p_t(i+1,j)

      p_ave = 0.5 * (p_b(i,j) + p_t(i,j))
      drho = (dRho_dT*T(i,j) + dRho_dS*S(i,j)) + dRho_dp * p_ave
      aaL = ((1.0 - Rho_T0_S0*alpha_ref) - drho*alpha_ref) / (Rho_T0_S0 + drho)

      p_ave = 0.5 * (p_b(i+1,j) + p_t(i+1,j))
      drho = (dRho_dT*T(i+1,j) + dRho_dS*S(i+1,j)) + dRho_dp * p_ave
      aaR = ((1.0 - Rho_T0_S0*alpha_ref) - drho*alpha_ref) / (Rho_T0_S0 + drho)

      intx_dza(i,j) = C1_6 * (2.0*((dpL*aaL) + (dpR*aaR)) + ((dpL*aaR) + (dpR*aaL)))
    else
      hL = (p_b(i,j) - p_t(i,j)) + dP_neglect
      hR = (p_b(i+1,j) - p_t(i+1,j)) + dP_neglect
      hWght = hWght * ( (hL-hR)/(hL+hR) )**2
      iDenom = 1.0 / ( hWght*(hR + hL) + hL*hR )
      hWt_LL = (hWght*hL + hR*hL) * iDenom ; hWt_LR = (hWght*hR) * iDenom
      hWt_RR = (hWght*hR + hR*hL) * iDenom ; hWt_RL = (hWght*hL) * iDenom

      intp(1) = dza(i,j) ; intp(5) = dza(i+1,j)
      do m=2,4
        wt_L = 0.25*real(5-m) ; wt_R = 1.0-wt_L
        wtT_L = (wt_L*hWt_LL) + (wt_R*hWt_RL) ; wtT_R = (wt_L*hWt_LR) + (wt_R*hWt_RR)

        ! T, S, and p are interpolated in the horizontal.  The p interpolation
        ! is linear, but for T and S it may be thickness weighted.
        dp = (wt_L*(p_b(i,j) - p_t(i,j))) + (wt_R*(p_b(i+1,j) - p_t(i+1,j)))
        p_ave = 0.5*((wt_L*(p_t(i,j)+p_b(i,j))) + (wt_R*(p_t(i+1,j)+p_b(i+1,j))))

        drho = (dRho_dT*((wtT_L*T(i,j)) + (wtT_R*T(i+1,j))) + &
                dRho_dS*((wtT_L*S(i,j)) + (wtT_R*S(i+1,j)))) + dRho_dp * p_ave
        ! alpha_anom = 1.0/(Rho_T0_S0  + drho)) - alpha_ref
        alpha_anom = ((1.0-Rho_T0_S0*alpha_ref) - drho*alpha_ref) / (Rho_T0_S0 + drho)
        intp(m) = alpha_anom*dp
      enddo
      ! Use Boole's rule to integrate the interface height anomaly values in y.
      intx_dza(i,j) = C1_90*(7.0*(intp(1)+intp(5)) + 32.0*(intp(2)+intp(4)) + &
                             12.0*intp(3))
    endif
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

    if (hWght <= 0.0) then
      dpL = p_b(i,j) - p_t(i,j) ; dpR = p_b(i,j+1) - p_t(i,j+1)

      p_ave = 0.5 * (p_b(i,j) + p_t(i,j)) + dRho_dp * p_ave
      drho = (dRho_dT*T(i,j) + dRho_dS*S(i,j)) + dRho_dp * p_ave
      aaL = ((1.0 - Rho_T0_S0*alpha_ref) - drho*alpha_ref) / (Rho_T0_S0 + drho)

      p_ave = 0.5 * (p_b(i,j+1) + p_t(i,j+1)) + dRho_dp * p_ave
      drho = (dRho_dT*T(i,j+1) + dRho_dS*S(i,j+1)) + dRho_dp * p_ave
      aaR = ((1.0 - Rho_T0_S0*alpha_ref) - drho*alpha_ref) / (Rho_T0_S0 + drho)

      inty_dza(i,j) = C1_6 * (2.0*((dpL*aaL) + (dpR*aaR)) + ((dpL*aaR) + (dpR*aaL)))
    else
      hL = (p_b(i,j) - p_t(i,j)) + dP_neglect
      hR = (p_b(i,j+1) - p_t(i,j+1)) + dP_neglect
      hWght = hWght * ( (hL-hR)/(hL+hR) )**2
      iDenom = 1.0 / ( hWght*(hR + hL) + hL*hR )
      hWt_LL = (hWght*hL + hR*hL) * iDenom ; hWt_LR = (hWght*hR) * iDenom
      hWt_RR = (hWght*hR + hR*hL) * iDenom ; hWt_RL = (hWght*hL) * iDenom

      intp(1) = dza(i,j) ; intp(5) = dza(i,j+1)
      do m=2,4
        wt_L = 0.25*real(5-m) ; wt_R = 1.0-wt_L
        wtT_L = (wt_L*hWt_LL) + (wt_R*hWt_RL) ; wtT_R = (wt_L*hWt_LR) + (wt_R*hWt_RR)

        ! T, S, and p are interpolated in the horizontal.  The p interpolation
        ! is linear, but for T and S it may be thickness weighted.
        dp = (wt_L*(p_b(i,j) - p_t(i,j))) + (wt_R*(p_b(i,j+1) - p_t(i,j+1)))
        p_ave = 0.5*((wt_L*(p_t(i,j)+p_b(i,j))) + (wt_R*(p_t(i,j+1)+p_b(i,j+1))))

        drho = (dRho_dT*((wtT_L*T(i,j)) + (wtT_R*T(i,j+1))) + &
                dRho_dS*((wtT_L*S(i,j)) + (wtT_R*S(i,j+1)))) + dRho_dp * p_ave
        ! alpha_anom = 1.0/(Rho_T0_S0  + drho)) - alpha_ref
        alpha_anom = ((1.0-Rho_T0_S0*alpha_ref) - drho*alpha_ref) / (Rho_T0_S0 + drho)
        intp(m) = alpha_anom*dp
      enddo
      ! Use Boole's rule to integrate the interface height anomaly values in y.
      inty_dza(i,j) = C1_90*(7.0*(intp(1)+intp(5)) + 32.0*(intp(2)+intp(4)) + &
                             12.0*intp(3))
    endif
  enddo ; enddo ; endif
end procedure int_spec_vol_dp_linear
module procedure calculate_density_array_linear
  integer :: j
  if (present(rho_ref)) then
    do j = start, start+npts-1
      rho(j) = density_anomaly_elem_linear(this, T(j), S(j), pressure(j), rho_ref)
    enddo
  else
    do j = start, start+npts-1
      rho(j) = density_elem_linear(this, T(j), S(j), pressure(j))
    enddo
  endif

end procedure calculate_density_array_linear
module procedure calculate_spec_vol_array_linear
  integer :: j
  if (present(spv_ref)) then
    do j = start, start+npts-1
      specvol(j) = spec_vol_anomaly_elem_linear(this, T(j), S(j), pressure(j), spv_ref)
    enddo
  else
    do j = start, start+npts-1
      specvol(j) = spec_vol_elem_linear(this, T(j), S(j), pressure(j) )
    enddo
  endif

end procedure calculate_spec_vol_array_linear
end submodule MOM_EOS_linear_s
