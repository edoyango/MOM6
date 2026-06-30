submodule (MOM_ice_shelf) MOM_ice_shelf_s
#include <MOM_memory.h>
#ifdef SYMMETRIC_MEMORY_
#  define GRID_SYM_ .true.
#else
#  define GRID_SYM_ .false.
#endif
  implicit none
contains
module procedure shelf_calc_flux
  type(ocean_grid_type), pointer :: G => NULL()  !< The grid structure used by the ice shelf.
  type(unit_scale_type), pointer :: US => NULL() !< Pointer to a structure containing
  type(ice_shelf_state), pointer :: ISS => NULL() !< A structure with elements that describe
  type(surface), pointer :: sfc_state => NULL()
  type(forcing), pointer :: fluxes => NULL()
  real, dimension(SZI_(CS%grid)) :: &
    Rhoml, &   !< Ocean mixed layer density [R ~> kg m-3].
    dR0_dT, &  !< Partial derivative of the mixed layer density
               !< with temperature [R C-1 ~> kg m-3 degC-1].
    dR0_dS, &  !< Partial derivative of the mixed layer density
               !< with salinity [R S-1 ~> kg m-3 ppt-1].
    p_int      !< The pressure at the ice-ocean interface [R L2 T-2 ~> Pa].

  real, dimension(SZI_(CS%grid),SZJ_(CS%grid)) :: &
    exch_vel_t, &   !< Sub-shelf thermal exchange velocity [Z T-1 ~> m s-1]
    exch_vel_s, &   !< Sub-shelf salt exchange velocity [Z T-1 ~> m s-1]
    dh_bdott, & !< Basal melt/accumulation over a time step, used for diagnostics [Z ~> m]
    dh_adott    !< Surface melt/accumulation over a time step, used for diagnostics [Z ~> m]
  real, dimension(SZDI_(CS%grid),SZDJ_(CS%grid)) :: &
    mass_flux  !< Total mass flux of freshwater across the ice-ocean interface. [R Z L2 T-1 ~> kg s-1]
  real, dimension(SZDI_(CS%grid),SZDJ_(CS%grid)) :: &
    haline_driving !< (SSS - S_boundary) ice-ocean
               !! interface, positive for melting and negative for freezing [S ~> ppt].
               !! This is computed as part of the ISOMIP diagnostics.
  real :: time_step !< Length of time over which these fluxes will be applied [T ~> s].
  real :: Itime_step !< Inverse of the length of time over which these fluxes will be applied [T-1 ~> s-1]
  real :: VK       !< Von Karman's constant [nondim]
  real :: ZETA_N   !< This is the stability constant xi_N = 0.052 from Holland & Jenkins '99
                   !! divided by the von Karman constant VK. Was 1/8. [nondim]
  real :: Rf_crit  !< critical flux Richardson number  [nondim]
  real :: I_2Zeta_N !< Half the inverse of Zeta_N [nondim].
  real :: I_LF     !< The inverse of the latent heat of fusion [Q-1 ~> kg J-1].
  real :: I_dt_LHF  ! The inverse of the timestep times the latent heat of fusion [Q-1 T-1 ~> kg J-1 s-1].
  real :: I_VK     !< The inverse of the Von Karman constant [nondim].
  real :: PR, SC   !< The Prandtl number and Schmidt number [nondim].

  ! 3 equations formulation variables
  real, dimension(SZDI_(CS%grid),SZDJ_(CS%grid)) :: &
    Sbdry     !< Salinities in the ocean at the interface with the ice shelf [S ~> ppt].
  real :: Sbdry_it ! The boundary salinity at an iteration [S ~> ppt]
  real :: S_a      ! A variable used to find salt roots [S-1 ~> ppt-1]
  real :: S_b      ! A variable used to find salt roots [nondim]
  real :: S_c      ! A variable used to find salt roots [S ~> ppt]
  real :: dS_it    !< The interface salinity change during an iteration [S ~> ppt].
  real :: hBL_neut !< The neutral boundary layer thickness [Z ~> m].
  real :: hBL_neut_h_molec !< The ratio of the neutral boundary layer thickness
                   !! to the molecular boundary layer thickness [nondim].
  real :: wT_flux !< The downward vertical flux of heat just inside the ocean [C Z T-1 ~> degC m s-1].
  real :: wB_flux !< The downward vertical flux of buoyancy just inside the ocean [Z2 T-3 ~> m2 s-3].
  real :: dB_dS   !< The derivative of buoyancy with salinity [Z T-2 S-1 ~> m s-2 ppt-1].
  real :: dB_dT   !< The derivative of buoyancy with temperature [Z T-2 C-1 ~> m s-2 degC-1].
  real :: I_n_star ! The inverse of the ratio of working boundary layer thickness
                   ! to the neutral thickness [nondim]
  real :: n_star_term ! A term in the expression for nstar [T3 Z-2 ~> s3 m-2]
  real :: absf     ! The absolute value of the Coriolis parameter [T-1 ~> s-1]
  real :: dIns_dwB !< The partial derivative of I_n_star with wB_flux, in [T3 Z-2 ~> s3 m-2]
  real :: dT_ustar ! The difference between the freezing point and the ocean boundary layer
                   ! temperature times the friction velocity [C Z T-1 ~> degC m s-1]
  real :: dS_ustar ! The difference between the salinity at the ice-ocean interface and the ocean
                   ! boundary layer salinity times the friction velocity [S Z T-1 ~> ppt m s-1]
  real :: ustar_h  ! The friction velocity in the water below the ice shelf [Z T-1 ~> m s-1]
  real :: Gam_turb ! A relative turbluent diffusivity [nondim]
  real :: Gam_mol_t, Gam_mol_s ! Relative coefficients of molecular diffusivities [nondim]
  real :: RhoCp     ! A typical ocean density times the heat capacity of water [Q R C-1 ~> J m-3 degC-1]
  real :: ln_neut   ! The log of the ratio of the neutral boundary layer thickness to the molecular
                    ! boundary layer thickness if it is greater than 1 or 0 otherwise [nondim]
  real :: mass_exch ! A mass exchange rate [R Z T-1 ~> kg m-2 s-1]
  real :: Sb_min, Sb_max ! Minimum and maximum boundary salinities [S ~> ppt]
  real :: dS_min, dS_max ! Minimum and maximum salinity changes [S ~> ppt]
  ! Variables used in iterating for wB_flux.
  real :: wB_flux_next ! The next interation's guess for wB_flux [Z2 T-3 ~> m2 s-3]
  real :: wB_flux_new  ! An updated value of wB_flux when Gam_turb is based on wB_flux [Z2 T-3 ~> m2 s-3]
  real :: wB_flux_max  ! The upper bound on wB_flux [Z2 T-3 ~> m2 s-3]
  real :: wB_flux_min  ! The lower bound on wB_flux [Z2 T-3 ~> m2 s-3]
  real :: dDwB_dwB     ! The slope of the change in wB_flux between iterations with wB_flux [nondim]
  real :: DwB_max      ! The change in wB_flux when it is wB_flux_max [Z2 T-3 ~> m2 s-3]
  real :: DwB_min      ! The change in wB_flux when it is wB_flux_min [Z2 T-3 ~> m2 s-3]
  real :: I_Gam_T, I_Gam_S  ! Terms that vary inversely with Gam_mol_T or Gam_mol_S and Gam_turb [nondim]
  real :: dG_dwB       ! The derivative of Gam_turb with wB [T3 Z-2 ~> s3 m-2]
  real :: taux2, tauy2 ! The squared surface stresses [R2 L2 Z2 T-4 ~> Pa2].
  real :: u2_av, v2_av ! The ice-area weighted average squared ocean velocities [L2 T-2 ~> m2 s-2]
  real :: asu1, asu2   ! Ocean areas covered by ice shelves at neighboring u-points [L2 ~> m2]
  real :: asv1, asv2   ! Ocean areas covered by ice shelves at neighboring v-points [L2 ~> m2]
  real :: I_au, I_av   ! The Adcroft reciprocals of the ice shelf areas at adjacent points [L-2 ~> m-2]
  real :: Irho0        ! The inverse of the mean density times a unit conversion factor [R-1 L Z-1 ~> m3 kg-1]
  logical :: Sb_min_set, Sb_max_set
  logical :: root_found
  logical :: update_ice_vel ! If true, it is time to update the ice shelf velocities.
  logical :: coupled_GL     ! If true, the grounding line position is determined based on
                            ! coupled ice-ocean dynamics.
  logical :: add_frazil ! If true, allow frazil formation to modify ice-shelf water flux
  real, parameter :: c2_3 = 2.0/3.0 ! Two thirds [nondim]
  character(len=320) :: mesg  ! The text of an error message
  integer, dimension(2) :: EOSdom ! The i-computational domain for the equation of state
  integer :: i, j, is, ie, js, je, ied, jed, it1, it3
  real :: vaf0, vaf0_A, vaf0_G ! The previous volumes above floatation [Z L2 ~> m3]
                               ! for all ice sheets, Antarctica only, or Greenland only

  if (.not. associated(CS)) call MOM_error(FATAL, "shelf_calc_flux: "// &
       "initialize_ice_shelf must be called before shelf_calc_flux.")
  call cpu_clock_begin(id_clock_shelf)

  G => CS%grid ; US => CS%US
  ISS => CS%ISS
  time_step = time_step_in
  Itime_step = 1./time_step

  dh_adott(:,:) = 0.0 ; dh_bdott(:,:) = 0.0

  if (CS%active_shelf_dynamics) then
    !calculate previous volumes above floatation
    if (CS%id_dvafdt     > 0) call volume_above_floatation(CS%dCS, G, ISS, vaf0)                 !all ice sheet
    if (CS%id_Ant_dvafdt > 0) call volume_above_floatation(CS%dCS, G, ISS, vaf0_A, hemisphere=0) !Antarctica only
    if (CS%id_Gr_dvafdt  > 0) call volume_above_floatation(CS%dCS, G, ISS, vaf0_G, hemisphere=1) !Greenland only
  endif

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; ied = G%ied ; jed = G%jed
  if (CS%data_override_shelf_fluxes .and. CS%active_shelf_dynamics) then
    call data_override(G%Domain, 'shelf_sfc_mass_flux', fluxes_in%shelf_sfc_mass_flux(is:ie,js:je), CS%Time, &
                       scale=US%kg_m2s_to_RZ_T)
    call pass_var(fluxes_in%shelf_sfc_mass_flux, G%domain, complete=.true.)
  endif

  if (CS%rotate_index) then
    allocate(sfc_state)
    call rotate_surface_state(sfc_state_in, sfc_state, CS%Grid, CS%turns)
    allocate(fluxes)
    call allocate_forcing_type(fluxes_in, G, fluxes, turns=CS%turns)
    call rotate_forcing(fluxes_in, fluxes, CS%turns)
  else
    sfc_state => sfc_state_in
    fluxes => fluxes_in
  endif
  ! useful parameters
  ZETA_N = CS%Zeta_N
  VK = CS%Vk
  Rf_crit = CS%Rc
  I_2Zeta_N = 0.5 / CS%Zeta_N
  I_LF = 1.0 / CS%Lat_fusion
  I_dt_LHF = 1.0 / (time_step * CS%Lat_fusion)
  SC = CS%kv_molec/CS%kd_molec_salt
  PR = CS%kv_molec/CS%kd_molec_temp
  I_VK = 1.0/VK
  RhoCp = CS%Rho_ocn * CS%Cp

  !first calculate molecular component
  Gam_mol_t = 12.5 * (PR**c2_3) - 6.0
  Gam_mol_s = 12.5 * (SC**c2_3) - 6.0

  ! GMM, zero some fields of the ice shelf structure (ice_shelf_CS)
  ! these fields are already set to zero during initialization
  ! However, they seem to be changed somewhere and, for diagnostic
  ! reasons, it is better to set them to zero again.
  exch_vel_t(:,:) = 0.0 ; exch_vel_s(:,:) = 0.0
  ISS%tflux_shelf(:,:) = 0.0 ; ISS%water_flux(:,:) = 0.0
  ISS%salt_flux(:,:) = 0.0 ; ISS%tflux_ocn(:,:) = 0.0 ; ISS%tfreeze(:,:) = 0.0
  ! define Sbdry to avoid Run-Time Check Failure, when melt is not computed.
  haline_driving(:,:) = 0.0
  Sbdry(:,:) = sfc_state%sss(:,:)

  !update time
  CS%Time = Time

  if (CS%override_shelf_movement) then
    CS%time_step = time_step
    ! update shelf mass
    if (CS%mass_from_file) call update_shelf_mass(G, US, CS, ISS, Time)
  endif

  if (CS%debug) then
    call hchksum(fluxes_in%frac_shelf_h, "frac_shelf_h before apply melting", CS%Grid_in%HI, haloshift=0)
    call hchksum(sfc_state_in%sst, "sst before apply melting", CS%Grid_in%HI, haloshift=0, unscale=US%C_to_degC)
    call hchksum(sfc_state_in%sss, "sss before apply melting", CS%Grid_in%HI, haloshift=0, unscale=US%S_to_ppt)
    call uvchksum("[uv]_ml before apply melting", sfc_state_in%u, sfc_state_in%v, &
                  CS%Grid_in%HI, haloshift=0, unscale=US%L_T_to_m_s)
    call hchksum(sfc_state_in%ocean_mass, "ocean_mass before apply melting", CS%Grid_in%HI, haloshift=0, &
                 unscale=US%RZ_to_kg_m2)
  endif

  ! Calculate the friction velocity under ice shelves, using taux_shelf and tauy_shelf if possible.
  if (allocated(sfc_state%taux_shelf) .and. allocated(sfc_state%tauy_shelf)) then
    call pass_vector(sfc_state%taux_shelf, sfc_state%tauy_shelf, G%domain, TO_ALL, CGRID_NE)
  endif
  Irho0 = US%Z_to_L / CS%Rho_ocn
  do j=js,je ; do i=is,ie ; if (fluxes%frac_shelf_h(i,j) > 0.0) then
    taux2 = 0.0 ; tauy2 = 0.0 ; u2_av = 0.0 ; v2_av = 0.0
    asu1 = (ISS%area_shelf_h(i-1,j) + ISS%area_shelf_h(i,j))
    asu2 = (ISS%area_shelf_h(i,j) + ISS%area_shelf_h(i+1,j))
    asv1 = (ISS%area_shelf_h(i,j-1) + ISS%area_shelf_h(i,j))
    asv2 = (ISS%area_shelf_h(i,j) + ISS%area_shelf_h(i,j+1))
    I_au = 0.0 ; if (asu1 + asu2 > 0.0) I_au = 1.0 / (asu1 + asu2)
    I_av = 0.0 ; if (asv1 + asv2 > 0.0) I_av = 1.0 / (asv1 + asv2)
    if (allocated(sfc_state%taux_shelf) .and. allocated(sfc_state%tauy_shelf)) then
      taux2 = (((asu1 * (sfc_state%taux_shelf(I-1,j)**2)) + (asu2 * (sfc_state%taux_shelf(I,j)**2))  ) * I_au)
      tauy2 = (((asv1 * (sfc_state%tauy_shelf(i,J-1)**2)) + (asv2 * (sfc_state%tauy_shelf(i,J)**2))  ) * I_av)
    endif
    u2_av = (((asu1 * (sfc_state%u(I-1,j)**2)) + (asu2 * sfc_state%u(I,j)**2)) * I_au)
    if (CS%ustar_from_vel_bugfix) then
      v2_av = (((asv1 * (sfc_state%v(i,J-1)**2)) + (asv2 * sfc_state%v(i,J)**2)) * I_av)
    else
      v2_av = (((asv1 * (sfc_state%v(i,J-1)**2)) + (asu2 * sfc_state%v(i,J)**2)) * I_av)
    endif

    if ((taux2 + tauy2 > 0.0) .and. .not.CS%ustar_shelf_from_vel) then
      if (CS%ustar_max >= 0.0) then
        fluxes%ustar_shelf(i,j) = MIN(CS%ustar_max, MAX(CS%ustar_bg, US%L_to_Z * &
            sqrt(Irho0 * sqrt(taux2 + tauy2) + CS%cdrag*CS%utide(i,j)**2)))
      else
        fluxes%ustar_shelf(i,j) = MAX(CS%ustar_bg, US%L_to_Z * &
            sqrt(Irho0 * sqrt(taux2 + tauy2) + CS%cdrag*CS%utide(i,j)**2))
      endif
    else   ! Take care of the cases when taux_shelf is not set or not allocated.
      fluxes%ustar_shelf(i,j) = MAX(CS%ustar_bg, US%L_TO_Z * &
          sqrt(CS%cdrag*((u2_av + v2_av) + CS%utide(i,j)**2)))
    endif
  else ! There is no shelf here.
    fluxes%ustar_shelf(i,j) = 0.0
  endif ; enddo ; enddo

  EOSdom(:) = EOS_domain(G%HI)
  do j=js,je
    ! Find the pressure at the ice-ocean interface, averaged only over the
    ! part of the cell covered by ice shelf.
    do i=is,ie ; p_int(i) = CS%g_Earth * ISS%mass_shelf(i,j) ; enddo

    ! Calculate insitu densities and expansion coefficients
    call calculate_density(sfc_state%sst(:,j), sfc_state%sss(:,j), p_int, Rhoml(:), &
                                 CS%eqn_of_state, EOSdom)
    call calculate_density_derivs(sfc_state%sst(:,j), sfc_state%sss(:,j), p_int, &
                                  dR0_dT, dR0_dS, CS%eqn_of_state, EOSdom)

    do i=is,ie
      if ((sfc_state%ocean_mass(i,j) > CS%col_mass_melt_threshold) .and. &
          (ISS%area_shelf_h(i,j) > 0.0) .and. CS%isthermo &
           .and. ISS%melt_mask(i,j)>0.0) then

        if (CS%threeeq) then
          !   Iteratively determine a self-consistent set of fluxes, with the ocean
          ! salinity just below the ice-shelf as the variable that is being
          ! iterated for.

          ustar_h = fluxes%ustar_shelf(i,j)

          ! Estimate the neutral ocean boundary layer thickness as the minimum of the
          ! reported ocean mixed layer thickness and the neutral Ekman depth.
          absf = 0.25*((abs(G%CoriolisBu(I,J)) + abs(G%CoriolisBu(I-1,J-1))) + &
                                 (abs(G%CoriolisBu(I,J-1)) + abs(G%CoriolisBu(I-1,J))))
          if (absf*sfc_state%Hml(i,j) <= VK*ustar_h) then ; hBL_neut = sfc_state%Hml(i,j)
          else ; hBL_neut = (VK*ustar_h) / absf ; endif
          hBL_neut_h_molec = ZETA_N * ((hBL_neut * ustar_h) / (5.0 * CS%kv_molec))
          ln_neut = 0.0 ; if (hBL_neut_h_molec > 1.0) ln_neut = log(hBL_neut_h_molec)
          n_star_term = (ZETA_N * hBL_neut * VK) / (Rf_crit * ustar_h**3)

          ! Determine the mixed layer buoyancy flux, wB_flux.
          dB_dS = (US%L_to_Z**2*CS%g_Earth / Rhoml(i)) * dR0_dS(i)
          dB_dT = (US%L_to_Z**2*CS%g_Earth / Rhoml(i)) * dR0_dT(i)

          if (CS%find_salt_root) then
            ! Solve for the skin salinity using the linearized liquidus parameters and
            ! balancing the turbulent fresh water flux in the near-boundary layer with
            ! the net fresh water or salt added by melting:
            ! (Cp/Lat_fusion)*Gamma_T_3Eq*(TFr_skin-T_ocn) = Gamma_S_3Eq*(S_skin-S_ocn)/S_skin

            ! S_a is always < 0.0 with a realistic expression for the freezing point.
            S_a = CS%dTFr_dS * CS%Gamma_T_3EQ * CS%Cp
            S_b = CS%Gamma_T_3EQ*CS%Cp*(CS%TFr_0_0 + CS%dTFr_dp*p_int(i) - sfc_state%sst(i,j)) - &
                  CS%Lat_fusion * CS%Gamma_S_3EQ    ! S_b Can take either sign, but is usually negative.
            S_c = CS%Lat_fusion * CS%Gamma_S_3EQ * sfc_state%sss(i,j) ! Always >= 0

            if (S_c == 0.0) then  ! The solution for fresh water.
              Sbdry(i,j) = 0.0
            elseif (S_a < 0.0) then ! This is the usual ocean case
              if (S_b < 0.0) then ! This is almost always the case
                Sbdry(i,j) = 2.0*S_c / (-S_b + SQRT(S_b*S_b - 4.*S_a*S_c))
              else
                Sbdry(i,j) = (S_b + SQRT(S_b*S_b - 4.*S_a*S_c)) / (-2.*S_a)
              endif
            elseif ((S_a == 0.0) .and. (S_b < 0.0)) then ! It should be the case that S_b < 0.
              Sbdry(i,j) = -S_c / S_b
            else
              call MOM_error(FATAL, "Impossible conditions found in 3-equation skin salinity calculation.")
            endif

            ! Safety check
            if (Sbdry(i,j) < 0.) then
              write(mesg,*) 'sfc_state%sss(i,j) = ',US%S_to_ppt*sfc_state%sss(i,j), &
                            'S_a, S_b, S_c', US%ppt_to_S*S_a, S_b, US%S_to_ppt*S_c
              call MOM_error(WARNING, mesg, .true.)
              call MOM_error(FATAL, "shelf_calc_flux: Negative salinity (Sbdry).")
            endif
          else
            ! Guess sss as the iteration starting point for the boundary salinity.
            Sbdry(i,j) = sfc_state%sss(i,j) ; Sb_max_set = .false.
            Sb_min_set = .false.
          endif !find_salt_root

          do it1 = 1,20
            ! Determine the potential temperature at the ice-ocean interface.
            ! The following two lines are equivalent:
            ! call calculate_TFreeze(Sbdry(i,j), p_int(i), ISS%tfreeze(i,j), CS%eqn_of_state, scale_from_EOS=.true.)
            call calculate_TFreeze(Sbdry(i:i,j), p_int(i:i), ISS%tfreeze(i:i,j), CS%eqn_of_state)

            dT_ustar = (ISS%tfreeze(i,j) - sfc_state%sst(i,j)) * ustar_h
            dS_ustar = (Sbdry(i,j) - sfc_state%sss(i,j)) * ustar_h

            if (CS%const_gamma) then
              ! If using a constant gamma_T, there are no effects of the buoyancy flux on the turbulence.
              I_Gam_T = CS%Gamma_T_3EQ
              I_Gam_S = CS%Gamma_S_3EQ
              wT_flux = dT_ustar * CS%Gamma_T_3EQ
              wB_flux = dB_dS * (dS_ustar * CS%Gamma_S_3EQ) + dB_dT * wT_flux
            elseif (.not.CS%buoy_flux_itt_bugfix) then
              ! Gamma_T and gamma_S are a function of the buoyancy flux, and there should have been
              ! iteration to find the root where wB_flux is consistent with the values of gamma with
              ! that flux, but it was omitted.
              Gam_turb = I_VK * (ln_neut + (I_2Zeta_N - 1.0))
              I_Gam_T = 1.0 / (Gam_mol_t + Gam_turb)
              I_Gam_S = 1.0 / (Gam_mol_s + Gam_turb)
              wB_flux = dB_dS * (dS_ustar * I_Gam_S) + dB_dT * (dT_ustar * I_Gam_T)

              if (wB_flux < 0.0) then  ! The stabilising buoyancy flux reduces the turbulent fluxes.
                I_n_star = sqrt(1.0 - n_star_term * wB_flux)
                if (hBL_neut_h_molec > I_n_star**2) then
                  Gam_turb = I_VK * ((ln_neut - 2.0*log(I_n_star)) + (I_2Zeta_N*I_n_star - 1.0))
                else ! The layer dominated by molecular viscosity is smaller than the boundary layer.
                  Gam_turb = I_VK * (I_2Zeta_N*I_n_star - 1.0)
                endif
                I_Gam_T = 1.0 / (Gam_mol_t + Gam_turb)
                I_Gam_S = 1.0 / (Gam_mol_s + Gam_turb)
              endif
              wT_flux = dT_ustar * I_Gam_T
            else  ! gamma_T and gamma_S are a function of the buoyancy flux with proper iteration.
              ! Find the root where wB_flux is consistent with the values of gamma with that flux.

              ! First, determine the buoyancy flux assuming no effects of stability
              ! on the turbulence.  Following H & J '99, this limit also applies
              ! when the buoyancy flux is destabilizing.
              Gam_turb = I_VK * (ln_neut + (I_2Zeta_N - 1.0))
              I_Gam_T = 1.0 / (Gam_mol_t + Gam_turb)
              I_Gam_S = 1.0 / (Gam_mol_s + Gam_turb)
              wB_flux = (dB_dS * dS_ustar) * I_Gam_S + (dB_dT * dT_ustar) * I_Gam_T

              if (wB_flux < 0.0) then
                ! The buoyancy flux is stabilizing and will reduce the turbulent
                ! fluxes, and iteration is required.

                ! n_star <= 1.0 is the ratio of working boundary layer thickness
                ! to the neutral thickness.  I_n_star is its inverse.
                I_n_star = sqrt(1.0 - n_star_term * wB_flux)
                if (hBL_neut_h_molec > I_n_star**2) then
                  Gam_turb = I_VK * ((ln_neut - 2.0*log(I_n_star)) + (I_2Zeta_N*I_n_star - 1.0))
                else !   The layer dominated by molecular viscosity is smaller than the boundary layer.
                  Gam_turb = I_VK * (I_2Zeta_N*I_n_star - 1.0)
                endif
                I_Gam_T = 1.0 / (Gam_mol_t + Gam_turb)
                I_Gam_S = 1.0 / (Gam_mol_s + Gam_turb)

                wB_flux_new = (dB_dS * dS_ustar) * I_Gam_S + (dB_dT * dT_ustar) * I_Gam_T
                root_found = (abs(wB_flux_new - wB_flux) < CS%buoy_flux_tol*(abs(wB_flux_new) + abs(wB_flux)))
                ! Do not update the flux if its maagnitude would be increased by the otherwise
                ! stabilizing buoyancy fluxes.  This can happen when the buoyancy flux
                ! is stabilizing when one of the heat or salt fluxes are destabilizing due
                ! to their different molecular properties.
                if (wB_flux_new <= wB_flux) root_found = .true.

                if (.not.root_found) then
                  wB_flux_max = 0.0 ; DwB_max = wB_flux
                  wB_flux_min = wB_flux ; DwB_min = wB_flux_new - wB_flux

                  if ((wB_flux_min*n_star_term < (1.0 - hBL_neut_h_molec)) .and. &
                      ((1.0 - hBL_neut_h_molec) < wB_flux_max*n_star_term)) then
                    ! The derivative of Gam_turb with wB_flux has a discontinuous change within the
                    ! bracketed range of values.  Take this discontinous slope value for a first
                    ! guess, because Newton's method and the false position method may not converge
                    ! quickly when this discontinuity is between a guess and the solution.
                    wB_flux = (1.0 - hBL_neut_h_molec) / n_star_term
                    I_n_star = sqrt(hBL_neut_h_molec)
                    Gam_turb = I_VK * (I_2Zeta_N*I_n_star - 1.0)
                    I_Gam_T = 1.0 / (Gam_mol_t + Gam_turb)
                    I_Gam_S = 1.0 / (Gam_mol_s + Gam_turb)
                    wB_flux_new = (dB_dS * dS_ustar) * I_Gam_S + (dB_dT * dT_ustar) * I_Gam_T

                    if (abs(wB_flux_new - wB_flux) <= CS%buoy_flux_tol*(abs(wB_flux_new) + abs(wB_flux))) then
                      ! The root has been found to within the tolerance at the kink.  This should be very rare.
                      root_found = .true.
                    elseif (wB_flux_new > wB_flux) then
                      ! The solution is in the limit where abs(wB_flux) is small and
                      ! Gam_turb = I_VK * ((ln_neut - 2.0*log(I_n_star)) + (I_2Zeta_N*I_n_star - 1.0))
                      wB_flux_min = wB_flux ; DwB_min = wB_flux_new - wB_flux
                    else
                      ! The solution is in the limt where abs(wB_flux) is large and
                      ! Gam_turb = I_VK * (I_2Zeta_N*I_n_star - 1.0)
                      wB_flux_max = wB_flux ; DwB_max = wB_flux_new - wB_flux
                    endif
                  endif
                endif

                if (.not.root_found) then
                  ! Use the false position for the next guess.
                  wB_flux = wB_flux_min + (wB_flux_max-wB_flux_min) * (DwB_min / (DwB_min - DwB_max))

                  do it3 = 1,30
                  ! Iterate using Newton's method with bounds or the false position method to find the root.

                    I_n_star = sqrt(1.0 - n_star_term * wB_flux)
                    dIns_dwB = -0.5 * n_star_term / I_n_star
                    if (hBL_neut_h_molec > I_n_star**2) then
                      Gam_turb = I_VK * ((ln_neut - 2.0*log(I_n_star)) + (I_2Zeta_N*I_n_star - 1.0))
                      dG_dwB =  I_VK * (( -2.0 / I_n_star + I_2Zeta_N) * dIns_dwB)
                    else
                      !   The layer dominated by molecular viscosity is smaller than the boundary layer.
                      Gam_turb = I_VK * (I_2Zeta_N*I_n_star - 1.0)
                      dG_dwB = I_VK * (I_2Zeta_N * dIns_dwB)
                    endif
                    I_Gam_T = 1.0 / (Gam_mol_t + Gam_turb)
                    I_Gam_S = 1.0 / (Gam_mol_s + Gam_turb)
                    wB_flux_new = (dB_dS * dS_ustar) * I_Gam_S + (dB_dT * dT_ustar) * I_Gam_T

                    ! Test for convergence to within tolerance at the point where wB_flux_new = wB_flux.
                    if (abs(wB_flux_new - wB_flux) <= CS%buoy_flux_tol*(abs(wB_flux_new) + abs(wB_flux))) &
                      root_found = .true.
                    if (root_found) exit

                    dDwB_dwB = -dG_dwB * ((dB_dS * dS_ustar) * I_Gam_S**2 + &
                                          (dB_dT * dT_ustar) * I_Gam_T**2) - 1.0
                    if ((dDwB_dwB >= 0.0) .or. &
                        ( wB_flux - wB_flux_new >= abs(dDwB_dwB)*(wB_flux_max - wB_flux)) .or. &
                        ( wB_flux - wB_flux_new <= abs(dDwB_dwB)*(wB_flux_min - wB_flux)) ) then
                      ! Use the False position method to determine the guess for the next iteration when
                      ! Newton's method would go out of bounds
                      wB_flux_next = wB_flux_min + (wB_flux_max-wB_flux_min) * (DwB_min / (DwB_min - DwB_max))
                    else
                      ! Use Newton's method for the next guess.
                      wB_flux_next = wB_flux - (wB_flux_new - wB_flux) / dDwB_dwB
                    endif

                    ! Reset one of the bounds inward.
                    if (wB_flux_new - wB_flux > 0) then
                      wB_flux_min = wB_flux ; DwB_min = wB_flux_new - wB_flux
                    else
                      wB_flux_max = wB_flux ; DwB_max = wB_flux_new - wB_flux
                    endif

                    ! Update wB_flux
                    wB_flux = wB_flux_next
                  enddo ! it3
                endif

              endif  ! End of test for first guess of wB_flux < 0.
              wT_flux = dT_ustar * I_Gam_T
            endif  ! End of test for CS%const_gamma

            ISS%tflux_ocn(i,j)  = RhoCp * wT_flux
            exch_vel_t(i,j) = ustar_h * I_Gam_T
            exch_vel_s(i,j) = ustar_h * I_Gam_S

            ! Calculate the heat flux inside the ice shelf.
            ! Vertical adv/diff as in H+J 1999, equations (26) & approx from (31).
            !   Q_ice = density_ice * CS%Cp_ice * K_ice * dT/dz (at interface)
            ! vertical adv/diff as in H+J 1999, equations (31) & (26)...
            !   dT/dz ~= min( (lprec/(density_ice*K_ice))*(CS%Temp_Ice-T_freeze) , 0.0 )
            ! If this approximation is not made, iterations are required... See H+J Fig 3.

            if (ISS%tflux_ocn(i,j) >= 0.0) then
              ! Freezing occurs due to downward ocean heat flux, so zero iout ce heat flux.
              ISS%water_flux(i,j) = -I_LF * ISS%tflux_ocn(i,j)
              ISS%tflux_shelf(i,j) = 0.0
            else
              if (CS%insulator) then
                !no conduction/perfect insulator
                ISS%tflux_shelf(i,j) = 0.0
                ISS%water_flux(i,j) = I_LF * (ISS%tflux_shelf(i,j) - ISS%tflux_ocn(i,j))

              else
                ! With melting, from H&J 1999, eqs (31) & (26)...
                !   Q_ice ~= Cp_ice * (CS%Temp_Ice-T_freeze) * lprec
                !   RhoLF*lprec = Q_ice - ISS%tflux_ocn(i,j)
                !   lprec = -(ISS%tflux_ocn(i,j)) / (CS%Lat_fusion + Cp_ice * (T_freeze-CS%Temp_Ice))
                ISS%water_flux(i,j) = -ISS%tflux_ocn(i,j) / &
                     (CS%Lat_fusion + CS%Cp_ice * (ISS%tfreeze(i,j) - CS%Temp_Ice))

                ISS%tflux_shelf(i,j) = ISS%tflux_ocn(i,j) + CS%Lat_fusion*ISS%water_flux(i,j)
              endif

            endif
            !other options: dTi/dz linear through shelf, with draft in [Z ~> m], KTI in [Z2 T-1 ~> m2 s-1]
            !    dTi_dz = (CS%Temp_Ice - ISS%tfreeze(i,j)) / draft(i,j)
            !    ISS%tflux_shelf(i,j) = Rho_Ice * CS%Cp_ice * KTI * dTi_dz


            if (CS%find_salt_root) then
              exit ! no need to do interaction, so exit loop
            else

              mass_exch = exch_vel_s(i,j) * CS%Rho_ocn
              Sbdry_it = (sfc_state%sss(i,j) * mass_exch + CS%Salin_ice * ISS%water_flux(i,j)) / &
                         (mass_exch + ISS%water_flux(i,j))
              dS_it = Sbdry_it - Sbdry(i,j)
              if (abs(dS_it) < 1.0e-4*(0.5*(sfc_state%sss(i,j) + Sbdry(i,j) + 1.0e-10*US%ppt_to_S))) exit

              if (dS_it < 0.0) then ! Sbdry is now the upper bound.
                if (Sb_max_set) then
                  if (Sbdry(i,j) > Sb_max) &
                    call MOM_error(FATAL,"shelf_calc_flux: Irregular iteration for Sbdry (max).")
                endif
                Sb_max = Sbdry(i,j) ; dS_max = dS_it ; Sb_max_set = .true.
              else ! Sbdry is now the lower bound.
                if (Sb_min_set) then
                  if (Sbdry(i,j) < Sb_min) &
                    call MOM_error(FATAL, "shelf_calc_flux: Irregular iteration for Sbdry (min).")
                endif
                Sb_min = Sbdry(i,j) ; dS_min = dS_it ; Sb_min_set = .true.
              endif ! dS_it < 0.0

              if (Sb_min_set .and. Sb_max_set) then
                ! Use the false position method for the next iteration.
                Sbdry(i,j) = Sb_min + (Sb_max-Sb_min) * (dS_min / (dS_min - dS_max))
              else
                Sbdry(i,j) = Sbdry_it
              endif ! Sb_min_set

              if (.not.CS%salt_flux_itt_bugfix) Sbdry(i,j) = Sbdry_it

            endif ! CS%find_salt_root

          enddo !it1
          ! Check for non-convergence and/or non-boundedness?

        else
          !   In the 2-equation form, the mixed layer turbulent exchange velocity
          ! is specified and large enough that the ocean salinity at the interface
          ! is about the same as the boundary layer salinity.
          ! The following two lines are equivalent:
          ! call calculate_TFreeze(Sbdry(i,j), p_int(i), ISS%tfreeze(i,j), CS%eqn_of_state, scale_from_EOS=.true.)
          call calculate_TFreeze(sfc_state%SSS(i:i,j), p_int(i:i), ISS%tfreeze(i:i,j), CS%eqn_of_state)

          exch_vel_t(i,j) = CS%gamma_t
          ISS%tflux_ocn(i,j) = RhoCp * exch_vel_t(i,j) * (ISS%tfreeze(i,j) - sfc_state%sst(i,j))
          ISS%tflux_shelf(i,j) = 0.0
          ISS%water_flux(i,j) = -I_LF * ISS%tflux_ocn(i,j)
          Sbdry(i,j) = 0.0
        endif
      elseif (ISS%area_shelf_h(i,j) > 0.0) then ! This is an ice-sheet, not a floating shelf.
        ISS%tflux_ocn(i,j) = 0.0
      else ! There is no ice shelf or sheet here.
        ISS%tflux_ocn(i,j) = 0.0
      endif

!      haline_driving(i,j) = sfc_state%sss(i,j) - Sbdry(i,j)

    enddo ! i-loop
  enddo ! j-loop

  if (allocated(sfc_state%frazil)) then
    add_frazil = .true.
  else
    add_frazil = .false.
  endif

  do j=js,je ; do i=is,ie
    ! ISS%water_flux = net liquid water into the ocean [R Z T-1 ~> kg m-2 s-1]
    if (CS%flux_factor/=1.0) then
      ISS%water_flux(i,j) = ISS%water_flux(i,j) * CS%flux_factor
      ISS%tflux_ocn(i,j) = ISS%tflux_ocn(i,j) * CS%flux_factor
      if (CS%threeeq .and. ISS%tflux_ocn(i,j) < 0.0 .and. (.not. CS%insulator)) &
        ISS%tflux_shelf(i,j)=ISS%tflux_ocn(i,j) + CS%Lat_fusion * ISS%water_flux(i,j)
    endif

    if ((sfc_state%ocean_mass(i,j) > CS%col_mass_melt_threshold) .and. &
        (ISS%area_shelf_h(i,j) > 0.0) .and.  (CS%isthermo)) then

      ! Set melt to zero above a cutoff pressure (CS%Rho_ocn*CS%cutoff_depth*CS%g_Earth).
      ! This is needed for the ISOMIP test case.
      if (ISS%mass_shelf(i,j) < CS%Rho_ocn*CS%cutoff_depth) then
        ISS%water_flux(i,j) = 0.0
      endif
      ! Compute haline driving, which is one of the diags. used in ISOMIP
      if (exch_vel_s(i,j)>0.) haline_driving(i,j) = (ISS%water_flux(i,j) * Sbdry(i,j)) / (CS%Rho_ocn * exch_vel_s(i,j))

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!Safety checks !!!!!!!!!!!!!!!!!!!!!!!!!
      !1)Check if haline_driving computed above is consistent with
      ! haline_driving = sfc_state%sss - Sbdry
      !if (ISS%water_flux(i,j) /= 0.0) then
      !   if (haline_driving(i,j) /= (sfc_state%sss(i,j) - Sbdry(i,j))) then
      !     write(mesg,*) 'at i,j=',i,j,' haline_driving, sss-Sbdry',US%S_to_ppt*haline_driving(i,j), &
      !                   US%S_to_ppt*(sfc_state%sss(i,j) - Sbdry(i,j))
      !     call MOM_error(FATAL, &
      !            "shelf_calc_flux: Inconsistency in melt and haline_driving"//trim(mesg))
      !   endif
      !endif

      ! 2) check if |melt| > 0 when ustar_shelf = 0.
      ! this should never happen
      if ((abs(ISS%water_flux(i,j))>0.0) .and. (fluxes%ustar_shelf(i,j) == 0.0)) then
        write(mesg,*) "|melt| = ",ISS%water_flux(i,j)," > 0 and ustar_shelf = 0. at i,j", i, j
        call MOM_error(FATAL, "shelf_calc_flux: "//trim(mesg))
      endif
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!End of safety checks !!!!!!!!!!!!!!!!!!!
    elseif (ISS%area_shelf_h(i,j) > 0.0) then
      ! This is grounded ice, that could be modified to melt if a geothermal heat flux were used.
      haline_driving(i,j) = 0.0
      ISS%water_flux(i,j) = 0.0
    endif ! area_shelf_h

    ! mass flux [R Z L2 T-1 ~> kg s-1], part of ISOMIP diags.
    mass_flux(i,j) = ISS%water_flux(i,j) * ISS%area_shelf_h(i,j)

    !Add frazil formation
    if (add_frazil .and. (ISS%hmask(i,j) == 1 .or. ISS%hmask(i,j) == 2)) &
      ISS%water_flux(i,j) = ISS%water_flux(i,j) - ISS%frazil(i,j) * I_dt_LHF
    fluxes%iceshelf_melt(i,j) = ISS%water_flux(i,j)
  enddo ; enddo ! i- and j-loops

  if (CS%active_shelf_dynamics .or. CS%override_shelf_movement) then
    call cpu_clock_begin(id_clock_pass)
    call pass_var(ISS%area_shelf_h, G%domain, complete=.false.)
    call pass_var(ISS%mass_shelf, G%domain)
    call cpu_clock_end(id_clock_pass)
  endif

  ! Melting has been computed, now is time to update thickness and mass
  if ( CS%override_shelf_movement .and. (.not.CS%mass_from_file)) then
    if (CS%bmb_diag) dh_bdott(is:ie,js:je) = ISS%h_shelf(is:ie,js:je)
    call change_thickness_using_melt(ISS, G, US, time_step, fluxes, CS%density_ice, CS%debug)
    if (CS%bmb_diag) dh_bdott(is:ie,js:je) = ISS%h_shelf(is:ie,js:je) - dh_bdott(is:ie,js:je)

    if (CS%debug) then
      call hchksum(ISS%h_shelf, "h_shelf after change thickness using melt", G%HI, haloshift=0, unscale=US%Z_to_m)
      call hchksum(ISS%mass_shelf, "mass_shelf after change thickness using melt", G%HI, haloshift=0, &
                   unscale=US%RZ_to_kg_m2)
    endif
  endif

  ! Melting has been computed, now is time to update thickness and mass with dynamic ice shelf
  if (CS%active_shelf_dynamics) then

    ISS%dhdt_shelf(:,:) = ISS%h_shelf(:,:)

    if (CS%bmb_diag) dh_bdott(is:ie,js:je) = ISS%h_shelf(is:ie,js:je)
    call change_thickness_using_melt(ISS, G, US, time_step, fluxes, CS%density_ice, CS%debug)
    if (CS%bmb_diag) dh_bdott(is:ie,js:je) = ISS%h_shelf(is:ie,js:je) - dh_bdott(is:ie,js:je)

    if (CS%debug) then
      call hchksum(ISS%h_shelf, "h_shelf after change thickness using melt", G%HI, haloshift=0, unscale=US%Z_to_m)
      call hchksum(ISS%mass_shelf, "mass_shelf after change thickness using melt", G%HI, haloshift=0, &
                   unscale=US%RZ_to_kg_m2)
    endif

    if (CS%smb_diag) dh_adott(is:ie,js:je) = ISS%h_shelf(is:ie,js:je)
    call change_thickness_using_precip(CS, ISS, G, US, fluxes, time_step, Time)
    if (CS%smb_diag) dh_adott(is:ie,js:je) = ISS%h_shelf(is:ie,js:je) - dh_adott(is:ie,js:je)

    if (CS%debug) then
      call hchksum(ISS%h_shelf, "h_shelf after change thickness using surf acc", G%HI, haloshift=0, unscale=US%Z_to_m)
      call hchksum(ISS%mass_shelf, "mass_shelf after change thickness using surf acc", G%HI, haloshift=0, &
                   unscale=US%RZ_to_kg_m2)
    endif

    update_ice_vel = .false.
    coupled_GL = (CS%GL_couple .and. .not. CS%solo_ice_sheet)

    ! advect the ice shelf, and advance the front. Calving will be in here somewhere as well..
    ! when we decide on how to do it
    call update_ice_shelf(CS%dCS, ISS, G, US, time_step, Time, CS%calve_ice_shelf_bergs, &
                          sfc_state%ocean_mass, coupled_GL)

    do j=js,je ; do i=is,ie
      ISS%dhdt_shelf(i,j) = (ISS%h_shelf(i,j) - ISS%dhdt_shelf(i,j))*Itime_step
    enddo ; enddo

    call IS_dynamics_post_data(time_step, Time, CS%dCS, ISS, G)
  endif

  if (CS%shelf_mass_is_dynamic) &
    call write_ice_shelf_energy(CS%dCS, G, US, ISS%mass_shelf, ISS%area_shelf_h, Time, &
                                time_step=real_to_time(time_step, unscale=US%T_to_s) )

  if (CS%debug) call MOM_forcing_chksum("Before add shelf flux", fluxes, G, CS%US, haloshift=0)

  ! pass on the updated ice sheet geometry (for pressure on ocean) and thermodynamic data
  call add_shelf_flux(G, US, CS, sfc_state, fluxes, time_step)

  call enable_averages(time_step, Time, CS%diag)
  if (CS%id_shelf_mass > 0) call post_data(CS%id_shelf_mass, ISS%mass_shelf, CS%diag)
  if (CS%id_area_shelf_h > 0) call post_data(CS%id_area_shelf_h, ISS%area_shelf_h, CS%diag)
  if (CS%id_ustar_shelf > 0) call post_data(CS%id_ustar_shelf, fluxes%ustar_shelf, CS%diag)
  if (CS%id_shelf_sfc_mass_flux > 0) call post_data(CS%id_shelf_sfc_mass_flux, fluxes%shelf_sfc_mass_flux, CS%diag)

  if (CS%id_melt > 0) call post_data(CS%id_melt, fluxes%iceshelf_melt, CS%diag)
  if (CS%id_thermal_driving > 0) call post_data(CS%id_thermal_driving, (sfc_state%sst-ISS%tfreeze), CS%diag)
  if (CS%id_Sbdry > 0) call post_data(CS%id_Sbdry, Sbdry, CS%diag)
  if (CS%id_haline_driving > 0) call post_data(CS%id_haline_driving, haline_driving, CS%diag)
  if (CS%id_mass_flux > 0) call post_data(CS%id_mass_flux, mass_flux, CS%diag)
  if (CS%id_u_ml > 0) call post_data(CS%id_u_ml, sfc_state%u, CS%diag)
  if (CS%id_v_ml > 0) call post_data(CS%id_v_ml, sfc_state%v, CS%diag)
  if (CS%id_tfreeze > 0) call post_data(CS%id_tfreeze, ISS%tfreeze, CS%diag)
  if (CS%id_tfl_shelf > 0) call post_data(CS%id_tfl_shelf, ISS%tflux_shelf, CS%diag)
  if (CS%id_exch_vel_t > 0) call post_data(CS%id_exch_vel_t, exch_vel_t, CS%diag)
  if (CS%id_exch_vel_s > 0) call post_data(CS%id_exch_vel_s, exch_vel_s, CS%diag)
  if (CS%id_h_shelf > 0) call post_data(CS%id_h_shelf, ISS%h_shelf, CS%diag)
  if (CS%id_dhdt_shelf > 0) call post_data(CS%id_dhdt_shelf, ISS%dhdt_shelf, CS%diag)
  if (CS%id_h_mask > 0) call post_data(CS%id_h_mask,ISS%hmask,CS%diag)
  if (CS%id_frazil > 0) call post_data(CS%id_frazil,ISS%frazil,CS%diag)
  if (CS%active_shelf_dynamics) &
      call process_and_post_scalar_data(CS, vaf0, vaf0_A, vaf0_G, Itime_step, dh_adott, dh_bdott)
  call disable_averaging(CS%diag)

  !reset used frazil
  if (add_frazil) ISS%frazil(:,:) = 0.0

  call cpu_clock_end(id_clock_shelf)

  if (CS%debug) call MOM_forcing_chksum("End of shelf calc flux", fluxes, G, CS%US, haloshift=0)

  if (CS%rotate_index) then
!   call rotate_surface_state(sfc_state, sfc_state_in, CS%Grid_in, -CS%turns)
    call rotate_forcing(fluxes, fluxes_in, -CS%turns)
    call deallocate_surface_state(sfc_state)
    deallocate(sfc_state)
    call deallocate_forcing_type(fluxes)
    deallocate(fluxes)
  endif

end procedure shelf_calc_flux
module procedure adjust_ice_sheet_frazil
  type(ocean_grid_type), pointer :: G => NULL()  !< The grid structure used by the ice shelf.
  type(ice_shelf_state), pointer :: ISS => NULL() !< A structure with elements that describe
  type(surface), pointer :: sfc_state => NULL()
  type(forcing), pointer :: fluxes => NULL()
  integer :: i,j,is,ie,js,je
  G => CS%grid ; ISS => CS%ISS

  if (CS%rotate_index) then
    allocate(sfc_state)
    call rotate_surface_state(sfc_state_in, sfc_state, G, CS%turns)
    allocate(fluxes)
    call allocate_forcing_type(fluxes_in, G, fluxes, turns=CS%turns)
    call rotate_forcing(fluxes_in, fluxes, CS%turns)
  else
    sfc_state => sfc_state_in
    fluxes => fluxes_in
  endif

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec

  do j=js,je ; do i=is,ie
    !Copy frazil to the ice sheet module where ice sheet is present.
    !No scaling to account for partial ice-sheet cells is necessary here, as
    !this is taken care of when applied to the ice sheet.
    if (fluxes%frac_shelf_h(i,j)>0.0) ISS%frazil(i,j) = sfc_state%frazil(i,j)
    !Remove the frazil that is used by the ice sheet from sfc_state%frazil
    !The sfc_state%frazil is sent to the sea-ice module
    sfc_state%frazil(i,j) = sfc_state%frazil(i,j) * (1.0-fluxes%frac_shelf_h(i,j))
  enddo ; enddo

  if (CS%rotate_index) then
    call rotate_surface_state(sfc_state, sfc_state_in, G, -CS%turns)
    ! call rotate_forcing(fluxes, fluxes_in, -CS%turns)
    call deallocate_surface_state(sfc_state)
    deallocate(sfc_state)
    call deallocate_forcing_type(fluxes)
    deallocate(fluxes)
  endif
end procedure adjust_ice_sheet_frazil
module procedure integrate_over_ice_sheet_area
  integer :: IS_ID ! local copy of hemisphere
  real, dimension(SZI_(G),SZJ_(G))  :: var_cell !< Variable integrated over the ice-sheet area of each cell
  integer, dimension(SZI_(G),SZJ_(G))  :: mask ! a mask for active cells depending on hemisphere indicated
  integer :: i, j
  if (present(hemisphere)) then
    IS_ID = hemisphere
  else
    IS_ID = -1
  endif

  mask(:,:) = 0
  if (IS_ID==0) then     !Antarctica (S. Hemisphere) only
    do j = G%jsc,G%jec ; do i = G%isc,G%iec
      if (ISS%hmask(i,j)>0 .and. G%geoLatT(i,j)<=0.0) mask(i,j)=1
    enddo ; enddo
  elseif (IS_ID==1) then !Greenland (N. Hemisphere) only
    do j = G%jsc,G%jec ; do i = G%isc,G%iec
      if (ISS%hmask(i,j)>0 .and. G%geoLatT(i,j)>0.0)  mask(i,j)=1
    enddo ; enddo
  else                   !All ice sheets
    mask(G%isc:G%iec,G%jsc:G%jec) = ISS%hmask(G%isc:G%iec,G%jsc:G%jec)
  endif

  var_cell(:,:) = 0.0
  do j = G%jsc,G%jec ; do i = G%isc,G%iec
    if (mask(i,j)>0) var_cell(i,j) = var(i,j) * ISS%area_shelf_h(i,j)
  enddo ; enddo

  var_out = reproducing_sum(var_cell, unscale=unscale*G%US%L_to_m**2)
end procedure integrate_over_ice_sheet_area
module procedure ice_sheet_calving_to_ocean_sfc
  type(ice_shelf_state), pointer :: ISS => NULL() !< A structure with elements that describe
  type(ocean_grid_type), pointer :: G => NULL()   !< A pointer to the ocean grid metric.
  integer :: is, ie, js, je
  G=>CS%Grid
  ISS => CS%ISS
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec

  calving = US%RZ_T_to_kg_m2s * ISS%calving(is:ie,js:je)
  calving_hflx = US%QRZ_T_to_W_m2 * ISS%calving_hflx(is:ie,js:je)

  !CS%calve_ice_shelf_bergs=.true.

end procedure ice_sheet_calving_to_ocean_sfc
module procedure change_thickness_using_melt
  real :: I_rho_ice ! Ice specific volume [R-1 ~> m3 kg-1]
  integer :: i, j
  I_rho_ice = 1.0 / density_ice


  do j=G%jsc,G%jec ; do i=G%isc,G%iec
    if ((ISS%hmask(i,j) == 1) .or. (ISS%hmask(i,j) == 2)) then
      ! first, zero out fluxes applied during previous time step
      if (associated(fluxes%lprec)) fluxes%lprec(i,j) = 0.0
      if (associated(fluxes%sens)) fluxes%sens(i,j) = 0.0
      if (associated(fluxes%frac_shelf_h)) fluxes%frac_shelf_h(i,j) = 0.0
      if (associated(fluxes%salt_flux)) fluxes%salt_flux(i,j) = 0.0

      if (ISS%water_flux(i,j) * time_step / density_ice < ISS%h_shelf(i,j)) then
        ISS%h_shelf(i,j) = ISS%h_shelf(i,j) - ISS%water_flux(i,j) * time_step / density_ice
      else
        ! the ice is about to melt away, so set thickness, area, and mask to zero
        ! NOTE: this is not mass conservative should maybe scale salt & heat flux for this cell
        ISS%h_shelf(i,j) = 0.0
        ISS%hmask(i,j) = 0.0
        ISS%area_shelf_h(i,j) = 0.0
      endif
      ISS%mass_shelf(i,j) = ISS%h_shelf(i,j) * density_ice
    endif
  enddo ; enddo

  call pass_var(ISS%area_shelf_h, G%domain, complete=.false.)
  call pass_var(ISS%h_shelf, G%domain, complete=.false.)
  call pass_var(ISS%hmask, G%domain, complete=.false.)
  call pass_var(ISS%mass_shelf, G%domain)

end procedure change_thickness_using_melt
module procedure add_shelf_forces
  type(ocean_grid_type), pointer :: G => NULL()   !< A pointer to the ocean grid metric.
  type(mech_forcing),    pointer :: forces     !< A structure with the driving mechanical forces
  real :: kv_rho_ice ! The viscosity of ice divided by its density [L4 T-1 R-1 Z-2 ~> m5 kg-1 s-1].
  real :: press_ice  ! The pressure of the ice shelf per unit area of ocean (not ice) [R L2 T-2 ~> Pa].
  logical :: find_area ! If true find the shelf areas at u & v points.
  logical :: rotate = .false.
  type(ice_shelf_state), pointer :: ISS => NULL() ! A structure with elements that describe
  integer :: i, j, is, ie, js, je, isd, ied, jsd, jed
  rotate = .false. ; if (present(external_call)) rotate = external_call

  if (CS%rotate_index .and. rotate) then
    if ((Ocn_grid%isc /= CS%Grid_in%isc) .or. (Ocn_grid%iec /= CS%Grid_in%iec) .or. &
        (Ocn_grid%jsc /= CS%Grid_in%jsc) .or. (Ocn_grid%jec /= CS%Grid_in%jec)) &
      call MOM_error(FATAL,"add_shelf_forces: Incompatible Ocean and external Ice shelf grids.")
    allocate(forces)
    call allocate_mech_forcing(forces_in, CS%Grid, forces)
    call rotate_mech_forcing(forces_in, CS%turns, forces)
  else
    if ((Ocn_grid%isc /= CS%Grid%isc) .or. (Ocn_grid%iec /= CS%Grid%iec) .or. &
        (Ocn_grid%jsc /= CS%Grid%jsc) .or. (Ocn_grid%jec /= CS%Grid%jec)) &
      call MOM_error(FATAL,"add_shelf_forces: Incompatible Ocean and internal Ice shelf grids.")

    forces=>forces_in
  endif

  G=>CS%Grid

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec
  isd = G%isd ; jsd = G%jsd ; ied = G%ied ; jed = G%jed

  ISS => CS%ISS

  find_area = .true. ; if (present(do_shelf_area)) find_area = do_shelf_area

  if (find_area) then
    ! The frac_shelf is set over the widest possible area. Could it be smaller?
    do j=jsd,jed ; do I=isd,ied-1
      forces%frac_shelf_u(I,j) = 0.0
      if ((G%areaT(i,j) + G%areaT(i+1,j) > 0.0)) & ! .and. (G%areaCu(I,j) > 0.0)) &
        forces%frac_shelf_u(I,j) = (ISS%area_shelf_h(i,j) + ISS%area_shelf_h(i+1,j)) / &
                                   (G%areaT(i,j) + G%areaT(i+1,j))
    enddo ; enddo
    do J=jsd,jed-1 ; do i=isd,ied
      forces%frac_shelf_v(i,J) = 0.0
      if ((G%areaT(i,j) + G%areaT(i,j+1) > 0.0)) & ! .and. (G%areaCv(i,J) > 0.0)) &
        forces%frac_shelf_v(i,J) = (ISS%area_shelf_h(i,j) + ISS%area_shelf_h(i,j+1)) / &
                                   (G%areaT(i,j) + G%areaT(i,j+1))
    enddo ; enddo
    call pass_vector(forces%frac_shelf_u, forces%frac_shelf_v, G%domain, TO_ALL, CGRID_NE)
  endif

  do j=js,je ; do i=is,ie
    press_ice = (ISS%area_shelf_h(i,j) * G%IareaT(i,j)) * (CS%g_Earth * ISS%mass_shelf(i,j))
    if (associated(forces%p_surf)) then
      if (.not.forces%accumulate_p_surf) forces%p_surf(i,j) = 0.0
      forces%p_surf(i,j) = forces%p_surf(i,j) + press_ice
    endif
    if (associated(forces%p_surf_full)) then
      if (.not.forces%accumulate_p_surf) forces%p_surf_full(i,j) = 0.0
      forces%p_surf_full(i,j) = forces%p_surf_full(i,j) + press_ice
    endif
  enddo ; enddo

  ! For various reasons, forces%rigidity_ice_[uv] is always updated here. Note
  ! that it may have been zeroed out where IOB is translated to forces and
  ! contributions from icebergs and the sea-ice pack added subsequently.
  !### THE RIGIDITY SHOULD ALSO INCORPORATE AREAL-COVERAGE INFORMATION.
  kv_rho_ice = CS%kv_ice / CS%density_ice
  do j=js,je ; do I=is-1,ie
    if (.not.forces%accumulate_rigidity) forces%rigidity_ice_u(I,j) = 0.0
    forces%rigidity_ice_u(I,j) = forces%rigidity_ice_u(I,j) + &
            kv_rho_ice * min(ISS%mass_shelf(i,j), ISS%mass_shelf(i+1,j))
  enddo ; enddo
  do J=js-1,je ; do i=is,ie
    if (.not.forces%accumulate_rigidity) forces%rigidity_ice_v(i,J) = 0.0
    forces%rigidity_ice_v(i,J) = forces%rigidity_ice_v(i,J) + &
            kv_rho_ice * min(ISS%mass_shelf(i,j), ISS%mass_shelf(i,j+1))
  enddo ; enddo

  if (CS%debug) then
    call uvchksum("rigidity_ice_[uv]", forces%rigidity_ice_u, &
        forces%rigidity_ice_v, CS%Grid%HI, symmetric=.true., &
        unscale=US%L_to_m**3*US%L_to_Z*US%s_to_T, scalar_pair=.true.)
    call uvchksum("frac_shelf_[uv]", forces%frac_shelf_u, &
        forces%frac_shelf_v, CS%Grid%HI, symmetric=.true., &
        scalar_pair=.true.)
  endif

  if (CS%rotate_index .and. rotate) then
    call rotate_mech_forcing(forces, -CS%turns, forces_in)
    call deallocate_mech_forcing(forces)
  endif

end procedure add_shelf_forces
module procedure add_shelf_pressure
  type(ocean_grid_type), pointer :: G => NULL()  ! A pointer to  ocean's grid structure.
  real :: press_ice       !< The pressure of the ice shelf per unit area of ocean (not ice) [R L2 T-2 ~> Pa].
  integer :: i, j, is, ie, js, je
  G=>CS%Grid
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec

  if ((CS%grid%isc /= G%isc) .or. (CS%grid%iec /= G%iec) .or. &
      (CS%grid%jsc /= G%jsc) .or. (CS%grid%jec /= G%jec)) &
    call MOM_error(FATAL,"add_shelf_pressure: Incompatible ocean and ice shelf grids.")

  do j=js,je ; do i=is,ie
    press_ice = (CS%ISS%area_shelf_h(i,j) * G%IareaT(i,j)) * (CS%g_Earth * CS%ISS%mass_shelf(i,j))
    if (associated(fluxes%p_surf)) then
      if (.not.fluxes%accumulate_p_surf) fluxes%p_surf(i,j) = 0.0
      fluxes%p_surf(i,j) = fluxes%p_surf(i,j) + press_ice
    endif
    if (associated(fluxes%p_surf_full)) then
      if (.not.fluxes%accumulate_p_surf) fluxes%p_surf_full(i,j) = 0.0
      fluxes%p_surf_full(i,j) = fluxes%p_surf_full(i,j) + press_ice
    endif
  enddo ; enddo

end procedure add_shelf_pressure
module procedure add_shelf_flux
  real :: frac_shelf       !< The fractional area covered by the ice shelf [nondim].
  real :: frac_open        !< The fractional area of the ocean that is not covered by the ice shelf [nondim].
  real :: delta_mass_shelf !< Change in ice shelf mass over one time step [R Z L2 T-1 ~> kg s-1]
  real :: balancing_flux   !< The fresh water flux that balances the integrated melt flux [R Z T-1 ~> kg m-2 s-1]
  real :: balancing_area   !< total area where the balancing flux is applied [L2 ~> m2]
  type(time_type) :: dTime !< The time step as a time_type
  type(time_type) :: Time0 !< The previous time (Time-dt)
  real, dimension(SZDI_(G),SZDJ_(G)) :: bal_frac  !< Fraction of the cell where the mass flux
  real, dimension(SZDI_(G),SZDJ_(G)) :: last_mass_shelf !< Ice shelf mass
  real, dimension(SZDI_(G),SZDJ_(G)) :: delta_float_mass   !< The change in the floating mass between
  real, dimension(SZDI_(G),SZDJ_(G))  :: last_h_shelf !< Ice shelf thickness [Z ~> m]
  real, dimension(SZDI_(G),SZDJ_(G))  :: last_hmask !< Ice shelf mask [nondim]
  real, dimension(SZDI_(G),SZDJ_(G))  :: last_area_shelf_h !< Ice shelf area [L2 ~> m2]
  real, dimension(SZDI_(G),SZDJ_(G))  :: delta_draft !< change in ice shelf draft thickness [L ~> m]
  type(ice_shelf_state), pointer :: ISS => NULL() !< A structure with elements that describe
  character(len=160) :: mesg  ! The text of an error message
  integer :: i, j, is, ie, js, je, isd, ied, jsd, jed
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec
  isd = G%isd ; jsd = G%jsd ; ied = G%ied ; jed = G%jed

  if ((CS%grid%isc /= G%isc) .or. (CS%grid%iec /= G%iec) .or. &
      (CS%grid%jsc /= G%jsc) .or. (CS%grid%jec /= G%jec)) &
    call MOM_error(FATAL,"add_shelf_flux: Incompatible ocean and ice shelf grids.")

  ISS => CS%ISS


  call add_shelf_pressure(G, US, CS, fluxes)

  ! Determine ustar and the square magnitude of the velocity in the
  ! bottom boundary layer. Together these give the TKE source and
  ! vertical decay scale.

  if (CS%debug) then
    if (allocated(sfc_state%taux_shelf) .and. allocated(sfc_state%tauy_shelf)) then
      call uvchksum("tau[xy]_shelf", sfc_state%taux_shelf, sfc_state%tauy_shelf, &
                    G%HI, haloshift=0, unscale=US%RZ_T_to_kg_m2s*US%L_T_to_m_s)
    endif
  endif

  if (CS%active_shelf_dynamics .or. CS%override_shelf_movement) then
    do j=jsd,jed ; do i=isd,ied
      if (G%areaT(i,j) > 0.0) &
        fluxes%frac_shelf_h(i,j) = min(1.0, ISS%area_shelf_h(i,j) * G%IareaT(i,j))
    enddo ; enddo
  endif

  if (CS%debug) then
    call MOM_forcing_chksum("Before adding shelf fluxes", fluxes, G, CS%US, haloshift=0)
  endif

  do j=js,je ; do i=is,ie ; if (ISS%area_shelf_h(i,j) > 0.0) then
    ! Replace fluxes intercepted by the ice shelf with fluxes from the ice shelf
    frac_shelf = min(1.0, ISS%area_shelf_h(i,j) * G%IareaT(i,j))
    frac_open = max(0.0, 1.0 - frac_shelf)

    if (associated(fluxes%sw)) fluxes%sw(i,j) = frac_open * fluxes%sw(i,j)
    if (associated(fluxes%sw_vis_dir)) fluxes%sw_vis_dir(i,j) = frac_open * fluxes%sw_vis_dir(i,j)
    if (associated(fluxes%sw_vis_dif)) fluxes%sw_vis_dif(i,j) = frac_open * fluxes%sw_vis_dif(i,j)
    if (associated(fluxes%sw_nir_dir)) fluxes%sw_nir_dir(i,j) = frac_open * fluxes%sw_nir_dir(i,j)
    if (associated(fluxes%sw_nir_dif)) fluxes%sw_nir_dif(i,j) = frac_open * fluxes%sw_nir_dif(i,j)
    if (associated(fluxes%lw)) fluxes%lw(i,j) = frac_open * fluxes%lw(i,j)
    if (associated(fluxes%latent)) fluxes%latent(i,j) = frac_open * fluxes%latent(i,j)
    if (associated(fluxes%evap)) fluxes%evap(i,j) = frac_open * fluxes%evap(i,j)
    if (associated(fluxes%lprec)) then
      if (ISS%water_flux(i,j) > 0.0) then
        fluxes%lprec(i,j) =  frac_shelf*ISS%water_flux(i,j) + frac_open * fluxes%lprec(i,j)
      else
        fluxes%lprec(i,j) = frac_open * fluxes%lprec(i,j)
        fluxes%evap(i,j) = fluxes%evap(i,j) + frac_shelf*ISS%water_flux(i,j)
      endif
    endif

    if (associated(fluxes%sens)) &
      fluxes%sens(i,j) = frac_shelf*ISS%tflux_ocn(i,j) + frac_open * fluxes%sens(i,j)
    ! The salt flux should be mostly from sea ice, so perhaps none should be intercepted and this should be changed.
    if (associated(fluxes%salt_flux)) &
      fluxes%salt_flux(i,j) = frac_shelf * ISS%salt_flux(i,j)*CS%flux_factor + frac_open * fluxes%salt_flux(i,j)
  endif ; enddo ; enddo

  if (CS%debug) then
    call hchksum(ISS%water_flux, "water_flux add shelf fluxes", G%HI, haloshift=0, unscale=US%RZ_T_to_kg_m2s)
    call hchksum(ISS%tflux_ocn, "tflux_ocn add shelf fluxes", G%HI, haloshift=0, unscale=US%QRZ_T_to_W_m2)
    call MOM_forcing_chksum("After adding shelf fluxes", fluxes, G, CS%US, haloshift=0)
  endif

  ! Keep sea level constant by removing mass via a balancing flux that might be applied
  ! in the open ocean or the sponge region (via virtual precip, vprec). Apply additional
  ! salt/heat fluxes so that the resultant surface buoyancy forcing is ~ 0.
  ! This is needed for some of the ISOMIP+ experiments.

  if (CS%constant_sea_level) then
    if (.not. associated(fluxes%salt_flux)) allocate(fluxes%salt_flux(ie,je))
    if (.not. associated(fluxes%vprec)) allocate(fluxes%vprec(ie,je))
    fluxes%salt_flux(:,:) = 0.0 ; fluxes%vprec(:,:) = 0.0

    ! take into account changes in mass (or thickness) when imposing ice shelf mass
    if (CS%override_shelf_movement .and. CS%mass_from_file) then
      dTime = real_to_time(CS%time_step, unscale=US%T_to_s)

      ! Compute changes in mass after at least one full time step
      if (CS%Time > dTime) then
        Time0 = CS%Time - dTime
        do j=js,je ; do i=is,ie
          last_hmask(i,j) = ISS%hmask(i,j) ; last_area_shelf_h(i,j) = ISS%area_shelf_h(i,j)
        enddo ; enddo
        call time_interp_external(CS%mass_handle, Time0, last_mass_shelf, scale=US%kg_m3_to_R*US%m_to_Z)
        do j=js,je ; do i=is,ie
          last_h_shelf(i,j) = last_mass_shelf(i,j) / CS%density_ice
        enddo ; enddo

        ! apply calving
        if (CS%min_thickness_simple_calve > 0.0) then
          call ice_shelf_min_thickness_calve(G, last_h_shelf, last_area_shelf_h, last_hmask, &
                                       CS%min_thickness_simple_calve, halo=0)
          ! convert to mass again
          do j=js,je ; do i=is,ie
            last_mass_shelf(i,j) = last_h_shelf(i,j) * CS%density_ice
          enddo ; enddo
        endif

        ! get total ice shelf mass at (Time-dt) and (Time), in kg
        do j=js,je ; do i=is,ie
          ! Just consider the change in the mass of the floating shelf.
          if ((sfc_state%ocean_mass(i,j) > CS%min_ocean_mass_float) .and. &
              (ISS%area_shelf_h(i,j) > 0.0)) then
            delta_float_mass(i,j) = ISS%mass_shelf(i,j) - last_mass_shelf(i,j)
          else
            delta_float_mass(i,j) = 0.0
          endif
        enddo ; enddo
        delta_mass_shelf = global_area_integral(delta_float_mass, G, tmp_scale=US%RZ_to_kg_m2, &
                                                area=ISS%area_shelf_h) / CS%time_step
      else! first time step
        delta_mass_shelf = 0.0
      endif
    else
      if (CS%active_shelf_dynamics) then ! change in ice_shelf draft
        do j=js,je ; do i=is,ie
          last_h_shelf(i,j) = ISS%h_shelf(i,j) - time_step * ISS%dhdt_shelf(i,j)
        enddo ; enddo
        call change_in_draft(CS%dCS, G, last_h_shelf, ISS%h_shelf, delta_draft)

        !this currently assumes area_shelf_h is constant over the time step
        delta_mass_shelf = global_area_integral(delta_draft, G, tmp_scale=US%RZ_to_kg_m2, &
                                                area=ISS%area_shelf_h) &
                                                * CS%Rho_ocn / CS%time_step
      else ! ice shelf mass does not change
        delta_mass_shelf = 0.0
      endif
    endif

    ! average total melt flux over sponge area (ISOMIP/MISOMIP only) or open ocean (general case)
    do j=js,je ; do i=is,ie
      if (CS%constant_sea_level_misomip) then !for ismip/misomip only
        if (G%geoLonT(i,j) >= 790.0) then
          bal_frac(i,j) = max(1.0 - ISS%area_shelf_h(i,j) * G%IareaT(i,j), 0.0)
        else
          bal_frac(i,j) = 0.0
        endif
      elseif ((G%mask2dT(i,j) > 0.0) .and. (ISS%area_shelf_h(i,j) * G%IareaT(i,j) < 1.0)) then !general case
        bal_frac(i,j) = max(1.0 - ISS%area_shelf_h(i,j) * G%IareaT(i,j), 0.0)
      else
        bal_frac(i,j) = 0.0
      endif
    enddo ; enddo

    balancing_area = global_area_integral(bal_frac, G, area=G%areaT, tmp_scale=1.0)
    if (balancing_area > 0.0) then
      balancing_flux = ( global_area_integral(ISS%water_flux, G, tmp_scale=US%RZ_T_to_kg_m2s, &
                                              area=ISS%area_shelf_h) + &
                         delta_mass_shelf ) / balancing_area
    else
      balancing_flux = 0.0
    endif

    ! apply fluxes
    do j=js,je ; do i=is,ie
      if (bal_frac(i,j) > 0.0) then
        ! evap is negative, and vprec has units of [R Z T-1 ~> kg m-2 s-1]
        fluxes%vprec(i,j) = -balancing_flux
        fluxes%sens(i,j) = fluxes%vprec(i,j) * CS%Cp * CS%T0 ! [Q R Z T-1 ~> W m-2]
        fluxes%salt_flux(i,j) = fluxes%vprec(i,j) * CS%S0*1.0e-3*US%S_to_ppt ! [1e-3 S R Z T-1 ~> kgSalt m-2 s-1]
      endif
    enddo ; enddo

    if (CS%debug) then
      write(mesg,*) 'Balancing flux (kg/(m^2 s)), dt = ', balancing_flux*US%RZ_T_to_kg_m2s, US%T_to_s*CS%time_step
      call MOM_mesg(mesg)
      call MOM_forcing_chksum("After constant sea level", fluxes, G, CS%US, haloshift=0)
    endif

  endif ! constant_sea_level

end procedure add_shelf_flux
module procedure initialize_ice_shelf
  type(ocean_grid_type), pointer :: G  => NULL(), OG  => NULL() ! Pointers to grids for convenience.
  type(unit_scale_type), pointer :: US => NULL() ! Pointer to a structure containing
  type(ice_shelf_state), pointer :: ISS => NULL() !< A structure with elements that describe
  type(directories)  :: dirs
  type(dyn_horgrid_type), pointer :: dG => NULL()
  type(dyn_horgrid_type), pointer :: dG_in => NULL()
  real :: meltrate_conversion ! The conversion factor to use for in the melt rate diagnostic
  real :: dz_ocean_min_float ! The minimum ocean thickness above which the ice shelf is considered
  real :: cdrag         ! The drag coefficient at the ice-ocean interface [nondim]
  real :: drag_bg_vel   ! A background velocity used in the quadratic drag [Z T-1 ~> m s-1]
  logical :: new_sim, save_IC
# include "version_variable.h"
  character(len=200) :: IC_file, inputdir  ! Input file names or paths
  character(len=40)  :: mdl = "MOM_ice_shelf"  ! This module's name.
  integer :: i, j, is, ie, js, je, isd, ied, jsd, jed, Isdq, Iedq, Jsdq, Jedq
  integer :: wd_halos(2)
  logical :: showCallTree
  logical :: read_TideAmp, debug
  logical :: global_indexing
  character(len=240) :: Tideamp_file  ! Input file names
  character(len=80)  :: tideamp_var ! Input file variable names
  real    :: utide  ! A tidal velocity [L T-1 ~> m s-1]
  real    :: col_thick_melt_thresh ! An ocean column thickness below which iceshelf melting
  real, allocatable, dimension(:,:) :: tmp2d ! Temporary array for ice shelf input data [L T-1 ~> m s-1]
  real, allocatable, dimension(:,:) :: maskT ! Temporary array for the tracer points masks [nondim]
  type(surface), pointer :: sfc_state => NULL()
  type(vardesc) :: u_desc, v_desc
  if (associated(CS)) then
    call MOM_error(FATAL, "MOM_ice_shelf.F90, initialize_ice_shelf: "// &
                          "called with an associated control structure.")
    return
  endif
  allocate(CS)

  !   Go through all of the infrastructure initialization calls, since this is
  ! being treated as an independent component that just happens to use the
  ! MOM's grid and infrastructure.
  call Get_MOM_Input(dirs=dirs)

  call MOM_IS_diag_mediator_infrastructure_init()

  ! Determining the internal unit scaling factors for this run.
  call unit_scaling_init(param_file, CS%US)

  call get_param(param_file, mdl, "ROTATE_INDEX", CS%rotate_index, &
      "Enable rotation of the horizontal indices.", default=.false., &
      debuggingParam=.true.)

  call get_param(param_file, "MOM", "GLOBAL_INDEXING", global_indexing, &
                 "If true, use a global lateral indexing convention, so "//&
                 "that corresponding points on different processors have "//&
                 "the same index. This does not work with static memory.", &
                 default=.false., layoutParam=.true.)

  ! Set up the ice-shelf domain and grid
  wd_halos(:)=0
  allocate(CS%Grid_in)
  call MOM_domains_init(CS%Grid_in%domain, param_file, min_halo=wd_halos, symmetric=GRID_SYM_,&
                        domain_name='MOM_Ice_Shelf_in', US=CS%US)
  !allocate(CS%Grid_in%HI)
  !call hor_index_init(CS%Grid%Domain, CS%Grid%HI, param_file, &
  !     local_indexing=.not.global_indexing)
  call MOM_grid_init(CS%Grid_in, param_file, CS%US)

  if (CS%rotate_index) then
  !   ! TODO: Index rotation currently only works when index rotation does not
  !   !   change the MPI rank of each domain.  Resolving this will require a
  !   !   modification to FMS PE assignment.
  !   !   For now, we only permit single-core runs.

    if (num_PEs() /= 1) call MOM_error(FATAL, "Index rotation is only supported on one PE.")

    call get_param(param_file, mdl, "INDEX_TURNS", CS%turns, &
         "Number of counterclockwise quarter-turn index rotations.", &
         default=1, debuggingParam=.true.)
    ! NOTE: If indices are rotated, then CS%Grid and CS%Grid_in must both be initialized.
    !   If not rotated, then CS%Grid_in and CS%Ggrid are the same grid.
    call create_dyn_horgrid(dG_in, CS%Grid_in%HI)
    call clone_MOM_domain(CS%Grid_in%Domain, dG_in%Domain)
    call set_grid_metrics(dG_in, param_file, CS%US)
    ! Set up the bottom depth, dG_in%bathyT, either analytically or from file
    call MOM_initialize_topography(dG_in%bathyT, CS%Grid_in%max_depth, dG_in, param_file, CS%US)

    ! The use of maskT here sets all ice shelf points to be unmasked.
    allocate(maskT(dG_in%isd:dG_in%ied,dG_in%jsd:dG_in%jed), source=1.0)
    call initialize_masks(dG_in, param_file, CS%US, maskT=maskT)
    deallocate(maskT)

    call copy_dyngrid_to_MOM_grid(dG_in, CS%Grid_in, CS%US)

    ! Now set up the rotated ice-shelf grid.
    allocate(CS%Grid)
    call clone_MOM_domain(CS%Grid_in%Domain, CS%Grid%Domain, turns=CS%turns)
    call rotate_hor_index(CS%Grid_in%HI, CS%turns, CS%Grid%HI)
    call MOM_grid_init(CS%Grid, param_file, CS%US, CS%Grid%HI)
    call create_dyn_horgrid(dG, CS%Grid%HI)
    call rotate_dyngrid(dG_in, dG, CS%US, CS%turns)
    call copy_dyngrid_to_MOM_grid(dG, CS%Grid, CS%US)

    call destroy_dyn_horgrid(dG_in)
    call destroy_dyn_horgrid(dG)
  else
    CS%Grid => CS%Grid_in
    dG => NULL()
    call create_dyn_horgrid(dG, CS%Grid%HI)
    call clone_MOM_domain(CS%Grid%Domain, dG%Domain)
    call set_grid_metrics(dG, param_file, CS%US)
    ! Set up the bottom depth, dG%bathyT, either analytically or from file
    call MOM_initialize_topography(dG%bathyT, CS%Grid%max_depth, dG, param_file, CS%US)

    ! The use of maskT here sets all ice shelf points to be unmasked.
    allocate(maskT(dG%isd:dG%ied,dG%jsd:dG%jed), source=1.0)
    call initialize_masks(dG, param_file, CS%US, maskT=maskT)
    deallocate(maskT)

    call copy_dyngrid_to_MOM_grid(dG, CS%Grid, CS%US)
    call destroy_dyn_horgrid(dG)
  endif
  G => CS%Grid

  allocate(CS%diag)
  call MOM_IS_diag_mediator_init(G, CS%US, param_file, CS%diag, component='MOM_IceShelf')
  ! This call sets up the diagnostic axes. These are needed,
  ! e.g. to generate the target grids below.
  call set_IS_axes_info(G, CS%diag)


  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec
  isd = G%isd ; jsd = G%jsd ; ied = G%ied ; jed = G%jed
  Isdq = G%IsdB ; Iedq = G%IedB ; Jsdq = G%JsdB ; Jedq = G%JedB

  ! The ocean grid possibly uses different symmetry.
  if (associated(ocn_grid)) then ; CS%ocn_grid => ocn_grid
  else ; CS%ocn_grid => CS%grid ; endif

  ! Convenience pointers
  OG => CS%ocn_grid
  US => CS%US

  ! Are we being called from the solo ice-sheet driver? When called by the ocean
  ! model solo_ice_sheet_in is not preset.
  CS%solo_ice_sheet = .false.
  if (present(solo_ice_sheet_in)) CS%solo_ice_sheet = solo_ice_sheet_in

  !if (present(Time_in)) Time = Time_in


  CS%override_shelf_movement = .false. ; CS%active_shelf_dynamics = .false.

  call log_version(param_file, mdl, version, "")
  call get_param(param_file, mdl, "DEBUG", debug, default=.false.)
  call get_param(param_file, mdl, "DEBUG_IS", CS%debug, &
                 "If true, write verbose debugging messages for the ice shelf.", &
                 default=debug)
  call get_param(param_file, mdl, "DYNAMIC_SHELF_MASS", CS%shelf_mass_is_dynamic, &
                 "If true, the ice sheet mass can evolve with time.", &
                 default=.false.)
  if (CS%shelf_mass_is_dynamic) then
    call get_param(param_file, mdl, "OVERRIDE_SHELF_MOVEMENT", CS%override_shelf_movement, &
                 "If true, user provided code specifies the ice-shelf "//&
                 "movement instead of the dynamic ice model.", default=.false.)
    CS%active_shelf_dynamics = .not.CS%override_shelf_movement
    call get_param(param_file, mdl, "DATA_OVERRIDE_SHELF_FLUXES", &
                  CS%data_override_shelf_fluxes, &
                 "If true, the data override feature is used to write "//&
                 "the surface mass flux deposition. This option is only "//&
                 "available for MOSAIC grid types.", default=.false.)
    call get_param(param_file, mdl, "GROUNDING_LINE_INTERPOLATE", CS%GL_regularize, &
                 "If true, regularize the floatation condition at the "//&
                 "grounding line as in Goldberg Holland Schoof 2009.", default=.false.)
    call get_param(param_file, mdl, "GROUNDING_LINE_COUPLE", CS%GL_couple, &
                 "If true, let the floatation condition be determined by "//&
                 "ocean column thickness. This means that update_OD_ffrac "//&
                 "will be called.  GL_REGULARIZE and GL_COUPLE are exclusive.", &
                 default=.false., do_not_log=CS%GL_regularize)
    if (CS%GL_regularize) CS%GL_couple = .false.
    if (CS%solo_ice_sheet) CS%GL_couple = .false.
  endif

  call get_param(param_file, mdl, "SHELF_THERMO", CS%isthermo, &
                 "If true, use a thermodynamically interactive ice shelf.", &
                 default=.false.)
  call get_param(param_file, mdl, "LATENT_HEAT_FUSION", CS%Lat_fusion, &
                 "The latent heat of fusion.", units="J/kg", default=hlf, scale=US%J_kg_to_Q)
  call get_param(param_file, mdl, "SHELF_THREE_EQN", CS%threeeq, &
                 "If true, use the three equation expression of "//&
                 "consistency to calculate the fluxes at the ice-ocean "//&
                 "interface.", default=.true.)
  call get_param(param_file, mdl, "SHELF_INSULATOR", CS%insulator, &
                 "If true, the ice shelf is a perfect insulator "//&
                 "(no conduction).", default=.false.)
  call get_param(param_file, mdl, "MELTING_CUTOFF_DEPTH", CS%cutoff_depth, &
                 "Depth above which the melt is set to zero (it must be >= 0) "//&
                 "Default value won't affect the solution.", units="m", default=0.0, scale=US%m_to_Z)
  if (CS%cutoff_depth < 0.) &
    call MOM_error(WARNING,"Initialize_ice_shelf: MELTING_CUTOFF_DEPTH must be >= 0.")

  call get_param(param_file, mdl, "CONST_SEA_LEVEL", CS%constant_sea_level, &
                 "If true, apply evaporative, heat and salt fluxes in "//&
                 "the sponge region. This will avoid a large increase "//&
                 "in sea level. This option is needed for some of the "//&
                 "ISOMIP+ experiments (Ocean3 and Ocean4). "//&
                 "IMPORTANT: it is not currently possible to do "//&
                 "prefect restarts using this flag.", default=.false.)
  call get_param(param_file, mdl, "CONST_SEA_LEVEL_MISOMIP", CS%constant_sea_level_misomip, &
                 "If true, constant_sea_level fluxes are applied only over "//&
                 "the surface sponge cells from the ISOMIP/MISOMIP configuration", default=.false.)
  call get_param(param_file, mdl, "MIN_OCEAN_FLOAT_THICK", dz_ocean_min_float, &
                 "The minimum ocean thickness above which the ice shelf is considered to be "//&
                 "floating when CONST_SEA_LEVEL = True.", &
                 default=0.1, units="m", scale=US%m_to_Z, do_not_log=.not.CS%constant_sea_level)

  call get_param(param_file, mdl, "ISOMIP_S_SUR_SPONGE", CS%S0, &
                 "Surface salinity in the restoring region.", &
                default=33.8, units='ppt', scale=US%ppt_to_S, do_not_log=.true.)

  call get_param(param_file, mdl, "ISOMIP_T_SUR_SPONGE", CS%T0, &
                "Surface temperature in the restoring region.", &
                default=-1.9, units='degC', scale=US%degC_to_C, do_not_log=.true.)

  call get_param(param_file, mdl, "SHELF_3EQ_GAMMA", CS%const_gamma, &
                 "If true, user specifies a constant nondimensional heat-transfer coefficient "//&
                 "(GAMMA_T_3EQ), from which the default salt-transfer coefficient is set "//&
                 "as GAMMA_T_3EQ/35. This is used with SHELF_THREE_EQN.", default=.false.)
  call get_param(param_file, mdl, "SHELF_S_ROOT", CS%find_salt_root, &
                 "If SHELF_S_ROOT = True, salinity at the ice/ocean interface (Sbdry) "//&
                 "is computed from a quadratic equation. Otherwise, the previous "//&
                 "interactive method to estimate Sbdry is used.", &
                 default=.false., do_not_log=.not.CS%threeeq)
  if (.not.CS%threeeq) then
    call get_param(param_file, mdl, "SHELF_2EQ_GAMMA_T", CS%gamma_t, &
                 "If SHELF_THREE_EQN is false, this the fixed turbulent "//&
                 "exchange velocity at the ice-ocean interface.", &
                 units="m s-1", scale=US%m_to_Z*US%T_to_s, fail_if_missing=.true.)
  endif
  if (CS%const_gamma .or. CS%find_salt_root) then
    call get_param(param_file, mdl, "SHELF_3EQ_GAMMA_T", CS%Gamma_T_3EQ, &
                 "Nondimensional heat-transfer coefficient.", &
                  units="nondim", default=2.2e-2)
    call get_param(param_file, mdl, "SHELF_3EQ_GAMMA_S", CS%Gamma_S_3EQ, &
                 "Nondimensional salt-transfer coefficient.", &
                 default=CS%Gamma_T_3EQ/35.0, units="nondim")
  endif

  call get_param(param_file, mdl, "ICE_SHELF_MASS_FROM_FILE", &
                 CS%mass_from_file, "Read the mass of the "//&
                 "ice shelf (every time step) from a file.", default=.false.)

  if (CS%find_salt_root) then ! read liquidus coeffs.
    call get_param(param_file, mdl, "TFREEZE_S0_P0", CS%TFr_0_0, &
                 "this is the freezing potential temperature at "//&
                 "S=0, P=0.", units="degC", default=0.0, scale=US%degC_to_C, do_not_log=.true.)
    call get_param(param_file, mdl, "DTFREEZE_DS", CS%dTFr_dS, &
                 "this is the derivative of the freezing potential temperature with salinity.", &
                 units="degC psu-1", default=-0.054, scale=US%degC_to_C*US%S_to_ppt, do_not_log=.true.)
    call get_param(param_file, mdl, "DTFREEZE_DP", CS%dTFr_dp, &
                 "this is the derivative of the freezing potential temperature with pressure.", &
                 units="degC Pa-1", default=0.0, scale=US%degC_to_C*US%RL2_T2_to_Pa, do_not_log=.true.)
  endif

  call get_param(param_file, mdl, "G_EARTH", CS%g_Earth, &
                 "The gravitational acceleration of the Earth.", &
                 units="m s-2", default=9.80, scale=US%m_s_to_L_T**2*US%Z_to_m)
  call get_param(param_file, mdl, "C_P", CS%Cp, &
                 "The heat capacity of sea water, approximated as a constant. "//&
                 "The default value is from the TEOS-10 definition of conservative temperature.", &
                 units="J kg-1 K-1", default=3991.86795711963, scale=US%J_kg_to_Q*US%C_to_degC)
  call get_param(param_file, mdl, "RHO_0", CS%Rho_ocn, &
                 "The mean ocean density used with BOUSSINESQ true to "//&
                 "calculate accelerations and the mass for conservation "//&
                 "properties, or with BOUSSINESQ false to convert some "//&
                 "parameters from vertical units of m to kg m-2.", &
                 units="kg m-3", default=1035.0, scale=US%kg_m3_to_R)
  call get_param(param_file, mdl, "C_P_ICE", CS%Cp_ice, &
                 "The heat capacity of ice.", units="J kg-1 K-1", scale=US%J_kg_to_Q*US%C_to_degC, &
                 default=2.10e3)
  if (CS%constant_sea_level) CS%min_ocean_mass_float = dz_ocean_min_float*CS%Rho_ocn

  call get_param(param_file, mdl, "ICE_SHELF_FLUX_FACTOR", CS%flux_factor, &
                 "Non-dimensional factor applied to shelf thermodynamic fluxes.", &
                 units="none", default=1.0)

  call get_param(param_file, mdl, "KV_ICE", CS%kv_ice, &
                 "The viscosity of the ice.", &
                 units="m2 s-1", default=1.0e10, scale=US%Z_to_L**2*US%m_to_L**2*US%T_to_s)
  call get_param(param_file, mdl, "KV_MOLECULAR", CS%kv_molec, &
                 "The molecular kinematic viscosity of sea water at the freezing temperature.", &
                 units="m2 s-1", default=1.95e-6, scale=US%m2_s_to_Z2_T)
  call get_param(param_file, mdl, "ICE_SHELF_SALINITY", CS%Salin_ice, &
                 "The salinity of the ice inside the ice shelf.", &
                 units="psu", default=0.0, scale=US%ppt_to_S)
  call get_param(param_file, mdl, "ICE_SHELF_TEMPERATURE", CS%Temp_ice, &
                 "The temperature at the center of the ice shelf.", &
                 units="degC", default=-15.0, scale=US%degC_to_C)
  call get_param(param_file, mdl, "KD_SALT_MOLECULAR", CS%kd_molec_salt, &
                 "The molecular diffusivity of salt in sea water at the "//&
                 "freezing point.", units="m2 s-1", default=8.02e-10, scale=US%m2_s_to_Z2_T)
  call get_param(param_file, mdl, "KD_TEMP_MOLECULAR", CS%kd_molec_temp, &
                 "The molecular diffusivity of heat in sea water at the "//&
                 "freezing point.", units="m2 s-1", default=1.41e-7, scale=US%m2_s_to_Z2_T)
  call get_param(param_file, mdl, "DT_FORCING", CS%time_step, &
                 "The time step for changing forcing, coupling with other "//&
                 "components, or potentially writing certain diagnostics. "//&
                 "The default value is given by DT.", units="s", default=0.0, scale=US%s_to_T)

  call get_param(param_file, mdl, "COL_THICK_MELT_THRESHOLD", col_thick_melt_thresh, &
                 "The minimum ocean column thickness where melting is allowed.", &
                 units="m", scale=US%m_to_Z, default=0.0)
  CS%col_mass_melt_threshold =  CS%Rho_ocn * col_thick_melt_thresh

  call get_param(param_file, mdl, "READ_TIDEAMP", read_TIDEAMP, &
                 "If true, read a file (given by TIDEAMP_FILE) containing "//&
                 "the tidal amplitude with INT_TIDE_DISSIPATION.", default=.false.)
  call get_param(param_file, mdl, "ICE_SHELF_LINEAR_SHELF_FRAC", CS%Zeta_N, &
                 "Ratio of HJ99 stability constant xi_N (ratio of maximum "//&
                 "mixing length to planetary boundary layer depth in "//&
                 "neutrally stable conditions) to the von Karman constant", &
                 units="nondim", default=0.13)
  call get_param(param_file, mdl, "ICE_SHELF_VK_CNST", CS%Vk, &
                 "Von Karman constant.", &
                 units="nondim", default=0.40)
  call get_param(param_file, mdl, "ICE_SHELF_RC", CS%Rc, &
                 "Critical flux Richardson number for ice melt ", &
                 units="nondim", default=0.20)
  call get_param(param_file, mdl, "ICE_SHELF_USTAR_FROM_VEL_BUGFIX", CS%ustar_from_vel_bugfix, &
                 "Bug fix for ice-area weighting of squared ocean velocities "//&
                 "used to calculate friction velocity under ice shelves", default=.false.)
  call get_param(param_file, mdl, "ICE_SHELF_BUOYANCY_FLUX_ITT_BUGFIX", CS%buoy_flux_itt_bugfix, &
                 "Bug fix of buoyancy iteration", default=.true., old_name="ICE_SHELF_BUOYANCY_FLUX_ITT_BUG")
  call get_param(param_file, mdl, "ICE_SHELF_SALT_FLUX_ITT_BUGFIX", CS%salt_flux_itt_bugfix, &
                 "Bug fix of salt iteration", default=.true., old_name="ICE_SHELF_SALT_FLUX_ITT_BUG")
  call get_param(param_file, mdl, "ICE_SHELF_BUOYANCY_FLUX_ITT_THRESHOLD", CS%buoy_flux_tol, &
                 "Convergence criterion of Newton's method for ice shelf "//&
                 "buoyancy iteration.", units="nondim", default=1.0e-4)

  if (PRESENT(sfc_state_in)) then
    ! assuming frazil is enabled in ocean. This could break some configurations?
    call allocate_surface_state(sfc_state_in, CS%Grid_in, use_temperature=.true., &
          do_integrals=.true., omit_frazil=.false., use_iceshelves=.true.)
    if (CS%rotate_index) then
      allocate(sfc_state)
      call rotate_surface_state(sfc_state_in, sfc_state, CS%Grid, CS%turns)
    else
      sfc_state => sfc_state_in
    endif
  endif


  call safe_alloc_ptr(CS%utide,isd,ied,jsd,jed) ; CS%utide(:,:) = 0.0

  if (read_TIDEAMP) then
    call get_param(param_file, mdl, "TIDEAMP_FILE", TideAmp_file, &
                 "The path to the file containing the spatially varying tidal amplitudes.", &
                 default="tideamp.nc")
    call get_param(param_file, mdl, "TIDEAMP_VARNAME", tideamp_var, &
                 "The name of the tidal amplitude variable in the input file.", &
                 default="tideamp")
     call get_param(param_file, mdl, "INPUTDIR", inputdir, default=".")
    inputdir = slasher(inputdir)
    TideAmp_file = trim(inputdir) // trim(TideAmp_file)
    if (CS%rotate_index) then
      allocate(tmp2d(CS%Grid_in%isd:CS%Grid_in%ied,CS%Grid_in%jsd:CS%Grid_in%jed), source=0.0)
      call MOM_read_data(TideAmp_file, tideamp_var, tmp2d, CS%Grid_in%domain, timelevel=1, scale=US%m_s_to_L_T)
      call rotate_array(tmp2d, CS%turns, CS%utide)
      deallocate(tmp2d)
    else
      call MOM_read_data(TideAmp_file, tideamp_var, CS%utide, CS%Grid%domain, timelevel=1, scale=US%m_s_to_L_T)
    endif
  else
    call get_param(param_file, mdl, "UTIDE", utide, &
                 "The constant tidal amplitude used with INT_TIDE_DISSIPATION.", &
                 units="m s-1", default=0.0 , scale=US%m_s_to_L_T)
    CS%utide(:,:) = utide
  endif

  call EOS_init(param_file, CS%eqn_of_state, US)

  !! new parameters that need to be in MOM_input

  if (CS%active_shelf_dynamics) then

    call get_param(param_file, mdl, "DENSITY_ICE", CS%density_ice, &
                 "A typical density of ice.", units="kg m-3", default=917.0, scale=US%kg_m3_to_R)

    call get_param(param_file, mdl, "INPUT_FLUX_ICE_SHELF", CS%input_flux, &
                 "volume flux at upstream boundary", units="m2 s-1", default=0., scale=US%m_to_Z*US%m_s_to_L_T)
    call get_param(param_file, mdl, "INPUT_THICK_ICE_SHELF", CS%input_thickness, &
                 "flux thickness at upstream boundary", units="m", default=1000., scale=US%m_to_Z)
  else
    ! This is here because of inconsistent defaults.  I don't know why.  RWH
    call get_param(param_file, mdl, "DENSITY_ICE", CS%density_ice, &
                 "A typical density of ice.", units="kg m-3", default=900.0, scale=US%kg_m3_to_R)
  endif
  call get_param(param_file, mdl, "MIN_THICKNESS_SIMPLE_CALVE", &
                CS%min_thickness_simple_calve, &
                 "Min thickness rule for the very simple calving law",&
                 units="m", default=0.0, scale=US%m_to_Z)

  call get_param(param_file, mdl, "USTAR_SHELF_BG", CS%ustar_bg, &
                 "The minimum value of ustar under ice shelves.", &
                 units="m s-1", default=0.0, scale=US%m_to_Z*US%T_to_s)
  call get_param(param_file, mdl, "CDRAG_SHELF", cdrag, &
       "CDRAG is the drag coefficient relating the magnitude of "//&
       "the velocity field to the surface stress.", units="nondim", &
       default=0.003)
  CS%cdrag = cdrag
  if (CS%ustar_bg <= 0.0) then
    call get_param(param_file, mdl, "DRAG_BG_VEL_SHELF", drag_bg_vel, &
                 "DRAG_BG_VEL is either the assumed bottom velocity (with "//&
                 "LINEAR_DRAG) or an unresolved  velocity that is "//&
                 "combined with the resolved velocity to estimate the "//&
                 "velocity magnitude.", units="m s-1", default=0.0, scale=US%m_to_Z*US%T_to_s)
    if (CS%cdrag*drag_bg_vel > 0.0) CS%ustar_bg = sqrt(CS%cdrag)*drag_bg_vel
  endif
  call get_param(param_file, mdl, "USTAR_SHELF_FROM_VEL", CS%ustar_shelf_from_vel, &
                 "If true, use the surface velocities to set the friction velocity under ice "//&
                 "shelves instead of using the previous values of the stresses.", &
                 default=.true.)
  call get_param(param_file, mdl, "USTAR_SHELF_MAX", CS%ustar_max, &
                 "The maximum value of ustar under ice shelves, or a negative value for no limit.", &
                 units="m s-1", default=-1.0, scale=US%m_to_Z*US%T_to_s, &
                 do_not_log=CS%ustar_shelf_from_vel)

  if (present(calve_ice_shelf_bergs)) CS%calve_ice_shelf_bergs=calve_ice_shelf_bergs

  ! Allocate and initialize state variables to default values
  call ice_shelf_state_init(CS%ISS, CS%grid)
  ISS => CS%ISS

  new_sim = .false.
  if ((dirs%input_filename(1:1) == 'n') .and. &
      (LEN_TRIM(dirs%input_filename) == 1)) new_sim = .true.

  ISS%area_shelf_h(:,:)=0.0
  ISS%h_shelf(:,:)=0.0
  ISS%hmask(:,:)=0.0
  ISS%mass_shelf(:,:)=0.0

  if (CS%override_shelf_movement .and. CS%mass_from_file) then

    ! initialize the ids for reading shelf mass from a netCDF
    call initialize_shelf_mass(G, param_file, CS, ISS)

    if (new_sim) then
      ! new simulation, initialize ice thickness as in the static case
      call initialize_ice_thickness(ISS%h_shelf, ISS%area_shelf_h, ISS%hmask, ISS%melt_mask, CS%Grid, CS%Grid_in, &
                                    US, param_file, CS%rotate_index, CS%turns)

    ! next make sure mass is consistent with thickness
      do j=G%jsd,G%jed ; do i=G%isd,G%ied
        if ((ISS%hmask(i,j) == 1) .or. (ISS%hmask(i,j) == 2) .or. (ISS%hmask(i,j)==3)) then
          ISS%mass_shelf(i,j) = ISS%h_shelf(i,j)*CS%density_ice
        endif
      enddo ; enddo

      if (CS%min_thickness_simple_calve > 0.0) &
        call ice_shelf_min_thickness_calve(G, ISS%h_shelf, ISS%area_shelf_h, ISS%hmask, &
                                           CS%min_thickness_simple_calve)
    endif
  endif

  if (CS%active_shelf_dynamics) then
    ! the only reason to initialize boundary conds is if the shelf is dynamic - MJH

    ! call initialize_ice_shelf_boundary ( CS%u_face_mask_bdry, CS%v_face_mask_bdry, &
    !                                      CS%u_flux_bdry_val, CS%v_flux_bdry_val, &
    !                                      CS%u_bdry_val, CS%v_bdry_val, CS%h_bdry_val, &
    !                                      ISS%hmask, G, param_file)

  endif

  ! Set up the restarts.

  call restart_init(param_file, CS%restart_CSp, "Shelf.res")
  call register_restart_field(ISS%mass_shelf, "shelf_mass", .true., CS%restart_CSp, &
                              "Ice shelf mass", "kg m-2", conversion=US%RZ_to_kg_m2)
  call register_restart_field(ISS%area_shelf_h, "shelf_area", .true., CS%restart_CSp, &
                              "Ice shelf area in cell", "m2", conversion=US%L_to_m**2)
  call register_restart_field(ISS%h_shelf, "h_shelf", .true., CS%restart_CSp, &
                              "ice sheet/shelf thickness", "m", conversion=US%Z_to_m)
  call register_restart_field(ISS%melt_mask, "melt_mask", .false., CS%restart_CSp, &
                              "Mask that is >0 where ice-shelf melting is allowed", "none")
  if (CS%calve_ice_shelf_bergs) then
    call register_restart_field(ISS%calving, "shelf_calving", .true., CS%restart_CSp, &
                                "Calving flux from ice shelf into icebergs", "kg m-2", conversion=US%RZ_to_kg_m2)
    call register_restart_field(ISS%calving_hflx, "shelf_calving_hflx", .true., CS%restart_CSp, &
                                "Calving heat flux from ice shelf into icebergs", "W m-2", conversion=US%QRZ_T_to_W_m2)
  endif

  if (PRESENT(sfc_state_in)) then
    if (allocated(sfc_state%taux_shelf) .and. allocated(sfc_state%tauy_shelf)) then
      u_desc = var_desc("taux_shelf", "Pa", "the zonal stress on the ocean under ice shelves", &
            hor_grid='Cu',z_grid='1')
      v_desc = var_desc("tauy_shelf", "Pa", "the meridional stress on the ocean under ice shelves", &
            hor_grid='Cv',z_grid='1')
      call register_restart_pair(sfc_state%taux_shelf, sfc_state%tauy_shelf, u_desc, v_desc, &
            .false., CS%restart_CSp, conversion=US%RLZ_T2_to_Pa)
    endif
  endif

  if (CS%active_shelf_dynamics) then
    call register_restart_field(ISS%hmask, "h_mask", .true., CS%restart_CSp, &
                                "ice sheet/shelf thickness mask" ,"none")
  endif

  if (CS%active_shelf_dynamics) then
    ! Allocate CS%dCS and specify additional restarts for ice shelf dynamics
    call register_ice_shelf_dyn_restarts(CS%Grid_in, US, param_file, CS%dCS, CS%restart_CSp)
  endif

  !GMM - I think we do not need to save ustar_shelf and iceshelf_melt in the restart file
  !if (.not. CS%solo_ice_sheet) then
  !  call register_restart_field(fluxes%ustar_shelf, "ustar_shelf", .false., CS%restart_CSp, &
  !                              "Friction velocity under ice shelves", "m s-1", conversion=US%Z_to_m*US%s_to_T)
  !endif

  CS%restart_output_dir = dirs%restart_output_dir

  if (present(fluxes_in)) then
     call initialize_ice_shelf_fluxes(CS, ocn_grid, US, fluxes_in)
     call register_restart_field(fluxes_in%shelf_sfc_mass_flux, "sfc_mass_flux", .true., CS%restart_CSp, &
        "ice shelf surface mass flux deposition from atmosphere", &
        'kg m-2 s-1', conversion=US%RZ_T_to_kg_m2s)
  endif

  if (new_sim .and. (.not. (CS%override_shelf_movement .and. CS%mass_from_file))) then
    ! This model is initialized internally or from a file.
    call initialize_ice_thickness(ISS%h_shelf, ISS%area_shelf_h, ISS%hmask, ISS%melt_mask, CS%Grid, CS%Grid_in, &
                                  US, param_file, CS%rotate_index, CS%turns)
    ! next make sure mass is consistent with thickness
    do j=G%jsd,G%jed ; do i=G%isd,G%ied
      if ((ISS%hmask(i,j) == 1) .or. (ISS%hmask(i,j) == 2) .or. (ISS%hmask(i,j) == 3)) then
        ISS%mass_shelf(i,j) = ISS%h_shelf(i,j)*CS%density_ice
      endif
    enddo ; enddo
    if (CS%debug) then
      call hchksum(ISS%mass_shelf, "IS init: mass_shelf", G%HI, haloshift=0, unscale=US%RZ_to_kg_m2)
      call hchksum(ISS%area_shelf_h, "IS init: area_shelf", G%HI, haloshift=0, unscale=US%L_to_m*US%L_to_m)
      call hchksum(ISS%hmask, "IS init: hmask", G%HI, haloshift=0)
    endif

  ! else ! Previous block for new_sim=.T., this block restores the state.
  elseif (.not.new_sim) then
    ! This line calls a subroutine that reads the initial conditions from a restart file.
    call MOM_mesg("MOM_ice_shelf.F90, initialize_ice_shelf: Restoring ice shelf from file.")
    call restore_state(dirs%input_filename, dirs%restart_input_dir, Time, G, CS%restart_CSp)

  endif ! .not. new_sim

!  do j=G%jsc,G%jec ; do i=G%isc,G%iec
!    ISS%area_shelf_h(i,j) = ISS%area_shelf_h(i,j)*G%mask2dT(i,j)
!  enddo ; enddo

  id_clock_shelf = cpu_clock_id('Ice shelf', grain=CLOCK_COMPONENT)
  id_clock_pass = cpu_clock_id(' Ice shelf halo updates', grain=CLOCK_ROUTINE)

  call cpu_clock_begin(id_clock_pass)
  call pass_var(ISS%area_shelf_h, G%domain, complete=.false.)
  call pass_var(ISS%h_shelf, G%domain, complete=.false.)
  call pass_var(ISS%mass_shelf, G%domain, complete=.false.)
  call pass_var(ISS%hmask, G%domain, complete=.false.)
  call pass_var(G%bathyT, G%domain)
  call cpu_clock_end(id_clock_pass)

  do j=jsd,jed ; do i=isd,ied
    if (ISS%area_shelf_h(i,j) > G%areaT(i,j)) then
      call MOM_error(WARNING,"Initialize_ice_shelf: area_shelf_h exceeds G%areaT.")
      ISS%area_shelf_h(i,j) = G%areaT(i,j)
    endif
  enddo ; enddo

  if (CS%debug) then
    call hchksum(ISS%area_shelf_h, "IS init: area_shelf_h", G%HI, haloshift=0, unscale=US%L_to_m*US%L_to_m)
  endif

  CS%Time = Time

  if (CS%active_shelf_dynamics .and. .not.CS%isthermo) then
    ISS%water_flux(:,:) = 0.0
  endif

  if (CS%shelf_mass_is_dynamic) &
    call initialize_ice_shelf_dyn(param_file, Time, ISS, CS%dCS, G, US, CS%diag, new_sim, CS%Cp_ice, &
    Time_init, directory, solo_ice_sheet_in)

  call fix_restart_unit_scaling(US, unscaled=.true.)

  call get_param(param_file, mdl, "SAVE_INITIAL_CONDS", save_IC, &
                 "If true, save the ice shelf initial conditions.", default=.false.)
  if (save_IC) call get_param(param_file, mdl, "SHELF_IC_OUTPUT_FILE", IC_file,&
                 "The name-root of the output file for the ice shelf initial conditions.", &
                 default="MOM_Shelf_IC")

  if (save_IC .and. .not.((dirs%input_filename(1:1) == 'r') .and. &
                          (LEN_TRIM(dirs%input_filename) == 1))) then
    showCallTree = callTree_showQuery()
    if (showCallTree) call callTree_waypoint("About to call save_restart (MOM_ice_shelf)")
    call save_restart(dirs%output_directory, CS%Time, CS%Grid_in, CS%restart_CSp, &
                      filename=IC_file, write_ic=.true.)
    if (showCallTree) call callTree_waypoint("Done with call to save_restart (MOM_ice_shelf)")
  endif

  CS%id_area_shelf_h = register_diag_field('ice_shelf_model', 'area_shelf_h', CS%diag%axesT1, CS%Time, &
      'Ice Shelf Area in cell', 'meter2', conversion=US%L_to_m**2)
  CS%id_shelf_mass = register_diag_field('ice_shelf_model', 'shelf_mass', CS%diag%axesT1, CS%Time, &
      'mass of shelf', 'kg/m^2', conversion=US%RZ_to_kg_m2)
  CS%id_h_shelf = register_diag_field('ice_shelf_model', 'h_shelf', CS%diag%axesT1, CS%Time, &
      'ice shelf thickness', 'm', conversion=US%Z_to_m)
  CS%id_dhdt_shelf = register_diag_field('ice_shelf_model', 'dhdt_shelf', CS%diag%axesT1, CS%Time, &
      'change in ice shelf thickness over time', 'm s-1', conversion=US%Z_to_m*US%s_to_T)
  CS%id_mass_flux = register_diag_field('ice_shelf_model', 'mass_flux', CS%diag%axesT1,&
      CS%Time, 'Total mass flux of freshwater across the ice-ocean interface.', &
      'kg/s', conversion=US%RZ_T_to_kg_m2s*US%L_to_m**2)

  if (CS%const_gamma) then ! use ISOMIP+ eq. with rho_fw = 1000. kg m-3
    meltrate_conversion = 86400.0*365.0*US%Z_to_m*US%s_to_T / (1000.0*US%kg_m3_to_R)
  else ! use original eq.
    meltrate_conversion = 86400.0*365.0*US%Z_to_m*US%s_to_T / CS%density_ice
  endif
  CS%id_melt = register_diag_field('ice_shelf_model', 'melt', CS%diag%axesT1, CS%Time, &
      'Ice Shelf Melt Rate', 'm yr-1', conversion=meltrate_conversion)
  CS%id_thermal_driving = register_diag_field('ice_shelf_model', 'thermal_driving', CS%diag%axesT1, CS%Time, &
      'pot. temp. in the boundary layer minus freezing pot. temp. at the ice-ocean interface.', &
      'Celsius', conversion=US%C_to_degC)
  CS%id_haline_driving = register_diag_field('ice_shelf_model', 'haline_driving', CS%diag%axesT1, CS%Time, &
      'salinity in the boundary layer minus salinity at the ice-ocean interface.', &
      'psu', conversion=US%S_to_ppt)
  CS%id_Sbdry = register_diag_field('ice_shelf_model', 'sbdry', CS%diag%axesT1, CS%Time, &
      'salinity at the ice-ocean interface.', 'psu', conversion=US%S_to_ppt)
  CS%id_u_ml = register_diag_field('ice_shelf_model', 'u_ml', CS%diag%axesCu1, CS%Time, &
      'Eastward vel. in the boundary layer (used to compute ustar)', 'm s-1', conversion=US%L_T_to_m_s)
  CS%id_v_ml = register_diag_field('ice_shelf_model', 'v_ml', CS%diag%axesCv1, CS%Time, &
      'Northward vel. in the boundary layer (used to compute ustar)', 'm s-1', conversion=US%L_T_to_m_s)
  CS%id_exch_vel_s = register_diag_field('ice_shelf_model', 'exch_vel_s', CS%diag%axesT1, CS%Time, &
      'Sub-shelf salinity exchange velocity', 'm s-1', conversion=US%Z_to_m*US%s_to_T)
  CS%id_exch_vel_t = register_diag_field('ice_shelf_model', 'exch_vel_t', CS%diag%axesT1, CS%Time, &
      'Sub-shelf thermal exchange velocity', 'm s-1' , conversion=US%Z_to_m*US%s_to_T)
  CS%id_tfreeze = register_diag_field('ice_shelf_model', 'tfreeze', CS%diag%axesT1, CS%Time, &
      'In Situ Freezing point at ice shelf interface', 'degC', conversion=US%C_to_degC)
  CS%id_tfl_shelf = register_diag_field('ice_shelf_model', 'tflux_shelf', CS%diag%axesT1, CS%Time, &
      'Heat conduction into ice shelf', 'W m-2', conversion=-US%QRZ_T_to_W_m2)
  CS%id_ustar_shelf = register_diag_field('ice_shelf_model', 'ustar_shelf', CS%diag%axesT1, CS%Time, &
      'Fric vel under shelf', 'm/s', conversion=US%Z_to_m*US%s_to_T)
  CS%id_frazil = register_diag_field('ice_shelf_model', 'frazil', CS%diag%axesT1, CS%Time, &
     'Frazil heat rejected by the ocean', 'J m-2', conversion=US%Q_to_J_kg*US%RZ_to_kg_m2)
  if (CS%active_shelf_dynamics) then
    CS%id_h_mask = register_diag_field('ice_shelf_model', 'h_mask', CS%diag%axesT1, CS%Time, &
       'ice shelf thickness mask', 'none', conversion=1.0)
  endif

  CS%id_shelf_sfc_mass_flux = register_diag_field('ice_shelf_model', 'sfc_mass_flux', CS%diag%axesT1, CS%Time, &
     'ice shelf surface mass flux deposition from atmosphere', &
     'kg m-2 s-1', conversion=US%RZ_T_to_kg_m2s)

  ! Scalars (area integrated over all ice sheets)
  CS%id_vaf = register_scalar_field('ice_shelf_model', 'int_vaf', CS%diag%axesT1, CS%Time, &
      'Area integrated ice sheet volume above floatation', 'm3', conversion=US%Z_to_m*US%L_to_m**2)
  CS%id_adott = register_scalar_field('ice_shelf_model', 'int_a', CS%diag%axesT1, CS%Time, &
      'Area integrated change in ice-sheet thickness ' //&
      'due to surface accum+melt during a DT_THERM time step', 'm3', conversion=US%Z_to_m*US%L_to_m**2)
  CS%id_g_adott = register_scalar_field('ice_shelf_model', 'int_a_ground', CS%diag%axesT1, CS%Time, &
      'Area integrated change in grounded ice-sheet thickness ' //&
      'due to surface accum+melt during a DT_THERM time step', 'm3', conversion=US%Z_to_m*US%L_to_m**2)
  CS%id_f_adott = register_scalar_field('ice_shelf_model', 'int_a_float', CS%diag%axesT1, CS%Time, &
      'Area integrated change in floating ice-shelf thickness ' //&
      'due to surface accum+melt during a DT_THERM time step', 'm3', conversion=US%Z_to_m*US%L_to_m**2)
  CS%id_bdott = register_scalar_field('ice_shelf_model', 'int_b', CS%diag%axesT1, CS%Time, &
      'Area integrated change in floating ice-shelf thickness '//&
      'due to basal accum+melt during a DT_THERM time step', 'm3', conversion=US%Z_to_m*US%L_to_m**2)
  CS%id_bdott_melt = register_scalar_field('ice_shelf_model', 'int_b_melt', CS%diag%axesT1, CS%Time, &
      'Area integrated basal melt over ice shelves during a DT_THERM time step', &
      units='m3', conversion=US%Z_to_m*US%L_to_m**2)
  CS%id_bdott_accum = register_scalar_field('ice_shelf_model', 'int_b_accum', CS%diag%axesT1, CS%Time, &
      'Area integrated basal accumulation over ice shelves during a DT_THERM a time step', &
      units='m3', conversion=US%Z_to_m*US%L_to_m**2)
  CS%id_t_area = register_scalar_field('ice_shelf_model', 'tot_area', CS%diag%axesT1, CS%Time, &
      'Total ice-sheet area', 'm2', conversion=US%L_to_m**2)
  CS%id_f_area = register_scalar_field('ice_shelf_model', 'tot_area_float', CS%diag%axesT1, CS%Time, &
      'Total area of floating ice shelves', 'm2', conversion=US%L_to_m**2)
  CS%id_g_area = register_scalar_field('ice_shelf_model', 'tot_area_ground', CS%diag%axesT1, CS%Time, &
      'Total area of grounded ice sheets', 'm2', conversion=US%L_to_m**2)
  !scalars (area integrated rates over all ice sheets)
  CS%id_dvafdt = register_scalar_field('ice_shelf_model', 'int_vafdot', CS%diag%axesT1, CS%Time, &
      'Area integrated rate of change in ice-sheet volume above floatation', &
      units='m3 s-1', conversion=US%Z_to_m*US%L_to_m**2*US%s_to_T)
  CS%id_adot = register_scalar_field('ice_shelf_model', 'int_adot', CS%diag%axesT1, CS%Time, &
      'Area integrated rate of change in ice-sheet thickness due to surface accum+melt', &
      units='m3 s-1', conversion=US%Z_to_m*US%L_to_m**2*US%s_to_T)
  CS%id_g_adot = register_scalar_field('ice_shelf_model', 'int_adot_ground', CS%diag%axesT1, CS%Time, &
      'Area integrated rate of change in grounded ice-sheet thickness due to surface accum+melt', &
      units='m3 s-1', conversion=US%Z_to_m*US%L_to_m**2*US%s_to_T)
  CS%id_f_adot = register_scalar_field('ice_shelf_model', 'int_adot_float', CS%diag%axesT1, CS%Time, &
      'Area integrated rate of change in floating ice-shelf thickness due to surface accum+melt', &
      units='m3 s-1', conversion=US%Z_to_m*US%L_to_m**2*US%s_to_T)
  CS%id_bdot = register_scalar_field('ice_shelf_model', 'int_bdot', CS%diag%axesT1, CS%Time, &
      'Area integrated rate of change in ice-shelf thickness due to basal accum+melt', &
      units='m3 s-1', conversion=US%Z_to_m*US%L_to_m**2*US%s_to_T)
  CS%id_bdot_melt = register_scalar_field('ice_shelf_model', 'int_bdot_melt', CS%diag%axesT1, CS%Time, &
      'Area integrated basal melt rate over ice shelves', &
      units='m3 s-1', conversion=US%Z_to_m*US%L_to_m**2*US%s_to_T)
  CS%id_bdot_accum = register_scalar_field('ice_shelf_model', 'int_bdot_accum', CS%diag%axesT1, CS%Time, &
      'Area integrated basal accumulation rate over ice shelves', &
      units='m3 s-1', conversion=US%Z_to_m*US%L_to_m**2*US%s_to_T)

  !scalars (area integrated over the Antarctic ice sheet)
  CS%id_Ant_vaf = register_scalar_field('ice_shelf_model', 'int_vaf_A', CS%diag%axesT1, CS%Time, &
      'Area integrated Antarctic ice sheet volume above floatation', 'm3', conversion=US%Z_to_m*US%L_to_m**2)
  CS%id_Ant_adott = register_scalar_field('ice_shelf_model', 'int_a_A', CS%diag%axesT1, CS%Time, &
      'Area integrated (Antarctic ice sheet) change in ice-sheet thickness ' //&
      'due to surface accum+melt during a DT_THERM time step', 'm3', conversion=US%Z_to_m*US%L_to_m**2)
  CS%id_Ant_g_adott = register_scalar_field('ice_shelf_model', 'int_a_ground_A', CS%diag%axesT1, CS%Time, &
      'Area integrated change in Antarctic grounded ice-sheet thickness ' //&
      'due to surface accum+melt during a DT_THERM time step', 'm3', conversion=US%Z_to_m*US%L_to_m**2)
  CS%id_Ant_f_adott = register_scalar_field('ice_shelf_model', 'int_a_float_A', CS%diag%axesT1, CS%Time, &
      'Area integrated change in Antarctic floating ice-shelf thickness ' //&
      'due to surface accum+melt during a DT_THERM time step', 'm3', conversion=US%Z_to_m*US%L_to_m**2)
  CS%id_Ant_bdott = register_scalar_field('ice_shelf_model', 'int_b_A', CS%diag%axesT1, CS%Time, &
      'Area integrated change in Antarctic floating ice-shelf thickness '//&
      'due to basal accum+melt during a DT_THERM time step', 'm3', conversion=US%Z_to_m*US%L_to_m**2)
  CS%id_Ant_bdott_melt = register_scalar_field('ice_shelf_model', 'int_b_melt_A', CS%diag%axesT1, CS%Time, &
      'Area integrated basal melt over Antarctic ice shelves during a DT_THERM time step', &
      units='m3', conversion=US%Z_to_m*US%L_to_m**2)
  CS%id_Ant_bdott_accum = register_scalar_field('ice_shelf_model', 'int_b_accum_A', CS%diag%axesT1, CS%Time, &
      'Area integrated basal accumulation over Antarctic ice shelves during a DT_THERM a time step', &
      units='m3', conversion=US%Z_to_m*US%L_to_m**2)
  CS%id_Ant_t_area = register_scalar_field('ice_shelf_model', 'tot_area_A', CS%diag%axesT1, CS%Time, &
      'Total area of Antarctic ice sheet', 'm2', conversion=US%L_to_m**2)
  CS%id_Ant_f_area = register_scalar_field('ice_shelf_model', 'tot_area_float_A', CS%diag%axesT1, CS%Time, &
      'Total area of Antarctic floating ice shelves', 'm2', conversion=US%L_to_m**2)
  CS%id_Ant_g_area = register_scalar_field('ice_shelf_model', 'tot_area_ground_A', CS%diag%axesT1, CS%Time, &
      'Total area of Antarctic grounded ice sheet', 'm2', conversion=US%L_to_m**2)
  !scalars (area integrated rates over the Antarctic ice sheet)
  CS%id_Ant_dvafdt = register_scalar_field('ice_shelf_model', 'int_vafdot_A', CS%diag%axesT1, CS%Time, &
      'Area integrated rate of change in Antarctic ice-sheet volume above floatation', &
      units='m3 s-1', conversion=US%Z_to_m*US%L_to_m**2*US%s_to_T)
  CS%id_Ant_adot = register_scalar_field('ice_shelf_model', 'int_adot_A', CS%diag%axesT1, CS%Time, &
      'Area integrated rate of change in Antarctic ice-sheet thickness due to surface accum+melt', &
      units='m3 s-1', conversion=US%Z_to_m*US%L_to_m**2*US%s_to_T)
  CS%id_Ant_g_adot = register_scalar_field('ice_shelf_model', 'int_adot_ground_A', CS%diag%axesT1, CS%Time, &
      'Area integrated rate of change in Antarctic grounded ice-sheet thickness due to surface accum+melt', &
      units='m3 s-1', conversion=US%Z_to_m*US%L_to_m**2*US%s_to_T)
  CS%id_Ant_f_adot = register_scalar_field('ice_shelf_model', 'int_adot_float_A', CS%diag%axesT1, CS%Time, &
      'Area integrated rate of change in Antarctic floating ice-shelf thickness due to surface accum+melt', &
      units='m3 s-1', conversion=US%Z_to_m*US%L_to_m**2*US%s_to_T)
  CS%id_Ant_bdot = register_scalar_field('ice_shelf_model', 'int_bdot_A', CS%diag%axesT1, CS%Time, &
      'Area integrated rate of change in Antarctic ice-shelf thickness due to basal accum+melt', &
      units='m3 s-1', conversion=US%Z_to_m*US%L_to_m**2*US%s_to_T)
  CS%id_Ant_bdot_melt = register_scalar_field('ice_shelf_model', 'int_bdot_melt_A', CS%diag%axesT1, CS%Time, &
      'Area integrated basal melt rate over Antarctic ice shelves', &
      units='m3 s-1', conversion=US%Z_to_m*US%L_to_m**2*US%s_to_T)
  CS%id_Ant_bdot_accum = register_scalar_field('ice_shelf_model', 'int_bdot_accum_A', CS%diag%axesT1, CS%Time, &
      'Area integrated basal accumulation rate over Antarctic ice shelves', &
      units='m3 s-1', conversion=US%Z_to_m*US%L_to_m**2*US%s_to_T)

  !scalars (area integrated over the Greenland ice sheet)
  CS%id_Gr_vaf = register_scalar_field('ice_shelf_model', 'int_vaf_G', CS%diag%axesT1, CS%Time, &
      'Area integrated Greenland ice sheet volume above floatation', 'm3', conversion=US%Z_to_m*US%L_to_m**2)
  CS%id_Gr_adott = register_scalar_field('ice_shelf_model', 'int_a_G', CS%diag%axesT1, CS%Time, &
      'Area integrated (Greenland ice sheet) change in ice-sheet thickness ' //&
      'due to surface accum+melt during a DT_THERM time step', 'm3', conversion=US%Z_to_m*US%L_to_m**2)
  CS%id_Gr_g_adott = register_scalar_field('ice_shelf_model', 'int_a_ground_G', CS%diag%axesT1, CS%Time, &
      'Area integrated change in Greenland grounded ice-sheet thickness ' //&
      'due to surface accum+melt during a DT_THERM time step', 'm3', conversion=US%Z_to_m*US%L_to_m**2)
  CS%id_Gr_f_adott = register_scalar_field('ice_shelf_model', 'int_a_float_G', CS%diag%axesT1, CS%Time, &
      'Area integrated change in Greenland floating ice-shelf thickness ' //&
      'due to surface accum+melt during a DT_THERM time step', 'm3', conversion=US%Z_to_m*US%L_to_m**2)
  CS%id_Gr_bdott = register_scalar_field('ice_shelf_model', 'int_b_G', CS%diag%axesT1, CS%Time, &
      'Area integrated change in Greenland floating ice-shelf thickness '//&
      'due to basal accum+melt during a DT_THERM time step', 'm3', conversion=US%Z_to_m*US%L_to_m**2)
  CS%id_Gr_bdott_melt = register_scalar_field('ice_shelf_model', 'int_b_melt_G', CS%diag%axesT1, CS%Time, &
      'Area integrated basal melt over Greenland ice shelves during a DT_THERM time step', &
      units='m3', conversion=US%Z_to_m*US%L_to_m**2)
  CS%id_Gr_bdott_accum = register_scalar_field('ice_shelf_model', 'int_b_accum_G', CS%diag%axesT1, CS%Time, &
      'Area integrated basal accumulation over Greenland ice shelves during a DT_THERM a time step', &
      units='m3', conversion=US%Z_to_m*US%L_to_m**2)
  CS%id_Gr_t_area = register_scalar_field('ice_shelf_model', 'tot_area_G', CS%diag%axesT1, CS%Time, &
      'Total area of Greenland ice sheet', 'm2', conversion=US%L_to_m**2)
  CS%id_Gr_f_area = register_scalar_field('ice_shelf_model', 'tot_area_float_G', CS%diag%axesT1, CS%Time, &
      'Total area of Greenland floating ice shelves', 'm2', conversion=US%L_to_m**2)
  CS%id_Gr_g_area = register_scalar_field('ice_shelf_model', 'tot_area_ground_G', CS%diag%axesT1, CS%Time, &
      'Total area of Greenland grounded ice sheet', 'm2', conversion=US%L_to_m**2)
  !scalars (area integrated rates over the Greenland ice sheet)
  CS%id_Gr_dvafdt = register_scalar_field('ice_shelf_model', 'int_vafdot_G', CS%diag%axesT1, CS%Time, &
      'Area integrated rate of change in Greenland ice-sheet volume above floatation', &
      units='m3 s-1', conversion=US%Z_to_m*US%L_to_m**2*US%s_to_T)
  CS%id_Gr_adot = register_scalar_field('ice_shelf_model', 'int_adot_G', CS%diag%axesT1, CS%Time, &
      'Area integrated rate of change in Greenland ice-sheet thickness due to surface accum+melt', &
      units='m3 s-1', conversion=US%Z_to_m*US%L_to_m**2*US%s_to_T)
  CS%id_Gr_g_adot = register_scalar_field('ice_shelf_model', 'int_adot_ground_G', CS%diag%axesT1, CS%Time, &
      'Area integrated rate of change in Greenland grounded ice-sheet thickness due to surface accum+melt', &
      units='m3 s-1', conversion=US%Z_to_m*US%L_to_m**2*US%s_to_T)
  CS%id_Gr_f_adot = register_scalar_field('ice_shelf_model', 'int_adot_float_G', CS%diag%axesT1, CS%Time, &
      'Area integrated rate of change in Greenland floating ice-shelf thickness due to surface accum+melt', &
      units='m3 s-1', conversion=US%Z_to_m*US%L_to_m**2*US%s_to_T)
  CS%id_Gr_bdot = register_scalar_field('ice_shelf_model', 'int_bdot_G', CS%diag%axesT1, CS%Time, &
      'Area integrated rate of change in Greenland ice-shelf thickness due to basal accum+melt', &
      units='m3 s-1', conversion=US%Z_to_m*US%L_to_m**2*US%s_to_T)
  CS%id_Gr_bdot_melt = register_scalar_field('ice_shelf_model', 'int_bdot_melt_G', CS%diag%axesT1, CS%Time, &
      'Area integrated basal melt rate over Greenland ice shelves', &
      units='m3 s-1', conversion=US%Z_to_m*US%L_to_m**2*US%s_to_T)
  CS%id_Gr_bdot_accum = register_scalar_field('ice_shelf_model', 'int_bdot_accum_G', CS%diag%axesT1, CS%Time, &
      'Area integrated basal accumulation rate over Greenland ice shelves', &
      units='m3 s-1', conversion=US%Z_to_m*US%L_to_m**2*US%s_to_T)

  !Flags to calculate diagnostics related to surface/basal mass balance
    if (CS%id_adott>0     .or. CS%id_g_adott>0     .or. CS%id_f_adott>0     .or. &
        CS%id_adot >0     .or. CS%id_g_adot >0     .or. CS%id_f_adot >0     .or. &
        CS%id_Ant_adott>0 .or. CS%id_Ant_g_adott>0 .or. CS%id_Ant_f_adott>0 .or. &
        CS%id_Ant_adot >0 .or. CS%id_Ant_g_adot >0 .or. CS%id_Ant_f_adot >0 .or. &
        CS%id_Gr_adott>0  .or. CS%id_Gr_g_adott>0  .or. CS%id_Gr_f_adott>0  .or. &
        CS%id_Gr_adot >0  .or. CS%id_Gr_g_adot >0  .or. CS%id_Gr_f_adot >0) then
      CS%smb_diag=.true.
    else
      CS%smb_diag=.false.
    endif

    if (CS%id_bdott>0     .or. CS%id_bdott_melt>0     .or. CS%id_bdott_accum>0     .or. &
        CS%id_bdot >0     .or. CS%id_bdot_melt >0     .or. CS%id_bdot_accum >0     .or. &
        CS%id_Ant_bdott>0 .or. CS%id_Ant_bdott_melt>0 .or. CS%id_Ant_bdott_accum>0 .or. &
        CS%id_Ant_bdot >0 .or. CS%id_Ant_bdot_melt >0 .or. CS%id_Ant_bdot_accum >0 .or. &
        CS%id_Gr_bdott>0  .or. CS%id_Gr_bdott_melt>0  .or. CS%id_Gr_bdott_accum>0  .or. &
        CS%id_Gr_bdot >0  .or. CS%id_Gr_bdot_melt >0  .or. CS%id_Gr_bdot_accum >0) then
      CS%bmb_diag=.true.
    else
      CS%bmb_diag=.false.
    endif

  call MOM_IS_diag_mediator_close_registration(CS%diag)

  if (present(forces_in)) call initialize_ice_shelf_forces(CS, ocn_grid, US, forces_in)

end procedure initialize_ice_shelf
module procedure initialize_ice_shelf_fluxes
  type(ocean_grid_type), pointer :: G  => NULL() ! Pointers to grids for convenience.
  type(forcing), pointer :: fluxes =>  NULL()
  integer :: i, j, isd, ied, jsd, jed
  G => CS%Grid
  isd = G%isd ; jsd = G%jsd ; ied = G%ied ; jed = G%jed

  ! Allocate the arrays for passing ice-shelf data through the forcing type.
  if (.not. CS%solo_ice_sheet) then
      call MOM_mesg("MOM_ice_shelf.F90, initialize_ice_shelf: allocating fluxes.")
   ! GMM: the following assures that water/heat fluxes are just allocated
   ! when SHELF_THERMO = True. These fluxes are necessary if one wants to
   ! use either ENERGETICS_SFC_PBL (ALE mode) or BULKMIXEDLAYER (layer mode).
    call allocate_forcing_type(CS%Grid_in, fluxes_in, ustar=.true., shelf=.true., &
         press=.true., water=CS%isthermo, heat=CS%isthermo, shelf_sfc_accumulation=CS%active_shelf_dynamics, &
         tau_mag=.true.)
  else
    call MOM_mesg("MOM_ice_shelf.F90, initialize_ice_shelf: allocating fluxes in solo mode.")
    call allocate_forcing_type(CS%Grid_in, fluxes_in, ustar=.true., shelf=.true., &
         press=.true., shelf_sfc_accumulation=CS%active_shelf_dynamics, tau_mag=.true.)
  endif
  if (CS%rotate_index) then
    allocate(fluxes)
    call allocate_forcing_type(fluxes_in, CS%Grid, fluxes, turns=CS%turns)
    call rotate_forcing(fluxes_in, fluxes, CS%turns)
  else
    fluxes => fluxes_in
  endif

  do j=jsd,jed ; do i=isd,ied
    if (G%areaT(i,j)>0.) fluxes%frac_shelf_h(i,j) = CS%ISS%area_shelf_h(i,j) / G%areaT(i,j)
  enddo ; enddo
  if (CS%debug) call hchksum(fluxes%frac_shelf_h, "IS init: frac_shelf_h", G%HI, haloshift=0)
  call add_shelf_pressure(ocn_grid, US, CS, fluxes)

  if (CS%rotate_index) then
    call rotate_forcing(fluxes, fluxes_in, -CS%turns)
    call deallocate_forcing_type(fluxes)
    deallocate(fluxes)
  endif

end procedure initialize_ice_shelf_fluxes
module procedure initialize_ice_shelf_forces
  type(mech_forcing), pointer :: forces => NULL()
  call MOM_mesg("MOM_ice_shelf.F90, initialize_ice_shelf: allocating forces.")

  if ((Ocn_grid%isc /= CS%Grid_in%isc) .or. (Ocn_grid%iec /= CS%Grid_in%iec) .or. &
      (Ocn_grid%jsc /= CS%Grid_in%jsc) .or. (Ocn_grid%jec /= CS%Grid_in%jec)) &
    call MOM_error(FATAL,"initialize_ice_shelf_forces: Incompatible ocean and external ice shelf grids.")

  call allocate_mech_forcing(CS%Grid_in, forces_in, ustar=.true., shelf=.true., press=.true., tau_mag=.true.)
  if (CS%rotate_index) then
    allocate(forces)
    call allocate_mech_forcing(forces_in, CS%Grid, forces)
    call rotate_mech_forcing(forces_in, CS%turns, forces)
  else
    forces=>forces_in
  endif

  call add_shelf_forces(CS%grid, US, CS, forces, do_shelf_area=.not.CS%solo_ice_sheet, &
                        external_call=.false.)

  if (CS%rotate_index) then
    call rotate_mech_forcing(forces, -CS%turns, forces_in)
    call deallocate_mech_forcing(forces)
  endif

end procedure initialize_ice_shelf_forces
module procedure initialize_shelf_mass
  integer :: i, j, is, ie, js, je
  logical :: read_shelf_area, new_sim_2
  character(len=240) :: config, inputdir, shelf_file, filename
  character(len=120) :: shelf_mass_var  ! The name of shelf mass in the file.
  character(len=120) :: shelf_area_var ! The name of shelf area in the file.
  character(len=40)  :: mdl = "MOM_ice_shelf"
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec

  new_sim_2 = .true. ; if (present(new_sim)) new_sim_2 = new_sim

  call get_param(param_file, mdl, "ICE_SHELF_CONFIG", config, &
                 "A string that specifies how the ice shelf is "//&
                 "initialized. Valid options include:\n"//&
                 " \tfile\t Read from a file.\n"//&
                 " \tzero\t Set shelf mass to 0 everywhere.\n"//&
                 " \tUSER\t Call USER_initialize_shelf_mass.\n", &
                 fail_if_missing=.true.)

  select case ( trim(config) )
    case ("file")

      call time_interp_external_init()

      call get_param(param_file, mdl, "INPUTDIR", inputdir, default=".")
      inputdir = slasher(inputdir)

      call get_param(param_file, mdl, "SHELF_FILE", shelf_file, &
              "If DYNAMIC_SHELF_MASS = True, OVERRIDE_SHELF_MOVEMENT = True "//&
              "and ICE_SHELF_MASS_FROM_FILE = True, this is the file from "//&
              "which to read the shelf mass and area.", &
               default="shelf_mass.nc")
      call get_param(param_file, mdl, "SHELF_MASS_VAR", shelf_mass_var, &
                 "The variable in SHELF_FILE with the shelf mass.", &
                 default="shelf_mass")
      call get_param(param_file, mdl, "READ_SHELF_AREA", read_shelf_area, &
                 "If true, also read the area covered by ice-shelf from SHELF_FILE.", &
                 default=.false.)

      filename = trim(slasher(inputdir))//trim(shelf_file)
      call log_param(param_file, mdl, "INPUTDIR/SHELF_FILE", filename)

      CS%mass_handle = init_external_field(filename, shelf_mass_var, &
                            MOM_domain=CS%Grid_in%Domain, verbose=CS%debug)

      if (read_shelf_area) then
         call get_param(param_file, mdl, "SHELF_AREA_VAR", shelf_area_var, &
                  "The variable in SHELF_FILE with the shelf area.", &
                  default="shelf_area")

         CS%area_handle = init_external_field(filename, shelf_area_var, &
                               MOM_domain=CS%Grid_in%Domain)
      endif

      if (.not.file_exists(filename, CS%Grid_in%Domain)) call MOM_error(FATAL, &
           " initialize_shelf_mass: Unable to open "//trim(filename))

    case ("zero")
      do j=js,je ; do i=is,ie
        ISS%mass_shelf(i,j) = 0.0
        ISS%area_shelf_h(i,j) = 0.0
      enddo ; enddo

    case ("USER")
      call USER_initialize_shelf_mass(ISS%mass_shelf, ISS%area_shelf_h, &
                   ISS%h_shelf, ISS%hmask, G, CS%US, CS%user_CS, param_file, new_sim_2)

    case default ;  call MOM_error(FATAL,"initialize_ice_shelf: "// &
      "Unrecognized ice shelf setup "//trim(config))
  end select

end procedure initialize_shelf_mass
module procedure change_thickness_using_precip
  integer :: i, j
  real :: I_rho_ice ! The specific volume of ice [R-1 ~> m3 kg-1]
  I_rho_ice = 1.0 / CS%density_ice

  !update time
!  CS%Time = Time

!    CS%time_step = time_step
    ! update surface mass flux  rate
!    if (CS%surf_mass_flux_from_file) call update_surf_mass_flux(G, US, CS, ISS, Time)

  do j=G%jsc,G%jec ; do i=G%isc,G%iec
    if ((ISS%hmask(i,j) == 1) .or. (ISS%hmask(i,j) == 2)) then

      if (-fluxes%shelf_sfc_mass_flux(i,j) * time_step * I_rho_ice  < ISS%h_shelf(i,j)) then
        ISS%h_shelf(i,j) = ISS%h_shelf(i,j) + fluxes%shelf_sfc_mass_flux(i,j) * time_step * I_rho_ice
      else
        ! the ice is about to ablate, so set thickness, area, and mask to zero
        ! NOTE: this is not mass conservative should maybe scale salt & heat flux for this cell
        ISS%h_shelf(i,j) = 0.0
        ISS%hmask(i,j) = 0.0
        ISS%area_shelf_h(i,j) = 0.0
      endif
      ISS%mass_shelf(i,j) = ISS%h_shelf(i,j) * CS%density_ice
    endif
  enddo ; enddo

  call pass_var(ISS%area_shelf_h, G%domain, complete=.false.)
  call pass_var(ISS%h_shelf, G%domain, complete=.false.)
  call pass_var(ISS%hmask, G%domain, complete=.false.)
  call pass_var(ISS%mass_shelf, G%domain, complete=.true.)

end procedure change_thickness_using_precip
module procedure update_shelf_mass
  integer :: i, j, is, ie, js, je
  real, allocatable, dimension(:,:) :: tmp2d ! Temporary array for storing ice shelf input data [R Z ~> kg m-2]
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec


  if (CS%rotate_index) then
    allocate(tmp2d(CS%Grid_in%isc:CS%Grid_in%iec,CS%Grid_in%jsc:CS%Grid_in%jec), source=0.0)
  else
    allocate(tmp2d(is:ie,js:je), source=0.0)
  endif

  call time_interp_external(CS%mass_handle, Time, tmp2d, scale=US%kg_m3_to_R*US%m_to_Z)
  call rotate_array(tmp2d, CS%turns, ISS%mass_shelf)
  deallocate(tmp2d)

  do j=js,je ; do i=is,ie
    ISS%area_shelf_h(i,j) = 0.0
    ISS%hmask(i,j) = 0.
    if (ISS%mass_shelf(i,j) > 0.0) then
      ISS%area_shelf_h(i,j) = G%areaT(i,j)
      ISS%h_shelf(i,j) = ISS%mass_shelf(i,j) / CS%density_ice
      ISS%hmask(i,j) = 1.
    endif
  enddo ; enddo

  !call USER_update_shelf_mass(ISS%mass_shelf, ISS%area_shelf_h, ISS%h_shelf, &
  !                            ISS%hmask, CS%grid, CS%user_CS, Time, .true.)

  if (CS%min_thickness_simple_calve > 0.0) then
    call ice_shelf_min_thickness_calve(G, ISS%h_shelf, ISS%area_shelf_h, ISS%hmask, &
                                       CS%min_thickness_simple_calve, halo=0)
  endif

  call pass_var(ISS%area_shelf_h, G%domain, complete=.false.)
  call pass_var(ISS%h_shelf, G%domain, complete=.false.)
  call pass_var(ISS%hmask, G%domain, complete=.false.)
  call pass_var(ISS%mass_shelf, G%domain, complete=.true.)

end procedure update_shelf_mass
module procedure ice_shelf_query
  integer :: i, j
  if (present(frac_shelf_h)) then
    do j=G%jsd,G%jed ; do i=G%isd,G%ied
      frac_shelf_h(i,j) = 0.0
      if (G%areaT(i,j)>0.) frac_shelf_h(i,j) = CS%ISS%area_shelf_h(i,j) / G%areaT(i,j)
    enddo ; enddo
  endif

  if (present(mass_shelf)) then
    do j=G%jsd,G%jed ; do i=G%isd,G%ied
      mass_shelf(i,j) = 0.0
      if (G%areaT(i,j)>0.) mass_shelf(i,j) = CS%ISS%mass_shelf(i,j)
    enddo ; enddo
  endif

  if (present(data_override_shelf_fluxes)) then
    data_override_shelf_fluxes=.false.
    if (CS%active_shelf_dynamics) data_override_shelf_fluxes = CS%data_override_shelf_fluxes
  endif

end procedure ice_shelf_query
module procedure ice_shelf_save_restart
  type(ocean_grid_type), pointer :: G => NULL()
  character(len=200) :: restart_dir
  G => CS%grid

  if (present(directory)) then ; restart_dir = directory
  else ; restart_dir = CS%restart_output_dir ; endif

  call save_restart(restart_dir, Time, CS%grid_in, CS%restart_CSp, time_stamped)

end procedure ice_shelf_save_restart
module procedure ice_shelf_end
  if (.not.associated(CS)) return

  call ice_shelf_state_end(CS%ISS)

  if (CS%active_shelf_dynamics) call ice_shelf_dyn_end(CS%dCS)

  call MOM_IS_diag_mediator_end(CS%diag)
  deallocate(CS)

end procedure ice_shelf_end
module procedure solo_step_ice_shelf
  type(ocean_grid_type), pointer :: G => NULL()  ! A pointer to the ocean's grid structure
  type(unit_scale_type), pointer :: US => NULL() ! Pointer to a structure containing
  type(ice_shelf_state), pointer :: ISS => NULL() !< A structure with elements that describe
  real :: remaining_time    ! The remaining time in this call [T ~> s]
  real :: time_step         ! The internal time step during this call [T ~> s]
  real :: full_time_step    ! The external time step (sum of internal time steps) during this call [T ~> s]
  real :: Ifull_time_step   ! The inverse of the external time step [T-1 ~> s-1]
  real :: min_time_step     ! The minimal required timestep that would indicate a fatal problem [T ~> s]
  character(len=240) :: mesg
  logical :: update_ice_vel ! If true, it is time to update the ice shelf velocities.
  logical :: coupled_GL     ! If true the grounding line position is determined based on
  integer :: is, ie, js, je, i, j
  real :: vaf0, vaf0_A, vaf0_G !The previous volumes above floatation
  real, dimension(SZI_(CS%grid),SZJ_(CS%grid)) :: &
    dh_adott_sum, &    ! Surface melt/accumulation over a full time step, used for diagnostics [Z ~> m]
    dh_adott           ! Surface melt/accumulation over a partial time step, used for diagnostics [Z ~> m]
  G => CS%grid
  US => CS%US
  ISS => CS%ISS
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec

  remaining_time = time_to_real(time_interval, scale=US%s_to_T)
  full_time_step = remaining_time
  Ifull_time_step = 1./full_time_step

  if (present (min_time_step_in)) then
    min_time_step = min_time_step_in
  else
    min_time_step = 1000.0*US%s_to_T ! At 1 km resolution this would imply ice is moving at ~1 meter per second
  endif

  write (mesg,*) "TIME in ice shelf call, yrs: ", time_to_real(Time)/(365. * 86400.)
  call MOM_mesg("solo_step_ice_shelf: "//mesg, 5)

  ISS%dhdt_shelf(:,:) = ISS%h_shelf(:,:)

  dh_adott(:,:)=0.0

  if (CS%smb_diag) dh_adott_sum(:,:) = 0.0

  !calculate previous volumes above floatation
  if (CS%id_dvafdt     > 0) call volume_above_floatation(CS%dCS, G, ISS, vaf0)                 !all ice sheet
  if (CS%id_Ant_dvafdt > 0) call volume_above_floatation(CS%dCS, G, ISS, vaf0_A, hemisphere=0) !Antarctica only
  if (CS%id_Gr_dvafdt  > 0) call volume_above_floatation(CS%dCS, G, ISS, vaf0_G, hemisphere=1) !Greenland only

  do while (remaining_time > 0.0)
    nsteps = nsteps+1

    ! If time_interval is not too long, this is unnecessary.
    time_step = min(ice_time_step_CFL(CS%dCS, ISS, G), remaining_time)

    write (mesg,*) "Ice model timestep = ", US%T_to_s*time_step, " seconds"
    if ((time_step < min_time_step) .and. (time_step < remaining_time))  then
      call MOM_error(FATAL, "MOM_ice_shelf:solo_step_ice_shelf: abnormally small timestep "//mesg)
    else
      call MOM_mesg("solo_step_ice_shelf: "//mesg, 5)
    endif

    if (CS%smb_diag) dh_adott(is:ie,js:je) = ISS%h_shelf(is:ie,js:je)
    call change_thickness_using_precip(CS, ISS, G, US, fluxes_in, time_step, Time)
    if (CS%smb_diag) dh_adott_sum(is:ie,js:je) = dh_adott_sum(is:ie,js:je) + &
                                             (ISS%h_shelf(is:ie,js:je) - dh_adott(is:ie,js:je))

    remaining_time = remaining_time - time_step

    ! If the last mini-timestep is a day or less, we cannot expect velocities to change by much.
    ! Do not update the velocities if the last step is very short.
    update_ice_vel = ((time_step > min_time_step) .or. (remaining_time > 0.0))
    coupled_GL = .false.

    call update_ice_shelf(CS%dCS, ISS, G, US, time_step, Time, CS%calve_ice_shelf_bergs, &
                          must_update_vel=update_ice_vel)

  enddo

  call write_ice_shelf_energy(CS%dCS, G, US, ISS%mass_shelf, ISS%area_shelf_h, Time, &
                              time_step=time_interval)
  do j=js,je ; do i=is,ie
    ISS%dhdt_shelf(i,j) = (ISS%h_shelf(i,j) - ISS%dhdt_shelf(i,j)) * Ifull_time_step
  enddo ; enddo

  call enable_averages(full_time_step, Time, CS%diag)
  if (CS%id_area_shelf_h > 0) call post_data(CS%id_area_shelf_h ,ISS%area_shelf_h,CS%diag)
  if (CS%id_h_shelf > 0)      call post_data(CS%id_h_shelf      ,ISS%h_shelf     ,CS%diag)
  if (CS%id_dhdt_shelf > 0)   call post_data(CS%id_dhdt_shelf   ,ISS%dhdt_shelf  ,CS%diag)
  if (CS%id_h_mask > 0)       call post_data(CS%id_h_mask       ,ISS%hmask       ,CS%diag)
  call process_and_post_scalar_data(CS, vaf0, vaf0_A, vaf0_G, Ifull_time_step, dh_adott, dh_adott*0.0)
  call disable_averaging(CS%diag)

  call IS_dynamics_post_data(full_time_step, Time, CS%dCS, ISS, G)
end procedure solo_step_ice_shelf
module procedure process_and_post_scalar_data
  real, dimension(SZI_(CS%grid),SZJ_(CS%grid)) :: tmp ! Temporary field used when calculating diagnostics [various]
  real, dimension(SZI_(CS%grid),SZJ_(CS%grid)) :: ones ! Temporary field used when calculating diagnostics [various]
  real :: vaf   ! The current ice-sheet volume above floatation [Z L2 ~> m3]
  real :: val   ! Temporary value when calculating scalar diagnostics [various]
  type(ocean_grid_type), pointer :: G => NULL()  ! A pointer to the ocean's grid structure
  type(unit_scale_type), pointer :: US => NULL() ! Pointer to a structure containing various unit conversion factors
  type(ice_shelf_state), pointer :: ISS => NULL() ! A structure with elements that describe the ice-shelf state
  integer :: is, ie, js, je, i, j
  G => CS%grid
  US => CS%US
  ISS => CS%ISS
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec

  !---ALL ICE SHEET---!
  if (CS%id_vaf > 0 .or. CS%id_dvafdt > 0) &  !calculate current volume above floatation (vaf)
    call volume_above_floatation(CS%dCS, G, ISS, vaf)
  if (CS%id_vaf    > 0) call post_scalar_data(CS%id_vaf   ,vaf                  ,CS%diag) !current vaf
  if (CS%id_dvafdt > 0) call post_scalar_data(CS%id_dvafdt,(vaf-vaf0)*Itime_step,CS%diag) !d(vaf)/dt
  if (CS%id_adott > 0 .or. CS%id_adot > 0) then !surface accumulation - surface melt
    val = integrate_over_ice_sheet_area(G, ISS, dh_adott, unscale=US%Z_to_m)
    if (CS%id_adott > 0) call post_scalar_data(CS%id_adott,val           ,CS%diag)
    if (CS%id_adot  > 0) call post_scalar_data(CS%id_adot ,val*Itime_step,CS%diag)
  endif
  if (CS%id_g_adott > 0 .or. CS%id_g_adot > 0) then !grounded only: surface accumulation - surface melt
    call masked_var_grounded(G,CS%dCS,dh_adott,tmp)
    val = integrate_over_ice_sheet_area(G, ISS, tmp, unscale=US%Z_to_m)
    if (CS%id_g_adott > 0) call post_scalar_data(CS%id_g_adott,val           ,CS%diag)
    if (CS%id_g_adot  > 0) call post_scalar_data(CS%id_g_adot ,val*Itime_step,CS%diag)
  endif
  if (CS%id_f_adott > 0 .or. CS%id_f_adot > 0) then !floating only: surface accumulation - surface melt
    call masked_var_grounded(G,CS%dCS,dh_adott,tmp)
    do j=js,je ; do i=is,ie
      tmp(i,j) = dh_adott(i,j) - tmp(i,j)
    enddo ; enddo
    val = integrate_over_ice_sheet_area(G, ISS, tmp, unscale=US%Z_to_m)
    if (CS%id_f_adott > 0) call post_scalar_data(CS%id_f_adott,val           ,CS%diag)
    if (CS%id_f_adot  > 0) call post_scalar_data(CS%id_f_adot ,val*Itime_step,CS%diag)
  endif
  if (CS%id_bdott > 0 .or. CS%id_bdot > 0) then !bottom accumulation - bottom melt
    val = integrate_over_ice_sheet_area(G, ISS, dh_bdott, unscale=US%Z_to_m)
    if (CS%id_bdott > 0) call post_scalar_data(CS%id_bdott,val           ,CS%diag)
    if (CS%id_bdot  > 0) call post_scalar_data(CS%id_bdot ,val*Itime_step,CS%diag)
  endif
  if (CS%id_bdott_melt > 0 .or. CS%id_bdot_melt > 0) then !bottom melt
    tmp(:,:)=0.0
    do j=js,je ; do i=is,ie
      if (dh_bdott(i,j) < 0) tmp(i,j) = -dh_bdott(i,j)
    enddo ; enddo
    val = integrate_over_ice_sheet_area(G, ISS, tmp, unscale=US%Z_to_m)
    if (CS%id_bdott_melt > 0) call post_scalar_data(CS%id_bdott_melt,val           ,CS%diag)
    if (CS%id_bdot_melt  > 0) call post_scalar_data(CS%id_bdot_melt ,val*Itime_step,CS%diag)
  endif
  if (CS%id_bdott_accum > 0 .or. CS%id_bdot_accum > 0) then !bottom accumulation
    tmp(:,:)=0.0
    do j=js,je ; do i=is,ie
      if (dh_bdott(i,j) > 0) tmp(i,j) = dh_bdott(i,j)
    enddo ; enddo
    val = integrate_over_ice_sheet_area(G, ISS, tmp, unscale=US%Z_to_m)
    if (CS%id_bdott_accum > 0) call post_scalar_data(CS%id_bdott_accum,val           ,CS%diag)
    if (CS%id_bdot_accum  > 0) call post_scalar_data(CS%id_bdot_accum ,val*Itime_step,CS%diag)
  endif
  if (CS%id_t_area > 0) then !ice sheet area
    tmp(:,:) = 1.0 ; val = integrate_over_ice_sheet_area(G, ISS, tmp, unscale=1.0)
    call post_scalar_data(CS%id_t_area,val,CS%diag)
  endif
  if (CS%id_g_area > 0 .or. CS%id_f_area > 0) then
    ones(:,:) = 1.0 ; call masked_var_grounded(G, CS%dCS, ones, tmp)
    if (CS%id_g_area > 0) then !grounded only ice sheet area
      val = integrate_over_ice_sheet_area(G, ISS,     tmp, unscale=1.0)
      call post_scalar_data(CS%id_g_area,val,CS%diag)
    endif
    if (CS%id_f_area > 0) then !floating only ice sheet area (ice shelf area)
      val = integrate_over_ice_sheet_area(G, ISS, 1.0-tmp, unscale=1.0)
      call post_scalar_data(CS%id_f_area,val,CS%diag)
    endif
  endif

  !---ANTARCTICA ONLY---!
  if (CS%id_Ant_vaf > 0 .or. CS%id_Ant_dvafdt > 0) &  !calculate current volume above floatation (vaf)
    call volume_above_floatation(CS%dCS, G, ISS, vaf, hemisphere=0)
  if (CS%id_Ant_vaf    > 0) call post_scalar_data(CS%id_Ant_vaf   ,vaf                  ,CS%diag) !current vaf
  if (CS%id_Ant_dvafdt > 0) call post_scalar_data(CS%id_Ant_dvafdt,(vaf-vaf0_A)*Itime_step,CS%diag) !d(vaf)/dt
  if (CS%id_Ant_adott > 0 .or. CS%id_Ant_adot > 0) then !surface accumulation - surface melt
    val = integrate_over_ice_sheet_area(G, ISS, dh_adott, unscale=US%Z_to_m, hemisphere=0)
    if (CS%id_Ant_adott > 0) call post_scalar_data(CS%id_Ant_adott,val           ,CS%diag)
    if (CS%id_Ant_adot  > 0) call post_scalar_data(CS%id_Ant_adot ,val*Itime_step,CS%diag)
  endif
  if (CS%id_Ant_g_adott > 0 .or. CS%id_Ant_g_adot > 0) then !grounded only: surface accumulation - surface melt
    call masked_var_grounded(G,CS%dCS,dh_adott,tmp)
    val = integrate_over_ice_sheet_area(G, ISS, tmp, unscale=US%Z_to_m, hemisphere=0)
    if (CS%id_Ant_g_adott > 0) call post_scalar_data(CS%id_Ant_g_adott,val           ,CS%diag)
    if (CS%id_Ant_g_adot  > 0) call post_scalar_data(CS%id_Ant_g_adot ,val*Itime_step,CS%diag)
  endif
  if (CS%id_Ant_f_adott > 0 .or. CS%id_Ant_f_adot > 0) then !floating only: surface accumulation - surface melt
    call masked_var_grounded(G,CS%dCS,dh_adott,tmp)
    do j=js,je ; do i=is,ie
      tmp(i,j) = dh_adott(i,j) - tmp(i,j)
    enddo ; enddo
    val = integrate_over_ice_sheet_area(G, ISS, tmp, unscale=US%Z_to_m, hemisphere=0)
    if (CS%id_Ant_f_adott > 0) call post_scalar_data(CS%id_Ant_f_adott,val           ,CS%diag)
    if (CS%id_Ant_f_adot  > 0) call post_scalar_data(CS%id_Ant_f_adot ,val*Itime_step,CS%diag)
  endif
  if (CS%id_Ant_bdott > 0 .or. CS%id_Ant_bdot > 0) then !bottom accumulation - bottom melt
    val = integrate_over_ice_sheet_area(G, ISS, dh_bdott, unscale=US%Z_to_m, hemisphere=0)
    if (CS%id_Ant_bdott > 0) call post_scalar_data(CS%id_Ant_bdott,val           ,CS%diag)
    if (CS%id_Ant_bdot  > 0) call post_scalar_data(CS%id_Ant_bdot ,val*Itime_step,CS%diag)
  endif
  if (CS%id_Ant_bdott_melt > 0 .or. CS%id_Ant_bdot_melt > 0) then !bottom melt
    tmp(:,:)=0.0
    do j=js,je ; do i=is,ie
      if (dh_bdott(i,j) < 0) tmp(i,j) = -dh_bdott(i,j)
    enddo ; enddo
    val = integrate_over_ice_sheet_area(G, ISS, tmp, unscale=US%Z_to_m, hemisphere=0)
    if (CS%id_Ant_bdott_melt > 0) call post_scalar_data(CS%id_Ant_bdott_melt,val           ,CS%diag)
    if (CS%id_Ant_bdot_melt  > 0) call post_scalar_data(CS%id_Ant_bdot_melt ,val*Itime_step,CS%diag)
  endif
  if (CS%id_Ant_bdott_accum > 0 .or. CS%id_Ant_bdot_accum > 0) then !bottom accumulation
    tmp(:,:)=0.0
    do j=js,je ; do i=is,ie
      if (dh_bdott(i,j) > 0) tmp(i,j) = dh_bdott(i,j)
    enddo ; enddo
    val = integrate_over_ice_sheet_area(G, ISS, tmp, unscale=US%Z_to_m, hemisphere=0)
    if (CS%id_Ant_bdott_accum > 0) call post_scalar_data(CS%id_Ant_bdott_accum,val           ,CS%diag)
    if (CS%id_Ant_bdot_accum  > 0) call post_scalar_data(CS%id_Ant_bdot_accum ,val*Itime_step,CS%diag)
  endif
  if (CS%id_Ant_t_area > 0) then !ice sheet area
    tmp(:,:) = 1.0 ; val = integrate_over_ice_sheet_area(G, ISS, tmp, unscale=1.0, hemisphere=0)
    call post_scalar_data(CS%id_Ant_t_area,val,CS%diag)
  endif
  if (CS%id_Ant_g_area > 0 .or. CS%id_Ant_f_area > 0) then
    ones(:,:) = 1.0 ; call masked_var_grounded(G, CS%dCS, ones, tmp)
    if (CS%id_Ant_g_area > 0) then !grounded only ice sheet area
      val = integrate_over_ice_sheet_area(G, ISS,     tmp, unscale=1.0, hemisphere=0)
      call post_scalar_data(CS%id_Ant_g_area,val,CS%diag)
    endif
    if (CS%id_Ant_f_area > 0) then !floating only ice sheet area (ice shelf area)
      val = integrate_over_ice_sheet_area(G, ISS, 1.0-tmp, unscale=1.0, hemisphere=0)
      call post_scalar_data(CS%id_Ant_f_area,val,CS%diag)
    endif
  endif

  !---GREENLAND ONLY---!
  if (CS%id_Gr_vaf > 0 .or. CS%id_Gr_dvafdt > 0) &  !calculate current volume above floatation (vaf)
    call volume_above_floatation(CS%dCS, G, ISS, vaf, hemisphere=1)
  if (CS%id_Gr_vaf    > 0) call post_scalar_data(CS%id_Gr_vaf   ,vaf                  ,CS%diag) !current vaf
  if (CS%id_Gr_dvafdt > 0) call post_scalar_data(CS%id_Gr_dvafdt,(vaf-vaf0_A)*Itime_step,CS%diag) !d(vaf)/dt
  if (CS%id_Gr_adott > 0 .or. CS%id_Gr_adot > 0) then !surface accumulation - surface melt
    val = integrate_over_ice_sheet_area(G, ISS, dh_adott, unscale=US%Z_to_m, hemisphere=1)
    if (CS%id_Gr_adott > 0) call post_scalar_data(CS%id_Gr_adott,val           ,CS%diag)
    if (CS%id_Gr_adot  > 0) call post_scalar_data(CS%id_Gr_adot ,val*Itime_step,CS%diag)
  endif
  if (CS%id_Gr_g_adott > 0 .or. CS%id_Gr_g_adot > 0) then !grounded only: surface accumulation - surface melt
    call masked_var_grounded(G,CS%dCS,dh_adott,tmp)
    val = integrate_over_ice_sheet_area(G, ISS, tmp, unscale=US%Z_to_m, hemisphere=1)
    if (CS%id_Gr_g_adott > 0) call post_scalar_data(CS%id_Gr_g_adott,val           ,CS%diag)
    if (CS%id_Gr_g_adot  > 0) call post_scalar_data(CS%id_Gr_g_adot ,val*Itime_step,CS%diag)
  endif
  if (CS%id_Gr_f_adott > 0 .or. CS%id_Gr_f_adot > 0) then !floating only: surface accumulation - surface melt
    call masked_var_grounded(G,CS%dCS,dh_adott,tmp)
    do j=js,je ; do i=is,ie
      tmp(i,j) = dh_adott(i,j) - tmp(i,j)
    enddo ; enddo
    val = integrate_over_ice_sheet_area(G, ISS, tmp, unscale=US%Z_to_m, hemisphere=1)
    if (CS%id_Gr_f_adott > 0) call post_scalar_data(CS%id_Gr_f_adott,val           ,CS%diag)
    if (CS%id_Gr_f_adot  > 0) call post_scalar_data(CS%id_Gr_f_adot ,val*Itime_step,CS%diag)
  endif
  if (CS%id_Gr_bdott > 0 .or. CS%id_Gr_bdot > 0) then !bottom accumulation - bottom melt
    val = integrate_over_ice_sheet_area(G, ISS, dh_bdott, unscale=US%Z_to_m, hemisphere=1)
    if (CS%id_Gr_bdott > 0) call post_scalar_data(CS%id_Gr_bdott,val           ,CS%diag)
    if (CS%id_Gr_bdot  > 0) call post_scalar_data(CS%id_Gr_bdot ,val*Itime_step,CS%diag)
  endif
  if (CS%id_Gr_bdott_melt > 0 .or. CS%id_Gr_bdot_melt > 0) then !bottom melt
    tmp(:,:)=0.0
    do j=js,je ; do i=is,ie
      if (dh_bdott(i,j) < 0) tmp(i,j) = -dh_bdott(i,j)
    enddo ; enddo
    val = integrate_over_ice_sheet_area(G, ISS, tmp, unscale=US%Z_to_m, hemisphere=1)
    if (CS%id_Gr_bdott_melt > 0) call post_scalar_data(CS%id_Gr_bdott_melt,val           ,CS%diag)
    if (CS%id_Gr_bdot_melt  > 0) call post_scalar_data(CS%id_Gr_bdot_melt ,val*Itime_step,CS%diag)
  endif
  if (CS%id_Gr_bdott_accum > 0 .or. CS%id_Gr_bdot_accum > 0) then !bottom accumulation
    tmp(:,:)=0.0
    do j=js,je ; do i=is,ie
      if (dh_bdott(i,j) > 0) tmp(i,j) = dh_bdott(i,j)
    enddo ; enddo
    val = integrate_over_ice_sheet_area(G, ISS, tmp, unscale=US%Z_to_m, hemisphere=1)
    if (CS%id_Gr_bdott_accum > 0) call post_scalar_data(CS%id_Gr_bdott_accum,val           ,CS%diag)
    if (CS%id_Gr_bdot_accum  > 0) call post_scalar_data(CS%id_Gr_bdot_accum ,val*Itime_step,CS%diag)
  endif
  if (CS%id_Gr_t_area > 0) then !ice sheet area
    tmp(:,:) = 1.0 ; val = integrate_over_ice_sheet_area(G, ISS, tmp, unscale=1.0, hemisphere=1)
    call post_scalar_data(CS%id_Gr_t_area,val,CS%diag)
  endif
  if (CS%id_Gr_g_area > 0 .or. CS%id_Gr_f_area > 0) then
    ones(:,:) = 1.0 ; call masked_var_grounded(G, CS%dCS, ones, tmp)
    if (CS%id_Gr_g_area > 0) then !grounded only ice sheet area
      val = integrate_over_ice_sheet_area(G, ISS,     tmp, unscale=1.0, hemisphere=1)
      call post_scalar_data(CS%id_Gr_g_area,val,CS%diag)
    endif
    if (CS%id_Gr_f_area > 0) then !floating only ice sheet area (ice shelf area)
      val = integrate_over_ice_sheet_area(G, ISS, 1.0-tmp, unscale=1.0, hemisphere=1)
      call post_scalar_data(CS%id_Gr_f_area,val,CS%diag)
    endif
  endif
end procedure process_and_post_scalar_data
end submodule MOM_ice_shelf_s
