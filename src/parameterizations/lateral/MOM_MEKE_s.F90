submodule (MOM_MEKE) MOM_MEKE_s
#include <MOM_memory.h>
  implicit none
contains
module procedure step_forward_MEKE
  real, dimension(SZI_(G),SZJ_(G)) :: &
    data_eke, &     ! EKE from file [L2 T-2 ~> m2 s-2]
    mass, &         ! The total mass of the water column [R Z ~> kg m-2].
    I_mass, &       ! The inverse of mass [R-1 Z-1 ~> m2 kg-1].
    depth_tot, &    ! The depth of the water column [H ~> m or kg m-2].
    src, &          ! The sum of all MEKE sources [L2 T-3 ~> W kg-1] (= m2 s-3).
    MEKE_decay, &   ! A diagnostic of the MEKE decay timescale [T-1 ~> s-1].
    src_adv, &      ! The MEKE source/tendency from the horizontal advection of MEKE [L2 T-3 ~> W kg-1] (= m2 s-3).
    src_mom_K4, &   ! The MEKE source/tendency from the bihamornic of MEKE [L2 T-3 ~> W kg-1] (= m2 s-3).
    src_btm_drag, & ! The MEKE source/tendency from the bottom drag acting on MEKE [L2 T-3 ~> W kg-1] (= m2 s-3).
    src_GM, &       ! The MEKE source/tendency from the thickness mixing (GM) [L2 T-3 ~> W kg-1] (= m2 s-3).
    src_mom_lp, &   ! The MEKE source/tendency from the Laplacian of the resolved flow [L2 T-3 ~> W kg-1] (= m2 s-3).
    src_mom_bh, &   ! The MEKE source/tendency from the biharmonic of the resolved flow [L2 T-3 ~> W kg-1] (= m2 s-3).
    damp_rate_s1, & ! The MEKE damping rate computed at the 1st Strang splitting stage [T-1 ~> s-1].
    MEKE_current, & ! A copy of MEKE for use in computing the MEKE damping [L2 T-2 ~> m2 s-2].
    drag_rate_visc, & ! Near-bottom velocity contribution to bottom drag [H T-1 ~> m s-1 or kg m-2 s-1]
    drag_rate, &    ! The MEKE spindown timescale due to bottom drag [T-1 ~> s-1].
    del2MEKE, &     ! Laplacian of MEKE, used for bi-harmonic diffusion [T-2 ~> s-2].
    del4MEKE, &     ! Time-integrated MEKE tendency arising from the biharmonic of MEKE [L2 T-2 ~> m2 s-2].
    LmixScale, &    ! Eddy mixing length [L ~> m].
    barotrFac2, &   ! Ratio of EKE_barotropic / EKE [nondim]
    bottomFac2, &   ! Ratio of EKE_bottom / EKE [nondim]
    tmp, &          ! Temporary variable for computation of diagnostic velocities [L T-1 ~> m s-1]
    equilibrium_value, & ! The equilibrium value of MEKE to be calculated at
                    ! each time step [L2 T-2 ~> m2 s-2]
    damp_rate, &    ! The MEKE damping rate [T-1 ~> s-1]
    damping         ! The net damping of a field after sdt_damp [nondim]

  real, dimension(SZIB_(G),SZJ_(G)) :: &
    MEKE_uflux, &   ! The zonal advective and diffusive flux of MEKE with units of [R Z L4 T-3 ~> kg m2 s-3].
                    ! In one place, MEKE_uflux is used as temporary work space with units of [L2 T-2 ~> m2 s-2].
    Kh_u, &         ! The zonal diffusivity that is actually used [L2 T-1 ~> m2 s-1].
    baroHu, &       ! Depth integrated accumulated zonal mass flux [R Z L2 ~> kg].
    drag_vel_u      ! A piston velocity associated with bottom drag at u-points [H T-1 ~> m s-1 or kg m-2 s-1]
  real, dimension(SZI_(G),SZJB_(G)) :: &
    MEKE_vflux, &   ! The meridional advective and diffusive flux of MEKE with units of [R Z L4 T-3 ~> kg m2 s-3].
                    ! In one place, MEKE_vflux is used as temporary work space with units of [L2 T-2 ~> m2 s-2].
    Kh_v, &         ! The meridional diffusivity that is actually used [L2 T-1 ~> m2 s-1].
    baroHv, &       ! Depth integrated accumulated meridional mass flux [R Z L2 ~> kg].
    drag_vel_v      ! A piston velocity associated with bottom drag at v-points [H T-1 ~> m s-1 or kg m-2 s-1]
  real :: bh_coeff  ! Biharmonic part of efficiency conversion in total MEKE [nondim]
  real :: Kh_here   ! The local horizontal viscosity [L2 T-1 ~> m2 s-1]
  real :: Inv_Kh_max ! The inverse of the local horizontal viscosity [T L-2 ~> s m-2]
  real :: K4_here   ! The local horizontal biharmonic viscosity [L4 T-1 ~> m4 s-1]
  real :: Inv_K4_max ! The inverse of the local horizontal biharmonic viscosity [T L-4 ~> s m-4]
  real :: cdrag2    ! The square of the drag coefficient times unit conversion factors [H2 L-2 ~> nondim or kg2 m-6]
  real :: advFac    ! The product of the advection scaling factor and 1/dt [T-1 ~> s-1]
  real :: mass_neglect ! A negligible mass [R Z ~> kg m-2].
  real :: sdt       ! dt to use locally [T ~> s] (could be scaled to accelerate)
  real :: sdt_damp  ! dt for damping [T ~> s] (sdt could be split).
  real :: damp_step ! Size of damping timestep relative to sdt [nondim]
  logical :: use_drag_rate ! Flag to indicate drag_rate is finite
  logical :: any_damping_diags_s1 ! True if any damped diagnostics are enabled in first stage
  logical :: any_damping_diags  ! True if any damped diagnostics are enabled
  integer :: i, j, k, is, ie, js, je, Isq, Ieq, Jsq, Jeq, nz
  real(kind=real32), dimension(size(MEKE%MEKE),NUM_FEATURES) :: features_array ! The array of features
                                        ! needed for the machine learning inference, with different
                                        ! units for the various subarrays [various]

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke
  Isq = G%IscB ; Ieq = G%IecB ; Jsq = G%JscB ; Jeq = G%JecB

  if (.not.CS%initialized) call MOM_error(FATAL, &
         "MOM_MEKE: Module must be initialized before it is used.")

  if ((CS%MEKE_Cd_scale > 0.0) .or. (CS%MEKE_Cb>0.) .or. CS%visc_drag) then
    use_drag_rate = .true.
  else
    use_drag_rate = .false.
  endif

  ! Only integrate the MEKE equations if MEKE is required.
  if (.not. allocated(MEKE%MEKE)) then
!   call MOM_error(FATAL, "MOM_MEKE: MEKE%MEKE is not associated!")
    return
  endif

  select case(CS%eke_src)
  case(EKE_PROG)
    if (CS%debug) then
      if (allocated(MEKE%mom_src)) &
        call hchksum(MEKE%mom_src, 'MEKE mom_src', G%HI, unscale=US%RZ3_T3_to_W_m2*US%L_to_Z**2)
      if (allocated(MEKE%mom_src_bh)) &
        call hchksum(MEKE%mom_src_bh, 'MEKE mom_src_bh', G%HI, unscale=US%RZ3_T3_to_W_m2*US%L_to_Z**2)
      if (allocated(MEKE%GME_snk)) &
        call hchksum(MEKE%GME_snk, 'MEKE GME_snk', G%HI, unscale=US%RZ3_T3_to_W_m2*US%L_to_Z**2)
      if (allocated(MEKE%GM_src)) &
        call hchksum(MEKE%GM_src, 'MEKE GM_src', G%HI, unscale=US%RZ3_T3_to_W_m2*US%L_to_Z**2)
      if (allocated(MEKE%MEKE)) &
        call hchksum(MEKE%MEKE, 'MEKE MEKE', G%HI, unscale=US%L_T_to_m_s**2)
      call uvchksum("MEKE SN_[uv]", SN_u, SN_v, G%HI, unscale=US%s_to_T, &
                    scalar_pair=.true.)
      call uvchksum("MEKE h[uv]", hu, hv, G%HI, haloshift=0, symmetric=.true., &
                    unscale=GV%H_to_m*US%L_to_m**2)
    endif

    sdt = dt*CS%MEKE_dtScale ! Scaled dt to use for time-stepping
    mass_neglect = GV%H_to_RZ * GV%H_subroundoff
    cdrag2 = CS%cdrag**2

    ! With a depth-dependent (and possibly strong) damping, it seems
    ! advisable to use Strang splitting between the damping and diffusion.
    damp_step = 1.
    if (CS%MEKE_KH >= 0. .or. CS%MEKE_K4 >= 0.) damp_step = 0.5
    sdt_damp = sdt * damp_step

    ! Calculate depth integrated mass exchange if doing advection [R Z L2 ~> kg]
    if (CS%MEKE_advection_factor>0.) then
      do j=js,je ; do I=is-1,ie
        baroHu(I,j) = 0.
      enddo ; enddo
      do k=1,nz
        do j=js,je ; do I=is-1,ie
          baroHu(I,j) = baroHu(I,j) + hu(I,j,k) * GV%H_to_RZ
        enddo ; enddo
      enddo
      do J=js-1,je ; do i=is,ie
        baroHv(i,J) = 0.
      enddo ; enddo
      do k=1,nz
        do J=js-1,je ; do i=is,ie
          baroHv(i,J) = baroHv(i,J) + hv(i,J,k) * GV%H_to_RZ
        enddo ; enddo
      enddo
      if (CS%MEKE_advection_bug) then
        ! This obviously incorrect code reproduces a bug in the original implementation of
        ! the MEKE advection.
        do j=js,je ; do I=is-1,ie
          baroHu(I,j) = hu(I,j,nz) * GV%H_to_RZ
        enddo ; enddo
        do J=js-1,je ; do i=is,ie
          baroHv(i,J) = hv(i,J,nz) * GV%H_to_RZ
        enddo ; enddo
      endif
    endif

    ! Calculate drag_rate_visc(i,j) which accounts for the model bottom mean flow
    if (CS%visc_drag .and. allocated(visc%Kv_bbl_u) .and. allocated(visc%Kv_bbl_v)) then
      !$OMP parallel do default(shared)
      do j=js,je ; do I=is-1,ie
        drag_vel_u(I,j) = 0.0
        if ((G%mask2dCu(I,j) > 0.0) .and. (visc%bbl_thick_u(I,j) > 0.0)) &
          drag_vel_u(I,j) = visc%Kv_bbl_u(I,j) / visc%bbl_thick_u(I,j)
      enddo ; enddo
      !$OMP parallel do default(shared)
      do J=js-1,je ; do i=is,ie
        drag_vel_v(i,J) = 0.0
        if ((G%mask2dCv(i,J) > 0.0) .and. (visc%bbl_thick_v(i,J) > 0.0)) &
          drag_vel_v(i,J) = visc%Kv_bbl_v(i,J) / visc%bbl_thick_v(i,J)
      enddo ; enddo

      !$OMP parallel do default(shared)
      do j=js,je ; do i=is,ie
        drag_rate_visc(i,j) = (0.25*G%IareaT(i,j) * &
                (((G%areaCu(I-1,j)*drag_vel_u(I-1,j)) + &
                  (G%areaCu(I,j)*drag_vel_u(I,j))) + &
                 ((G%areaCv(i,J-1)*drag_vel_v(i,J-1)) + &
                  (G%areaCv(i,J)*drag_vel_v(i,J))) ) )
      enddo ; enddo
    else
      !$OMP parallel do default(shared)
      do j=js,je ; do i=is,ie
        drag_rate_visc(i,j) = 0.
      enddo ; enddo
    endif

    !$OMP parallel do default(shared)
    do j=js-1,je+1
      do i=is-1,ie+1 ; mass(i,j) = 0.0 ; enddo
      do k=1,nz ; do i=is-1,ie+1
        mass(i,j) = mass(i,j) + G%mask2dT(i,j) * (GV%H_to_RZ * h(i,j,k)) ! [R Z ~> kg m-2]
      enddo ; enddo
      do i=is-1,ie+1
        I_mass(i,j) = 0.0
        if (mass(i,j) > 0.0) I_mass(i,j) = 1.0 / mass(i,j) ! [R-1 Z-1 ~> m2 kg-1]
      enddo
    enddo

    if (CS%fixed_total_depth) then
      if (GV%Boussinesq) then
        !$OMP parallel do default(shared)
        do j=js-1,je+1 ; do i=is-1,ie+1
          depth_tot(i,j) = max(G%meanSL(i,j) + G%bathyT(i,j), 0.0) * GV%Z_to_H
        enddo ; enddo
      else
        !$OMP parallel do default(shared)
        do j=js-1,je+1 ; do i=is-1,ie+1
          depth_tot(i,j) = max(G%meanSL(i,j) + G%bathyT(i,j), 0.0) * CS%rho_fixed_total_depth * GV%RZ_to_H
        enddo ; enddo
      endif
    else
      !$OMP parallel do default(shared)
      do j=js-1,je+1 ; do i=is-1,ie+1
        depth_tot(i,j) = mass(i,j) * GV%RZ_to_H
      enddo ; enddo
    endif

    if (CS%initialize) then
      call MEKE_equilibrium(CS, MEKE, G, GV, US, SN_u, SN_v, drag_rate_visc, I_mass, depth_tot)
      CS%initialize = .false.
    endif

    ! Calculates bottomFac2, barotrFac2 and LmixScale
    call MEKE_lengthScales(CS, MEKE, G, GV, US, SN_u, SN_v, MEKE%MEKE, depth_tot, bottomFac2, barotrFac2, LmixScale)
    if (CS%debug) then
      if (CS%visc_drag) &
        call uvchksum("MEKE drag_vel_[uv]", drag_vel_u, drag_vel_v, G%HI, &
                      unscale=GV%H_to_mks*US%s_to_T, scalar_pair=.true.)
      call hchksum(mass, 'MEKE mass',G%HI,haloshift=1, unscale=US%RZ_to_kg_m2)
      call hchksum(drag_rate_visc, 'MEKE drag_rate_visc', G%HI, unscale=GV%H_to_mks*US%s_to_T)
      call hchksum(bottomFac2, 'MEKE bottomFac2', G%HI)
      call hchksum(barotrFac2, 'MEKE barotrFac2', G%HI)
      call hchksum(LmixScale, 'MEKE LmixScale', G%HI, unscale=US%L_to_m)
    endif

    if (allocated(MEKE%Le)) then
      !$OMP parallel do default(shared)
      do j=js,je ; do i=is,ie
        MEKE%Le(i,j) = LmixScale(i,j)
      enddo ; enddo
    endif

    ! Aggregate sources of MEKE (background, frictional and GM)
    !$OMP parallel do default(shared)
    do j=js,je ; do i=is,ie
      src(i,j) = CS%MEKE_BGsrc
    enddo ; enddo

    ! Initialize diagnostics
    if (CS%id_src_adv > 0) src_adv(is:ie, js:je) = 0.
    if (CS%id_src_GM > 0) src_GM(is:ie, js:je) = 0.
    if (CS%id_src_mom_lp > 0) src_mom_lp(is:ie, js:je) = 0.
    if (CS%id_src_mom_bh > 0) src_mom_bh(is:ie, js:je) = 0.
    if (CS%id_src_mom_K4 > 0) src_mom_K4(is:ie, js:je) = 0.
    if (CS%id_src_btm_drag > 0) src_btm_drag(is:ie, js:je) = 0.

    ! Identify any damped diagnostics in first stage of Strang splitting
    any_damping_diags_s1 = any([ &
        CS%id_src_GM > 0, &
        CS%id_src_mom_lp > 0, &
        CS%id_src_mom_bh > 0, &
        CS%id_src_btm_drag > 0 &
    ])

    ! Identify any damped diagnostics
    any_damping_diags = any([ &
        any_damping_diags_s1, &
        CS%id_src_adv > 0, &
        CS%id_src_mom_K4 > 0 &
    ])

    if (CS%MEKE_FrCoeff > 0.) then
      !$OMP parallel do default(shared)
      do j=js,je ; do i=is,ie
        src(i,j) = src(i,j) - CS%MEKE_FrCoeff * I_mass(i,j) * MEKE%mom_src(i,j)
      enddo ; enddo
    endif

    if (allocated(MEKE%mom_src_bh)) then
      if (CS%MEKE_bhFrCoeff > 0. .and. CS%MEKE_FrCoeff > 0.) then
        bh_coeff = CS%MEKE_bhFrCoeff - CS%MEKE_FrCoeff
      else
        bh_coeff = CS%MEKE_bhFrCoeff
      endif

      !$OMP parallel do default(shared)
      do j=js,je ; do i=is,ie
        src(i,j) = src(i,j) - bh_coeff * I_mass(i,j) * MEKE%mom_src_bh(i,j)
      enddo ; enddo

      if (CS%id_src_mom_lp > 0) then
        !$OMP parallel do default(shared)
        do j=js,je ; do i=is,ie
          src_mom_lp(i,j) = -CS%MEKE_FrCoeff * I_mass(i,j) &
              * (MEKE%mom_src(i,j) - MEKE%mom_src_bh(i,j))
        enddo ; enddo
      endif

      if (CS%id_src_mom_bh > 0) then
        !$OMP parallel do default(shared)
        do j=js,je ; do i=is,ie
          src_mom_bh(i,j) = -CS%MEKE_bhFrCoeff * I_mass(i,j) * MEKE%mom_src_bh(i,j)
        enddo ; enddo
      endif
    endif

    if (allocated(MEKE%GME_snk)) then
      !$OMP parallel do default(shared)
      do j=js,je ; do i=is,ie
        src(i,j) = src(i,j) - CS%MEKE_GMECoeff*I_mass(i,j)*MEKE%GME_snk(i,j)
      enddo ; enddo
    endif

    if (allocated(MEKE%GM_src)) then
      if (CS%GM_src_alt) then
        !$OMP parallel do default(shared)
        do j=js,je ; do i=is,ie
          src(i,j) = src(i,j) - CS%MEKE_GMcoeff*MEKE%GM_src(i,j) / &
                     (GV%H_to_RZ * MAX(CS%MEKE_min_depth_tot, depth_tot(i,j)))
        enddo ; enddo
      else
        !$OMP parallel do default(shared)
        do j=js,je ; do i=is,ie
          src(i,j) = src(i,j) - CS%MEKE_GMcoeff*I_mass(i,j)*MEKE%GM_src(i,j)
        enddo ; enddo

        do j=js,je ; do i=is,ie
          src_GM(i,j) = -CS%MEKE_GMcoeff*I_mass(i,j)*MEKE%GM_src(i,j)
        enddo ; enddo
      endif
    endif

    if (CS%MEKE_equilibrium_restoring) then
      call MEKE_equilibrium_restoring(CS, G, GV, US, SN_u, SN_v, depth_tot, &
                                      equilibrium_value)
      do j=js,je ; do i=is,ie
        src(i,j) = src(i,j) - CS%MEKE_restoring_rate*(MEKE%MEKE(i,j) - equilibrium_value(i,j))
      enddo ; enddo
    endif

    if (CS%debug) then
      call hchksum(src, "MEKE src", G%HI, haloshift=0, unscale=US%L_to_m**2*US%s_to_T**3)
    endif

    ! Increase EKE by a full time-steps worth of source
    !$OMP parallel do default(shared)
    do j=js,je ; do i=is,ie
      MEKE_current(i,j) = MEKE%MEKE(i,j)
      MEKE%MEKE(i,j) = (MEKE%MEKE(i,j) + sdt*src(i,j))*G%mask2dT(i,j)
    enddo ; enddo

    if (use_drag_rate) then
      ! Calculate a viscous drag rate (includes BBL contributions from mean flow and eddies)
      !$OMP parallel do default(shared)
      do j=js,je ; do i=is,ie
        drag_rate(i,j) = (GV%H_to_RZ * I_mass(i,j)) * sqrt( drag_rate_visc(i,j)**2 + &
                 cdrag2 * ( max(0.0, 2.0*bottomFac2(i,j)*MEKE%MEKE(i,j)) + CS%MEKE_Uscale**2 ) )
      enddo ; enddo
    else
      !$OMP parallel do default(shared)
      do j=js,je ; do i=is,ie
        drag_rate(i,j) = 0.
      enddo ; enddo
    endif

    ! First stage of Strang splitting

    !$OMP parallel do default(shared)
    do j=js,je ; do i=is,ie
      damp_rate(i,j) = CS%MEKE_damping + drag_rate(i,j) * bottomFac2(i,j)

      if (MEKE%MEKE(i,j) < 0.) damp_rate(i,j) = 0.
      ! notice that the above line ensures a damping only if MEKE is positive,
      ! while leaving MEKE unchanged if it is negative
    enddo ; enddo

    ! NOTE: MEKE%MEKE cannot use `damping` since we must preserve the existing
    !   bit-reproducible solution.
    !$OMP parallel do default(shared)
    do j=js,je ; do i=is,ie
      MEKE%MEKE(i,j) =  MEKE%MEKE(i,j) / (1. + sdt_damp * damp_rate(i,j))
    enddo ; enddo

    if (any_damping_diags_s1) then
      !$OMP parallel do default(shared)
      do j=js,je ; do i=is,ie
        damping(i,j) = 1. / (1. + sdt_damp * damp_rate(i,j))
      enddo ; enddo

      if (CS%id_decay > 0) then
        !$OMP parallel do default(shared)
        do j=js,je ; do i=is,ie
          MEKE_decay(i,j) = damp_rate(i,j) * G%mask2dT(i,j)
        enddo ; enddo
      endif

      if (CS%id_src_GM > 0) then
        !$OMP parallel do default(shared)
        do j=js,je ; do i=is,ie
          src_GM(i,j) = src_GM(i,j) * damping(i,j)
        enddo ; enddo
      endif

      if (CS%id_src_mom_lp > 0) then
        !$OMP parallel do default(shared)
        do j=js,je ; do i=is,ie
          src_mom_lp(i,j) = src_mom_lp(i,j) * damping(i,j)
        enddo ; enddo
      endif

      if (CS%id_src_mom_bh > 0) then
        !$OMP parallel do default(shared)
        do j=js,je ; do i=is,ie
          src_mom_bh(i,j) = src_mom_bh(i,j) * damping(i,j)
        enddo ; enddo
      endif

      if (CS%id_src_btm_drag > 0) then
        !$OMP parallel do default(shared)
        do j=js,je ; do i=is,ie
          src_btm_drag(i,j) = -MEKE_current(i,j) * ( &
              damp_step * (damp_rate(i,j) * damping(i,j)) &
          )
        enddo ; enddo

        ! Store the effective damping rate if sdt is split
        if (CS%MEKE_KH >= 0. .or. CS%MEKE_K4 >= 0.) then
          !$OMP parallel do default(shared)
          do j=js,je ; do i=is,ie
            damp_rate_s1(i,j) = damp_rate(i,j) * damping(i,j)
          enddo ; enddo
        endif
      endif
    endif

    if (CS%kh_flux_enabled .or. CS%MEKE_K4 >= 0.0) then
      ! Update MEKE in the halos for lateral or bi-harmonic diffusion
      call cpu_clock_begin(CS%id_clock_pass)
      call do_group_pass(CS%pass_MEKE, G%Domain)
      call cpu_clock_end(CS%id_clock_pass)
    endif

    if (CS%MEKE_K4 >= 0.0) then
      ! Calculate Laplacian of MEKE using MEKE_uflux and MEKE_vflux as temporary work space.
      !$OMP parallel do default(shared)
      do j=js-1,je+1 ; do I=is-2,ie+1
        ! MEKE_uflux is used here as workspace with units of [L2 T-2 ~> m2 s-2].
        MEKE_uflux(I,j) = (G%dy_Cu(I,j)*G%IdxCu_OBCmask(I,j)) * &
            (MEKE%MEKE(i+1,j) - MEKE%MEKE(i,j))
      ! This would have units of [R Z L2 T-2 ~> kg s-2]
      ! MEKE_uflux(I,j) = ((G%dy_Cu(I,j)*G%IdxCu(I,j)) * &
      !     ((2.0*mass(i,j)*mass(i+1,j)) / ((mass(i,j)+mass(i+1,j)) + mass_neglect)) ) * &
      !     (MEKE%MEKE(i+1,j) - MEKE%MEKE(i,j))
      enddo ; enddo
      !$OMP parallel do default(shared)
      do J=js-2,je+1 ; do i=is-1,ie+1
        ! MEKE_vflux is used here as workspace with units of [L2 T-2 ~> m2 s-2].
        MEKE_vflux(i,J) = (G%dx_Cv(i,J)*G%IdyCv_OBCmask(i,J)) * &
            (MEKE%MEKE(i,j+1) - MEKE%MEKE(i,j))
      ! This would have units of [R Z L2 T-2 ~> kg s-2]
      ! MEKE_vflux(i,J) = ((G%dx_Cv(i,J)*G%IdyCv(i,J)) * &
      !     ((2.0*mass(i,j)*mass(i,j+1)) / ((mass(i,j)+mass(i,j+1)) + mass_neglect)) ) * &
      !     (MEKE%MEKE(i,j+1) - MEKE%MEKE(i,j))
      enddo ; enddo

      !$OMP parallel do default(shared)
      do j=js-1,je+1 ; do i=is-1,ie+1 ! del2MEKE has units [T-2 ~> s-2].
        del2MEKE(i,j) = G%IareaT(i,j) * &
            ((MEKE_uflux(I,j) - MEKE_uflux(I-1,j)) + (MEKE_vflux(i,J) - MEKE_vflux(i,J-1)))
      enddo ; enddo

      ! Bi-harmonic diffusion of MEKE
      !$OMP parallel do default(shared) private(K4_here,Inv_K4_max)
      do j=js,je ; do I=is-1,ie
        K4_here = CS%MEKE_K4 ! [L4 T-1 ~> m4 s-1]
        ! Limit Kh to avoid CFL violations.
        Inv_K4_max = 64.0 * sdt * ((G%dy_Cu(I,j)*G%IdxCu(I,j)) * &
                     max(G%IareaT(i,j), G%IareaT(i+1,j)))**2
        if (K4_here*Inv_K4_max > 0.3) K4_here = 0.3 / Inv_K4_max

        ! Here the units of MEKE_uflux are [R Z L4 T-3 ~> kg m2 s-3].
        MEKE_uflux(I,j) = ((K4_here * (G%dy_Cu(I,j)*G%IdxCu(I,j))) * &
            ((2.0*mass(i,j)*mass(i+1,j)) / ((mass(i,j)+mass(i+1,j)) + mass_neglect)) ) * &
            (del2MEKE(i+1,j) - del2MEKE(i,j))
      enddo ; enddo
      !$OMP parallel do default(shared) private(K4_here,Inv_K4_max)
      do J=js-1,je ; do i=is,ie
        K4_here = CS%MEKE_K4 ! [L4 T-1 ~> m4 s-1]
        Inv_K4_max = 64.0 * sdt * ((G%dx_Cv(i,J)*G%IdyCv(i,J)) * max(G%IareaT(i,j), G%IareaT(i,j+1)))**2
        if (K4_here*Inv_K4_max > 0.3) K4_here = 0.3 / Inv_K4_max

        ! Here the units of MEKE_vflux are [R Z L4 T-3 ~> kg m2 s-3].
        MEKE_vflux(i,J) = ((K4_here * (G%dx_Cv(i,J)*G%IdyCv(i,J))) * &
            ((2.0*mass(i,j)*mass(i,j+1)) / ((mass(i,j)+mass(i,j+1)) + mass_neglect)) ) * &
            (del2MEKE(i,j+1) - del2MEKE(i,j))
      enddo ; enddo
      ! Store change in MEKE arising from the bi-harmonic in del4MEKE [L2 T-2 ~> m2 s-2].
      !$OMP parallel do default(shared)
      do j=js,je ; do i=is,ie
        del4MEKE(i,j) = (sdt*(G%IareaT(i,j)*I_mass(i,j))) * &
            ((MEKE_uflux(I-1,j) - MEKE_uflux(I,j)) + &
             (MEKE_vflux(i,J-1) - MEKE_vflux(i,J)))
        src_mom_K4(i,j) = (G%IareaT(i,j)*I_mass(i,j))  * &
            ((MEKE_uflux(I-1,j) - MEKE_uflux(I,j)) + &
             (MEKE_vflux(i,J-1) - MEKE_vflux(i,J)))
      enddo ; enddo
    endif !

    if (CS%kh_flux_enabled) then
      ! Lateral diffusion of MEKE
      Kh_here = max(0., CS%MEKE_Kh)
      !$OMP parallel do default(shared) firstprivate(Kh_here) private(Inv_Kh_max)
      do j=js,je ; do I=is-1,ie
        ! Limit Kh to avoid CFL violations.
        if (allocated(MEKE%Kh)) &
          Kh_here = max(0., CS%MEKE_Kh) + &
              CS%KhMEKE_Fac*0.5*(MEKE%Kh(i,j)+MEKE%Kh(i+1,j))
        if (allocated(MEKE%Kh_diff)) &
          Kh_here = max(0.,CS%MEKE_Kh) + &
              CS%KhMEKE_Fac*0.5*(MEKE%Kh_diff(i,j)+MEKE%Kh_diff(i+1,j))
        Inv_Kh_max = 2.0*sdt * ((G%dy_Cu(I,j)*G%IdxCu(I,j)) * &
                     max(G%IareaT(i,j),G%IareaT(i+1,j)))
        if (Kh_here*Inv_Kh_max > 0.25) Kh_here = 0.25 / Inv_Kh_max
        Kh_u(I,j) = Kh_here

        ! Here the units of MEKE_uflux and MEKE_vflux are [R Z L4 T-3 ~> kg m2 s-3].
        MEKE_uflux(I,j) = ((Kh_here * (G%dy_Cu(I,j)*G%IdxCu(I,j))) * &
            ((2.0*mass(i,j)*mass(i+1,j)) / ((mass(i,j)+mass(i+1,j)) + mass_neglect)) ) * &
            (MEKE%MEKE(i,j) - MEKE%MEKE(i+1,j))
      enddo ; enddo
      !$OMP parallel do default(shared) firstprivate(Kh_here) private(Inv_Kh_max)
      do J=js-1,je ; do i=is,ie
        if (allocated(MEKE%Kh)) &
          Kh_here = max(0.,CS%MEKE_Kh) + CS%KhMEKE_Fac * 0.5*(MEKE%Kh(i,j)+MEKE%Kh(i,j+1))
        if (allocated(MEKE%Kh_diff)) &
          Kh_here = max(0.,CS%MEKE_Kh) + CS%KhMEKE_Fac * 0.5*(MEKE%Kh_diff(i,j)+MEKE%Kh_diff(i,j+1))
        Inv_Kh_max = 2.0*sdt * ((G%dx_Cv(i,J)*G%IdyCv(i,J)) * max(G%IareaT(i,j),G%IareaT(i,j+1)))
        if (Kh_here*Inv_Kh_max > 0.25) Kh_here = 0.25 / Inv_Kh_max
        Kh_v(i,J) = Kh_here

        ! Here the units of MEKE_uflux and MEKE_vflux are [R Z L4 T-3 ~> kg m2 s-3].
        MEKE_vflux(i,J) = ((Kh_here * (G%dx_Cv(i,J)*G%IdyCv(i,J))) * &
            ((2.0*mass(i,j)*mass(i,j+1)) / ((mass(i,j)+mass(i,j+1)) + mass_neglect)) ) * &
            (MEKE%MEKE(i,j) - MEKE%MEKE(i,j+1))
      enddo ; enddo
      if (CS%MEKE_advection_factor>0.) then
        advFac = CS%MEKE_advection_factor / sdt ! [T-1 ~> s-1]
        !$OMP parallel do default(shared)
        do j=js,je ; do I=is-1,ie
          ! Here the units of the quantities added to MEKE_uflux are [R Z L4 T-3 ~> kg m2 s-3].
          if (baroHu(I,j)>0.) then
            MEKE_uflux(I,j) = MEKE_uflux(I,j) + baroHu(I,j)*MEKE%MEKE(i,j)*advFac
          elseif (baroHu(I,j)<0.) then
            MEKE_uflux(I,j) = MEKE_uflux(I,j) + baroHu(I,j)*MEKE%MEKE(i+1,j)*advFac
          endif
        enddo ; enddo
        !$OMP parallel do default(shared)
        do J=js-1,je ; do i=is,ie
          ! Here the units of the quantities added to MEKE_vflux are [R Z L4 T-3 ~> kg m2 s-3].
          if (baroHv(i,J)>0.) then
            MEKE_vflux(i,J) = MEKE_vflux(i,J) + baroHv(i,J)*MEKE%MEKE(i,j)*advFac
          elseif (baroHv(i,J)<0.) then
            MEKE_vflux(i,J) = MEKE_vflux(i,J) + baroHv(i,J)*MEKE%MEKE(i,j+1)*advFac
          endif
        enddo ; enddo
      endif

      !$OMP parallel do default(shared)
      do j=js,je ; do i=is,ie
        MEKE%MEKE(i,j) = MEKE%MEKE(i,j) + (sdt*(G%IareaT(i,j)*I_mass(i,j))) * &
            ((MEKE_uflux(I-1,j) - MEKE_uflux(I,j)) + &
             (MEKE_vflux(i,J-1) - MEKE_vflux(i,J)))
      enddo ; enddo

      if (CS%id_src_adv > 0) then
        !$OMP parallel do default(shared)
        do j=js,je ; do i=is,ie
          src_adv(i,j) = (G%IareaT(i,j)*I_mass(i,j)) * &
              ((MEKE_uflux(I-1,j) - MEKE_uflux(I,j)) + &
               (MEKE_vflux(i,J-1) - MEKE_vflux(i,J)))
        enddo ; enddo
      endif
    endif ! MEKE_KH>0

    ! Add on bi-harmonic tendency
    if (CS%MEKE_K4 >= 0.0) then
      !$OMP parallel do default(shared)
      do j=js,je ; do i=is,ie
        MEKE%MEKE(i,j) = MEKE%MEKE(i,j) + del4MEKE(i,j)
      enddo ; enddo
    endif

    ! Second stage of Strang splitting
    if (CS%MEKE_KH >= 0.0 .or. CS%MEKE_K4 >= 0.0) then
      ! Recalculate the drag rate, since MEKE has changed.
      if (use_drag_rate) then
        !$OMP parallel do default(shared)
        do j=js,je ; do i=is,ie
          drag_rate(i,j) = (GV%H_to_RZ * I_mass(i,j)) * sqrt( drag_rate_visc(i,j)**2 + &
                 cdrag2 * ( max(0.0, 2.0*bottomFac2(i,j)*MEKE%MEKE(i,j)) + CS%MEKE_Uscale**2 ) )
        enddo ; enddo
      endif

      !$OMP parallel do default(shared)
      do j=js,je ; do i=is,ie
        damp_rate(i,j) = CS%MEKE_damping + drag_rate(i,j) * bottomFac2(i,j)

        if (MEKE%MEKE(i,j) < 0.) damp_rate(i,j) = 0.
        ! notice that the above line ensures a damping only if MEKE is positive,
        ! while leaving MEKE unchanged if it is negative
      enddo ; enddo

      ! NOTE: MEKE%MEKE cannot use `damping` since we must preserve the
      !   existing bit-reproducible solution.
      !$OMP parallel do default(shared)
      do j=js,je ; do i=is,ie
        MEKE%MEKE(i,j) =  MEKE%MEKE(i,j) / (1. + sdt_damp * damp_rate(i,j))
      enddo ; enddo

      if (any_damping_diags) then
        !$OMP parallel do default(shared)
        do j=js,je ; do i=is,ie
          damping(i,j) = 1. / (1. + sdt_damp * damp_rate(i,j))
        enddo ; enddo

        if (CS%id_decay > 0) then
          !$OMP parallel do default(shared)
          do j=js,je ; do i=is,ie
            MEKE_decay(i,j) = damp_rate(i,j) * G%mask2dT(i,j)
          enddo ; enddo
        endif

        if (CS%id_src_GM > 0) then
          !$OMP parallel do default(shared)
          do j=js,je ; do i=is,ie
            src_GM(i,j) = src_GM(i,j) * damping(i,j)
          enddo ; enddo
        endif

        if (CS%id_src_mom_lp > 0) then
          !$OMP parallel do default(shared)
          do j=js,je ; do i=is,ie
            src_mom_lp(i,j) = src_mom_lp(i,j) * damping(i,j)
          enddo ; enddo
        endif

        if (CS%id_src_mom_bh > 0) then
          !$OMP parallel do default(shared)
          do j=js,je ; do i=is,ie
            src_mom_bh(i,j) = src_mom_bh(i,j) * damping(i,j)
          enddo ; enddo
        endif

        if (CS%id_src_adv > 0) then
          !$OMP parallel do default(shared)
          do j=js,je ; do i=is,ie
            src_adv(i,j) = src_adv(i,j) * damping(i,j)
          enddo ; enddo
        endif

        if (CS%id_src_mom_K4 > 0) then
          !$OMP parallel do default(shared)
          do j=js,je ; do i=is,ie
            src_mom_K4(i,j) = src_mom_K4(i,j) * damping(i,j)
          enddo ; enddo
        endif

        if (CS%id_src_btm_drag > 0) then
          !$OMP parallel do default(shared)
          do j=js,je ; do i=is,ie
            src_btm_drag(i,j) = -MEKE_current(i,j) * (damp_step &
                * ((damp_rate(i,j) + damp_rate_s1(i,j)) * damping(i,j)) &
            )
          enddo ; enddo
        endif
      endif
    endif ! MEKE_KH>=0

    if (CS%debug) then
      call hchksum(MEKE%MEKE, "MEKE post-update MEKE", G%HI, haloshift=0, unscale=US%L_T_to_m_s**2)
    endif

  case(EKE_FILE)
    call time_interp_external(CS%eke_handle, Time, data_eke, scale=US%m_s_to_L_T**2)
    do j=js,je ; do i=is,ie
      MEKE%MEKE(i,j) = data_eke(i,j) * G%mask2dT(i,j)
    enddo ; enddo
    call MEKE_lengthScales(CS, MEKE, G, GV, US, SN_u, SN_v, MEKE%MEKE, depth_tot, bottomFac2, barotrFac2, LmixScale)
  case(EKE_DBCLIENT)
    call pass_vector(u, v, G%Domain)
    call MEKE_lengthScales(CS, MEKE, G, GV, US, SN_u, SN_v, MEKE%MEKE, depth_tot, bottomFac2, barotrFac2, LmixScale)
    call ML_MEKE_calculate_features(G, GV, US, CS, MEKE%Rd_dx_h, u, v, tv, h, dt, features_array)
    call predict_MEKE(G, US, CS, SIZE(h), Time, features_array, MEKE%MEKE)
  case default
    call MOM_error(FATAL,"Invalid method specified for calculating EKE")
  end select

  if (CS%MEKE_positive) then
    !$OMP parallel do default(shared)
    do j=js,je ; do i=is,ie
      MEKE%MEKE(i,j) = MAX(0., MEKE%MEKE(i,j))
    enddo ; enddo
  endif

  call cpu_clock_begin(CS%id_clock_pass)
  call do_group_pass(CS%pass_MEKE, G%Domain)
  call cpu_clock_end(CS%id_clock_pass)

  ! Calculate diffusivity for main model to use
  if (CS%MEKE_KhCoeff>0.) then
    if (.not.CS%MEKE_GEOMETRIC) then
      if (CS%use_old_lscale) then
        if (CS%Rd_as_max_scale) then
          !$OMP parallel do default(shared)
          do j=js,je ; do i=is,ie
            MEKE%Kh(i,j) = (CS%MEKE_KhCoeff * &
                       sqrt(2.*max(0.,barotrFac2(i,j)*MEKE%MEKE(i,j))*G%areaT(i,j)) ) * &
                       min(MEKE%Rd_dx_h(i,j), 1.0)
          enddo ; enddo
        else
          !$OMP parallel do default(shared)
          do j=js,je ; do i=is,ie
            MEKE%Kh(i,j) = CS%MEKE_KhCoeff * &
                sqrt(2.*max(0., barotrFac2(i,j)*MEKE%MEKE(i,j))*G%areaT(i,j))
          enddo ; enddo
        endif
      else
        !$OMP parallel do default(shared)
        do j=js,je ; do i=is,ie
          MEKE%Kh(i,j) = CS%MEKE_KhCoeff * &
              sqrt(2.*max(0., barotrFac2(i,j)*MEKE%MEKE(i,j))) * LmixScale(i,j)
        enddo ; enddo
      endif
    endif
  endif

  ! Calculate viscosity for the main model to use
  if (CS%viscosity_coeff_Ku /=0.) then
    do j=js,je ; do i=is,ie
      MEKE%Ku(i,j) = CS%viscosity_coeff_Ku * sqrt(2.*max(0.,MEKE%MEKE(i,j))) * LmixScale(i,j)
    enddo ; enddo
  endif

  if (CS%viscosity_coeff_Au /=0.) then
    do j=js,je ; do i=is,ie
      MEKE%Au(i,j) = CS%viscosity_coeff_Au * sqrt(2.*max(0.,MEKE%MEKE(i,j))) * LmixScale(i,j)**3
    enddo ; enddo
  endif

  if (allocated(MEKE%Kh) .or. allocated(MEKE%Ku) .or. allocated(MEKE%Au) &
      .or. allocated(MEKE%Le)) then
    call cpu_clock_begin(CS%id_clock_pass)
    call do_group_pass(CS%pass_Kh, G%Domain)
    call cpu_clock_end(CS%id_clock_pass)
  endif

  ! Offer fields for averaging.
  if (any([CS%id_Ue, CS%id_Ub, CS%id_Ut] > 0)) &
    tmp(:,:) = 0.
  if (CS%id_MEKE>0) call post_data(CS%id_MEKE, MEKE%MEKE, CS%diag)
  if (CS%id_Ue>0) then
    do j=js,je ; do i=is,ie
      tmp(i,j) = sqrt(max(0., 2. * MEKE%MEKE(i,j)))
    enddo ; enddo
    call post_data(CS%id_Ue, tmp, CS%diag)
  endif
  if (CS%id_Ub>0) then
    do j=js,je ; do i=is,ie
      tmp(i,j) = sqrt(max(0., 2. * MEKE%MEKE(i,j) * bottomFac2(i,j)))
    enddo ; enddo
    call post_data(CS%id_Ub, tmp, CS%diag)
  endif
  if (CS%id_Ut>0) then
    do j=js,je ; do i=is,ie
      tmp(i,j) = sqrt(max(0., 2. * MEKE%MEKE(i,j) * barotrFac2(i,j)))
    enddo ; enddo
    call post_data(CS%id_Ut, tmp, CS%diag)
  endif
  if (CS%id_Kh>0) call post_data(CS%id_Kh, MEKE%Kh, CS%diag)
  if (CS%id_Ku>0) call post_data(CS%id_Ku, MEKE%Ku, CS%diag)
  if (CS%id_Au>0) call post_data(CS%id_Au, MEKE%Au, CS%diag)
  if (CS%id_KhMEKE_u>0) call post_data(CS%id_KhMEKE_u, Kh_u, CS%diag)
  if (CS%id_KhMEKE_v>0) call post_data(CS%id_KhMEKE_v, Kh_v, CS%diag)
  if (CS%id_src>0) call post_data(CS%id_src, src, CS%diag)
  if (CS%id_src_adv>0) call post_data(CS%id_src_adv, src_adv, CS%diag)
  if (CS%id_src_mom_K4>0) call post_data(CS%id_src_mom_K4, src_mom_K4, CS%diag)
  if (CS%id_src_btm_drag>0) call post_data(CS%id_src_btm_drag, src_btm_drag, CS%diag)
  if (CS%id_src_GM>0) call post_data(CS%id_src_GM, src_GM, CS%diag)
  if (CS%id_src_mom_lp>0) call post_data(CS%id_src_mom_lp, src_mom_lp, CS%diag)
  if (CS%id_src_mom_bh>0) call post_data(CS%id_src_mom_bh, src_mom_bh, CS%diag)
  if (CS%id_decay>0) call post_data(CS%id_decay, MEKE_decay, CS%diag)
  if (CS%id_GM_src>0) call post_data(CS%id_GM_src, MEKE%GM_src, CS%diag)
  if (CS%id_mom_src>0) call post_data(CS%id_mom_src, MEKE%mom_src, CS%diag)
  if (CS%id_mom_src_bh>0) call post_data(CS%id_mom_src_bh, MEKE%mom_src_bh, CS%diag)
  if (CS%id_GME_snk>0) call post_data(CS%id_GME_snk, MEKE%GME_snk, CS%diag)
  if (CS%id_Le>0) call post_data(CS%id_Le, LmixScale, CS%diag)
  if (CS%id_gamma_b>0) then
    do j=js,je ; do i=is,ie
      bottomFac2(i,j) = sqrt(bottomFac2(i,j))
    enddo ; enddo
    call post_data(CS%id_gamma_b, bottomFac2, CS%diag)
  endif
  if (CS%id_gamma_t>0) then
    do j=js,je ; do i=is,ie
      barotrFac2(i,j) = sqrt(barotrFac2(i,j))
    enddo ; enddo
    call post_data(CS%id_gamma_t, barotrFac2, CS%diag)
  endif

end procedure step_forward_MEKE
module procedure MEKE_equilibrium
  real :: beta ! Combined topographic and planetary vorticity gradient [T-1 L-1 ~> s-1 m-1]
  real :: SN   ! The local Eady growth rate [T-1 ~> s-1]
  real :: bottomFac2, barotrFac2    ! Vertical structure factors [nondim]
  real :: LmixScale, LRhines, LEady ! Various mixing length scales [L ~> m]
  real :: KhCoeff ! A copy of MEKE_KhCoeff from the control structure [nondim]
  real :: Kh    ! A lateral diffusivity [L2 T-1 ~> m2 s-1]
  real :: Ubg2  ! Background (tidal?) velocity squared [L2 T-2 ~> m2 s-2]
  real :: cd2   ! The square of the drag coefficient times unit conversion factors [H2 L-2 ~> nondim or kg2 m-6]
  real :: drag_rate ! The MEKE spindown timescale due to bottom drag [T-1 ~> s-1].
  real :: src   ! The sum of MEKE sources [L2 T-3 ~> W kg-1]
  real :: ldamping  ! The MEKE damping rate [T-1 ~> s-1].
  real :: EKE, EKEmin, EKEmax, EKEerr ! [L2 T-2 ~> m2 s-2]
  real :: resid, ResMin, ResMax ! Residuals [L2 T-3 ~> W kg-1]
  real :: FatH    ! Coriolis parameter at h points, used to compute topographic beta [T-1 ~> s-1]
  real :: beta_topo_x, beta_topo_y    ! Topographic PV gradients in x and y [T-1 L-1 ~> s-1 m-1]
  real :: h_neglect ! A negligible thickness [H ~> m or kg m-2]
  integer :: i, j, is, ie, js, je, n1, n2
  real :: tolerance ! Width of EKE bracket [L2 T-2 ~> m2 s-2].
  logical :: useSecant, debugIteration
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec

  debugIteration = .false.
  KhCoeff = CS%MEKE_KhCoeff
  Ubg2 = CS%MEKE_Uscale**2
  cd2 = CS%cdrag**2
  tolerance = 1.0e-12*US%m_s_to_L_T**2
  h_neglect = GV%H_subroundoff

!$OMP do
  do j=js,je ; do i=is,ie
    ! SN = 0.25*max( (SN_u(I,j) + SN_u(I-1,j)) + (SN_v(i,J) + SN_v(i,J-1)), 0.)
    ! This avoids extremes values in equilibrium solution due to bad values in SN_u, SN_v
    SN = min(SN_u(I,j), SN_u(I-1,j), SN_v(i,J), SN_v(i,J-1))

    if (CS%MEKE_equilibrium_alt) then
      MEKE%MEKE(i,j) = (CS%MEKE_GEOMETRIC_alpha * SN * depth_tot(i,j))**2 / cd2
    else
      FatH = 0.25*((G%CoriolisBu(I,J) + G%CoriolisBu(I-1,J-1)) + &
                   (G%CoriolisBu(I-1,J) + G%CoriolisBu(I,J-1))) ! Coriolis parameter at h points

      ! Since zero-bathymetry cells are masked, this avoids calculations on land
      if (CS%MEKE_topographic_beta == 0. .or. (depth_tot(i,j) == 0.0)) then
        beta_topo_x = 0. ; beta_topo_y = 0.
      else
        !### Consider different combinations of these estimates of topographic beta.
        beta_topo_x = -CS%MEKE_topographic_beta * FatH * 0.5 * ( &
                      (depth_tot(i+1,j)-depth_tot(i,j)) * G%IdxCu(I,j)  &
                  / max(depth_tot(i+1,j), depth_tot(i,j), h_neglect) &
              +       (depth_tot(i,j)-depth_tot(i-1,j)) * G%IdxCu(I-1,j) &
                  / max(depth_tot(i,j), depth_tot(i-1,j), h_neglect) )
        beta_topo_y = -CS%MEKE_topographic_beta * FatH * 0.5 * ( &
                      (depth_tot(i,j+1)-depth_tot(i,j)) * G%IdyCv(i,J)  &
                  / max(depth_tot(i,j+1), depth_tot(i,j), h_neglect) + &
                      (depth_tot(i,j)-depth_tot(i,j-1)) * G%IdyCv(i,J-1) &
                  / max(depth_tot(i,j), depth_tot(i,j-1), h_neglect) )
      endif
      beta =  sqrt(((G%dF_dx(i,j) + beta_topo_x)**2) + &
                   ((G%dF_dy(i,j) + beta_topo_y)**2) )

      if (KhCoeff*SN*I_mass(i,j)>0.) then
        ! Solve resid(E) = 0, where resid = Kh(E) * (SN)^2 - damp_rate(E) E
        EKEmin = 0.   ! Use the trivial root as the left bracket
        ResMin = 0.   ! Need to detect direction of left residual
        EKEmax = 0.01*US%m_s_to_L_T**2 ! First guess at right bracket
        useSecant = .false. ! Start using a bisection method

        ! First find right bracket for which resid<0
        resid = 1.0*US%m_to_L**2*US%T_to_s**3 ; n1 = 0
        do while (resid>0.)
          n1 = n1 + 1
          EKE = EKEmax
          call MEKE_lengthScales_0d(CS, US, G%areaT(i,j), beta, depth_tot(i,j), &
                                    MEKE%Rd_dx_h(i,j), SN, EKE, &
                                    bottomFac2, barotrFac2, LmixScale, LRhines, LEady)
          ! TODO: Should include resolution function in Kh
          Kh = (KhCoeff * sqrt(2.*barotrFac2*EKE) * LmixScale)
          src = Kh * (SN * SN)
          drag_rate = (GV%H_to_RZ * I_mass(i,j)) * sqrt(drag_rate_visc(i,j)**2 + cd2 * ( 2.0*bottomFac2*EKE + Ubg2 ) )
          ldamping = CS%MEKE_damping + drag_rate * bottomFac2
          resid = src - ldamping * EKE
          ! if (debugIteration) then
          !   write(0,*) n1, 'EKE=',EKE,'resid=',resid
          !   write(0,*) 'EKEmin=',EKEmin,'ResMin=',ResMin
          !   write(0,*) 'src=',src,'ldamping=',ldamping
          !   write(0,*) 'gamma-b=',bottomFac2,'gamma-t=',barotrFac2
          !   write(0,*) 'drag_visc=',drag_rate_visc(i,j),'Ubg2=',Ubg2
          ! endif
          if (resid>0.) then    ! EKE is to the left of the root
            EKEmin = EKE        ! so we move the left bracket here
            EKEmax = 10. * EKE  ! and guess again for the right bracket
            if (resid<ResMin) useSecant = .true.
            ResMin = resid
            if (EKEmax > 2.e17*US%m_s_to_L_T**2) then
              if (debugIteration) stop 'Something has gone very wrong'
              debugIteration = .true.
              resid = 1. ; n1 = 0
              EKEmin = 0. ; ResMin = 0.
              EKEmax = 0.01*US%m_s_to_L_T**2
              useSecant = .false.
            endif
          endif
        enddo ! while(resid>0.) searching for right bracket
        ResMax = resid

        ! Bisect the bracket
        n2 = 0 ; EKEerr = EKEmax - EKEmin
        do while (EKEerr > tolerance)
          n2 = n2 + 1
          if (useSecant) then
            EKE = EKEmin + (EKEmax - EKEmin) * (ResMin / (ResMin - ResMax))
          else
            EKE = 0.5 * (EKEmin + EKEmax)
          endif
          EKEerr = min( EKE-EKEmin, EKEmax-EKE )
          ! TODO: Should include resolution function in Kh
          Kh = (KhCoeff * sqrt(2.*barotrFac2*EKE) * LmixScale)
          src = Kh * (SN * SN)
          drag_rate = (GV%H_to_RZ * I_mass(i,j)) * sqrt( drag_rate_visc(i,j)**2 + cd2 * ( 2.0*bottomFac2*EKE + Ubg2 ) )
          ldamping = CS%MEKE_damping + drag_rate * bottomFac2
          resid = src - ldamping * EKE
          if (useSecant .and. resid>ResMin) useSecant = .false.
          if (resid>0.) then              ! EKE is to the left of the root
            EKEmin = EKE                  ! so we move the left bracket here
            if (resid<ResMin) useSecant = .true.
            ResMin = resid                ! Save this for the secant method
          elseif (resid<0.) then          ! EKE is to the right of the root
            EKEmax = EKE                  ! so we move the right bracket here
            ResMax = resid                ! Save this for the secant method
          else
            exit                          ! resid=0 => EKE is exactly at the root
          endif
          if (n2>200) stop 'Failing to converge?'
        enddo ! while(EKEmax-EKEmin>tolerance)

      else
        EKE = 0.
      endif
      MEKE%MEKE(i,j) = EKE
    endif
  enddo ; enddo

end procedure MEKE_equilibrium
module procedure MEKE_equilibrium_restoring
  real :: SN                      ! The local Eady growth rate [T-1 ~> s-1]
  integer :: i, j, is, ie, js, je ! local indices
  real :: cd2                     ! The square of the drag coefficient [nondim]
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec
  cd2 = CS%cdrag**2
  equilibrium_value(:,:) = 0.0

!$OMP do
  do j=js,je ; do i=is,ie
    ! SN = 0.25*max( (SN_u(I,j) + SN_u(I-1,j)) + (SN_v(i,J) + SN_v(i,J-1)), 0.)
    ! This avoids extremes values in equilibrium solution due to bad values in SN_u, SN_v
    SN = min(SN_u(I,j), SN_u(I-1,j), SN_v(i,J), SN_v(i,J-1))
    equilibrium_value(i,j) = (CS%MEKE_GEOMETRIC_alpha * SN * depth_tot(i,j))**2 / cd2
  enddo ; enddo

  if (CS%id_MEKE_equilibrium>0) call post_data(CS%id_MEKE_equilibrium, equilibrium_value, CS%diag)
end procedure MEKE_equilibrium_restoring
module procedure MEKE_lengthScales
  real, dimension(SZI_(G),SZJ_(G)) :: LRhines, LEady  ! Possible mixing length scales [L ~> m]
  real :: beta ! Combined topographic and planetary vorticity gradient [T-1 L-1 ~> s-1 m-1]
  real :: SN   ! The local Eady growth rate [T-1 ~> s-1]
  real :: FatH ! Coriolis parameter at h points [T-1 ~> s-1]
  real :: beta_topo_x, beta_topo_y  ! Topographic PV gradients in x and y [T-1 L-1 ~> s-1 m-1]
  real :: h_neglect ! A negligible thickness [H ~> m or kg m-2]
  integer :: i, j, is, ie, js, je
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec
  h_neglect = GV%H_subroundoff

!$OMP do
  do j=js,je ; do i=is,ie
    if (.not.CS%use_old_lscale) then
      if (CS%aEady > 0.) then
        SN = 0.25 * ( (SN_u(I,j) + SN_u(I-1,j)) + (SN_v(i,J) + SN_v(i,J-1)) )
      else
        SN = 0.
      endif
      FatH = 0.25* ( ( G%CoriolisBu(I,J) + G%CoriolisBu(I-1,J-1) ) + &
                     ( G%CoriolisBu(I-1,J) + G%CoriolisBu(I,J-1) ) )  ! Coriolis parameter at h points

      ! If depth_tot is zero, then a division by zero FPE will be raised.  In this
      ! case, we apply Adcroft's rule of reciprocals and set the term to zero.
      ! Since zero-bathymetry cells are masked, this should not affect values.
      if (CS%MEKE_topographic_beta == 0. .or. (depth_tot(i,j) == 0.0)) then
        beta_topo_x = 0. ; beta_topo_y = 0.
      else
        !### Consider different combinations of these estimates of topographic beta.
        beta_topo_x = -CS%MEKE_topographic_beta * FatH * 0.5 * ( &
                      (depth_tot(i+1,j)-depth_tot(i,j)) * G%IdxCu(I,j)  &
                 / max(depth_tot(i+1,j), depth_tot(i,j), h_neglect) &
              +       (depth_tot(i,j)-depth_tot(i-1,j)) * G%IdxCu(I-1,j) &
                 / max(depth_tot(i,j), depth_tot(i-1,j), h_neglect) )
        beta_topo_y = -CS%MEKE_topographic_beta * FatH * 0.5 * ( &
                      (depth_tot(i,j+1)-depth_tot(i,j)) * G%IdyCv(i,J)  &
                 / max(depth_tot(i,j+1), depth_tot(i,j), h_neglect) + &
                      (depth_tot(i,j)-depth_tot(i,j-1)) * G%IdyCv(i,J-1) &
                 / max(depth_tot(i,j), depth_tot(i,j-1), h_neglect) )
      endif
      beta =  sqrt(((G%dF_dx(i,j) + beta_topo_x)**2) + &
                   ((G%dF_dy(i,j) + beta_topo_y)**2) )

    else
      beta = 0.
    endif
    ! Returns bottomFac2, barotrFac2 and LmixScale
    call MEKE_lengthScales_0d(CS, US, G%areaT(i,j), beta, depth_tot(i,j),  &
                              MEKE%Rd_dx_h(i,j), SN, MEKE%MEKE(i,j), &
                              bottomFac2(i,j), barotrFac2(i,j), LmixScale(i,j), &
                              LRhines(i,j), LEady(i,j))
  enddo ; enddo
  if (CS%id_Lrhines>0) call post_data(CS%id_LRhines, LRhines, CS%diag)
  if (CS%id_Leady>0) call post_data(CS%id_LEady, LEady, CS%diag)

end procedure MEKE_lengthScales
module procedure MEKE_lengthScales_0d
  real :: Lgrid, Ldeform, Lfrict ! Length scales [L ~> m]
  real :: Ue  ! An eddy velocity [L T-1 ~> m s-1]
  Lgrid = sqrt(area)               ! Grid scale
  Ldeform = Lgrid * Rd_dx          ! Deformation scale
  Lfrict = depth_tot / CS%cdrag    ! Frictional arrest scale
  ! gamma_b^2 is the ratio of bottom eddy energy to mean column eddy energy
  ! used in calculating bottom drag
  bottomFac2 = CS%MEKE_CD_SCALE**2
  if (Lfrict*CS%MEKE_Cb>0.) bottomFac2 = bottomFac2 + 1./( 1. + CS%MEKE_Cb*(Ldeform/Lfrict) )**0.8
  bottomFac2 = max(bottomFac2, CS%MEKE_min_gamma)
  ! gamma_t^2 is the ratio of barotropic eddy energy to mean column eddy energy
  ! used in the velocity scale for diffusivity
  barotrFac2 = 1.
  if (Lfrict*CS%MEKE_Ct>0.) barotrFac2 = 1. / ( 1. + CS%MEKE_Ct*(Ldeform/Lfrict) )**0.25
  barotrFac2 = max(barotrFac2, CS%MEKE_min_gamma)
  if (CS%use_old_lscale) then
    if (CS%Rd_as_max_scale) then
      LmixScale = min(Ldeform, Lgrid) ! The smaller of Ld or dx
    else
      LmixScale = Lgrid
    endif
  else
    Ue = sqrt( 2.0 * max( 0., barotrFac2*EKE ) ) ! Barotropic eddy flow scale
    Lrhines = sqrt( Ue / max( beta, 1.e-30*US%T_to_s*US%L_to_m ) )       ! Rhines scale
    if (CS%aEady > 0.) then
      Leady = Ue / max( SN, 1.e-15*US%T_to_s ) ! Bound Eady time-scale < 1e15 seconds
    else
      Leady = 0.
    endif
    if (CS%use_min_lscale) then
      LmixScale = CS%lscale_maxval
      if (CS%aDeform*Ldeform > 0.) LmixScale = min(LmixScale,CS%aDeform*Ldeform)
      if (CS%aFrict *Lfrict  > 0.) LmixScale = min(LmixScale,CS%aFrict *Lfrict)
      if (CS%aRhines*Lrhines > 0.) LmixScale = min(LmixScale,CS%aRhines*Lrhines)
      if (CS%aEady  *Leady   > 0.) LmixScale = min(LmixScale,CS%aEady  *Leady)
      if (CS%aGrid  *Lgrid   > 0.) LmixScale = min(LmixScale,CS%aGrid  *Lgrid)
      if (CS%Lfixed          > 0.) LmixScale = min(LmixScale,CS%Lfixed)
    else
      LmixScale = 0.
      if (CS%aDeform*Ldeform > 0.) LmixScale = LmixScale + 1./(CS%aDeform*Ldeform)
      if (CS%aFrict *Lfrict  > 0.) LmixScale = LmixScale + 1./(CS%aFrict *Lfrict)
      if (CS%aRhines*Lrhines > 0.) LmixScale = LmixScale + 1./(CS%aRhines*Lrhines)
      if (CS%aEady  *Leady   > 0.) LmixScale = LmixScale + 1./(CS%aEady  *Leady)
      if (CS%aGrid  *Lgrid   > 0.) LmixScale = LmixScale + 1./(CS%aGrid  *Lgrid)
      if (CS%Lfixed          > 0.) LmixScale = LmixScale + 1./CS%Lfixed
      if (LmixScale > 0.) LmixScale = 1. / LmixScale
    endif
  endif

end procedure MEKE_lengthScales_0d
module procedure MEKE_init
  real :: MEKE_restoring_timescale ! The timescale used to nudge MEKE toward its equilibrium value [T ~> s]
  real :: cdrag            ! The default bottom drag coefficient [nondim].
  character(len=200) :: eke_filename, eke_varname, inputdir
  character(len=16) :: eke_source_str
  integer :: i, j, is, ie, js, je, isd, ied, jsd, jed
  logical :: laplacian, biharmonic, coldStart
# include "version_variable.h"
  character(len=40)  :: mdl = "MOM_MEKE" ! This module's name.
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed

  ! Determine whether this module will be used
  call get_param(param_file, mdl, "USE_MEKE", MEKE_init, default=.false., do_not_log=.true.)
  call log_version(param_file, mdl, version, "", all_default=.not.MEKE_init)
  call get_param(param_file, mdl, "USE_MEKE", MEKE_init, &
                 "If true, turns on the MEKE scheme which calculates "// &
                 "a sub-grid mesoscale eddy kinetic energy budget.", &
                 default=.false.)
  if (.not. MEKE_init) return
  CS%initialized = .true.
  call get_param(param_file, mdl, "MEKE_IN_DYNAMICS", meke_in_dynamics, &
                 "If true, step MEKE forward with the dynamics "// &
                 "otherwise with the tracer timestep.", &
                 default=.true.)

  call get_param(param_file, mdl, "EKE_SOURCE", eke_source_str, &
                 "Determine the where EKE comes from:\n" // &
                 "  'prog': Calculated solving EKE equation\n"// &
                 "  'file': Read in from a file\n"            // &
                 "  'dbclient': Retrieved from ML-database", default='prog')

  call MOM_mesg("MEKE_init: reading parameters ", 5)

  select case (lowercase(eke_source_str))
  case("file")
    CS%eke_src = EKE_FILE
    call time_interp_external_init
    call get_param(param_file, mdl, "EKE_FILE", eke_filename, &
                 "A file in which to find the eddy kineteic energy variable.", &
                 default="eke_file.nc")
    call get_param(param_file, mdl, "EKE_VARIABLE", eke_varname, &
                 "The name of the eddy kinetic energy variable to read from "//&
                 "EKE_FILE to use in MEKE.", &
                 default="eke")
    call get_param(param_file, mdl, "INPUTDIR", inputdir, &
                 "The directory in which all input files are found.", &
                 default=".", do_not_log=.true.)
    inputdir = slasher(inputdir)

    eke_filename = trim(inputdir) // trim(eke_filename)
    CS%eke_handle = init_external_field(eke_filename, eke_varname, domain=G%Domain%mpp_domain)
  case("prog")
    CS%eke_src = EKE_PROG
    ! Read all relevant parameters and write them to the model log.
    call get_param(param_file, mdl, "MEKE_DAMPING", CS%MEKE_damping, &
                   "The local depth-independent MEKE dissipation rate.", &
                   units="s-1", default=0.0, scale=US%T_to_s)
    call get_param(param_file, mdl, "MEKE_CD_SCALE", CS%MEKE_Cd_scale, &
                   "The ratio of the bottom eddy velocity to the column mean "//&
                   "eddy velocity, i.e. sqrt(2*MEKE). This should be less than 1 "//&
                   "to account for the surface intensification of MEKE.", &
                   units="nondim", default=0.)
    call get_param(param_file, mdl, "MEKE_CB", CS%MEKE_Cb, &
                   "A coefficient in the expression for the ratio of bottom projected "//&
                   "eddy energy and mean column energy (see Jansen et al. 2015).",&
                   units="nondim", default=25.)
    call get_param(param_file, mdl, "MEKE_MIN_GAMMA2", CS%MEKE_min_gamma, &
                   "The minimum allowed value of gamma_b^2.",&
                   units="nondim", default=0.0001)
    call get_param(param_file, mdl, "MEKE_CT", CS%MEKE_Ct, &
                   "A coefficient in the expression for the ratio of barotropic "//&
                   "eddy energy and mean column energy (see Jansen et al. 2015).",&
                   units="nondim", default=50.)
    call get_param(param_file, mdl, "MEKE_GMCOEFF", CS%MEKE_GMcoeff, &
                   "The efficiency of the conversion of potential energy "//&
                   "into MEKE by the thickness mixing parameterization. "//&
                   "If MEKE_GMCOEFF is negative, this conversion is not "//&
                   "used or calculated.", units="nondim", default=-1.0)
    call get_param(param_file, mdl, "MEKE_GEOMETRIC", CS%MEKE_GEOMETRIC, &
                   "If MEKE_GEOMETRIC is true, uses the GM coefficient formulation "//&
                   "from the GEOMETRIC framework (Marshall et al., 2012).", default=.false.)
    call get_param(param_file, mdl, "MEKE_GEOMETRIC_ALPHA", CS%MEKE_GEOMETRIC_alpha, &
                   "The nondimensional coefficient governing the efficiency of the GEOMETRIC \n"//&
                   "thickness diffusion.", units="nondim", default=0.05)
    call get_param(param_file, mdl, "MEKE_EQUILIBRIUM_ALT", CS%MEKE_equilibrium_alt, &
                   "If true, use an alternative formula for computing the (equilibrium) "//&
                   "initial value of MEKE.", default=.false.)
    call get_param(param_file, mdl, "MEKE_EQUILIBRIUM_RESTORING", CS%MEKE_equilibrium_restoring, &
                   "If true, restore MEKE back to its equilibrium value, which is calculated at "//&
                   "each time step.", default=.false.)
    if (CS%MEKE_equilibrium_restoring) then
      call get_param(param_file, mdl, "MEKE_RESTORING_TIMESCALE", MEKE_restoring_timescale, &
                     "The timescale used to nudge MEKE toward its equilibrium value.", &
                     units="s", default=1e6, scale=US%s_to_T)
      CS%MEKE_restoring_rate = 1.0 / MEKE_restoring_timescale
    endif

    call get_param(param_file, mdl, "MEKE_FRCOEFF", CS%MEKE_FrCoeff, &
                   "The efficiency of the conversion of mean energy into "//&
                   "MEKE.  If MEKE_FRCOEFF is negative, this conversion "//&
                   "is not used or calculated.", units="nondim", default=-1.0)
    call get_param(param_file, mdl, "MEKE_BHFRCOEFF", CS%MEKE_bhFrCoeff, &
                 "The efficiency of the conversion of mean energy into "//&
                 "MEKE by the biharmonic dissipation.  If MEKE_bhFRCOEFF is negative, this conversion "//&
                 "is not used or calculated.", units="nondim", default=-1.0)
    call get_param(param_file, mdl, "MEKE_GMECOEFF", CS%MEKE_GMECoeff, &
                   "The efficiency of the conversion of MEKE into mean energy "//&
                   "by GME.  If MEKE_GMECOEFF is negative, this conversion "//&
                   "is not used or calculated.", units="nondim", default=-1.0)
    call get_param(param_file, mdl, "MEKE_BGSRC", CS%MEKE_BGsrc, &
                   "A background energy source for MEKE.", &
                   units="W kg-1", default=0.0, scale=US%m_to_L**2*US%T_to_s**3)
    call get_param(param_file, mdl, "MEKE_KH", CS%MEKE_Kh, &
                   "A background lateral diffusivity of MEKE. "//&
                   "Use a negative value to not apply lateral diffusion to MEKE.", &
                   units="m2 s-1", default=-1.0, scale=US%m_to_L**2*US%T_to_s)
    call get_param(param_file, mdl, "MEKE_K4", CS%MEKE_K4, &
                   "A lateral bi-harmonic diffusivity of MEKE. "//&
                   "Use a negative value to not apply bi-harmonic diffusion to MEKE.", &
                   units="m4 s-1", default=-1.0, scale=US%m_to_L**4*US%T_to_s)
    call get_param(param_file, mdl, "MEKE_DTSCALE", CS%MEKE_dtScale, &
                   "A scaling factor to accelerate the time evolution of MEKE.", &
                   units="nondim", default=1.0)
    call get_param(param_file, mdl, "MEKE_POSITIVE", CS%MEKE_positive, &
                   "If true, it guarantees that MEKE will always be >= 0.", &
                   default=.false.)
  case("dbclient")
    CS%eke_src = EKE_DBCLIENT
    call ML_MEKE_init(diag, G, US, Time, param_file, dbcomms_CS, CS)
  case default
    call MOM_error(FATAL, "Invalid method selected for calculating EKE")
  end select
  ! GMM, make sure all parameters used to calculated MEKE are within the above if

  call get_param(param_file, mdl, "MEKE_KHCOEFF", CS%MEKE_KhCoeff, &
                 "A scaling factor in the expression for eddy diffusivity "//&
                 "which is otherwise proportional to the MEKE velocity- "//&
                 "scale times an eddy mixing-length. This factor "//&
                 "must be >0 for MEKE to contribute to the thickness/ "//&
                 "and tracer diffusivity in the rest of the model.", &
                 units="nondim", default=1.0)
  call get_param(param_file, mdl, "MEKE_USCALE", CS%MEKE_Uscale, &
                 "The background velocity that is combined with MEKE to "//&
                 "calculate the bottom drag.", units="m s-1", default=0.0, scale=US%m_s_to_L_T)
  call get_param(param_file, mdl, "MEKE_GM_SRC_ALT", CS%GM_src_alt, &
                 "If true, use the GM energy conversion form S^2*N^2*kappa rather "//&
                 "than the streamfunction for the MEKE GM source term.", default=.false.)
  call get_param(param_file, mdl, "MEKE_MIN_DEPTH_TOT", CS%MEKE_min_depth_tot, &
                 "The minimum total depth over which to distribute MEKE energy sources.  "//&
                 "When the total depth is less than this, the sources are scaled away.", &
                 units="m", default=1.0, scale=GV%m_to_H, do_not_log=.not.CS%GM_src_alt)
  call get_param(param_file, mdl, "MEKE_VISC_DRAG", CS%visc_drag, &
                 "If true, use the vertvisc_type to calculate the bottom "//&
                 "drag acting on MEKE.", default=.true.)
  call get_param(param_file, mdl, "MEKE_KHTH_FAC", MEKE%KhTh_fac, &
                 "A factor that maps MEKE%Kh to KhTh.", units="nondim", default=0.0)
  call get_param(param_file, mdl, "MEKE_KHTR_FAC", MEKE%KhTr_fac, &
                 "A factor that maps MEKE%Kh to KhTr.", units="nondim", default=0.0)
  call get_param(param_file, mdl, "MEKE_KHMEKE_FAC", CS%KhMEKE_Fac, &
                 "A factor that maps MEKE%Kh to Kh for MEKE itself.", &
                 units="nondim", default=0.0)
  call get_param(param_file, mdl, "MEKE_OLD_LSCALE", CS%use_old_lscale, &
                 "If true, use the old formula for length scale which is "//&
                 "a function of grid spacing and deformation radius.",  &
                 default=.false.)
  call get_param(param_file, mdl, "MEKE_MIN_LSCALE", CS%use_min_lscale, &
                 "If true, use a strict minimum of provided length scales "//&
                 "rather than harmonic mean.",  &
                 default=.false.)
  call get_param(param_file, mdl, "MEKE_LSCALE_MAX_VAL", CS%lscale_maxval, &
                 "The ceiling on the value of the MEKE length scale when MEKE_MIN_LSCALE=True.  "//&
                 "The default is the distance from the equator to the pole on Earth, as "//&
                 "estimated by enlightenment era scientists, but should probably scale with RAD_EARTH.", &
                 units="m", default=1.0e7, scale=US%m_to_L, do_not_log=.not.CS%use_min_lscale)
  call get_param(param_file, mdl, "MEKE_RD_MAX_SCALE", CS%Rd_as_max_scale, &
                 "If true, the length scale used by MEKE is the minimum of "//&
                 "the deformation radius or grid-spacing. Only used if "//&
                 "MEKE_OLD_LSCALE=True", default=.false.)
  call get_param(param_file, mdl, "MEKE_VISCOSITY_COEFF_KU", CS%viscosity_coeff_Ku, &
                 "If non-zero, is the scaling coefficient in the expression for "//&
                 "viscosity used to parameterize harmonic lateral momentum mixing by "//&
                 "unresolved eddies represented by MEKE. Can be negative to "//&
                 "represent backscatter from the unresolved eddies.", &
                 units="nondim", default=0.0)
  call get_param(param_file, mdl, "MEKE_VISCOSITY_COEFF_AU", CS%viscosity_coeff_Au, &
                 "If non-zero, is the scaling coefficient in the expression for "//&
                 "viscosity used to parameterize biharmonic lateral momentum mixing by "//&
                 "unresolved eddies represented by MEKE. Can be negative to "//&
                 "represent backscatter from the unresolved eddies.", &
                 units="nondim", default=0.0)
  call get_param(param_file, mdl, "MEKE_FIXED_MIXING_LENGTH", CS%Lfixed, &
                 "If positive, is a fixed length contribution to the expression "//&
                 "for mixing length used in MEKE-derived diffusivity.", &
                 units="m", default=0.0, scale=US%m_to_L)
  call get_param(param_file, mdl, "MEKE_FIXED_TOTAL_DEPTH", CS%fixed_total_depth, &
                 "If true, use the nominal bathymetric depth as the estimate of the "//&
                 "time-varying ocean depth.  Otherwise base the depth on the total ocean mass "//&
                 "per unit area.", default=.true.)
  call get_param(param_file, mdl, "MEKE_TOTAL_DEPTH_RHO", CS%rho_fixed_total_depth, &
                 "A density used to translate the nominal bathymetric depth into an estimate "//&
                 "of the total ocean mass per unit area when MEKE_FIXED_TOTAL_DEPTH is true.", &
                 units="kg m-3", default=GV%Rho0*US%R_to_kg_m3, scale=US%kg_m3_to_R, &
                 do_not_log=(GV%Boussinesq.or.(.not.CS%fixed_total_depth)))

  call get_param(param_file, mdl, "MEKE_ALPHA_DEFORM", CS%aDeform, &
                 "If positive, is a coefficient weighting the deformation scale "//&
                 "in the expression for mixing length used in MEKE-derived diffusivity.", &
                 units="nondim", default=0.0)
  call get_param(param_file, mdl, "MEKE_ALPHA_RHINES", CS%aRhines, &
                 "If positive, is a coefficient weighting the Rhines scale "//&
                 "in the expression for mixing length used in MEKE-derived diffusivity.", &
                 units="nondim", default=0.0)
  call get_param(param_file, mdl, "MEKE_ALPHA_EADY", CS%aEady, &
                 "If positive, is a coefficient weighting the Eady length scale "//&
                 "in the expression for mixing length used in MEKE-derived diffusivity.", &
                 units="nondim", default=0.0)
  call get_param(param_file, mdl, "MEKE_ALPHA_FRICT", CS%aFrict, &
                 "If positive, is a coefficient weighting the frictional arrest scale "//&
                 "in the expression for mixing length used in MEKE-derived diffusivity.", &
                 units="nondim", default=0.0)
  call get_param(param_file, mdl, "MEKE_ALPHA_GRID", CS%aGrid, &
                 "If positive, is a coefficient weighting the grid-spacing as a scale "//&
                 "in the expression for mixing length used in MEKE-derived diffusivity.", &
                 units="nondim", default=0.0)
  call get_param(param_file, mdl, "MEKE_COLD_START", coldStart, &
                 "If true, initialize EKE to zero. Otherwise a local equilibrium solution "//&
                 "is used as an initial condition for EKE.", default=.false.)
  call get_param(param_file, mdl, "MEKE_BACKSCAT_RO_C", MEKE%backscatter_Ro_c, &
                 "The coefficient in the Rossby number function for scaling the biharmonic "//&
                 "frictional energy source. Setting to non-zero enables the Rossby number function.", &
                 units="nondim", default=0.0)
  call get_param(param_file, mdl, "MEKE_BACKSCAT_RO_POW", MEKE%backscatter_Ro_pow, &
                 "The power in the Rossby number function for scaling the biharmonic "//&
                 "frictional energy source.", units="nondim", default=0.0)
  call get_param(param_file, mdl, "MEKE_ADVECTION_FACTOR", CS%MEKE_advection_factor, &
                 "A scale factor in front of advection of eddy energy. Zero turns advection off. "//&
                 "Using unity would be normal but other values could accommodate a mismatch "//&
                 "between the advecting barotropic flow and the vertical structure of MEKE.", &
                 units="nondim", default=0.0)
  call get_param(param_file, mdl, "MEKE_ADVECTION_BUG", CS%MEKE_advection_bug, &
                 "If true, recover a bug in the calculation of the barotropic transport for "//&
                 "the advection of MEKE.  With the bug, only the transports in the deepest "//&
                 "layer are used.", default=.false., do_not_log=(CS%MEKE_advection_factor<=0.))
  call get_param(param_file, mdl, "MEKE_TOPOGRAPHIC_BETA", CS%MEKE_topographic_beta, &
                 "A scale factor to determine how much topographic beta is weighed in " //&
                 "computing beta in the expression of Rhines scale. Use 1 if full "//&
                 "topographic beta effect is considered; use 0 if it's completely ignored.", &
                 units="nondim", default=0.0)
  call get_param(param_file, mdl, "SQG_USE_MEKE", CS%sqg_use_MEKE, &
                 "If true, the eddy scale of MEKE is used for the SQG vertical structure ",&
                 default=.false.)

  ! Nonlocal module parameters
  call get_param(param_file, mdl, "CDRAG", cdrag, &
                 "CDRAG is the drag coefficient relating the magnitude of the velocity "//&
                 "field to the bottom stress.", units="nondim", default=0.003)
  call get_param(param_file, mdl, "MEKE_CDRAG", CS%cdrag, &
                 "Drag coefficient relating the magnitude of the velocity "//&
                 "field to the bottom stress in MEKE.", units="nondim", default=cdrag, scale=US%L_to_m*GV%m_to_H)
  call get_param(param_file, mdl, "LAPLACIAN", laplacian, default=.false., do_not_log=.true.)
  call get_param(param_file, mdl, "BIHARMONIC", biharmonic, default=.false., do_not_log=.true.)

  if (CS%viscosity_coeff_Ku/=0. .and. .not. laplacian) call MOM_error(FATAL, &
                 "LAPLACIAN must be true if MEKE_VISCOSITY_COEFF_KU is true.")

  if (CS%viscosity_coeff_Au/=0. .and. .not. biharmonic) call MOM_error(FATAL, &
                 "BIHARMONIC must be true if MEKE_VISCOSITY_COEFF_AU is true.")

  call get_param(param_file, mdl, "DEBUG", CS%debug, default=.false., do_not_log=.true.)

  ! Identify if any lateral diffusive processes are active
  CS%kh_flux_enabled = .false.
  if ((CS%MEKE_KH >= 0.0)  .or. (CS%KhMEKE_FAC > 0.0) .or. (CS%MEKE_advection_factor > 0.0)) &
    CS%kh_flux_enabled = .true.

! Register fields for output from this module.
  CS%diag => diag
  CS%id_MEKE = register_diag_field('ocean_model', 'MEKE', diag%axesT1, Time, &
     'Mesoscale Eddy Kinetic Energy', 'm2 s-2', conversion=US%L_T_to_m_s**2)
  if (.not. allocated(MEKE%MEKE)) CS%id_MEKE = -1
  CS%id_Kh = register_diag_field('ocean_model', 'MEKE_KH', diag%axesT1, Time, &
     'MEKE derived diffusivity', 'm2 s-1', conversion=US%L_to_m**2*US%s_to_T)
  if (.not. allocated(MEKE%Kh)) CS%id_Kh = -1
  CS%id_Ku = register_diag_field('ocean_model', 'MEKE_KU', diag%axesT1, Time, &
     'MEKE derived lateral viscosity', 'm2 s-1', conversion=US%L_to_m**2*US%s_to_T)
  if (.not. allocated(MEKE%Ku)) CS%id_Ku = -1
  CS%id_Au = register_diag_field('ocean_model', 'MEKE_AU', diag%axesT1, Time, &
     'MEKE derived lateral biharmonic viscosity', 'm4 s-1', conversion=US%L_to_m**4*US%s_to_T)
  if (.not. allocated(MEKE%Au)) CS%id_Au = -1
  CS%id_Ue = register_diag_field('ocean_model', 'MEKE_Ue', diag%axesT1, Time, &
     'MEKE derived eddy-velocity scale', 'm s-1', conversion=US%L_T_to_m_s)
  if (.not. allocated(MEKE%MEKE)) CS%id_Ue = -1
  CS%id_Ub = register_diag_field('ocean_model', 'MEKE_Ub', diag%axesT1, Time, &
     'MEKE derived bottom eddy-velocity scale', 'm s-1', conversion=US%L_T_to_m_s)
  if (.not. allocated(MEKE%MEKE)) CS%id_Ub = -1
  CS%id_Ut = register_diag_field('ocean_model', 'MEKE_Ut', diag%axesT1, Time, &
     'MEKE derived barotropic eddy-velocity scale', 'm s-1', conversion=US%L_T_to_m_s)
  if (.not. allocated(MEKE%MEKE)) CS%id_Ut = -1
  CS%id_src = register_diag_field('ocean_model', 'MEKE_src', diag%axesT1, Time, &
     'MEKE energy source', 'm2 s-3', conversion=(US%L_T_to_m_s**2)*US%s_to_T)

  CS%id_src_adv = register_diag_field('ocean_model', 'MEKE_src_adv', diag%axesT1, Time, &
     'MEKE energy source from the horizontal advection of MEKE', 'm2 s-3', conversion=(US%L_T_to_m_s**2)*US%s_to_T)

  CS%id_src_btm_drag = register_diag_field('ocean_model', 'MEKE_src_btm_drag', diag%axesT1, Time, &
     'MEKE energy source from the bottom drag acting on MEKE', 'm2 s-3', conversion=(US%L_T_to_m_s**2)*US%s_to_T)

  if (CS%MEKE_K4 >= 0.) &
    CS%id_src_mom_K4 = register_diag_field('ocean_model', 'MEKE_src_mom_K4', &
        diag%axesT1, Time, 'MEKE energy source from the biharmonic of MEKE', &
        'm2 s-3', conversion=(US%L_T_to_m_s**2)*US%s_to_T)

  if (CS%MEKE_GMcoeff >= 0.) &
    CS%id_src_GM = register_diag_field('ocean_model', 'MEKE_src_GM', &
        diag%axesT1, Time, 'MEKE energy source from the thickness mixing (GM scheme)', &
        'm2 s-3', conversion=(US%L_T_to_m_s**2)*US%s_to_T)

  if (CS%MEKE_FrCoeff >= 0.) &
    CS%id_src_mom_lp = register_diag_field('ocean_model', 'MEKE_src_mom_lp', &
        diag%axesT1, Time, 'MEKE energy source from the Laplacian of resolved flows', &
        'm2 s-3', conversion=(US%L_T_to_m_s**2)*US%s_to_T)

  if (CS%MEKE_bhFrCoeff >= 0.) &
    CS%id_src_mom_bh = register_diag_field('ocean_model', 'MEKE_src_mom_bh', &
        diag%axesT1, Time, 'MEKE energy source from the biharmonic of resolved flows', &
        'm2 s-3', conversion=(US%L_T_to_m_s**2)*US%s_to_T)

  CS%id_decay = register_diag_field('ocean_model', 'MEKE_decay', diag%axesT1, Time, &
     'MEKE decay rate', 's-1', conversion=US%s_to_T)
  CS%id_GM_src = register_diag_field('ocean_model', 'MEKE_GM_src', diag%axesT1, Time, &
     'MEKE energy available from thickness mixing', &
     'W m-2', conversion=US%RZ3_T3_to_W_m2*US%L_to_Z**2)
  if (.not. allocated(MEKE%GM_src)) CS%id_GM_src = -1
  CS%id_mom_src = register_diag_field('ocean_model', 'MEKE_mom_src',diag%axesT1, Time, &
     'MEKE energy available from momentum', &
     'W m-2', conversion=US%RZ3_T3_to_W_m2*US%L_to_Z**2)
  if (.not. allocated(MEKE%mom_src)) CS%id_mom_src = -1
  CS%id_mom_src_bh = register_diag_field('ocean_model', 'MEKE_mom_src_bh',diag%axesT1, Time, &
     'MEKE energy available from the biharmonic dissipation of momentum', &
     'W m-2', conversion=US%RZ3_T3_to_W_m2*US%L_to_Z**2)
  if (.not. allocated(MEKE%mom_src_bh)) CS%id_mom_src_bh = -1
  CS%id_GME_snk = register_diag_field('ocean_model', 'MEKE_GME_snk',diag%axesT1, Time, &
     'MEKE energy lost to GME backscatter', &
     'W m-2', conversion=US%RZ3_T3_to_W_m2*US%L_to_Z**2)
  if (.not. allocated(MEKE%GME_snk)) CS%id_GME_snk = -1
  CS%id_Le = register_diag_field('ocean_model', 'MEKE_Le', diag%axesT1, Time, &
     'Eddy mixing length used in the MEKE derived eddy diffusivity', 'm', conversion=US%L_to_m)
  CS%id_Lrhines = register_diag_field('ocean_model', 'MEKE_Lrhines', diag%axesT1, Time, &
     'Rhines length scale used in the MEKE derived eddy diffusivity', 'm', conversion=US%L_to_m)
  CS%id_Leady = register_diag_field('ocean_model', 'MEKE_Leady', diag%axesT1, Time, &
     'Eady length scale used in the MEKE derived eddy diffusivity', 'm', conversion=US%L_to_m)
  CS%id_gamma_b = register_diag_field('ocean_model', 'MEKE_gamma_b', diag%axesT1, Time, &
     'Ratio of bottom-projected eddy velocity to column-mean eddy velocity', 'nondim')
  CS%id_gamma_t = register_diag_field('ocean_model', 'MEKE_gamma_t', diag%axesT1, Time, &
     'Ratio of barotropic eddy velocity to column-mean eddy velocity', 'nondim')

  if (CS%kh_flux_enabled) then
    CS%id_KhMEKE_u = register_diag_field('ocean_model', 'KHMEKE_u', diag%axesCu1, Time, &
     'Zonal diffusivity of MEKE', 'm2 s-1', conversion=US%L_to_m**2*US%s_to_T)
    CS%id_KhMEKE_v = register_diag_field('ocean_model', 'KHMEKE_v', diag%axesCv1, Time, &
     'Meridional diffusivity of MEKE', 'm2 s-1', conversion=US%L_to_m**2*US%s_to_T)
  endif

  if (CS%MEKE_equilibrium_restoring) then
    CS%id_MEKE_equilibrium = register_diag_field('ocean_model', 'MEKE_equilibrium', diag%axesT1, Time, &
     'Equilibrated Mesoscale Eddy Kinetic Energy', 'm2 s-2', conversion=US%L_T_to_m_s**2)
  endif

  CS%id_clock_pass = cpu_clock_id('(Ocean continuity halo updates)', grain=CLOCK_ROUTINE)


  ! Detect whether this instance of MEKE_init() is at the beginning of a run
  ! or after a restart. If at the beginning, we will initialize MEKE to a local
  ! equilibrium.
  CS%initialize = .not.query_initialized(MEKE%MEKE, "MEKE", restart_CS)
  if (coldStart) CS%initialize = .false.
  if (CS%initialize) call MOM_error(WARNING, &
                       "MEKE_init: Initializing MEKE with a local equilibrium balance.")
  if (allocated(MEKE%Le)) then
    if (.not.query_initialized(MEKE%Le, "MEKE_Le", restart_CS)) then
      !$OMP parallel do default(shared)
      do j=js,je ; do i=is,ie
        MEKE%Le(i,j) = sqrt(G%areaT(i,j))
      enddo ; enddo
    endif
  endif

  ! Set up group passes.  In the case of a restart, these fields need a halo update now.
  if (allocated(MEKE%MEKE)) then
    call create_group_pass(CS%pass_MEKE, MEKE%MEKE, G%Domain)
    if (allocated(MEKE%Kh_diff)) call create_group_pass(CS%pass_MEKE, MEKE%Kh_diff, G%Domain)
    if (.not.CS%initialize) call do_group_pass(CS%pass_MEKE, G%Domain)
  endif
  if (allocated(MEKE%Kh)) call create_group_pass(CS%pass_Kh, MEKE%Kh, G%Domain)
  if (allocated(MEKE%Ku)) call create_group_pass(CS%pass_Kh, MEKE%Ku, G%Domain)
  if (allocated(MEKE%Au)) call create_group_pass(CS%pass_Kh, MEKE%Au, G%Domain)
  if (allocated(MEKE%Le)) call create_group_pass(CS%pass_Kh, MEKE%Le, G%Domain)

  if (allocated(MEKE%Kh) .or. allocated(MEKE%Ku) .or. allocated(MEKE%Au) &
      .or. allocated(MEKE%Le)) &
    call do_group_pass(CS%pass_Kh, G%Domain)

end procedure MEKE_init
module procedure ML_MEKE_init
  character(len=200)  :: inputdir, backend, model_filename
  integer :: db_return_code, batch_size
  character(len=40) :: mdl = "MOM_ML_MEKE"
  write(CS%key_suffix, '(A,I6.6)') '_', PE_here()
  ! Put some basic information into the database
  db_return_code = 0
  db_return_code = CS%client%put_tensor("meta"//CS%key_suffix, &
    REAL([G%isd_global, G%idg_offset, G%jsd_global, G%jdg_offset]),[4]) + db_return_code
  db_return_code = CS%client%put_tensor("geolat"//CS%key_suffix, G%geoLatT, shape(G%geoLatT)) + db_return_code
  db_return_code = CS%client%put_tensor("geolon"//CS%key_suffix, G%geoLonT, shape(G%geoLonT)) + db_return_code
  db_return_code = CS%client%put_tensor("EKE_shape"//CS%key_suffix, shape(G%geolonT), [2]) + db_return_code

  if (CS%client%SR_error_parser(db_return_code)) call MOM_error(FATAL, "Putting metadata into the database failed")

  call read_param(param_file, "INPUTDIR", inputdir)
  inputdir = slasher(inputdir)

  call get_param(param_file, mdl, "BATCH_SIZE", batch_size, "Batch size to use for inference", default=1)
  call get_param(param_file, mdl, "EKE_BACKEND", backend, &
                 "The computational backend to use for EKE inference (CPU or GPU)", default="GPU")
  call get_param(param_file, mdl, "EKE_MODEL", model_filename, &
                 "Filename of the a saved pyTorch model to use", fail_if_missing = .true.)
  call get_param(param_file, mdl, "EKE_MAX", CS%eke_max, &
                 "Maximum value of EKE allowed when inferring EKE", &
                 units="m2 s-2", default=2., scale=US%m_s_to_L_T**2)

  ! Set the machine learning model
  if (dbcomms_CS%colocated) then
    if (modulo(PE_here(),dbcomms_CS%colocated_stride) == 0) then
      db_return_code = CS%client%set_model_from_file(CS%model_key, trim(inputdir)//trim(model_filename), &
                                                  "TORCH", backend, batch_size=batch_size)
    endif
  else
    if (is_root_pe()) then
      db_return_code = CS%client%set_model_from_file(CS%model_key, trim(inputdir)//trim(model_filename), &
                                                  "TORCH", backend, batch_size=batch_size)
    endif
  endif
  if (CS%client%SR_error_parser(db_return_code)) then
    call MOM_error(FATAL, "MEKE: set_model failed")
  endif

  call get_param(param_file, mdl, "ONLINE_ANALYSIS", CS%online_analysis, &
               "If true, post EKE used in MOM6 to the database for analysis", default=.true.)

  ! Set various clock ids
  CS%id_client_init   = cpu_clock_id('(ML_MEKE client init)', grain=CLOCK_ROUTINE)
  CS%id_put_tensor    = cpu_clock_id('(ML_MEKE put tensor)', grain=CLOCK_ROUTINE)
  CS%id_run_model     = cpu_clock_id('(ML_MEKE run model)', grain=CLOCK_ROUTINE)
  CS%id_unpack_tensor = cpu_clock_id('(ML_MEKE unpack tensor )', grain=CLOCK_ROUTINE)

  ! Diagnostics for ML_MEKE
  CS%id_mke = register_diag_field('ocean_model', 'MEKE_MKE', diag%axesT1, Time, &
     'Surface mean (resolved) kinetic energy used in MEKE', 'm2 s-2', conversion=US%L_T_to_m_s**2)
  CS%id_slope_z= register_diag_field('ocean_model', 'MEKE_slope_z', diag%axesT1, Time, &
     'Vertically averaged isopyncal slope magnitude used in MEKE', 'nondim', conversion=US%Z_to_L)
  CS%id_slope_x= register_diag_field('ocean_model', 'MEKE_slope_x', diag%axesCui, Time, &
     'Isopycnal slope in the x-direction used in MEKE', 'nondim', conversion=US%Z_to_L)
  CS%id_slope_y= register_diag_field('ocean_model', 'MEKE_slope_y', diag%axesCvi, Time, &
     'Isopycnal slope in the y-direction used in MEKE', 'nondim', conversion=US%Z_to_L)
  CS%id_rv = register_diag_field('ocean_model', 'MEKE_RV', diag%axesT1, Time, &
     'Surface relative vorticity used in MEKE', 's-1', conversion=US%s_to_T)

end procedure ML_MEKE_init
module procedure ML_MEKE_calculate_features
  real, dimension(SZI_(G),SZJ_(G)) :: mke      ! Surface kinetic energy per unit mass [L2 T-2 ~> m2 s-2]
  real, dimension(SZI_(G),SZJ_(G)) :: slope_z  ! Vertically averaged isoneutral slopes [Z L-1 ~> nondim]
  real, dimension(SZIB_(G),SZJB_(G)) :: rv_z   ! Surface relative vorticity [T-1 ~> s-1]
  real, dimension(SZIB_(G),SZJB_(G)) :: rv_z_t ! Surface relative vorticity interpolated to tracer points [T-1 ~> s-1]
  real, dimension(SZIB_(G),SZJ_(G), SZK_(G)) :: h_u ! Thickness at u point [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJB_(G), SZK_(G)) :: h_v ! Thickness at v point [H ~> m or kg m-2]
  real, dimension(SZIB_(G),SZJ_(G),SZK_(G)+1) :: slope_x ! Isoneutral slope at U point [Z L-1 ~> nondim]
  real, dimension(SZI_(G),SZJB_(G),SZK_(G)+1) :: slope_y ! Isoneutral slope at V point [Z L-1 ~> nondim]
  real, dimension(SZIB_(G),SZJ_(G)) :: slope_x_vert_avg ! Isoneutral slope at U point [Z L-1 ~> nondim]
  real, dimension(SZI_(G),SZJB_(G)) :: slope_y_vert_avg ! Isoneutral slope at V point [Z L-1 ~> nondim]
  real, dimension(SZI_(G), SZJ_(G), SZK_(G)+1) ::  e ! The interface heights relative to mean sea level [Z ~> m].
  real :: slope_t  ! Slope interpolated to thickness points [Z L-1 ~> nondim]
  real :: u_t, v_t ! u and v interpolated to thickness points [L T-1 ~> m s-1]
  real :: dvdx, dudy ! Components of relative vorticity [T-1 ~> s-1]
  real :: a_e, a_w, a_n, a_s ! Fractional areas of neighboring cells for interpolating velocities [nondim]
  real :: Idenom    ! A normalizing factor in calculating weighted averages of areas [L-2 ~> m-2]
  real :: sum_area  ! A sum of adjacent cell areas [L2 ~> m2]
  integer :: i, j, k, is, ie, js, je, Isq, Ieq, Jsq, Jeq, nz
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke
  Isq = G%IscB ; Ieq = G%IecB ; Jsq = G%JscB ; Jeq = G%JecB

  ! Calculate various features for used to infer eddy kinetic energy
  ! Linear interpolation to estimate thickness at a velocity points
  do k=1,nz ; do j=js-1,je+1 ; do i=is-1,ie+1
    h_u(I,j,k) = 0.5*(h(i,j,k)*G%mask2dT(i,j) + h(i+1,j,k)*G%mask2dT(i+1,j)) + GV%Angstrom_H
    h_v(i,J,k) = 0.5*(h(i,j,k)*G%mask2dT(i,j) + h(i,j+1,k)*G%mask2dT(i,j+1)) + GV%Angstrom_H
  enddo ; enddo ; enddo
  call find_eta(h, tv, G, GV, US, e, halo_size=2)
  ! Note the hard-coded dimenisional constant in the following line.
  call calc_isoneutral_slopes(G, GV, US, h, e, tv, dt*1.e-7*GV%m2_s_to_HZ_T, .false., slope_x, slope_y)
  call pass_vector(slope_x, slope_y, G%Domain)
  do j=js-1,je+1 ; do i=is-1,ie+1
    slope_x_vert_avg(I,j) = vertical_average_interface(slope_x(i,j,:), h_u(i,j,:), GV%H_subroundoff)
    slope_y_vert_avg(i,J) = vertical_average_interface(slope_y(i,j,:), h_v(i,j,:), GV%H_subroundoff)
  enddo ; enddo
  slope_z(:,:) = 0.

  call pass_vector(slope_x_vert_avg, slope_y_vert_avg, G%Domain)
  do j=js,je ; do i=is,ie
    ! Calculate weights for interpolation from velocity points to h points
    sum_area = G%areaCu(I-1,j) + G%areaCu(I,j)
    if (sum_area>0.0) then
      Idenom = sqrt(0.5*G%IareaT(i,j) / sum_area)
      a_w = G%areaCu(I-1,j) * Idenom
      a_e = G%areaCu(I,j) * Idenom
    else
      a_w = 0.0 ; a_e = 0.0
    endif

    sum_area = G%areaCv(i,J-1) + G%areaCv(i,J)
    if (sum_area>0.0) then
      Idenom = sqrt(0.5*G%IareaT(i,j) / sum_area)
      a_s = G%areaCv(i,J-1) * Idenom
      a_n = G%areaCv(i,J) * Idenom
    else
      a_s = 0.0 ; a_n = 0.0
    endif

    ! Calculate mean kinetic energy
    u_t = (a_e*u(I,j,1)) + (a_w*u(I-1,j,1))
    v_t = (a_n*v(i,J,1)) + (a_s*v(i,J-1,1))
    mke(i,j) = 0.5*( (u_t*u_t) + (v_t*v_t) )

    ! Calculate the magnitude of the slope
    slope_t = slope_x_vert_avg(I,j)*a_e+slope_x_vert_avg(I-1,j)*a_w
    slope_z(i,j) = sqrt(slope_t*slope_t)
    slope_t = slope_y_vert_avg(i,J)*a_n+slope_y_vert_avg(i,J-1)*a_s
    slope_z(i,j) = 0.5*(slope_z(i,j) + sqrt(slope_t*slope_t))*G%mask2dT(i,j)
  enddo ; enddo
  call pass_var(slope_z, G%Domain)

  ! Calculate relative vorticity
  do J=Jsq-1,Jeq+1 ; do I=Isq-1,Ieq+1
    dvdx = ((v(i+1,J,1)*G%dyCv(i+1,J)) - (v(i,J,1)*G%dyCv(i,J)))
    dudy = ((u(I,j+1,1)*G%dxCu(I,j+1)) - (u(I,j,1)*G%dxCu(I,j)))
    ! Assumed no slip
    rv_z(I,J) = (2.0-G%mask2dBu(I,J)) * (dvdx - dudy) * G%IareaBu(I,J)
  enddo ; enddo
  ! Interpolate RV to t-point, revisit this calculation to include metrics
  do j=js,je ; do i=is,ie
    rv_z_t(i,j) = 0.25*(rv_z(i-1,j) + rv_z(i,j) + rv_z(i-1,j-1) + rv_z(i,j-1))
  enddo ; enddo


  ! Construct the feature array
  features_array(:,mke_idx) = pack(mke,.true.)
  features_array(:,slope_z_idx) = pack(slope_z,.true.)
  features_array(:,rd_dx_z_idx) = pack(Rd_dx_h,.true.)
  features_array(:,rv_idx) = pack(rv_z_t,.true.)

  if (CS%id_rv>0) call post_data(CS%id_rv, rv_z, CS%diag)
  if (CS%id_mke>0) call post_data(CS%id_mke, mke, CS%diag)
  if (CS%id_slope_z>0) call post_data(CS%id_slope_z, slope_z, CS%diag)
  if (CS%id_slope_x>0) call post_data(CS%id_slope_x, slope_x, CS%diag)
  if (CS%id_slope_y>0) call post_data(CS%id_slope_y, slope_y, CS%diag)
end procedure ML_MEKE_calculate_features
module procedure predict_MEKE
  integer :: db_return_code
  character(len=255), dimension(1) :: model_out, model_in
  character(len=255) :: time_suffix
  real(kind=real32), dimension(SIZE(MEKE)) :: MEKE_vec ! A one-dimensional array of the natural log of eddy kinetic
  real, dimension(size(MEKE,1),size(MEKE,2)) :: ln_MEKE  ! the natural log of eddy kinetic energy
  real, dimension(size(MEKE,1),size(MEKE,2)) :: MEKE_mks  ! The eddy kinetic energy in mks units [m2 s-2]
  integer :: i, j, is, ie, js, je
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec
!> Use the database client to call a machine learning model to predict eddy kinetic energy
  call cpu_clock_begin(CS%id_put_tensor)
  db_return_code = CS%client%put_tensor("features"//CS%key_suffix, features_array, shape(features_array))
  call cpu_clock_end(CS%id_put_tensor)

  ! Run the ML model to predict EKE and return the result
  model_out(1) = "EKE"//CS%key_suffix
  model_in(1) = "features"//CS%key_suffix
  call cpu_clock_begin(CS%id_run_model)
  db_return_code = CS%client%run_model(CS%model_key, model_in, model_out)
  call cpu_clock_end(CS%id_run_model)
  if (CS%client%SR_error_parser(db_return_code)) then
    call MOM_error(FATAL, "MEKE: run_model failed")
  endif
  call cpu_clock_begin(CS%id_unpack_tensor)
  db_return_code = CS%client%unpack_tensor( model_out(1), MEKE_vec, shape(MEKE_vec) )
  call cpu_clock_end(CS%id_unpack_tensor)

  ln_MEKE = reshape(MEKE_vec, shape(MEKE))
  ! Zero out the halos.  These will usually be reset by the pass_var in a few lines.
  MEKE_mks(:,:) = 0.0
  do j=js,je ; do i=is,ie
    MEKE_mks(i,j) = MIN(exp(ln_MEKE(i,j)), US%L_T_to_m_s**2*CS%eke_max)
  enddo ; enddo
  call pass_var(MEKE_mks, G%Domain, halo=1)

  if (CS%online_analysis) then
    write(time_suffix,"(F16.0)") time_type_to_real(Time)
    db_return_code = CS%client%put_tensor(trim("EKE_")//trim(adjustl(time_suffix))//CS%key_suffix, &
                                          MEKE_mks, shape(MEKE))
  endif

  ! Copy MEKE_mks into the argument in rescaled units.
  ! MEKE(:,:) = 0.0  ! This would fill in the wider halos of this intent(out) array.
  do j=js-1,je+1 ; do i=is-1,ie+1
    MEKE(i,j) = US%m_s_to_L_T**2 * MEKE_mks(i,j)
  enddo ; enddo

end procedure predict_MEKE
module procedure vertical_average_interface
  real :: htot  ! Twice the sum of the layer thicknesses interpolated to interior interfaces [H ~> m or kg m-2]
  real :: inv_htot ! The inverse of htot  [H-1 ~> m-1 or m2 kg-1]
  integer :: k, nk
  nk = size(h)
  htot = h_min
  do k=2,nk
    htot = htot + (h(k-1)+h(k))
  enddo
  inv_htot = 1./htot

  vertical_average_interface = 0.
  do K=2,nk
    vertical_average_interface = vertical_average_interface + (w(k)*(h(k-1)+h(k)))*inv_htot
  enddo
end procedure vertical_average_interface
module procedure MEKE_alloc_register_restart
  real :: MEKE_GMcoeff, MEKE_FrCoeff, MEKE_bhFrCoeff, MEKE_GMECoeff  ! Coefficients for various terms [nondim]
  real :: MEKE_KHCoeff, MEKE_viscCoeff_Ku, MEKE_viscCoeff_Au  ! Coefficients for various terms [nondim]
  logical :: Use_KH_in_MEKE
  logical :: useMEKE
  logical :: sqg_use_MEKE
  integer :: isd, ied, jsd, jed
  useMEKE = .false. ; call read_param(param_file,"USE_MEKE",useMEKE)

! Read these parameters to determine what should be in the restarts
  MEKE_GMcoeff = -1. ; call read_param(param_file,"MEKE_GMCOEFF",MEKE_GMcoeff)
  MEKE_FrCoeff = -1. ; call read_param(param_file,"MEKE_FRCOEFF",MEKE_FrCoeff)
  MEKE_bhFrCoeff = -1. ; call read_param(param_file,"MEKE_bhFRCOEFF",MEKE_bhFrCoeff)
  MEKE_GMEcoeff = -1. ; call read_param(param_file,"MEKE_GMECOEFF",MEKE_GMEcoeff)
  MEKE_KhCoeff = 1. ; call read_param(param_file,"MEKE_KHCOEFF",MEKE_KhCoeff)
  MEKE_viscCoeff_Ku = 0. ; call read_param(param_file,"MEKE_VISCOSITY_COEFF_KU",MEKE_viscCoeff_Ku)
  MEKE_viscCoeff_Au = 0. ; call read_param(param_file,"MEKE_VISCOSITY_COEFF_AU",MEKE_viscCoeff_Au)
  Use_KH_in_MEKE = .false. ; call read_param(param_file,"USE_KH_IN_MEKE", Use_KH_in_MEKE)
  sqg_use_MEKE = .false. ; call read_param(param_file,"SQG_USE_MEKE", sqg_use_MEKE)

  if (.not. useMEKE) return

! Allocate memory
  call MOM_mesg("MEKE_alloc_register_restart: allocating and registering", 5)
  isd = HI%isd ; ied = HI%ied ; jsd = HI%jsd ; jed = HI%jed
  allocate(MEKE%MEKE(isd:ied,jsd:jed), source=0.0)
  call register_restart_field(MEKE%MEKE, "MEKE", .false., restart_CS, &
           longname="Mesoscale Eddy Kinetic Energy", units="m2 s-2", conversion=US%L_T_to_m_s**2)

  if (MEKE_GMcoeff>=0.) allocate(MEKE%GM_src(isd:ied,jsd:jed), source=0.0)
  if (MEKE_FrCoeff>=0. .or. MEKE_bhFrCoeff>=0. .or. MEKE_GMECoeff>=0.) &
    allocate(MEKE%mom_src(isd:ied,jsd:jed), source=0.0)
  if (MEKE_bhFrCoeff >= 0.) &
    allocate(MEKE%mom_src_bh(isd:ied,jsd:jed), source=0.0)
  if (MEKE_FrCoeff<0.) MEKE_FrCoeff = 0.
  if (MEKE_bhFrCoeff<0.) MEKE_bhFrCoeff = 0.
  if (MEKE_GMECoeff>=0.) allocate(MEKE%GME_snk(isd:ied,jsd:jed), source=0.0)
  if (MEKE_KhCoeff>=0.) then
    allocate(MEKE%Kh(isd:ied,jsd:jed), source=0.0)
    call register_restart_field(MEKE%Kh, "MEKE_Kh", .false., restart_CS, &
             longname="Lateral diffusivity from Mesoscale Eddy Kinetic Energy", &
             units="m2 s-1", conversion=US%L_to_m**2*US%s_to_T)
  endif
  allocate(MEKE%Rd_dx_h(isd:ied,jsd:jed), source=0.0)
  if (MEKE_viscCoeff_Ku/=0.) then
    allocate(MEKE%Ku(isd:ied,jsd:jed), source=0.0)
    call register_restart_field(MEKE%Ku, "MEKE_Ku", .false., restart_CS, &
             longname="Lateral viscosity from Mesoscale Eddy Kinetic Energy", &
             units="m2 s-1", conversion=US%L_to_m**2*US%s_to_T)
  endif
  if (sqg_use_MEKE) then
    allocate(MEKE%Le(isd:ied,jsd:jed), source=0.0)
    call register_restart_field(MEKE%Le, "MEKE_Le", .false., restart_CS, &
             longname="Eddy length scale from Mesoscale Eddy Kinetic Energy", &
             units="m", conversion=US%L_to_m)
  endif
  if (Use_Kh_in_MEKE) then
    allocate(MEKE%Kh_diff(isd:ied,jsd:jed), source=0.0)
    call register_restart_field(MEKE%Kh_diff, "MEKE_Kh_diff", .false., restart_CS, &
             longname="Copy of thickness diffusivity for diffusing MEKE", &
             units="m2 s-1", conversion=US%L_to_m**2*US%s_to_T)
  endif

  if (MEKE_viscCoeff_Au/=0.) then
    allocate(MEKE%Au(isd:ied,jsd:jed), source=0.0)
    call register_restart_field(MEKE%Au, "MEKE_Au", .false., restart_CS, &
             longname="Lateral biharmonic viscosity from Mesoscale Eddy Kinetic Energy", &
             units="m4 s-1", conversion=US%L_to_m**4*US%s_to_T)
  endif

end procedure MEKE_alloc_register_restart
module procedure MEKE_end
  if (allocated(MEKE%Au)) deallocate(MEKE%Au)
  if (allocated(MEKE%Kh_diff)) deallocate(MEKE%Kh_diff)
  if (allocated(MEKE%Ku)) deallocate(MEKE%Ku)
  if (allocated(MEKE%Rd_dx_h)) deallocate(MEKE%Rd_dx_h)
  if (allocated(MEKE%Kh)) deallocate(MEKE%Kh)
  if (allocated(MEKE%GME_snk)) deallocate(MEKE%GME_snk)
  if (allocated(MEKE%mom_src)) deallocate(MEKE%mom_src)
  if (allocated(MEKE%mom_src_bh)) deallocate(MEKE%mom_src_bh)
  if (allocated(MEKE%GM_src)) deallocate(MEKE%GM_src)
  if (allocated(MEKE%MEKE)) deallocate(MEKE%MEKE)
  if (allocated(MEKE%Le)) deallocate(MEKE%Le)
end procedure MEKE_end
end submodule MOM_MEKE_s
