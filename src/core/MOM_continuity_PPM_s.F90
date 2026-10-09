! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> The zonal mass flux routines of MOM_continuity_PPM, kept in a submodule of their own so that
!! changing them does not force the modules that use MOM_continuity_PPM to be recompiled.
submodule (MOM_continuity_PPM) MOM_continuity_PPM_s

#include <MOM_memory.h>

implicit none

contains

module procedure zonal_mass_flux
  ! Local variables
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)) :: &
    duhdu, &      ! Partial derivative of uh with u [H L ~> m2 or kg m-1].
    visc_rem      ! A copy of visc_rem_u or an array of 1's [nondim].
  real, dimension(SZIB_(G),SZJ_(G)) :: &
    du, &         ! Corrective barotropic change in the velocity to give uhbt [L T-1 ~> m s-1].
    du_min_CFL, & ! Lower limit on du correction to avoid CFL violations [L T-1 ~> m s-1]
    du_max_CFL, & ! Upper limit on du correction to avoid CFL violations [L T-1 ~> m s-1]
    duhdu_tot_0, & ! Summed partial derivative of uh with u [H L ~> m2 or kg m-1].
    uh_tot_0, &   ! Summed transport with no barotropic correction [H L2 T-1 ~> m3 s-1 or kg s-1].
    visc_rem_max, & ! The column maximum of visc_rem [nondim].
    FAuI          ! The sums of zonal face areas [H L ~> m2 or kg m-1].
  logical, dimension(SZIB_(G),SZJ_(G)) :: &
    do_I, &       ! Indicates the points where the barotropic and baroclinic transports are reconciled
    simple_OBC_pt ! Indicates points with specified transport OBCs
  integer, dimension(SZIB_(G),SZJ_(G)) :: &
    open_dir      ! 1 or -1 at faces on open boundary segments where the flow is taken from the
                  ! cell to the west or east of the face, or 0 elsewhere [nondim]
  real :: FA_u    ! A sum of zonal face areas [H L ~> m2 or kg m-1].
  real :: I_vrm   ! 1.0 / visc_rem_max [nondim]
  real :: CFL_dt  ! The maximum CFL ratio of the adjusted velocities divided by
                  ! the time step [T-1 ~> s-1].
  real :: I_dt    ! 1.0 / dt [T-1 ~> s-1].
  real :: du_lim  ! The velocity change that give a relative CFL of 1 [L T-1 ~> m s-1].
  real :: dx_E, dx_W ! Effective x-grid spacings to the east and west [L ~> m].
  real :: h_marg_min ! A copy of CS%h_marg_min for use on the device [H ~> m or kg m-2]
  real :: H_subroundoff ! A copy of GV%H_subroundoff for use on the device [H ~> m or kg m-2]
  type(cont_loop_bounds_type) :: LB
  integer :: i, j, k, ish, ieh, jsh, jeh, n, nz
  integer :: IsdB, IedB, jsd, jed ! The data domain bounds at u points
  integer :: l_seg ! The OBC segment number
  logical :: use_visc_rem, set_BT_cont, set_h_u
  logical :: vol_CFL, aggress_adjust, use_visc_rem_max ! Copies of CS fields for use on the device
  logical :: local_specified_BC, local_Flather_OBC, local_open_BC, any_simple_OBC  ! OBC-related logicals

  use_visc_rem = present(visc_rem_u)

  set_BT_cont = .false. ; if (present(BT_cont)) set_BT_cont = (associated(BT_cont))
  set_h_u = .false. ; if (set_BT_cont) set_h_u = allocated(BT_cont%h_u)

  !   The arrays that are used on the device are copied in before the clock is started.  The output
  ! arrays that are only partly set here are copied in too, so that the copies back to the host
  ! leave the rest of them unchanged.
  !$omp target enter data map(to: G)
  !$omp target enter data map(to: G%dy_Cu, G%dxCu, G%IareaT, G%areaT, G%IdxT, G%dxT, G%mask2dCu)
  !$omp target enter data map(to: u, h_in, h_W, h_E, por_face_areaU, uh)
  if (use_visc_rem) then
    !$omp target enter data map(to: visc_rem_u)
  endif
  if (present(uhbt)) then
    !$omp target enter data map(to: uhbt)
  endif
  if (present(u_cor)) then
    !$omp target enter data map(to: u_cor)
  endif
  if (present(du_cor)) then
    !$omp target enter data map(alloc: du_cor)
  endif
  if (set_BT_cont) call zonal_BT_cont_to_device(BT_cont, set_h_u)

  call cpu_clock_begin(id_clock_correct)

  !$omp target enter data map(alloc: duhdu, visc_rem, du, du_min_CFL, du_max_CFL, duhdu_tot_0, &
  !$omp                              uh_tot_0, visc_rem_max, FAuI, do_I, simple_OBC_pt, open_dir)

  local_specified_BC = .false. ; local_Flather_OBC = .false. ; local_open_BC = .false.
  if (associated(OBC)) then ; if (OBC%OBC_pe) then
    local_specified_BC = OBC%specified_u_BCs_exist_globally
    local_Flather_OBC = OBC%Flather_u_BCs_exist_globally
    local_open_BC = OBC%open_u_BCs_exist_globally
  endif ; endif

  if (present(LB_in)) then
    LB = LB_in
  else
    LB%ish = G%isc ; LB%ieh = G%iec ; LB%jsh = G%jsc ; LB%jeh = G%jec
  endif
  ish = LB%ish ; ieh = LB%ieh ; jsh = LB%jsh ; jeh = LB%jeh ; nz = GV%ke
  IsdB = G%IsdB ; IedB = G%IedB ; jsd = G%jsd ; jed = G%jed

  vol_CFL = CS%vol_CFL ; aggress_adjust = CS%aggress_adjust ; use_visc_rem_max = CS%use_visc_rem_max
  h_marg_min = CS%h_marg_min ; H_subroundoff = GV%H_subroundoff

  if (present(du_cor)) then
    !$omp target teams distribute parallel do collapse(2)
    do j=jsd,jed ; do I=IsdB,IedB
      du_cor(I,j) = 0.0
    enddo ; enddo
  endif

  CFL_dt = CS%CFL_limit_adjust / dt
  I_dt = 1.0 / dt
  if (aggress_adjust) CFL_dt = I_dt

  if (local_open_BC) then
    ! Note which faces are on open boundary segments, so that OBC need not be used on the device.
    do j=jsh,jeh ; do I=ish-1,ieh
      open_dir(I,j) = 0
      if (OBC%segnum_u(I,j) /= 0) then
        if (OBC%segment(abs(OBC%segnum_u(I,j)))%open) open_dir(I,j) = sign(1, OBC%segnum_u(I,j))
      endif
    enddo ; enddo
    !$omp target update to(open_dir)
  endif

  ! Set uh and duhdu.
  if (use_visc_rem) then
    !$omp target teams distribute parallel do collapse(3)
    do k=1,nz ; do j=jsh,jeh ; do I=ish-1,ieh
      visc_rem(I,j,k) = visc_rem_u(I,j,k)
    enddo ; enddo ; enddo
  else
    !$omp target teams distribute parallel do collapse(3)
    do k=1,nz ; do j=jsh,jeh ; do I=ish-1,ieh
      visc_rem(I,j,k) = 1.0
    enddo ; enddo ; enddo
  endif
  !$omp target teams distribute parallel do collapse(3)
  do k=1,nz ; do j=jsh,jeh ; do I=ish-1,ieh
    call flux_elem(u(I,j,k), h_in(i,j,k), h_in(i+1,j,k), h_W(i,j,k), h_W(i+1,j,k), h_E(i,j,k), &
                   h_E(i+1,j,k), uh(I,j,k), duhdu(I,j,k), visc_rem(I,j,k), G%dy_Cu(I,j), &
                   G%IareaT(i,j), G%IareaT(i+1,j), G%IdxT(i,j), G%IdxT(i+1,j), dt, vol_CFL, &
                   por_face_areaU(I,j,k), h_marg_min)
    if (local_open_BC) &
      call flux_elem_OBC(u(I,j,k), h_in(i,j,k), h_in(i+1,j,k), uh(I,j,k), duhdu(I,j,k), &
                         visc_rem(I,j,k), por_face_areaU(I,j,k), G%dy_Cu(I,j), h_marg_min, &
                         open_dir(I,j))
  enddo ; enddo ; enddo
  if (local_specified_BC) then
    !$omp target teams distribute parallel do collapse(3) private(l_seg)
    do k=1,nz ; do j=jsh,jeh ; do I=ish-1,ieh ; if (OBC%segnum_u(I,j) /= 0) then
      l_seg = abs(OBC%segnum_u(I,j))
      if (OBC%segment(l_seg)%specified) uh(I,j,k) = OBC%segment(l_seg)%normal_trans(I,j,k)
    endif ; enddo ; enddo ; enddo
  endif

  if (present(uhbt) .or. set_BT_cont) then
    !$omp target teams
    if (use_visc_rem .and. use_visc_rem_max) then
      !$omp distribute parallel do collapse(2)
      do j=jsh,jeh ; do I=ish-1,ieh
        visc_rem_max(I,j) = 0.0
      enddo ; enddo
      do k=1,nz
        !$omp distribute parallel do collapse(2)
        do j=jsh,jeh ; do I=ish-1,ieh
          visc_rem_max(I,j) = max(visc_rem_max(I,j), visc_rem(I,j,k))
        enddo ; enddo
      enddo
    else
      !$omp distribute parallel do collapse(2)
      do j=jsh,jeh ; do I=ish-1,ieh
        visc_rem_max(I,j) = 1.0
      enddo ; enddo
    endif
    !   Set limits on du that will keep the CFL number between -1 and 1.
    ! This should be adequate to keep the root bracketed in all cases.
    !$omp distribute parallel do collapse(2) private(I_vrm, dx_W, dx_E)
    do j=jsh,jeh ; do I=ish-1,ieh
      I_vrm = 0.0
      if (visc_rem_max(I,j) > 0.0) I_vrm = 1.0 / visc_rem_max(I,j)
      if (vol_CFL) then
        dx_W = ratio_max(G%areaT(i,j), G%dy_Cu(I,j), 1000.0*G%dxT(i,j))
        dx_E = ratio_max(G%areaT(i+1,j), G%dy_Cu(I,j), 1000.0*G%dxT(i+1,j))
      else ; dx_W = G%dxT(i,j) ; dx_E = G%dxT(i+1,j) ; endif
      du_max_CFL(I,j) = 2.0* (CFL_dt * dx_W) * I_vrm
      du_min_CFL(I,j) = -2.0 * (CFL_dt * dx_E) * I_vrm
      uh_tot_0(I,j) = 0.0 ; duhdu_tot_0(I,j) = 0.0
    enddo ; enddo
    do k=1,nz
      !$omp distribute parallel do collapse(2)
      do j=jsh,jeh ; do I=ish-1,ieh
        duhdu_tot_0(I,j) = duhdu_tot_0(I,j) + duhdu(I,j,k)
        uh_tot_0(I,j) = uh_tot_0(I,j) + uh(I,j,k)
      enddo ; enddo
      if (use_visc_rem) then
        if (aggress_adjust) then
          !$omp distribute parallel do collapse(2) private(dx_W, dx_E, du_lim)
          do j=jsh,jeh ; do I=ish-1,ieh
            if (vol_CFL) then
              dx_W = ratio_max(G%areaT(i,j), G%dy_Cu(I,j), 1000.0*G%dxT(i,j))
              dx_E = ratio_max(G%areaT(i+1,j), G%dy_Cu(I,j), 1000.0*G%dxT(i+1,j))
            else ; dx_W = G%dxT(i,j) ; dx_E = G%dxT(i+1,j) ; endif

            du_lim = 0.499*((dx_W*I_dt - u(I,j,k)) + MIN(0.0,u(I-1,j,k)))
            if (du_max_CFL(I,j) * visc_rem(I,j,k) > du_lim) &
              du_max_CFL(I,j) = du_lim / visc_rem(I,j,k)

            du_lim = 0.499*((-dx_E*I_dt - u(I,j,k)) + MAX(0.0,u(I+1,j,k)))
            if (du_min_CFL(I,j) * visc_rem(I,j,k) < du_lim) &
              du_min_CFL(I,j) = du_lim / visc_rem(I,j,k)
          enddo ; enddo
        else
          !$omp distribute parallel do collapse(2) private(dx_W, dx_E)
          do j=jsh,jeh ; do I=ish-1,ieh
            if (vol_CFL) then
              dx_W = ratio_max(G%areaT(i,j), G%dy_Cu(I,j), 1000.0*G%dxT(i,j))
              dx_E = ratio_max(G%areaT(i+1,j), G%dy_Cu(I,j), 1000.0*G%dxT(i+1,j))
            else ; dx_W = G%dxT(i,j) ; dx_E = G%dxT(i+1,j) ; endif

            if (du_max_CFL(I,j) * visc_rem(I,j,k) > dx_W*CFL_dt - u(I,j,k)*G%mask2dCu(I,j)) &
              du_max_CFL(I,j) = (dx_W*CFL_dt - u(I,j,k)) / visc_rem(I,j,k)
            if (du_min_CFL(I,j) * visc_rem(I,j,k) < -dx_E*CFL_dt - u(I,j,k)*G%mask2dCu(I,j)) &
              du_min_CFL(I,j) = -(dx_E*CFL_dt + u(I,j,k)) / visc_rem(I,j,k)
          enddo ; enddo
        endif
      else
        if (aggress_adjust) then
          !$omp distribute parallel do collapse(2) private(dx_W, dx_E)
          do j=jsh,jeh ; do I=ish-1,ieh
            if (vol_CFL) then
              dx_W = ratio_max(G%areaT(i,j), G%dy_Cu(I,j), 1000.0*G%dxT(i,j))
              dx_E = ratio_max(G%areaT(i+1,j), G%dy_Cu(I,j), 1000.0*G%dxT(i+1,j))
            else ; dx_W = G%dxT(i,j) ; dx_E = G%dxT(i+1,j) ; endif

            du_max_CFL(I,j) = MIN(du_max_CFL(I,j), 0.499 * &
                        ((dx_W*I_dt - u(I,j,k)) + MIN(0.0,u(I-1,j,k))) )
            du_min_CFL(I,j) = MAX(du_min_CFL(I,j), 0.499 * &
                        ((-dx_E*I_dt - u(I,j,k)) + MAX(0.0,u(I+1,j,k))) )
          enddo ; enddo
        else
          !$omp distribute parallel do collapse(2) private(dx_W, dx_E)
          do j=jsh,jeh ; do I=ish-1,ieh
            if (vol_CFL) then
              dx_W = ratio_max(G%areaT(i,j), G%dy_Cu(I,j), 1000.0*G%dxT(i,j))
              dx_E = ratio_max(G%areaT(i+1,j), G%dy_Cu(I,j), 1000.0*G%dxT(i+1,j))
            else ; dx_W = G%dxT(i,j) ; dx_E = G%dxT(i+1,j) ; endif

            du_max_CFL(I,j) = MIN(du_max_CFL(I,j), dx_W*CFL_dt - u(I,j,k))
            du_min_CFL(I,j) = MAX(du_min_CFL(I,j), -(dx_E*CFL_dt + u(I,j,k)))
          enddo ; enddo
        endif
      endif
    enddo
    !$omp distribute parallel do collapse(2)
    do j=jsh,jeh ; do I=ish-1,ieh
      du_max_CFL(I,j) = max(du_max_CFL(I,j),0.0)
      du_min_CFL(I,j) = min(du_min_CFL(I,j),0.0)
    enddo ; enddo
    !$omp end target teams

    any_simple_OBC = .false.
    if (local_specified_BC .or. local_Flather_OBC) then
      !$omp target teams distribute parallel do collapse(2) private(l_seg) &
      !$omp   reduction(.or.: any_simple_OBC) map(tofrom: any_simple_OBC)
      do j=jsh,jeh ; do I=ish-1,ieh
        l_seg = abs(OBC%segnum_u(I,j))

        ! Avoid reconciling barotropic/baroclinic transports if transport is specified
        simple_OBC_pt(I,j) = .false.
        if (l_seg /= OBC_NONE) simple_OBC_pt(I,j) = OBC%segment(l_seg)%specified
        do_I(I,j) = .not.simple_OBC_pt(I,j)
        any_simple_OBC = any_simple_OBC .or. simple_OBC_pt(I,j)
      enddo ; enddo
    else
      !$omp target teams distribute parallel do collapse(2)
      do j=jsh,jeh ; do I=ish-1,ieh
        do_I(I,j) = .true.
      enddo ; enddo
    endif

    if (present(uhbt)) then
      ! Find du and uh.
      call zonal_flux_adjust(u, h_in, h_W, h_E, uhbt, uh_tot_0, duhdu_tot_0, du, &
                             du_max_CFL, du_min_CFL, dt, G, GV, US, CS, visc_rem, &
                             ish-1, ieh, jsh, jeh, do_I, por_face_areaU, uh, &
                             local_open_BC, open_dir)

      if (present(u_cor)) then
        !$omp target teams distribute parallel do collapse(3)
        do k=1,nz ; do j=jsh,jeh ; do I=ish-1,ieh
          u_cor(I,j,k) = u(I,j,k) + du(I,j) * visc_rem(I,j,k)
        enddo ; enddo ; enddo
        if (any_simple_OBC) then
          !$omp target teams distribute parallel do collapse(3)
          do k=1,nz ; do j=jsh,jeh ; do I=ish-1,ieh ; if (simple_OBC_pt(I,j)) then
            u_cor(I,j,k) = OBC%segment(abs(OBC%segnum_u(I,j)))%normal_vel(I,j,k)
          endif ; enddo ; enddo ; enddo
        endif
      endif ! u-corrected

      if (present(du_cor)) then
        !$omp target teams distribute parallel do collapse(2)
        do j=jsh,jeh ; do I=ish-1,ieh
          du_cor(I,j) = du(I,j)
        enddo ; enddo
      endif

    endif

    if (set_BT_cont) then
      call set_zonal_BT_cont(u, h_in, h_W, h_E, BT_cont, uh_tot_0, duhdu_tot_0, &
                             du_max_CFL, du_min_CFL, dt, G, GV, US, CS, visc_rem, &
                             visc_rem_max, ish-1, ieh, jsh, jeh, do_I, por_face_areaU, open_dir)
      if (any_simple_OBC) then
        !$omp target teams
        !$omp distribute parallel do collapse(2)
        do j=jsh,jeh ; do I=ish-1,ieh
          if (simple_OBC_pt(I,j)) FAuI(I,j) = H_subroundoff*G%dy_Cu(I,j)
        enddo ; enddo
        ! NOTE: simple_OBC_pt should prevent access to segment OBC_NONE
        do k=1,nz
          !$omp distribute parallel do collapse(2) private(l_seg)
          do j=jsh,jeh ; do I=ish-1,ieh ; if (simple_OBC_pt(I,j)) then
            l_seg = abs(OBC%segnum_u(I,j))
            if ((abs(OBC%segment(l_seg)%normal_vel(I,j,k)) > 0.0) .and. (OBC%segment(l_seg)%specified)) &
              FAuI(I,j) = FAuI(I,j) + OBC%segment(l_seg)%normal_trans(I,j,k) / OBC%segment(l_seg)%normal_vel(I,j,k)
          endif ; enddo ; enddo
        enddo
        !$omp distribute parallel do collapse(2)
        do j=jsh,jeh ; do I=ish-1,ieh ; if (simple_OBC_pt(I,j)) then
          BT_cont%FA_u_W0(I,j) = FAuI(I,j) ; BT_cont%FA_u_E0(I,j) = FAuI(I,j)
          BT_cont%FA_u_WW(I,j) = FAuI(I,j) ; BT_cont%FA_u_EE(I,j) = FAuI(I,j)
          BT_cont%uBT_WW(I,j) = 0.0 ; BT_cont%uBT_EE(I,j) = 0.0
        endif ; enddo ; enddo
        !$omp end target teams
      endif
    endif ! set_BT_cont

  endif ! present(uhbt) or set_BT_cont

  if (local_open_BC .and. set_BT_cont) then
    ! This is done on the host, between copies of the face areas from and back to the device.
    !$omp target update from(BT_cont%FA_u_W0, BT_cont%FA_u_E0, BT_cont%FA_u_WW, &
    !$omp                    BT_cont%FA_u_EE, BT_cont%uBT_WW, BT_cont%uBT_EE)
    do n = 1, OBC%number_of_segments
      if (OBC%segment(n)%open .and. OBC%segment(n)%is_E_or_W) then
        I = OBC%segment(n)%HI%IsdB
        if (OBC%segment(n)%direction == OBC_DIRECTION_E) then
          do j = OBC%segment(n)%HI%Jsd, OBC%segment(n)%HI%Jed
            FA_u = 0.0
            do k=1,nz ; FA_u = FA_u + h_in(i,j,k)*(G%dy_Cu(I,j)*por_face_areaU(I,j,k)) ; enddo
            BT_cont%FA_u_W0(I,j) = FA_u ; BT_cont%FA_u_E0(I,j) = FA_u
            BT_cont%FA_u_WW(I,j) = FA_u ; BT_cont%FA_u_EE(I,j) = FA_u
            BT_cont%uBT_WW(I,j) = 0.0 ; BT_cont%uBT_EE(I,j) = 0.0
          enddo
        else
          do j = OBC%segment(n)%HI%Jsd, OBC%segment(n)%HI%Jed
            FA_u = 0.0
            do k=1,nz ; FA_u = FA_u + h_in(i+1,j,k)*(G%dy_Cu(I,j)*por_face_areaU(I,j,k)) ; enddo
            BT_cont%FA_u_W0(I,j) = FA_u ; BT_cont%FA_u_E0(I,j) = FA_u
            BT_cont%FA_u_WW(I,j) = FA_u ; BT_cont%FA_u_EE(I,j) = FA_u
            BT_cont%uBT_WW(I,j) = 0.0 ; BT_cont%uBT_EE(I,j) = 0.0
          enddo
        endif
      endif
    enddo
    !$omp target update to(BT_cont%FA_u_W0, BT_cont%FA_u_E0, BT_cont%FA_u_WW, &
    !$omp                  BT_cont%FA_u_EE, BT_cont%uBT_WW, BT_cont%uBT_EE)
  endif

  if (set_h_u) then
    if (present(u_cor)) then
      call zonal_flux_thickness(u_cor, h_in, h_W, h_E, BT_cont%h_u, dt, G, GV, US, LB, &
                                vol_CFL, CS%marginal_faces, OBC, por_face_areaU, visc_rem_u)
    else
      call zonal_flux_thickness(u, h_in, h_W, h_E, BT_cont%h_u, dt, G, GV, US, LB, &
                                vol_CFL, CS%marginal_faces, OBC, por_face_areaU, visc_rem_u)
    endif
  endif

  !$omp target exit data map(release: duhdu, visc_rem, du, du_min_CFL, du_max_CFL, duhdu_tot_0, &
  !$omp                               uh_tot_0, visc_rem_max, FAuI, do_I, simple_OBC_pt, open_dir)

  call cpu_clock_end(id_clock_correct)

  ! The results are copied back to the host, and the inputs released, after the clock is stopped.
  if (set_BT_cont) call zonal_BT_cont_from_device(BT_cont, set_h_u)
  if (present(du_cor)) then
    !$omp target exit data map(from: du_cor)
  endif
  if (present(u_cor)) then
    !$omp target exit data map(from: u_cor)
  endif
  if (present(uhbt)) then
    !$omp target exit data map(release: uhbt)
  endif
  if (use_visc_rem) then
    !$omp target exit data map(release: visc_rem_u)
  endif
  !$omp target exit data map(from: uh) map(release: u, h_in, h_W, h_E, por_face_areaU)
  !$omp target exit data map(release: G%dy_Cu, G%dxCu, G%IareaT, G%areaT, G%IdxT, G%dxT, G%mask2dCu)
  !$omp target exit data map(release: G)

end procedure zonal_mass_flux

module procedure zonal_flux_thickness
  ! Local variables
  real :: CFL  ! The CFL number based on the local velocity and grid spacing [nondim]
  real :: curv_3 ! A measure of the thickness curvature over a grid length [H ~> m or kg m-2]
  real :: h_avg  ! The average thickness of a flux [H ~> m or kg m-2].
  real :: h_marg ! The marginal thickness of a flux [H ~> m or kg m-2].
  logical :: local_open_BC
  integer :: i, j, k, ish, ieh, jsh, jeh, nz, n
  ish = LB%ish ; ieh = LB%ieh ; jsh = LB%jsh ; jeh = LB%jeh ; nz = GV%ke

  !   This works on the device, with all of the arrays already there.
  !$omp target teams distribute parallel do collapse(3) private(CFL, curv_3, h_marg, h_avg)
  do k=1,nz ; do j=jsh,jeh ; do I=ish-1,ieh
    if (u(I,j,k) > 0.0) then
      if (vol_CFL) then ; CFL = (u(I,j,k) * dt) * (G%dy_Cu(I,j) * G%IareaT(i,j))
      else ; CFL = u(I,j,k) * dt * G%IdxT(i,j) ; endif
      curv_3 = (h_W(i,j,k) + h_E(i,j,k)) - 2.0*h(i,j,k)
      h_avg = h_E(i,j,k) + CFL * (0.5*(h_W(i,j,k) - h_E(i,j,k)) + curv_3*(CFL - 1.5))
      h_marg = h_E(i,j,k) + CFL * ((h_W(i,j,k) - h_E(i,j,k)) + 3.0*curv_3*(CFL - 1.0))
    elseif (u(I,j,k) < 0.0) then
      if (vol_CFL) then ; CFL = (-u(I,j,k)*dt) * (G%dy_Cu(I,j) * G%IareaT(i+1,j))
      else ; CFL = -u(I,j,k) * dt * G%IdxT(i+1,j) ; endif
      curv_3 = (h_W(i+1,j,k) + h_E(i+1,j,k)) - 2.0*h(i+1,j,k)
      h_avg = h_W(i+1,j,k) + CFL * (0.5*(h_E(i+1,j,k)-h_W(i+1,j,k)) + curv_3*(CFL - 1.5))
      h_marg = h_W(i+1,j,k) + CFL * ((h_E(i+1,j,k)-h_W(i+1,j,k)) + &
                                    3.0*curv_3*(CFL - 1.0))
    else
      h_avg = 0.5 * (h_W(i+1,j,k) + h_E(i,j,k))
      !   The choice to use the arithmetic mean here is somewhat arbitrarily, but
      ! it should be noted that h_W(i+1,j,k) and h_E(i,j,k) are usually the same.
      h_marg = 0.5 * (h_W(i+1,j,k) + h_E(i,j,k))
 !    h_marg = (2.0 * h_W(i+1,j,k) * h_E(i,j,k)) / &
 !             (h_W(i+1,j,k) + h_E(i,j,k) + GV%H_subroundoff)
    endif

    if (marginal) then ; h_u(I,j,k) = h_marg
    else ; h_u(I,j,k) = h_avg ; endif
  enddo ; enddo ; enddo
  if (present(visc_rem_u)) then
    ! Scale back the thickness to account for the effects of viscosity and the fractional open
    ! thickness to give an appropriate non-normalized weight for each layer in determining the
    ! barotropic acceleration.
    !$omp target teams distribute parallel do collapse(3)
    do k=1,nz ; do j=jsh,jeh ; do I=ish-1,ieh
      h_u(I,j,k) = h_u(I,j,k) * (visc_rem_u(I,j,k) * por_face_areaU(I,j,k))
    enddo ; enddo ; enddo
  else
    !$omp target teams distribute parallel do collapse(3)
    do k=1,nz ; do j=jsh,jeh ; do I=ish-1,ieh
      h_u(I,j,k) = h_u(I,j,k) * por_face_areaU(I,j,k)
    enddo ; enddo ; enddo
  endif

  local_open_BC = .false.
  if (associated(OBC)) local_open_BC = OBC%open_u_BCs_exist_globally
  if (local_open_BC) then
    ! This is done on the host, between copies of h_u from and back to the device.
    !$omp target update from(h_u)
    do n = 1, OBC%number_of_segments
      if (OBC%segment(n)%open .and. OBC%segment(n)%is_E_or_W) then
        I = OBC%segment(n)%HI%IsdB
        if (OBC%segment(n)%direction == OBC_DIRECTION_E) then
          if (present(visc_rem_u)) then ; do k=1,nz
            do j = OBC%segment(n)%HI%jsd, OBC%segment(n)%HI%jed
              h_u(I,j,k) = h(i,j,k) * (visc_rem_u(I,j,k) * por_face_areaU(I,j,k))
            enddo
          enddo ; else ; do k=1,nz
            do j = OBC%segment(n)%HI%jsd, OBC%segment(n)%HI%jed
              h_u(I,j,k) = h(i,j,k) * por_face_areaU(I,j,k)
            enddo
          enddo ; endif
        else
          if (present(visc_rem_u)) then ; do k=1,nz
            do j = OBC%segment(n)%HI%jsd, OBC%segment(n)%HI%jed
              h_u(I,j,k) = h(i+1,j,k) * (visc_rem_u(I,j,k) * por_face_areaU(I,j,k))
            enddo
          enddo ; else ; do k=1,nz
            do j = OBC%segment(n)%HI%jsd, OBC%segment(n)%HI%jed
              h_u(I,j,k) = h(i+1,j,k) * por_face_areaU(I,j,k)
            enddo
          enddo ; endif
        endif
      endif
    enddo
    !$omp target update to(h_u)
  endif

end procedure zonal_flux_thickness


!> Copies BT_cont and the components that zonal_mass_flux sets to the device.  This is done here,
!! with BT_cont passed as a plain argument, because when amdflang maps the components of a structure
!! that is reached through a pointer argument, it takes the structure to be as many times its size
!! as there are elements in a component, and the mapping fails.
subroutine zonal_BT_cont_to_device(BT_cont, set_h_u)
  type(BT_cont_type), intent(inout) :: BT_cont !< A structure with elements that describe the
                                               !! effective open face areas as a function of barotropic flow.
  logical,            intent(in)    :: set_h_u !< If true, BT_cont%h_u is copied as well

  ! BT_cont is mapped apart from, and before, its components.
  !$omp target enter data map(to: BT_cont)
  !$omp target enter data map(to: BT_cont%FA_u_EE, BT_cont%FA_u_E0, BT_cont%FA_u_W0, &
  !$omp                           BT_cont%FA_u_WW, BT_cont%uBT_WW, BT_cont%uBT_EE)
  if (set_h_u) then
    !$omp target enter data map(to: BT_cont%h_u)
  endif
end subroutine zonal_BT_cont_to_device

!> Copies the components of BT_cont that zonal_mass_flux sets back from the device, and releases
!! BT_cont there.
subroutine zonal_BT_cont_from_device(BT_cont, set_h_u)
  type(BT_cont_type), intent(inout) :: BT_cont !< A structure with elements that describe the
                                               !! effective open face areas as a function of barotropic flow.
  logical,            intent(in)    :: set_h_u !< If true, BT_cont%h_u is copied as well

  if (set_h_u) then
    !$omp target exit data map(from: BT_cont%h_u)
  endif
  !$omp target exit data map(from: BT_cont%FA_u_EE, BT_cont%FA_u_E0, BT_cont%FA_u_W0, &
  !$omp                          BT_cont%FA_u_WW, BT_cont%uBT_WW, BT_cont%uBT_EE)
  ! Released after, and apart from, its components.
  !$omp target exit data map(release: BT_cont)
end subroutine zonal_BT_cont_from_device

!> Returns the barotropic velocity adjustment that gives the
!! desired barotropic (layer-summed) transport.
subroutine zonal_flux_adjust(u, h_in, h_W, h_E, uhbt, uh_tot_0, duhdu_tot_0, &
                             du, du_max_CFL, du_min_CFL, dt, G, GV, US, CS, visc_rem, &
                             i_start, i_end, j_start, j_end, do_I_in, por_face_areaU, uh_3d, &
                             local_open_BC, open_dir)
  type(ocean_grid_type),                     intent(in)    :: G    !< Ocean's grid structure.
  type(verticalGrid_type),                   intent(in)    :: GV   !< Ocean's vertical grid structure.
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), intent(in)   :: u    !< Zonal velocity [L T-1 ~> m s-1].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in)    :: h_in !< Layer thickness used to
                                                                   !! calculate fluxes [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in)    :: h_W  !< West edge thickness in the
                                                                   !! reconstruction [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in)    :: h_E  !< East edge thickness in the
                                                                   !! reconstruction [H ~> m or kg m-2].
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), intent(in)   :: visc_rem !< Both the fraction of the
                       !! momentum originally in a layer that remains after a time-step of viscosity, and
                       !! the fraction of a time-step's worth of a barotropic acceleration that a layer
                       !! experiences after viscosity is applied [nondim].
                       !! Visc_rem is between 0 (at the bottom) and 1 (far above the bottom).
  real, dimension(SZIB_(G),SZJ_(G)),         intent(in)    :: uhbt !< The summed volume flux
                       !! through zonal faces [H L2 T-1 ~> m3 s-1 or kg s-1].
  real, dimension(SZIB_(G),SZJ_(G)),         intent(in)    :: du_max_CFL  !< Maximum acceptable
                       !! value of du [L T-1 ~> m s-1].
  real, dimension(SZIB_(G),SZJ_(G)),         intent(in)    :: du_min_CFL  !< Minimum acceptable
                       !! value of du [L T-1 ~> m s-1].
  real, dimension(SZIB_(G),SZJ_(G)),         intent(in)    :: uh_tot_0    !< The summed transport
                       !! with 0 adjustment [H L2 T-1 ~> m3 s-1 or kg s-1].
  real, dimension(SZIB_(G),SZJ_(G)),         intent(in)    :: duhdu_tot_0 !< The partial derivative
                       !! of du_err with du at 0 adjustment [H L ~> m2 or kg m-1].
  real, dimension(SZIB_(G),SZJ_(G)),         intent(out)   :: du !<
                       !! The barotropic velocity adjustment [L T-1 ~> m s-1].
  real,                                      intent(in)    :: dt   !< Time increment [T ~> s].
  type(unit_scale_type),                     intent(in)    :: US   !< A dimensional unit scaling type
  type(continuity_PPM_CS),                   intent(in)    :: CS   !< This module's control structure.
  integer,                                   intent(in)    :: i_start !< Start of the I index range
  integer,                                   intent(in)    :: i_end   !< End of the I index range
  integer,                                   intent(in)    :: j_start !< Start of the j index range
  integer,                                   intent(in)    :: j_end   !< End of the j index range
  logical, dimension(SZIB_(G),SZJ_(G)),      intent(in)    :: do_I_in !< A logical flag indicating
                       !! which points to work on.
  real, dimension(SZIB_(G),SZJ_(G),SZK_(G)), intent(in)    :: por_face_areaU !< fractional open area
                       !! of U-faces [nondim]
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), intent(inout) :: uh_3d !< Volume flux through zonal
                       !! faces = u*h*dy [H L2 T-1 ~> m3 s-1 or kg s-1], updated at the points that
                       !! are adjusted.
  logical,                                   intent(in)    :: local_open_BC !< True if there are open
                       !! boundary faces on this PE whose fluxes are set from open_dir.
  integer, dimension(SZIB_(G),SZJ_(G)),      intent(in)    :: open_dir !< 1 or -1 at faces on open
                       !! boundary segments where the flow is taken from the cell to the west or east
                       !! of the face, or 0 elsewhere [nondim]
  ! Local variables
  real, dimension(SZIB_(G),SZJ_(G)) :: &
    uh_err, &  ! Difference between uhbt and the summed uh [H L2 T-1 ~> m3 s-1 or kg s-1].
    uh_err_best, & ! The smallest value of uh_err found so far [H L2 T-1 ~> m3 s-1 or kg s-1].
    duhdu_tot,&! Summed partial derivative of uh with u [H L ~> m2 or kg m-1].
    du_min, &  ! Lower limit on du correction based on CFL limits and previous iterations [L T-1 ~> m s-1]
    du_max     ! Upper limit on du correction based on CFL limits and previous iterations [L T-1 ~> m s-1]
  logical, dimension(SZIB_(G),SZJ_(G)) :: &
    do_I       ! Indicates the points that are still being adjusted
  real :: u_new   ! The velocity with the correction added [L T-1 ~> m s-1].
  real :: duhdu   ! Partial derivative of uh with u [H L ~> m2 or kg m-1].
  real :: du_prev ! The previous value of du [L T-1 ~> m s-1].
  real :: ddu     ! The change in du from the previous iteration [L T-1 ~> m s-1].
  real :: tol_eta ! The tolerance for the current iteration [H ~> m or kg m-2].
  real :: tol_vel ! The tolerance for velocity in the current iteration [L T-1 ~> m s-1].
  real :: tol_eta_ref ! A copy of CS%tol_eta for use on the device [H ~> m or kg m-2].
  real :: h_marg_min  ! A copy of CS%h_marg_min for use on the device [H ~> m or kg m-2]
  logical :: vol_CFL, better_iter ! Copies of CS fields for use on the device
  integer :: i, j, k, nz, itt
#ifndef _OPENMP
  logical :: domore ! True if any point still needs to be adjusted
#endif
  integer, parameter :: max_itts = 20

  nz = GV%ke

  tol_vel = CS%tol_vel ; tol_eta_ref = CS%tol_eta ; better_iter = CS%better_iter
  vol_CFL = CS%vol_CFL ; h_marg_min = CS%h_marg_min

  !$omp target enter data map(alloc: uh_err, uh_err_best, duhdu_tot, du_min, du_max, do_I)

  !$omp target teams private(tol_eta)
  !$omp distribute parallel do collapse(2)
  do j=j_start,j_end ; do I=i_start,i_end
    du(I,j) = 0.0 ; do_I(I,j) = do_I_in(I,j)
    du_max(I,j) = du_max_CFL(I,j) ; du_min(I,j) = du_min_CFL(I,j)
    uh_err(I,j) = uh_tot_0(I,j) - uhbt(I,j) ; duhdu_tot(I,j) = duhdu_tot_0(I,j)
    uh_err_best(I,j) = abs(uh_err(I,j))
  enddo ; enddo

  do itt=1,max_itts
    select case (itt)
      case (:1) ; tol_eta = 1e-6 * tol_eta_ref
      case (2)  ; tol_eta = 1e-4 * tol_eta_ref
      case (3)  ; tol_eta = 1e-2 * tol_eta_ref
      case default ; tol_eta = tol_eta_ref
    end select

    !$omp distribute parallel do collapse(2)
    do j=j_start,j_end ; do I=i_start,i_end
      if (uh_err(I,j) > 0.0) then ; du_max(I,j) = du(I,j)
      elseif (uh_err(I,j) < 0.0) then ; du_min(I,j) = du(I,j)
      else ; do_I(I,j) = .false. ; endif
    enddo ; enddo
#ifndef _OPENMP
    domore = .false.
#endif
    !$omp distribute parallel do collapse(2) private(ddu, du_prev)
    do j=j_start,j_end ; do I=i_start,i_end ; if (do_I(I,j)) then
      if ((dt * min(G%IareaT(i,j),G%IareaT(i+1,j))*abs(uh_err(I,j)) > tol_eta) .or. &
          (better_iter .and. ((abs(uh_err(I,j)) > tol_vel * duhdu_tot(I,j)) .or. &
                                 (abs(uh_err(I,j)) > uh_err_best(I,j))) )) then
      !   Use Newton's method, provided it stays bounded.  Otherwise bisect
      ! the value with the appropriate bound.
        ddu = -uh_err(I,j) / duhdu_tot(I,j)
        du_prev = du(I,j)
        du(I,j) = du(I,j) + ddu
        if (abs(ddu) < 1.0e-15*abs(du(I,j))) then
          do_I(I,j) = .false. ! ddu is small enough to quit.
        elseif (ddu > 0.0) then
          if (du(I,j) >= du_max(I,j)) then
            du(I,j) = 0.5*(du_prev + du_max(I,j))
            if (du_max(I,j) - du_prev < 1.0e-15*abs(du(I,j))) do_I(I,j) = .false.
          endif
        else ! ddu < 0.0
          if (du(I,j) <= du_min(I,j)) then
            du(I,j) = 0.5*(du_prev + du_min(I,j))
            if (du_prev - du_min(I,j) < 1.0e-15*abs(du(I,j))) do_I(I,j) = .false.
          endif
        endif
#ifndef _OPENMP
        if (do_I(I,j)) domore = .true.
#endif
      else
        do_I(I,j) = .false.
      endif
    endif ; enddo ; enddo
#ifndef _OPENMP
    ! Without OpenMP, stop as soon as every point has converged.  Iterations after that would not
    ! change any of the results, so with OpenMP, where the points are spread across teams that
    ! cannot share this flag, they are simply carried out.
    if (.not.domore) exit
#endif

    !$omp distribute parallel do collapse(2)
    do j=j_start,j_end ; do I=i_start,i_end
      uh_err(I,j) = -uhbt(I,j) ; duhdu_tot(I,j) = 0.0
    enddo ; enddo
    do k=1,nz
      !$omp distribute parallel do collapse(2) private(u_new, duhdu)
      do j=j_start,j_end ; do I=i_start,i_end ; if (do_I(I,j)) then
        u_new = u(I,j,k) + du(I,j) * visc_rem(I,j,k)
        call flux_elem(u_new, h_in(i,j,k), h_in(i+1,j,k), h_W(i,j,k), h_W(i+1,j,k), h_E(i,j,k), &
                       h_E(i+1,j,k), uh_3d(I,j,k), duhdu, visc_rem(I,j,k), G%dy_Cu(I,j), &
                       G%IareaT(i,j), G%IareaT(i+1,j), G%IdxT(i,j), G%IdxT(i+1,j), dt, vol_CFL, &
                       por_face_areaU(I,j,k), h_marg_min)
        if (local_open_BC) &
          call flux_elem_OBC(u_new, h_in(i,j,k), h_in(i+1,j,k), uh_3d(I,j,k), duhdu, &
                             visc_rem(I,j,k), por_face_areaU(I,j,k), G%dy_Cu(I,j), h_marg_min, &
                             open_dir(I,j))
        uh_err(I,j) = uh_err(I,j) + uh_3d(I,j,k)
        duhdu_tot(I,j) = duhdu_tot(I,j) + duhdu
      endif ; enddo ; enddo
    enddo
    !$omp distribute parallel do collapse(2)
    do j=j_start,j_end ; do I=i_start,i_end
      uh_err_best(I,j) = min(uh_err_best(I,j), abs(uh_err(I,j)))
    enddo ; enddo
  enddo ! itt-loop
  ! If there are any faces which have not converged to within the tolerance,
  ! so-be-it, or else use a final upwind correction?
  ! This never seems to happen with 20 iterations as max_itt.
  !$omp end target teams

  !$omp target exit data map(release: uh_err, uh_err_best, duhdu_tot, du_min, du_max, do_I)

end subroutine zonal_flux_adjust

!> Sets a structure that describes the zonal barotropic volume or mass fluxes as a
!! function of barotropic flow to agree closely with the sum of the layer's transports.
subroutine set_zonal_BT_cont(u, h_in, h_W, h_E, BT_cont, uh_tot_0, duhdu_tot_0, &
                             du_max_CFL, du_min_CFL, dt, G, GV, US, CS, visc_rem, &
                             visc_rem_max, i_start, i_end, j_start, j_end, do_I, por_face_areaU, &
                             open_dir)
  type(ocean_grid_type),                     intent(in)    :: G    !< Ocean's grid structure.
  type(verticalGrid_type),                   intent(in)    :: GV   !< Ocean's vertical grid structure.
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), intent(in)   :: u    !< Zonal velocity [L T-1 ~> m s-1].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in)    :: h_in !< Layer thickness used to
                                                                   !! calculate fluxes [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in)    :: h_W  !< West edge thickness in the
                                                                   !! reconstruction [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in)    :: h_E  !< East edge thickness in the
                                                                   !! reconstruction [H ~> m or kg m-2].
  type(BT_cont_type),                        intent(inout) :: BT_cont !< A structure with elements
                       !! that describe the effective open face areas as a function of barotropic flow.
  real, dimension(SZIB_(G),SZJ_(G)),         intent(in)    :: uh_tot_0    !< The summed transport
                       !! with 0 adjustment [H L2 T-1 ~> m3 s-1 or kg s-1].
  real, dimension(SZIB_(G),SZJ_(G)),         intent(in)    :: duhdu_tot_0 !< The partial derivative
                       !! of du_err with du at 0 adjustment [H L ~> m2 or kg m-1].
  real, dimension(SZIB_(G),SZJ_(G)),         intent(in)    :: du_max_CFL  !< Maximum acceptable
                       !! value of du [L T-1 ~> m s-1].
  real, dimension(SZIB_(G),SZJ_(G)),         intent(in)    :: du_min_CFL  !< Minimum acceptable
                       !! value of du [L T-1 ~> m s-1].
  real,                                      intent(in)    :: dt   !< Time increment [T ~> s].
  type(unit_scale_type),                     intent(in)    :: US   !< A dimensional unit scaling type
  type(continuity_PPM_CS),                   intent(in)    :: CS   !< This module's control structure.
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), intent(in)   :: visc_rem !< Both the fraction of the
                       !! momentum originally in a layer that remains after a time-step of viscosity, and
                       !! the fraction of a time-step's worth of a barotropic acceleration that a layer
                       !! experiences after viscosity is applied [nondim].
                       !! Visc_rem is between 0 (at the bottom) and 1 (far above the bottom).
  real, dimension(SZIB_(G),SZJ_(G)),         intent(in)    :: visc_rem_max !< Maximum allowable
                       !! visc_rem [nondim].
  integer,                                   intent(in)    :: i_start !< Start of the I index range
  integer,                                   intent(in)    :: i_end   !< End of the I index range
  integer,                                   intent(in)    :: j_start !< Start of the j index range
  integer,                                   intent(in)    :: j_end   !< End of the j index range
  logical, dimension(SZIB_(G),SZJ_(G)),      intent(in)    :: do_I    !< A logical flag indicating
                       !! which points to work on.
  real, dimension(SZIB_(G),SZJ_(G),SZK_(G)), intent(in)    :: por_face_areaU !< fractional open area
                       !! of U-faces [nondim]
  integer, dimension(SZIB_(G),SZJ_(G)),      intent(in)    :: open_dir !< An array that is passed
                       !! on to zonal_flux_adjust, but is not used there by this call [nondim]
  ! Local variables
  real, dimension(SZIB_(G),SZJ_(G)) :: &
    du0, &        ! The barotropic velocity increment that gives 0 transport [L T-1 ~> m s-1].
    zeros, &      ! An array of full of 0 transports [H L2 T-1 ~> m3 s-1 or kg s-1]
    duL, duR, &   ! The barotropic velocity increments that give the westerly
                  ! (duL) and easterly (duR) test velocities [L T-1 ~> m s-1].
    du_CFL, &     ! The velocity increment that corresponds to CFL_min [L T-1 ~> m s-1].
    FAmt_L, FAmt_R, & ! The summed effective marginal face areas for the 3
    FAmt_0, &     ! test velocities [H L ~> m2 or kg m-1].
    uhtot_L, &    ! The summed transport with the westerly (uhtot_L) and
    uhtot_R       ! and easterly (uhtot_R) test velocities [H L2 T-1 ~> m3 s-1 or kg s-1].
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)) :: &
    uh_tmp        ! The layer transports found while diagnosing du0, which are not used
                  ! [H L2 T-1 ~> m3 s-1 or kg s-1].
  real :: &
    u_L, u_R, &   ! The westerly (u_L), easterly (u_R), and zero-barotropic
    u_0, &        ! transport (u_0) layer test velocities [L T-1 ~> m s-1].
    duhdu_L, &    ! The effective layer marginal face areas with the westerly
    duhdu_R, &    ! (_L), easterly (_R), and zero-barotropic (_0) test
    duhdu_0, &    ! velocities [H L ~> m2 or kg m-1].
    uh_L, uh_R, & ! The layer transports with the westerly (_L), easterly (_R),
    uh_0          ! and zero-barotropic (_0) test velocities [H L2 T-1 ~> m3 s-1 or kg s-1].
  real :: FA_0    ! The effective face area with 0 barotropic transport [L H ~> m2 or kg m-1].
  real :: FA_avg  ! The average effective face area [L H ~> m2 or kg m-1], nominally given by
                  ! the realized transport divided by the barotropic velocity.
  real :: visc_rem_lim ! The larger of visc_rem and min_visc_rem [nondim]. This
                       ! limiting is necessary to keep the inverse of visc_rem
                       ! from leading to large CFL numbers.
  real :: min_visc_rem ! The smallest permitted value for visc_rem that is used
                       ! in finding the barotropic velocity that changes the
                       ! flow direction [nondim].  This is necessary to keep the inverse
                       ! of visc_rem from leading to large CFL numbers.
  real :: CFL_min ! A minimal increment in the CFL to try to ensure that the
                  ! flow is truly upwind [nondim]
  real :: Idt     ! The inverse of the time step [T-1 ~> s-1].
  real :: h_marg_min ! A copy of CS%h_marg_min for use on the device [H ~> m or kg m-2]
  logical :: vol_CFL ! A copy of CS%vol_CFL for use on the device
  integer :: i, j, k, nz

  nz = GV%ke ; Idt = 1.0 / dt
  min_visc_rem = 0.1 ; CFL_min = 1e-6
  vol_CFL = CS%vol_CFL ; h_marg_min = CS%h_marg_min

  !$omp target enter data map(alloc: du0, zeros, duL, duR, du_CFL, FAmt_L, FAmt_R, FAmt_0, &
  !$omp                              uhtot_L, uhtot_R, uh_tmp)

  ! Diagnose the zero-transport correction, du0.
  !$omp target teams distribute parallel do collapse(2)
  do j=j_start,j_end ; do I=i_start,i_end
    zeros(I,j) = 0.0
  enddo ; enddo
  call zonal_flux_adjust(u, h_in, h_W, h_E, zeros, uh_tot_0, duhdu_tot_0, du0, &
                         du_max_CFL, du_min_CFL, dt, G, GV, US, CS, visc_rem, &
                         i_start, i_end, j_start, j_end, do_I, por_face_areaU, uh_tmp, &
                         .false., open_dir)

  ! Determine the westerly- and easterly- fluxes.  Choose a sufficiently
  ! negative velocity correction for the easterly-flux, and a sufficiently
  ! positive correction for the westerly-flux.
  !$omp target teams
  !$omp distribute parallel do collapse(2)
  do j=j_start,j_end ; do I=i_start,i_end
    du_CFL(I,j) = (CFL_min * Idt) * G%dxCu(I,j)
    duR(I,j) = min(0.0,du0(I,j) - du_CFL(I,j))
    duL(I,j) = max(0.0,du0(I,j) + du_CFL(I,j))
    FAmt_L(I,j) = 0.0 ; FAmt_R(I,j) = 0.0 ; FAmt_0(I,j) = 0.0
    uhtot_L(I,j) = 0.0 ; uhtot_R(I,j) = 0.0
  enddo ; enddo

  do k=1,nz
    !$omp distribute parallel do collapse(2) private(visc_rem_lim)
    do j=j_start,j_end ; do I=i_start,i_end ; if (do_I(I,j)) then
      visc_rem_lim = max(visc_rem(I,j,k), min_visc_rem*visc_rem_max(I,j))
      if (visc_rem_lim > 0.0) then ! This is almost always true for ocean points.
        if (u(I,j,k) + duR(I,j)*visc_rem_lim > -du_CFL(I,j)*visc_rem(I,j,k)) &
          duR(I,j) = -(u(I,j,k) + du_CFL(I,j)*visc_rem(I,j,k)) / visc_rem_lim
        if (u(I,j,k) + duL(I,j)*visc_rem_lim < du_CFL(I,j)*visc_rem(I,j,k)) &
          duL(I,j) = -(u(I,j,k) - du_CFL(I,j)*visc_rem(I,j,k)) / visc_rem_lim
      endif
    endif ; enddo ; enddo
  enddo

  do k=1,nz
    !$omp distribute parallel do collapse(2) &
    !$omp   private(u_L, u_R, u_0, duhdu_0, duhdu_L, duhdu_R, uh_0, uh_L, uh_R)
    do j=j_start,j_end ; do I=i_start,i_end ; if (do_I(I,j)) then
      u_L = u(I,j,k) + duL(I,j) * visc_rem(I,j,k)
      u_R = u(I,j,k) + duR(I,j) * visc_rem(I,j,k)
      u_0 = u(I,j,k) + du0(I,j) * visc_rem(I,j,k)
      call flux_elem(u_0, h_in(i,j,k), h_in(i+1,j,k), h_W(i,j,k), h_W(i+1,j,k), h_E(i,j,k), &
                     h_E(i+1,j,k), uh_0, duhdu_0, visc_rem(I,j,k), G%dy_Cu(I,j), &
                     G%IareaT(i,j), G%IareaT(i+1,j), G%IdxT(i,j), G%IdxT(i+1,j), dt, vol_CFL, &
                     por_face_areaU(I,j,k), h_marg_min)
      call flux_elem(u_L, h_in(i,j,k), h_in(i+1,j,k), h_W(i,j,k), h_W(i+1,j,k), h_E(i,j,k), &
                     h_E(i+1,j,k), uh_L, duhdu_L, visc_rem(I,j,k), G%dy_Cu(I,j), &
                     G%IareaT(i,j), G%IareaT(i+1,j), G%IdxT(i,j), G%IdxT(i+1,j), dt, vol_CFL, &
                     por_face_areaU(I,j,k), h_marg_min)
      call flux_elem(u_R, h_in(i,j,k), h_in(i+1,j,k), h_W(i,j,k), h_W(i+1,j,k), h_E(i,j,k), &
                     h_E(i+1,j,k), uh_R, duhdu_R, visc_rem(I,j,k), G%dy_Cu(I,j), &
                     G%IareaT(i,j), G%IareaT(i+1,j), G%IdxT(i,j), G%IdxT(i+1,j), dt, vol_CFL, &
                     por_face_areaU(I,j,k), h_marg_min)
      FAmt_0(I,j) = FAmt_0(I,j) + duhdu_0
      FAmt_L(I,j) = FAmt_L(I,j) + duhdu_L
      FAmt_R(I,j) = FAmt_R(I,j) + duhdu_R
      uhtot_L(I,j) = uhtot_L(I,j) + uh_L
      uhtot_R(I,j) = uhtot_R(I,j) + uh_R
    endif ; enddo ; enddo
  enddo

  !$omp distribute parallel do collapse(2) private(FA_0, FA_avg)
  do j=j_start,j_end ; do I=i_start,i_end
    if (do_I(I,j)) then
      FA_0 = FAmt_0(I,j) ; FA_avg = FAmt_0(I,j)
      if ((duL(I,j) - du0(I,j)) /= 0.0) &
        FA_avg = uhtot_L(I,j) / (duL(I,j) - du0(I,j))
      if (FA_avg > max(FA_0, FAmt_L(I,j))) then ; FA_avg = max(FA_0, FAmt_L(I,j))
      elseif (FA_avg < min(FA_0, FAmt_L(I,j))) then ; FA_0 = FA_avg ; endif

      BT_cont%FA_u_W0(I,j) = FA_0 ; BT_cont%FA_u_WW(I,j) = FAmt_L(I,j)
      if (abs(FA_0-FAmt_L(I,j)) <= 1e-12*FA_0) then ; BT_cont%uBT_WW(I,j) = 0.0 ; else
        BT_cont%uBT_WW(I,j) = (1.5 * (duL(I,j) - du0(I,j))) * &
                              ((FAmt_L(I,j) - FA_avg) / (FAmt_L(I,j) - FA_0))
      endif

      FA_0 = FAmt_0(I,j) ; FA_avg = FAmt_0(I,j)
      if ((duR(I,j) - du0(I,j)) /= 0.0) &
        FA_avg = uhtot_R(I,j) / (duR(I,j) - du0(I,j))
      if (FA_avg > max(FA_0, FAmt_R(I,j))) then ; FA_avg = max(FA_0, FAmt_R(I,j))
      elseif (FA_avg < min(FA_0, FAmt_R(I,j))) then ; FA_0 = FA_avg ; endif

      BT_cont%FA_u_E0(I,j) = FA_0 ; BT_cont%FA_u_EE(I,j) = FAmt_R(I,j)
      if (abs(FAmt_R(I,j) - FA_0) <= 1e-12*FA_0) then ; BT_cont%uBT_EE(I,j) = 0.0 ; else
        BT_cont%uBT_EE(I,j) = (1.5 * (duR(I,j) - du0(I,j))) * &
                              ((FAmt_R(I,j) - FA_avg) / (FAmt_R(I,j) - FA_0))
      endif
    else
      BT_cont%FA_u_W0(I,j) = 0.0 ; BT_cont%FA_u_WW(I,j) = 0.0
      BT_cont%FA_u_E0(I,j) = 0.0 ; BT_cont%FA_u_EE(I,j) = 0.0
      BT_cont%uBT_WW(I,j) = 0.0 ; BT_cont%uBT_EE(I,j) = 0.0
    endif
  enddo ; enddo
  !$omp end target teams

  !$omp target exit data map(release: du0, zeros, duL, duR, du_CFL, FAmt_L, FAmt_R, FAmt_0, &
  !$omp                               uhtot_L, uhtot_R, uh_tmp)

end subroutine set_zonal_BT_cont

!> Calculates the mass or volume flux through a single zonal or meridional face and its partial
!! derivative with the face velocity, from a PPM reconstruction of the thicknesses on either side.
subroutine flux_elem(u, h, h_p1, h_L, h_L_p1, h_R, h_R_p1, uh, duhdu, visc_rem, &
                     G_dy_Cu, G_IareaT, G_IareaT_p1, G_IdxT, G_IdxT_p1, dt, &
                     vol_CFL, por_face_area, h_marg_min)
  !$omp declare target
  real,    intent(in)  :: u        !< Zonal or meridional velocity [L T-1 ~> m s-1].
  real,    intent(in)  :: h        !< Layer thickness [H ~> m or kg m-2].
  real,    intent(in)  :: h_p1     !< Layer thickness, offset by 1 [H ~> m or kg m-2].
  real,    intent(in)  :: h_L      !< West or south edge thickness [H ~> m or kg m-2].
  real,    intent(in)  :: h_L_p1   !< West or south edge thickness, offset by 1 [H ~> m or kg m-2].
  real,    intent(in)  :: h_R      !< East or north edge thickness [H ~> m or kg m-2].
  real,    intent(in)  :: h_R_p1   !< East or north edge thickness, offset by 1 [H ~> m or kg m-2].
  real,    intent(out) :: uh       !< Zonal or meridional mass or volume transport
                                   !! [H L2 T-1 ~> m3 s-1 or kg s-1].
  real,    intent(out) :: duhdu    !< Partial derivative of uh with u [H L ~> m2 or kg m-1].
  real,    intent(in)  :: visc_rem !< Both the fraction of the momentum originally in a layer that
                        !! remains after a time-step of viscosity, and the fraction of a time-step's
                        !! worth of a barotropic acceleration that a layer experiences after viscosity
                        !! is applied [nondim].  Visc_rem is between 0 (at the bottom) and 1 (far above).
  real,    intent(in)  :: G_dy_Cu  !< The unblocked length of the face [L ~> m].
  real,    intent(in)  :: G_IareaT !< The inverse of the area of the cell on the low side
                                   !! of the face [L-2 ~> m-2].
  real,    intent(in)  :: G_IareaT_p1 !< The inverse of the area of the cell on the high side
                                   !! of the face [L-2 ~> m-2].
  real,    intent(in)  :: G_IdxT   !< The inverse of the grid spacing of the cell on the low side
                                   !! of the face [L-1 ~> m-1].
  real,    intent(in)  :: G_IdxT_p1 !< The inverse of the grid spacing of the cell on the high side
                                   !! of the face [L-1 ~> m-1].
  real,    intent(in)  :: dt       !< Time increment [T ~> s].
  logical, intent(in)  :: vol_CFL  !< If true, rescale the ratio of face areas to the cell areas
                                   !! when estimating the CFL number.
  real,    intent(in)  :: por_face_area !< fractional open area of the face [nondim].
  real,    intent(in)  :: h_marg_min !< Negligible floor on h_marg [H ~> m or kg m-2]
  ! Local variables
  real :: CFL    ! The CFL number based on the local velocity and grid spacing [nondim]
  real :: curv_3 ! A measure of the thickness curvature over a grid length [H ~> m or kg m-2]
  real :: h_marg ! The marginal thickness of a flux [H ~> m or kg m-2].
  real :: dy_por ! The open length of the face [L ~> m].
  real :: dh     ! The difference between the edge thicknesses in the upwind cell [H ~> m or kg m-2].

  dy_por = G_dy_Cu * por_face_area
  if (u > 0.0) then
    if (vol_CFL) then ; CFL = (u * dt) * (G_dy_Cu * G_IareaT)
    else ; CFL = u * dt * G_IdxT ; endif
    curv_3 = (h_L + h_R) - 2.0*h
    dh = h_L - h_R
    uh = dy_por * u * (h_R + CFL * (0.5*dh + curv_3*(CFL - 1.5)))
    h_marg = h_R + CFL * (dh + 3.0*curv_3*(CFL - 1.0))
  elseif (u < 0.0) then
    if (vol_CFL) then ; CFL = (-u * dt) * (G_dy_Cu * G_IareaT_p1)
    else ; CFL = -u * dt * G_IdxT_p1 ; endif
    curv_3 = (h_L_p1 + h_R_p1) - 2.0*h_p1
    dh = h_R_p1 - h_L_p1
    uh = dy_por * u * (h_L_p1 + CFL * (0.5*dh + curv_3*(CFL - 1.5)))
    h_marg = h_L_p1 + CFL * (dh + 3.0*curv_3*(CFL - 1.0))
  else
    uh = 0.0
    h_marg = 0.5 * (h_L_p1 + h_R)
  endif
  h_marg = max(h_marg, h_marg_min)
  duhdu = dy_por * h_marg * visc_rem

end subroutine flux_elem

!> Replaces the transport through a single zonal or meridional face and its partial derivative
!! with the face velocity with simple upwind estimates if the face is on an open boundary.
subroutine flux_elem_OBC(u, h, h_p1, uh, duhdu, visc_rem, por_face_area, G_dy_Cu, h_marg_min, open_dir)
  !$omp declare target
  real,                 intent(in)    :: u        !< Zonal or meridional velocity [L T-1 ~> m s-1].
  real,                 intent(in)    :: h        !< Layer thickness [H ~> m or kg m-2].
  real,                 intent(in)    :: h_p1     !< Layer thickness, offset by 1 [H ~> m or kg m-2].
  real,                 intent(inout) :: uh       !< Zonal or meridional mass or volume transport
                                                  !! [H L2 T-1 ~> m3 s-1 or kg s-1].
  real,                 intent(inout) :: duhdu    !< Partial derivative of uh with u [H L ~> m2 or kg m-1].
  real,                 intent(in)    :: visc_rem !< Both the fraction of the momentum originally in a
                        !! layer that remains after a time-step of viscosity, and the fraction of a
                        !! time-step's worth of a barotropic acceleration that a layer experiences after
                        !! viscosity is applied [nondim].
  real,                 intent(in)    :: por_face_area !< fractional open area of the face [nondim].
  real,                 intent(in)    :: G_dy_Cu  !< The unblocked length of the face [L ~> m].
  real,                 intent(in)    :: h_marg_min !< Negligible floor on h_marg [H ~> m or kg m-2]
  integer,              intent(in)    :: open_dir !< 1 or -1 if the face is on an open boundary
                                                  !! segment where the flow is taken from the cell on
                                                  !! the low or high side of the face, or 0 [nondim]

  if (open_dir > 0) then !  OBC_DIRECTION_E or OBC_DIRECTION_N
    uh = (G_dy_Cu * por_face_area) * u * h
    duhdu = (G_dy_Cu * por_face_area) * max(h, h_marg_min) * visc_rem
  elseif (open_dir < 0) then !  OBC_DIRECTION_W or OBC_DIRECTION_S
    uh = (G_dy_Cu * por_face_area) * u * h_p1
    duhdu = (G_dy_Cu * por_face_area) * max(h_p1, h_marg_min) * visc_rem
  endif

end subroutine flux_elem_OBC

module procedure meridional_mass_flux
  ! Local variables
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)) :: &
    dvhdv, &      ! Partial derivative of vh with v [H L ~> m2 or kg m-1].
    visc_rem      ! A copy of visc_rem_v or an array of 1's [nondim].
  real, dimension(SZI_(G),SZJB_(G)) :: &
    dv, &         ! Corrective barotropic change in the velocity to give vhbt [L T-1 ~> m s-1].
    dv_min_CFL, & ! Lower limit on dv correction to avoid CFL violations [L T-1 ~> m s-1]
    dv_max_CFL, & ! Upper limit on dv correction to avoid CFL violations [L T-1 ~> m s-1]
    dvhdv_tot_0, & ! Summed partial derivative of vh with v [H L ~> m2 or kg m-1].
    vh_tot_0, &   ! Summed transport with no barotropic correction [H L2 T-1 ~> m3 s-1 or kg s-1].
    visc_rem_max, & ! The column maximum of visc_rem [nondim]
    FAvi          ! The sums of meridional face areas [H L ~> m2 or kg m-1].
  logical, dimension(SZI_(G),SZJB_(G)) :: &
    do_I, &       ! Indicates the points where the barotropic and baroclinic transports are reconciled
    simple_OBC_pt ! Indicates points with specified transport OBCs
  integer, dimension(SZI_(G),SZJB_(G)) :: &
    open_dir      ! 1 or -1 at faces on open boundary segments where the flow is taken from the
                  ! cell to the south or north of the face, or 0 elsewhere [nondim]
  real :: FA_v    ! A sum of meridional face areas [H L ~> m2 or kg m-1].
  real :: I_vrm   ! 1.0 / visc_rem_max [nondim]
  real :: CFL_dt  ! The maximum CFL ratio of the adjusted velocities divided by
                  ! the time step [T-1 ~> s-1].
  real :: I_dt    ! 1.0 / dt [T-1 ~> s-1].
  real :: dv_lim  ! The velocity change that give a relative CFL of 1 [L T-1 ~> m s-1].
  real :: dy_N, dy_S ! Effective y-grid spacings to the north and south [L ~> m].
  real :: h_marg_min ! A copy of CS%h_marg_min for use on the device [H ~> m or kg m-2]
  real :: H_subroundoff ! A copy of GV%H_subroundoff for use on the device [H ~> m or kg m-2]
  type(cont_loop_bounds_type) :: LB
  integer :: i, j, k, ish, ieh, jsh, jeh, n, nz
  integer :: isd, ied, JsdB, JedB ! The data domain bounds at v points
  integer :: l_seg ! The OBC segment number
  logical :: use_visc_rem, set_BT_cont, set_h_v
  logical :: vol_CFL, aggress_adjust, use_visc_rem_max ! Copies of CS fields for use on the device
  logical :: local_specified_BC, local_Flather_OBC, local_open_BC, any_simple_OBC  ! OBC-related logicals

  use_visc_rem = present(visc_rem_v)

  set_BT_cont = .false. ; if (present(BT_cont)) set_BT_cont = (associated(BT_cont))
  set_h_v = .false. ; if (set_BT_cont) set_h_v = allocated(BT_cont%h_v)

  !   The arrays that are used on the device are copied in before the clock is started.  The output
  ! arrays that are only partly set here are copied in too, so that the copies back to the host
  ! leave the rest of them unchanged.
  !$omp target enter data map(to: G)
  !$omp target enter data map(to: G%dx_Cv, G%dyCv, G%IareaT, G%areaT, G%IdyT, G%dyT, G%mask2dCv)
  !$omp target enter data map(to: v, h_in, h_S, h_N, por_face_areaV, vh)
  if (use_visc_rem) then
    !$omp target enter data map(to: visc_rem_v)
  endif
  if (present(vhbt)) then
    !$omp target enter data map(to: vhbt)
  endif
  if (present(v_cor)) then
    !$omp target enter data map(to: v_cor)
  endif
  if (present(dv_cor)) then
    !$omp target enter data map(alloc: dv_cor)
  endif
  if (set_BT_cont) call merid_BT_cont_to_device(BT_cont, set_h_v)

  call cpu_clock_begin(id_clock_correct)

  !$omp target enter data map(alloc: dvhdv, visc_rem, dv, dv_min_CFL, dv_max_CFL, dvhdv_tot_0, &
  !$omp                              vh_tot_0, visc_rem_max, FAvi, do_I, simple_OBC_pt, open_dir)

  local_specified_BC = .false. ; local_Flather_OBC = .false. ; local_open_BC = .false.
  if (associated(OBC)) then ; if (OBC%OBC_pe) then
    local_specified_BC = OBC%specified_v_BCs_exist_globally
    local_Flather_OBC = OBC%Flather_v_BCs_exist_globally
    local_open_BC = OBC%open_v_BCs_exist_globally
  endif ; endif

  if (present(LB_in)) then
    LB = LB_in
  else
    LB%ish = G%isc ; LB%ieh = G%iec ; LB%jsh = G%jsc ; LB%jeh = G%jec
  endif
  ish = LB%ish ; ieh = LB%ieh ; jsh = LB%jsh ; jeh = LB%jeh ; nz = GV%ke
  isd = G%isd ; ied = G%ied ; JsdB = G%JsdB ; JedB = G%JedB

  vol_CFL = CS%vol_CFL ; aggress_adjust = CS%aggress_adjust ; use_visc_rem_max = CS%use_visc_rem_max
  h_marg_min = CS%h_marg_min ; H_subroundoff = GV%H_subroundoff

  if (present(dv_cor)) then
    !$omp target teams distribute parallel do collapse(2)
    do J=JsdB,JedB ; do i=isd,ied
      dv_cor(i,J) = 0.0
    enddo ; enddo
  endif

  CFL_dt = CS%CFL_limit_adjust / dt
  I_dt = 1.0 / dt
  if (aggress_adjust) CFL_dt = I_dt

  if (local_open_BC) then
    ! Note which faces are on open boundary segments, so that OBC need not be used on the device.
    do J=jsh-1,jeh ; do i=ish,ieh
      open_dir(i,J) = 0
      if (OBC%segnum_v(i,J) /= 0) then
        if (OBC%segment(abs(OBC%segnum_v(i,J)))%open) open_dir(i,J) = sign(1, OBC%segnum_v(i,J))
      endif
    enddo ; enddo
    !$omp target update to(open_dir)
  endif

  ! This sets vh and dvhdv.
  if (use_visc_rem) then
    !$omp target teams distribute parallel do collapse(3)
    do k=1,nz ; do J=jsh-1,jeh ; do i=ish,ieh
      visc_rem(i,J,k) = visc_rem_v(i,J,k)
    enddo ; enddo ; enddo
  else
    !$omp target teams distribute parallel do collapse(3)
    do k=1,nz ; do J=jsh-1,jeh ; do i=ish,ieh
      visc_rem(i,J,k) = 1.0
    enddo ; enddo ; enddo
  endif
  !$omp target teams distribute parallel do collapse(3)
  do k=1,nz ; do J=jsh-1,jeh ; do i=ish,ieh
    call flux_elem(v(i,J,k), h_in(i,j,k), h_in(i,j+1,k), h_S(i,j,k), h_S(i,j+1,k), h_N(i,j,k), &
                   h_N(i,j+1,k), vh(i,J,k), dvhdv(i,J,k), visc_rem(i,J,k), G%dx_Cv(i,J), &
                   G%IareaT(i,j), G%IareaT(i,j+1), G%IdyT(i,j), G%IdyT(i,j+1), dt, vol_CFL, &
                   por_face_areaV(i,J,k), h_marg_min)
    if (local_open_BC) &
      call flux_elem_OBC(v(i,J,k), h_in(i,j,k), h_in(i,j+1,k), vh(i,J,k), dvhdv(i,J,k), &
                         visc_rem(i,J,k), por_face_areaV(i,J,k), G%dx_Cv(i,J), h_marg_min, &
                         open_dir(i,J))
  enddo ; enddo ; enddo
  if (local_specified_BC) then
    !$omp target teams distribute parallel do collapse(3) private(l_seg)
    do k=1,nz ; do J=jsh-1,jeh ; do i=ish,ieh ; if (OBC%segnum_v(i,J) /= 0) then
      l_seg = abs(OBC%segnum_v(i,J))
      if (OBC%segment(l_seg)%specified) vh(i,J,k) = OBC%segment(l_seg)%normal_trans(i,J,k)
    endif ; enddo ; enddo ; enddo
  endif

  if (present(vhbt) .or. set_BT_cont) then
    !$omp target teams
    if (use_visc_rem .and. use_visc_rem_max) then
      !$omp distribute parallel do collapse(2)
      do J=jsh-1,jeh ; do i=ish,ieh
        visc_rem_max(i,J) = 0.0
      enddo ; enddo
      do k=1,nz
        !$omp distribute parallel do collapse(2)
        do J=jsh-1,jeh ; do i=ish,ieh
          visc_rem_max(i,J) = max(visc_rem_max(i,J), visc_rem(i,J,k))
        enddo ; enddo
      enddo
    else
      !$omp distribute parallel do collapse(2)
      do J=jsh-1,jeh ; do i=ish,ieh
        visc_rem_max(i,J) = 1.0
      enddo ; enddo
    endif
    !   Set limits on dv that will keep the CFL number between -1 and 1.
    ! This should be adequate to keep the root bracketed in all cases.
    !$omp distribute parallel do collapse(2) private(I_vrm, dy_S, dy_N)
    do J=jsh-1,jeh ; do i=ish,ieh
      I_vrm = 0.0
      if (visc_rem_max(i,J) > 0.0) I_vrm = 1.0 / visc_rem_max(i,J)
      if (vol_CFL) then
        dy_S = ratio_max(G%areaT(i,j), G%dx_Cv(i,J), 1000.0*G%dyT(i,j))
        dy_N = ratio_max(G%areaT(i,j+1), G%dx_Cv(i,J), 1000.0*G%dyT(i,j+1))
      else ; dy_S = G%dyT(i,j) ; dy_N = G%dyT(i,j+1) ; endif
      dv_max_CFL(i,J) = 2.0 * (CFL_dt * dy_S) * I_vrm
      dv_min_CFL(i,J) = -2.0 * (CFL_dt * dy_N) * I_vrm
      vh_tot_0(i,J) = 0.0 ; dvhdv_tot_0(i,J) = 0.0
    enddo ; enddo
    do k=1,nz
      !$omp distribute parallel do collapse(2)
      do J=jsh-1,jeh ; do i=ish,ieh
        dvhdv_tot_0(i,J) = dvhdv_tot_0(i,J) + dvhdv(i,J,k)
        vh_tot_0(i,J) = vh_tot_0(i,J) + vh(i,J,k)
      enddo ; enddo
      if (use_visc_rem) then
        if (aggress_adjust) then
          !$omp distribute parallel do collapse(2) private(dy_S, dy_N, dv_lim)
          do J=jsh-1,jeh ; do i=ish,ieh
            if (vol_CFL) then
              dy_S = ratio_max(G%areaT(i,j), G%dx_Cv(i,J), 1000.0*G%dyT(i,j))
              dy_N = ratio_max(G%areaT(i,j+1), G%dx_Cv(i,J), 1000.0*G%dyT(i,j+1))
            else ; dy_S = G%dyT(i,j) ; dy_N = G%dyT(i,j+1) ; endif
            dv_lim = 0.499*((dy_S*I_dt - v(i,J,k)) + MIN(0.0,v(i,J-1,k)))
            if (dv_max_CFL(i,J) * visc_rem(i,J,k) > dv_lim) &
              dv_max_CFL(i,J) = dv_lim / visc_rem(i,J,k)

            dv_lim = 0.499*((-dy_N*CFL_dt - v(i,J,k)) + MAX(0.0,v(i,J+1,k)))
            if (dv_min_CFL(i,J) * visc_rem(i,J,k) < dv_lim) &
              dv_min_CFL(i,J) = dv_lim / visc_rem(i,J,k)
          enddo ; enddo
        else
          !$omp distribute parallel do collapse(2) private(dy_S, dy_N)
          do J=jsh-1,jeh ; do i=ish,ieh
            if (vol_CFL) then
              dy_S = ratio_max(G%areaT(i,j), G%dx_Cv(i,J), 1000.0*G%dyT(i,j))
              dy_N = ratio_max(G%areaT(i,j+1), G%dx_Cv(i,J), 1000.0*G%dyT(i,j+1))
            else ; dy_S = G%dyT(i,j) ; dy_N = G%dyT(i,j+1) ; endif
            if (dv_max_CFL(i,J) * visc_rem(i,J,k) > dy_S*CFL_dt - v(i,J,k)*G%mask2dCv(i,J)) &
              dv_max_CFL(i,J) = (dy_S*CFL_dt - v(i,J,k)) / visc_rem(i,J,k)
            if (dv_min_CFL(i,J) * visc_rem(i,J,k) < -dy_N*CFL_dt - v(i,J,k)*G%mask2dCv(i,J)) &
              dv_min_CFL(i,J) = -(dy_N*CFL_dt + v(i,J,k)) / visc_rem(i,J,k)
          enddo ; enddo
        endif
      else
        if (aggress_adjust) then
          !$omp distribute parallel do collapse(2) private(dy_S, dy_N)
          do J=jsh-1,jeh ; do i=ish,ieh
            if (vol_CFL) then
              dy_S = ratio_max(G%areaT(i,j), G%dx_Cv(i,J), 1000.0*G%dyT(i,j))
              dy_N = ratio_max(G%areaT(i,j+1), G%dx_Cv(i,J), 1000.0*G%dyT(i,j+1))
            else ; dy_S = G%dyT(i,j) ; dy_N = G%dyT(i,j+1) ; endif
            dv_max_CFL(i,J) = min(dv_max_CFL(i,J), 0.499 * &
                        ((dy_S*I_dt - v(i,J,k)) + MIN(0.0,v(i,J-1,k))) )
            dv_min_CFL(i,J) = max(dv_min_CFL(i,J), 0.499 * &
                        ((-dy_N*I_dt - v(i,J,k)) + MAX(0.0,v(i,J+1,k))) )
          enddo ; enddo
        else
          !$omp distribute parallel do collapse(2) private(dy_S, dy_N)
          do J=jsh-1,jeh ; do i=ish,ieh
            if (vol_CFL) then
              dy_S = ratio_max(G%areaT(i,j), G%dx_Cv(i,J), 1000.0*G%dyT(i,j))
              dy_N = ratio_max(G%areaT(i,j+1), G%dx_Cv(i,J), 1000.0*G%dyT(i,j+1))
            else ; dy_S = G%dyT(i,j) ; dy_N = G%dyT(i,j+1) ; endif
            dv_max_CFL(i,J) = min(dv_max_CFL(i,J), dy_S*CFL_dt - v(i,J,k))
            dv_min_CFL(i,J) = max(dv_min_CFL(i,J), -(dy_N*CFL_dt + v(i,J,k)))
          enddo ; enddo
        endif
      endif
    enddo
    !$omp distribute parallel do collapse(2)
    do J=jsh-1,jeh ; do i=ish,ieh
      dv_max_CFL(i,J) = max(dv_max_CFL(i,J),0.0)
      dv_min_CFL(i,J) = min(dv_min_CFL(i,J),0.0)
    enddo ; enddo
    !$omp end target teams

    any_simple_OBC = .false.
    if (local_specified_BC .or. local_Flather_OBC) then
      !$omp target teams distribute parallel do collapse(2) private(l_seg) &
      !$omp   reduction(.or.: any_simple_OBC) map(tofrom: any_simple_OBC)
      do J=jsh-1,jeh ; do i=ish,ieh
        l_seg = abs(OBC%segnum_v(i,J))

        ! Avoid reconciling barotropic/baroclinic transports if transport is specified
        simple_OBC_pt(i,J) = .false.
        if (l_seg /= 0) simple_OBC_pt(i,J) = OBC%segment(l_seg)%specified
        do_I(i,J) = .not.simple_OBC_pt(i,J)
        any_simple_OBC = any_simple_OBC .or. simple_OBC_pt(i,J)
      enddo ; enddo
    else
      !$omp target teams distribute parallel do collapse(2)
      do J=jsh-1,jeh ; do i=ish,ieh
        do_I(i,J) = .true.
      enddo ; enddo
    endif

    if (present(vhbt)) then
      ! Find dv and vh.
      call meridional_flux_adjust(v, h_in, h_S, h_N, vhbt, vh_tot_0, dvhdv_tot_0, dv, &
                                  dv_max_CFL, dv_min_CFL, dt, G, GV, US, CS, visc_rem, &
                                  ish, ieh, jsh-1, jeh, do_I, por_face_areaV, vh, &
                                  local_open_BC, open_dir)

      if (present(v_cor)) then
        !$omp target teams distribute parallel do collapse(3)
        do k=1,nz ; do J=jsh-1,jeh ; do i=ish,ieh
          v_cor(i,J,k) = v(i,J,k) + dv(i,J) * visc_rem(i,J,k)
        enddo ; enddo ; enddo
        if (any_simple_OBC) then
          !$omp target teams distribute parallel do collapse(3)
          do k=1,nz ; do J=jsh-1,jeh ; do i=ish,ieh ; if (simple_OBC_pt(i,J)) then
            v_cor(i,J,k) = OBC%segment(abs(OBC%segnum_v(i,J)))%normal_vel(i,J,k)
          endif ; enddo ; enddo ; enddo
        endif
      endif ! v-corrected

      if (present(dv_cor)) then
        !$omp target teams distribute parallel do collapse(2)
        do J=jsh-1,jeh ; do i=ish,ieh
          dv_cor(i,J) = dv(i,J)
        enddo ; enddo
      endif

    endif

    if (set_BT_cont) then
      call set_merid_BT_cont(v, h_in, h_S, h_N, BT_cont, vh_tot_0, dvhdv_tot_0, &
                             dv_max_CFL, dv_min_CFL, dt, G, GV, US, CS, visc_rem, &
                             visc_rem_max, ish, ieh, jsh-1, jeh, do_I, por_face_areaV, open_dir)
      if (any_simple_OBC) then
        !$omp target teams
        !$omp distribute parallel do collapse(2)
        do J=jsh-1,jeh ; do i=ish,ieh
          if (simple_OBC_pt(i,J)) FAvi(i,J) = H_subroundoff*G%dx_Cv(i,J)
        enddo ; enddo
        ! NOTE: simple_OBC_pt should prevent access to segment OBC_NONE
        do k=1,nz
          !$omp distribute parallel do collapse(2) private(l_seg)
          do J=jsh-1,jeh ; do i=ish,ieh ; if (simple_OBC_pt(i,J)) then
            l_seg = abs(OBC%segnum_v(i,J))
            if ((abs(OBC%segment(l_seg)%normal_vel(i,J,k)) > 0.0) .and. (OBC%segment(l_seg)%specified)) &
              FAvi(i,J) = FAvi(i,J) + OBC%segment(l_seg)%normal_trans(i,J,k) / OBC%segment(l_seg)%normal_vel(i,J,k)
          endif ; enddo ; enddo
        enddo
        !$omp distribute parallel do collapse(2)
        do J=jsh-1,jeh ; do i=ish,ieh ; if (simple_OBC_pt(i,J)) then
          BT_cont%FA_v_S0(i,J) = FAvi(i,J) ; BT_cont%FA_v_N0(i,J) = FAvi(i,J)
          BT_cont%FA_v_SS(i,J) = FAvi(i,J) ; BT_cont%FA_v_NN(i,J) = FAvi(i,J)
          BT_cont%vBT_SS(i,J) = 0.0 ; BT_cont%vBT_NN(i,J) = 0.0
        endif ; enddo ; enddo
        !$omp end target teams
      endif
    endif ! set_BT_cont

  endif ! present(vhbt) or set_BT_cont

  if (local_open_BC .and. set_BT_cont) then
    ! This is done on the host, between copies of the face areas from and back to the device.
    !$omp target update from(BT_cont%FA_v_S0, BT_cont%FA_v_N0, BT_cont%FA_v_SS, &
    !$omp                    BT_cont%FA_v_NN, BT_cont%vBT_SS, BT_cont%vBT_NN)
    do n = 1, OBC%number_of_segments
      if (OBC%segment(n)%open .and. OBC%segment(n)%is_N_or_S) then
        J = OBC%segment(n)%HI%JsdB
        if (OBC%segment(n)%direction == OBC_DIRECTION_N) then
          do i = OBC%segment(n)%HI%Isd, OBC%segment(n)%HI%Ied
            FA_v = 0.0
            do k=1,nz ; FA_v = FA_v + h_in(i,j,k)*(G%dx_Cv(i,J)*por_face_areaV(i,J,k)) ; enddo
            BT_cont%FA_v_S0(i,J) = FA_v ; BT_cont%FA_v_N0(i,J) = FA_v
            BT_cont%FA_v_SS(i,J) = FA_v ; BT_cont%FA_v_NN(i,J) = FA_v
            BT_cont%vBT_SS(i,J) = 0.0 ; BT_cont%vBT_NN(i,J) = 0.0
          enddo
        else
          do i = OBC%segment(n)%HI%Isd, OBC%segment(n)%HI%Ied
            FA_v = 0.0
            do k=1,nz ; FA_v = FA_v + h_in(i,j+1,k)*(G%dx_Cv(i,J)*por_face_areaV(i,J,k)) ; enddo
            BT_cont%FA_v_S0(i,J) = FA_v ; BT_cont%FA_v_N0(i,J) = FA_v
            BT_cont%FA_v_SS(i,J) = FA_v ; BT_cont%FA_v_NN(i,J) = FA_v
            BT_cont%vBT_SS(i,J) = 0.0 ; BT_cont%vBT_NN(i,J) = 0.0
          enddo
        endif
      endif
    enddo
    !$omp target update to(BT_cont%FA_v_S0, BT_cont%FA_v_N0, BT_cont%FA_v_SS, &
    !$omp                  BT_cont%FA_v_NN, BT_cont%vBT_SS, BT_cont%vBT_NN)
  endif

  if (set_h_v) then
    if (present(v_cor)) then
      call meridional_flux_thickness(v_cor, h_in, h_S, h_N, BT_cont%h_v, dt, G, GV, US, LB, &
                                    vol_CFL, CS%marginal_faces, OBC, por_face_areaV, visc_rem_v)
    else
      call meridional_flux_thickness(v, h_in, h_S, h_N, BT_cont%h_v, dt, G, GV, US, LB, &
                                    vol_CFL, CS%marginal_faces, OBC, por_face_areaV, visc_rem_v)
    endif
  endif

  !$omp target exit data map(release: dvhdv, visc_rem, dv, dv_min_CFL, dv_max_CFL, dvhdv_tot_0, &
  !$omp                               vh_tot_0, visc_rem_max, FAvi, do_I, simple_OBC_pt, open_dir)

  call cpu_clock_end(id_clock_correct)

  ! The results are copied back to the host, and the inputs released, after the clock is stopped.
  if (set_BT_cont) call merid_BT_cont_from_device(BT_cont, set_h_v)
  if (present(dv_cor)) then
    !$omp target exit data map(from: dv_cor)
  endif
  if (present(v_cor)) then
    !$omp target exit data map(from: v_cor)
  endif
  if (present(vhbt)) then
    !$omp target exit data map(release: vhbt)
  endif
  if (use_visc_rem) then
    !$omp target exit data map(release: visc_rem_v)
  endif
  !$omp target exit data map(from: vh) map(release: v, h_in, h_S, h_N, por_face_areaV)
  !$omp target exit data map(release: G%dx_Cv, G%dyCv, G%IareaT, G%areaT, G%IdyT, G%dyT, G%mask2dCv)
  !$omp target exit data map(release: G)

end procedure meridional_mass_flux

module procedure meridional_flux_thickness
  ! Local variables
  real :: CFL ! The CFL number based on the local velocity and grid spacing [nondim]
  real :: curv_3 ! A measure of the thickness curvature over a grid length,
                 ! with the same units as h [H ~> m or kg m-2] .
  real :: h_avg  ! The average thickness of a flux [H ~> m or kg m-2].
  real :: h_marg ! The marginal thickness of a flux [H ~> m or kg m-2].
  logical :: local_open_BC
  integer :: i, j, k, ish, ieh, jsh, jeh, n, nz
  ish = LB%ish ; ieh = LB%ieh ; jsh = LB%jsh ; jeh = LB%jeh ; nz = GV%ke

  !   This works on the device, with all of the arrays already there.
  !$omp target teams distribute parallel do collapse(3) private(CFL, curv_3, h_marg, h_avg)
  do k=1,nz ; do J=jsh-1,jeh ; do i=ish,ieh
    if (v(i,J,k) > 0.0) then
      if (vol_CFL) then ; CFL = (v(i,J,k) * dt) * (G%dx_Cv(i,J) * G%IareaT(i,j))
      else ; CFL = v(i,J,k) * dt * G%IdyT(i,j) ; endif
      curv_3 = (h_S(i,j,k) + h_N(i,j,k)) - 2.0*h(i,j,k)
      h_avg = h_N(i,j,k) + CFL * (0.5*(h_S(i,j,k) - h_N(i,j,k)) + curv_3*(CFL - 1.5))
      h_marg = h_N(i,j,k) + CFL * ((h_S(i,j,k) - h_N(i,j,k)) + &
                                3.0*curv_3*(CFL - 1.0))
    elseif (v(i,J,k) < 0.0) then
      if (vol_CFL) then ; CFL = (-v(i,J,k)*dt) * (G%dx_Cv(i,J) * G%IareaT(i,j+1))
      else ; CFL = -v(i,J,k) * dt * G%IdyT(i,j+1) ; endif
      curv_3 = (h_S(i,j+1,k) + h_N(i,j+1,k)) - 2.0*h(i,j+1,k)
      h_avg = h_S(i,j+1,k) + CFL * (0.5*(h_N(i,j+1,k)-h_S(i,j+1,k)) + curv_3*(CFL - 1.5))
      h_marg = h_S(i,j+1,k) + CFL * ((h_N(i,j+1,k)-h_S(i,j+1,k)) + &
                                    3.0*curv_3*(CFL - 1.0))
    else
      h_avg = 0.5 * (h_S(i,j+1,k) + h_N(i,j,k))
      !   The choice to use the arithmetic mean here is somewhat arbitrarily, but
      ! it should be noted that h_S(i+1,j,k) and h_N(i,j,k) are usually the same.
      h_marg = 0.5 * (h_S(i,j+1,k) + h_N(i,j,k))
 !    h_marg = (2.0 * h_S(i,j+1,k) * h_N(i,j,k)) / &
 !             (h_S(i,j+1,k) + h_N(i,j,k) + GV%H_subroundoff)
    endif

    if (marginal) then ; h_v(i,J,k) = h_marg
    else ; h_v(i,J,k) = h_avg ; endif
  enddo ; enddo ; enddo

  if (present(visc_rem_v)) then
    ! Scale back the thickness to account for the effects of viscosity and the fractional open
    ! thickness to give an appropriate non-normalized weight for each layer in determining the
    ! barotropic acceleration.
    !$omp target teams distribute parallel do collapse(3)
    do k=1,nz ; do J=jsh-1,jeh ; do i=ish,ieh
      h_v(i,J,k) = h_v(i,J,k) * (visc_rem_v(i,J,k) * por_face_areaV(i,J,k))
    enddo ; enddo ; enddo
  else
    !$omp target teams distribute parallel do collapse(3)
    do k=1,nz ; do J=jsh-1,jeh ; do i=ish,ieh
      h_v(i,J,k) = h_v(i,J,k) * por_face_areaV(i,J,k)
    enddo ; enddo ; enddo
  endif

  local_open_BC = .false.
  if (associated(OBC)) local_open_BC = OBC%open_v_BCs_exist_globally
  if (local_open_BC) then
    ! This is done on the host, between copies of h_v from and back to the device.
    !$omp target update from(h_v)
    do n = 1, OBC%number_of_segments
      if (OBC%segment(n)%open .and. OBC%segment(n)%is_N_or_S) then
        J = OBC%segment(n)%HI%JsdB
        if (OBC%segment(n)%direction == OBC_DIRECTION_N) then
          if (present(visc_rem_v)) then ; do k=1,nz
            do i = OBC%segment(n)%HI%isd, OBC%segment(n)%HI%ied
              h_v(i,J,k) = h(i,j,k) * (visc_rem_v(i,J,k) * por_face_areaV(i,J,k))
            enddo
          enddo ; else ; do k=1,nz
            do i = OBC%segment(n)%HI%isd, OBC%segment(n)%HI%ied
              h_v(i,J,k) = h(i,j,k) * por_face_areaV(i,J,k)
            enddo
          enddo ; endif
        else
          if (present(visc_rem_v)) then ; do k=1,nz
            do i = OBC%segment(n)%HI%isd, OBC%segment(n)%HI%ied
              h_v(i,J,k) = h(i,j+1,k) * (visc_rem_v(i,J,k) * por_face_areaV(i,J,k))
            enddo
          enddo ; else ; do k=1,nz
            do i = OBC%segment(n)%HI%isd, OBC%segment(n)%HI%ied
              h_v(i,J,k) = h(i,j+1,k) * por_face_areaV(i,J,k)
            enddo
          enddo ; endif
        endif
      endif
    enddo
    !$omp target update to(h_v)
  endif

end procedure meridional_flux_thickness

!> Copies BT_cont and the components that meridional_mass_flux sets to the device.  As with
!! zonal_BT_cont_to_device, BT_cont is passed here as a plain argument to avoid an amdflang mapping
!! failure with its pointer argument in meridional_mass_flux.
subroutine merid_BT_cont_to_device(BT_cont, set_h_v)
  type(BT_cont_type), intent(inout) :: BT_cont !< A structure with elements that describe the
                                               !! effective open face areas as a function of barotropic flow.
  logical,            intent(in)    :: set_h_v !< If true, BT_cont%h_v is copied as well

  ! BT_cont is mapped apart from, and before, its components.
  !$omp target enter data map(to: BT_cont)
  !$omp target enter data map(to: BT_cont%FA_v_NN, BT_cont%FA_v_N0, BT_cont%FA_v_S0, &
  !$omp                           BT_cont%FA_v_SS, BT_cont%vBT_SS, BT_cont%vBT_NN)
  if (set_h_v) then
    !$omp target enter data map(to: BT_cont%h_v)
  endif
end subroutine merid_BT_cont_to_device

!> Copies the components of BT_cont that meridional_mass_flux sets back from the device, and
!! releases BT_cont there.
subroutine merid_BT_cont_from_device(BT_cont, set_h_v)
  type(BT_cont_type), intent(inout) :: BT_cont !< A structure with elements that describe the
                                               !! effective open face areas as a function of barotropic flow.
  logical,            intent(in)    :: set_h_v !< If true, BT_cont%h_v is copied as well

  if (set_h_v) then
    !$omp target exit data map(from: BT_cont%h_v)
  endif
  !$omp target exit data map(from: BT_cont%FA_v_NN, BT_cont%FA_v_N0, BT_cont%FA_v_S0, &
  !$omp                          BT_cont%FA_v_SS, BT_cont%vBT_SS, BT_cont%vBT_NN)
  ! Released after, and apart from, its components.
  !$omp target exit data map(release: BT_cont)
end subroutine merid_BT_cont_from_device

!> Returns the barotropic velocity adjustment that gives the desired barotropic (layer-summed) transport.
subroutine meridional_flux_adjust(v, h_in, h_S, h_N, vhbt, vh_tot_0, dvhdv_tot_0, &
                                  dv, dv_max_CFL, dv_min_CFL, dt, G, GV, US, CS, visc_rem, &
                                  i_start, i_end, j_start, j_end, do_I_in, por_face_areaV, vh_3d, &
                                  local_open_BC, open_dir)
  type(ocean_grid_type),   intent(in)    :: G    !< Ocean's grid structure.
  type(verticalGrid_type), intent(in)    :: GV   !< Ocean's vertical grid structure.
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                           intent(in)    :: v    !< Meridional velocity [L T-1 ~> m s-1].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in)    :: h_in !< Layer thickness used to calculate fluxes [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),&
                           intent(in)    :: h_S  !< South edge thickness in the reconstruction [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in)    :: h_N  !< North edge thickness in the reconstruction [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                           intent(in)    :: visc_rem
                             !< Both the fraction of the momentum originally
                             !! in a layer that remains after a time-step of viscosity, and the
                             !! fraction of a time-step's worth of a barotropic acceleration that
                             !! a layer experiences after viscosity is applied [nondim].
                             !! Visc_rem is between 0 (at the bottom) and 1 (far above the bottom).
  real, dimension(SZI_(G),SZJB_(G)), intent(in) :: vhbt !< The summed volume flux through meridional faces
                                                        !! [H L2 T-1 ~> m3 s-1 or kg s-1].
  real, dimension(SZI_(G),SZJB_(G)), intent(in) :: dv_max_CFL !< Maximum acceptable value of dv [L T-1 ~> m s-1].
  real, dimension(SZI_(G),SZJB_(G)), intent(in) :: dv_min_CFL !< Minimum acceptable value of dv [L T-1 ~> m s-1].
  real, dimension(SZI_(G),SZJB_(G)), intent(in) :: vh_tot_0   !< The summed transport with 0 adjustment
                                                        !! [H L2 T-1 ~> m3 s-1 or kg s-1].
  real, dimension(SZI_(G),SZJB_(G)), intent(in) :: dvhdv_tot_0 !< The partial derivative of dv_err with
                                                        !! dv at 0 adjustment [H L ~> m2 or kg m-1].
  real, dimension(SZI_(G),SZJB_(G)), intent(out) :: dv  !< The barotropic velocity adjustment [L T-1 ~> m s-1].
  real,                     intent(in)    :: dt   !< Time increment [T ~> s].
  type(unit_scale_type),    intent(in)    :: US   !< A dimensional unit scaling type
  type(continuity_PPM_CS),  intent(in)    :: CS   !< This module's control structure.
  integer,                  intent(in)    :: i_start !< Start of the i index range
  integer,                  intent(in)    :: i_end   !< End of the i index range
  integer,                  intent(in)    :: j_start !< Start of the J index range
  integer,                  intent(in)    :: j_end   !< End of the J index range
  logical, dimension(SZI_(G),SZJB_(G)), &
                            intent(in)    :: do_I_in  !< A flag indicating which points to work on.
  real, dimension(SZI_(G),SZJB_(G),SZK_(G)), &
                            intent(in)    :: por_face_areaV !< fractional open area of V-faces [nondim]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                            intent(inout) :: vh_3d !< Volume flux through meridional faces = v*h*dx
                                                   !! [H L2 T-1 ~> m3 s-1 or kg s-1], updated at the
                                                   !! points that are adjusted.
  logical,                  intent(in)    :: local_open_BC !< True if there are open boundary faces
                                                   !! on this PE whose fluxes are set from open_dir.
  integer, dimension(SZI_(G),SZJB_(G)), &
                            intent(in)    :: open_dir !< 1 or -1 at faces on open boundary segments
                                                   !! where the flow is taken from the cell to the south
                                                   !! or north of the face, or 0 elsewhere [nondim]
  ! Local variables
  real, dimension(SZI_(G),SZJB_(G)) :: &
    vh_err, &  ! Difference between vhbt and the summed vh [H L2 T-1 ~> m3 s-1 or kg s-1].
    vh_err_best, & ! The smallest value of vh_err found so far [H L2 T-1 ~> m3 s-1 or kg s-1].
    dvhdv_tot,&! Summed partial derivative of vh with u [H L ~> m2 or kg m-1].
    dv_min, &  ! Lower limit on dv correction based on CFL limits and previous iterations [L T-1 ~> m s-1]
    dv_max     ! Upper limit on dv correction based on CFL limits and previous iterations [L T-1 ~> m s-1]
  logical, dimension(SZI_(G),SZJB_(G)) :: &
    do_I       ! Indicates the points that are still being adjusted
  real :: v_new   ! The velocity with the correction added [L T-1 ~> m s-1].
  real :: dvhdv   ! Partial derivative of vh with v [H L ~> m2 or kg m-1].
  real :: dv_prev ! The previous value of dv [L T-1 ~> m s-1].
  real :: ddv     ! The change in dv from the previous iteration [L T-1 ~> m s-1].
  real :: tol_eta ! The tolerance for the current iteration [H ~> m or kg m-2].
  real :: tol_vel ! The tolerance for velocity in the current iteration [L T-1 ~> m s-1].
  real :: tol_eta_ref ! A copy of CS%tol_eta for use on the device [H ~> m or kg m-2].
  real :: h_marg_min  ! A copy of CS%h_marg_min for use on the device [H ~> m or kg m-2]
  logical :: vol_CFL, better_iter ! Copies of CS fields for use on the device
  integer :: i, j, k, nz, itt
#ifndef _OPENMP
  logical :: domore ! True if any point still needs to be adjusted
#endif
  integer, parameter :: max_itts = 20

  nz = GV%ke

  tol_vel = CS%tol_vel ; tol_eta_ref = CS%tol_eta ; better_iter = CS%better_iter
  vol_CFL = CS%vol_CFL ; h_marg_min = CS%h_marg_min

  !$omp target enter data map(alloc: vh_err, vh_err_best, dvhdv_tot, dv_min, dv_max, do_I)

  !$omp target teams private(tol_eta)
  !$omp distribute parallel do collapse(2)
  do J=j_start,j_end ; do i=i_start,i_end
    dv(i,J) = 0.0 ; do_I(i,J) = do_I_in(i,J)
    dv_max(i,J) = dv_max_CFL(i,J) ; dv_min(i,J) = dv_min_CFL(i,J)
    vh_err(i,J) = vh_tot_0(i,J) - vhbt(i,J) ; dvhdv_tot(i,J) = dvhdv_tot_0(i,J)
    vh_err_best(i,J) = abs(vh_err(i,J))
  enddo ; enddo

  do itt=1,max_itts
    select case (itt)
      case (:1) ; tol_eta = 1e-6 * tol_eta_ref
      case (2)  ; tol_eta = 1e-4 * tol_eta_ref
      case (3)  ; tol_eta = 1e-2 * tol_eta_ref
      case default ; tol_eta = tol_eta_ref
    end select

    !$omp distribute parallel do collapse(2)
    do J=j_start,j_end ; do i=i_start,i_end
      if (vh_err(i,J) > 0.0) then ; dv_max(i,J) = dv(i,J)
      elseif (vh_err(i,J) < 0.0) then ; dv_min(i,J) = dv(i,J)
      else ; do_I(i,J) = .false. ; endif
    enddo ; enddo
#ifndef _OPENMP
    domore = .false.
#endif
    !$omp distribute parallel do collapse(2) private(ddv, dv_prev)
    do J=j_start,j_end ; do i=i_start,i_end ; if (do_I(i,J)) then
      if ((dt * min(G%IareaT(i,j),G%IareaT(i,j+1))*abs(vh_err(i,J)) > tol_eta) .or. &
          (better_iter .and. ((abs(vh_err(i,J)) > tol_vel * dvhdv_tot(i,J)) .or. &
                                 (abs(vh_err(i,J)) > vh_err_best(i,J))) )) then
        !   Use Newton's method, provided it stays bounded.  Otherwise bisect
        ! the value with the appropriate bound.
        ddv = -vh_err(i,J) / dvhdv_tot(i,J)
        dv_prev = dv(i,J)
        dv(i,J) = dv(i,J) + ddv
        if (abs(ddv) < 1.0e-15*abs(dv(i,J))) then
          do_I(i,J) = .false. ! ddv is small enough to quit.
        elseif (ddv > 0.0) then
          if (dv(i,J) >= dv_max(i,J)) then
            dv(i,J) = 0.5*(dv_prev + dv_max(i,J))
            if (dv_max(i,J) - dv_prev < 1.0e-15*abs(dv(i,J))) do_I(i,J) = .false.
          endif
        else ! dvv(i,J) < 0.0
          if (dv(i,J) <= dv_min(i,J)) then
            dv(i,J) = 0.5*(dv_prev + dv_min(i,J))
            if (dv_prev - dv_min(i,J) < 1.0e-15*abs(dv(i,J))) do_I(i,J) = .false.
          endif
        endif
#ifndef _OPENMP
        if (do_I(i,J)) domore = .true.
#endif
      else
        do_I(i,J) = .false.
      endif
    endif ; enddo ; enddo
#ifndef _OPENMP
    ! Without OpenMP, stop as soon as every point has converged.  Iterations after that would not
    ! change any of the results, so with OpenMP, where the points are spread across teams that
    ! cannot share this flag, they are simply carried out.
    if (.not.domore) exit
#endif

    !$omp distribute parallel do collapse(2)
    do J=j_start,j_end ; do i=i_start,i_end
      vh_err(i,J) = -vhbt(i,J) ; dvhdv_tot(i,J) = 0.0
    enddo ; enddo
    do k=1,nz
      !$omp distribute parallel do collapse(2) private(v_new, dvhdv)
      do J=j_start,j_end ; do i=i_start,i_end ; if (do_I(i,J)) then
        v_new = v(i,J,k) + dv(i,J) * visc_rem(i,J,k)
        call flux_elem(v_new, h_in(i,j,k), h_in(i,j+1,k), h_S(i,j,k), h_S(i,j+1,k), h_N(i,j,k), &
                       h_N(i,j+1,k), vh_3d(i,J,k), dvhdv, visc_rem(i,J,k), G%dx_Cv(i,J), &
                       G%IareaT(i,j), G%IareaT(i,j+1), G%IdyT(i,j), G%IdyT(i,j+1), dt, vol_CFL, &
                       por_face_areaV(i,J,k), h_marg_min)
        if (local_open_BC) &
          call flux_elem_OBC(v_new, h_in(i,j,k), h_in(i,j+1,k), vh_3d(i,J,k), dvhdv, &
                             visc_rem(i,J,k), por_face_areaV(i,J,k), G%dx_Cv(i,J), h_marg_min, &
                             open_dir(i,J))
        vh_err(i,J) = vh_err(i,J) + vh_3d(i,J,k)
        dvhdv_tot(i,J) = dvhdv_tot(i,J) + dvhdv
      endif ; enddo ; enddo
    enddo
    !$omp distribute parallel do collapse(2)
    do J=j_start,j_end ; do i=i_start,i_end
      vh_err_best(i,J) = min(vh_err_best(i,J), abs(vh_err(i,J)))
    enddo ; enddo
  enddo ! itt-loop
  ! If there are any faces which have not converged to within the tolerance,
  ! so-be-it, or else use a final upwind correction?
  ! This never seems to happen with 20 iterations as max_itt.
  !$omp end target teams

  !$omp target exit data map(release: vh_err, vh_err_best, dvhdv_tot, dv_min, dv_max, do_I)

end subroutine meridional_flux_adjust

!> Sets of a structure that describes the meridional barotropic volume or mass fluxes as a
!! function of barotropic flow to agree closely with the sum of the layer's transports.
subroutine set_merid_BT_cont(v, h_in, h_S, h_N, BT_cont, vh_tot_0, dvhdv_tot_0, &
                             dv_max_CFL, dv_min_CFL, dt, G, GV, US, CS, visc_rem, &
                             visc_rem_max, i_start, i_end, j_start, j_end, do_I, por_face_areaV, &
                             open_dir)
  type(ocean_grid_type),                     intent(in)    :: G    !< Ocean's grid structure.
  type(verticalGrid_type),                   intent(in)    :: GV   !< Ocean's vertical grid structure.
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), intent(in)   :: v    !< Meridional velocity [L T-1 ~> m s-1].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in)    :: h_in !< Layer thickness used to calculate fluxes,
                                                                   !! [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in)    :: h_S  !< South edge thickness in the reconstruction,
                                                                   !! [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in)    :: h_N  !< North edge thickness in the reconstruction,
                                                                   !! [H ~> m or kg m-2].
  type(BT_cont_type),                        intent(inout) :: BT_cont !< A structure with elements
                       !! that describe the effective open face areas as a function of barotropic flow.
  real, dimension(SZI_(G),SZJB_(G)),         intent(in)    :: vh_tot_0    !< The summed transport
                       !! with 0 adjustment [H L2 T-1 ~> m3 s-1 or kg s-1].
  real, dimension(SZI_(G),SZJB_(G)),         intent(in)    :: dvhdv_tot_0 !< The partial derivative
                       !! of du_err with dv at 0 adjustment [H L ~> m2 or kg m-1].
  real, dimension(SZI_(G),SZJB_(G)),         intent(in)    :: dv_max_CFL !< Maximum acceptable value
                       !! of dv [L T-1 ~> m s-1].
  real, dimension(SZI_(G),SZJB_(G)),         intent(in)    :: dv_min_CFL !< Minimum acceptable value
                       !! of dv [L T-1 ~> m s-1].
  real,                                      intent(in)    :: dt   !< Time increment [T ~> s].
  type(unit_scale_type),                     intent(in)    :: US   !< A dimensional unit scaling type
  type(continuity_PPM_CS),                   intent(in)    :: CS   !< This module's control structure.
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), intent(in)   :: visc_rem !< Both the fraction of the
                       !! momentum originally in a layer that remains after a time-step
                       !! of viscosity, and the fraction of a time-step's worth of a barotropic
                       !! acceleration that a layer experiences after viscosity is applied [nondim].
                       !! Visc_rem is between 0 (at the bottom) and 1 (far above the bottom).
  real, dimension(SZI_(G),SZJB_(G)),         intent(in)    :: visc_rem_max !< Maximum allowable
                       !! visc_rem [nondim]
  integer,                                   intent(in)    :: i_start !< Start of the i index range
  integer,                                   intent(in)    :: i_end   !< End of the i index range
  integer,                                   intent(in)    :: j_start !< Start of the J index range
  integer,                                   intent(in)    :: j_end   !< End of the J index range
  logical, dimension(SZI_(G),SZJB_(G)),      intent(in)    :: do_I    !< A logical flag indicating
                       !! which points to work on.
  real, dimension(SZI_(G),SZJB_(G),SZK_(G)), intent(in)    :: por_face_areaV !< fractional open area
                       !! of V-faces [nondim]
  integer, dimension(SZI_(G),SZJB_(G)),      intent(in)    :: open_dir !< An array that is passed
                       !! on to meridional_flux_adjust, but is not used there by this call [nondim]
  ! Local variables
  real, dimension(SZI_(G),SZJB_(G)) :: &
    dv0, &        ! The barotropic velocity increment that gives 0 transport [L T-1 ~> m s-1].
    dvL, dvR, &   ! The barotropic velocity increments that give the southerly
                  ! (dvL) and northerly (dvR) test velocities [L T-1 ~> m s-1].
    zeros, &      ! An array of full of 0 transports [H L2 T-1 ~> m3 s-1 or kg s-1]
    dv_CFL, &     ! The velocity increment that corresponds to CFL_min [L T-1 ~> m s-1].
    FAmt_L, FAmt_R, & ! The summed effective marginal face areas for the 3
    FAmt_0, &     ! test velocities [H L ~> m2 or kg m-1].
    vhtot_L, &    ! The summed transport with the southerly (vhtot_L) and
    vhtot_R       ! and northerly (vhtot_R) test velocities [H L2 T-1 ~> m3 s-1 or kg s-1].
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)) :: &
    vh_tmp        ! The layer transports found while diagnosing dv0, which are not used
                  ! [H L2 T-1 ~> m3 s-1 or kg s-1].
  real :: &
    v_L, v_R, &   ! The southerly (v_L), northerly (v_R), and zero-barotropic
    v_0, &        ! transport (v_0) layer test velocities [L T-1 ~> m s-1].
    dvhdv_L, &    ! The effective layer marginal face areas with the southerly
    dvhdv_R, &    ! (_L), northerly (_R), and zero-barotropic (_0) test
    dvhdv_0, &    ! velocities [H L ~> m2 or kg m-1].
    vh_L, vh_R, & ! The layer transports with the southerly (_L), northerly (_R)
    vh_0          ! and zero-barotropic (_0) test velocities [H L2 T-1 ~> m3 s-1 or kg s-1].
  real :: FA_0    ! The effective face area with 0 barotropic transport [H L ~> m2 or kg m-1].
  real :: FA_avg  ! The average effective face area [H L ~> m2 or kg m-1], nominally given by
                  ! the realized transport divided by the barotropic velocity.
  real :: visc_rem_lim ! The larger of visc_rem and min_visc_rem [nondim]  This
                       ! limiting is necessary to keep the inverse of visc_rem
                       ! from leading to large CFL numbers.
  real :: min_visc_rem ! The smallest permitted value for visc_rem that is used
                       ! in finding the barotropic velocity that changes the
                       ! flow direction [nondim].  This is necessary to keep the inverse
                       ! of visc_rem from leading to large CFL numbers.
  real :: CFL_min ! A minimal increment in the CFL to try to ensure that the
                  ! flow is truly upwind [nondim]
  real :: Idt     ! The inverse of the time step [T-1 ~> s-1].
  real :: h_marg_min ! A copy of CS%h_marg_min for use on the device [H ~> m or kg m-2]
  logical :: vol_CFL ! A copy of CS%vol_CFL for use on the device
  integer :: i, j, k, nz

  nz = GV%ke ; Idt = 1.0 / dt
  min_visc_rem = 0.1 ; CFL_min = 1e-6
  vol_CFL = CS%vol_CFL ; h_marg_min = CS%h_marg_min

  !$omp target enter data map(alloc: dv0, zeros, dvL, dvR, dv_CFL, FAmt_L, FAmt_R, FAmt_0, &
  !$omp                              vhtot_L, vhtot_R, vh_tmp)

  ! Diagnose the zero-transport correction, dv0.
  !$omp target teams distribute parallel do collapse(2)
  do J=j_start,j_end ; do i=i_start,i_end
    zeros(i,J) = 0.0
  enddo ; enddo
  call meridional_flux_adjust(v, h_in, h_S, h_N, zeros, vh_tot_0, dvhdv_tot_0, dv0, &
                              dv_max_CFL, dv_min_CFL, dt, G, GV, US, CS, visc_rem, &
                              i_start, i_end, j_start, j_end, do_I, por_face_areaV, vh_tmp, &
                              .false., open_dir)

  !   Determine the southerly- and northerly- fluxes.  Choose a sufficiently
  ! negative velocity correction for the northerly-flux, and a sufficiently
  ! positive correction for the southerly-flux.
  !$omp target teams
  !$omp distribute parallel do collapse(2)
  do J=j_start,j_end ; do i=i_start,i_end ; if (do_I(i,J)) then
    dv_CFL(i,J) = (CFL_min * Idt) * G%dyCv(i,J)
    dvR(i,J) = min(0.0,dv0(i,J) - dv_CFL(i,J))
    dvL(i,J) = max(0.0,dv0(i,J) + dv_CFL(i,J))
    FAmt_L(i,J) = 0.0 ; FAmt_R(i,J) = 0.0 ; FAmt_0(i,J) = 0.0
    vhtot_L(i,J) = 0.0 ; vhtot_R(i,J) = 0.0
  endif ; enddo ; enddo

  do k=1,nz
    !$omp distribute parallel do collapse(2) private(visc_rem_lim)
    do J=j_start,j_end ; do i=i_start,i_end ; if (do_I(i,J)) then
      visc_rem_lim = max(visc_rem(i,J,k), min_visc_rem*visc_rem_max(i,J))
      if (visc_rem_lim > 0.0) then ! This is almost always true for ocean points.
        if (v(i,J,k) + dvR(i,J)*visc_rem_lim > -dv_CFL(i,J)*visc_rem(i,J,k)) &
          dvR(i,J) = -(v(i,J,k) + dv_CFL(i,J)*visc_rem(i,J,k)) / visc_rem_lim
        if (v(i,J,k) + dvL(i,J)*visc_rem_lim < dv_CFL(i,J)*visc_rem(i,J,k)) &
          dvL(i,J) = -(v(i,J,k) - dv_CFL(i,J)*visc_rem(i,J,k)) / visc_rem_lim
      endif
    endif ; enddo ; enddo
  enddo

  do k=1,nz
    !$omp distribute parallel do collapse(2) &
    !$omp   private(v_L, v_R, v_0, dvhdv_0, dvhdv_L, dvhdv_R, vh_0, vh_L, vh_R)
    do J=j_start,j_end ; do i=i_start,i_end ; if (do_I(i,J)) then
      v_L = v(i,J,k) + dvL(i,J) * visc_rem(i,J,k)
      v_R = v(i,J,k) + dvR(i,J) * visc_rem(i,J,k)
      v_0 = v(i,J,k) + dv0(i,J) * visc_rem(i,J,k)
      call flux_elem(v_0, h_in(i,j,k), h_in(i,j+1,k), h_S(i,j,k), h_S(i,j+1,k), h_N(i,j,k), &
                     h_N(i,j+1,k), vh_0, dvhdv_0, visc_rem(i,J,k), G%dx_Cv(i,J), &
                     G%IareaT(i,j), G%IareaT(i,j+1), G%IdyT(i,j), G%IdyT(i,j+1), dt, vol_CFL, &
                     por_face_areaV(i,J,k), h_marg_min)
      call flux_elem(v_L, h_in(i,j,k), h_in(i,j+1,k), h_S(i,j,k), h_S(i,j+1,k), h_N(i,j,k), &
                     h_N(i,j+1,k), vh_L, dvhdv_L, visc_rem(i,J,k), G%dx_Cv(i,J), &
                     G%IareaT(i,j), G%IareaT(i,j+1), G%IdyT(i,j), G%IdyT(i,j+1), dt, vol_CFL, &
                     por_face_areaV(i,J,k), h_marg_min)
      call flux_elem(v_R, h_in(i,j,k), h_in(i,j+1,k), h_S(i,j,k), h_S(i,j+1,k), h_N(i,j,k), &
                     h_N(i,j+1,k), vh_R, dvhdv_R, visc_rem(i,J,k), G%dx_Cv(i,J), &
                     G%IareaT(i,j), G%IareaT(i,j+1), G%IdyT(i,j), G%IdyT(i,j+1), dt, vol_CFL, &
                     por_face_areaV(i,J,k), h_marg_min)
      FAmt_0(i,J) = FAmt_0(i,J) + dvhdv_0
      FAmt_L(i,J) = FAmt_L(i,J) + dvhdv_L
      FAmt_R(i,J) = FAmt_R(i,J) + dvhdv_R
      vhtot_L(i,J) = vhtot_L(i,J) + vh_L
      vhtot_R(i,J) = vhtot_R(i,J) + vh_R
    endif ; enddo ; enddo
  enddo

  !$omp distribute parallel do collapse(2) private(FA_0, FA_avg)
  do J=j_start,j_end ; do i=i_start,i_end
    if (do_I(i,J)) then
      FA_0 = FAmt_0(i,J) ; FA_avg = FAmt_0(i,J)
      if ((dvL(i,J) - dv0(i,J)) /= 0.0) &
        FA_avg = vhtot_L(i,J) / (dvL(i,J) - dv0(i,J))
      if (FA_avg > max(FA_0, FAmt_L(i,J))) then ; FA_avg = max(FA_0, FAmt_L(i,J))
      elseif (FA_avg < min(FA_0, FAmt_L(i,J))) then ; FA_0 = FA_avg ; endif
      BT_cont%FA_v_S0(i,J) = FA_0 ; BT_cont%FA_v_SS(i,J) = FAmt_L(i,J)
      if (abs(FA_0-FAmt_L(i,J)) <= 1e-12*FA_0) then ; BT_cont%vBT_SS(i,J) = 0.0 ; else
        BT_cont%vBT_SS(i,J) = (1.5 * (dvL(i,J) - dv0(i,J))) * &
                     ((FAmt_L(i,J) - FA_avg) / (FAmt_L(i,J) - FA_0))
      endif

      FA_0 = FAmt_0(i,J) ; FA_avg = FAmt_0(i,J)
      if ((dvR(i,J) - dv0(i,J)) /= 0.0) &
        FA_avg = vhtot_R(i,J) / (dvR(i,J) - dv0(i,J))
      if (FA_avg > max(FA_0, FAmt_R(i,J))) then ; FA_avg = max(FA_0, FAmt_R(i,J))
      elseif (FA_avg < min(FA_0, FAmt_R(i,J))) then ; FA_0 = FA_avg ; endif
      BT_cont%FA_v_N0(i,J) = FA_0 ; BT_cont%FA_v_NN(i,J) = FAmt_R(i,J)
      if (abs(FAmt_R(i,J) - FA_0) <= 1e-12*FA_0) then ; BT_cont%vBT_NN(i,J) = 0.0 ; else
        BT_cont%vBT_NN(i,J) = (1.5 * (dvR(i,J) - dv0(i,J))) * &
                     ((FAmt_R(i,J) - FA_avg) / (FAmt_R(i,J) - FA_0))
      endif
    else
      BT_cont%FA_v_S0(i,J) = 0.0 ; BT_cont%FA_v_SS(i,J) = 0.0
      BT_cont%FA_v_N0(i,J) = 0.0 ; BT_cont%FA_v_NN(i,J) = 0.0
      BT_cont%vBT_SS(i,J) = 0.0 ; BT_cont%vBT_NN(i,J) = 0.0
    endif
  enddo ; enddo
  !$omp end target teams

  !$omp target exit data map(release: dv0, zeros, dvL, dvR, dv_CFL, FAmt_L, FAmt_R, FAmt_0, &
  !$omp                               vhtot_L, vhtot_R, vh_tmp)

end subroutine set_merid_BT_cont

module procedure ratio_max
  !$omp declare target
  if (abs(a) > abs(maxrat*b)) then
    ratio = maxrat
  else
    ratio = a / b
  endif
end procedure ratio_max

end submodule MOM_continuity_PPM_s
