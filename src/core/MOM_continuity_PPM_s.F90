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
  real :: FA_u    ! A sum of zonal face areas [H L ~> m2 or kg m-1].
  real :: I_vrm   ! 1.0 / visc_rem_max [nondim]
  real :: CFL_dt  ! The maximum CFL ratio of the adjusted velocities divided by
                  ! the time step [T-1 ~> s-1].
  real :: I_dt    ! 1.0 / dt [T-1 ~> s-1].
  real :: du_lim  ! The velocity change that give a relative CFL of 1 [L T-1 ~> m s-1].
  real :: dx_E, dx_W ! Effective x-grid spacings to the east and west [L ~> m].
  type(cont_loop_bounds_type) :: LB
  integer :: i, j, k, ish, ieh, jsh, jeh, n, nz
  integer :: l_seg ! The OBC segment number
  logical :: use_visc_rem, set_BT_cont
  logical :: local_specified_BC, local_Flather_OBC, local_open_BC, any_simple_OBC  ! OBC-related logicals

  call cpu_clock_begin(id_clock_correct)

  use_visc_rem = present(visc_rem_u)

  set_BT_cont = .false. ; if (present(BT_cont)) set_BT_cont = (associated(BT_cont))

  local_specified_BC = .false. ; local_Flather_OBC = .false. ; local_open_BC = .false.
  if (associated(OBC)) then ; if (OBC%OBC_pe) then
    local_specified_BC = OBC%specified_u_BCs_exist_globally
    local_Flather_OBC = OBC%Flather_u_BCs_exist_globally
    local_open_BC = OBC%open_u_BCs_exist_globally
  endif ; endif

  if (present(du_cor)) then
    do j=G%jsd,G%jed ; do I=G%IsdB,G%IedB
      du_cor(I,j) = 0.0
    enddo ; enddo
  endif

  if (present(LB_in)) then
    LB = LB_in
  else
    LB%ish = G%isc ; LB%ieh = G%iec ; LB%jsh = G%jsc ; LB%jeh = G%jec
  endif
  ish = LB%ish ; ieh = LB%ieh ; jsh = LB%jsh ; jeh = LB%jeh ; nz = GV%ke

  CFL_dt = CS%CFL_limit_adjust / dt
  I_dt = 1.0 / dt
  if (CS%aggress_adjust) CFL_dt = I_dt

  ! Set uh and duhdu.
  do k=1,nz
    if (use_visc_rem) then
      do j=jsh,jeh ; do I=ish-1,ieh
        visc_rem(I,j,k) = visc_rem_u(I,j,k)
      enddo ; enddo
    else
      do j=jsh,jeh ; do I=ish-1,ieh
        visc_rem(I,j,k) = 1.0
      enddo ; enddo
    endif
    do j=jsh,jeh ; do I=ish-1,ieh
      call flux_elem(u(I,j,k), h_in(i,j,k), h_in(i+1,j,k), h_W(i,j,k), h_W(i+1,j,k), h_E(i,j,k), &
                     h_E(i+1,j,k), uh(I,j,k), duhdu(I,j,k), visc_rem(I,j,k), G%dy_Cu(I,j), &
                     G%IareaT(i,j), G%IareaT(i+1,j), G%IdxT(i,j), G%IdxT(i+1,j), dt, CS%vol_CFL, &
                     por_face_areaU(I,j,k), CS%h_marg_min)
      if (local_open_BC) &
        call flux_elem_OBC(u(I,j,k), h_in(i,j,k), h_in(i+1,j,k), uh(I,j,k), duhdu(I,j,k), &
                           visc_rem(I,j,k), por_face_areaU(I,j,k), G%dy_Cu(I,j), CS%h_marg_min, &
                           OBC, OBC%segnum_u(I,j))
    enddo ; enddo
    if (local_specified_BC) then
      do j=jsh,jeh ; do I=ish-1,ieh ; if (OBC%segnum_u(I,j) /= 0) then
        l_seg = abs(OBC%segnum_u(I,j))
        if (OBC%segment(l_seg)%specified) uh(I,j,k) = OBC%segment(l_seg)%normal_trans(I,j,k)
      endif ; enddo ; enddo
    endif
  enddo

  if (present(uhbt) .or. set_BT_cont) then
    if (use_visc_rem .and. CS%use_visc_rem_max) then
      do j=jsh,jeh ; do I=ish-1,ieh
        visc_rem_max(I,j) = 0.0
      enddo ; enddo
      do k=1,nz
        do j=jsh,jeh ; do I=ish-1,ieh
          visc_rem_max(I,j) = max(visc_rem_max(I,j), visc_rem(I,j,k))
        enddo ; enddo
      enddo
    else
      do j=jsh,jeh ; do I=ish-1,ieh
        visc_rem_max(I,j) = 1.0
      enddo ; enddo
    endif
    !   Set limits on du that will keep the CFL number between -1 and 1.
    ! This should be adequate to keep the root bracketed in all cases.
    do j=jsh,jeh ; do I=ish-1,ieh
      I_vrm = 0.0
      if (visc_rem_max(I,j) > 0.0) I_vrm = 1.0 / visc_rem_max(I,j)
      if (CS%vol_CFL) then
        dx_W = ratio_max(G%areaT(i,j), G%dy_Cu(I,j), 1000.0*G%dxT(i,j))
        dx_E = ratio_max(G%areaT(i+1,j), G%dy_Cu(I,j), 1000.0*G%dxT(i+1,j))
      else ; dx_W = G%dxT(i,j) ; dx_E = G%dxT(i+1,j) ; endif
      du_max_CFL(I,j) = 2.0* (CFL_dt * dx_W) * I_vrm
      du_min_CFL(I,j) = -2.0 * (CFL_dt * dx_E) * I_vrm
      uh_tot_0(I,j) = 0.0 ; duhdu_tot_0(I,j) = 0.0
    enddo ; enddo
    do k=1,nz
      do j=jsh,jeh ; do I=ish-1,ieh
        duhdu_tot_0(I,j) = duhdu_tot_0(I,j) + duhdu(I,j,k)
        uh_tot_0(I,j) = uh_tot_0(I,j) + uh(I,j,k)
      enddo ; enddo
    enddo
    if (use_visc_rem) then
      if (CS%aggress_adjust) then
        do k=1,nz
          do j=jsh,jeh ; do I=ish-1,ieh
            if (CS%vol_CFL) then
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
        enddo
      else
        do k=1,nz
          do j=jsh,jeh ; do I=ish-1,ieh
            if (CS%vol_CFL) then
              dx_W = ratio_max(G%areaT(i,j), G%dy_Cu(I,j), 1000.0*G%dxT(i,j))
              dx_E = ratio_max(G%areaT(i+1,j), G%dy_Cu(I,j), 1000.0*G%dxT(i+1,j))
            else ; dx_W = G%dxT(i,j) ; dx_E = G%dxT(i+1,j) ; endif

            if (du_max_CFL(I,j) * visc_rem(I,j,k) > dx_W*CFL_dt - u(I,j,k)*G%mask2dCu(I,j)) &
              du_max_CFL(I,j) = (dx_W*CFL_dt - u(I,j,k)) / visc_rem(I,j,k)
            if (du_min_CFL(I,j) * visc_rem(I,j,k) < -dx_E*CFL_dt - u(I,j,k)*G%mask2dCu(I,j)) &
              du_min_CFL(I,j) = -(dx_E*CFL_dt + u(I,j,k)) / visc_rem(I,j,k)
          enddo ; enddo
        enddo
      endif
    else
      if (CS%aggress_adjust) then
        do k=1,nz
          do j=jsh,jeh ; do I=ish-1,ieh
            if (CS%vol_CFL) then
              dx_W = ratio_max(G%areaT(i,j), G%dy_Cu(I,j), 1000.0*G%dxT(i,j))
              dx_E = ratio_max(G%areaT(i+1,j), G%dy_Cu(I,j), 1000.0*G%dxT(i+1,j))
            else ; dx_W = G%dxT(i,j) ; dx_E = G%dxT(i+1,j) ; endif

            du_max_CFL(I,j) = MIN(du_max_CFL(I,j), 0.499 * &
                        ((dx_W*I_dt - u(I,j,k)) + MIN(0.0,u(I-1,j,k))) )
            du_min_CFL(I,j) = MAX(du_min_CFL(I,j), 0.499 * &
                        ((-dx_E*I_dt - u(I,j,k)) + MAX(0.0,u(I+1,j,k))) )
          enddo ; enddo
        enddo
      else
        do k=1,nz
          do j=jsh,jeh ; do I=ish-1,ieh
            if (CS%vol_CFL) then
              dx_W = ratio_max(G%areaT(i,j), G%dy_Cu(I,j), 1000.0*G%dxT(i,j))
              dx_E = ratio_max(G%areaT(i+1,j), G%dy_Cu(I,j), 1000.0*G%dxT(i+1,j))
            else ; dx_W = G%dxT(i,j) ; dx_E = G%dxT(i+1,j) ; endif

            du_max_CFL(I,j) = MIN(du_max_CFL(I,j), dx_W*CFL_dt - u(I,j,k))
            du_min_CFL(I,j) = MAX(du_min_CFL(I,j), -(dx_E*CFL_dt + u(I,j,k)))
          enddo ; enddo
        enddo
      endif
    endif
    do j=jsh,jeh ; do I=ish-1,ieh
      du_max_CFL(I,j) = max(du_max_CFL(I,j),0.0)
      du_min_CFL(I,j) = min(du_min_CFL(I,j),0.0)
    enddo ; enddo

    any_simple_OBC = .false.
    if (local_specified_BC .or. local_Flather_OBC) then
      do j=jsh,jeh ; do I=ish-1,ieh
        l_seg = abs(OBC%segnum_u(I,j))

        ! Avoid reconciling barotropic/baroclinic transports if transport is specified
        simple_OBC_pt(I,j) = .false.
        if (l_seg /= OBC_NONE) simple_OBC_pt(I,j) = OBC%segment(l_seg)%specified
        do_I(I,j) = .not.simple_OBC_pt(I,j)
        any_simple_OBC = any_simple_OBC .or. simple_OBC_pt(I,j)
      enddo ; enddo
    else
      do j=jsh,jeh ; do I=ish-1,ieh
        do_I(I,j) = .true.
      enddo ; enddo
    endif

    if (present(uhbt)) then
      ! Find du and uh.
      call zonal_flux_adjust(u, h_in, h_W, h_E, uhbt, uh_tot_0, duhdu_tot_0, du, &
                             du_max_CFL, du_min_CFL, dt, G, GV, US, CS, visc_rem, &
                             ish-1, ieh, jsh, jeh, do_I, por_face_areaU, uh, OBC=OBC)

      if (present(u_cor)) then
        do k=1,nz
          do j=jsh,jeh ; do I=ish-1,ieh
            u_cor(I,j,k) = u(I,j,k) + du(I,j) * visc_rem(I,j,k)
          enddo ; enddo
          if (any_simple_OBC) then
            do j=jsh,jeh ; do I=ish-1,ieh ; if (simple_OBC_pt(I,j)) then
              u_cor(I,j,k) = OBC%segment(abs(OBC%segnum_u(I,j)))%normal_vel(I,j,k)
            endif ; enddo ; enddo
          endif
        enddo
      endif ! u-corrected

      if (present(du_cor)) then
        do j=jsh,jeh ; do I=ish-1,ieh
          du_cor(I,j) = du(I,j)
        enddo ; enddo
      endif

    endif

    if (set_BT_cont) then
      call set_zonal_BT_cont(u, h_in, h_W, h_E, BT_cont, uh_tot_0, duhdu_tot_0, &
                             du_max_CFL, du_min_CFL, dt, G, GV, US, CS, visc_rem, &
                             visc_rem_max, ish-1, ieh, jsh, jeh, do_I, por_face_areaU)
      if (any_simple_OBC) then
        do j=jsh,jeh ; do I=ish-1,ieh
          if (simple_OBC_pt(I,j)) FAuI(I,j) = GV%H_subroundoff*G%dy_Cu(I,j)
        enddo ; enddo
        ! NOTE: simple_OBC_pt should prevent access to segment OBC_NONE
        do k=1,nz
          do j=jsh,jeh ; do I=ish-1,ieh ; if (simple_OBC_pt(I,j)) then
            l_seg = abs(OBC%segnum_u(I,j))
            if ((abs(OBC%segment(l_seg)%normal_vel(I,j,k)) > 0.0) .and. (OBC%segment(l_seg)%specified)) &
              FAuI(I,j) = FAuI(I,j) + OBC%segment(l_seg)%normal_trans(I,j,k) / OBC%segment(l_seg)%normal_vel(I,j,k)
          endif ; enddo ; enddo
        enddo
        do j=jsh,jeh ; do I=ish-1,ieh ; if (simple_OBC_pt(I,j)) then
          BT_cont%FA_u_W0(I,j) = FAuI(I,j) ; BT_cont%FA_u_E0(I,j) = FAuI(I,j)
          BT_cont%FA_u_WW(I,j) = FAuI(I,j) ; BT_cont%FA_u_EE(I,j) = FAuI(I,j)
          BT_cont%uBT_WW(I,j) = 0.0 ; BT_cont%uBT_EE(I,j) = 0.0
        endif ; enddo ; enddo
      endif
    endif ! set_BT_cont

  endif ! present(uhbt) or set_BT_cont

  if (local_open_BC .and. set_BT_cont) then
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
  endif

  if  (set_BT_cont) then ; if (allocated(BT_cont%h_u)) then
    if (present(u_cor)) then
      call zonal_flux_thickness(u_cor, h_in, h_W, h_E, BT_cont%h_u, dt, G, GV, US, LB, &
                                CS%vol_CFL, CS%marginal_faces, OBC, por_face_areaU, visc_rem_u)
    else
      call zonal_flux_thickness(u, h_in, h_W, h_E, BT_cont%h_u, dt, G, GV, US, LB, &
                                CS%vol_CFL, CS%marginal_faces, OBC, por_face_areaU, visc_rem_u)
    endif
  endif ; endif

  call cpu_clock_end(id_clock_correct)

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

  !$OMP parallel do default(shared) private(CFL,curv_3,h_marg,h_avg)
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
    !$OMP parallel do default(shared)
    do k=1,nz ; do j=jsh,jeh ; do I=ish-1,ieh
      h_u(I,j,k) = h_u(I,j,k) * (visc_rem_u(I,j,k) * por_face_areaU(I,j,k))
    enddo ; enddo ; enddo
  else
    !$OMP parallel do default(shared)
    do k=1,nz ; do j=jsh,jeh ; do I=ish-1,ieh
      h_u(I,j,k) = h_u(I,j,k) * por_face_areaU(I,j,k)
    enddo ; enddo ; enddo
  endif

  local_open_BC = .false.
  if (associated(OBC)) local_open_BC = OBC%open_u_BCs_exist_globally
  if (local_open_BC) then
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
  endif

end procedure zonal_flux_thickness

!> Returns the barotropic velocity adjustment that gives the
!! desired barotropic (layer-summed) transport.
subroutine zonal_flux_adjust(u, h_in, h_W, h_E, uhbt, uh_tot_0, duhdu_tot_0, &
                             du, du_max_CFL, du_min_CFL, dt, G, GV, US, CS, visc_rem, &
                             i_start, i_end, j_start, j_end, do_I_in, por_face_areaU, uh_3d, OBC)
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
  type(ocean_OBC_type),            optional, pointer       :: OBC !< Open boundaries control structure.
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
  integer :: i, j, k, nz, itt
  logical :: local_open_BC ! True if there are open OBC points on this PE
#ifndef _OPENMP
  logical :: domore ! True if any point still needs to be adjusted
#endif
  integer, parameter :: max_itts = 20

  nz = GV%ke
  local_open_BC = .false.
  if (present(OBC)) then ; if (associated(OBC)) local_open_BC = OBC%open_u_BCs_exist_globally ; endif

  tol_vel = CS%tol_vel

  do j=j_start,j_end ; do I=i_start,i_end
    du(I,j) = 0.0 ; do_I(I,j) = do_I_in(I,j)
    du_max(I,j) = du_max_CFL(I,j) ; du_min(I,j) = du_min_CFL(I,j)
    uh_err(I,j) = uh_tot_0(I,j) - uhbt(I,j) ; duhdu_tot(I,j) = duhdu_tot_0(I,j)
    uh_err_best(I,j) = abs(uh_err(I,j))
  enddo ; enddo

  do itt=1,max_itts
    select case (itt)
      case (:1) ; tol_eta = 1e-6 * CS%tol_eta
      case (2)  ; tol_eta = 1e-4 * CS%tol_eta
      case (3)  ; tol_eta = 1e-2 * CS%tol_eta
      case default ; tol_eta = CS%tol_eta
    end select

    do j=j_start,j_end ; do I=i_start,i_end
      if (uh_err(I,j) > 0.0) then ; du_max(I,j) = du(I,j)
      elseif (uh_err(I,j) < 0.0) then ; du_min(I,j) = du(I,j)
      else ; do_I(I,j) = .false. ; endif
    enddo ; enddo
#ifndef _OPENMP
    domore = .false.
#endif
    do j=j_start,j_end ; do I=i_start,i_end ; if (do_I(I,j)) then
      if ((dt * min(G%IareaT(i,j),G%IareaT(i+1,j))*abs(uh_err(I,j)) > tol_eta) .or. &
          (CS%better_iter .and. ((abs(uh_err(I,j)) > tol_vel * duhdu_tot(I,j)) .or. &
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

    do j=j_start,j_end ; do I=i_start,i_end
      uh_err(I,j) = -uhbt(I,j) ; duhdu_tot(I,j) = 0.0
    enddo ; enddo
    do k=1,nz
      do j=j_start,j_end ; do I=i_start,i_end ; if (do_I(I,j)) then
        u_new = u(I,j,k) + du(I,j) * visc_rem(I,j,k)
        call flux_elem(u_new, h_in(i,j,k), h_in(i+1,j,k), h_W(i,j,k), h_W(i+1,j,k), h_E(i,j,k), &
                       h_E(i+1,j,k), uh_3d(I,j,k), duhdu, visc_rem(I,j,k), G%dy_Cu(I,j), &
                       G%IareaT(i,j), G%IareaT(i+1,j), G%IdxT(i,j), G%IdxT(i+1,j), dt, CS%vol_CFL, &
                       por_face_areaU(I,j,k), CS%h_marg_min)
        if (local_open_BC) &
          call flux_elem_OBC(u_new, h_in(i,j,k), h_in(i+1,j,k), uh_3d(I,j,k), duhdu, &
                             visc_rem(I,j,k), por_face_areaU(I,j,k), G%dy_Cu(I,j), CS%h_marg_min, &
                             OBC, OBC%segnum_u(I,j))
        uh_err(I,j) = uh_err(I,j) + uh_3d(I,j,k)
        duhdu_tot(I,j) = duhdu_tot(I,j) + duhdu
      endif ; enddo ; enddo
    enddo
    do j=j_start,j_end ; do I=i_start,i_end
      uh_err_best(I,j) = min(uh_err_best(I,j), abs(uh_err(I,j)))
    enddo ; enddo
  enddo ! itt-loop
  ! If there are any faces which have not converged to within the tolerance,
  ! so-be-it, or else use a final upwind correction?
  ! This never seems to happen with 20 iterations as max_itt.

end subroutine zonal_flux_adjust

!> Sets a structure that describes the zonal barotropic volume or mass fluxes as a
!! function of barotropic flow to agree closely with the sum of the layer's transports.
subroutine set_zonal_BT_cont(u, h_in, h_W, h_E, BT_cont, uh_tot_0, duhdu_tot_0, &
                             du_max_CFL, du_min_CFL, dt, G, GV, US, CS, visc_rem, &
                             visc_rem_max, i_start, i_end, j_start, j_end, do_I, por_face_areaU)
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
  integer :: i, j, k, nz

  nz = GV%ke ; Idt = 1.0 / dt
  min_visc_rem = 0.1 ; CFL_min = 1e-6

  ! Diagnose the zero-transport correction, du0.
  do j=j_start,j_end ; do I=i_start,i_end
    zeros(I,j) = 0.0
  enddo ; enddo
  call zonal_flux_adjust(u, h_in, h_W, h_E, zeros, uh_tot_0, duhdu_tot_0, du0, &
                         du_max_CFL, du_min_CFL, dt, G, GV, US, CS, visc_rem, &
                         i_start, i_end, j_start, j_end, do_I, por_face_areaU, uh_tmp)

  ! Determine the westerly- and easterly- fluxes.  Choose a sufficiently
  ! negative velocity correction for the easterly-flux, and a sufficiently
  ! positive correction for the westerly-flux.
  do j=j_start,j_end ; do I=i_start,i_end
    du_CFL(I,j) = (CFL_min * Idt) * G%dxCu(I,j)
    duR(I,j) = min(0.0,du0(I,j) - du_CFL(I,j))
    duL(I,j) = max(0.0,du0(I,j) + du_CFL(I,j))
    FAmt_L(I,j) = 0.0 ; FAmt_R(I,j) = 0.0 ; FAmt_0(I,j) = 0.0
    uhtot_L(I,j) = 0.0 ; uhtot_R(I,j) = 0.0
  enddo ; enddo

  do k=1,nz
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
    do j=j_start,j_end ; do I=i_start,i_end ; if (do_I(I,j)) then
      u_L = u(I,j,k) + duL(I,j) * visc_rem(I,j,k)
      u_R = u(I,j,k) + duR(I,j) * visc_rem(I,j,k)
      u_0 = u(I,j,k) + du0(I,j) * visc_rem(I,j,k)
      call flux_elem(u_0, h_in(i,j,k), h_in(i+1,j,k), h_W(i,j,k), h_W(i+1,j,k), h_E(i,j,k), &
                     h_E(i+1,j,k), uh_0, duhdu_0, visc_rem(I,j,k), G%dy_Cu(I,j), &
                     G%IareaT(i,j), G%IareaT(i+1,j), G%IdxT(i,j), G%IdxT(i+1,j), dt, CS%vol_CFL, &
                     por_face_areaU(I,j,k), CS%h_marg_min)
      call flux_elem(u_L, h_in(i,j,k), h_in(i+1,j,k), h_W(i,j,k), h_W(i+1,j,k), h_E(i,j,k), &
                     h_E(i+1,j,k), uh_L, duhdu_L, visc_rem(I,j,k), G%dy_Cu(I,j), &
                     G%IareaT(i,j), G%IareaT(i+1,j), G%IdxT(i,j), G%IdxT(i+1,j), dt, CS%vol_CFL, &
                     por_face_areaU(I,j,k), CS%h_marg_min)
      call flux_elem(u_R, h_in(i,j,k), h_in(i+1,j,k), h_W(i,j,k), h_W(i+1,j,k), h_E(i,j,k), &
                     h_E(i+1,j,k), uh_R, duhdu_R, visc_rem(I,j,k), G%dy_Cu(I,j), &
                     G%IareaT(i,j), G%IareaT(i+1,j), G%IdxT(i,j), G%IdxT(i+1,j), dt, CS%vol_CFL, &
                     por_face_areaU(I,j,k), CS%h_marg_min)
      FAmt_0(I,j) = FAmt_0(I,j) + duhdu_0
      FAmt_L(I,j) = FAmt_L(I,j) + duhdu_L
      FAmt_R(I,j) = FAmt_R(I,j) + duhdu_R
      uhtot_L(I,j) = uhtot_L(I,j) + uh_L
      uhtot_R(I,j) = uhtot_R(I,j) + uh_R
    endif ; enddo ; enddo
  enddo

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

end subroutine set_zonal_BT_cont

!> Calculates the mass or volume flux through a single zonal or meridional face and its partial
!! derivative with the face velocity, from a PPM reconstruction of the thicknesses on either side.
subroutine flux_elem(u, h, h_p1, h_L, h_L_p1, h_R, h_R_p1, uh, duhdu, visc_rem, &
                     G_dy_Cu, G_IareaT, G_IareaT_p1, G_IdxT, G_IdxT_p1, dt, &
                     vol_CFL, por_face_area, h_marg_min)
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
subroutine flux_elem_OBC(u, h, h_p1, uh, duhdu, visc_rem, por_face_area, G_dy_Cu, h_marg_min, OBC, l_seg)
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
  type(ocean_OBC_type), intent(in)    :: OBC      !< Open boundaries control structure.
  integer,              intent(in)    :: l_seg    !< The signed segment number of the face, or 0.

  if (l_seg /= 0) then
    if (OBC%segment(abs(l_seg))%open) then
      if (l_seg > 0) then !  OBC_DIRECTION_E or OBC_DIRECTION_N
        uh = (G_dy_Cu * por_face_area) * u * h
        duhdu = (G_dy_Cu * por_face_area) * max(h, h_marg_min) * visc_rem
      else !  OBC_DIRECTION_W or OBC_DIRECTION_S
        uh = (G_dy_Cu * por_face_area) * u * h_p1
        duhdu = (G_dy_Cu * por_face_area) * max(h_p1, h_marg_min) * visc_rem
      endif
    endif
  endif

end subroutine flux_elem_OBC

module procedure ratio_max
  if (abs(a) > abs(maxrat*b)) then
    ratio = maxrat
  else
    ratio = a / b
  endif
end procedure ratio_max

end submodule MOM_continuity_PPM_s
