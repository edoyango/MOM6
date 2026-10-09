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
  real, dimension(SZIB_(G),SZK_(GV)) :: duhdu ! Partial derivative of uh with u [H L ~> m2 or kg m-1].
  real, dimension(SZIB_(G)) :: &
    du, &         ! Corrective barotropic change in the velocity to give uhbt [L T-1 ~> m s-1].
    du_min_CFL, & ! Lower limit on du correction to avoid CFL violations [L T-1 ~> m s-1]
    du_max_CFL, & ! Upper limit on du correction to avoid CFL violations [L T-1 ~> m s-1]
    duhdu_tot_0, & ! Summed partial derivative of uh with u [H L ~> m2 or kg m-1].
    uh_tot_0, &   ! Summed transport with no barotropic correction [H L2 T-1 ~> m3 s-1 or kg s-1].
    visc_rem_max  ! The column maximum of visc_rem [nondim].
  logical, dimension(SZIB_(G)) :: do_I
  real, dimension(SZIB_(G),SZK_(GV)) :: &
    visc_rem      ! A 2-D copy of visc_rem_u or an array of 1's [nondim].
  real, dimension(SZIB_(G)) :: FAuI  ! A list of sums of zonal face areas [H L ~> m2 or kg m-1].
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
  logical :: simple_OBC_pt(SZIB_(G))  ! Indicates points in a row with specified transport OBCs

  call cpu_clock_begin(id_clock_correct)

  use_visc_rem = present(visc_rem_u)

  set_BT_cont = .false. ; if (present(BT_cont)) set_BT_cont = (associated(BT_cont))

  local_specified_BC = .false. ; local_Flather_OBC = .false. ; local_open_BC = .false.
  if (associated(OBC)) then ; if (OBC%OBC_pe) then
    local_specified_BC = OBC%specified_u_BCs_exist_globally
    local_Flather_OBC = OBC%Flather_u_BCs_exist_globally
    local_open_BC = OBC%open_u_BCs_exist_globally
  endif ; endif

  if (present(du_cor)) du_cor(:,:) = 0.0

  if (present(LB_in)) then
    LB = LB_in
  else
    LB%ish = G%isc ; LB%ieh = G%iec ; LB%jsh = G%jsc ; LB%jeh = G%jec
  endif
  ish = LB%ish ; ieh = LB%ieh ; jsh = LB%jsh ; jeh = LB%jeh ; nz = GV%ke

  CFL_dt = CS%CFL_limit_adjust / dt
  I_dt = 1.0 / dt
  if (CS%aggress_adjust) CFL_dt = I_dt

  if (.not.use_visc_rem) visc_rem(:,:) = 1.0
  !$OMP parallel do default(shared) private(do_I,duhdu,du,du_max_CFL,du_min_CFL,uh_tot_0, &
  !$OMP                                     duhdu_tot_0,FAuI,visc_rem_max,I_vrm,du_lim,dx_E,dx_W, &
  !$OMP                                     simple_OBC_pt,any_simple_OBC,l_seg) &
  !$OMP                        firstprivate(visc_rem)
  do j=jsh,jeh
    do I=ish-1,ieh ; do_I(I) = .true. ; enddo
    ! Set uh and duhdu.
    do k=1,nz
      if (use_visc_rem) then ; do I=ish-1,ieh
        visc_rem(I,k) = visc_rem_u(I,j,k)
      enddo ; endif
      call zonal_flux_layer(u(:,j,k), h_in(:,j,k), h_W(:,j,k), h_E(:,j,k), &
                            uh(:,j,k), duhdu(:,k), visc_rem(:,k), &
                            dt, G, US, j, ish, ieh, do_I, CS%vol_CFL, por_face_areaU(:,j,k), &
                            CS%h_marg_min, OBC)
      if (local_specified_BC) then
        do I=ish-1,ieh ; if (OBC%segnum_u(I,j) /= 0) then
          l_seg = abs(OBC%segnum_u(I,j))
          if (OBC%segment(l_seg)%specified) uh(I,j,k) = OBC%segment(l_seg)%normal_trans(I,j,k)
        endif ; enddo
      endif
    enddo

    if (present(uhbt) .or. set_BT_cont) then
      if (use_visc_rem .and. CS%use_visc_rem_max) then
        visc_rem_max(:) = 0.0
        do k=1,nz ; do I=ish-1,ieh
          visc_rem_max(I) = max(visc_rem_max(I), visc_rem(I,k))
        enddo ; enddo
      else
        visc_rem_max(:) = 1.0
      endif
      !   Set limits on du that will keep the CFL number between -1 and 1.
      ! This should be adequate to keep the root bracketed in all cases.
      do I=ish-1,ieh
        I_vrm = 0.0
        if (visc_rem_max(I) > 0.0) I_vrm = 1.0 / visc_rem_max(I)
        if (CS%vol_CFL) then
          dx_W = ratio_max(G%areaT(i,j), G%dy_Cu(I,j), 1000.0*G%dxT(i,j))
          dx_E = ratio_max(G%areaT(i+1,j), G%dy_Cu(I,j), 1000.0*G%dxT(i+1,j))
        else ; dx_W = G%dxT(i,j) ; dx_E = G%dxT(i+1,j) ; endif
        du_max_CFL(I) = 2.0* (CFL_dt * dx_W) * I_vrm
        du_min_CFL(I) = -2.0 * (CFL_dt * dx_E) * I_vrm
        uh_tot_0(I) = 0.0 ; duhdu_tot_0(I) = 0.0
      enddo
      do k=1,nz ; do I=ish-1,ieh
        duhdu_tot_0(I) = duhdu_tot_0(I) + duhdu(I,k)
        uh_tot_0(I) = uh_tot_0(I) + uh(I,j,k)
      enddo ; enddo
      if (use_visc_rem) then
        if (CS%aggress_adjust) then
          do k=1,nz ; do I=ish-1,ieh
            if (CS%vol_CFL) then
              dx_W = ratio_max(G%areaT(i,j), G%dy_Cu(I,j), 1000.0*G%dxT(i,j))
              dx_E = ratio_max(G%areaT(i+1,j), G%dy_Cu(I,j), 1000.0*G%dxT(i+1,j))
            else ; dx_W = G%dxT(i,j) ; dx_E = G%dxT(i+1,j) ; endif

            du_lim = 0.499*((dx_W*I_dt - u(I,j,k)) + MIN(0.0,u(I-1,j,k)))
            if (du_max_CFL(I) * visc_rem(I,k) > du_lim) &
              du_max_CFL(I) = du_lim / visc_rem(I,k)

            du_lim = 0.499*((-dx_E*I_dt - u(I,j,k)) + MAX(0.0,u(I+1,j,k)))
            if (du_min_CFL(I) * visc_rem(I,k) < du_lim) &
              du_min_CFL(I) = du_lim / visc_rem(I,k)
          enddo ; enddo
        else
          do k=1,nz ; do I=ish-1,ieh
            if (CS%vol_CFL) then
              dx_W = ratio_max(G%areaT(i,j), G%dy_Cu(I,j), 1000.0*G%dxT(i,j))
              dx_E = ratio_max(G%areaT(i+1,j), G%dy_Cu(I,j), 1000.0*G%dxT(i+1,j))
            else ; dx_W = G%dxT(i,j) ; dx_E = G%dxT(i+1,j) ; endif

            if (du_max_CFL(I) * visc_rem(I,k) > dx_W*CFL_dt - u(I,j,k)*G%mask2dCu(I,j)) &
              du_max_CFL(I) = (dx_W*CFL_dt - u(I,j,k)) / visc_rem(I,k)
            if (du_min_CFL(I) * visc_rem(I,k) < -dx_E*CFL_dt - u(I,j,k)*G%mask2dCu(I,j)) &
              du_min_CFL(I) = -(dx_E*CFL_dt + u(I,j,k)) / visc_rem(I,k)
          enddo ; enddo
        endif
      else
        if (CS%aggress_adjust) then
          do k=1,nz ; do I=ish-1,ieh
            if (CS%vol_CFL) then
              dx_W = ratio_max(G%areaT(i,j), G%dy_Cu(I,j), 1000.0*G%dxT(i,j))
              dx_E = ratio_max(G%areaT(i+1,j), G%dy_Cu(I,j), 1000.0*G%dxT(i+1,j))
            else ; dx_W = G%dxT(i,j) ; dx_E = G%dxT(i+1,j) ; endif

            du_max_CFL(I) = MIN(du_max_CFL(I), 0.499 * &
                        ((dx_W*I_dt - u(I,j,k)) + MIN(0.0,u(I-1,j,k))) )
            du_min_CFL(I) = MAX(du_min_CFL(I), 0.499 * &
                        ((-dx_E*I_dt - u(I,j,k)) + MAX(0.0,u(I+1,j,k))) )
          enddo ; enddo
        else
          do k=1,nz ; do I=ish-1,ieh
            if (CS%vol_CFL) then
              dx_W = ratio_max(G%areaT(i,j), G%dy_Cu(I,j), 1000.0*G%dxT(i,j))
              dx_E = ratio_max(G%areaT(i+1,j), G%dy_Cu(I,j), 1000.0*G%dxT(i+1,j))
            else ; dx_W = G%dxT(i,j) ; dx_E = G%dxT(i+1,j) ; endif

            du_max_CFL(I) = MIN(du_max_CFL(I), dx_W*CFL_dt - u(I,j,k))
            du_min_CFL(I) = MAX(du_min_CFL(I), -(dx_E*CFL_dt + u(I,j,k)))
          enddo ; enddo
        endif
      endif
      do I=ish-1,ieh
        du_max_CFL(I) = max(du_max_CFL(I),0.0)
        du_min_CFL(I) = min(du_min_CFL(I),0.0)
      enddo

      any_simple_OBC = .false.
      if (present(uhbt) .or. set_BT_cont) then
        if (local_specified_BC .or. local_Flather_OBC) then ; do I=ish-1,ieh
          l_seg = abs(OBC%segnum_u(I,j))

          ! Avoid reconciling barotropic/baroclinic transports if transport is specified
          simple_OBC_pt(I) = .false.
          if (l_seg /= OBC_NONE) simple_OBC_pt(I) = OBC%segment(l_seg)%specified
          do_I(I) = .not.simple_OBC_pt(I)
          any_simple_OBC = any_simple_OBC .or. simple_OBC_pt(I)
        enddo ; else ; do I=ish-1,ieh
          do_I(I) = .true.
        enddo ; endif
      endif

      if (present(uhbt)) then
        ! Find du and uh.
        call zonal_flux_adjust(u, h_in, h_W, h_E, uhbt(:,j), uh_tot_0, duhdu_tot_0, du, &
                               du_max_CFL, du_min_CFL, dt, G, GV, US, CS, visc_rem, &
                               j, ish, ieh, do_I, por_face_areaU, uh, OBC=OBC)

        if (present(u_cor)) then ; do k=1,nz
          do I=ish-1,ieh ; u_cor(I,j,k) = u(I,j,k) + du(I) * visc_rem(I,k) ; enddo
          if (any_simple_OBC) then ; do I=ish-1,ieh ; if (simple_OBC_pt(I)) then
            u_cor(I,j,k) = OBC%segment(abs(OBC%segnum_u(I,j)))%normal_vel(I,j,k)
          endif ; enddo ; endif
        enddo ; endif ! u-corrected

        if (present(du_cor)) then
          do I=ish-1,ieh ; du_cor(I,j) = du(I) ; enddo
        endif

      endif

      if (set_BT_cont) then
        call set_zonal_BT_cont(u, h_in, h_W, h_E, BT_cont, uh_tot_0, duhdu_tot_0,&
                               du_max_CFL, du_min_CFL, dt, G, GV, US, CS, visc_rem, &
                               visc_rem_max, j, ish, ieh, do_I, por_face_areaU)
        if (any_simple_OBC) then
          do I=ish-1,ieh
            if (simple_OBC_pt(I)) FAuI(I) = GV%H_subroundoff*G%dy_Cu(I,j)
          enddo
          ! NOTE: simple_OBC_pt(I) should prevent access to segment OBC_NONE
          do k=1,nz ; do I=ish-1,ieh ; if (simple_OBC_pt(I)) then
            l_seg = abs(OBC%segnum_u(I,j))
            if ((abs(OBC%segment(l_seg)%normal_vel(I,j,k)) > 0.0) .and. (OBC%segment(l_seg)%specified)) &
              FAuI(I) = FAuI(I) + OBC%segment(l_seg)%normal_trans(I,j,k) / OBC%segment(l_seg)%normal_vel(I,j,k)
          endif ; enddo ; enddo
          do I=ish-1,ieh ; if (simple_OBC_pt(I)) then
            BT_cont%FA_u_W0(I,j) = FAuI(I) ; BT_cont%FA_u_E0(I,j) = FAuI(I)
            BT_cont%FA_u_WW(I,j) = FAuI(I) ; BT_cont%FA_u_EE(I,j) = FAuI(I)
            BT_cont%uBT_WW(I,j) = 0.0 ; BT_cont%uBT_EE(I,j) = 0.0
          endif ; enddo
        endif
      endif ! set_BT_cont

    endif ! present(uhbt) or set_BT_cont

  enddo ! j-loop

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

module procedure zonal_flux_adjust
  ! Local variables
  real, dimension(SZIB_(G),SZK_(GV)) :: &
    uh_aux, &  ! An auxiliary zonal volume flux [H L2 T-1 ~> m3 s-1 or kg s-1].
    duhdu      ! Partial derivative of uh with u [H L ~> m2 or kg m-1].
  real, dimension(SZIB_(G)) :: &
    uh_err, &  ! Difference between uhbt and the summed uh [H L2 T-1 ~> m3 s-1 or kg s-1].
    uh_err_best, & ! The smallest value of uh_err found so far [H L2 T-1 ~> m3 s-1 or kg s-1].
    u_new, &   ! The velocity with the correction added [L T-1 ~> m s-1].
    duhdu_tot,&! Summed partial derivative of uh with u [H L ~> m2 or kg m-1].
    du_min, &  ! Lower limit on du correction based on CFL limits and previous iterations [L T-1 ~> m s-1]
    du_max     ! Upper limit on du correction based on CFL limits and previous iterations [L T-1 ~> m s-1]
  real :: du_prev ! The previous value of du [L T-1 ~> m s-1].
  real :: ddu     ! The change in du from the previous iteration [L T-1 ~> m s-1].
  real :: tol_eta ! The tolerance for the current iteration [H ~> m or kg m-2].
  real :: tol_vel ! The tolerance for velocity in the current iteration [L T-1 ~> m s-1].
  integer :: i, k, nz, itt, max_itts = 20
  logical :: domore, do_I(SZIB_(G))

  nz = GV%ke

  uh_aux(:,:) = 0.0 ; duhdu(:,:) = 0.0

  if (present(uh_3d)) then ; do k=1,nz ; do I=ish-1,ieh
    uh_aux(i,k) = uh_3d(I,j,k)
  enddo ; enddo ; endif

  do I=ish-1,ieh
    du(I) = 0.0 ; do_I(I) = do_I_in(I)
    du_max(I) = du_max_CFL(I) ; du_min(I) = du_min_CFL(I)
    uh_err(I) = uh_tot_0(I) - uhbt(I) ; duhdu_tot(I) = duhdu_tot_0(I)
    uh_err_best(I) = abs(uh_err(I))
  enddo

  do itt=1,max_itts
    select case (itt)
      case (:1) ; tol_eta = 1e-6 * CS%tol_eta
      case (2)  ; tol_eta = 1e-4 * CS%tol_eta
      case (3)  ; tol_eta = 1e-2 * CS%tol_eta
      case default ; tol_eta = CS%tol_eta
    end select
    tol_vel = CS%tol_vel

    do I=ish-1,ieh
      if (uh_err(I) > 0.0) then ; du_max(I) = du(I)
      elseif (uh_err(I) < 0.0) then ; du_min(I) = du(I)
      else ; do_I(I) = .false. ; endif
    enddo
    domore = .false.
    do I=ish-1,ieh ; if (do_I(I)) then
      if ((dt * min(G%IareaT(i,j),G%IareaT(i+1,j))*abs(uh_err(I)) > tol_eta) .or. &
          (CS%better_iter .and. ((abs(uh_err(I)) > tol_vel * duhdu_tot(I)) .or. &
                                 (abs(uh_err(I)) > uh_err_best(I))) )) then
      !   Use Newton's method, provided it stays bounded.  Otherwise bisect
      ! the value with the appropriate bound.
        ddu = -uh_err(I) / duhdu_tot(I)
        du_prev = du(I)
        du(I) = du(I) + ddu
        if (abs(ddu) < 1.0e-15*abs(du(I))) then
          do_I(I) = .false. ! ddu is small enough to quit.
        elseif (ddu > 0.0) then
          if (du(I) >= du_max(I)) then
            du(I) = 0.5*(du_prev + du_max(I))
            if (du_max(I) - du_prev < 1.0e-15*abs(du(I))) do_I(I) = .false.
          endif
        else ! ddu < 0.0
          if (du(I) <= du_min(I)) then
            du(I) = 0.5*(du_prev + du_min(I))
            if (du_prev - du_min(I) < 1.0e-15*abs(du(I))) do_I(I) = .false.
          endif
        endif
        if (do_I(I)) domore = .true.
      else
        do_I(I) = .false.
      endif
    endif ; enddo
    if (.not.domore) exit

    if ((itt < max_itts) .or. present(uh_3d)) then ; do k=1,nz
      do I=ish-1,ieh ; u_new(I) = u(I,j,k) + du(I) * visc_rem(I,k) ; enddo
      call zonal_flux_layer(u_new, h_in(:,j,k), h_W(:,j,k), h_E(:,j,k), &
                            uh_aux(:,k), duhdu(:,k), visc_rem(:,k), &
                            dt, G, US, j, ish, ieh, do_I, CS%vol_CFL, por_face_areaU(:,j,k), &
                            CS%h_marg_min, OBC)
    enddo ; endif

    if (itt < max_itts) then
      do I=ish-1,ieh
        uh_err(I) = -uhbt(I) ; duhdu_tot(I) = 0.0
      enddo
      do k=1,nz ; do I=ish-1,ieh
        uh_err(I) = uh_err(I) + uh_aux(I,k)
        duhdu_tot(I) = duhdu_tot(I) + duhdu(I,k)
      enddo ; enddo
      do I=ish-1,ieh
        uh_err_best(I) = min(uh_err_best(I), abs(uh_err(I)))
      enddo
    endif
  enddo ! itt-loop
  ! If there are any faces which have not converged to within the tolerance,
  ! so-be-it, or else use a final upwind correction?
  ! This never seems to happen with 20 iterations as max_itt.

  if (present(uh_3d)) then ; do k=1,nz ; do I=ish-1,ieh
    uh_3d(I,j,k) = uh_aux(I,k)
  enddo ; enddo ; endif

end procedure zonal_flux_adjust

module procedure set_zonal_BT_cont
  ! Local variables
  real, dimension(SZIB_(G)) :: &
    du0, &        ! The barotropic velocity increment that gives 0 transport [L T-1 ~> m s-1].
    duL, duR, &   ! The barotropic velocity increments that give the westerly
                  ! (duL) and easterly (duR) test velocities [L T-1 ~> m s-1].
    zeros, &      ! An array of full of 0 transports [H L2 T-1 ~> m3 s-1 or kg s-1]
    du_CFL, &     ! The velocity increment that corresponds to CFL_min [L T-1 ~> m s-1].
    u_L, u_R, &   ! The westerly (u_L), easterly (u_R), and zero-barotropic
    u_0, &        ! transport (u_0) layer test velocities [L T-1 ~> m s-1].
    duhdu_L, &    ! The effective layer marginal face areas with the westerly
    duhdu_R, &    ! (_L), easterly (_R), and zero-barotropic (_0) test
    duhdu_0, &    ! velocities [H L ~> m2 or kg m-1].
    uh_L, uh_R, & ! The layer transports with the westerly (_L), easterly (_R),
    uh_0, &       ! and zero-barotropic (_0) test velocities [H L2 T-1 ~> m3 s-1 or kg s-1].
    FAmt_L, FAmt_R, & ! The summed effective marginal face areas for the 3
    FAmt_0, &     ! test velocities [H L ~> m2 or kg m-1].
    uhtot_L, &    ! The summed transport with the westerly (uhtot_L) and
    uhtot_R       ! and easterly (uhtot_R) test velocities [H L2 T-1 ~> m3 s-1 or kg s-1].
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
  logical :: domore
  integer :: i, k, nz

  nz = GV%ke ; Idt = 1.0 / dt
  min_visc_rem = 0.1 ; CFL_min = 1e-6

 ! Diagnose the zero-transport correction, du0.
  do I=ish-1,ieh ; zeros(I) = 0.0 ; enddo
  call zonal_flux_adjust(u, h_in, h_W, h_E, zeros, uh_tot_0, duhdu_tot_0, du0, &
                         du_max_CFL, du_min_CFL, dt, G, GV, US, CS, visc_rem, &
                         j, ish, ieh, do_I, por_face_areaU)

  ! Determine the westerly- and easterly- fluxes.  Choose a sufficiently
  ! negative velocity correction for the easterly-flux, and a sufficiently
  ! positive correction for the westerly-flux.
  domore = .false.
  do I=ish-1,ieh
    if (do_I(I)) domore = .true.
    du_CFL(I) = (CFL_min * Idt) * G%dxCu(I,j)
    duR(I) = min(0.0,du0(I) - du_CFL(I))
    duL(I) = max(0.0,du0(I) + du_CFL(I))
    FAmt_L(I) = 0.0 ; FAmt_R(I) = 0.0 ; FAmt_0(I) = 0.0
    uhtot_L(I) = 0.0 ; uhtot_R(I) = 0.0
  enddo

  if (.not.domore) then
    do k=1,nz ; do I=ish-1,ieh
      BT_cont%FA_u_W0(I,j) = 0.0 ; BT_cont%FA_u_WW(I,j) = 0.0
      BT_cont%FA_u_E0(I,j) = 0.0 ; BT_cont%FA_u_EE(I,j) = 0.0
      BT_cont%uBT_WW(I,j) = 0.0 ; BT_cont%uBT_EE(I,j) = 0.0
    enddo ; enddo
    return
  endif

  do k=1,nz ; do I=ish-1,ieh ; if (do_I(I)) then
    visc_rem_lim = max(visc_rem(I,k), min_visc_rem*visc_rem_max(I))
    if (visc_rem_lim > 0.0) then ! This is almost always true for ocean points.
      if (u(I,j,k) + duR(I)*visc_rem_lim > -du_CFL(I)*visc_rem(I,k)) &
        duR(I) = -(u(I,j,k) + du_CFL(I)*visc_rem(I,k)) / visc_rem_lim
      if (u(I,j,k) + duL(I)*visc_rem_lim < du_CFL(I)*visc_rem(I,k)) &
        duL(I) = -(u(I,j,k) - du_CFL(I)*visc_rem(I,k)) / visc_rem_lim
    endif
  endif ; enddo ; enddo

  do k=1,nz
    do I=ish-1,ieh ; if (do_I(I)) then
      u_L(I) = u(I,j,k) + duL(I) * visc_rem(I,k)
      u_R(I) = u(I,j,k) + duR(I) * visc_rem(I,k)
      u_0(I) = u(I,j,k) + du0(I) * visc_rem(I,k)
    endif ; enddo
    call zonal_flux_layer(u_0, h_in(:,j,k), h_W(:,j,k), h_E(:,j,k), uh_0, duhdu_0, &
                          visc_rem(:,k), dt, G, US, j, ish, ieh, do_I, CS%vol_CFL, &
                          por_face_areaU(:,j,k), CS%h_marg_min)
    call zonal_flux_layer(u_L, h_in(:,j,k), h_W(:,j,k), h_E(:,j,k), uh_L, duhdu_L, &
                          visc_rem(:,k), dt, G, US, j, ish, ieh, do_I, CS%vol_CFL, &
                          por_face_areaU(:,j,k), CS%h_marg_min)
    call zonal_flux_layer(u_R, h_in(:,j,k), h_W(:,j,k), h_E(:,j,k), uh_R, duhdu_R, &
                          visc_rem(:,k), dt, G, US, j, ish, ieh, do_I, CS%vol_CFL, &
                          por_face_areaU(:,j,k), CS%h_marg_min)
    do I=ish-1,ieh ; if (do_I(I)) then
      FAmt_0(I) = FAmt_0(I) + duhdu_0(I)
      FAmt_L(I) = FAmt_L(I) + duhdu_L(I)
      FAmt_R(I) = FAmt_R(I) + duhdu_R(I)
      uhtot_L(I) = uhtot_L(I) + uh_L(I)
      uhtot_R(I) = uhtot_R(I) + uh_R(I)
    endif ; enddo
  enddo
  do I=ish-1,ieh ; if (do_I(I)) then
    FA_0 = FAmt_0(I) ; FA_avg = FAmt_0(I)
    if ((duL(I) - du0(I)) /= 0.0) &
      FA_avg = uhtot_L(I) / (duL(I) - du0(I))
    if (FA_avg > max(FA_0, FAmt_L(I))) then ; FA_avg = max(FA_0, FAmt_L(I))
    elseif (FA_avg < min(FA_0, FAmt_L(I))) then ; FA_0 = FA_avg ; endif

    BT_cont%FA_u_W0(I,j) = FA_0 ; BT_cont%FA_u_WW(I,j) = FAmt_L(I)
    if (abs(FA_0-FAmt_L(I)) <= 1e-12*FA_0) then ; BT_cont%uBT_WW(I,j) = 0.0 ; else
      BT_cont%uBT_WW(I,j) = (1.5 * (duL(I) - du0(I))) * &
                            ((FAmt_L(I) - FA_avg) / (FAmt_L(I) - FA_0))
    endif

    FA_0 = FAmt_0(I) ; FA_avg = FAmt_0(I)
    if ((duR(I) - du0(I)) /= 0.0) &
      FA_avg = uhtot_R(I) / (duR(I) - du0(I))
    if (FA_avg > max(FA_0, FAmt_R(I))) then ; FA_avg = max(FA_0, FAmt_R(I))
    elseif (FA_avg < min(FA_0, FAmt_R(I))) then ; FA_0 = FA_avg ; endif

    BT_cont%FA_u_E0(I,j) = FA_0 ; BT_cont%FA_u_EE(I,j) = FAmt_R(I)
    if (abs(FAmt_R(I) - FA_0) <= 1e-12*FA_0) then ; BT_cont%uBT_EE(I,j) = 0.0 ; else
      BT_cont%uBT_EE(I,j) = (1.5 * (duR(I) - du0(I))) * &
                            ((FAmt_R(I) - FA_avg) / (FAmt_R(I) - FA_0))
    endif
  else
    BT_cont%FA_u_W0(I,j) = 0.0 ; BT_cont%FA_u_WW(I,j) = 0.0
    BT_cont%FA_u_E0(I,j) = 0.0 ; BT_cont%FA_u_EE(I,j) = 0.0
    BT_cont%uBT_WW(I,j) = 0.0 ; BT_cont%uBT_EE(I,j) = 0.0
  endif ; enddo

end procedure set_zonal_BT_cont

module procedure ratio_max
  if (abs(a) > abs(maxrat*b)) then
    ratio = maxrat
  else
    ratio = a / b
  endif
end procedure ratio_max

end submodule MOM_continuity_PPM_s
