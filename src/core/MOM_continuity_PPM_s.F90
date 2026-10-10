! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> The zonal and meridional mass flux routines of MOM_continuity_PPM and the routines they call, kept in
!! a submodule of their own so that changing them does not force the modules that use
!! MOM_continuity_PPM to be recompiled.
submodule (MOM_continuity_PPM) MOM_continuity_PPM_s

#include <MOM_memory.h>

implicit none

contains

module procedure zonal_mass_flux
  ! Local variables
  ! The block-sized work arrays below are dimensioned (nIIB,njj[,nz]).  In this zonal
  ! routine their first (nIIB) block index is staggered at the zonal velocity (u, C-grid I)
  ! points, while the second (njj) block index is at the tracer-row (j) points.
  real, dimension(nIIB,njj,SZK_(GV)) :: &
    duhdu,      & ! Partial derivative of uh with u, at u points [H L ~> m2 or kg m-1].
    visc_rem      ! A 2-D copy of visc_rem_u or an array of 1's, at u points [nondim].
  real, dimension(nIIB,njj) :: &
    du, &         ! Corrective barotropic change in the velocity to give uhbt [L T-1 ~> m s-1].
    du_min_CFL, & ! Lower limit on du correction to avoid CFL violations [L T-1 ~> m s-1]
    du_max_CFL, & ! Upper limit on du correction to avoid CFL violations [L T-1 ~> m s-1]
    duhdu_tot_0, & ! Summed partial derivative of uh with u [H L ~> m2 or kg m-1].
    uh_tot_0, &   ! Summed transport with no barotropic correction [H L2 T-1 ~> m3 s-1 or kg s-1].
    visc_rem_max  ! The column maximum of visc_rem [nondim].
  logical, dimension(nIIB,njj) :: do_ij
  real, dimension(nIIB,njj,SZK_(GV)) :: &
    uh_b          ! A block-sized copy of uh [H L2 T-1 ~> m3 s-1 or kg s-1].
  real, dimension(nIIB,njj) :: &
    uhbt_b        ! A block-sized copy of uhbt [H L2 T-1 ~> m3 s-1 or kg s-1].
  real, dimension(nIIB,njj) :: FAuI  ! A list of sums of zonal face areas [H L ~> m2 or kg m-1].
  real :: FA_u    ! A sum of zonal face areas [H L ~> m2 or kg m-1].
  real :: I_vrm   ! 1.0 / visc_rem_max [nondim]
  real :: CFL_dt  ! The maximum CFL ratio of the adjusted velocities divided by
                  ! the time step [T-1 ~> s-1].
  real :: I_dt    ! 1.0 / dt [T-1 ~> s-1].
  real :: du_lim  ! The velocity change that give a relative CFL of 1 [L T-1 ~> m s-1].
  real :: dx_E, dx_W ! Effective x-grid spacings to the east and west [L ~> m].
  type(cont_loop_bounds_type) :: LB
  integer :: i, j, k, ish, ieh, jsh, jeh, n, nz, II, jj
  integer :: l_seg ! The OBC segment number
  logical :: use_visc_rem, set_BT_cont
  logical :: local_specified_BC, local_Flather_OBC, local_open_BC, any_simple_OBC  ! OBC-related logicals
  logical :: simple_OBC_pt(nIIB,njj)  ! Indicates points in a row with specified transport OBCs
  integer :: jsb, jeb !< The j-index bounds of the current block.
  integer :: IsbB, IebB !< The I-index bounds of the current block.
  integer :: jje, IIe !< The block-local j- and I-index end bounds.
  integer :: nteams
  logical :: set_h_u
  real :: h_marg_min ! A copy of CS%h_marg_min for use in kernels [H ~> m or kg m-2].
  logical :: vol_CFL ! A copy of CS%vol_CFL for use in kernels.

  use_visc_rem = present(visc_rem_u)

  set_BT_cont = .false. ; if (present(BT_cont)) set_BT_cont = (associated(BT_cont))
  set_h_u = .false. ; if (set_BT_cont) set_h_u = allocated(BT_cont%h_u)

  !   The arrays that are used on the device are copied in before the clock is started.  The output
  ! arrays that are only partly set here are copied in too, so that the copies back to the host
  ! leave the rest of them unchanged.
  !$omp target enter data map(to: G)
  !$omp target enter data map(to: G%dy_Cu, G%dxCu, G%IareaT, G%areaT, G%IdxT, G%dxT, G%mask2dCu)
  !$omp target enter data map(to: CS)
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
  !$omp                              uh_tot_0, visc_rem_max, do_ij, uh_b, uhbt_b)

  local_specified_BC = .false. ; local_Flather_OBC = .false. ; local_open_BC = .false.
  if (associated(OBC)) then ; if (OBC%OBC_pe) then
    local_specified_BC = OBC%specified_u_BCs_exist_globally
    local_Flather_OBC = OBC%Flather_u_BCs_exist_globally
    local_open_BC = OBC%open_u_BCs_exist_globally
  endif ; endif

  if (present(du_cor)) then
    !$omp target teams distribute parallel do collapse(2)
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
  h_marg_min = CS%h_marg_min ; vol_CFL = CS%vol_CFL
  I_dt = 1.0 / dt
  if (CS%aggress_adjust) CFL_dt = I_dt

  if (.not.use_visc_rem) then
    !$omp target teams distribute parallel do collapse(3)
    do k=1,nz ; do jj=1,njj ; do II=1,nIIB
      visc_rem(II,jj,k) = 1.0
    enddo ; enddo ; enddo
  endif

  do jsb=jsh,jeh,njj ; do IsbB=ish-1,ieh,nIIB
    IebB = min(IsbB + nIIB - 1, ieh)
    jeb = min(jsb + njj - 1, jeh)
    IIe = IebB - IsbB + 1
    jje = jeb - jsb + 1
    !$omp target teams distribute parallel do collapse(2)
    do jj=1,jje ; do II=1,IIe
      do_ij(II,jj) = .true.
    enddo ; enddo

    ! visc_rem_u is copied in a kernel of its own, as a kernel that names an absent argument fails.
    if (use_visc_rem) then
      !$omp target teams distribute parallel do collapse(3) private(I,j)
      do k=1,nz ; do jj=1,jje ; do II=1,IIe
        I = IsbB + II - 1
        j = jsb + jj - 1
        visc_rem(II,jj,k) = visc_rem_u(I,j,k)
      enddo ; enddo ; enddo
    endif

    ! Set uh and duhdu.
    !$omp target teams distribute parallel do collapse(3) private(I,j,l_seg) &
    !$omp   firstprivate(IsbB, jsb, dt, local_open_BC, local_specified_BC, vol_CFL, h_marg_min)
    do k=1,nz ; do jj=1,jje ; do II=1,IIe
      I = IsbB + II - 1
      j = jsb + jj - 1

        call flux_elem(u(I,j,k), h_in(i,j,k), h_in(i+1,j,k), &
            h_W(i,j,k), h_W(i+1,j,k), h_E(i,j,k), h_E(i+1,j,k), &
            uh_b(II,jj,k), duhdu(II,jj,k), visc_rem(II,jj,k), G%dy_Cu(I,j), &
            G%IareaT(i,j), G%IareaT(i+1,j), G%IdxT(i,j), G%IdxT(i+1,j), dt, &
            vol_CFL, por_face_areaU(I,j,k), h_marg_min)

        if (local_open_BC) &
          call flux_elem_OBC(u(I,j,k), h_in(i,j,k), h_in(i+1,j,k), &
              uh_b(II,jj,k), duhdu(II,jj,k), visc_rem(II,jj,k), &
              por_face_areaU(I,j,k), G%dy_Cu(I,j), OBC, OBC%segnum_u(I,j), &
              h_marg_min)

      if (local_specified_BC) then
          if (OBC%segnum_u(I,j) /= 0) then
            l_seg = abs(OBC%segnum_u(I,j))
            if (OBC%segment(l_seg)%specified) &
              uh_b(II,jj,k) = OBC%segment(l_seg)%normal_trans(I,j,k)
          endif
      endif
    enddo ; enddo ; enddo

    if (present(uhbt) .or. set_BT_cont) then
      if (use_visc_rem .and. CS%use_visc_rem_max) then
        !$omp target teams distribute parallel do collapse(2)
        do jj=1,jje ; do II=1,IIe
          visc_rem_max(II,jj) = 0.0

        do k=1,nz
          visc_rem_max(II,jj) = max(visc_rem_max(II,jj), visc_rem(II,jj,k))
        enddo ; enddo ; enddo
      else
        !$omp target teams distribute parallel do collapse(2)
        do jj=1,jje ; do II=1,IIe
          visc_rem_max(II,jj) = 1.0
        enddo ; enddo
      endif
      !   Set limits on du that will keep the CFL number between -1 and 1.
      ! This should be adequate to keep the root bracketed in all cases.
      !$omp target teams distribute parallel do collapse(2) private(I,j,I_vrm,dx_W,dx_E)
      do jj=1,jje ; do II=1,IIe
        I = IsbB + II - 1
        j = jsb + jj - 1

        I_vrm = 0.0
        if (visc_rem_max(II,jj) > 0.0) I_vrm = 1.0 / visc_rem_max(II,jj)
        if (CS%vol_CFL) then
          dx_W = ratio_max(G%areaT(i,j), G%dy_Cu(I,j), 1000.0*G%dxT(i,j))
          dx_E = ratio_max(G%areaT(i+1,j), G%dy_Cu(I,j), 1000.0*G%dxT(i+1,j))
        else ; dx_W = G%dxT(i,j) ; dx_E = G%dxT(i+1,j) ; endif
        du_max_CFL(II,jj) = 2.0* (CFL_dt * dx_W) * I_vrm
        du_min_CFL(II,jj) = -2.0 * (CFL_dt * dx_E) * I_vrm
        uh_tot_0(II,jj) = 0.0 ; duhdu_tot_0(II,jj) = 0.0

      do k=1,nz
          duhdu_tot_0(II,jj) = duhdu_tot_0(II,jj) + duhdu(II,jj,k)
          uh_tot_0(II,jj) = uh_tot_0(II,jj) + uh_b(II,jj,k)
      enddo
      enddo ; enddo

      if (use_visc_rem) then
        if (CS%aggress_adjust) then
          !$omp target teams distribute parallel do collapse(2) private(I,j,dx_E,Dx_W,du_lim)
            do jj=1,jje ; do II=1,IIe
              I = IsbB + II - 1
              j = jsb + jj - 1
          do k=1,nz

              if (CS%vol_CFL) then
                dx_W = ratio_max(G%areaT(i,j), G%dy_Cu(I,j), 1000.0*G%dxT(i,j))
                dx_E = ratio_max(G%areaT(i+1,j), G%dy_Cu(I,j), 1000.0*G%dxT(i+1,j))
              else ; dx_W = G%dxT(i,j) ; dx_E = G%dxT(i+1,j) ; endif

              du_lim = 0.499*((dx_W*I_dt - u(I,j,k)) + MIN(0.0,u(I-1,j,k)))
              if (du_max_CFL(II,jj) * visc_rem(II,jj,k) > du_lim) &
                du_max_CFL(II,jj) = du_lim / visc_rem(II,jj,k)

              du_lim = 0.499*((-dx_E*I_dt - u(I,j,k)) + MAX(0.0,u(I+1,j,k)))
              if (du_min_CFL(II,jj) * visc_rem(II,jj,k) < du_lim) &
                du_min_CFL(II,jj) = du_lim / visc_rem(II,jj,k)
          enddo
            enddo ; enddo
        else
          !$omp target teams distribute parallel do collapse(2) private(I,j,dx_E,dx_W)
            do jj=1,jje ; do II=1,IIe
              I = IsbB + II - 1
              j = jsb + jj - 1
          do k=1,nz

              if (CS%vol_CFL) then
                dx_W = ratio_max(G%areaT(i,j), G%dy_Cu(I,j), 1000.0*G%dxT(i,j))
                dx_E = ratio_max(G%areaT(i+1,j), G%dy_Cu(I,j), 1000.0*G%dxT(i+1,j))
              else ; dx_W = G%dxT(i,j) ; dx_E = G%dxT(i+1,j) ; endif

              if (du_max_CFL(II,jj) * visc_rem(II,jj,k) > dx_W*CFL_dt - u(I,j,k)*G%mask2dCu(I,j)) &
                du_max_CFL(II,jj) = (dx_W*CFL_dt - u(I,j,k)) / visc_rem(II,jj,k)
              if (du_min_CFL(II,jj) * visc_rem(II,jj,k) < -dx_E*CFL_dt - u(I,j,k)*G%mask2dCu(I,j)) &
                du_min_CFL(II,jj) = -(dx_E*CFL_dt + u(I,j,k)) / visc_rem(II,jj,k)
          enddo
            enddo ; enddo
        endif
      else
        if (CS%aggress_adjust) then
          !$omp target teams distribute parallel do collapse(2) private(I,j,dx_E,dx_W)
            do jj=1,jje ; do II=1,IIe
              I = IsbB + II - 1
              j = jsb + jj - 1
          do k=1,nz

              if (CS%vol_CFL) then
                dx_W = ratio_max(G%areaT(i,j), G%dy_Cu(I,j), 1000.0*G%dxT(i,j))
                dx_E = ratio_max(G%areaT(i+1,j), G%dy_Cu(I,j), 1000.0*G%dxT(i+1,j))
              else ; dx_W = G%dxT(i,j) ; dx_E = G%dxT(i+1,j) ; endif

              du_max_CFL(II,jj) = MIN(du_max_CFL(II,jj), 0.499 * &
                          ((dx_W*I_dt - u(I,j,k)) + MIN(0.0,u(I-1,j,k))) )
              du_min_CFL(II,jj) = MAX(du_min_CFL(II,jj), 0.499 * &
                          ((-dx_E*I_dt - u(I,j,k)) + MAX(0.0,u(I+1,j,k))) )
          enddo
            enddo ; enddo
        else
          !$omp target teams distribute parallel do collapse(2) private(I,j,dx_E,dx_W)
            do jj=1,jje ; do II=1,IIe
              I = IsbB + II - 1
              j = jsb + jj - 1
          do k=1,nz

              if (CS%vol_CFL) then
                dx_W = ratio_max(G%areaT(i,j), G%dy_Cu(I,j), 1000.0*G%dxT(i,j))
                dx_E = ratio_max(G%areaT(i+1,j), G%dy_Cu(I,j), 1000.0*G%dxT(i+1,j))
              else ; dx_W = G%dxT(i,j) ; dx_E = G%dxT(i+1,j) ; endif

              du_max_CFL(II,jj) = MIN(du_max_CFL(II,jj), dx_W*CFL_dt - u(I,j,k))
              du_min_CFL(II,jj) = MAX(du_min_CFL(II,jj), -(dx_E*CFL_dt + u(I,j,k)))
          enddo
            enddo ; enddo
        endif
      endif
      !$omp target teams distribute parallel do collapse(2)
      do jj=1,jje ; do II=1,IIe
        du_max_CFL(II,jj) = max(du_max_CFL(II,jj),0.0)
        du_min_CFL(II,jj) = min(du_min_CFL(II,jj),0.0)
      enddo ; enddo

      any_simple_OBC = .false.
      if (present(uhbt) .or. set_BT_cont) then
        if (local_specified_BC .or. local_Flather_OBC) then
          do jj=1,jje ; do II=1,IIe
            I = IsbB + II - 1
            j = jsb + jj - 1

            l_seg = abs(OBC%segnum_u(I,j))

            ! Avoid reconciling barotropic/baroclinic transports if transport is specified
            simple_OBC_pt(II,jj) = .false.
            if (l_seg /= OBC_NONE) simple_OBC_pt(II,jj) = OBC%segment(l_seg)%specified
            do_ij(II,jj) = .not.simple_OBC_pt(II,jj)
            any_simple_OBC = any_simple_OBC .or. simple_OBC_pt(II,jj)
          enddo ; enddo
        else
          !$omp target teams distribute parallel do collapse(2)
          do jj=1,jje ; do II=1,IIe
            do_ij(II,jj) = .true.
          enddo ; enddo
        endif
      endif

      if (present(uhbt)) then
        !$omp target teams distribute parallel do collapse(2) private(I,j)
        do jj=1,jje ; do II=1,IIe
          I = IsbB + II - 1
          j = jsb + jj - 1
          uhbt_b(II,jj) = uhbt(I,j)
        enddo ; enddo

        ! Find du and uh.
        call zonal_flux_adjust(u, h_in, h_W, h_E, uhbt_b, uh_tot_0, duhdu_tot_0, du, &
                              du_max_CFL, du_min_CFL, dt, G, GV, US, CS, visc_rem, &
                              IsbB, IebB, jsb, jeb, do_ij, por_face_areaU, nIIB, njj, &
                              uh_b, OBC=OBC)

        if (present(u_cor)) then
          !$omp target teams distribute parallel do collapse(3) private(I,j)
          do k=1,nz ; do jj=1,jje ; do II=1,IIe
            I = IsbB + II - 1
            j = jsb + jj - 1
            u_cor(I,j,k) = u(I,j,k) + du(II,jj) * visc_rem(II,jj,k)
          enddo ; enddo ; enddo

          if (any_simple_OBC) then
            do k=1,nz ; do jj=1,jje ; do II=1,IIe
              if (simple_OBC_pt(II,jj)) then
                I = IsbB + II - 1
                j = jsb + jj - 1
                u_cor(I,j,k) = OBC%segment(abs(OBC%segnum_u(I,j)))%normal_vel(I,j,k)
              endif
            enddo ; enddo ; enddo
          endif
        endif ! u-corrected

        if (present(du_cor)) then
          !$omp target teams distribute parallel do collapse(2) private(I,j)
          do jj=1,jje ; do II=1,IIe
            I = IsbB + II - 1
            j = jsb + jj - 1
            du_cor(I,j) = du(II,jj)
          enddo ; enddo
        endif

      endif

      if (set_BT_cont) then
        call set_zonal_BT_cont(u, h_in, h_W, h_E, BT_cont, uh_tot_0, duhdu_tot_0,&
                              du_max_CFL, du_min_CFL, dt, G, GV, US, CS, visc_rem, &
                              visc_rem_max, IsbB, IebB, jsb, jeb, do_ij, &
                              por_face_areaU, nIIB, njj)
        if (any_simple_OBC) then
          do jj=1,jje ; do II=1,IIe
            if (simple_OBC_pt(II,jj)) then
              I = IsbB + II - 1
              j = jsb + jj - 1
              FAuI(II,jj) = GV%H_subroundoff*G%dy_Cu(I,j)
            endif
          enddo ; enddo

          ! NOTE: simple_OBC_pt should prevent access to segment OBC_NONE
          do k=1,nz
            do jj=1,jje ; do II=1,IIe
              if (simple_OBC_pt(II,jj)) then
                I = IsbB + II - 1
                j = jsb + jj - 1

                l_seg = abs(OBC%segnum_u(I,j))
                if ((abs(OBC%segment(l_seg)%normal_vel(I,j,k)) > 0.0) .and. &
                    (OBC%segment(l_seg)%specified)) &
                  FAuI(II,jj) = FAuI(II,jj) + &
                    OBC%segment(l_seg)%normal_trans(I,j,k) / OBC%segment(l_seg)%normal_vel(I,j,k)
              endif
            enddo ; enddo
          enddo

          do jj=1,jje ; do II=1,IIe
            if (simple_OBC_pt(II,jj)) then
              I = IsbB + II - 1
              j = jsb + jj - 1
              BT_cont%FA_u_W0(I,j) = FAuI(II,jj) ; BT_cont%FA_u_E0(I,j) = FAuI(II,jj)
              BT_cont%FA_u_WW(I,j) = FAuI(II,jj) ; BT_cont%FA_u_EE(I,j) = FAuI(II,jj)
              BT_cont%uBT_WW(I,j) = 0.0 ; BT_cont%uBT_EE(I,j) = 0.0
            endif
          enddo ; enddo
        endif
      endif ! set_BT_cont
    endif ! present(uhbt) or set_BT_cont

    !$omp target teams distribute parallel do collapse(3) private(I,j)
    do k=1,nz ; do jj=1,jje ; do II=1,IIe
      I = IsbB + II - 1
      j = jsb + jj - 1
      uh(I,j,k) = uh_b(II,jj,k)
    enddo ; enddo ; enddo
  enddo ; enddo ! ij block loop

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

  !$omp target exit data map(release: duhdu, visc_rem, du, du_min_CFL, du_max_CFL, duhdu_tot_0, &
  !$omp                              uh_tot_0, visc_rem_max, do_ij, uh_b, uhbt_b)

  call cpu_clock_end(id_clock_correct)

  !   The outputs are copied back to the host, and the inputs released, after the clock is stopped.
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
  !$omp target exit data map(release: CS)
  !$omp target exit data map(release: G%dy_Cu, G%dxCu, G%IareaT, G%areaT, G%IdxT, G%dxT, G%mask2dCu)
  !$omp target exit data map(release: G)

end procedure zonal_mass_flux

module procedure flux_elem
  !DIR$ ATTRIBUTES FORCEINLINE :: flux_elem
  !$omp declare target
  real :: CFL
    ! CFL number based on the local velocity and grid spacing [nondim]
  real :: curv_3
    ! Measure of thickness curvature over a grid length [H ~> m or kg m-2]
  real :: h_marg
    ! Marginal thickness of a flux [H ~> m or kg m-2].
  real :: tmp
    ! temporary variable to store precalculated values
  real :: dh
    ! h differential between E/W

  ! Set new values of uh and duhdu.
  tmp = G_dy_Cu * por_face_area ! precalculate things

  if (u > 0.) then
    if (vol_CFL) then
      CFL = (u * dt) * (G_dy_Cu * G_IareaT)
    else
      CFL = u * dt * G_IdxT
    endif

    curv_3 = (h_L + h_R) - 2. * h
    dh = h_L - h_R
    uh = tmp * u * (h_R + CFL * (0.5 * dh + curv_3 * (CFL - 1.5)))
    h_marg = h_R + CFL * (dh + 3. * curv_3 * (CFL - 1.))
  elseif (u < 0.) then
    if (vol_CFL) then
      CFL = (-u * dt) * (G_dy_Cu * G_IareaT_p1)
    else
      CFL = -u * dt * G_IdxT_p1
    endif

    curv_3 = (h_L_p1 + h_R_p1) - 2. * h_p1
    dh = h_R_p1 - h_L_p1
    uh = tmp * u * (h_L_p1 + CFL * (0.5 * dh + curv_3 * (CFL - 1.5)))
    h_marg = h_L_p1 + CFL * (dh + 3. * curv_3 * (CFL - 1.))
  else
    uh = 0.
    h_marg = 0.5 * (h_L_p1 + h_R)
  endif

  h_marg = max(h_marg, h_marg_min)
  duhdu = tmp * h_marg * visc_rem
end procedure flux_elem

module procedure flux_elem_OBC
  !$omp declare target
  if (l_seg /= 0) then
    if (OBC%segment(abs(l_seg))%open) then
      if (l_seg > 0) then
        ! OBC_DIRECTION_E or OBC_DIRECTION_N
        uh = (G_dy_Cu * por_face_area) * u * h
        duhdu = (G_dy_Cu * por_face_area) * max(h, h_marg_min) * visc_rem
      else
        ! OBC_DIRECTION_W or OBC_DIRECTION_S
        uh = (G_dy_Cu * por_face_area) * u * h_p1
        duhdu = (G_dy_Cu * por_face_area) * max(h_p1, h_marg_min) * visc_rem
      endif
    endif
  endif
end procedure flux_elem_OBC

module procedure zonal_flux_thickness
  ! Local variables
  real :: CFL  ! The CFL number based on the local velocity and grid spacing [nondim]
  real :: curv_3 ! A measure of the thickness curvature over a grid length [H ~> m or kg m-2]
  logical :: local_open_BC
  integer :: i, j, k, ish, ieh, jsh, jeh, nz, n
  real :: dh
  ish = LB%ish ; ieh = LB%ieh ; jsh = LB%jsh ; jeh = LB%jeh ; nz = GV%ke

  !$omp target teams distribute parallel do collapse(3) private(CFL, curv_3, dh)
  do k=1,nz ; do j=jsh,jeh ; do I=ish-1,ieh
    if (u(I,j,k) > 0.0) then
      if (vol_CFL) then ; CFL = (u(I,j,k) * dt) * (G%dy_Cu(I,j) * G%IareaT(i,j))
      else ; CFL = u(I,j,k) * dt * G%IdxT(i,j) ; endif
      curv_3 = (h_W(i,j,k) + h_E(i,j,k)) - 2.0*h(i,j,k)
      dh = h_W(i,j,k) - h_E(i,j,k)
      if (marginal) then
        h_u(I,j,k) = h_E(i,j,k) + CFL * (dh + 3.0*curv_3*(CFL - 1.0))
      else
        h_u(I,j,k) = h_E(i,j,k) + CFL * (0.5*dh + curv_3*(CFL - 1.5))
      endif
    elseif (u(I,j,k) < 0.0) then
      if (vol_CFL) then ; CFL = (-u(I,j,k)*dt) * (G%dy_Cu(I,j) * G%IareaT(i+1,j))
      else ; CFL = -u(I,j,k) * dt * G%IdxT(i+1,j) ; endif
      curv_3 = (h_W(i+1,j,k) + h_E(i+1,j,k)) - 2.0*h(i+1,j,k)
      dh = h_E(i+1,j,k)-h_W(i+1,j,k)
      if (marginal) then
        h_u(I,j,k) = h_W(i+1,j,k) + CFL * (dh + 3.0*curv_3*(CFL - 1.0))
      else
        h_u(I,j,k) = h_W(i+1,j,k) + CFL * (0.5*dh + curv_3*(CFL - 1.5))
      endif
    else
      !   The choice to use the arithmetic mean here is somewhat arbitrarily, but
      ! it should be noted that h_W(i+1,j,k) and h_E(i,j,k) are usually the same.
      h_u(I,j,k) = 0.5 * (h_W(i+1,j,k) + h_E(i,j,k))
 !    h_marg = (2.0 * h_W(i+1,j,k) * h_E(i,j,k)) / &
 !             (h_W(i+1,j,k) + h_E(i,j,k) + GV%H_subroundoff)
    endif

    if (present(visc_rem_u)) then
      ! Scale back the thickness to account for the effects of viscosity and the fractional open
      ! thickness to give an appropriate non-normalized weight for each layer in determining the
      ! barotropic acceleration.
      h_u(I,j,k) = h_u(I,j,k) * (visc_rem_u(I,j,k) * por_face_areaU(I,j,k))
    else
      h_u(I,j,k) = h_u(I,j,k) * por_face_areaU(I,j,k)
    endif
  enddo ; enddo ; enddo

  local_open_BC = .false.
  if (associated(OBC)) local_open_BC = OBC%open_u_BCs_exist_globally
  if (local_open_BC) then
    do n = 1, OBC%number_of_segments
      if (OBC%segment(n)%open .and. OBC%segment(n)%is_E_or_W) then
        I = OBC%segment(n)%HI%IsdB
        if (OBC%segment(n)%direction == OBC_DIRECTION_E) then
          if (present(visc_rem_u)) then
            do k=1,nz ; do j=OBC%segment(n)%HI%jsd,OBC%segment(n)%HI%jed
              h_u(I,j,k) = h(i,j,k) * (visc_rem_u(I,j,k) * por_face_areaU(I,j,k))
            enddo ; enddo
          else
            do k=1,nz ; do j=OBC%segment(n)%HI%jsd,OBC%segment(n)%HI%jed
              h_u(I,j,k) = h(i,j,k) * por_face_areaU(I,j,k)
            enddo ; enddo
          endif
        else
          if (present(visc_rem_u)) then
            do k=1,nz ; do j=OBC%segment(n)%HI%jsd,OBC%segment(n)%HI%jed
              h_u(I,j,k) = h(i+1,j,k) * (visc_rem_u(I,j,k) * por_face_areaU(I,j,k))
            enddo ; enddo
          else
            do k=1,nz ; do j=OBC%segment(n)%HI%jsd,OBC%segment(n)%HI%jed
              h_u(I,j,k) = h(i+1,j,k) * por_face_areaU(I,j,k)
            enddo ; enddo
          endif
        endif
      endif
    enddo
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

module procedure zonal_flux_adjust
  ! Local variables
  ! The per-column iteration state is held in private scalars rather than in (nIIB,njj) arrays,
  ! so that it stays in registers instead of being loaded and stored at every level.
  real :: &
    u_new, &   ! The velocity with the correction added [L T-1 ~> m s-1].
    uh, &      ! The volume flux through a zonal face [H L2 T-1 ~> m3 s-1 or kg s-1].
    duhdu      ! Partial derivative of uh with u [H L ~> m2 or kg m-1].
  real :: du_c       ! The barotropic velocity correction for this column [L T-1 ~> m s-1].
  real :: uh_err     ! Difference between uhbt and the summed uh [H L2 T-1 ~> m3 s-1 or kg s-1].
  real :: uh_err_best ! The smallest value of uh_err found so far [H L2 T-1 ~> m3 s-1 or kg s-1].
  real :: duhdu_tot  ! Summed partial derivative of uh with u [H L ~> m2 or kg m-1].
  real :: du_min     ! Lower limit on du correction based on CFL limits and previous iterations [L T-1 ~> m s-1]
  real :: du_max     ! Upper limit on du correction based on CFL limits and previous iterations [L T-1 ~> m s-1]
  real :: dy_Cu      ! The unblocked length of this u-face [L ~> m].
  real :: IareaT_W, IareaT_E ! The inverse areas of the tracer cells west and east of this face [L-2 ~> m-2].
  real :: IdxT_W, IdxT_E     ! The inverse zonal grid spacings west and east of this face [L-1 ~> m-1].
  real :: du_prev ! The previous value of du [L T-1 ~> m s-1].
  real :: ddu     ! The change in du from the previous iteration [L T-1 ~> m s-1].
  real :: tol_eta ! The tolerance for the current iteration [H ~> m or kg m-2].
  real :: tol_eta_CS ! A copy of CS%tol_eta for use in the kernel [H ~> m or kg m-2].
  real :: tol_vel ! The tolerance for velocity in the current iteration [L T-1 ~> m s-1].
  real :: h_marg_min ! A copy of CS%h_marg_min for use in the kernel [H ~> m or kg m-2].
  logical :: vol_CFL ! A copy of CS%vol_CFL for use in the kernel.
  logical :: better_iter ! A copy of CS%better_iter for use in the kernel.
  logical :: do_c ! If true, this column is still being iterated.
  integer :: i, j, k, nz, itt, II, jj, IIe, jje
  logical :: local_open_BC ! True if there are open OBC points on this PE [nondim].
  type(ocean_OBC_type), pointer :: OBC_loc ! OBC, or null if it is absent, for use in kernels,
                                           ! as a kernel that names an absent argument fails.
  integer, parameter :: max_itts = 20

  nz = GV%ke
  IIe = IebB - IsbB + 1
  jje = jeb - jsb + 1
  local_open_BC = .false.
  OBC_loc => null()
  if (present(OBC)) then
    OBC_loc => OBC
    if (associated(OBC)) local_open_BC = OBC%open_u_BCs_exist_globally
  endif

  tol_vel = CS%tol_vel ; tol_eta_CS = CS%tol_eta ; better_iter = CS%better_iter
  h_marg_min = CS%h_marg_min ; vol_CFL = CS%vol_CFL

  !$omp target teams distribute parallel do collapse(2) &
  !$omp   private(I, j, tol_eta, ddu, du_prev, u_new, uh, duhdu, du_c, uh_err, uh_err_best, &
  !$omp           duhdu_tot, du_min, du_max, do_c, dy_Cu, IareaT_W, IareaT_E, IdxT_W, IdxT_E)
  do jj=1,jje ; do II=1,IIe
    I = IsbB + II - 1
    j = jsb + jj - 1
    dy_Cu = G%dy_Cu(I,j)
    IareaT_W = G%IareaT(i,j) ; IareaT_E = G%IareaT(i+1,j)
    IdxT_W = G%IdxT(i,j) ; IdxT_E = G%IdxT(i+1,j)
    du_c = 0.0 ; do_c = do_ij_in(II,jj)
    du_max = du_max_CFL(II,jj) ; du_min = du_min_CFL(II,jj)
    uh_err = uh_tot_0(II,jj) - uhbt(II,jj) ; duhdu_tot = duhdu_tot_0(II,jj)
    uh_err_best = abs(uh_err)

  do itt=1,max_itts
    select case (itt)
      case (:1) ; tol_eta = 1e-6 * tol_eta_CS
      case (2)  ; tol_eta = 1e-4 * tol_eta_CS
      case (3)  ; tol_eta = 1e-2 * tol_eta_CS
      case default ; tol_eta = tol_eta_CS
    end select

      if (uh_err > 0.0) then ; du_max = du_c
      elseif (uh_err < 0.0) then ; du_min = du_c
      else ; do_c = .false. ; endif

      if (do_c) then
        if ((dt * min(IareaT_W,IareaT_E)*abs(uh_err) > tol_eta) .or. &
            (better_iter .and. ((abs(uh_err) > tol_vel * duhdu_tot) .or. &
                                (abs(uh_err) > uh_err_best)) )) then
        !   Use Newton's method, provided it stays bounded.  Otherwise bisect
        ! the value with the appropriate bound.
          ddu = -uh_err / duhdu_tot
          du_prev = du_c
          du_c = du_c + ddu
          if (abs(ddu) < 1.0e-15*abs(du_c)) then
            do_c = .false. ! ddu is small enough to quit.
          elseif (ddu > 0.0) then
            if (du_c >= du_max) then
              du_c = 0.5*(du_prev + du_max)
              if (du_max - du_prev < 1.0e-15*abs(du_c)) do_c = .false.
            endif
          else ! ddu < 0.0
            if (du_c <= du_min) then
              du_c = 0.5*(du_prev + du_min)
              if (du_prev - du_min < 1.0e-15*abs(du_c)) do_c = .false.
            endif
          endif
        else
          do_c = .false.
        endif
      endif

      if (.not.do_c) exit

      uh_err = -uhbt(II,jj) ; duhdu_tot = 0.0

      do k=1,nz
        u_new = u(I,j,k) + du_c * visc_rem(II,jj,k)

        call flux_elem(u_new, h_in(i,j,k), h_in(i+1,j,k), &
            h_W(i,j,k), h_W(i+1,j,k), h_E(i,j,k), h_E(i+1,j,k), &
            uh, duhdu, visc_rem(II,jj,k), dy_Cu, &
            IareaT_W, IareaT_E, IdxT_W, IdxT_E, &
            dt, vol_CFL, por_face_areaU(I,j,k), h_marg_min)

        if (local_open_BC) &
          call flux_elem_OBC(u_new, h_in(i,j,k), h_in(i+1,j,k), &
              uh, duhdu, visc_rem(II,jj,k), &
              por_face_areaU(I,j,k), dy_Cu, OBC_loc, OBC_loc%segnum_u(I,j), &
              h_marg_min)

        uh_3d(II,jj,k) = uh
        uh_err = uh_err + uh
        duhdu_tot = duhdu_tot + duhdu
      enddo

      uh_err_best = min(uh_err_best, abs(uh_err))
    enddo ! itt-loop
    du(II,jj) = du_c
  enddo ; enddo
  ! If there are any faces which have not converged to within the tolerance,
  ! so-be-it, or else use a final upwind correction?
  ! This never seems to happen with 20 iterations as max_itt.

end procedure zonal_flux_adjust

module procedure set_zonal_BT_cont
  ! Local variables
  ! The block-sized work arrays below are dimensioned (nIIB,njj[,nz]).  In this zonal
  ! routine their first (nIIB) block index is staggered at the zonal velocity (u, C-grid I)
  ! points, while the second (njj) block index is at the tracer-row (j) points.
  real, dimension(nIIB,njj) :: &
    du0           ! The barotropic velocity increment that gives 0 transport [L T-1 ~> m s-1].
  real, dimension(nIIB,njj) :: &
    zeros         ! An array of full of 0 transports [H L2 T-1 ~> m3 s-1 or kg s-1]
  real, dimension(nIIB,njj,SZK_(GV)) :: &
    uh_tmp        ! A temporary array for the unused layer transports [H L2 T-1 ~> m3 s-1 or kg s-1]
  ! The per-column values below are private scalars in the kernel, so that they stay in registers.
  real :: &
    duL, duR, &   ! The barotropic velocity increments that give the westerly
                  ! (duL) and easterly (duR) test velocities [L T-1 ~> m s-1].
    du_CFL, &     ! The velocity increment that corresponds to CFL_min [L T-1 ~> m s-1].
    du0_c, &      ! A copy of du0 for this column [L T-1 ~> m s-1].
    FAmt_L, FAmt_R, & ! The summed effective marginal face areas for the 3
    FAmt_0, &     ! test velocities [H L ~> m2 or kg m-1].
    uhtot_L, &    ! The summed transport with the westerly (uhtot_L) and
    uhtot_R       ! and easterly (uhtot_R) test velocities [H L2 T-1 ~> m3 s-1 or kg s-1].
  real :: dy_Cu      ! The unblocked length of this u-face [L ~> m].
  real :: IareaT_W, IareaT_E ! The inverse areas of the tracer cells west and east of this face [L-2 ~> m-2].
  real :: IdxT_W, IdxT_E     ! The inverse zonal grid spacings west and east of this face [L-1 ~> m-1].
  real :: visc_rem_max_c ! A copy of visc_rem_max for this column [nondim].
  real :: h_marg_min ! A copy of CS%h_marg_min for use in the kernel [H ~> m or kg m-2].
  logical :: vol_CFL ! A copy of CS%vol_CFL for use in the kernel.
  logical :: do_c    ! A copy of do_ij for this column.
  real :: &
    u_L, u_R, &   ! The westerly (u_L), easterly (u_R), and zero-barotropic
    u_0, &        ! transport (u_0) layer test velocities [L T-1 ~> m s-1].
    duhdu_L, &    ! The effective layer marginal face areas with the westerly
    duhdu_R, &    ! (_L), easterly (_R), and zero-barotropic (_0) test
    duhdu_0, &    ! velocities [H L ~> m2 or kg m-1].
    uh_L, uh_R, & ! The layer transports with the westerly (_L), easterly (_R),
    uh_0       ! and zero-barotropic (_0) test velocities [H L2 T-1 ~> m3 s-1 or kg s-1].
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
  integer :: i, j, k, nz, II, jj, IIe, jje !< Tile loop indices.

  nz = GV%ke ; Idt = 1.0 / dt
  IIe = IebB - IsbB + 1 ; jje = jeb - jsb + 1
  min_visc_rem = 0.1 ; CFL_min = 1e-6
  h_marg_min = CS%h_marg_min ; vol_CFL = CS%vol_CFL

  !$omp target enter data map(alloc: du0, zeros, uh_tmp)

 ! Diagnose the zero-transport correction, du0.
 !$omp target teams distribute parallel do collapse(2)
  do jj=1,jje ; do II=1,IIe
    zeros(II,jj) = 0.0
  enddo ; enddo
  call zonal_flux_adjust(u, h_in, h_W, h_E, zeros, uh_tot_0, duhdu_tot_0, du0, &
                         du_max_CFL, du_min_CFL, dt, G, GV, US, CS, visc_rem, &
                         IsbB, IebB, jsb, jeb, do_ij, por_face_areaU, &
                         nIIB, njj, uh_tmp)

  ! Determine the westerly- and easterly- fluxes.  Choose a sufficiently
  ! negative velocity correction for the easterly-flux, and a sufficiently
  ! positive correction for the westerly-flux.
  !$omp target teams distribute parallel do collapse(2) &
  !$omp   private(I, j, visc_rem_lim, u_L, u_R, u_0, uh_L, uh_R, uh_0, duhdu_L, duhdu_R, duhdu_0, &
  !$omp           FA_0, FA_avg, duL, duR, du_CFL, du0_c, FAmt_L, FAmt_R, FAmt_0, uhtot_L, uhtot_R, &
  !$omp           dy_Cu, IareaT_W, IareaT_E, IdxT_W, IdxT_E, visc_rem_max_c, do_c)
  do jj=1,jje ; do II=1,IIe
    I = IsbB + II - 1
    j = jsb + jj - 1
    dy_Cu = G%dy_Cu(I,j)
    IareaT_W = G%IareaT(i,j) ; IareaT_E = G%IareaT(i+1,j)
    IdxT_W = G%IdxT(i,j) ; IdxT_E = G%IdxT(i+1,j)
    visc_rem_max_c = visc_rem_max(II,jj) ; do_c = do_ij(II,jj) ; du0_c = du0(II,jj)
    du_CFL = (CFL_min * Idt) * G%dxCu(I,j)
    duR = min(0.0,du0_c - du_CFL)
    duL = max(0.0,du0_c + du_CFL)
    FAmt_L = 0.0 ; FAmt_R = 0.0 ; FAmt_0 = 0.0
    uhtot_L = 0.0 ; uhtot_R = 0.0

  do k=1,nz

      if (do_c) then
        visc_rem_lim = max(visc_rem(II,jj,k), min_visc_rem*visc_rem_max_c)
        if (visc_rem_lim > 0.0) then ! This is almost always true for ocean points.
          if (u(I,j,k) + duR*visc_rem_lim > -du_CFL*visc_rem(II,jj,k)) &
            duR = -(u(I,j,k) + du_CFL*visc_rem(II,jj,k)) / visc_rem_lim
          if (u(I,j,k) + duL*visc_rem_lim < du_CFL*visc_rem(II,jj,k)) &
            duL = -(u(I,j,k) - du_CFL*visc_rem(II,jj,k)) / visc_rem_lim
        endif
      endif
  enddo

      if (do_c) then
  do k=1,nz
        u_L = u(I,j,k) + duL * visc_rem(II,jj,k)
        u_R = u(I,j,k) + duR * visc_rem(II,jj,k)
        u_0 = u(I,j,k) + du0_c * visc_rem(II,jj,k)

        call flux_elem(u_0, h_in(i,j,k), h_in(i+1,j,k), &
            h_W(i,j,k), h_W(i+1,j,k), h_E(i,j,k), h_E(i+1,j,k), &
            uh_0, duhdu_0, visc_rem(II,jj,k), dy_Cu, &
            IareaT_W, IareaT_E, &
            IdxT_W, IdxT_E, dt, vol_CFL, &
            por_face_areaU(I,j,k), h_marg_min)

        call flux_elem(u_L, h_in(i,j,k), h_in(i+1,j,k), &
            h_W(i,j,k), h_W(i+1,j,k), h_E(i,j,k), h_E(i+1,j,k), &
            uh_L, duhdu_L, visc_rem(II,jj,k), dy_Cu, &
            IareaT_W, IareaT_E, &
            IdxT_W, IdxT_E, dt, vol_CFL, &
            por_face_areaU(I,j,k), h_marg_min)

        call flux_elem(u_R, h_in(i,j,k), h_in(i+1,j,k), &
            h_W(i,j,k), h_W(i+1,j,k), h_E(i,j,k), h_E(i+1,j,k), &
            uh_R, duhdu_R, visc_rem(II,jj,k), dy_Cu, &
            IareaT_W, IareaT_E, IdxT_W, IdxT_E, dt, &
            vol_CFL, por_face_areaU(I,j,k), h_marg_min)

        ! XXX: Add zero to force rounding and prevent FMA optimization.
        FAmt_0 = FAmt_0 + (duhdu_0 + 0.)
        FAmt_L = FAmt_L + (duhdu_L + 0.)
        FAmt_R = FAmt_R + (duhdu_R + 0.)
        uhtot_L = uhtot_L + uh_L
        uhtot_R = uhtot_R + uh_R
  enddo

      FA_0 = FAmt_0
      FA_avg = FAmt_0

      if (duL - du0_c /= 0.) &
        FA_avg = uhtot_L / (duL - du0_c)

      if (FA_avg > max(FA_0, FAmt_L)) then
        FA_avg = max(FA_0, FAmt_L)
      elseif (FA_avg < min(FA_0, FAmt_L)) then
        FA_0 = FA_avg
      endif

      BT_cont%FA_u_W0(I,j) = FA_0
      BT_cont%FA_u_WW(I,j) = FAmt_L

      if (abs(FA_0-FAmt_L) <= 1e-12 * FA_0) then
        BT_cont%uBT_WW(I,j) = 0.
      else
        BT_cont%uBT_WW(I,j) = (1.5 * (duL - du0_c)) &
            * ((FAmt_L - FA_avg) / (FAmt_L - FA_0))
      endif

      FA_0 = FAmt_0
      FA_avg = FAmt_0

      if (duR - du0_c /= 0.) &
        FA_avg = uhtot_R / (duR - du0_c)

      if (FA_avg > max(FA_0, FAmt_R)) then
        FA_avg = max(FA_0, FAmt_R)
      elseif (FA_avg < min(FA_0, FAmt_R)) then
        FA_0 = FA_avg
      endif

      BT_cont%FA_u_E0(I,j) = FA_0 ; BT_cont%FA_u_EE(I,j) = FAmt_R

      if (abs(FAmt_R - FA_0) <= 1e-12 * FA_0) then
        BT_cont%uBT_EE(I,j) = 0.
      else
        BT_cont%uBT_EE(I,j) = (1.5 * (duR - du0_c)) &
            * ((FAmt_R - FA_avg) / (FAmt_R - FA_0))
      endif
    else
      BT_cont%FA_u_W0(I,j) = 0.
      BT_cont%FA_u_WW(I,j) = 0.
      BT_cont%FA_u_E0(I,j) = 0.
      BT_cont%FA_u_EE(I,j) = 0.
      BT_cont%uBT_WW(I,j) = 0.
      BT_cont%uBT_EE(I,j) = 0.
    endif
  enddo ; enddo

  !$omp target exit data map(release: du0, zeros, uh_tmp)
end procedure set_zonal_BT_cont

module procedure meridional_mass_flux
  ! Local variables
  ! The block-sized work arrays below are dimensioned (nii,nJJB[,nz]).  In this meridional
  ! routine their first (nii) block index is at the tracer-column (i) points, while the second
  ! (nJJB) block index is staggered at the meridional velocity (v, C-grid J) points.
  real, dimension(nii,nJJB) :: &
    vhbt_b
  real, dimension(nii,nJJB,SZK_(GV)) :: &
    vh_b
  real, dimension(nii,nJJB,SZK_(GV)) :: &
    dvhdv         ! Partial derivative of vh with v, at v points [H L ~> m2 or kg m-1].
  real, dimension(nii,nJJB) :: &
    dv, &         ! Corrective barotropic change in the velocity to give vhbt [L T-1 ~> m s-1].
    dv_min_CFL, & ! Lower limit on dv correction to avoid CFL violations [L T-1 ~> m s-1]
    dv_max_CFL, & ! Upper limit on dv correction to avoid CFL violations [L T-1 ~> m s-1]
    dvhdv_tot_0, & ! Summed partial derivative of vh with v [H L ~> m2 or kg m-1].
    vh_tot_0, &   ! Summed transport with no barotropic correction [H L2 T-1 ~> m3 s-1 or kg s-1].
    visc_rem_max  ! The column maximum of visc_rem [nondim]
  logical, dimension(nii,nJJB) :: do_ij
  real, dimension(nii,nJJB) :: FAvi  ! A list of sums of meridional face areas
                                                   ! [H L ~> m2 or kg m-1].
  real :: FA_v    ! A sum of meridional face areas [H L ~> m2 or kg m-1].
  real, dimension(nii,nJJB,SZK_(GV)) :: &
    visc_rem      ! A 2-D copy of visc_rem_v or an array of 1's [nondim]
  real :: I_vrm   ! 1.0 / visc_rem_max [nondim]
  real :: CFL_dt  ! The maximum CFL ratio of the adjusted velocities divided by
                  ! the time step [T-1 ~> s-1].
  real :: I_dt    ! 1.0 / dt [T-1 ~> s-1].
  real :: dv_lim  ! The velocity change that give a relative CFL of 1 [L T-1 ~> m s-1].
  real :: dy_N, dy_S ! Effective y-grid spacings to the north and south [L ~> m].
  type(cont_loop_bounds_type) :: LB
  integer :: i, j, k, ish, ieh, jsh, jeh, n, nz
  integer :: l_seg ! The OBC segment number
  logical :: use_visc_rem, set_BT_cont
  logical :: local_specified_BC, local_Flather_OBC ! OBC-related logicals
  logical :: local_open_BC, any_simple_OBC          ! OBC-related logicals
  logical :: simple_OBC_pt(nii,nJJB) ! Indicates points with specified transport OBCs
  type(OBC_segment_type), pointer :: segment => NULL()
  integer :: isb, ieb, JsbB, JebB, ii, JJ, iie, JJe
  integer :: nteams
  logical :: set_h_v

  use_visc_rem = present(visc_rem_v)

  set_BT_cont = .false. ; if (present(BT_cont)) set_BT_cont = (associated(BT_cont))
  set_h_v = .false. ; if (set_BT_cont) set_h_v = allocated(BT_cont%h_v)

  !   The arrays that are used on the device are copied in before the clock is started.  The output
  ! arrays that are only partly set here are copied in too, so that the copies back to the host
  ! leave the rest of them unchanged.
  !$omp target enter data map(to: G)
  !$omp target enter data map(to: G%dx_Cv, G%dyCv, G%IareaT, G%areaT, G%IdyT, G%dyT, G%mask2dCv)
  !$omp target enter data map(to: CS)
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
  !$omp                              vh_tot_0, visc_rem_max, do_ij, vh_b, vhbt_b)

  local_specified_BC = .false. ; local_Flather_OBC = .false. ; local_open_BC = .false.
  if (associated(OBC)) then ; if (OBC%OBC_pe) then
    local_specified_BC = OBC%specified_v_BCs_exist_globally
    local_Flather_OBC = OBC%Flather_v_BCs_exist_globally
    local_open_BC = OBC%open_v_BCs_exist_globally
  endif ; endif

  if (present(dv_cor)) then
    !$omp target teams distribute parallel do collapse(2)
    do J=G%jsdB,G%jedB ; do i=G%isd,G%ied
      dv_cor(i,J) = 0.0
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

  if (.not.use_visc_rem) then
    !$omp target teams distribute parallel do collapse(3)
    do k=1,nz ; do JJ=1,nJJB ; do ii=1,nii
      visc_rem(ii,JJ,k) = 1.0
    enddo ; enddo ; enddo
  endif

  do JsbB = jsh-1,jeh,nJJB ; do isb = ish,ieh,nii
    JebB = min(JsbB + nJJB - 1, jeh)
    ieb = min(isb + nii - 1, ieh)
    JJe = JebB - JsbB + 1
    iie = ieb - isb + 1

    !$omp target teams distribute parallel do collapse(2)
    do JJ=1,JJe ; do ii=1,iie
      do_ij(ii,JJ) = .true.
    enddo ; enddo

    ! visc_rem_v is copied in a kernel of its own, as a kernel that names an absent argument fails.
    if (use_visc_rem) then
      !$omp target teams distribute parallel do collapse(3) private(i,J)
      do k=1,nz ; do JJ=1,JJe ; do ii=1,iie
        i = isb + ii - 1
        J = JsbB + JJ - 1
        visc_rem(ii,JJ,k) = visc_rem_v(i,J,k)
      enddo ; enddo ; enddo
    endif

    ! This sets vh and dvhdv.
    !$omp target teams distribute parallel do collapse(3) private(i,J,l_seg)
    do k=1,nz ; do JJ=1,JJe ; do ii=1,iie
      i = isb + ii - 1
      J = JsbB + JJ - 1

      call flux_elem(v(i,J,k), h_in(i,J,k), h_in(i,J+1,k), &
          h_S(i,J,k), h_S(i,J+1,k), h_N(i,J,k), h_N(i,J+1,k), &
          vh_b(ii,JJ,k), dvhdv(ii,JJ,k), visc_rem(ii,JJ,k), G%dx_Cv(i,J), &
          G%IareaT(i,J), G%IareaT(i,J+1), G%IdyT(i,J), G%IdyT(i,J+1), dt, &
          CS%vol_CFL, por_face_areaV(i,J,k), CS%h_marg_min)

      if (local_open_BC) &
        call flux_elem_OBC(v(i,J,k), h_in(i,J,k), h_in(i,J+1,k), &
            vh_b(ii,JJ,k), dvhdv(ii,JJ,k), visc_rem(ii,JJ,k), &
            por_face_areaV(i,J,k), G%dx_Cv(i,J), OBC, OBC%segnum_v(i,J), &
            CS%h_marg_min)

      if (local_specified_BC) then
        if (OBC%segnum_v(i,J) /= 0) then
          l_seg = abs(OBC%segnum_v(i,J))
          if (OBC%segment(l_seg)%specified) &
            vh_b(ii,JJ,k) = OBC%segment(l_seg)%normal_trans(i,J,k)
        endif
      endif
    enddo ; enddo ; enddo

    if (present(vhbt) .or. set_BT_cont) then
      if (use_visc_rem .and. CS%use_visc_rem_max) then
        !$omp target teams distribute parallel do collapse(2)
        do JJ=1,JJe ; do ii=1,iie
          visc_rem_max(ii,JJ) = 0.0
          do k=1,nz
            visc_rem_max(ii,JJ) = max(visc_rem_max(ii,JJ), visc_rem(ii,JJ,k))
          enddo
        enddo ; enddo
      else
        !$omp target teams distribute parallel do collapse(2)
        do JJ=1,JJe ; do ii=1,iie
          visc_rem_max(ii,JJ) = 1.0
        enddo ; enddo
      endif

      !   Set limits on dv that will keep the CFL number between -1 and 1.
      ! This should be adequate to keep the root bracketed in all cases.
      !$omp target teams distribute parallel do collapse(2) private(i,J,I_vrm,dy_S,dy_N)
      do JJ=1,JJe ; do ii=1,iie
        i = isb + ii - 1
        J = JsbB + JJ - 1

        I_vrm = 0.0
        if (visc_rem_max(ii,JJ) > 0.0) I_vrm = 1.0 / visc_rem_max(ii,JJ)
        if (CS%vol_CFL) then
          dy_S = ratio_max(G%areaT(i,j), G%dx_Cv(i,J), 1000.0*G%dyT(i,j))
          dy_N = ratio_max(G%areaT(i,j+1), G%dx_Cv(i,J), 1000.0*G%dyT(i,j+1))
        else ; dy_S = G%dyT(i,j) ; dy_N = G%dyT(i,j+1) ; endif
        dv_max_CFL(ii,JJ) = 2.0 * (CFL_dt * dy_S) * I_vrm
        dv_min_CFL(ii,JJ) = -2.0 * (CFL_dt * dy_N) * I_vrm
        vh_tot_0(ii,JJ) = 0.0 ; dvhdv_tot_0(ii,JJ) = 0.0

        do k=1,nz
          dvhdv_tot_0(ii,JJ) = dvhdv_tot_0(ii,JJ) + dvhdv(ii,JJ,k)
          vh_tot_0(ii,JJ) = vh_tot_0(ii,JJ) + vh_b(ii,JJ,k)
        enddo
      enddo ; enddo

      if (use_visc_rem) then
        if (CS%aggress_adjust) then
          !$omp target teams distribute parallel do collapse(2) private(i,J,dy_S,dy_N,dv_lim)
          do JJ=1,JJe ; do ii=1,iie
            i = isb + ii - 1
            J = JsbB + JJ - 1
            do k=1,nz
              if (CS%vol_CFL) then
                dy_S = ratio_max(G%areaT(i,j), G%dx_Cv(i,J), 1000.0*G%dyT(i,j))
                dy_N = ratio_max(G%areaT(i,j+1), G%dx_Cv(i,J), 1000.0*G%dyT(i,j+1))
              else ; dy_S = G%dyT(i,j) ; dy_N = G%dyT(i,j+1) ; endif
              dv_lim = 0.499*((dy_S*I_dt - v(i,J,k)) + MIN(0.0,v(i,J-1,k)))
              if (dv_max_CFL(ii,JJ) * visc_rem(ii,JJ,k) > dv_lim) &
                dv_max_CFL(ii,JJ) = dv_lim / visc_rem(ii,JJ,k)

              dv_lim = 0.499*((-dy_N*CFL_dt - v(i,J,k)) + MAX(0.0,v(i,J+1,k)))
              if (dv_min_CFL(ii,JJ) * visc_rem(ii,JJ,k) < dv_lim) &
                dv_min_CFL(ii,JJ) = dv_lim / visc_rem(ii,JJ,k)
            enddo
          enddo ; enddo
        else
          !$omp target teams distribute parallel do collapse(2) private(i,J,dy_S,dy_N)
          do JJ=1,JJe ; do ii=1,iie
            i = isb + ii - 1
            J = JsbB + JJ - 1
            do k=1,nz
              if (CS%vol_CFL) then
                dy_S = ratio_max(G%areaT(i,j), G%dx_Cv(i,J), 1000.0*G%dyT(i,j))
                dy_N = ratio_max(G%areaT(i,j+1), G%dx_Cv(i,J), 1000.0*G%dyT(i,j+1))
              else ; dy_S = G%dyT(i,j) ; dy_N = G%dyT(i,j+1) ; endif
              if (dv_max_CFL(ii,JJ) * visc_rem(ii,JJ,k) > &
                  dy_S*CFL_dt - v(i,J,k)*G%mask2dCv(i,J)) &
                dv_max_CFL(ii,JJ) = (dy_S*CFL_dt - v(i,J,k)) / visc_rem(ii,JJ,k)
              if (dv_min_CFL(ii,JJ) * visc_rem(ii,JJ,k) < &
                  -dy_N*CFL_dt - v(i,J,k)*G%mask2dCv(i,J)) &
                dv_min_CFL(ii,JJ) = -(dy_N*CFL_dt + v(i,J,k)) / visc_rem(ii,JJ,k)
            enddo
          enddo ; enddo
        endif
      else
        if (CS%aggress_adjust) then
          !$omp target teams distribute parallel do collapse(2) private(i,J,dy_S,dy_N)
          do JJ=1,JJe ; do ii=1,iie
            i = isb + ii - 1
            J = JsbB + JJ - 1
            do k=1,nz
              if (CS%vol_CFL) then
                dy_S = ratio_max(G%areaT(i,j), G%dx_Cv(i,J), 1000.0*G%dyT(i,j))
                dy_N = ratio_max(G%areaT(i,j+1), G%dx_Cv(i,J), 1000.0*G%dyT(i,j+1))
              else ; dy_S = G%dyT(i,j) ; dy_N = G%dyT(i,j+1) ; endif
              dv_max_CFL(ii,JJ) = min(dv_max_CFL(ii,JJ), 0.499 * &
                          ((dy_S*I_dt - v(i,J,k)) + MIN(0.0,v(i,J-1,k))) )
              dv_min_CFL(ii,JJ) = max(dv_min_CFL(ii,JJ), 0.499 * &
                          ((-dy_N*I_dt - v(i,J,k)) + MAX(0.0,v(i,J+1,k))) )
            enddo
          enddo ; enddo
        else
          !$omp target teams distribute parallel do collapse(2) private(i,J,dy_S,dy_N)
          do JJ=1,JJe ; do ii=1,iie
            i = isb + ii - 1
            J = JsbB + JJ - 1
            do k=1,nz
              if (CS%vol_CFL) then
                dy_S = ratio_max(G%areaT(i,j), G%dx_Cv(i,J), 1000.0*G%dyT(i,j))
                dy_N = ratio_max(G%areaT(i,j+1), G%dx_Cv(i,J), 1000.0*G%dyT(i,j+1))
              else ; dy_S = G%dyT(i,j) ; dy_N = G%dyT(i,j+1) ; endif
              dv_max_CFL(ii,JJ) = min(dv_max_CFL(ii,JJ), dy_S*CFL_dt - v(i,J,k))
              dv_min_CFL(ii,JJ) = max(dv_min_CFL(ii,JJ), -(dy_N*CFL_dt + v(i,J,k)))
            enddo
          enddo ; enddo
        endif
      endif
      !$omp target teams distribute parallel do collapse(2)
      do JJ=1,JJe ; do ii=1,iie
        dv_max_CFL(ii,JJ) = max(dv_max_CFL(ii,JJ),0.0)
        dv_min_CFL(ii,JJ) = min(dv_min_CFL(ii,JJ),0.0)
      enddo ; enddo

      any_simple_OBC = .false.
      if (present(vhbt) .or. set_BT_cont) then
        if (local_specified_BC .or. local_Flather_OBC) then
          do JJ=1,JJe ; do ii=1,iie
            i = isb + ii - 1
            J = JsbB + JJ - 1

            l_seg = abs(OBC%segnum_v(i,J))

            ! Avoid reconciling barotropic/baroclinic transports if transport is specified
            simple_OBC_pt(ii,JJ) = .false.
            if (l_seg /= 0) simple_OBC_pt(ii,JJ) = OBC%segment(l_seg)%specified
            do_ij(ii,JJ) = .not.simple_OBC_pt(ii,JJ)
            any_simple_OBC = any_simple_OBC .or. simple_OBC_pt(ii,JJ)
          enddo ; enddo
        else
          !$omp target teams distribute parallel do collapse(2)
          do JJ=1,JJe ; do ii=1,iie
            do_ij(ii,JJ) = .true.
          enddo ; enddo
        endif
      endif

      if (present(vhbt)) then
        !$omp target teams distribute parallel do collapse(2) private(i,J)
        do JJ=1,JJe ; do ii=1,iie
          i = isb + ii - 1
          J = JsbB + JJ - 1
          vhbt_b(ii,JJ) = vhbt(i,J)
        enddo ; enddo
        ! Find dv and vh.
        call meridional_flux_adjust(v, h_in, h_S, h_N, vhbt_b, vh_tot_0, dvhdv_tot_0, dv, &
                                    dv_max_CFL, dv_min_CFL, dt, G, GV, US, CS, visc_rem, &
                                    isb, ieb, JsbB, JebB, do_ij, por_face_areaV, &
                                    nii, nJJB, vh_b, OBC)

        if (present(v_cor)) then
          !$omp target teams distribute parallel do collapse(3) private(i,J)
          do k=1,nz ; do JJ=1,JJe ; do ii=1,iie
            i = isb + ii - 1
            J = JsbB + JJ - 1
            v_cor(i,J,k) = v(i,J,k) + dv(ii,JJ) * visc_rem(ii,JJ,k)
          enddo ; enddo ; enddo

          if (any_simple_OBC) then
            do k=1,nz ; do JJ=1,JJe ; do ii=1,iie
              if (simple_OBC_pt(ii,JJ)) then
                i = isb + ii - 1
                J = JsbB + JJ - 1
                v_cor(i,J,k) = OBC%segment(abs(OBC%segnum_v(i,J)))%normal_vel(i,J,k)
              endif
            enddo ; enddo ; enddo
          endif
        endif

        if (present(dv_cor)) then
          !$omp target teams distribute parallel do collapse(2) private(i,J)
          do JJ=1,JJe ; do ii=1,iie
            i = isb + ii - 1
            J = JsbB + JJ - 1
            dv_cor(i,J) = dv(ii,JJ)
          enddo ; enddo
        endif
      endif

      if (set_BT_cont) then
        call set_merid_BT_cont(v, h_in, h_S, h_N, BT_cont, vh_tot_0, dvhdv_tot_0, &
                               dv_max_CFL, dv_min_CFL, dt, G, GV, US, CS, visc_rem, &
                               visc_rem_max, isb, ieb, JsbB, JebB, do_ij, &
                               por_face_areaV, nii, nJJB)
        if (any_simple_OBC) then
          do JJ=1,JJe ; do ii=1,iie
            if (simple_OBC_pt(ii,JJ)) then
              i = isb + ii - 1
              J = JsbB + JJ - 1
              FAvi(ii,JJ) = GV%H_subroundoff*G%dx_Cv(i,J)
            endif
          enddo ; enddo

          ! NOTE: simple_OBC_pt should prevent access to segment OBC_NONE
          do k=1,nz
            do JJ=1,JJe ; do ii=1,iie
              if (simple_OBC_pt(ii,JJ)) then
                i = isb + ii - 1
                J = JsbB + JJ - 1

                l_seg = abs(OBC%segnum_v(i,J))
                if ((abs(OBC%segment(l_seg)%normal_vel(i,J,k)) > 0.0) .and. &
                    (OBC%segment(l_seg)%specified)) &
                  FAvi(ii,JJ) = FAvi(ii,JJ) + &
                    OBC%segment(l_seg)%normal_trans(i,J,k) / OBC%segment(l_seg)%normal_vel(i,J,k)
              endif
            enddo ; enddo
          enddo

          do JJ=1,JJe ; do ii=1,iie
            if (simple_OBC_pt(ii,JJ)) then
              i = isb + ii - 1
              J = JsbB + JJ - 1
              BT_cont%FA_v_S0(i,J) = FAvi(ii,JJ) ; BT_cont%FA_v_N0(i,J) = FAvi(ii,JJ)
              BT_cont%FA_v_SS(i,J) = FAvi(ii,JJ) ; BT_cont%FA_v_NN(i,J) = FAvi(ii,JJ)
              BT_cont%vBT_SS(i,J) = 0.0 ; BT_cont%vBT_NN(i,J) = 0.0
            endif
          enddo ; enddo
        endif
      endif ! set_BT_cont

    endif ! present(vhbt) or set_BT_cont

    !$omp target teams distribute parallel do collapse(3) private(i,J)
    do k=1,nz ; do JJ=1,JJe ; do ii=1,iie
      i = isb + ii - 1
      J = JsbB + JJ - 1
      vh(i,J,k) = vh_b(ii,JJ,k)
    enddo ; enddo ; enddo

  enddo ; enddo ! ij block loops

  if (local_open_BC .and. set_BT_cont) then
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
  endif

  if (set_BT_cont) then ; if (allocated(BT_cont%h_v)) then
    if (present(v_cor)) then
      call meridional_flux_thickness(v_cor, h_in, h_S, h_N, BT_cont%h_v, dt, G, GV, US, LB, &
                                    CS%vol_CFL, CS%marginal_faces, OBC, por_face_areaV, visc_rem_v)
    else
      call meridional_flux_thickness(v, h_in, h_S, h_N, BT_cont%h_v, dt, G, GV, US, LB, &
                                    CS%vol_CFL, CS%marginal_faces, OBC, por_face_areaV, visc_rem_v)
    endif
  endif ; endif

  !$omp target exit data map(release: dvhdv, visc_rem, dv, dv_min_CFL, dv_max_CFL, dvhdv_tot_0, &
  !$omp                              vh_tot_0, visc_rem_max, do_ij, vh_b, vhbt_b)

  call cpu_clock_end(id_clock_correct)

  !   The outputs are copied back to the host, and the inputs released, after the clock is stopped.
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
  !$omp target exit data map(release: CS)
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
  real :: dh
  ish = LB%ish ; ieh = LB%ieh ; jsh = LB%jsh ; jeh = LB%jeh ; nz = GV%ke

  !$omp target teams distribute parallel do collapse(3) private(CFL, curv_3, dh)
  do k=1,nz ; do J=jsh-1,jeh ; do i=ish,ieh
    if (v(i,J,k) > 0.0) then
      if (vol_CFL) then ; CFL = (v(i,J,k) * dt) * (G%dx_Cv(i,J) * G%IareaT(i,j))
      else ; CFL = v(i,J,k) * dt * G%IdyT(i,j) ; endif
      curv_3 = (h_S(i,j,k) + h_N(i,j,k)) - 2.0*h(i,j,k)
      dh = h_S(i,J,k) - h_N(i,J,k)
      if (marginal) then
        h_v(i,J,k) = h_N(i,j,k) + CFL * (dh + 3.0*curv_3*(CFL - 1.0))
      else
        h_v(i,J,k) = h_N(i,j,k) + CFL * (0.5*dh + curv_3*(CFL - 1.5))
      endif
    elseif (v(i,J,k) < 0.0) then
      if (vol_CFL) then ; CFL = (-v(i,J,k)*dt) * (G%dx_Cv(i,J) * G%IareaT(i,j+1))
      else ; CFL = -v(i,J,k) * dt * G%IdyT(i,j+1) ; endif
      curv_3 = (h_S(i,j+1,k) + h_N(i,j+1,k)) - 2.0*h(i,j+1,k)
      dh = h_N(i,j+1,k)-h_S(i,j+1,k)
      if (marginal) then
        h_v(i,J,k) = h_S(i,j+1,k) + CFL * (dh + 3.0*curv_3*(CFL - 1.0))
      else
        h_v(i,J,k) = h_S(i,j+1,k) + CFL * (0.5*dh + curv_3*(CFL - 1.5))
      endif
    else
      !   The choice to use the arithmetic mean here is somewhat arbitrarily, but
      ! it should be noted that h_S(i+1,j,k) and h_N(i,j,k) are usually the same.
      h_v(i,J,k) = 0.5 * (h_S(i,j+1,k) + h_N(i,j,k))
 !    h_marg = (2.0 * h_S(i,j+1,k) * h_N(i,j,k)) / &
 !             (h_S(i,j+1,k) + h_N(i,j,k) + GV%H_subroundoff)
    endif

    if (present(visc_rem_v)) then
      ! Scale back the thickness to account for the effects of viscosity and the fractional open
      ! thickness to give an appropriate non-normalized weight for each layer in determining the
      ! barotropic acceleration.
      h_v(i,J,k) = h_v(i,J,k) * (visc_rem_v(i,J,k) * por_face_areaV(i,J,k))
    else
      h_v(i,J,k) = h_v(i,J,k) * por_face_areaV(i,J,k)
    endif
  enddo ; enddo ; enddo

  local_open_BC = .false.
  if (associated(OBC)) local_open_BC = OBC%open_v_BCs_exist_globally

  if (local_open_BC) then
    do n = 1, OBC%number_of_segments
      if (OBC%segment(n)%open .and. OBC%segment(n)%is_N_or_S) then
        J = OBC%segment(n)%HI%JsdB
        if (OBC%segment(n)%direction == OBC_DIRECTION_N) then
          if (present(visc_rem_v)) then
            do k=1,nz ; do i=OBC%segment(n)%HI%isd,OBC%segment(n)%HI%ied
              h_v(i,J,k) = h(i,J,k) * (visc_rem_v(i,J,k) * por_face_areaV(i,J,k))
            enddo ; enddo
          else
            do k=1,nz ; do i=OBC%segment(n)%HI%isd,OBC%segment(n)%HI%ied
              h_v(i,J,k) = h(i,J,k) * por_face_areaV(i,J,k)
            enddo ; enddo
          endif
        else
          if (present(visc_rem_v)) then
            do k=1,nz ; do i=OBC%segment(n)%HI%isd,OBC%segment(n)%HI%ied
              h_v(i,J,k) = h(i,J+1,k) * (visc_rem_v(i,J,k) * por_face_areaV(i,J,k))
            enddo ; enddo
          else
            do k=1,nz ; do i=OBC%segment(n)%HI%isd,OBC%segment(n)%HI%ied
              h_v(i,J,k) = h(i,J+1,k) * por_face_areaV(i,J,k)
            enddo ; enddo
          endif
        endif
      endif
    enddo
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

module procedure meridional_flux_adjust
  ! Local variables
  ! The per-column iteration state is held in private scalars rather than in (nii,nJJB) arrays,
  ! so that it stays in registers instead of being loaded and stored at every level.
  real :: &
    dvhdv      ! Partial derivative of vh with v [H L ~> m2 or kg m-1].
  real :: dv_c       ! The barotropic velocity correction for this column [L T-1 ~> m s-1].
  real :: vh_err     ! Difference between vhbt and the summed vh [H L2 T-1 ~> m3 s-1 or kg s-1].
  real :: vh_err_best ! The smallest value of vh_err found so far [H L2 T-1 ~> m3 s-1 or kg s-1].
  real :: dvhdv_tot  ! Summed partial derivative of vh with u [H L ~> m2 or kg m-1].
  real :: dv_min     ! Lower limit on dv correction based on CFL limits and previous iterations [L T-1 ~> m s-1]
  real :: dv_max     ! Upper limit on dv correction based on CFL limits and previous iterations [L T-1 ~> m s-1]
  real :: dx_Cv      ! The unblocked length of this v-face [L ~> m].
  real :: IareaT_S, IareaT_N ! The inverse areas of the tracer cells south and north of this face [L-2 ~> m-2].
  real :: IdyT_S, IdyT_N     ! The inverse meridional grid spacings south and north of this face [L-1 ~> m-1].
  real :: tol_eta_CS ! A copy of CS%tol_eta for use in the kernel [H ~> m or kg m-2].
  real :: h_marg_min ! A copy of CS%h_marg_min for use in the kernel [H ~> m or kg m-2].
  logical :: vol_CFL ! A copy of CS%vol_CFL for use in the kernel.
  logical :: better_iter ! A copy of CS%better_iter for use in the kernel.
  logical :: do_c ! If true, this column is still being iterated.
  real :: &
    vh, &      ! The volume flux through a meridional face [H L2 T-1 ~> m3 s-1 or kg s-1].
    v_new   ! The velocity with the correction added [L T-1 ~> m s-1].
  real :: dv_prev ! The previous value of dv [L T-1 ~> m s-1].
  real :: ddv     ! The change in dv from the previous iteration [L T-1 ~> m s-1].
  real :: tol_eta ! The tolerance for the current iteration [H ~> m or kg m-2].
  real :: tol_vel ! The tolerance for velocity in the current iteration [L T-1 ~> m s-1].
  integer :: i, j, k, nz, itt, ii, JJ, iie, JJe
  logical :: local_open_BC
  type(ocean_OBC_type), pointer :: OBC_loc ! OBC, or null if it is absent, for use in kernels,
                                           ! as a kernel that names an absent argument fails.
  integer, parameter :: max_itts = 20

  local_open_BC = .false.
  OBC_loc => null()
  if (present(OBC)) then
    OBC_loc => OBC
    if (associated(OBC)) local_open_BC = OBC%open_v_BCs_exist_globally
  endif

  nz = GV%ke
  iie = ieb - isb + 1
  JJe = JebB - JsbB + 1
  tol_vel = CS%tol_vel ; tol_eta_CS = CS%tol_eta ; better_iter = CS%better_iter
  h_marg_min = CS%h_marg_min ; vol_CFL = CS%vol_CFL

  !$omp target teams distribute parallel do collapse(2) &
  !$omp   private(i, J, tol_eta, ddv, dv_prev, v_new, vh, dvhdv, dv_c, vh_err, vh_err_best, &
  !$omp           dvhdv_tot, dv_min, dv_max, do_c, dx_Cv, IareaT_S, IareaT_N, IdyT_S, IdyT_N)
  do JJ=1,JJe ; do ii=1,iie
    i = isb + ii - 1
    J = JsbB + JJ - 1
    dx_Cv = G%dx_Cv(i,J)
    IareaT_S = G%IareaT(i,J) ; IareaT_N = G%IareaT(i,J+1)
    IdyT_S = G%IdyT(i,J) ; IdyT_N = G%IdyT(i,J+1)
    dv_c = 0.0 ; do_c = do_ij_in(ii,JJ)
    dv_max = dv_max_CFL(ii,JJ) ; dv_min = dv_min_CFL(ii,JJ)
    vh_err = vh_tot_0(ii,JJ) - vhbt(ii,JJ) ; dvhdv_tot = dvhdv_tot_0(ii,JJ)
    vh_err_best = abs(vh_err)

    do itt=1,max_itts
      select case (itt)
        case (:1) ; tol_eta = 1e-6 * tol_eta_CS
        case (2)  ; tol_eta = 1e-4 * tol_eta_CS
        case (3)  ; tol_eta = 1e-2 * tol_eta_CS
        case default ; tol_eta = tol_eta_CS
      end select

      if (vh_err > 0.0) then ; dv_max = dv_c
      elseif (vh_err < 0.0) then ; dv_min = dv_c
      else ; do_c = .false. ; endif

      if (do_c) then
        if ((dt * min(IareaT_S,IareaT_N)*abs(vh_err) > tol_eta) .or. &
            (better_iter .and. ((abs(vh_err) > tol_vel * dvhdv_tot) .or. &
                                   (abs(vh_err) > vh_err_best)) )) then
          !   Use Newton's method, provided it stays bounded.  Otherwise bisect
          ! the value with the appropriate bound.
          ddv = -vh_err / dvhdv_tot
          dv_prev = dv_c
          dv_c = dv_c + ddv
          if (abs(ddv) < 1.0e-15*abs(dv_c)) then
            do_c = .false. ! ddv is small enough to quit.
          elseif (ddv > 0.0) then
            if (dv_c >= dv_max) then
              dv_c = 0.5*(dv_prev + dv_max)
              if (dv_max - dv_prev < 1.0e-15*abs(dv_c)) do_c = .false.
            endif
          else ! dvv(i,J) < 0.0
            if (dv_c <= dv_min) then
              dv_c = 0.5*(dv_prev + dv_min)
              if (dv_prev - dv_min < 1.0e-15*abs(dv_c)) do_c = .false.
            endif
          endif
        else
          do_c = .false.
        endif
      endif

      if (.not.do_c) exit

      vh_err = -vhbt(ii,JJ) ; dvhdv_tot = 0.0

      do k=1,nz
        v_new = v(i,J,k) + dv_c * visc_rem(ii,JJ,k)
        call flux_elem(v_new, h_in(i,J,k), h_in(i,J+1,k), &
            h_S(i,J,k), h_S(i,J+1,k), h_N(i,J,k), h_N(i,J+1,k), &
            vh, dvhdv, visc_rem(ii,JJ,k), dx_Cv, &
            IareaT_S, IareaT_N, IdyT_S, IdyT_N, dt, &
            vol_CFL, por_face_areaV(i,J,k), h_marg_min)

        if (local_open_BC) &
          call flux_elem_OBC(v_new, h_in(i,J,k), h_in(i,J+1,k), &
              vh, dvhdv, visc_rem(ii,JJ,k), &
              por_face_areaV(i,J,k), dx_Cv, OBC_loc, OBC_loc%segnum_v(i,J), &
              h_marg_min)

        vh_3d(ii,JJ,k) = vh
        vh_err = vh_err + vh
        dvhdv_tot = dvhdv_tot + dvhdv
      enddo

      vh_err_best = min(vh_err_best, abs(vh_err))
    enddo ! itt-loop
    dv(ii,JJ) = dv_c
  enddo ; enddo
  ! If there are any faces which have not converged to within the tolerance,
  ! so-be-it, or else use a final upwind correction?
  ! This never seems to happen with 20 iterations as max_itt.

end procedure meridional_flux_adjust

module procedure set_merid_BT_cont
  ! Local variables
  ! The block-sized work arrays below are dimensioned (nii,nJJB[,nz]).  In this meridional
  ! routine their first (nii) block index is at the tracer-column (i) points, while the second
  ! (nJJB) block index is staggered at the meridional velocity (v, C-grid J) points.
  real, dimension(nii,nJJB) :: &
    dv0           ! The barotropic velocity increment that gives 0 transport [L T-1 ~> m s-1].
  ! The per-column values below are private scalars in the kernel, so that they stay in registers.
  real :: &
    dv0_c, &      ! A copy of dv0 for this column [L T-1 ~> m s-1].
    dvL, dvR, &   ! The barotropic velocity increments that give the southerly
                  ! (dvL) and northerly (dvR) test velocities [L T-1 ~> m s-1].
    dv_CFL, &     ! The velocity increment that corresponds to CFL_min [L T-1 ~> m s-1].
    FAmt_L, FAmt_R, & ! The summed effective marginal face areas for the 3
    FAmt_0, &     ! test velocities [H L ~> m2 or kg m-1].
    vhtot_L, &    ! The summed transport with the southerly (vhtot_L) and
    vhtot_R       ! and northerly (vhtot_R) test velocities [H L2 T-1 ~> m3 s-1 or kg s-1].
  real, dimension(nii,nJJB) :: &
    zeros      ! An array of full of 0 transports [H L2 T-1 ~> m3 s-1 or kg s-1]
  real, dimension(nii,nJJB,SZK_(GV)) :: &
    vh_tmp     ! A temporary array for the unused layer transports [H L2 T-1 ~> m3 s-1 or kg s-1]
  real :: &
    vh_L, vh_R, & ! The layer transports with the southerly (_L), northerly (_R)
    vh_0, &       ! and zero-barotropic (_0) test velocities [H L2 T-1 ~> m3 s-1 or kg s-1].
    dvhdv_L, &    ! The effective layer marginal face areas with the southerly
    dvhdv_R, &    ! (_L), northerly (_R), and zero-barotropic (_0) test
    dvhdv_0, &    ! velocities [H L ~> m2 or kg m-1].
    v_L, v_R, &   ! The southerly (v_L), northerly (v_R), and zero-barotropic
    v_0        ! transport (v_0) layer test velocities [L T-1 ~> m s-1].
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
  real :: dx_Cv      ! The unblocked length of this v-face [L ~> m].
  real :: IareaT_S, IareaT_N ! The inverse areas of the tracer cells south and north of this face [L-2 ~> m-2].
  real :: IdyT_S, IdyT_N     ! The inverse meridional grid spacings south and north of this face [L-1 ~> m-1].
  real :: visc_rem_max_c ! A copy of visc_rem_max for this column [nondim].
  real :: h_marg_min ! A copy of CS%h_marg_min for use in the kernel [H ~> m or kg m-2].
  logical :: vol_CFL ! A copy of CS%vol_CFL for use in the kernel.
  logical :: do_c    ! A copy of do_ij for this column.
  integer :: i, j, k, nz, ii, JJ, iie, JJe

  nz = GV%ke ; Idt = 1.0 / dt
  iie = ieb - isb + 1 ; JJe = JebB - JsbB + 1
  min_visc_rem = 0.1 ; CFL_min = 1e-6
  h_marg_min = CS%h_marg_min ; vol_CFL = CS%vol_CFL

  !$omp target enter data map(alloc: dv0, zeros, vh_tmp)

 ! Diagnose the zero-transport correction, dv0.
  !$omp target teams distribute parallel do collapse(2)
  do JJ=1,JJe ; do ii=1,iie
    zeros(ii,JJ) = 0.0
  enddo ; enddo

  call meridional_flux_adjust(v, h_in, h_S, h_N, zeros, vh_tot_0, dvhdv_tot_0, dv0, &
                              dv_max_CFL, dv_min_CFL, dt, G, GV, US, CS, visc_rem, &
                              isb, ieb, JsbB, JebB, do_ij, por_face_areaV, &
                              nii, nJJB, vh_tmp)

  !   Determine the southerly- and northerly- fluxes.  Choose a sufficiently
  ! negative velocity correction for the northerly-flux, and a sufficiently
  ! positive correction for the southerly-flux.
  !$omp target teams distribute parallel do collapse(2) &
  !$omp   private(i, J, visc_rem_lim, v_L, v_R, v_0, vh_L, vh_R, vh_0, dvhdv_L, dvhdv_R, dvhdv_0, &
  !$omp           FA_0, FA_avg, dvL, dvR, dv_CFL, dv0_c, FAmt_L, FAmt_R, FAmt_0, vhtot_L, vhtot_R, &
  !$omp           dx_Cv, IareaT_S, IareaT_N, IdyT_S, IdyT_N, visc_rem_max_c, do_c)
  do JJ=1,JJe ; do ii=1,iie
    i = isb + ii - 1
    J = JsbB + JJ - 1
    dx_Cv = G%dx_Cv(i,J)
    IareaT_S = G%IareaT(i,J) ; IareaT_N = G%IareaT(i,J+1)
    IdyT_S = G%IdyT(i,J) ; IdyT_N = G%IdyT(i,J+1)
    visc_rem_max_c = visc_rem_max(ii,JJ) ; do_c = do_ij(ii,JJ) ; dv0_c = dv0(ii,JJ)
    dv_CFL = (CFL_min * Idt) * G%dyCv(i,J)
    dvR = min(0.0,dv0_c - dv_CFL)
    dvL = max(0.0,dv0_c + dv_CFL)
    FAmt_L = 0.0 ; FAmt_R = 0.0 ; FAmt_0 = 0.0
    vhtot_L = 0.0 ; vhtot_R = 0.0

    if (do_c) then
      do k=1,nz
        visc_rem_lim = max(visc_rem(ii,JJ,k), min_visc_rem*visc_rem_max_c)
        if (visc_rem_lim > 0.0) then ! This is almost always true for ocean points.
          if (v(i,J,k) + dvR*visc_rem_lim > -dv_CFL*visc_rem(ii,JJ,k)) &
            dvR = -(v(i,J,k) + dv_CFL*visc_rem(ii,JJ,k)) / visc_rem_lim
          if (v(i,J,k) + dvL*visc_rem_lim < dv_CFL*visc_rem(ii,JJ,k)) &
            dvL = -(v(i,J,k) - dv_CFL*visc_rem(ii,JJ,k)) / visc_rem_lim
        endif
      enddo

      do k=1,nz
        v_L = v(i,J,k) + dvL * visc_rem(ii,JJ,k)
        v_R = v(i,J,k) + dvR * visc_rem(ii,JJ,k)
        v_0 = v(i,J,k) + dv0_c * visc_rem(ii,JJ,k)

        call flux_elem(v_0, h_in(i,J,k), h_in(i,J+1,k), &
            h_S(i,J,k), h_S(i,J+1,k), h_N(i,J,k), h_N(i,J+1,k), &
            vh_0, dvhdv_0, visc_rem(ii,JJ,k), dx_Cv, &
            IareaT_S, IareaT_N, IdyT_S, IdyT_N, dt, &
            vol_CFL, por_face_areaV(i,J,k), h_marg_min)

        call flux_elem(v_L, h_in(i,J,k), h_in(i,J+1,k), &
            h_S(i,J,k), h_S(i,J+1,k), h_N(i,J,k), h_N(i,J+1,k), &
            vh_L, dvhdv_L, visc_rem(ii,JJ,k), dx_Cv, &
            IareaT_S, IareaT_N, IdyT_S, IdyT_N, dt, &
            vol_CFL, por_face_areaV(i,J,k), h_marg_min)

        call flux_elem(v_R, h_in(i,J,k), h_in(i,J+1,k), &
            h_S(i,J,k), h_S(i,J+1,k), h_N(i,J,k), h_N(i,J+1,k), &
            vh_R, dvhdv_R, visc_rem(ii,JJ,k), dx_Cv, &
            IareaT_S, IareaT_N, IdyT_S, IdyT_N, dt, &
            vol_CFL, por_face_areaV(i,J,k), h_marg_min)

        ! XXX: Add zero to force rounding and prevent FMA optimization.
        FAmt_0 = FAmt_0 + (dvhdv_0 + 0.)
        FAmt_L = FAmt_L + (dvhdv_L + 0.)
        FAmt_R = FAmt_R + (dvhdv_R + 0.)
        vhtot_L = vhtot_L + vh_L
        vhtot_R = vhtot_R + vh_R
      enddo

      FA_0 = FAmt_0
      FA_avg = FAmt_0

      if (dvL - dv0_c /= 0.) &
        FA_avg = vhtot_L / (dvL - dv0_c)

      if (FA_avg > max(FA_0, FAmt_L)) then
        FA_avg = max(FA_0, FAmt_L)
      elseif (FA_avg < min(FA_0, FAmt_L)) then
        FA_0 = FA_avg
      endif

      BT_cont%FA_v_S0(i,J) = FA_0
      BT_cont%FA_v_SS(i,J) = FAmt_L

      if (abs(FA_0-FAmt_L) <= 1e-12 * FA_0) then
        BT_cont%vBT_SS(i,J) = 0.
      else
        BT_cont%vBT_SS(i,J) = (1.5 * (dvL - dv0_c)) &
            * ((FAmt_L - FA_avg) / (FAmt_L - FA_0))
      endif

      FA_0 = FAmt_0
      FA_avg = FAmt_0

      if (dvR - dv0_c /= 0.) &
        FA_avg = vhtot_R / (dvR - dv0_c)

      if (FA_avg > max(FA_0, FAmt_R)) then
        FA_avg = max(FA_0, FAmt_R)
      elseif (FA_avg < min(FA_0, FAmt_R)) then
        FA_0 = FA_avg
      endif

      BT_cont%FA_v_N0(i,J) = FA_0 ; BT_cont%FA_v_NN(i,J) = FAmt_R

      if (abs(FAmt_R - FA_0) <= 1e-12 * FA_0) then
        BT_cont%vBT_NN(i,J) = 0.
      else
        BT_cont%vBT_NN(i,J) = (1.5 * (dvR - dv0_c)) &
            * ((FAmt_R - FA_avg) / (FAmt_R - FA_0))
      endif
    else
      BT_cont%FA_v_S0(i,J) = 0.
      BT_cont%FA_v_SS(i,J) = 0.
      BT_cont%FA_v_N0(i,J) = 0.
      BT_cont%FA_v_NN(i,J) = 0.
      BT_cont%vBT_SS(i,J) = 0.
      BT_cont%vBT_NN(i,J) = 0.
    endif
  enddo ; enddo

  !$omp target exit data map(release: dv0, zeros, vh_tmp)
end procedure set_merid_BT_cont

module procedure ratio_max
  !$omp declare target
  if (abs(a) > abs(maxrat*b)) then
    ratio = maxrat
  else
    ratio = a / b
  endif
end procedure ratio_max

end submodule MOM_continuity_PPM_s
