! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

#include "do_concurrent_compat.h"

!> Does full convective adjustment of unstable regions via a strong diffusivity.
module MOM_full_convection

use MOM_grid,              only : ocean_grid_type
use MOM_interface_heights, only : thickness_to_dz
use MOM_unit_scaling,      only : unit_scale_type
use MOM_variables,         only : thermo_var_ptrs
use MOM_verticalGrid,      only : verticalGrid_type
use MOM_EOS,               only : calculate_density_derivs

implicit none ; private

#include <MOM_memory.h>

public full_convection

contains

!> Calculate new temperatures and salinities that have been subject to full convective mixing.
subroutine full_convection(G, GV, US, h, tv, T_adj, S_adj, p_surf, Kddt_smooth, halo, nii, njj)
  type(ocean_grid_type),   intent(in)    :: G     !< The ocean's grid structure
  type(verticalGrid_type), intent(in)    :: GV    !< The ocean's vertical grid structure
  type(unit_scale_type),   intent(in)    :: US    !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in)    :: h     !< Layer thicknesses [H ~> m or kg m-2]
  type(thermo_var_ptrs),   intent(in)    :: tv    !< A structure pointing to various
                                                  !! thermodynamic variables
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(out)   :: T_adj !< Adjusted potential temperature [C ~> degC].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(out)   :: S_adj !< Adjusted salinity [S ~> ppt].
  real, dimension(:,:),    pointer       :: p_surf !< The pressure at the ocean surface [R L2 T-2 ~> Pa] (or NULL).
  real,                    intent(in)    :: Kddt_smooth  !< A smoothing vertical diffusivity
                                                  !! times a timestep [H Z ~> m2 or kg m-1].
  integer,                 intent(in)    :: halo  !< Halo width over which to compute
  integer,                 intent(in)    :: nii   !< i block size for array calculations.
  integer,                 intent(in)    :: njj   !< j block size for array calculations.

  ! Local variables
  real, dimension(nii,njj,SZK_(GV)+1) :: &
    dRho_dT, &  ! The derivative of density with temperature [R C-1 ~> kg m-3 degC-1]
    dRho_dS     ! The derivative of density with salinity [R S-1 ~> kg m-3 ppt-1].
  real :: dz(nii,njj,SZK_(GV)) ! Height change across layers [Z ~> m]
  real :: h_neglect     ! A thickness that is so small it is usually lost
                        ! in roundoff and can be neglected [H ~> m or kg m-2].
! logical :: use_EOS    ! If true, density is calculated from T & S using an equation of state.
  real, dimension(nii,njj,SZK0_(G)) :: &
    Te_a, & ! A partially updated temperature estimate including the influence from
            ! mixing with layers above rescaled by a factor of d_a [C ~> degC].
            ! This array is discretized on tracer cells, but contains an extra
            ! layer at the top for algorithmic convenience.
    Se_a    ! A partially updated salinity estimate including the influence from
            ! mixing with layers above rescaled by a factor of d_a [S ~> ppt].
            ! This array is discretized on tracer cells, but contains an extra
            ! layer at the top for algorithmic convenience.
  real, dimension(nii,njj,SZK_(GV)+1) :: &
    Te_b, & ! A partially updated temperature estimate including the influence from
            ! mixing with layers below rescaled by a factor of d_b [C ~> degC].
            ! This array is discretized on tracer cells, but contains an extra
            ! layer at the bottom for algorithmic convenience.
    Se_b    ! A partially updated salinity estimate including the influence from
            ! mixing with layers below rescaled by a factor of d_b [S ~> ppt].
            ! This array is discretized on tracer cells, but contains an extra
            ! layer at the bottom for algorithmic convenience.
  real, dimension(nii,njj,SZK_(GV)+1) :: &
    c_a, &  ! The fractional influence of the properties of the layer below
            ! in the final properties with a downward-first solver [nondim]
    d_a, &  ! The fractional influence of the properties of the layer in question
            ! and layers above in the final properties with a downward-first solver [nondim]
            ! d_a = 1.0 - c_a
    c_b, &  ! The fractional influence of the properties of the layer above
            ! in the final properties with a upward-first solver [nondim]
    d_b     ! The fractional influence of the properties of the layer in question
            ! and layers below in the final properties with a upward-first solver [nondim]
            ! d_b = 1.0 - c_b
  real, dimension(nii,njj,SZK_(GV)+1) :: &
    mix     !< The amount of mixing across the interface between layers [H ~> m or kg m-2].
  real :: mix_len  ! The length-scale of mixing, when it is active [H ~> m or kg m-2]
  real :: h_b, h_a ! The thicknesses of the layers above and below an interface [H ~> m or kg m-2]
  real :: b_b, b_a ! Inverse pivots used by the tridiagonal solver [H-1 ~> m-1 or m2 kg-1].

  logical, dimension(nii,njj) :: do_ij ! Do more work on this column.
  logical, dimension(nii,njj) :: last_down ! The last setup pass was downward.
  integer, dimension(nii,njj) :: change_ct ! The number of interfaces where the
                         ! mixing has changed this iteration.
  integer :: changed_col ! The number of columns whose mixing changed.
  integer :: i, j, k, is, ie, js, je, nz, itt, ii, jj, isb, ieb, jsb, jeb, iie, jje

  is = G%isc-halo ; ie = G%iec+halo ; js = G%jsc-halo ; je = G%jec+halo
  nz = GV%ke

  if (.not.associated(tv%eqn_of_state)) return

  h_neglect = GV%H_subroundoff
  mix_len = (1.0e20 * nz) * (G%max_depth * US%Z_to_m * GV%m_to_H)

  !$omp target enter data &
  !$omp   map(alloc: dRho_dT, dRho_dS, dz, Te_a, Se_a, Te_b, Se_b, c_a, d_a, c_b, d_b, mix, do_ij, &
  !$omp     last_down, change_ct)

  do jsb=js,je,njj ; do isb=is,ie,nii
    ieb = min(isb + nii - 1, ie)
    jeb = min(jsb + njj - 1, je)
    iie = ieb - isb + 1
    jje = jeb - jsb + 1

    do concurrent (k=1:nz+1, jj=1:jje, ii=1:iie)
      mix(ii,jj,k) = 0.0 ; d_b(ii,jj,k) = 1.0
      ! These would be Te_b(ii,jj,k) = tv%T(i,j,k), etc., but the values are not used
      Te_b(ii,jj,k) = 0.0 ; Se_b(ii,jj,k) = 0.0
    enddo

    ! Find the vertical distances across layers.
    call thickness_to_dz(h, tv, dz, G, GV, US, niB=nii, njB=njj, do_offload=.true., &
                         i_lo=isb, i_hi=ieb, j_lo=jsb, j_hi=jeb)

    call smoothed_dRdT_dRdS(h, dz, tv, Kddt_smooth, dRho_dT, dRho_dS, G, GV, US, nii, njj, &
                            isb, ieb, jsb, jeb, p_surf)

    do concurrent (jj=1:jje, ii=1:iie) DO_LOCALITY(local(i, j))
      i = isb + ii - 1
      j = jsb + jj - 1

      do_ij(ii,jj) = (G%mask2dT(i,j) > 0.0)

      d_a(ii,jj,1) = 1.0
      last_down(ii,jj) = .true. ! This is set for debuggers.
      ! These are extra values are used for convenience in the stability test
      Te_a(ii,jj,0) = 0.0 ; Se_a(ii,jj,0) = 0.0
    enddo

    do itt=1,nz ! At least 2 interfaces will change with each full pass, or the
                ! iterations stop, so the maximum count of nz is very conservative.

      do concurrent (jj=1:jje, ii=1:iie)
        change_ct(ii,jj) = 0
      enddo
      ! Move down the water column, finding unstable interfaces, and building up the
      ! temporary arrays for the tridiagonal solver.
      do concurrent (jj=1:jje) DO_LOCALITY(local(i, j, h_a, h_b, b_a))
        do K=2,nz
          do concurrent (ii=1:iie)
            i = isb + ii - 1
            j = jsb + jj - 1

            if (do_ij(ii,jj)) then
              h_a = h(i,j,k-1) + h_neglect ; h_b = h(i,j,k) + h_neglect
              if (mix(ii,jj,K) <= 0.0) then
                if (is_unstable(dRho_dT(ii,jj,K), dRho_dS(ii,jj,K), h_a, h_b, mix(ii,jj,K-1), mix(ii,jj,K+1), &
                                tv%T(i,j,k-1), tv%T(i,j,k), tv%S(i,j,k-1), tv%S(i,j,k), &
                                Te_a(ii,jj,k-2), Te_b(ii,jj,k+1), Se_a(ii,jj,k-2), Se_b(ii,jj,k+1), &
                                d_a(ii,jj,K-1), d_b(ii,jj,K+1))) then
                  mix(ii,jj,K) = mix_len
                  change_ct(ii,jj) = change_ct(ii,jj) + 1
                endif
              endif

              b_a = 1.0 / ((h_a + d_a(ii,jj,K-1)*mix(ii,jj,K-1)) + mix(ii,jj,K))
              if (mix(ii,jj,K) <= 0.0) then
                c_a(ii,jj,K) = 0.0 ; d_a(ii,jj,K) = 1.0
              else
                d_a(ii,jj,K) = b_a * (h_a + d_a(ii,jj,K-1)*mix(ii,jj,K-1)) ! = 1.0-c_a(ii,jj,K)
                c_a(ii,jj,K) = 1.0 ; if (d_a(ii,jj,K) > epsilon(b_a)) c_a(ii,jj,K) = b_a * mix(ii,jj,K)
              endif

              if (K>2) then
                Te_a(ii,jj,k-1) = b_a * (h_a*tv%T(i,j,k-1) + mix(ii,jj,K-1)*Te_a(ii,jj,k-2))
                Se_a(ii,jj,k-1) = b_a * (h_a*tv%S(i,j,k-1) + mix(ii,jj,K-1)*Se_a(ii,jj,k-2))
              else
                Te_a(ii,jj,k-1) = b_a * (h_a*tv%T(i,j,k-1))
                Se_a(ii,jj,k-1) = b_a * (h_a*tv%S(i,j,k-1))
              endif
            endif
          enddo
        enddo
      enddo

      ! Determine which columns might have further instabilities.
      changed_col = 0
      do concurrent (jj=1:jje, ii=1:iie, do_ij(ii,jj)) DO_LOCALITY(reduce(+:changed_col))
        if (change_ct(ii,jj) == 0) then
          last_down(ii,jj) = .true. ; do_ij(ii,jj) = .false.
        else
          changed_col = changed_col + 1 ; change_ct(ii,jj) = 0
        endif
      enddo
      if (changed_col == 0) exit ! No more columns are unstable.

      ! This is the same as above, but with the direction reversed (bottom to top)
      do concurrent (jj=1:jje) DO_LOCALITY(local(i, j, h_a, h_b, b_b))
        do K=nz,2,-1
          do concurrent (ii=1:iie)
            i = isb + ii - 1
            j = jsb + jj - 1

            if (do_ij(ii,jj)) then
              h_a = h(i,j,k-1) + h_neglect ; h_b = h(i,j,k) + h_neglect
              if (mix(ii,jj,K) <= 0.0) then
                if (is_unstable(dRho_dT(ii,jj,K), dRho_dS(ii,jj,K), h_a, h_b, mix(ii,jj,K-1), mix(ii,jj,K+1), &
                                tv%T(i,j,k-1), tv%T(i,j,k), tv%S(i,j,k-1), tv%S(i,j,k), &
                                Te_a(ii,jj,k-2), Te_b(ii,jj,k+1), Se_a(ii,jj,k-2), Se_b(ii,jj,k+1), &
                                d_a(ii,jj,K-1), d_b(ii,jj,K+1))) then
                  mix(ii,jj,K) = mix_len
                  change_ct(ii,jj) = change_ct(ii,jj) + 1
                endif
              endif

              b_b = 1.0 / ((h_b + d_b(ii,jj,K+1)*mix(ii,jj,K+1)) + mix(ii,jj,K))
              if (mix(ii,jj,K) <= 0.0) then
                c_b(ii,jj,K) = 0.0 ; d_b(ii,jj,K) = 1.0
              else
                d_b(ii,jj,K) = b_b * (h_b + d_b(ii,jj,K+1)*mix(ii,jj,K+1)) ! = 1.0-c_b(ii,jj,K)
                c_b(ii,jj,K) = 1.0 ; if (d_b(ii,jj,K) > epsilon(b_b)) c_b(ii,jj,K) = b_b * mix(ii,jj,K)
              endif

              if (k<nz) then
                Te_b(ii,jj,k) = b_b * (h_b*tv%T(i,j,k) + mix(ii,jj,K+1)*Te_b(ii,jj,k+1))
                Se_b(ii,jj,k) = b_b * (h_b*tv%S(i,j,k) + mix(ii,jj,K+1)*Se_b(ii,jj,k+1))
              else
                Te_b(ii,jj,k) = b_b * (h_b*tv%T(i,j,k))
                Se_b(ii,jj,k) = b_b * (h_b*tv%S(i,j,k))
              endif
            endif
          enddo
        enddo
      enddo

      ! Determine which columns might have further instabilities.
      changed_col = 0
      do concurrent (jj=1:jje, ii=1:iie, do_ij(ii,jj)) DO_LOCALITY(reduce(+:changed_col))
        if (change_ct(ii,jj) == 0) then
          last_down(ii,jj) = .false. ; do_ij(ii,jj) = .false.
        else
          changed_col = changed_col + 1 ; change_ct(ii,jj) = 0
        endif
      enddo
      if (changed_col == 0) exit ! No more columns are unstable.

    enddo  ! End of iterations, all columns are now stable.

    ! Do the final return pass on the columns where the penultimate pass was downward.
    do concurrent (jj=1:jje, ii=1:iie) DO_LOCALITY(local(i, j))
      i = isb + ii - 1
      j = jsb + jj - 1
      do_ij(ii,jj) = ((G%mask2dT(i,j) > 0.0) .and. last_down(ii,jj))
    enddo
    do concurrent (jj=1:jje, ii=1:iie, do_ij(ii,jj)) DO_LOCALITY(local(i, j, h_a, b_a))
      i = isb + ii - 1
      j = jsb + jj - 1
      h_a = h(i,j,nz) + h_neglect
      b_a = 1.0 / (h_a + d_a(ii,jj,nz)*mix(ii,jj,nz))
      T_adj(i,j,nz) = b_a * (h_a*tv%T(i,j,nz) + mix(ii,jj,nz)*Te_a(ii,jj,nz-1))
      S_adj(i,j,nz) = b_a * (h_a*tv%S(i,j,nz) + mix(ii,jj,nz)*Se_a(ii,jj,nz-1))
    enddo
    do concurrent (jj=1:jje) DO_LOCALITY(local(i, j))
      do k=nz-1,1,-1
        do concurrent (ii=1:iie, do_ij(ii,jj))
          i = isb + ii - 1
          j = jsb + jj - 1
          T_adj(i,j,k) = Te_a(ii,jj,k) + c_a(ii,jj,K+1)*T_adj(i,j,k+1)
          S_adj(i,j,k) = Se_a(ii,jj,k) + c_a(ii,jj,K+1)*S_adj(i,j,k+1)
        enddo
      enddo
    enddo

    ! Keeping this confuses NVHPC - it thinks k is needed after the loop on 144
    ! do concurrent (jj=1:jje, ii=1:iie, do_ij(ii,jj)) DO_LOCALITY(reduce(max:k))
    !   k = 1 ! A hook for debugging.
    ! enddo

    ! Do the final return pass on the columns where the penultimate pass was upward.
    ! Also do a simple copy of T & S values on land points.
    do concurrent (jj=1:jje, ii=1:iie) DO_LOCALITY(local(i, j, h_b, b_b))
      i = isb + ii - 1
      j = jsb + jj - 1
      do_ij(ii,jj) = ((G%mask2dT(i,j) > 0.0) .and. .not.last_down(ii,jj))
      if (do_ij(ii,jj)) then
        h_b = h(i,j,1) + h_neglect
        b_b = 1.0 / (h_b + d_b(ii,jj,2)*mix(ii,jj,2))
        T_adj(i,j,1) = b_b * (h_b*tv%T(i,j,1) + mix(ii,jj,2)*Te_b(ii,jj,2))
        S_adj(i,j,1) = b_b * (h_b*tv%S(i,j,1) + mix(ii,jj,2)*Se_b(ii,jj,2))
      elseif (G%mask2dT(i,j) <= 0.0) then
        T_adj(i,j,1) = tv%T(i,j,1) ; S_adj(i,j,1) = tv%S(i,j,1)
      endif
    enddo
    do concurrent (jj=1:jje) DO_LOCALITY(local(i, j))
      do k=2,nz
        do concurrent (ii=1:iie)
          i = isb + ii - 1
          j = jsb + jj - 1
          if (do_ij(ii,jj)) then
            T_adj(i,j,k) = Te_b(ii,jj,k) + c_b(ii,jj,K)*T_adj(i,j,k-1)
            S_adj(i,j,k) = Se_b(ii,jj,k) + c_b(ii,jj,K)*S_adj(i,j,k-1)
          elseif (G%mask2dT(i,j) <= 0.0) then
            T_adj(i,j,k) = tv%T(i,j,k) ; S_adj(i,j,k) = tv%S(i,j,k)
          endif
        enddo
      enddo
    enddo

    ! For some reason causes present table errors on GPU
    ! do concurrent (jj=1:jje, ii=1:iie, do_ij(ii,jj)) DO_LOCALITY(reduce(max:k))
    !   k = 1 ! A hook for debugging.
    ! enddo

  enddo ; enddo ! ij block loops

  !$omp target exit data &
  !$omp   map(release: dRho_dT, dRho_dS, dz, Te_a, Se_a, Te_b, Se_b, c_a, d_a, c_b, d_b, mix, do_ij, &
  !$omp     last_down, change_ct)

  k = 1 ! A hook for debugging.

  ! The following set of expressions for the final values are derived from the partial
  ! updates for the estimated temperatures and salinities around an interface, then directly
  ! solving for the final temperatures and salinities.  They are here for later reference
  ! and to document an intermediate step in the stability calculation.
    ! hp_a = (h_a + d_a(i,j,K-1)*mix(i,j,K-1))
    ! hp_b = (h_b + d_b(i,j,K+1)*mix(i,j,K+1))
    ! b2_c = 1.0 / (hp_a*hp_b + (hp_a + hp_b) * mix(i,j,K))
    ! Th_a = h_a*tv%T(i,j,k-1) + mix(i,j,K-1)*Te_a(i,j,k-2)
    ! Th_b = h_b*tv%T(i,j,k)   + mix(i,j,K+1)*Te_b(i,j,k+1)
    ! T_fin(i,k)   = ( (hp_a + mix(i,j,K)) * Th_b  + Th_a * mix(i,j,K) ) * b2_c
    ! T_fin(i,k-1) = ( (hp_b + mix(i,j,K)) * Th_a  + Th_b * mix(i,j,K) ) * b2_c
    ! Sh_a = h_a*tv%S(i,j,k-1) + mix(i,j,K-1)*Se_a(i,j,k-2)
    ! Sh_b = h_b*tv%S(i,j,k)   + mix(i,j,K+1)*Se_b(i,j,k+1)
    ! S_fin(i,k)   = ( (hp_a + mix(i,j,K)) * Sh_b  + Sh_a * mix(i,j,K) ) * b2_c
    ! S_fin(i,k-1) = ( (hp_b + mix(i,j,K)) * Sh_a  + Sh_b * mix(i,j,K) ) * b2_c

end subroutine full_convection

!> This function returns True if the profiles around the given interface will be
!! statically unstable after mixing above below.  The arguments are the ocean state
!! above and below, including partial calculations from a tridiagonal solver.
pure function is_unstable(dRho_dT, dRho_dS, h_a, h_b, mix_A, mix_B, T_a, T_b, S_a, S_b, &
                          Te_aa, Te_bb, Se_aa, Se_bb, d_A, d_B)
  real, intent(in) :: dRho_dT !< The derivative of in situ density with temperature [R C-1 ~> kg m-3 degC-1]
  real, intent(in) :: dRho_dS !< The derivative of in situ density with salinity [R S-1 ~> kg m-3 ppt-1]
  real, intent(in) :: h_a     !< The thickness of the layer above [H ~> m or kg m-2]
  real, intent(in) :: h_b     !< The thickness of the layer below [H ~> m or kg m-2]
  real, intent(in) :: mix_A   !< The time integrated mixing rate of the interface above [H ~> m or kg m-2]
  real, intent(in) :: mix_B   !< The time integrated mixing rate of the interface below [H ~> m or kg m-2]
  real, intent(in) :: T_a     !< The initial temperature of the layer above [C ~> degC]
  real, intent(in) :: T_b     !< The initial temperature of the layer below [C ~> degC]
  real, intent(in) :: S_a     !< The initial salinity of the layer below [S ~> ppt]
  real, intent(in) :: S_b     !< The initial salinity of the layer below [S ~> ppt]
  real, intent(in) :: Te_aa   !< The estimated temperature two layers above rescaled by d_A [C ~> degC]
  real, intent(in) :: Te_bb   !< The estimated temperature two layers below rescaled by d_B [C ~> degC]
  real, intent(in) :: Se_aa   !< The estimated salinity two layers above rescaled by d_A [S ~> ppt]
  real, intent(in) :: Se_bb   !< The estimated salinity two layers below rescaled by d_B [S ~> ppt]
  real, intent(in) :: d_A     !< The rescaling dependency across the interface above [nondim]
  real, intent(in) :: d_B     !< The rescaling dependency across the interface below [nondim]
  logical :: is_unstable !< The return value, true if the profile is statically unstable
                         !! around the interface in question.

  ! These expressions for the local stability are long, but they have been carefully
  ! grouped for accuracy even when the mixing rates are huge or tiny, and common
  ! positive definite factors that would appear in the final expression for the
  ! locally referenced potential density difference across an interface have been omitted.
  is_unstable = (dRho_dT * ((h_a * h_b * (T_b - T_a) + &
                             mix_A*mix_B * (d_A*Te_bb - d_B*Te_aa)) + &
                            (h_a*mix_B * (Te_bb - d_B*T_a) + &
                             h_b*mix_A * (d_A*T_b - Te_aa)) ) + &
                 dRho_dS * ((h_a * h_b * (S_b - S_a) + &
                             mix_A*mix_B * (d_A*Se_bb - d_B*Se_aa)) + &
                            (h_a*mix_B * (Se_bb - d_B*S_a) + &
                             h_b*mix_A * (d_A*S_b - Se_aa)) ) < 0.0)
end function is_unstable

!> Returns the partial derivatives of locally referenced potential density with
!! temperature and salinity after the properties have been smoothed with a small
!! constant diffusivity.
subroutine smoothed_dRdT_dRdS(h, dz, tv, Kddt, dR_dT, dR_dS, G, GV, US, nii, njj, isb, ieb, jsb, jeb, p_surf)
  type(ocean_grid_type),   intent(in)  :: G    !< The ocean's grid structure
  type(verticalGrid_type), intent(in)  :: GV   !< The ocean's vertical grid structure
  integer,                 intent(in)  :: nii  !< i block size for array calculations.
  integer,                 intent(in)  :: njj  !< j block size for array calculations.
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in)  :: h    !< Layer thicknesses [H ~> m or kg m-2]
  real, dimension(nii,njj,SZK_(GV)), &
                           intent(in)  :: dz   !< Height change across layers [Z ~> m]
  type(thermo_var_ptrs),   intent(in)  :: tv   !< A structure pointing to various
                                               !! thermodynamic variables
  real,                    intent(in)  :: Kddt !< A diffusivity times a time increment [H Z ~> m2 or kg m-1].
  real, dimension(nii,njj,SZK_(GV)+1), &
                           intent(out) :: dR_dT !< Derivative of locally referenced
                                               !! potential density with temperature [R C-1 ~> kg m-3 degC-1]
  real, dimension(nii,njj,SZK_(GV)+1), &
                           intent(out) :: dR_dS !< Derivative of locally referenced
                                               !! potential density with salinity [R S-1 ~> kg m-3 ppt-1]
  type(unit_scale_type),   intent(in)  :: US   !< A dimensional unit scaling type
  integer,                 intent(in)  :: isb  !< Start of i index range.
  integer,                 intent(in)  :: ieb  !< End of i index range.
  integer,                 intent(in)  :: jsb  !< Start of j index range.
  integer,                 intent(in)  :: jeb  !< End of j index range.
  real, dimension(:,:),    pointer     :: p_surf !< The pressure at the ocean surface [R L2 T-2 ~> Pa].

  ! Local variables
  real :: mix(nii,njj,SZK_(GV)+1)  ! The diffusive mixing length (kappa*dt)/dz
                                   ! between layers within in a timestep [H ~> m or kg m-2].
  real :: b1(nii,njj)              ! A tridiagonal solver variable [H-1 ~> m-1 or m2 kg-1]
  real :: d1(nii,njj)              ! A tridiagonal solver variable [nondim]
  real :: c1(nii,njj,SZK_(GV))     ! A tridiagonal solver variable [nondim]
  real :: T_f(nii,njj,SZK_(GV))    ! Filtered temperatures [C ~> degC]
  real :: S_f(nii,njj,SZK_(GV))    ! Filtered salinities [S ~> ppt]
  real :: pres(nii,njj)            ! Interface pressures [R L2 T-2 ~> Pa].
  real :: T_EOS(nii,njj)           ! Filtered and vertically averaged temperatures [C ~> degC]
  real :: S_EOS(nii,njj)           ! Filtered and vertically averaged salinities [S ~> ppt]
  real :: kap_dt_x2                ! The product of 2*kappa*dt [H Z ~> m2 or kg m-1].
  real :: dz_neglect, h0           ! A negligible vertical distances [Z ~> m]
  real :: h_neglect                ! A negligible thickness to allow for zero thicknesses
                                   ! [H ~> m or kg m-2].
  real :: h_tr                     ! The thickness at tracer points, plus h_neglect [H ~> m or kg m-2].
  integer, dimension(2,2) :: EOSdom  ! The block-local i- and j-domain for the equation of state
  integer :: i, j, k, nz, ii, jj, iie, jje

  nz = GV%ke
  iie = ieb - isb + 1 ; jje = jeb - jsb + 1

  h_neglect = GV%H_subroundoff
  dz_neglect = GV%dz_subroundoff
  kap_dt_x2 = 2.0*Kddt

  !$omp target enter data map(alloc: mix, b1, d1, c1, T_f, S_f, pres, T_EOS, S_EOS)

  if (Kddt <= 0.0) then
    do concurrent (k=1:nz, jj=1:jje, ii=1:iie) DO_LOCALITY(local(i, j))
      i = isb + ii - 1
      j = jsb + jj - 1
      T_f(ii,jj,k) = tv%T(i,j,k) ; S_f(ii,jj,k) = tv%S(i,j,k)
    enddo
  else
    h0 = 1.0e-16*sqrt(GV%H_to_m*US%m_to_Z*Kddt) + dz_neglect
    do concurrent (jj=1:jje, ii=1:iie) DO_LOCALITY(local(i, j, h_tr))
      i = isb + ii - 1
      j = jsb + jj - 1

      mix(ii,jj,2) = kap_dt_x2 / ((dz(ii,jj,1)+dz(ii,jj,2)) + h0)

      h_tr = h(i,j,1) + h_neglect
      b1(ii,jj) = 1.0 / (h_tr + mix(ii,jj,2))
      d1(ii,jj) = b1(ii,jj) * h(i,j,1)
      T_f(ii,jj,1) = (b1(ii,jj)*h_tr)*tv%T(i,j,1)
      S_f(ii,jj,1) = (b1(ii,jj)*h_tr)*tv%S(i,j,1)
    enddo
    do concurrent (jj=1:jje) DO_LOCALITY(local(i, j, h_tr))
      do k=2,nz-1
        do concurrent (ii=1:iie)
          i = isb + ii - 1
          j = jsb + jj - 1

          mix(ii,jj,K+1) = kap_dt_x2 / ((dz(ii,jj,k)+dz(ii,jj,k+1)) + h0)

          c1(ii,jj,k) = mix(ii,jj,K) * b1(ii,jj)
          h_tr = h(i,j,k) + h_neglect
          b1(ii,jj) = 1.0 / ((h_tr + d1(ii,jj)*mix(ii,jj,K)) + mix(ii,jj,K+1))
          d1(ii,jj) = b1(ii,jj) * (h_tr + d1(ii,jj)*mix(ii,jj,K))
          T_f(ii,jj,k) = b1(ii,jj) * (h_tr*tv%T(i,j,k) + mix(ii,jj,K)*T_f(ii,jj,k-1))
          S_f(ii,jj,k) = b1(ii,jj) * (h_tr*tv%S(i,j,k) + mix(ii,jj,K)*S_f(ii,jj,k-1))
        enddo
      enddo
    enddo
    do concurrent (jj=1:jje, ii=1:iie) DO_LOCALITY(local(i, j, h_tr))
      i = isb + ii - 1
      j = jsb + jj - 1

      c1(ii,jj,nz) = mix(ii,jj,nz) * b1(ii,jj)
      h_tr = h(i,j,nz) + h_neglect
      b1(ii,jj) = 1.0 / (h_tr + d1(ii,jj)*mix(ii,jj,nz))
      T_f(ii,jj,nz) = b1(ii,jj) * (h_tr*tv%T(i,j,nz) + mix(ii,jj,nz)*T_f(ii,jj,nz-1))
      S_f(ii,jj,nz) = b1(ii,jj) * (h_tr*tv%S(i,j,nz) + mix(ii,jj,nz)*S_f(ii,jj,nz-1))
    enddo
    do concurrent (jj=1:jje)
      do k=nz-1,1,-1
        do concurrent (ii=1:iie)
          T_f(ii,jj,k) = T_f(ii,jj,k) + c1(ii,jj,k+1)*T_f(ii,jj,k+1)
          S_f(ii,jj,k) = S_f(ii,jj,k) + c1(ii,jj,k+1)*S_f(ii,jj,k+1)
        enddo
      enddo
    enddo
  endif

  if (associated(p_surf)) then
    do concurrent (jj=1:jje, ii=1:iie) DO_LOCALITY(local(i, j))
      i = isb + ii - 1
      j = jsb + jj - 1
      pres(ii,jj) = p_surf(i,j)
    enddo
  else
    do concurrent (jj=1:jje, ii=1:iie) ; pres(ii,jj) = 0.0 ; enddo
  endif

  ! TODO: implement k-blocking so full k can be passed to calculate_density_derivs
  !       to maximize occupancy for GPU.
  EOSdom(1,1) = 1 ; EOSdom(1,2) = iie
  EOSdom(2,1) = 1 ; EOSdom(2,2) = jje
  call calculate_density_derivs(T_f(:,:,1), S_f(:,:,1), pres, dR_dT(:,:,1), dR_dS(:,:,1), tv%eqn_of_state, EOSdom)
  do concurrent (jj=1:jje, ii=1:iie) DO_LOCALITY(local(i, j))
    i = isb + ii - 1
    j = jsb + jj - 1
    pres(ii,jj) = pres(ii,jj) + h(i,j,1)*(GV%H_to_RZ*GV%g_Earth)
  enddo
  do K=2,nz
    do concurrent (jj=1:jje, ii=1:iie)
      T_EOS(ii,jj) = 0.5*(T_f(ii,jj,k-1) + T_f(ii,jj,k))
      S_EOS(ii,jj) = 0.5*(S_f(ii,jj,k-1) + S_f(ii,jj,k))
    enddo
    call calculate_density_derivs(T_EOS, S_EOS, pres, dR_dT(:,:,K), dR_dS(:,:,K), tv%eqn_of_state, EOSdom)
    do concurrent (jj=1:jje, ii=1:iie) DO_LOCALITY(local(i, j))
      i = isb + ii - 1
      j = jsb + jj - 1

      pres(ii,jj) = pres(ii,jj) + h(i,j,k)*(GV%H_to_RZ*GV%g_Earth)
    enddo
  enddo
  call calculate_density_derivs(T_f(:,:,nz), S_f(:,:,nz), pres, dR_dT(:,:,nz+1), dR_dS(:,:,nz+1), &
                                tv%eqn_of_state, EOSdom)

  !$omp target exit data map(release: mix, b1, d1, c1, T_f, S_f, pres, T_EOS, S_EOS)

end subroutine smoothed_dRdT_dRdS

end module MOM_full_convection
