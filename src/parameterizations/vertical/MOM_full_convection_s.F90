submodule (MOM_full_convection) MOM_full_convection_s
#include <MOM_memory.h>
  implicit none
contains
module procedure full_convection
  real, dimension(SZI_(G),SZK_(GV)+1) :: &
    dRho_dT, &  ! The derivative of density with temperature [R C-1 ~> kg m-3 degC-1]
    dRho_dS     ! The derivative of density with salinity [R S-1 ~> kg m-3 ppt-1].
  real :: dz(SZI_(G),SZK_(GV)) ! Height change across layers [Z ~> m]
  real :: h_neglect     ! A thickness that is so small it is usually lost
  real, dimension(SZI_(G),SZK0_(G)) :: &
    Te_a, & ! A partially updated temperature estimate including the influence from
            ! mixing with layers above rescaled by a factor of d_a [C ~> degC].
    Se_a    ! A partially updated salinity estimate including the influence from
            ! mixing with layers above rescaled by a factor of d_a [S ~> ppt].
            ! This array is discretized on tracer cells, but contains an extra
            ! layer at the top for algorithmic convenience.
  real, dimension(SZI_(G),SZK_(GV)+1) :: &
    Te_b, & ! A partially updated temperature estimate including the influence from
            ! mixing with layers below rescaled by a factor of d_b [C ~> degC].
            ! This array is discretized on tracer cells, but contains an extra
            ! layer at the bottom for algorithmic convenience.
    Se_b    ! A partially updated salinity estimate including the influence from
            ! mixing with layers below rescaled by a factor of d_b [S ~> ppt].
            ! This array is discretized on tracer cells, but contains an extra
            ! layer at the bottom for algorithmic convenience.
  real, dimension(SZI_(G),SZK_(GV)+1) :: &
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
  real, dimension(SZI_(G),SZK_(GV)+1) :: &
    mix     !< The amount of mixing across the interface between layers [H ~> m or kg m-2].
  real :: mix_len  ! The length-scale of mixing, when it is active [H ~> m or kg m-2]
  real :: h_b, h_a ! The thicknesses of the layers above and below an interface [H ~> m or kg m-2]
  real :: b_b, b_a ! Inverse pivots used by the tridiagonal solver [H-1 ~> m-1 or m2 kg-1].

  logical, dimension(SZI_(G)) :: do_i ! Do more work on this column.
  logical, dimension(SZI_(G)) :: last_down ! The last setup pass was downward.
  integer, dimension(SZI_(G)) :: change_ct ! The number of interfaces where the
                         ! mixing has changed this iteration.
  integer :: changed_col ! The number of columns whose mixing changed.
  integer :: i, j, k, is, ie, js, je, nz, itt

  is = G%isc-halo ; ie = G%iec+halo ; js = G%jsc-halo ; je = G%jec+halo
  nz = GV%ke

  if (.not.associated(tv%eqn_of_state)) return

  h_neglect = GV%H_subroundoff
  mix_len = (1.0e20 * nz) * (G%max_depth * US%Z_to_m * GV%m_to_H)

  do j=js,je
    mix(:,:) = 0.0 ; d_b(:,:) = 1.0
    ! These would be Te_b(:,:) = tv%T(:,j,:), etc., but the values are not used
    Te_b(:,:) = 0.0 ; Se_b(:,:) = 0.0

    ! Find the vertical distances across layers.
    call thickness_to_dz(h, tv, dz, j, G, GV, halo_size=halo)

    call smoothed_dRdT_dRdS(h, dz, tv, Kddt_smooth, dRho_dT, dRho_dS, G, GV, US, j, p_surf, halo)

    do i=is,ie
      do_i(i) = (G%mask2dT(i,j) > 0.0)

      d_a(i,1) = 1.0
      last_down(i) = .true. ! This is set for debuggers.
      ! These are extra values are used for convenience in the stability test
      Te_a(i,0) = 0.0 ; Se_a(i,0) = 0.0
    enddo

    do itt=1,nz ! At least 2 interfaces will change with each full pass, or the
                ! iterations stop, so the maximum count of nz is very conservative.

      do i=is,ie ; change_ct(i) = 0 ; enddo
      ! Move down the water column, finding unstable interfaces, and building up the
      ! temporary arrays for the tridiagonal solver.
      do K=2,nz ; do i=is,ie ; if (do_i(i)) then

        h_a = h(i,j,k-1) + h_neglect ; h_b = h(i,j,k) + h_neglect
        if (mix(i,K) <= 0.0) then
          if (is_unstable(dRho_dT(i,K), dRho_dS(i,K), h_a, h_b, mix(i,K-1), mix(i,K+1), &
                          tv%T(i,j,k-1), tv%T(i,j,k), tv%S(i,j,k-1), tv%S(i,j,k), &
                          Te_a(i,k-2), Te_b(i,k+1), Se_a(i,k-2), Se_b(i,k+1), &
                          d_a(i,K-1), d_b(i,K+1))) then
            mix(i,K) = mix_len
            change_ct(i) = change_ct(i) + 1
          endif
        endif

        b_a = 1.0 / ((h_a + d_a(i,K-1)*mix(i,K-1)) + mix(i,K))
        if (mix(i,K) <= 0.0) then
          c_a(i,K) = 0.0 ; d_a(i,K) = 1.0
        else
          d_a(i,K) = b_a * (h_a + d_a(i,K-1)*mix(i,K-1)) ! = 1.0-c_a(i,K)
          c_a(i,K) = 1.0 ; if (d_a(i,K) > epsilon(b_a)) c_a(i,K) = b_a * mix(i,K)
        endif

        if (K>2) then
          Te_a(i,k-1) = b_a * (h_a*tv%T(i,j,k-1) + mix(i,K-1)*Te_a(i,k-2))
          Se_a(i,k-1) = b_a * (h_a*tv%S(i,j,k-1) + mix(i,K-1)*Se_a(i,k-2))
        else
          Te_a(i,k-1) = b_a * (h_a*tv%T(i,j,k-1))
          Se_a(i,k-1) = b_a * (h_a*tv%S(i,j,k-1))
        endif
      endif ; enddo ; enddo

      ! Determine which columns might have further instabilities.
      changed_col = 0
      do i=is,ie ; if (do_i(i)) then
        if (change_ct(i) == 0) then
          last_down(i) = .true. ; do_i(i) = .false.
        else
          changed_col = changed_col + 1 ; change_ct(i) = 0
        endif
      endif ; enddo
      if (changed_col == 0) exit ! No more columns are unstable.

      ! This is the same as above, but with the direction reversed (bottom to top)
      do K=nz,2,-1 ; do i=is,ie ; if (do_i(i)) then

        h_a = h(i,j,k-1) + h_neglect ; h_b = h(i,j,k) + h_neglect
        if (mix(i,K) <= 0.0) then
          if (is_unstable(dRho_dT(i,K), dRho_dS(i,K), h_a, h_b, mix(i,K-1), mix(i,K+1), &
                          tv%T(i,j,k-1), tv%T(i,j,k), tv%S(i,j,k-1), tv%S(i,j,k), &
                          Te_a(i,k-2), Te_b(i,k+1), Se_a(i,k-2), Se_b(i,k+1), &
                          d_a(i,K-1), d_b(i,K+1))) then
            mix(i,K) = mix_len
            change_ct(i) = change_ct(i) + 1
          endif
        endif

        b_b = 1.0 / ((h_b + d_b(i,K+1)*mix(i,K+1)) + mix(i,K))
        if (mix(i,K) <= 0.0) then
          c_b(i,K) = 0.0 ; d_b(i,K) = 1.0
        else
          d_b(i,K) = b_b * (h_b + d_b(i,K+1)*mix(i,K+1)) ! = 1.0-c_b(i,K)
          c_b(i,K) = 1.0 ; if (d_b(i,K) > epsilon(b_b)) c_b(i,K) = b_b * mix(i,K)
        endif

        if (k<nz) then
          Te_b(i,k) = b_b * (h_b*tv%T(i,j,k) + mix(i,K+1)*Te_b(i,k+1))
          Se_b(i,k) = b_b * (h_b*tv%S(i,j,k) + mix(i,K+1)*Se_b(i,k+1))
        else
          Te_b(i,k) = b_b * (h_b*tv%T(i,j,k))
          Se_b(i,k) = b_b * (h_b*tv%S(i,j,k))
        endif
      endif ; enddo ; enddo

      ! Determine which columns might have further instabilities.
      changed_col = 0
      do i=is,ie ; if (do_i(i)) then
        if (change_ct(i) == 0) then
          last_down(i) = .false. ; do_i(i) = .false.
        else
          changed_col = changed_col + 1 ; change_ct(i) = 0
        endif
      endif ; enddo
      if (changed_col == 0) exit ! No more columns are unstable.

    enddo  ! End of iterations, all columns are now stable.

    ! Do the final return pass on the columns where the penultimate pass was downward.
    do i=is,ie ; do_i(i) = ((G%mask2dT(i,j) > 0.0) .and. last_down(i)) ; enddo
    do i=is,ie ; if (do_i(i)) then
      h_a = h(i,j,nz) + h_neglect
      b_a = 1.0 / (h_a + d_a(i,nz)*mix(i,nz))
      T_adj(i,j,nz) = b_a * (h_a*tv%T(i,j,nz) + mix(i,nz)*Te_a(i,nz-1))
      S_adj(i,j,nz) = b_a * (h_a*tv%S(i,j,nz) + mix(i,nz)*Se_a(i,nz-1))
    endif ; enddo
    do k=nz-1,1,-1 ; do i=is,ie ; if (do_i(i)) then
      T_adj(i,j,k) = Te_a(i,k) + c_a(i,K+1)*T_adj(i,j,k+1)
      S_adj(i,j,k) = Se_a(i,k) + c_a(i,K+1)*S_adj(i,j,k+1)
    endif ; enddo ; enddo

    do i=is,ie ; if (do_i(i)) then
      k = 1 ! A hook for debugging.
    endif ; enddo

    ! Do the final return pass on the columns where the penultimate pass was upward.
    ! Also do a simple copy of T & S values on land points.
    do i=is,ie
      do_i(i) = ((G%mask2dT(i,j) > 0.0) .and. .not.last_down(i))
      if (do_i(i)) then
        h_b = h(i,j,1) + h_neglect
        b_b = 1.0 / (h_b + d_b(i,2)*mix(i,2))
        T_adj(i,j,1) = b_b * (h_b*tv%T(i,j,1) + mix(i,2)*Te_b(i,2))
        S_adj(i,j,1) = b_b * (h_b*tv%S(i,j,1) + mix(i,2)*Se_b(i,2))
      elseif (G%mask2dT(i,j) <= 0.0) then
        T_adj(i,j,1) = tv%T(i,j,1) ; S_adj(i,j,1) = tv%S(i,j,1)
      endif
    enddo
    do k=2,nz ; do i=is,ie
      if (do_i(i)) then
        T_adj(i,j,k) = Te_b(i,k) + c_b(i,K)*T_adj(i,j,k-1)
        S_adj(i,j,k) = Se_b(i,k) + c_b(i,K)*S_adj(i,j,k-1)
      elseif (G%mask2dT(i,j) <= 0.0) then
        T_adj(i,j,k) = tv%T(i,j,k) ; S_adj(i,j,k) = tv%S(i,j,k)
      endif
    enddo ; enddo

    do i=is,ie ; if (do_i(i)) then
      k = 1 ! A hook for debugging.
    endif ; enddo

  enddo ! j-loop

  k = 1 ! A hook for debugging.

  ! The following set of expressions for the final values are derived from the partial
  ! updates for the estimated temperatures and salinities around an interface, then directly
  ! solving for the final temperatures and salinities.  They are here for later reference
  ! and to document an intermediate step in the stability calculation.
    ! hp_a = (h_a + d_a(i,K-1)*mix(i,K-1))
    ! hp_b = (h_b + d_b(i,K+1)*mix(i,K+1))
    ! b2_c = 1.0 / (hp_a*hp_b + (hp_a + hp_b) * mix(i,K))
    ! Th_a = h_a*tv%T(i,j,k-1) + mix(i,K-1)*Te_a(i,k-2)
    ! Th_b = h_b*tv%T(i,j,k)   + mix(i,K+1)*Te_b(i,k+1)
    ! T_fin(i,k)   = ( (hp_a + mix(i,K)) * Th_b  + Th_a * mix(i,K) ) * b2_c
    ! T_fin(i,k-1) = ( (hp_b + mix(i,K)) * Th_a  + Th_b * mix(i,K) ) * b2_c
    ! Sh_a = h_a*tv%S(i,j,k-1) + mix(i,K-1)*Se_a(i,k-2)
    ! Sh_b = h_b*tv%S(i,j,k)   + mix(i,K+1)*Se_b(i,k+1)
    ! S_fin(i,k)   = ( (hp_a + mix(i,K)) * Sh_b  + Sh_a * mix(i,K) ) * b2_c
    ! S_fin(i,k-1) = ( (hp_b + mix(i,K)) * Sh_a  + Sh_b * mix(i,K) ) * b2_c

end procedure full_convection
module procedure is_unstable
  is_unstable = (dRho_dT * ((h_a * h_b * (T_b - T_a) + &
                             mix_A*mix_B * (d_A*Te_bb - d_B*Te_aa)) + &
                            (h_a*mix_B * (Te_bb - d_B*T_a) + &
                             h_b*mix_A * (d_A*T_b - Te_aa)) ) + &
                 dRho_dS * ((h_a * h_b * (S_b - S_a) + &
                             mix_A*mix_B * (d_A*Se_bb - d_B*Se_aa)) + &
                            (h_a*mix_B * (Se_bb - d_B*S_a) + &
                             h_b*mix_A * (d_A*S_b - Se_aa)) ) < 0.0)
end procedure is_unstable
module procedure smoothed_dRdT_dRdS
  real :: mix(SZI_(G),SZK_(GV)+1)  ! The diffusive mixing length (kappa*dt)/dz
  real :: b1(SZI_(G))              ! A tridiagonal solver variable [H-1 ~> m-1 or m2 kg-1]
  real :: d1(SZI_(G))              ! A tridiagonal solver variable [nondim]
  real :: c1(SZI_(G),SZK_(GV))     ! A tridiagonal solver variable [nondim]
  real :: T_f(SZI_(G),SZK_(GV))    ! Filtered temperatures [C ~> degC]
  real :: S_f(SZI_(G),SZK_(GV))    ! Filtered salinities [S ~> ppt]
  real :: pres(SZI_(G))            ! Interface pressures [R L2 T-2 ~> Pa].
  real :: T_EOS(SZI_(G))           ! Filtered and vertically averaged temperatures [C ~> degC]
  real :: S_EOS(SZI_(G))           ! Filtered and vertically averaged salinities [S ~> ppt]
  real :: kap_dt_x2                ! The product of 2*kappa*dt [H Z ~> m2 or kg m-1].
  real :: dz_neglect, h0           ! A negligible vertical distances [Z ~> m]
  real :: h_neglect                ! A negligible thickness to allow for zero thicknesses
  real :: h_tr                     ! The thickness at tracer points, plus h_neglect [H ~> m or kg m-2].
  integer, dimension(2) :: EOSdom  ! The i-computational domain for the equation of state
  integer :: i, k, is, ie, nz
  is = G%isc-halo ; ie = G%iec+halo
  nz = GV%ke

  h_neglect = GV%H_subroundoff
  dz_neglect = GV%dz_subroundoff
  kap_dt_x2 = 2.0*Kddt

  if (Kddt <= 0.0) then
    do k=1,nz ; do i=is,ie
      T_f(i,k) = tv%T(i,j,k) ; S_f(i,k) = tv%S(i,j,k)
    enddo ; enddo
  else
    h0 = 1.0e-16*sqrt(GV%H_to_m*US%m_to_Z*Kddt) + dz_neglect
    do i=is,ie
      mix(i,2) = kap_dt_x2 / ((dz(i,1)+dz(i,2)) + h0)

      h_tr = h(i,j,1) + h_neglect
      b1(i) = 1.0 / (h_tr + mix(i,2))
      d1(i) = b1(i) * h(i,j,1)
      T_f(i,1) = (b1(i)*h_tr)*tv%T(i,j,1)
      S_f(i,1) = (b1(i)*h_tr)*tv%S(i,j,1)
    enddo
    do k=2,nz-1 ; do i=is,ie
      mix(i,K+1) = kap_dt_x2 / ((dz(i,k)+dz(i,k+1)) + h0)

      c1(i,k) = mix(i,K) * b1(i)
      h_tr = h(i,j,k) + h_neglect
      b1(i) = 1.0 / ((h_tr + d1(i)*mix(i,K)) + mix(i,K+1))
      d1(i) = b1(i) * (h_tr + d1(i)*mix(i,K))
      T_f(i,k) = b1(i) * (h_tr*tv%T(i,j,k) + mix(i,K)*T_f(i,k-1))
      S_f(i,k) = b1(i) * (h_tr*tv%S(i,j,k) + mix(i,K)*S_f(i,k-1))
    enddo ; enddo
    do i=is,ie
      c1(i,nz) = mix(i,nz) * b1(i)
      h_tr = h(i,j,nz) + h_neglect
      b1(i) = 1.0 / (h_tr + d1(i)*mix(i,nz))
      T_f(i,nz) = b1(i) * (h_tr*tv%T(i,j,nz) + mix(i,nz)*T_f(i,nz-1))
      S_f(i,nz) = b1(i) * (h_tr*tv%S(i,j,nz) + mix(i,nz)*S_f(i,nz-1))
    enddo
    do k=nz-1,1,-1 ; do i=is,ie
      T_f(i,k) = T_f(i,k) + c1(i,k+1)*T_f(i,k+1)
      S_f(i,k) = S_f(i,k) + c1(i,k+1)*S_f(i,k+1)
    enddo ; enddo
  endif

  if (associated(p_surf)) then
    do i=is,ie ; pres(i) = p_surf(i,j) ; enddo
  else
    do i=is,ie ; pres(i) = 0.0 ; enddo
  endif
  EOSdom(:) = EOS_domain(G%HI, halo)
  call calculate_density_derivs(T_f(:,1), S_f(:,1), pres, dR_dT(:,1), dR_dS(:,1), tv%eqn_of_state, EOSdom)
  do i=is,ie ; pres(i) = pres(i) + h(i,j,1)*(GV%H_to_RZ*GV%g_Earth) ; enddo
  do K=2,nz
    do i=is,ie
      T_EOS(i) = 0.5*(T_f(i,k-1) + T_f(i,k))
      S_EOS(i) = 0.5*(S_f(i,k-1) + S_f(i,k))
    enddo
    call calculate_density_derivs(T_EOS, S_EOS, pres, dR_dT(:,K), dR_dS(:,K), tv%eqn_of_state, EOSdom)
    do i=is,ie ; pres(i) = pres(i) + h(i,j,k)*(GV%H_to_RZ*GV%g_Earth) ; enddo
  enddo
  call calculate_density_derivs(T_f(:,nz), S_f(:,nz), pres, dR_dT(:,nz+1), dR_dS(:,nz+1), &
                                tv%eqn_of_state, EOSdom)

end procedure smoothed_dRdT_dRdS
end submodule MOM_full_convection_s
