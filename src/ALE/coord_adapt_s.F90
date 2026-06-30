submodule (coord_adapt) coord_adapt_s
#include <MOM_memory.h>
  implicit none
contains
module procedure init_coord_adapt
  if (associated(CS)) call MOM_error(FATAL, "init_coord_adapt: CS already associated")
  allocate(CS)
  allocate(CS%coordinateResolution(nk))

  CS%nk = nk
  CS%coordinateResolution(:) = coordinateResolution(:)

  ! Set real parameter default values
  CS%adaptTimeRatio = 1e-1 ! Nondim.
  CS%adaptAlpha     = 1.0  ! Nondim.
  CS%adaptZoom      = 200.0 * m_to_H    ! [H ~> m or kg m-2]
  CS%adaptZoomCoeff = 0.0  ! Nondim.
  CS%adaptBuoyCoeff = 0.0  ! Nondim.
  CS%adaptDrho0     = 0.5 * kg_m3_to_R  ! [R ~> kg m-3]

end procedure init_coord_adapt
module procedure end_coord_adapt
  if (.not. associated(CS)) return
  deallocate(CS%coordinateResolution)
  deallocate(CS)
end procedure end_coord_adapt
module procedure set_adapt_params
  if (.not. associated(CS)) call MOM_error(FATAL, "set_adapt_params: CS not associated")

  if (present(adaptTimeRatio)) CS%adaptTimeRatio = adaptTimeRatio
  if (present(adaptAlpha)) CS%adaptAlpha = adaptAlpha
  if (present(adaptZoom)) CS%adaptZoom = adaptZoom
  if (present(adaptZoomCoeff)) CS%adaptZoomCoeff = adaptZoomCoeff
  if (present(adaptBuoyCoeff)) CS%adaptBuoyCoeff = adaptBuoyCoeff
  if (present(adaptDrho0)) CS%adaptDrho0 = adaptDrho0
  if (present(adaptDoMin)) CS%adaptDoMin = adaptDoMin
end procedure set_adapt_params
module procedure build_adapt_column
  integer :: k, nz
  real :: h_up        ! The upwind source grid thickness based on the direction of the
  real :: b1          ! The inverse of the tridiagonal denominator [nondim]
  real :: b_denom_1   ! The leading term in the tridiagonal denominator [nondim]
  real :: d1          ! A term in the tridiagonal expressions [nondim]
  real :: depth       ! Depth in thickness units [H ~> m or kg m-2]
  real :: nominal_z   ! A nominal interface position in thickness units [H ~> m or kg m-2]
  real :: stretching  ! A stretching factor for the water column [nondim]
  real :: drdz  ! The vertical density gradient [R H-1 ~> kg m-4 or m-1]
  real, dimension(SZK_(GV)+1) :: alpha ! drho/dT [R C-1 ~> kg m-3 degC-1]
  real, dimension(SZK_(GV)+1) :: beta  ! drho/dS [R S-1 ~> kg m-3 ppt-1]
  real, dimension(SZK_(GV)+1) :: del2sigma ! Laplacian of in situ density times grid spacing [R ~> kg m-3]
  real, dimension(SZK_(GV)+1) :: dh_d2s ! Thickness change in response to del2sigma [H ~> m or kg m-2]
  real, dimension(SZK_(GV)) :: kGrid ! grid diffusivity on layers [nondim]
  real, dimension(SZK_(GV)) :: c1 ! A tridiagonal work array [nondim]
  nz = CS%nk

  ! set bottom and surface of zNext
  zNext(1) = 0.
  zNext(nz+1) = zInt(i,j,nz+1)

  ! local depth for scaling diffusivity
  depth = nom_depth_H(i,j)

  ! initialize del2sigma and the thickness change response to it zero
  del2sigma(:) = 0.0 ; dh_d2s(:) = 0.0

  ! calculate del-squared of neutral density by a
  ! stencilled finite difference
  ! TODO: this needs to be adjusted to account for vanished layers near topography

  ! up (j-1)
  if (G%mask2dT(i,j-1) > 0.0) then
    call calculate_density_derivs( &
         0.5 * (tInt(i,j,2:nz) + tInt(i,j-1,2:nz)), &
         0.5 * (sInt(i,j,2:nz) + sInt(i,j-1,2:nz)), &
         0.5 * (zInt(i,j,2:nz) + zInt(i,j-1,2:nz)) * (GV%H_to_RZ * GV%g_Earth), &
         alpha, beta, tv%eqn_of_state, (/2,nz/) )

    del2sigma(2:nz) = del2sigma(2:nz) + &
         (alpha(2:nz) * (tInt(i,j-1,2:nz) - tInt(i,j,2:nz)) + &
          beta(2:nz)  * (sInt(i,j-1,2:nz) - sInt(i,j,2:nz)))
  endif
  ! down (j+1)
  if (G%mask2dT(i,j+1) > 0.0) then
    call calculate_density_derivs( &
         0.5 * (tInt(i,j,2:nz) + tInt(i,j+1,2:nz)), &
         0.5 * (sInt(i,j,2:nz) + sInt(i,j+1,2:nz)), &
         0.5 * (zInt(i,j,2:nz) + zInt(i,j+1,2:nz)) * (GV%H_to_RZ * GV%g_Earth), &
         alpha, beta, tv%eqn_of_state, (/2,nz/) )

    del2sigma(2:nz) = del2sigma(2:nz) + &
         (alpha(2:nz) * (tInt(i,j+1,2:nz) - tInt(i,j,2:nz)) + &
          beta(2:nz)  * (sInt(i,j+1,2:nz) - sInt(i,j,2:nz)))
  endif
  ! left (i-1)
  if (G%mask2dT(i-1,j) > 0.0) then
    call calculate_density_derivs( &
         0.5 * (tInt(i,j,2:nz) + tInt(i-1,j,2:nz)), &
         0.5 * (sInt(i,j,2:nz) + sInt(i-1,j,2:nz)), &
         0.5 * (zInt(i,j,2:nz) + zInt(i-1,j,2:nz)) * (GV%H_to_RZ * GV%g_Earth), &
         alpha, beta, tv%eqn_of_state, (/2,nz/) )

    del2sigma(2:nz) = del2sigma(2:nz) + &
         (alpha(2:nz) * (tInt(i-1,j,2:nz) - tInt(i,j,2:nz)) + &
          beta(2:nz)  * (sInt(i-1,j,2:nz) - sInt(i,j,2:nz)))
  endif
  ! right (i+1)
  if (G%mask2dT(i+1,j) > 0.0) then
    call calculate_density_derivs( &
         0.5 * (tInt(i,j,2:nz) + tInt(i+1,j,2:nz)), &
         0.5 * (sInt(i,j,2:nz) + sInt(i+1,j,2:nz)), &
         0.5 * (zInt(i,j,2:nz) + zInt(i+1,j,2:nz)) * (GV%H_to_RZ * GV%g_Earth), &
         alpha, beta, tv%eqn_of_state, (/2,nz/) )

    del2sigma(2:nz) = del2sigma(2:nz) + &
         (alpha(2:nz) * (tInt(i+1,j,2:nz) - tInt(i,j,2:nz)) + &
          beta(2:nz)  * (sInt(i+1,j,2:nz) - sInt(i,j,2:nz)))
  endif

  ! at this point, del2sigma contains the local neutral density curvature at
  ! h-points, on interfaces
  ! we need to divide by drho/dz to give an interfacial displacement
  !
  ! a positive curvature means we're too light relative to adjacent columns,
  ! so del2sigma needs to be positive too (push the interface deeper)
  call calculate_density_derivs(tInt(i,j,:), sInt(i,j,:), zInt(i,j,:) * (GV%H_to_RZ * GV%g_Earth), &
       alpha, beta, tv%eqn_of_state, (/1,nz+1/) )
  do K = 2, nz
    ! TODO make lower bound here configurable
    dh_d2s(K) = del2sigma(K) * (0.5 * (h(i,j,k-1) + h(i,j,k))) / &
         max(alpha(K) * (tv%T(i,j,k) - tv%T(i,j,k-1)) + &
             beta(K)  * (tv%S(i,j,k) - tv%S(i,j,k-1)), 1e-20*US%kg_m3_to_R)

    ! don't move the interface so far that it would tangle with another
    ! interface in the direction we're moving (or exceed a Nyquist limit
    ! that could cause oscillations of the interface)
    h_up = merge(h(i,j,k), h(i,j,k-1), dh_d2s(K) > 0.)
    dh_d2s(K) = 0.5 * CS%adaptAlpha * &
         sign(min(abs(del2sigma(K)), 0.5 * h_up), dh_d2s(K))

    ! update interface positions so we can diffuse them
    zNext(K) = zInt(i,j,K) + dh_d2s(K)
  enddo

  ! solve diffusivity equation to smooth grid
  ! upper diagonal coefficients: -kGrid(2:nz)
  ! lower diagonal coefficients: -kGrid(1:nz-1)
  ! diagonal coefficients:       1 + (kGrid(1:nz-1) + kGrid(2:nz))
  !
  ! first, calculate the diffusivities within layers
  do k = 1, nz
    ! calculate the dr bit of drdz
    drdz = 0.5 * (alpha(K) + alpha(K+1)) * (tInt(i,j,K+1) - tInt(i,j,K)) + &
           0.5 * (beta(K)  + beta(K+1))  * (sInt(i,j,K+1) - sInt(i,j,K))
    ! divide by dz from the new interface positions
    drdz = drdz / (zNext(K) - zNext(K+1) + GV%H_subroundoff)
    ! don't do weird stuff in unstably-stratified regions
    drdz = max(drdz, 0.)

    ! set vertical grid diffusivity
    kGrid(k) = (CS%adaptTimeRatio * nz**2 * depth) * &
         ( CS%adaptZoomCoeff / (CS%adaptZoom + 0.5*(zNext(K) + zNext(K+1))) + &
           (CS%adaptBuoyCoeff * drdz / CS%adaptDrho0) + &
           max(1.0 - CS%adaptZoomCoeff - CS%adaptBuoyCoeff, 0.0) / depth)
  enddo

  ! initial denominator (first diagonal element)
  b1 = 1.0
  ! initial Q_1 = 1 - q_1 = 1 - 0/1
  d1 = 1.0
  ! work on all interior interfaces
  do K = 2, nz
    ! calculate numerator of Q_k
    b_denom_1 = 1. + d1 * kGrid(k-1)
    ! update denominator for k
    b1 = 1.0 / (b_denom_1 + kGrid(k))

    ! calculate q_k
    c1(K) = kGrid(k) * b1
    ! update Q_k = 1 - q_k
    d1 = b_denom_1 * b1

    ! update RHS
    zNext(K) = b1 * (zNext(K) + kGrid(k-1)*zNext(K-1))
  enddo
  ! final substitution
  do K = nz, 2, -1
    zNext(K) = zNext(K) + c1(K)*zNext(K+1)
  enddo

  if (CS%adaptDoMin) then
    nominal_z = 0.
    stretching = zInt(i,j,nz+1) / depth

    do k = 2, nz+1
      nominal_z = nominal_z + CS%coordinateResolution(k-1) * stretching
      ! take the deeper of the calculated and nominal positions
      zNext(K) = max(zNext(K), nominal_z)
      ! interface can't go below topography
      zNext(K) = min(zNext(K), zInt(i,j,nz+1))
    enddo
  endif
end procedure build_adapt_column
end submodule coord_adapt_s
