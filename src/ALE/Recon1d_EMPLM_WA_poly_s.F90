submodule (Recon1d_EMPLM_WA_poly) Recon1d_EMPLM_WA_poly_s
  implicit none
contains
module procedure reconstruct
  integer :: n
  real :: slope ! Difference of u across cell [A]
  call this%reconstruct_parent(h, u)

  n = this%n

  ! Fix reconstruction for first cell
  slope = - PLM_extrapolate_slope( h(2), h(1), this%h_neglect, u(2), u(1) )
  this%ul(1) = u(1) - 0.5 * slope
  this%ur(1) = u(1) + 0.5 * slope
  this%poly_coef(1,1) = this%ul(1)
  this%poly_coef(1,2) = this%ur(1) - this%ul(1)

  ! Fix reconstruction for last cell
  slope = PLM_extrapolate_slope( h(n-1), h(n), this%h_neglect, u(n-1), u(n) )
  this%ul(n) = u(n) - 0.5 * slope
  this%ur(n) = u(n) + 0.5 * slope
  this%poly_coef(n,1) = this%ul(n)
  this%poly_coef(n,2) = this%ur(n) - this%ul(n)

end procedure reconstruct
module procedure PLM_extrapolate_slope
  real :: left_edge ! Left edge value [A]
  real :: hl, hc ! Left and central cell thicknesses [H]
  hl = h_l + h_neglect
  hc = h_c + h_neglect

  ! The h2 scheme is used to compute the left edge value
  left_edge = (u_l*hc + u_c*hl) / (hl + hc)

  PLM_extrapolate_slope = 2.0 * ( u_c - left_edge )

end procedure PLM_extrapolate_slope
module procedure check_reconstruction
  integer :: k
  check_reconstruction = .false.

  do k = 1, this%n
    if ( abs( this%u_mean(k) - u(k) ) > 0. ) check_reconstruction = .true.
  enddo

  ! Check implied curvature
  do k = 1, this%n
    if ( ( this%u_mean(k) - this%ul(k) ) * ( this%ur(k) - this%u_mean(k) ) < 0. ) check_reconstruction = .true.
  enddo

  ! These two checks fail MOM_remapping:test_recon_consistency in the presence of vanished layers
  ! e.g. intel/2023.2.0 on gaea at iter=26

! ! Check bounding of right edges, w.r.t. the cell means
! do K = 1, this%n-1
!   if ( ( this%ur(k) - this%u_mean(k) ) * ( this%u_mean(k+1) - this%ur(k) ) < 0. ) check_reconstruction = .true.
! enddo

! ! Check bounding of left edges, w.r.t. the cell means
! do K = 2, this%n
!   if ( ( this%u_mean(k) - this%ul(k) ) * ( this%ul(k) - this%u_mean(k-1) ) < 0. ) check_reconstruction = .true.
! enddo

  ! Check order of u, ur, ul
  ! Note that in the OM4-era implementation, we were not consistent for top and bottom layers due
  ! extrapolation using cell means rather than edge values, hence reduced range for K
  do K = 2, this%n-2
    if ( ( this%ur(k) - this%u_mean(k) ) * ( this%ul(k+1) - this%ur(k) ) < 0. ) check_reconstruction = .true.
  enddo

  ! Check bounding of left edges, w.r.t. this cell mean and the previous cell right edge
  do K = 3, this%n-1
    if ( ( this%u_mean(k) - this%ul(k) ) * ( this%ul(k) - this%ur(k-1) ) < 0. ) check_reconstruction = .true.
  enddo

end procedure check_reconstruction
module procedure unit_tests
  real, allocatable :: ul(:), ur(:), um(:) ! test values [A]
  real, allocatable :: ull(:), urr(:) ! test values [A]
  type(testing) :: test ! convenience functions
  integer :: k
  call test%set( stdout=stdout ) ! Sets the stdout channel in test
  call test%set( stderr=stderr ) ! Sets the stderr channel in test
  call test%set( verbose=verbose ) ! Sets the verbosity flag in test

  call this%init(3)
  call test%test( this%n /= 3, 'Setting number of levels')
  allocate( um(3), ul(3), ur(3), ull(3), urr(3) )

  call this%reconstruct( (/2.,2.,2./), (/1.,3.,5./) )
  call test%real_arr(3, this%u_mean, (/1.,3.,5./), 'Setting cell values')

  do k = 1, 3
    ul(k) = this%f(k, 0.)
    um(k) = this%f(k, 0.5)
    ur(k) = this%f(k, 1.)
  enddo
  call test%real_arr(3, ul, (/0.,2.,4./), 'Evaluation on left edge')
  call test%real_arr(3, um, (/1.,3.,5./), 'Evaluation in center')
  call test%real_arr(3, ur, (/2.,4.,6./), 'Evaluation on right edge')

  do k = 1, 3
    ul(k) = this%dfdx(k, 0.)
    um(k) = this%dfdx(k, 0.5)
    ur(k) = this%dfdx(k, 1.)
  enddo
  call test%real_arr(3, ul, (/2.,2.,2./), 'dfdx on left edge')
  call test%real_arr(3, um, (/2.,2.,2./), 'dfdx in center')
  call test%real_arr(3, ur, (/2.,2.,2./), 'dfdx on right edge')

  do k = 1, 3
    um(k) = this%average(k, 0.5, 0.75) ! Average from x=0.25 to 0.75 in each cell
  enddo
  call test%real_arr(3, um, (/1.25,3.25,5.25/), 'Return interval average')

  unit_tests = test%summarize('EMPLM_WA_poly:unit_tests')

end procedure unit_tests
end submodule Recon1d_EMPLM_WA_poly_s
