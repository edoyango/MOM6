submodule (Recon1d_EMPLM_CWK) Recon1d_EMPLM_CWK_s
  implicit none
contains
module procedure reconstruct
  real :: slp ! The PLM slopes (difference across cell) [A]
  real :: sigma_l, sigma_c, sigma_r ! Left, central and right slope estimates as
  real :: u_min, u_max ! Minimum and maximum value across cell [A]
  real :: u_l, u_r, u_c ! Left, right, and center values [A]
  real :: u_e(this%n+1) ! Average of edge values [A]
  integer :: k, n
  n = this%n

  call this%reconstruct_parent(h, u)

  this%ur(1) = this%ul(2)
  this%ul(1) = u(1) + ( u(1) - this%ur(1) )

  this%ul(n) = this%ur(n-1)
  this%ur(n) = u(n) + ( u(n) - this%ul(n) )

end procedure reconstruct
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
    um(k) = this%average(k, 0.25, 0.75) ! Average from x=0.25 to 0.75 in each cell
  enddo
  call test%real_arr(3, um, (/1.,3.,5./), 'Return interval average')

  call this%destroy()
  deallocate( um, ul, ur, ull, urr )

  allocate( um(4), ul(4), ur(4) )
  call this%init(4)

  ! These values lead to non-monotonic reconstuctions which are
  ! valid for transport problems but not always appropriate for
  ! remapping to arbitrary resolution grids.
  ! The O(h^2) slopes are -, 2, 2, - and the limited
  ! slopes are 0, 1, 1, 0 so the everywhere the reconstructions
  ! are bounded by neighbors but ur(2) and ul(3) are out-of-order.
  call this%reconstruct( (/1.,1.,1.,1./), (/0.,3.,4.,7./) )
  do k = 1, 4
    ul(k) = this%f(k, 0.)
    ur(k) = this%f(k, 1.)
  enddo
  call test%real_arr(4, ul, (/-2.5,2.5,3.5,4.5/), 'Evaluation on left edge')
  call test%real_arr(4, ur, (/2.5,3.5,4.5,9.5/), 'Evaluation on right edge')

  deallocate( um, ul, ur )

  unit_tests = test%summarize('EMPLM_CWK:unit_tests')

end procedure unit_tests
end submodule Recon1d_EMPLM_CWK_s
