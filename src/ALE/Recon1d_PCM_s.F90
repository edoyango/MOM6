submodule (Recon1d_PCM) Recon1d_PCM_s
  implicit none
contains
module procedure init
  if (present(h_neglect)) this%n = n ! no-op to avoid compiler warning about unused dummy argument
  if (present(check)) this%check = check

  this%n = n

  allocate( this%u_mean(n) )

end procedure init
module procedure reconstruct
  integer :: k
  this%u_mean(1) = h(1) ! no-op to avoid compiler warning about unused dummy argument

  do k = 1, this%n
    this%u_mean(k) = u(k)
  enddo

end procedure reconstruct
module procedure f
  f = this%u_mean(k)

end procedure f
module procedure dfdx
  dfdx = 0.

end procedure dfdx
module procedure x
  real :: slp  ! Difference across cell [A]
  slp = this%u_mean(min(k+1,this%n)) - this%u_mean(max(k-1,1))
  if ( abs(slp) > 0. ) slp = sign(1., slp)
  x = 0.5 ! Fall back if t==u_mean
  ! if t>u_mean & slp=1 then x=1
  ! if t<u_mean & slp=1 then x=0
  ! if t>u_mean & slp=-1 then x=0
  ! if t<u_mean & slp=-1 then x=1
  ! if slp=0 then x=0.5
  if ( abs(t - this%u_mean(k)) > 0. ) x = 0.5 + slp * sign(0.5, t - this%u_mean(k))
end procedure x
module procedure average
  average = xb + xa ! no-op to avoid compiler warnings about unused dummy argument
  average = this%u_mean(k)

end procedure average
module procedure destroy
  deallocate( this%u_mean )

end procedure destroy
module procedure check_reconstruction
  integer :: k
  check_reconstruction = .false.

  do k = 1, this%n
    if ( abs( this%u_mean(k) - u(k) ) > 0. ) check_reconstruction = .true.
  enddo

end procedure check_reconstruction
module procedure unit_tests
  real, allocatable :: ul(:), ur(:), um(:) ! test values [A]
  type(testing) :: test ! convenience functions
  integer :: k
  call test%set( stdout=stdout ) ! Sets the stdout channel in test
  call test%set( stderr=stderr ) ! Sets the stderr channel in test
  call test%set( verbose=verbose ) ! Sets the verbosity flag in test

  call this%init(3)
  call test%test( this%n /= 3, 'Setting number of levels')
  allocate( um(3), ul(3), ur(3) )

  call this%reconstruct( (/2.,2.,2./), (/1.,3.,5./) )
  call test%real_arr(3, this%u_mean, (/1.,3.,5./), 'Setting cell values')

  do k = 1, 3
    ul(k) = this%f(k, 0.)
    um(k) = this%f(k, 0.5)
    ur(k) = this%f(k, 1.)
  enddo
  call test%real_arr(3, ul, (/1.,3.,5./), 'Evaluation on left edge')
  call test%real_arr(3, um, (/1.,3.,5./), 'Evaluation in center')
  call test%real_arr(3, ur, (/1.,3.,5./), 'Evaluation on right edge')

  do k = 1, 3
    ul(k) = this%dfdx(k, 0.)
    um(k) = this%dfdx(k, 0.5)
    ur(k) = this%dfdx(k, 1.)
  enddo
  call test%real_arr(3, ul, (/0.,0.,0./), 'dfdx on left edge')
  call test%real_arr(3, um, (/0.,0.,0./), 'dfdx in center')
  call test%real_arr(3, ur, (/0.,0.,0./), 'dfdx on right edge')

  call test%real_scalar( this%x(1,0.), 0., 'f-1(1,0)=0')
  call test%real_scalar( this%x(1,1.), 0.5, 'f-1(1,1)=0.5')
  call test%real_scalar( this%x(1,3.), 1., 'f-1(1,3)=1')
  call test%real_scalar( this%x(2,1.), 0., 'f-1(2,1)=0')
  call test%real_scalar( this%x(2,3.), 0.5, 'f-1(2,3)=0.5')
  call test%real_scalar( this%x(2,5.), 1., 'f-1(2,5)=1')
  call test%real_scalar( this%x(3,3.), 0., 'f-1(3,3)=0')
  call test%real_scalar( this%x(3,5.), 0.5, 'f-1(3,5)=0.5')
  call test%real_scalar( this%x(3,7.), 1., 'f-1(3,7)=1')

  do k = 1, 3
    um(k) = this%average(k, 0.5, 0.75)
  enddo
  call test%real_arr(3, um, (/1.,3.,5./), 'Return interval average')

  unit_tests = test%summarize('PCM:unit_tests')

end procedure unit_tests
end submodule Recon1d_PCM_s
