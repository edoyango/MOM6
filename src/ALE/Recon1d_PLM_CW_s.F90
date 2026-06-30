submodule (Recon1d_PLM_CW) Recon1d_PLM_CW_s
  implicit none
contains
module procedure init
  this%n = n

  allocate( this%u_mean(n) )
  allocate( this%ul(n) )
  allocate( this%ur(n) )

  this%h_neglect = tiny( this%u_mean(1) )
  if (present(h_neglect)) this%h_neglect = h_neglect
  this%check = .false.
  if (present(check)) this%check = check

end procedure init
module procedure reconstruct
  real :: slp ! The PLM slopes (difference across cell) [A]
  real :: sigma_l, sigma_c, sigma_r ! Left, central and right slope estimates as
  real :: u_min, u_max ! Minimum and maximum value across cell [A]
  real :: u_l, u_r, u_c ! Left, right, and center values [A]
  real :: h_l, h_c, h_r ! Thickness of left, center and right cells [H]
  real :: h_c0 ! Thickness of center with h_neglect added [H]
  integer :: k, n
  n = this%n

  ! Loop over all cells
  do k = 1, n
    this%u_mean(k) = u(k)
  enddo

  ! Boundary cells use PCM
  this%ul(1) = u(1)
  this%ur(1) = u(1)

  ! Loop over interior cells
  do k = 2, n-1
    u_l = u(k-1)
    u_c = u(k)
    u_r = u(k+1)

    ! Side differences
    sigma_r = u_r - u_c
    sigma_l = u_c - u_l

    h_l = h(k-1)
    h_c = h(k)
    h_r = h(k+1)
    ! Avoids division by zero
    h_c0 = h_c + this%h_neglect

    ! This is the second order slope given by equation 1.7 of
    ! Piecewise Parabolic Method, Colella and Woodward (1984),
    ! http://dx.doi.org/10.1016/0021-991(84)90143-8.
    ! For uniform resolution it simplifies to ( u_r - u_l )/2 .
    sigma_c = ( h_c / ( h_c0 + ( h_l + h_r ) ) ) * ( &
                  ( 2.*h_l + h_c ) / ( h_r + h_c0 ) * sigma_r &
                + ( 2.*h_r + h_c ) / ( h_l + h_c0 ) * sigma_l )

    ! Limit slope so that reconstructions are bounded by neighbors
    u_min = min( u_l, u_c, u_r )
    u_max = max( u_l, u_c, u_r )

    if ( (sigma_l * sigma_r) > 0.0 ) then
      ! This limits the slope so that the edge values are bounded by the two cell averages spanning the edge
      slp = sign( min( abs(sigma_c), 2.*min( u_c - u_min, u_max - u_c ) ), sigma_c )
    else
      ! Extrema in the mean values require a PCM reconstruction
      slp = 0.0
    endif

    ! Left edge
    u_min = min( u_c, u_l )
    u_max = max( u_c, u_l )
    u_l = u_c - 0.5 * slp
    this%ul(k) = max( min( u_l, u_max), u_min )

    ! Right edge
    u_min = min( u_c, u_r )
    u_max = max( u_c, u_r )
    u_r = u_c + 0.5 * slp
    this%ur(k) = max( min( u_r, u_max), u_min )
  enddo

  ! Boundary cells use PCM
  this%ul(n) = u(n)
  this%ur(n) = u(n)

end procedure reconstruct
module procedure f
  real :: xc ! Bounded version of x [nondim]
  real :: du ! Difference across cell [A]
  real :: u_a, u_b ! Two estimate of f [A]
  du = this%ur(k) - this%ul(k)
  xc = max( 0., min( 1., x ) )

  ! This expression for u_a can overshoot u_r but is good for x<<1
  u_a = this%ul(k) + du * xc
  ! This expression for u_b can overshoot u_l but is good for 1-x<<1
  u_b = this%ur(k) + du * ( xc - 1. )

  ! Since u_a and u_b are both bounded, this will perserve uniformity
  f = 0.5 * ( u_a + u_b )

end procedure f
module procedure dfdx
  dfdx = this%ur(k) - this%ul(k)

end procedure dfdx
module procedure x
  real :: slp ! Difference across cell [A]
  slp = this%ur(k) - this%ul(k)
  if ( abs(slp) > 0. ) then
    x = ( t - this%ul(k) ) / slp
    x = max( 0., min( x, 1. ) )
  else
    slp = this%ul(min(k+1,this%n)) - this%ur(max(k-1,1))
    if ( abs(slp) > 0. ) slp = sign(1., slp)
    x = 0.5 ! Fall back if t==u_mean
    ! if t>u_mean & slp=1 then x=1
    ! if t<u_mean & slp=1 then x=0
    ! if t>u_mean & slp=-1 then x=0
    ! if t<u_mean & slp=-1 then x=1
    ! if t=u_mean then slp=0
    ! if slp=0 then x=0.5
    if ( abs(t - this%u_mean(k)) > 0. ) x = 0.5 + slp * sign(0.5, t - this%u_mean(k))
  endif
end procedure x
module procedure average
  real :: xmab ! Mid-point between xa and xb (0 to 1)
  real :: u_a, u_b ! Values at xa and xb [A]
  xmab = 0.5 * ( xa + xb )

  ! The following expression is exact at xmab=0 and xmab=1,
  ! i.e. gives the numerically correct values.
  ! It is not obvious that the expression is monotonic but according to
  ! https://math.stackexchange.com/questions/907329/accurate-floating-point-linear-interpolation
  ! it will be for the default rounding behavior. Otherwise is it
  ! then possible this expression can be outside the range of ul and ur?
! average = this%ul(k) * ( 1. - xmab ) + this%ur(k) * xmab
  ! Emperically it fails the uniform value test

  ! The following is more complicated but seems to ensure being within bounds.
  ! This expression for u_a can overshoot u_r but is good for xmab<<1
  u_a = this%ul(k) + ( this%ur(k)  - this%ul(k) ) * xmab
  ! This expression for u_b can overshoot u_l but is good for 1-xmab<<1
  u_b = this%ur(k) + ( this%ul(k)  - this%ur(k) ) * ( 1. - xmab )
  ! Replace xmab with -1 for xmab<0.5, 1 for xmab>=0.5
! xmab = sign(1., xmab-0.5)
  ! Select either u_a or u_b, depending whether mid-point of xa, xb is smaller/larger than 0.5
! average = xmab * u_b + ( 1. - xmab ) * u_a

  ! Since u_a and u_b are both bounded, this will perserve uniformity but will the
  ! sum be bounded? Emperically it seems to work...
  average = 0.5 * ( u_a + u_b )

end procedure average
module procedure destroy
  deallocate( this%u_mean, this%ul, this%ur )

end procedure destroy
module procedure check_reconstruction
  integer :: k
  check_reconstruction = .false.

  do k = 1, this%n
    if ( abs( this%u_mean(k) - u(k) ) > 0. ) check_reconstruction = .true.
  enddo

  ! Check the cell reconstruction is monotonic within each cell (it should be as a straight line)
  do k = 1, this%n
    if ( ( this%u_mean(k) - this%ul(k) ) * ( this%ur(k) - this%u_mean(k) ) < 0. ) check_reconstruction = .true.
  enddo

  ! Check the cell is a straight line (to within machine precision)
  do k = 1, this%n
    if ( abs(2. * this%u_mean(k) - ( this%ul(k) + this%ur(k) )) > epsilon(this%u_mean(1)) * &
         max(abs(2. * this%u_mean(k)), abs(this%ul(k)), abs(this%ur(k))) ) check_reconstruction = .true.
  enddo

  ! Check bounding of right edges, w.r.t. the cell means
  do K = 1, this%n-1
    if ( ( this%ur(k) - this%u_mean(k) ) * ( this%u_mean(k+1) - this%ur(k) ) < 0. ) check_reconstruction = .true.
  enddo

  ! Check bounding of left edges, w.r.t. the cell means
  do K = 2, this%n
    if ( ( this%u_mean(k) - this%ul(k) ) * ( this%ul(k) - this%u_mean(k-1) ) < 0. ) check_reconstruction = .true.
  enddo

  ! PLM is not globally monotonic (expected)

! ! Check bounding of right edges, w.r.t. this cell mean and the next cell left edge
! do K = 1, this%n-1
!   if ( ( this%ur(k) - this%u_mean(k) ) * ( this%ul(k+1) - this%ur(k) ) < 0. ) check_reconstruction = .true.
! enddo

! ! Check bounding of left edges, w.r.t. this cell mean and the previous cell right edge
! do K = 2, this%n
!   if ( ( this%u_mean(k) - this%ul(k) ) * ( this%ul(k) - this%ur(k-1) ) < 0. ) check_reconstruction = .true.
! enddo

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
  call test%real_arr(3, ul, (/1.,2.,5./), 'Evaluation on left edge')
  call test%real_arr(3, um, (/1.,3.,5./), 'Evaluation in center')
  call test%real_arr(3, ur, (/1.,4.,5./), 'Evaluation on right edge')

  do k = 1, 3
    ul(k) = this%dfdx(k, 0.)
    um(k) = this%dfdx(k, 0.5)
    ur(k) = this%dfdx(k, 1.)
  enddo
  call test%real_arr(3, ul, (/0.,2.,0./), 'dfdx on left edge')
  call test%real_arr(3, um, (/0.,2.,0./), 'dfdx in center')
  call test%real_arr(3, ur, (/0.,2.,0./), 'dfdx on right edge')

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
    um(k) = this%average(k, 0.5, 0.75) ! Average from x=0.25 to 0.75 in each cell
  enddo
  call test%real_arr(3, um, (/1.,3.25,5./), 'Return interval average')

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
  call test%real_arr(4, ul, (/0.,2.,3.,7./), 'Evaluation on left edge')
  call test%real_arr(4, ur, (/0.,4.,5.,7./), 'Evaluation on right edge')

  deallocate( um, ul, ur )

  unit_tests = test%summarize('PLM_CW:unit_tests')

end procedure unit_tests
end submodule Recon1d_PLM_CW_s
