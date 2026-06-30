submodule (Recon1d_MPLM_WA_poly) Recon1d_MPLM_WA_poly_s
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

  this%degree = 2
  allocate( this%poly_coef(n,2) )

end procedure init
module procedure reconstruct
  real :: slp(this%n) ! The PLM slopes (difference across cell) [A]
  real :: mslp(this%n) ! The monotonized PLM slopes [A]
  real :: e_r, edge ! Edge values [A]
  real :: almost_one  ! A value that is slightly smaller than 1 [nondim]
  integer :: k, n
  n = this%n

  ! Loop over all cells
  do k = 1, n
    this%u_mean(k) = u(k)
  enddo

  ! Loop on interior cells
  do k = 2, n-1
    slp(k) = PLM_slope_wa(h(k-1), h(k), h(k+1), this%h_neglect, u(k-1), u(k), u(k+1))
  enddo ! end loop on interior cells

  ! Boundary cells use PCM. Extrapolation is handled after monotonization.
  slp(1) = 0.
  slp(n) = 0.

  ! This loop adjusts the slope so that edge values are monotonic.
  do k = 2, n-1
    mslp(k) = PLM_monotonized_slope( u(k-1), u(k), u(k+1), slp(k-1), slp(k), slp(k+1) )
  enddo ! end loop on interior cells
  mslp(1) = 0.
  mslp(n) = 0.

  ! Store and return edge values and polynomial coefficients.
  almost_one = 1. - epsilon(e_r)
  this%ul(1) = u(1)
  this%ur(1) = u(1)
  this%poly_coef(1,1) = u(1)
  this%poly_coef(1,2) = 0.
  do k = 2, n-1
    this%ul(k) = u(k) - 0.5 * mslp(k) ! Left edge value of cell k
    this%ur(k) = u(k) + 0.5 * mslp(k) ! Right edge value of cell k

    this%poly_coef(k,1) = this%ul(k)
    this%poly_coef(k,2) = this%ur(k) - this%ul(k)
    ! Check to see if this evaluation of the polynomial at x=1 would be
    ! monotonic w.r.t. the next cell's edge value. If not, scale back!
    edge = this%poly_coef(k,2) + this%poly_coef(k,1)
    e_r = u(k+1) - 0.5 * sign( mslp(k+1), slp(k+1) )
    if ( (edge-u(k))*(e_r-edge)<0.) then
      this%poly_coef(k,2) = this%poly_coef(k,2) * almost_one
    endif
  enddo
  this%ul(n) = u(n)
  this%ur(n) = u(n)
  this%poly_coef(n,1) = u(n)
  this%poly_coef(n,2) = 0.

end procedure reconstruct
module procedure PLM_slope_wa
  real :: sigma_l, sigma_c, sigma_r ! Left, central and right slope estimates as
  real :: u_min, u_max ! Minimum and maximum value across cell [A]
  sigma_r = u_r - u_c
  sigma_l = u_c - u_l

  ! Quasi-second order difference
  sigma_c = 2.0 * ( u_r - u_l ) * ( h_c / ( h_l + 2.0*h_c + h_r + h_neglect) )

  ! Limit slope so that reconstructions are bounded by neighbors
  u_min = min( u_l, u_c, u_r )
  u_max = max( u_l, u_c, u_r )
  if ( (sigma_l * sigma_r) > 0.0 ) then
    ! This limits the slope so that the edge values are bounded by the
    ! two cell averages spanning the edge.
    PLM_slope_wa = sign( min( abs(sigma_c), 2.*min( u_c - u_min, u_max - u_c ) ), sigma_c )
  else
    ! Extrema in the mean values require a PCM reconstruction avoid generating
    ! larger extreme values.
    PLM_slope_wa = 0.0
  endif

  ! This block tests to see if roundoff causes edge values to be out of bounds
  if (u_c - 0.5*abs(PLM_slope_wa) < u_min .or.  u_c + 0.5*abs(PLM_slope_wa) > u_max) then
    PLM_slope_wa = PLM_slope_wa * ( 1. - epsilon(PLM_slope_wa) )
  endif

  ! An attempt to avoid inconsistency when the values become unrepresentable.
  ! ### The following 1.E-140 is dimensionally inconsistent. A newer version of
  ! PLM is progress that will avoid the need for such rounding.
  if (abs(PLM_slope_wa) < 1.E-140) PLM_slope_wa = 0.

end procedure PLM_slope_wa
module procedure PLM_monotonized_slope
  real :: e_r, e_l, edge ! Right, left and temporary edge values [A]
  real :: almost_two ! The number 2, almost [nondim]
  real :: slp ! Magnitude of PLM central slope [A]
  almost_two = 2. * ( 1. - epsilon(s_c) )

  ! Edge values of neighbors abutting this cell
  e_r = u_l + 0.5*s_l
  e_l = u_r - 0.5*s_r
  slp = abs(s_c)

  ! Check that left edge is between right edge of cell to the left and this cell mean
  edge = u_c - 0.5 * s_c
  if ( ( edge - e_r ) * ( u_c - edge ) < 0. ) then
    edge = 0.5 * ( edge + e_r )
    slp = min( slp, abs( edge - u_c ) * almost_two )
  endif

  ! Check that right edge is between left edge of cell to the right and this cell mean
  edge = u_c + 0.5 * s_c
  if ( ( edge - u_c ) * ( e_l - edge ) < 0. ) then
    edge = 0.5 * ( edge + e_l )
    slp = min( slp, abs( edge - u_c ) * almost_two )
  endif

  PLM_monotonized_slope = sign( slp, s_c )

end procedure PLM_monotonized_slope
module procedure average
  average = this%poly_coef(k,1) &
          + this%poly_coef(k,2) * 0.5 * ( xb + xa )

end procedure average
module procedure remap_to_sub_grid
  integer :: i_sub ! Index of sub-cell
  integer :: i0 ! Index into h0(1:n0), source column
  integer :: i_max ! Used to record which sub-cell is the largest contribution of a source cell
  real :: dh_max ! Used to record which sub-cell is the largest contribution of a source cell [H]
  real :: xa, xb ! Non-dimensional position within a source cell (0..1) [nondim]
  real :: dh ! The width of the sub-cell [H]
  real :: duh ! The total amount of accumulated stuff (u*h) [A H]
  real :: dh0_eff ! Running sum of source cell thickness [H]
  integer :: i0_last_thick_cell, n0
  real :: u0_min(this%n), u0_max(this%n) ! Min/max of u0 for each source cell [A]
  real :: ul, ur ! left/right edge values of cell i0
  n0 = this%n

  i0_last_thick_cell = 0
  do i0 = 1, n0
    ul = this%ul(i0)
    ur = this%ur(i0)
    u0_min(i0) = min(ul, ur)
    u0_max(i0) = max(ul, ur)
    if (h0(i0)>0.) i0_last_thick_cell = i0
  enddo

  ! Loop over each sub-cell to calculate average/integral values within each sub-cell.
  ! Uses: h_sub, isub_src, h0_eff
  ! Sets: u_sub, uh_sub
  xa = 0.
  dh0_eff = 0.
  u02_err = 0.
  do i_sub = 1, n0+n1

    ! Sub-cell thickness from loop above
    dh = h_sub(i_sub)

    ! Source cell
    i0 = isub_src(i_sub)

    ! Evaluate average and integral for sub-cell i_sub.
    ! Integral is over distance dh but expressed in terms of non-dimensional
    ! positions with source cell from xa to xb  (0 <= xa <= xb <= 1).
    dh0_eff = dh0_eff + dh ! Cumulative thickness within the source cell
    if (h0(i0)>0.) then
      xb = dh0_eff / h0(i0) ! This expression yields xa <= xb <= 1.0
      xb = min(1., xb) ! This is only needed when the total target column is wider than the source column
      u_sub(i_sub) = this%average( i0, xa, xb )
    else ! Vanished cell
      xb = 1.
      u_sub(i_sub) = u0(i0)
    endif
    uh_sub(i_sub) = dh * u_sub(i_sub)
    u_sub(i_sub) = max( u_sub(i_sub), u0_min(i0) )
    u_sub(i_sub) = min( u_sub(i_sub), u0_max(i0) )

    if (isub_src(i_sub+1) /= i0) then
      ! If the next sub-cell is in a different source cell, reset the position counters
      dh0_eff = 0.
      xa = 0.
    else
      xa = xb ! Next integral will start at end of last
    endif

  enddo
  i_sub = n0+n1+1
  ! Sub-cell thickness from loop above
  dh = h_sub(i_sub)
  ! Source cell
  i0 = isub_src(i_sub)

  ! Evaluate average and integral for sub-cell i_sub.
  ! Integral is over distance dh but expressed in terms of non-dimensional
  ! positions with source cell from xa to xb  (0 <= xa <= xb <= 1).
  dh0_eff = dh0_eff + dh ! Cumulative thickness within the source cell
  if (h0(i0)>0.) then
    xb = dh0_eff / h0(i0) ! This expression yields xa <= xb <= 1.0
    xb = min(1., xb) ! This is only needed when the total target column is wider than the source column
    u_sub(i_sub) = this%average( i0, xa, xb )
  else ! Vanished cell
    xb = 1.
    u_sub(i_sub) = u0(i0)
  endif
  u_sub(i_sub) = max( u_sub(i_sub), u0_min(i0) )
  u_sub(i_sub) = min( u_sub(i_sub), u0_max(i0) )
  uh_sub(i_sub) = dh * u_sub(i_sub)

  ! Loop over each source cell substituting the integral/average for the thickest sub-cell (within
  ! the source cell) with the residual of the source cell integral minus the other sub-cell integrals
  ! aka a genius algorithm for accurate conservation when remapping from Robert Hallberg (\@Hallberg-NOAA).
  ! Uses: i0_last_thick_cell, isrc_max, h_sub, isrc_start, isrc_end, uh_sub, u0, h0
  ! Updates: uh_sub
  do i0 = 1, i0_last_thick_cell
    i_max = isrc_max(i0)
    dh_max = h_sub(i_max)
    if (dh_max > 0.) then
      ! duh will be the sum of sub-cell integrals within the source cell except for the thickest sub-cell.
      duh = 0.
      do i_sub = isrc_start(i0), isrc_end(i0)
        if (i_sub /= i_max) duh = duh + uh_sub(i_sub)
      enddo
      uh_sub(i_max) = u0(i0)*h0(i0) - duh
      u02_err = u02_err + max( abs(uh_sub(i_max)), abs(u0(i0)*h0(i0)), abs(duh) )
    endif
  enddo

  ! This should not generally be used
  if (this%check) then
    if ( this%check_reconstruction(h0, u0) ) stop 912 ! A debugger is required to understand why this failed
  endif

end procedure remap_to_sub_grid
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

  ! Check bounding of left edges, w.r.t. this cell mean and the previous cell right edge
  ! Check order of u, ur, ul
  ! Note that in OM4 implementation, we were not consistent for top and bottom layers due
  ! extrapolation using cell means rather than edge values
  do K = 2, this%n-2
    if ( ( this%ur(k) - this%u_mean(k) ) * ( this%ul(k+1) - this%ur(k) ) < 0. ) check_reconstruction = .true.
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

  do k = 1, 3
    um(k) = this%average(k, 0.5, 0.75) ! Average from x=0.25 to 0.75 in each cell
  enddo
  call test%real_arr(3, um, (/1.,3.25,5./), 'Return interval average')

  unit_tests = test%summarize('MPLM_WA_poly:unit_tests')

end procedure unit_tests
end submodule Recon1d_MPLM_WA_poly_s
