submodule (PLM_functions) PLM_functions_s
  implicit none
contains
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
module procedure PLM_slope_cw
  real :: sigma_l, sigma_c, sigma_r ! Left, central and right slope estimates as
  real :: u_min, u_max ! Minimum and maximum value across cell [A]
  real :: h_cn ! Thickness of center cell [H]
  h_cn = h_c + h_neglect

  ! Side differences
  sigma_r = u_r - u_c
  sigma_l = u_c - u_l

  ! This is the second order slope given by equation 1.7 of
  ! Piecewise Parabolic Method, Colella and Woodward (1984),
  ! http://dx.doi.org/10.1016/0021-991(84)90143-8.
  ! For uniform resolution it simplifies to ( u_r - u_l )/2 .
  sigma_c = ( h_c / ( h_cn + ( h_l + h_r ) ) ) * ( &
                ( 2.*h_l + h_c ) / ( h_r + h_cn ) * sigma_r &
              + ( 2.*h_r + h_c ) / ( h_l + h_cn ) * sigma_l )

  ! Limit slope so that reconstructions are bounded by neighbors
  u_min = min( u_l, u_c, u_r )
  u_max = max( u_l, u_c, u_r )
  if ( (sigma_l * sigma_r) > 0.0 ) then
    ! This limits the slope so that the edge values are bounded by the
    ! two cell averages spanning the edge.
    PLM_slope_cw = sign( min( abs(sigma_c), 2.*min( u_c - u_min, u_max - u_c ) ), sigma_c )
  else
    ! Extrema in the mean values require a PCM reconstruction avoid generating
    ! larger extreme values.
    PLM_slope_cw = 0.0
  endif

  ! This block tests to see if roundoff causes edge values to be out of bounds
  if (u_c - 0.5*abs(PLM_slope_cw) < u_min .or.  u_c + 0.5*abs(PLM_slope_cw) > u_max) then
    PLM_slope_cw = PLM_slope_cw * ( 1. - epsilon(PLM_slope_cw) )
  endif

  ! An attempt to avoid inconsistency when the values become unrepresentable.
  ! ### The following 1.E-140 is dimensionally inconsistent. A newer version of
  ! PLM is progress that will avoid the need for such rounding.
  if (abs(PLM_slope_cw) < 1.E-140) PLM_slope_cw = 0.

end procedure PLM_slope_cw
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
module procedure PLM_extrapolate_slope
  real :: left_edge ! Left edge value [A]
  real :: hl, hc ! Left and central cell thicknesses [H]
  hl = h_l + h_neglect
  hc = h_c + h_neglect

  ! The h2 scheme is used to compute the left edge value
  left_edge = (u_l*hc + u_c*hl) / (hl + hc)

  PLM_extrapolate_slope = 2.0 * ( u_c - left_edge )

end procedure PLM_extrapolate_slope
module procedure PLM_reconstruction
  integer       :: k           ! loop index
  real          :: u_l, u_r    ! left and right cell averages [A]
  real          :: slope       ! retained PLM slope for a normalized cell width [A]
  real          :: e_r         ! The edge value in the neighboring cell [A]
  real          :: edge        ! The projected edge value in the cell [A]
  real          :: almost_one  ! A value that is slightly smaller than 1 [nondim]
  real, dimension(N) :: slp    ! The first guess at the normalized tracer slopes [A]
  real, dimension(N) :: mslp   ! The monotonized normalized tracer slopes [A]
  almost_one = 1. - epsilon(slope)

  ! Loop on interior cells
  do k = 2,N-1
    slp(k) = PLM_slope_wa(h(k-1), h(k), h(k+1), h_neglect, u(k-1), u(k), u(k+1))
  enddo ! end loop on interior cells

  ! Boundary cells use PCM. Extrapolation is handled after monotonization.
  slp(1) = 0.
  slp(N) = 0.

  ! This loop adjusts the slope so that edge values are monotonic.
  do K = 2, N-1
    mslp(k) = PLM_monotonized_slope( u(k-1), u(k), u(k+1), slp(k-1), slp(k), slp(k+1) )
  enddo ! end loop on interior cells
  mslp(1) = 0.
  mslp(N) = 0.

  ! Store and return edge values and polynomial coefficients.
  edge_values(1,1) = u(1)
  edge_values(1,2) = u(1)
  ppoly_coef(1,1) = u(1)
  ppoly_coef(1,2) = 0.
  do k = 2, N-1
    slope = mslp(k)
    u_l = u(k) - 0.5 * slope ! Left edge value of cell k
    u_r = u(k) + 0.5 * slope ! Right edge value of cell k

    edge_values(k,1) = u_l
    edge_values(k,2) = u_r
    ppoly_coef(k,1) = u_l
    ppoly_coef(k,2) = ( u_r - u_l )
    ! Check to see if this evaluation of the polynomial at x=1 would be
    ! monotonic w.r.t. the next cell's edge value. If not, scale back!
    edge = ppoly_coef(k,2) + ppoly_coef(k,1)
    e_r = u(k+1) - 0.5 * sign( mslp(k+1), slp(k+1) )
    if ( (edge-u(k))*(e_r-edge)<0.) then
      ppoly_coef(k,2) = ppoly_coef(k,2) * almost_one
    endif
  enddo
  edge_values(N,1) = u(N)
  edge_values(N,2) = u(N)
  ppoly_coef(N,1) = u(N)
  ppoly_coef(N,2) = 0.

end procedure PLM_reconstruction
module procedure PLM_boundary_extrapolation
  real    :: slope     ! retained PLM slope for a normalized cell width [A]
  slope = - PLM_extrapolate_slope( h(2), h(1), h_neglect, u(2), u(1) )

  edge_values(1,1) = u(1) - 0.5 * slope
  edge_values(1,2) = u(1) + 0.5 * slope

  ppoly_coef(1,1) = edge_values(1,1)
  ppoly_coef(1,2) = edge_values(1,2) - edge_values(1,1)

  ! Extrapolate from N-1 to N to estimate slope
  slope = PLM_extrapolate_slope( h(N-1), h(N), h_neglect, u(N-1), u(N) )

  edge_values(N,1) = u(N) - 0.5 * slope
  edge_values(N,2) = u(N) + 0.5 * slope

  ppoly_coef(N,1) = edge_values(N,1)
  ppoly_coef(N,2) = edge_values(N,2) - edge_values(N,1)

end procedure PLM_boundary_extrapolation
end submodule PLM_functions_s
