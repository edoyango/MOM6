submodule (P3M_functions) P3M_functions_s
  implicit none
contains
module procedure P3M_interpolation
  call P3M_limiter( N, h, u, edge_values, ppoly_S, ppoly_coef, h_neglect, &
                    answer_date=answer_date )

end procedure P3M_interpolation
module procedure P3M_limiter
  integer :: k            ! loop index
  logical :: monotonic    ! boolean indicating whether the cubic is monotonic
  real    :: u0_l, u0_r   ! edge values [A]
  real    :: u1_l, u1_r   ! edge slopes [A H-1]
  real    :: u_l, u_c, u_r        ! left, center and right cell averages [A]
  real    :: h_l, h_c, h_r        ! left, center and right cell widths [H]
  real    :: sigma_l, sigma_c, sigma_r  ! left, center and right van Leer slopes [A H-1]
  real    :: slope        ! retained PLM slope [A H-1]
  call bound_edge_values( N, h, u, edge_values, h_neglect, answer_date=answer_date )

  ! 2. Systematically average discontinuous edge values
  call average_discontinuous_edge_values( N, edge_values )


  ! 3. Loop on cells and do the following
  !     (a) Build cubic curve
  !     (b) Check if cubic curve is monotonic
  !     (c) If not, monotonize cubic curve and rebuild it
  do k = 1,N

    ! Get edge values, edge slopes and cell width
    u0_l = edge_values(k,1)
    u0_r = edge_values(k,2)
    u1_l = ppoly_S(k,1)
    u1_r = ppoly_S(k,2)

    ! Get cell widths and cell averages (boundary cells are assumed to
    ! be local extrema for the sake of slopes)
    u_c = u(k)
    h_c = h(k)

    if ( k == 1 ) then
      h_l = h(k)
      u_l = u(k)
    else
      h_l = h(k-1)
      u_l = u(k-1)
    endif

    if ( k == N ) then
      h_r = h(k)
      u_r = u(k)
    else
      h_r = h(k+1)
      u_r = u(k+1)
    endif

    ! Compute limited slope
    sigma_l = 2.0 * ( u_c - u_l ) / ( h_c + h_neglect )
    sigma_c = 2.0 * ( u_r - u_l ) / ( h_l + 2.0*h_c + h_r + h_neglect )
    sigma_r = 2.0 * ( u_r - u_c ) / ( h_c + h_neglect )

    if ( (sigma_l * sigma_r) > 0.0 ) then
      slope = sign( min(abs(sigma_l),abs(sigma_c),abs(sigma_r)), sigma_c )
    else
      slope = 0.0
    endif

    ! If the slopes are small, set them to zero to prevent asymmetric representation near extrema.
    if ( abs(u1_l*h_c) < epsilon(u_c)*abs(u_c) ) u1_l = 0.0
    if ( abs(u1_r*h_c) < epsilon(u_c)*abs(u_c) ) u1_r = 0.0

    ! The edge slopes are limited from above by the respective
    ! one-sided slopes
    if ( abs(u1_l) > abs(sigma_l) ) then
      u1_l = sigma_l
    endif

    if ( abs(u1_r) > abs(sigma_r) ) then
      u1_r = sigma_r
    endif

    ! Build cubic interpolant (compute the coefficients)
    call build_cubic_interpolant( h, k, edge_values, ppoly_S, ppoly_coef )

    ! Check whether cubic is monotonic
    monotonic = is_cubic_monotonic( ppoly_coef, k )

    ! If cubic is not monotonic, monotonize it by modifiying the
    ! edge slopes, store the new edge slopes and recompute the
    ! cubic coefficients
    if ( .not.monotonic ) then
      call monotonize_cubic( h_c, u0_l, u0_r, sigma_l, sigma_r, slope, u1_l, u1_r )
    endif

    ! Store edge slopes
    ppoly_S(k,1) = u1_l
    ppoly_S(k,2) = u1_r

    ! Recompute coefficients of cubic
    call build_cubic_interpolant( h, k, edge_values, ppoly_S, ppoly_coef )

  enddo ! loop on cells

end procedure P3M_limiter
module procedure P3M_boundary_extrapolation
  integer :: i0, i1
  logical :: monotonic    ! boolean indicating whether the cubic is monotonic
  real    :: u0, u1  ! Values of u in two adjacent cells [A]
  real    :: h0, h1  ! Values of h in two adjacent cells, plus a smal increment [H]
  real    :: b, c, d ! Temporary variables [A]
  real    :: u0_l, u0_r ! Left and right edge values [A]
  real    :: u1_l, u1_r ! Left and right edge slopes [A H-1]
  real    :: slope   ! The cell center slope [A H-1]
  real    :: hNeglect_edge ! Negligibly small thickness [H]
  hNeglect_edge = h_neglect ; if (present(h_neglect_edge)) hNeglect_edge = h_neglect_edge

  ! ----- Left boundary -----
  i0 = 1
  i1 = 2
  h0 = h(i0) + hNeglect_edge
  h1 = h(i1) + hNeglect_edge
  u0 = u(i0)
  u1 = u(i1)

  ! Compute the left edge slope in neighboring cell and express it in
  ! the global coordinate system
  b = ppoly_coef(i1,2)
  u1_r = b / h1     ! derivative evaluated at xi = 0.0, expressed w.r.t. x

  ! Limit the right slope by the PLM limited slope
  slope = 2.0 * ( u1 - u0 ) / ( h0 + h_neglect )
  if ( abs(u1_r) > abs(slope) ) then
    u1_r = slope
  endif

  ! The right edge value in the boundary cell is taken to be the left
  ! edge value in the neighboring cell
  u0_r = edge_values(i1,1)

  ! Given the right edge value and slope, we determine the left
  ! edge value and slope by computing the parabola as determined by
  ! the right edge value and slope and the boundary cell average
  u0_l = 3.0 * u0 + 0.5 * h0*u1_r - 2.0 * u0_r
  u1_l = ( - 6.0 * u0 - 2.0 * h0*u1_r + 6.0 * u0_r) / ( h0 + h_neglect )

  ! Check whether the edge values are monotonic. For example, if the left edge
  ! value is larger than the right edge value while the slope is positive, the
  ! edge values are inconsistent and we need to modify the left edge value
  if ( (u0_r-u0_l) * slope < 0.0 ) then
    u0_l = u0_r
    u1_l = 0.0
    u1_r = 0.0
  endif

  ! Store edge values and slope, build cubic and check monotonicity
  edge_values(i0,1) = u0_l
  edge_values(i0,2) = u0_r
  ppoly_S(i0,1) = u1_l
  ppoly_S(i0,2) = u1_r

  ! Store edge values and slope, build cubic and check monotonicity
  call build_cubic_interpolant( h, i0, edge_values, ppoly_S, ppoly_coef )
  monotonic = is_cubic_monotonic( ppoly_coef, i0 )

  if ( .not.monotonic ) then
    call monotonize_cubic( h0, u0_l, u0_r, 0.0, slope, slope, u1_l, u1_r )

    ! Rebuild cubic after monotonization
    ppoly_S(i0,1) = u1_l
    ppoly_S(i0,2) = u1_r
    call build_cubic_interpolant( h, i0, edge_values, ppoly_S, ppoly_coef )

  endif

  ! ----- Right boundary -----
  i0 = N-1
  i1 = N
  h0 = h(i0) + hNeglect_edge
  h1 = h(i1) + hNeglect_edge
  u0 = u(i0)
  u1 = u(i1)

  ! Compute the right edge slope in neighboring cell and express it in
  ! the global coordinate system
  b = ppoly_coef(i0,2)
  c = ppoly_coef(i0,3)
  d = ppoly_coef(i0,4)
  u1_l = (b + 2*c + 3*d) / ( h0 + h_neglect ) ! derivative evaluated at xi = 1.0

  ! Limit the left slope by the PLM limited slope
  slope = 2.0 * ( u1 - u0 ) / ( h1 + h_neglect )
  if ( abs(u1_l) > abs(slope) ) then
    u1_l = slope
  endif

  ! The left edge value in the boundary cell is taken to be the right
  ! edge value in the neighboring cell
  u0_l = edge_values(i0,2)

  ! Given the left edge value and slope, we determine the right
  ! edge value and slope by computing the parabola as determined by
  ! the left edge value and slope and the boundary cell average
  u0_r = 3.0 * u1 - 0.5 * h1*u1_l - 2.0 * u0_l
  u1_r = ( 6.0 * u1 - 2.0 * h1*u1_l - 6.0 * u0_l) / ( h1 + h_neglect )

  ! Check whether the edge values are monotonic. For example, if the right edge
  ! value is smaller than the left edge value while the slope is positive, the
  ! edge values are inconsistent and we need to modify the right edge value
  if ( (u0_r-u0_l) * slope < 0.0 ) then
    u0_r = u0_l
    u1_l = 0.0
    u1_r = 0.0
  endif

  ! Store edge values and slope, build cubic and check monotonicity
  edge_values(i1,1) = u0_l
  edge_values(i1,2) = u0_r
  ppoly_S(i1,1) = u1_l
  ppoly_S(i1,2) = u1_r

  call build_cubic_interpolant( h, i1, edge_values, ppoly_S, ppoly_coef )
  monotonic = is_cubic_monotonic( ppoly_coef, i1 )

  if ( .not.monotonic ) then
    call monotonize_cubic( h1, u0_l, u0_r, slope, 0.0, slope, u1_l, u1_r )

    ! Rebuild cubic after monotonization
    ppoly_S(i1,1) = u1_l
    ppoly_S(i1,2) = u1_r
    call build_cubic_interpolant( h, i1, edge_values, ppoly_S, ppoly_coef )

  endif

end procedure P3M_boundary_extrapolation
module procedure build_cubic_interpolant
  real          :: u0_l, u0_r       ! edge values [A]
  real          :: u1_l, u1_r       ! edge slopes times the cell width [A]
  real          :: h_c              ! cell width  [H]
  real          :: a0, a1, a2, a3   ! cubic coefficients [A]
  h_c = h(k)

  u0_l = edge_values(k,1)
  u0_r = edge_values(k,2)

  u1_l = ppoly_S(k,1) * h_c
  u1_r = ppoly_S(k,2) * h_c

  a0 = u0_l
  a1 = u1_l
  a2 = 3.0 * ( u0_r - u0_l ) - u1_r - 2.0 * u1_l
  a3 = u1_r + u1_l + 2.0 * ( u0_l - u0_r )

  ppoly_coef(k,1) = a0
  ppoly_coef(k,2) = a1
  ppoly_coef(k,3) = a2
  ppoly_coef(k,4) = a3

end procedure build_cubic_interpolant
module procedure is_cubic_monotonic
  real :: a, b, c   ! Coefficients of the first derivative of the cubic [A]
  a = ppoly_coef(k,2)
  b = 2.0 * ppoly_coef(k,3)
  c = 3.0 * ppoly_coef(k,4)

  ! Look for real roots of the quadratic derivative equation, c*x**2 + b*x + a = 0, in (0, 1)
  if (b*b - 4.0*a*c <= 0.0) then  ! The cubic is monotonic everywhere.
    is_cubic_monotonic = .true.
  elseif (a * (a + (b + c)) < 0.0) then ! The derivative changes sign between the endpoints of (0, 1)
    is_cubic_monotonic = .false.
  elseif (b * (b + 2.0*c) < 0.0) then ! The second derivative changes sign inside of (0, 1)
    is_cubic_monotonic = .false.
  else
    is_cubic_monotonic = .true.
  endif

end procedure is_cubic_monotonic
module procedure monotonize_cubic
  logical       :: found_ip
  logical       :: inflexion_l  ! bool telling if inflex. pt must be on left
  logical       :: inflexion_r  ! bool telling if inflex. pt must be on right
  real          :: a1, a2, a3   ! Temporary slopes times the cell width [A]
  real          :: u1_l_tmp     ! trial left edge slope [A H-1]
  real          :: u1_r_tmp     ! trial right edge slope [A H-1]
  real          :: xi_ip        ! location of inflexion point in cell coordinates (0,1) [nondim]
  real          :: slope_ip     ! slope at inflexion point times cell width [A]
  found_ip = .false.
  inflexion_l = .false.
  inflexion_r = .false.

  ! If the edge slopes are inconsistent w.r.t. the limited PLM slope,
  ! set them to zero
  if ( u1_l*slope <= 0.0 ) then
    u1_l = 0.0
  endif

  if ( u1_r*slope <= 0.0 ) then
    u1_r = 0.0
  endif

  ! Compute the location of the inflexion point, which is the root
  ! of the second derivative
  a1 = h * u1_l
  a2 = 3.0 * ( u0_r - u0_l ) - h*(u1_r + 2.0*u1_l)
  a3 = h*(u1_r + u1_l) + 2.0*(u0_l - u0_r)

  ! There is a possible root (and inflexion point) only if a3 is nonzero.
  ! When a3 is zero, the second derivative of the cubic is constant (the
  ! cubic degenerates into a parabola) and no inflexion point exists.
  if ( a3 /= 0.0 ) then
    ! Location of inflexion point
    xi_ip = - a2 / (3.0 * a3)

    ! If the inflexion point lies in [0,1], change boolean value
    if ( (xi_ip >= 0.0) .AND. (xi_ip <= 1.0) ) then
      found_ip = .true.
    endif
  endif

  ! When there is an inflexion point within [0,1], check the slope
  ! to see if it is consistent with the limited PLM slope. If not,
  ! decide on which side we want to collapse the inflexion point.
  ! If the inflexion point lies on one of the edges, the cubic is
  ! guaranteed to be monotonic
  if ( found_ip ) then
    slope_ip = a1 + 2.0*a2*xi_ip + 3.0*a3*xi_ip*xi_ip

    ! Check whether slope is consistent
    if ( slope_ip*slope < 0.0 ) then
      if ( abs(sigma_l) < abs(sigma_r)  ) then
        inflexion_l = .true.
      else
        inflexion_r = .true.
      endif
    endif
  endif ! found_ip

  ! At this point, if the cubic is not monotonic, we know where the
  ! inflexion point should lie. When the cubic is monotonic, both
  ! 'inflexion_l' and 'inflexion_r' are false and nothing is to be done.

  ! Move inflexion point on the left
  if ( inflexion_l ) then

    u1_l_tmp = 1.5*(u0_r-u0_l)/h - 0.5*u1_r
    u1_r_tmp = 3.0*(u0_r-u0_l)/h - 2.0*u1_l

    if ( (u1_l_tmp*slope < 0.0) .AND. (u1_r_tmp*slope < 0.0) ) then

      u1_l = 0.0
      u1_r = 3.0 * (u0_r - u0_l) / h

    elseif (u1_l_tmp*slope < 0.0) then

      u1_r = u1_r_tmp
      u1_l = 1.5*(u0_r - u0_l)/h - 0.5*u1_r

    elseif (u1_r_tmp*slope < 0.0) then

      u1_l = u1_l_tmp
      u1_r = 3.0*(u0_r - u0_l)/h - 2.0*u1_l

    else

      u1_l = u1_l_tmp
      u1_r = u1_r_tmp

    endif

  endif ! end treating case with inflexion point on the left

  ! Move inflexion point on the right
  if ( inflexion_r ) then

    u1_l_tmp = 3.0*(u0_r-u0_l)/h - 2.0*u1_r
    u1_r_tmp = 1.5*(u0_r-u0_l)/h - 0.5*u1_l

    if ( (u1_l_tmp*slope < 0.0) .AND. (u1_r_tmp*slope < 0.0) ) then

      u1_l = 3.0 * (u0_r - u0_l) / h
      u1_r = 0.0

    elseif (u1_l_tmp*slope < 0.0) then

      u1_r = u1_r_tmp
      u1_l = 3.0*(u0_r - u0_l)/h - 2.0*u1_r

    elseif (u1_r_tmp*slope < 0.0) then

      u1_l = u1_l_tmp
      u1_r = 1.5*(u0_r - u0_l)/h - 0.5*u1_l

    else

      u1_l = u1_l_tmp
      u1_r = u1_r_tmp

    endif

  endif ! end treating case with inflexion point on the right

  ! Zero out negligibly small slopes.
  if ( abs(u1_l*h) < epsilon(u0_l) * (abs(u0_l) + abs(u0_r)) ) u1_l = 0.0
  if ( abs(u1_r*h) < epsilon(u0_l) * (abs(u0_l) + abs(u0_r)) ) u1_r = 0.0

end procedure monotonize_cubic
end submodule P3M_functions_s
