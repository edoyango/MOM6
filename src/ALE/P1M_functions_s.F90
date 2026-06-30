submodule (P1M_functions) P1M_functions_s
  implicit none
contains
module procedure P1M_interpolation
  integer   :: k            ! loop index
  real      :: u0_l, u0_r   ! edge values (left and right) [A]
  call bound_edge_values( N, h, u, edge_values, h_neglect, answer_date=answer_date )

  ! Systematically average discontinuous edge values (routine found in
  ! 'edge_values.F90')
  call average_discontinuous_edge_values( N, edge_values )

  ! Loop on interior cells to build interpolants
  do k = 1,N

    u0_l = edge_values(k,1)
    u0_r = edge_values(k,2)

    ppoly_coef(k,1) = u0_l
    ppoly_coef(k,2) = u0_r - u0_l

  enddo ! end loop on interior cells

end procedure P1M_interpolation
module procedure P1M_boundary_extrapolation
  real          :: u0, u1               ! cell averages [A]
  real          :: h0, h1               ! corresponding cell widths [H]
  real          :: slope                ! retained PLM slope [A]
  real          :: u0_l, u0_r           ! edge values [A]
  h0 = h(1)
  h1 = h(2)

  u0 = u(1)
  u1 = u(2)

  ! The standard PLM slope is computed as a first estimate for the
  ! interpolation within the cell
  slope = 2.0 * ( u1 - u0 )

  ! The right edge value is then computed and we check whether this
  ! right edge value is consistent: it cannot be larger than the edge
  ! value in the neighboring cell if the data set is increasing.
  ! If the right value is found to too large, the slope is further limited
  ! by using the edge value in the neighboring cell.
  u0_r = u0 + 0.5 * slope

  if ( (u1 - u0) * (edge_values(2,1) - u0_r) < 0.0 ) then
    slope = 2.0 * ( edge_values(2,1) - u0 )
  endif

  ! Using the limited slope, the left edge value is reevaluated and
  ! the interpolant coefficients recomputed
  if ( h0 /= 0.0 ) then
    edge_values(1,1) = u0 - 0.5 * slope
  else
    edge_values(1,1) = u0
  endif

  ppoly_coef(1,1) = edge_values(1,1)
  ppoly_coef(1,2) = edge_values(1,2) - edge_values(1,1)

  ! ------------------------------------------
  ! Right edge value in the left boundary cell
  ! ------------------------------------------
  h0 = h(N-1)
  h1 = h(N)

  u0 = u(N-1)
  u1 = u(N)

  slope = 2.0 * ( u1 - u0 )

  u0_l = u1 - 0.5 * slope

  if ( (u1 - u0) * (u0_l - edge_values(N-1,2)) < 0.0 ) then
    slope = 2.0 * ( u1 - edge_values(N-1,2) )
  endif

  if ( h1 /= 0.0 ) then
    edge_values(N,2) = u1 + 0.5 * slope
  else
    edge_values(N,2) = u1
  endif

  ppoly_coef(N,1) = edge_values(N,1)
  ppoly_coef(N,2) = edge_values(N,2) - edge_values(N,1)

end procedure P1M_boundary_extrapolation
end submodule P1M_functions_s
