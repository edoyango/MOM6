submodule (regrid_interp) regrid_interp_s
  implicit none
contains
module procedure regridding_set_ppolys
  real :: h_neg_edge  ! A negligibly small width for the purpose of edge value
  logical :: extrapolate
  h_neg_edge = h_neglect ; if (present(h_neglect_edge)) h_neg_edge = h_neglect_edge

  ! Reset piecewise polynomials
  ppoly0_E(:,:) = 0.0
  ppoly0_S(:,:) = 0.0
  ppoly0_coefs(:,:) = 0.0

  extrapolate = CS%boundary_extrapolation

  ! Compute the interpolated profile of the density field and build grid
  select case (CS%interpolation_scheme)

    case ( INTERPOLATION_P1M_H2 )
      degree = DEGREE_1
      call edge_values_explicit_h2( n0, h0, densities, ppoly0_E )
      call P1M_interpolation( n0, h0, densities, ppoly0_E, ppoly0_coefs, h_neglect, answer_date=CS%answer_date )
      if (extrapolate) then
        call P1M_boundary_extrapolation( n0, h0, densities, ppoly0_E, ppoly0_coefs )
      endif

    case ( INTERPOLATION_P1M_H4 )
      degree = DEGREE_1
      if ( n0 >= 4 ) then
        call edge_values_explicit_h4( n0, h0, densities, ppoly0_E, h_neg_edge, answer_date=CS%answer_date )
      else
        call edge_values_explicit_h2( n0, h0, densities, ppoly0_E )
      endif
      call P1M_interpolation( n0, h0, densities, ppoly0_E, ppoly0_coefs, h_neglect, answer_date=CS%answer_date )
      if (extrapolate) then
        call P1M_boundary_extrapolation( n0, h0, densities, ppoly0_E, ppoly0_coefs )
      endif

    case ( INTERPOLATION_P1M_IH4 )
      degree = DEGREE_1
      if ( n0 >= 4 ) then
        call edge_values_implicit_h4( n0, h0, densities, ppoly0_E, h_neg_edge, answer_date=CS%answer_date )
      else
        call edge_values_explicit_h2( n0, h0, densities, ppoly0_E )
      endif
      call P1M_interpolation( n0, h0, densities, ppoly0_E, ppoly0_coefs, h_neglect, answer_date=CS%answer_date )
      if (extrapolate) then
        call P1M_boundary_extrapolation( n0, h0, densities, ppoly0_E, ppoly0_coefs )
      endif

    case ( INTERPOLATION_PLM )
      degree = DEGREE_1
      call PLM_reconstruction( n0, h0, densities, ppoly0_E, ppoly0_coefs, h_neglect )
      if (extrapolate) then
        call PLM_boundary_extrapolation( n0, h0, densities, ppoly0_E, ppoly0_coefs, h_neglect )
      endif

    case ( INTERPOLATION_PPM_CW )
      if ( n0 >= 4 ) then
        degree = DEGREE_2
        call edge_values_explicit_h4cw( n0, h0, densities, ppoly0_E, h_neg_edge )
        call PPM_monotonicity(   n0,     densities, ppoly0_E )
        call PPM_reconstruction( n0, h0, densities, ppoly0_E, ppoly0_coefs, h_neglect, answer_date=CS%answer_date )
        if (extrapolate) then
          call PPM_boundary_extrapolation( n0, h0, densities, ppoly0_E, &
                                           ppoly0_coefs, h_neglect )
        endif
      else
        degree = DEGREE_1
        call edge_values_explicit_h2( n0, h0, densities, ppoly0_E )
        call P1M_interpolation( n0, h0, densities, ppoly0_E, ppoly0_coefs, h_neglect, answer_date=CS%answer_date )
        if (extrapolate) then
          call P1M_boundary_extrapolation( n0, h0, densities, ppoly0_E, ppoly0_coefs )
        endif
      endif

    case ( INTERPOLATION_PPM_H4 )
      if ( n0 >= 4 ) then
        degree = DEGREE_2
        call edge_values_explicit_h4( n0, h0, densities, ppoly0_E, h_neg_edge, answer_date=CS%answer_date )
        call PPM_reconstruction( n0, h0, densities, ppoly0_E, ppoly0_coefs, h_neglect, answer_date=CS%answer_date )
        if (extrapolate) then
          call PPM_boundary_extrapolation( n0, h0, densities, ppoly0_E, &
                                           ppoly0_coefs, h_neglect )
        endif
      else
        degree = DEGREE_1
        call edge_values_explicit_h2( n0, h0, densities, ppoly0_E )
        call P1M_interpolation( n0, h0, densities, ppoly0_E, ppoly0_coefs, h_neglect, answer_date=CS%answer_date )
        if (extrapolate) then
          call P1M_boundary_extrapolation( n0, h0, densities, ppoly0_E, ppoly0_coefs )
        endif
      endif

    case ( INTERPOLATION_PPM_IH4 )
      if ( n0 >= 4 ) then
        degree = DEGREE_2
        call edge_values_implicit_h4( n0, h0, densities, ppoly0_E, h_neg_edge, answer_date=CS%answer_date )
        call PPM_reconstruction( n0, h0, densities, ppoly0_E, ppoly0_coefs, h_neglect, answer_date=CS%answer_date )
        if (extrapolate) then
          call PPM_boundary_extrapolation( n0, h0, densities, ppoly0_E, &
                                           ppoly0_coefs, h_neglect )
        endif
      else
        degree = DEGREE_1
        call edge_values_explicit_h2( n0, h0, densities, ppoly0_E )
        call P1M_interpolation( n0, h0, densities, ppoly0_E, ppoly0_coefs, h_neglect, answer_date=CS%answer_date )
        if (extrapolate) then
          call P1M_boundary_extrapolation( n0, h0, densities, ppoly0_E, ppoly0_coefs )
        endif
      endif

    case ( INTERPOLATION_P3M_IH4IH3 )
      if ( n0 >= 4 ) then
        degree = DEGREE_3
        call edge_values_implicit_h4( n0, h0, densities, ppoly0_E, h_neg_edge, answer_date=CS%answer_date )
        call edge_slopes_implicit_h3( n0, h0, densities, ppoly0_S, h_neglect, answer_date=CS%answer_date )
        call P3M_interpolation( n0, h0, densities, ppoly0_E, ppoly0_S, &
                                ppoly0_coefs, h_neglect, answer_date=CS%answer_date )
        if (extrapolate) then
          call P3M_boundary_extrapolation( n0, h0, densities, ppoly0_E, ppoly0_S, &
                                           ppoly0_coefs, h_neglect, h_neg_edge )
        endif
      else
        degree = DEGREE_1
        call edge_values_explicit_h2( n0, h0, densities, ppoly0_E )
        call P1M_interpolation( n0, h0, densities, ppoly0_E, ppoly0_coefs, h_neglect, answer_date=CS%answer_date )
        if (extrapolate) then
          call P1M_boundary_extrapolation( n0, h0, densities, ppoly0_E, ppoly0_coefs )
        endif
      endif

    case ( INTERPOLATION_P3M_IH6IH5 )
      if ( n0 >= 6 ) then
        degree = DEGREE_3
        call edge_values_implicit_h6( n0, h0, densities, ppoly0_E, h_neg_edge, answer_date=CS%answer_date )
        call edge_slopes_implicit_h5( n0, h0, densities, ppoly0_S, h_neglect, answer_date=CS%answer_date )
        call P3M_interpolation( n0, h0, densities, ppoly0_E, ppoly0_S, &
                                ppoly0_coefs, h_neglect, answer_date=CS%answer_date )
        if (extrapolate) then
          call P3M_boundary_extrapolation( n0, h0, densities, ppoly0_E, ppoly0_S, &
                   ppoly0_coefs, h_neglect, h_neglect_edge )
        endif
      else
        degree = DEGREE_1
        call edge_values_explicit_h2( n0, h0, densities, ppoly0_E )
        call P1M_interpolation( n0, h0, densities, ppoly0_E, ppoly0_coefs, h_neglect, answer_date=CS%answer_date )
        if (extrapolate) then
          call P1M_boundary_extrapolation( n0, h0, densities, ppoly0_E, ppoly0_coefs )
        endif
      endif

    case ( INTERPOLATION_PQM_IH4IH3 )
      if ( n0 >= 4 ) then
        degree = DEGREE_4
        call edge_values_implicit_h4( n0, h0, densities, ppoly0_E, h_neg_edge, answer_date=CS%answer_date )
        call edge_slopes_implicit_h3( n0, h0, densities, ppoly0_S, h_neglect, answer_date=CS%answer_date )
        call PQM_reconstruction( n0, h0, densities, ppoly0_E, ppoly0_S, &
                                 ppoly0_coefs, h_neglect, answer_date=CS%answer_date )
        if (extrapolate) then
          call PQM_boundary_extrapolation_v1( n0, h0, densities, ppoly0_E, ppoly0_S, &
                                 ppoly0_coefs, h_neglect )
        endif
      else
        degree = DEGREE_1
        call edge_values_explicit_h2( n0, h0, densities, ppoly0_E )
        call P1M_interpolation( n0, h0, densities, ppoly0_E, ppoly0_coefs, h_neglect, answer_date=CS%answer_date )
        if (extrapolate) then
          call P1M_boundary_extrapolation( n0, h0, densities, ppoly0_E, ppoly0_coefs )
        endif
      endif

    case ( INTERPOLATION_PQM_IH6IH5 )
      if ( n0 >= 6 ) then
        degree = DEGREE_4
        call edge_values_implicit_h6( n0, h0, densities, ppoly0_E, h_neg_edge, answer_date=CS%answer_date )
        call edge_slopes_implicit_h5( n0, h0, densities, ppoly0_S, h_neglect, answer_date=CS%answer_date )
        call PQM_reconstruction( n0, h0, densities, ppoly0_E, ppoly0_S, &
                                 ppoly0_coefs, h_neglect, answer_date=CS%answer_date )
        if (extrapolate) then
          call PQM_boundary_extrapolation_v1( n0, h0, densities, ppoly0_E, ppoly0_S, &
                                 ppoly0_coefs, h_neglect )
        endif
      else
        degree = DEGREE_1
        call edge_values_explicit_h2( n0, h0, densities, ppoly0_E )
        call P1M_interpolation( n0, h0, densities, ppoly0_E, ppoly0_coefs, h_neglect, answer_date=CS%answer_date )
        if (extrapolate) then
          call P1M_boundary_extrapolation( n0, h0, densities, ppoly0_E, ppoly0_coefs )
        endif
      endif
  end select

end procedure regridding_set_ppolys
module procedure interpolate_grid
  integer        :: k ! loop index
  real           :: t ! current interface target density [A]
  x1(1) = x0(1)
  x1(n1+1) = x0(n0+1)

  ! Find coordinates for interior target values
  do k = 2,n1
    t = target_values(k)
    x1(k) = get_polynomial_coordinate ( n0, h0, x0, ppoly0_E, ppoly0_coefs, t, degree, &
                                        answer_date=answer_date )
    h1(k-1) = x1(k) - x1(k-1)
  enddo
  h1(n1) = x1(n1+1) - x1(n1)

end procedure interpolate_grid
module procedure build_and_interpolate_grid
  real, dimension(n0,2) :: ppoly0_E   ! Polynomial edge values [R ~> kg m-3]
  real, dimension(n0,2) :: ppoly0_S   ! Polynomial edge slopes [R H-1 ~> kg m-4 or m-1] or [R Z-1 ~> kg m-4]
  real, dimension(n0,DEGREE_MAX+1) :: ppoly0_C  ! Polynomial interpolant coeficients on the local 0-1 grid [R ~> kg m-3]
  integer :: degree
  call regridding_set_ppolys(CS, densities, n0, h0, ppoly0_E, ppoly0_S, ppoly0_C, &
       degree, h_neglect, h_neglect_edge)
  call interpolate_grid(n0, h0, x0, ppoly0_E, ppoly0_C, target_values, degree, &
       n1, h1, x1, answer_date=CS%answer_date)
end procedure build_and_interpolate_grid
module procedure get_polynomial_coordinate
  real                        :: xi0         ! normalized target coordinate [nondim]
  real, dimension(DEGREE_MAX) :: a           ! polynomial coefficients [A]
  real                        :: numerator   ! The numerator of an expression [A]
  real                        :: denominator ! The denominator of an expression [A]
  real                        :: delta       ! Newton-Raphson increment [nondim]
  real                        :: eps         ! offset used to get away from boundaries [nondim]
  real                        :: grad        ! gradient during N-R iterations [A]
  integer :: i, k, iter  ! loop indices
  integer :: k_found     ! index of target cell
  character(len=320) :: mesg
  logical :: use_2018_answers  ! If true use older, less accurate expressions.
  eps = NR_OFFSET
  k_found = -1
  use_2018_answers = (answer_date < 20190101)

  ! If the target value is outside the range of all values, we
  ! force the target coordinate to be equal to the lowest or
  ! largest value, depending on which bound is overtaken
  if ( target_value <= edge_values(1,1) ) then
    x_tgt = x_g(1)
    return  ! return because there is no need to look further
  endif

  ! Since discontinuous edge values are allowed, we check whether the target
  ! value lies between two discontinuous edge values at interior interfaces
  do k = 2,N
    if ( ( target_value >= edge_values(k-1,2) ) .AND. ( target_value <= edge_values(k,1) ) ) then
      x_tgt = x_g(k)
      return   ! return because there is no need to look further
    endif
  enddo

  ! If the target value is outside the range of all values, we
  ! force the target coordinate to be equal to the lowest or
  ! largest value, depending on which bound is overtaken
  if ( target_value >= edge_values(N,2) ) then
    x_tgt = x_g(N+1)
    return  ! return because there is no need to look further
  endif

  ! At this point, we know that the target value is bounded and does not
  ! lie between discontinuous, monotonic edge values. Therefore,
  ! there is a unique solution. We loop on all cells and find which one
  ! contains the target value. The variable k_found holds the index value
  ! of the cell where the taregt value lies.
  do k = 1,N
    if ( ( target_value > edge_values(k,1) ) .AND. ( target_value < edge_values(k,2) ) ) then
      k_found = k
      exit
    endif
  enddo

  ! At this point, 'k_found' should be strictly positive. If not, this is
  ! a major failure because it means we could not find any target cell
  ! despite the fact that the target value lies between the extremes. It
  ! means there is a major problem with the interpolant. This needs to be
  ! reported.
  if ( k_found == -1 ) then
    write(mesg,*) 'Could not find target coordinate', target_value, 'in get_polynomial_coordinate. This is '//&
                  'caused by an inconsistent interpolant (perhaps not monotonically increasing):', &
                  target_value, edge_values(1,1), edge_values(N,2)
    call MOM_error( FATAL, mesg )
  endif

  ! Reset all polynomial coefficients to 0 and copy those pertaining to
  ! the found cell
  a(:) = 0.0
  do i = 1,degree+1
    a(i) = ppoly_coefs(k_found,i)
  enddo

  ! Guess the middle of the cell to start Newton-Raphson iterations
  xi0 = 0.5

  ! Newton-Raphson iterations
  do iter = 1,NR_ITERATIONS

    if (use_2018_answers) then
      numerator = a(1) + a(2)*xi0 + a(3)*xi0*xi0 + a(4)*xi0*xi0*xi0 + &
                  a(5)*xi0*xi0*xi0*xi0 - target_value
      denominator = a(2) + 2*a(3)*xi0 + 3*a(4)*xi0*xi0 + 4*a(5)*xi0*xi0*xi0
    else  ! These expressions are mathematicaly equivalent but more accurate.
      numerator = (a(1) - target_value) + xi0*(a(2) + xi0*(a(3) + xi0*(a(4) + a(5)*xi0)))
      denominator = a(2) + xi0*(2.*a(3) + xi0*(3.*a(4) + 4.*a(5)*xi0))
    endif

    delta = -numerator / denominator

    xi0 = xi0 + delta

    ! Check whether new estimate is out of bounds. If the new estimate is
    ! indeed out of bounds, we manually set it to be equal to the overtaken
    ! bound with a small offset towards the interior when the gradient of
    ! the function at the boundary is zero (in which case, the Newton-Raphson
    ! algorithm does not converge).
    if ( xi0 < 0.0 ) then
      xi0 = 0.0
      grad = a(2)
      if ( grad == 0.0 ) xi0 = xi0 + eps
    endif

    if ( xi0 > 1.0 ) then
      xi0 = 1.0
      if (use_2018_answers) then
        grad = a(2) + 2*a(3) + 3*a(4) + 4*a(5)
      else  ! These expressions are mathematicaly equivalent but more accurate.
        grad = a(2) + (2.*a(3) + (3.*a(4) + 4.*a(5)))
      endif
      if ( grad == 0.0 ) xi0 = xi0 - eps
    endif

    ! break if converged or too many iterations taken
    if ( abs(delta) < NR_TOLERANCE ) exit
  enddo ! end Newton-Raphson iterations

  x_tgt = x_g(k_found) + xi0 * h(k_found)
end procedure get_polynomial_coordinate
module procedure interpolation_scheme
  select case ( uppercase(trim(interp_scheme)) )
    case ("P1M_H2");     interpolation_scheme = INTERPOLATION_P1M_H2
    case ("P1M_H4");     interpolation_scheme = INTERPOLATION_P1M_H4
    case ("P1M_IH2");    interpolation_scheme = INTERPOLATION_P1M_IH4
    case ("PLM");        interpolation_scheme = INTERPOLATION_PLM
    case ("PPM_CW");     interpolation_scheme = INTERPOLATION_PPM_CW
    case ("PPM_H4");     interpolation_scheme = INTERPOLATION_PPM_H4
    case ("PPM_IH4");    interpolation_scheme = INTERPOLATION_PPM_IH4
    case ("P3M_IH4IH3"); interpolation_scheme = INTERPOLATION_P3M_IH4IH3
    case ("P3M_IH6IH5"); interpolation_scheme = INTERPOLATION_P3M_IH6IH5
    case ("PQM_IH4IH3"); interpolation_scheme = INTERPOLATION_PQM_IH4IH3
    case ("PQM_IH6IH5"); interpolation_scheme = INTERPOLATION_PQM_IH6IH5
    case default ; call MOM_error(FATAL, "regrid_interp: "//&
     "Unrecognized choice for INTERPOLATION_SCHEME ("//trim(interp_scheme)//").")
  end select
end procedure interpolation_scheme
module procedure set_interp_scheme
  CS%interpolation_scheme = interpolation_scheme(interp_scheme)
end procedure set_interp_scheme
module procedure set_interp_extrap
  CS%boundary_extrapolation = extrap
end procedure set_interp_extrap
module procedure set_interp_answer_date
  CS%answer_date = answer_date
end procedure set_interp_answer_date
end submodule regrid_interp_s
