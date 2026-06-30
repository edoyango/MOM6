! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> Cubic interpolation functions
module P3M_functions

use regrid_edge_values, only : bound_edge_values, average_discontinuous_edge_values

implicit none ; private

public P3M_interpolation
public P3M_boundary_extrapolation


  interface
module subroutine P3M_interpolation( N, h, u, edge_values, ppoly_S, ppoly_coef, h_neglect, answer_date )
  integer,              intent(in)    :: N !< Number of cells
  real, dimension(:),   intent(in)    :: h !< cell widths (size N) [H]
  real, dimension(:),   intent(in)    :: u !< cell averages (size N) in arbitrary units [A]
  real, dimension(:,:), intent(inout) :: edge_values   !< Edge value of polynomial [A]
  real, dimension(:,:), intent(inout) :: ppoly_S   !< Edge slope of polynomial [A H-1].
  real, dimension(:,:), intent(inout) :: ppoly_coef !< Coefficients of polynomial [A]
  real,                 intent(in)    :: h_neglect !< A negligibly small width for the
                                          !! purpose of cell reconstructions [H]
  integer,    optional, intent(in)    :: answer_date  !< The vintage of the expressions to use

  ! Call the limiter for p3m, which takes care of everything from
  ! computing the coefficients of the cubic to monotonizing it.
  ! This routine could be called directly instead of having to call
  ! 'P3M_interpolation' first but we do that to provide an homogeneous
  ! interface.
end subroutine P3M_interpolation
module subroutine P3M_limiter( N, h, u, edge_values, ppoly_S, ppoly_coef, h_neglect, answer_date )
  integer,              intent(in)    :: N !< Number of cells
  real, dimension(:),   intent(in)    :: h !< cell widths (size N) [H]
  real, dimension(:),   intent(in)    :: u !< cell averages (size N) in arbitrary units [A]
  real, dimension(:,:), intent(inout) :: edge_values !< Edge value of polynomial [A]
  real, dimension(:,:), intent(inout) :: ppoly_S  !< Edge slope of polynomial [A H-1]
  real, dimension(:,:), intent(inout) :: ppoly_coef !< Coefficients of polynomial [A]
  real,                 intent(in)    :: h_neglect !< A negligibly small width for
                                           !! the purpose of cell reconstructions [H]
  integer,    optional, intent(in)    :: answer_date  !< The vintage of the expressions to use

  ! Local variables

  ! 1. Bound edge values (boundary cells are assumed to be local extrema)
end subroutine P3M_limiter
module subroutine P3M_boundary_extrapolation( N, h, u, edge_values, ppoly_S, ppoly_coef, &
                                       h_neglect, h_neglect_edge )
  integer,              intent(in)    :: N !< Number of cells
  real, dimension(:),   intent(in)    :: h !< cell widths (size N) [H]
  real, dimension(:),   intent(in)    :: u !< cell averages (size N) in arbitrary units [A]
  real, dimension(:,:), intent(inout) :: edge_values !< Edge value of polynomial [A]
  real, dimension(:,:), intent(inout) :: ppoly_S !< Edge slope of polynomial [A H-1]
  real, dimension(:,:), intent(inout) :: ppoly_coef !< Coefficients of polynomial [A]
  real,                 intent(in)    :: h_neglect !< A negligibly small width for the
                                          !! purpose of cell reconstructions [H]
  real,       optional, intent(in)    :: h_neglect_edge !< A negligibly small width for the purpose
                                          !! of finding edge values [H].  The default is h_neglect.
  ! Local variables

end subroutine P3M_boundary_extrapolation
module subroutine build_cubic_interpolant( h, k, edge_values, ppoly_S, ppoly_coef )
  real, dimension(:),   intent(in)    :: h !< cell widths (size N) [H]
  integer,              intent(in)    :: k !< The index of the cell to work on
  real, dimension(:,:), intent(in)    :: edge_values !< Edge value of polynomial in arbitrary units [A]
  real, dimension(:,:), intent(in)    :: ppoly_S    !< Edge slope of polynomial [A H-1]
  real, dimension(:,:), intent(inout) :: ppoly_coef !< Coefficients of polynomial [A]

  ! Local variables

end subroutine build_cubic_interpolant
logical module function is_cubic_monotonic( ppoly_coef, k )
  real, dimension(:,:), intent(in) :: ppoly_coef !< Coefficients of cubic polynomial in arbitrary units [A]
  integer,              intent(in) :: k  !< The index of the cell to work on
  ! Local variables

end function is_cubic_monotonic
module subroutine monotonize_cubic( h, u0_l, u0_r, sigma_l, sigma_r, slope, u1_l, u1_r )
  real, intent(in)      :: h       !< cell width [H]
  real, intent(in)      :: u0_l    !< left edge value in arbitrary units [A]
  real, intent(in)      :: u0_r    !< right edge value [A]
  real, intent(in)      :: sigma_l !< left 2nd-order slopes [A H-1]
  real, intent(in)      :: sigma_r !< right 2nd-order slopes [A H-1]
  real, intent(in)      :: slope   !< limited PLM slope [A H-1]
  real, intent(inout)   :: u1_l    !< left edge slopes [A H-1]
  real, intent(inout)   :: u1_r    !< right edge slopes [A H-1]
  ! Local variables

end subroutine monotonize_cubic
  end interface

end module P3M_functions
