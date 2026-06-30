! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> Provides functions used with the Piecewise-Parabolic-Method in the vertical ALE algorithm.
module PPM_functions

! First version was created by Laurent White, June 2008.
! Substantially re-factored January 2016.

!! @todo Re-factor PPM_boundary_extrapolation to give round-off safe and
!!       optimization independent results.

use regrid_edge_values, only : bound_edge_values, check_discontinuous_edge_values

implicit none ; private

public PPM_reconstruction, PPM_boundary_extrapolation, PPM_monotonicity


  interface
module subroutine PPM_reconstruction( N, h, u, edge_values, ppoly_coef, h_neglect, answer_date)
  integer,              intent(in)    :: N !< Number of cells
  real, dimension(N),   intent(in)    :: h !< Cell widths [H]
  real, dimension(N),   intent(in)    :: u !< Cell averages in arbitrary coordinates [A]
  real, dimension(N,2), intent(inout) :: edge_values !< Edge values [A]
  real, dimension(N,3), intent(inout) :: ppoly_coef !< Polynomial coefficients, mainly [A]
  real,                 intent(in)    :: h_neglect !< A negligibly small width [H]
  integer,    optional, intent(in)    :: answer_date  !< The vintage of the expressions to use

  ! Local variables

  ! PPM limiter
end subroutine PPM_reconstruction
module subroutine PPM_limiter_standard( N, h, u, edge_values, h_neglect, answer_date )
  integer,              intent(in)    :: N !< Number of cells
  real, dimension(:),   intent(in)    :: h !< cell widths (size N) [H]
  real, dimension(:),   intent(in)    :: u !< cell average properties (size N) [A]
  real, dimension(:,:), intent(inout) :: edge_values !< Potentially modified edge values [A]
  real,                 intent(in)    :: h_neglect !< A negligibly small width [H]
  integer,    optional, intent(in)    :: answer_date  !< The vintage of the expressions to use

  ! Local variables

  ! Bound edge values
end subroutine PPM_limiter_standard
module subroutine PPM_monotonicity( N, u, edge_values )
  integer,              intent(in)    :: N !< Number of cells
  real, dimension(:),   intent(in)    :: u !< cell average properties (size N) [A]
  real, dimension(:,:), intent(inout) :: edge_values !< Potentially modified edge values [A]

  ! Local variables

  ! Loop on interior cells to impose monotonicity
  ! Eq. 1.10 of (Colella & Woodward, JCP 84)
end subroutine PPM_monotonicity
module subroutine PPM_boundary_extrapolation( N, h, u, edge_values, ppoly_coef, h_neglect)
!------------------------------------------------------------------------------
! Reconstruction by parabolas within boundary cells.
!
! The following explanations apply to the left boundary cell. The same
! reasoning holds for the right boundary cell.
!
! A parabola needs to be built in the cell and requires three degrees of
! freedom, which are the right edge value and slope and the cell average.
! The right edge values and slopes are taken to be that of the neighboring
! cell (i.e., the left edge value and slope of the neighboring cell).
! The resulting parabola is not necessarily monotonic and the traditional
! PPM limiter is used to modify one of the edge values in order to yield
! a monotonic parabola.
!
! N:     number of cells in grid
! h:     thicknesses of grid cells
! u:     cell averages to use in constructing piecewise polynomials
! edge_values : edge values of piecewise polynomials
! ppoly_coef : coefficients of piecewise polynomials
!
! It is assumed that the size of the array 'u' is equal to the number of cells
! defining 'grid' and 'ppoly'. No consistency check is performed here.
!------------------------------------------------------------------------------

  ! Arguments
  integer,              intent(in)    :: N !< Number of cells
  real, dimension(:),   intent(in)    :: h !< cell widths (size N) [H]
  real, dimension(:),   intent(in)    :: u !< cell averages (size N) [A]
  real, dimension(:,:), intent(inout) :: edge_values    !< edge values of piecewise polynomials [A]
  real, dimension(:,:), intent(inout) :: ppoly_coef !< coefficients of piecewise polynomials, mainly [A]
  real,                 intent(in)    :: h_neglect  !< A negligibly small width for
                                           !! the purpose of cell reconstructions [H]

  ! Local variables
                        ! of a reconstructed distribution [A]
                        ! the cell being worked on [A]

  ! ----- Left boundary -----
end subroutine PPM_boundary_extrapolation
  end interface

end module PPM_functions
