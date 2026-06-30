! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> Linear interpolation functions
module P1M_functions

use regrid_edge_values, only : bound_edge_values, average_discontinuous_edge_values

implicit none ; private

! The following routines are visible to the outside world
public P1M_interpolation, P1M_boundary_extrapolation


  interface
module subroutine P1M_interpolation( N, h, u, edge_values, ppoly_coef, h_neglect, answer_date )
  integer,              intent(in)    :: N !< Number of cells
  real, dimension(:),   intent(in)    :: h !< cell widths (size N) [H]
  real, dimension(:),   intent(in)    :: u !< cell average properties (size N) [A]
  real, dimension(:,:), intent(inout) :: edge_values !< Potentially modified edge values [A]
  real, dimension(:,:), intent(inout) :: ppoly_coef !< Potentially modified
                                           !! piecewise polynomial coefficients, mainly [A]
  real,                 intent(in)    :: h_neglect !< A negligibly small width [H]
  integer,    optional, intent(in)    :: answer_date  !< The vintage of the expressions to use

  ! Local variables

  ! Bound edge values (routine found in 'edge_values.F90')
end subroutine P1M_interpolation
module subroutine P1M_boundary_extrapolation( N, h, u, edge_values, ppoly_coef )
  ! Arguments
  integer,              intent(in)    :: N !< Number of cells
  real, dimension(:),   intent(in)    :: h !< cell widths (size N) [H]
  real, dimension(:),   intent(in)    :: u !< cell averages (size N) [A]
  real, dimension(:,:), intent(inout) :: edge_values !< edge values of piecewise polynomials [A]
  real, dimension(:,:), intent(inout) :: ppoly_coef !< coefficients of piecewise polynomials, mainly [A]

  ! Local variables

  ! -----------------------------------------
  ! Left edge value in the left boundary cell
  ! -----------------------------------------
end subroutine P1M_boundary_extrapolation
  end interface

end module P1M_functions
