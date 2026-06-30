! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> Piecewise quartic reconstruction functions
module PQM_functions

use regrid_edge_values, only : bound_edge_values, check_discontinuous_edge_values

implicit none ; private

public PQM_reconstruction, PQM_boundary_extrapolation, PQM_boundary_extrapolation_v1


  interface
module subroutine PQM_reconstruction( N, h, u, edge_values, edge_slopes, ppoly_coef, h_neglect, answer_date )
  integer,              intent(in)    :: N !< Number of cells
  real, dimension(:),   intent(in)    :: h !< cell widths (size N) [H]
  real, dimension(:),   intent(in)    :: u !< cell averages (size N) [A]
  real, dimension(:,:), intent(inout) :: edge_values    !< Edge value of polynomial [A]
  real, dimension(:,:), intent(inout) :: edge_slopes    !< Edge slope of polynomial [A H-1]
  real, dimension(:,:), intent(inout) :: ppoly_coef !< Coefficients of polynomial, mainly [A]
  real,                 intent(in)    :: h_neglect  !< A negligibly small width for
                                           !! the purpose of cell reconstructions [H]
  integer,    optional, intent(in)    :: answer_date  !< The vintage of the expressions to use

  ! Local variables

  ! PQM limiter
end subroutine PQM_reconstruction
module subroutine PQM_limiter( N, h, u, edge_values, edge_slopes, h_neglect, answer_date )
  integer,              intent(in)    :: N !< Number of cells
  real, dimension(:),   intent(in)    :: h !< cell widths (size N) [H]
  real, dimension(:),   intent(in)    :: u !< cell average properties (size N) [A]
  real, dimension(:,:), intent(inout) :: edge_values !< Potentially modified edge values [A]
  real, dimension(:,:), intent(inout) :: edge_slopes !< Potentially modified edge slopes [A H-1]
  real,                 intent(in)    :: h_neglect !< A negligibly small width for
                                           !! the purpose of cell reconstructions [H]
  integer,    optional, intent(in)    :: answer_date  !< The vintage of the expressions to use

  ! Local variables

  ! Bound edge values
end subroutine PQM_limiter
module subroutine PQM_boundary_extrapolation( N, h, u, edge_values, ppoly_coef )
  integer,              intent(in)    :: N !< Number of cells
  real, dimension(:),   intent(in)    :: h !< cell widths (size N) [H]
  real, dimension(:),   intent(in)    :: u !< cell averages (size N) [A]
  real, dimension(:,:), intent(inout) :: edge_values    !< Edge value of polynomial [A]
  real, dimension(:,:), intent(inout) :: ppoly_coef !< Coefficients of polynomial, mainly [A]
  ! Local variables

  ! ----- Left boundary -----
end subroutine PQM_boundary_extrapolation
module subroutine PQM_boundary_extrapolation_v1( N, h, u, edge_values, edge_slopes, ppoly_coef, h_neglect )
  integer,              intent(in)    :: N !< Number of cells
  real, dimension(:),   intent(in)    :: h !< cell widths (size N) [H]
  real, dimension(:),   intent(in)    :: u !< cell averages (size N) [A]
  real, dimension(:,:), intent(inout) :: edge_values    !< Edge value of polynomial [A]
  real, dimension(:,:), intent(inout) :: edge_slopes    !< Edge slope of polynomial [A H-1]
  real, dimension(:,:), intent(inout) :: ppoly_coef !< Coefficients of polynomial, mainly [A]
  real,                 intent(in)    :: h_neglect  !< A negligibly small width for
                                           !! the purpose of cell reconstructions [H]
  ! Local variables

  ! ----- Left boundary (TOP) -----
end subroutine PQM_boundary_extrapolation_v1
  end interface

end module PQM_functions
