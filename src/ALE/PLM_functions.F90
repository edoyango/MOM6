! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> Piecewise linear reconstruction functions
module PLM_functions

implicit none ; private

public PLM_boundary_extrapolation
public PLM_extrapolate_slope
public PLM_monotonized_slope
public PLM_reconstruction
public PLM_slope_wa
public PLM_slope_cw


  interface
real elemental pure module function PLM_slope_wa(h_l, h_c, h_r, h_neglect, u_l, u_c, u_r)
  real, intent(in) :: h_l !< Thickness of left cell in arbitrary grid thickness units [H]
  real, intent(in) :: h_c !< Thickness of center cell in arbitrary grid thickness units [H]
  real, intent(in) :: h_r !< Thickness of right cell in arbitrary grid thickness units [H]
  real, intent(in) :: h_neglect !< A negligible thickness [H]
  real, intent(in) :: u_l !< Value of left cell in arbitrary units [A]
  real, intent(in) :: u_c !< Value of center cell in arbitrary units [A]
  real, intent(in) :: u_r !< Value of right cell in arbitrary units [A]
  ! Local variables
                                    ! differences across the cell [A]

  ! Side differences
end function PLM_slope_wa
real elemental pure module function PLM_slope_cw(h_l, h_c, h_r, h_neglect, u_l, u_c, u_r)
  real, intent(in) :: h_l !< Thickness of left cell in arbitrary grid thickness units [H]
  real, intent(in) :: h_c !< Thickness of center cell in arbitrary grid thickness units [H]
  real, intent(in) :: h_r !< Thickness of right cell in arbitrary grid thickness units [H]
  real, intent(in) :: h_neglect !< A negligible thickness [H]
  real, intent(in) :: u_l !< Value of left cell in arbitrary units [A]
  real, intent(in) :: u_c !< Value of center cell in arbitrary units [A]
  real, intent(in) :: u_r !< Value of right cell in arbitrary units [A]
  ! Local variables
                                    ! differences across the cell [A]

end function PLM_slope_cw
real elemental pure module function PLM_monotonized_slope(u_l, u_c, u_r, s_l, s_c, s_r)
  real, intent(in) :: u_l !< Value of left cell in arbitrary units [A]
  real, intent(in) :: u_c !< Value of center cell in arbitrary units [A]
  real, intent(in) :: u_r !< Value of right cell in arbitrary units [A]
  real, intent(in) :: s_l !< PLM slope of left cell [A]
  real, intent(in) :: s_c !< PLM slope of center cell [A]
  real, intent(in) :: s_r !< PLM slope of right cell [A]
  ! Local variables

end function PLM_monotonized_slope
real elemental pure module function PLM_extrapolate_slope(h_l, h_c, h_neglect, u_l, u_c)
  real, intent(in) :: h_l !< Thickness of left cell in arbitrary grid thickness units [H]
  real, intent(in) :: h_c !< Thickness of center cell in arbitrary grid thickness units [H]
  real, intent(in) :: h_neglect !< A negligible thickness [H]
  real, intent(in) :: u_l !< Value of left cell in arbitrary units [A]
  real, intent(in) :: u_c !< Value of center cell in arbitrary units [A]
  ! Local variables

  ! Avoid division by zero for vanished cells
end function PLM_extrapolate_slope
module subroutine PLM_reconstruction( N, h, u, edge_values, ppoly_coef, h_neglect )
  integer,              intent(in)    :: N !< Number of cells
  real, dimension(:),   intent(in)    :: h !< cell widths (size N) [H]
  real, dimension(:),   intent(in)    :: u !< cell averages (size N) in arbitrary units [A]
  real, dimension(:,:), intent(inout) :: edge_values !< edge values of piecewise polynomials,
                                           !! with the same units as u [A].
  real, dimension(:,:), intent(inout) :: ppoly_coef !< coefficients of piecewise polynomials, mainly
                                           !! with the same units as u [A].
  real,                 intent(in)    :: h_neglect !< A negligibly small width for
                                           !! the purpose of cell reconstructions
                                           !! in the same units as h [H]

  ! Local variables

end subroutine PLM_reconstruction
module subroutine PLM_boundary_extrapolation( N, h, u, edge_values, ppoly_coef, h_neglect )
  integer,              intent(in)    :: N !< Number of cells
  real, dimension(:),   intent(in)    :: h !< cell widths (size N) [H]
  real, dimension(:),   intent(in)    :: u !< cell averages (size N) in arbitrary units [A]
  real, dimension(:,:), intent(inout) :: edge_values !< edge values of piecewise polynomials,
                                           !! with the same units as u [A].
  real, dimension(:,:), intent(inout) :: ppoly_coef !< coefficients of piecewise polynomials, mainly
                                           !! with the same units as u [A].
  real,                 intent(in)    :: h_neglect !< A negligibly small width for
                                           !! the purpose of cell reconstructions
                                           !! in the same units as h [H]
  ! Local variables

  ! Extrapolate from 2 to 1 to estimate slope
end subroutine PLM_boundary_extrapolation
  end interface

end module PLM_functions
