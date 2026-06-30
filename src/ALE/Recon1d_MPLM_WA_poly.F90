! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> Monotonized Piecewise Linear Method 1D reconstruction using polynomial representation
!!
!! This implementation of PLM follows White and Adcroft, 2008 \cite white2008.
!! The PLM slopes are first limited following Colella and Woodward, 1984, but are then
!! further limited to ensure the edge values moving across cell boundaries are monotone.
!! The first and last cells are always limited to PCM.
!!
!! This stores and evaluates the reconstruction using a polynomial representation which is
!! not preferred but was the form used in OM4.
module Recon1d_MPLM_WA_poly

use Recon1d_MPLM_WA, only : MPLM_WA, testing

implicit none ; private

public MPLM_WA_poly, testing

!> Limited Monotonic PLM reconstruction following White and Adcroft, 2008
!!
!! The source for the methods ultimately used by this class are:
!! - init()                    *locally defined
!! - reconstruct()             *locally defined
!! - average()                 *locally defined
!! - f()                    -> recon1d_mplm_wa   -> recon1d_plm_cw.f()
!! - dfdx()                 -> recon1d_mplm_wa   -> recon1d_plm_cw.dfdx()
!! - check_reconstruction()    *locally defined
!! - unit_tests()              *locally defined
!! - destroy()              -> recon1d_mplm_wa   -> recon1d_plm_cw.destroy()
!! - remap_to_sub_grid()       *locally defined
!! - init_parent()          -> init()
!! - reconstruct_parent()   -> reconstruct()
type, extends (MPLM_WA) :: MPLM_WA_poly

  ! Legacy representation
  integer :: degree !< Degree of polynomial used in legacy representation
  real, allocatable, dimension(:,:) :: poly_coef !< Polynomial coefficients in legacy representation

contains
  !> Implementation of the MPLM_WA_poly initialization
  procedure :: init => init
  !> Implementation of the MPLM_WA_poly reconstruction
  procedure :: reconstruct => reconstruct
  !> Implementation of the MPLM_WA_poly average over an interval [A]
  procedure :: average => average
  !> Implementation of check reconstruction for the MPLM_WA_poly reconstruction
  procedure :: check_reconstruction => check_reconstruction
  !> Implementation of unit tests for the MPLM_WA_poly reconstruction
  procedure :: unit_tests => unit_tests

  !> Duplicate interface to init()
  procedure :: init_parent => init
  !> Duplicate interface to reconstruct()
  procedure :: reconstruct_parent => reconstruct

#undef USE_BASE_CLASS_REMAP
#ifndef USE_BASE_CLASS_REMAP
! This block is here to test whether the compiler can do better if we have local copies of
! the remapping functions.
  !> Remaps the column to subgrid h_sub
  procedure :: remap_to_sub_grid => remap_to_sub_grid
#endif

end type MPLM_WA_poly


  interface
module subroutine init(this, n, h_neglect, check)
  class(MPLM_WA_poly), intent(out) :: this      !< This reconstruction
  integer,             intent(in)  :: n         !< Number of cells in this column
  real, optional,      intent(in)  :: h_neglect !< A negligibly small width used in cell reconstructions [H]
  logical, optional,   intent(in)  :: check     !< If true, enable some consistency checking

end subroutine init
module subroutine reconstruct(this, h, u)
  class(MPLM_WA_poly), intent(inout) :: this !< This reconstruction
  real,           intent(in)    :: h(*) !< Grid spacing (thickness) [typically H]
  real,           intent(in)    :: u(*) !< Cell mean values [A]
  ! Local variables

end subroutine reconstruct
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
real elemental pure module function PLM_monotonized_slope(u_l, u_c, u_r, s_l, s_c, s_r)
  real, intent(in) :: u_l !< Value of left cell in arbitrary units [A]
  real, intent(in) :: u_c !< Value of center cell in arbitrary units [A]
  real, intent(in) :: u_r !< Value of right cell in arbitrary units [A]
  real, intent(in) :: s_l !< PLM slope of left cell [A]
  real, intent(in) :: s_c !< PLM slope of center cell [A]
  real, intent(in) :: s_r !< PLM slope of right cell [A]
  ! Local variables

end function PLM_monotonized_slope
real module function average(this, k, xa, xb)
  class(MPLM_WA_poly), intent(in) :: this !< This reconstruction
  integer,        intent(in) :: k    !< Cell number
  real,           intent(in) :: xa   !< Start of averaging interval on element (0 to 1)
  real,           intent(in) :: xb   !< End of averaging interval on element (0 to 1)

end function average
module subroutine remap_to_sub_grid(this, h0, u0, n1, h_sub, &
                                   isrc_start, isrc_end, isrc_max, isub_src, &
                                   u_sub, uh_sub, u02_err)
  class(MPLM_WA_poly), intent(in) :: this !< 1-D reconstruction type
  real,    intent(in)  :: h0(*)  !< Source grid widths (size n0) [H]
  real,    intent(in)  :: u0(*)  !< Source grid widths (size n0) [H]
  integer, intent(in)  :: n1      !< Number of cells in target grid
  real,    intent(in)  :: h_sub(*) !< Overlapping sub-cell thicknesses, h_sub [H]
  integer, intent(in)  :: isrc_start(*) !< Index of first sub-cell within each source cell
  integer, intent(in)  :: isrc_end(*) !< Index of last sub-cell within each source cell
  integer, intent(in)  :: isrc_max(*) !< Index of thickest sub-cell within each source cell
  integer, intent(in)  :: isub_src(*) !< Index of source cell for each sub-cell
  real,    intent(out) :: u_sub(*) !< Sub-cell cell averages (size n1) [A]
  real,    intent(out) :: uh_sub(*) !< Sub-cell cell integrals (size n1) [A H]
  real,    intent(out) :: u02_err !< Integrated reconstruction error estimates [A H]
  ! Local variables

end subroutine remap_to_sub_grid
logical module function check_reconstruction(this, h, u)
  class(MPLM_WA_poly), intent(in) :: this !< This reconstruction
  real,                intent(in) :: h(*) !< Grid spacing (thickness) [typically H]
  real,                intent(in) :: u(*) !< Cell mean values [A]
  ! Local variables

end function check_reconstruction
logical module function unit_tests(this, verbose, stdout, stderr)
  class(MPLM_WA_poly), intent(inout) :: this    !< This reconstruction
  logical,             intent(in)    :: verbose !< True, if verbose
  integer,             intent(in)    :: stdout  !< I/O channel for stdout
  integer,             intent(in)    :: stderr  !< I/O channel for stderr
  ! Local variables

end function unit_tests
  end interface

end module Recon1d_MPLM_WA_poly
