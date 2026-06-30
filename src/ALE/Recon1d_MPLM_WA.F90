! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> Monotonized Piecewise Linear Method 1D reconstruction
!!
!! This implementation of PLM follows White and Adcroft, 2008 \cite white2008.
!! The PLM slopes are first limited following Colella and Woodward, 1984, but are then
!! further limited to ensure the edge values moving across cell boundaries are monotone.
!! The first and last cells are always limited to PCM.
!!
!! This differs from recon1d_mplm_wa_poly in the internally not polynomial representations
!! are referred to.
module Recon1d_MPLM_WA

use Recon1d_PLM_CW, only : PLM_CW, testing

implicit none ; private

public MPLM_WA, testing

!> Limited Monotonic PLM reconstruction following White and Adcroft, 2008
!!
!! The source for the methods ultimately used by this class are:
!! - init()                 -> recon1d_plm_cw.init()
!! - reconstruct()             *locally defined
!! - average()              -> recon1d_plm_cw.average()
!! - f()                    -> recon1d_plm_cw.f()
!! - dfdx()                 -> recon1d_plm_cw.dfdx()
!! - check_reconstruction()    *locally defined
!! - unit_tests()              *locally defined
!! - destroy()              -> recon1d_plm_cw.destroy()
!! - remap_to_sub_grid()    -> recon1d_plm_cw           -> recon1d_type.remap_to_sub_grid()
!! - init_parent()          -> recon1d_plm_cw.init()
!! - reconstruct_parent()   -> reconstruct()
type, extends (PLM_CW) :: MPLM_WA

contains
  !> Implementation of the MPLM_WA reconstruction
  procedure :: reconstruct => reconstruct
  !> Implementation of check reconstruction for the MPLM_WA reconstruction
  procedure :: check_reconstruction => check_reconstruction
  !> Implementation of unit tests for the MPLM_WA reconstruction
  procedure :: unit_tests => unit_tests

  !> Duplicate interface to reconstruct()
  procedure :: reconstruct_parent => reconstruct

end type MPLM_WA


  interface
module subroutine reconstruct(this, h, u)
  class(MPLM_WA), intent(inout) :: this !< This reconstruction
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

  ! Comparison are made assuming +ve slopes
end function PLM_monotonized_slope
logical module function check_reconstruction(this, h, u)
  class(MPLM_WA), intent(in) :: this !< This reconstruction
  real,           intent(in) :: h(*) !< Grid spacing (thickness) [typically H]
  real,           intent(in) :: u(*) !< Cell mean values [A]
  ! Local variables

end function check_reconstruction
logical module function unit_tests(this, verbose, stdout, stderr)
  class(MPLM_WA), intent(inout) :: this    !< This reconstruction
  logical,        intent(in)    :: verbose !< True, if verbose
  integer,        intent(in)    :: stdout  !< I/O channel for stdout
  integer,        intent(in)    :: stderr  !< I/O channel for stderr
  ! Local variables

end function unit_tests
  end interface

end module Recon1d_MPLM_WA
