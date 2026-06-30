! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> Extrapolated-Monotonized Piecewise Linear Method 1D reconstruction
!!
!! This extends MPLM_WA, following White and Adcroft, 2008 \cite white2008, by extrapolating for the slopes of the
!! first and last cells. This extrapolation is used by White et al., 2009, during grid-generation.
module Recon1d_EMPLM_WA

use Recon1d_MPLM_WA, only : MPLM_WA, testing

implicit none ; private

public EMPLM_WA

!> Extraplated Monotonic PLM reconstruction of White and Adcroft, 2008
!!
!! The source for the methods ultimately used by this class are:
!! - init()                 -> recon1d_mplm_wa           -> recon1d_plm_cw.init()
!! - reconstruct()             *locally defined
!! - average()              -> recon1d_mplm_wa           -> recon1d_plm_cw.average()
!! - f()                    -> recon1d_mplm_wa           -> recon1d_plm_cw.f()
!! - dfdx()                 -> recon1d_mplm_wa           -> recon1d_plm_cw.dfdx()
!! - check_reconstruction() -> recon1d_mplm_wa.check_reconstruction()
!! - unit_tests()              *locally defined
!! - destroy()              -> recon1d_mplm_wa           -> recon1d_plm_cw.destroy()
!! - remap_to_sub_grid()    -> recon1d_mplm_wa           -> recon1d_plm_cw        -> recon1d_type.remap_to_sub_grid()
!! - init_parent()          -> recon1d_mplm_wa           -> recon1d_plm_cw.init()
!! - reconstruct_parent()   -> recon1d_mplm_wa.reconstruct()
type, extends (MPLM_WA) :: EMPLM_WA

contains
  !> Implementation of the EMPLM_WA reconstruction with boundary extrapolation
  procedure :: reconstruct => reconstruct
  !> Implementation of unit tests for the EMPLM_WA reconstruction
  procedure :: unit_tests => unit_tests

end type EMPLM_WA


  interface
module subroutine reconstruct(this, h, u)
  class(EMPLM_WA), intent(inout) :: this !< This reconstruction
  real,            intent(in)    :: h(*) !< Grid spacing (thickness) [typically H]
  real,            intent(in)    :: u(*) !< Cell mean values [A]
  ! Local variables

  ! Use parent (MPLM_WA) reconstruction
end subroutine reconstruct
real elemental pure module function PLM_extrapolate_slope(h_l, h_c, h_neglect, u_l, u_c)
  real, intent(in) :: h_l !< Thickness of left cell in arbitrary grid thickness units [H]
  real, intent(in) :: h_c !< Thickness of center cell in arbitrary grid thickness units [H]
  real, intent(in) :: h_neglect !< A negligible thickness [H]
  real, intent(in) :: u_l !< Value of left cell in arbitrary units [A]
  real, intent(in) :: u_c !< Value of center cell in arbitrary units [A]
  ! Local variables

  ! Avoid division by zero for vanished cells
end function PLM_extrapolate_slope
logical module function unit_tests(this, verbose, stdout, stderr)
  class(EMPLM_WA), intent(inout) :: this    !< This reconstruction
  logical,         intent(in)    :: verbose !< True, if verbose
  integer,         intent(in)    :: stdout  !< I/O channel for stdout
  integer,         intent(in)    :: stderr  !< I/O channel for stderr
  ! Local variables

end function unit_tests
  end interface

end module Recon1d_EMPLM_WA
