! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> Extrapolated-Monotonized Piecewise Linear Method 1D reconstruction
!!
!! This extends MPLM_poly, following White and Adcroft, 2008 \cite white2008, by extraplating for the slopes of the
!! first and last cells. This extrapolation is used by White et al., 2009, during grid-generation.
!!
!! This stores and evaluates the reconstruction using a polynomial representation which is not preferred
!! but was the form used in OM4.
module Recon1d_EMPLM_WA_poly

use Recon1d_MPLM_WA_poly, only : MPLM_WA_poly, testing

implicit none ; private

public EMPLM_WA_poly

!> Extrapolation Limited Monotonic PLM reconstruction following White and Adcroft, 2008
!!
!! The source for the methods ultimately used by this class are:
!! - init()                 -> recon1d_mplm_wa_poly.init()
!! - reconstruct()          -> recon1d_mplm_wa_poly.reconstruct()
!! - average()              -> recon1d_mplm_wa_poly.average()
!! - f()                    -> recon1d_mplm_wa_poly           -> recon1d_mplm_wa -> recon1d_plm_cw.f()
!! - dfdx()                 -> recon1d_mplm_wa_poly           -> recon1d_mplm_wa -> recon1d_plm_cw.dfdx()
!! - check_reconstruction()    *locally defined
!! - unit_tests()              *locally defined
!! - destroy()              -> recon1d_mplm_wa_poly           -> recon1d_mplm_wa -> recon1d_plm_cw.destroy()
!! - remap_to_sub_grid()       *locally defined
!! - init_parent()          -> recon1d_mplm_wa_poly           -> recon1d_mplm_wa.init()
!! - reconstruct_parent()   -> recon1d_mplm_wa_poly           -> recon1d_mplm_wa.reconstruct()
type, extends (MPLM_WA_poly) :: EMPLM_WA_poly

contains
  !> Implementation of the EMPLM_WA_poly reconstruction with boundary extrapolation
  procedure :: reconstruct => reconstruct
  !> Implementation of check reconstruction for the EMPLM_WA_poly reconstruction
  procedure :: check_reconstruction => check_reconstruction
  !> Implementation of unit tests for the EMPLM_WA_poly reconstruction
  procedure :: unit_tests => unit_tests

end type EMPLM_WA_poly


  interface
module subroutine reconstruct(this, h, u)
  class(EMPLM_WA_poly), intent(inout) :: this !< This reconstruction
  real,                 intent(in)    :: h(*) !< Grid spacing (thickness) [typically H]
  real,                 intent(in)    :: u(*) !< Cell mean values [A]
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
logical module function check_reconstruction(this, h, u)
  class(EMPLM_WA_poly), intent(in) :: this !< This reconstruction
  real,                 intent(in) :: h(*) !< Grid spacing (thickness) [typically H]
  real,                 intent(in) :: u(*) !< Cell mean values [A]
  ! Local variables

end function check_reconstruction
logical module function unit_tests(this, verbose, stdout, stderr)
  class(EMPLM_WA_poly), intent(inout) :: this    !< This reconstruction
  logical,              intent(in)    :: verbose !< True, if verbose
  integer,              intent(in)    :: stdout  !< I/O channel for stdout
  integer,              intent(in)    :: stderr  !< I/O channel for stderr
  ! Local variables

end function unit_tests
  end interface

end module Recon1d_EMPLM_WA_poly
