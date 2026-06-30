! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> Piecewise Linear Method 1D reconstruction in index space and boundary extrapolation
!!
!! This implementation of PLM follows Colella and Woodward, 1984 \cite colella1984, except for assuming
!! uniform resolution so that the method is independent of grid spacing. The cell-wise reconstructions
!! are limited so that the edge values (which are also the extrema in a cell) are bounded by the neighbors.
!! The slope of the first and last cells are set so that the first interior edge values match the interior
!! cell (i.e. extrapolates from the interior).
module Recon1d_EMPLM_CWK

use Recon1d_type, only : testing
use Recon1d_MPLM_CWK, only : MPLM_CWK

implicit none ; private

public EMPLM_CWK, testing

!> PLM reconstruction following Colella and Woodward, 1984
!!
!! Implemented by extending recon1d_mplm_cwk.
!!
!! The source for the methods ultimately used by this class are:
!! - init()                 -> recon1d_mplm_cwk -> recon1d_plm_cw.init()
!! - reconstruct()             *locally defined
!! - average()              -> recon1d_mplm_cwk -> recon1d_plm_cw.average()
!! - f()                    -> recon1d_mplm_cwk -> recon1d_plm_cw.f()
!! - dfdx()                 -> recon1d_mplm_cwk -> recon1d_plm_cw.dfdx()
!! - check_reconstruction() -> recon1d_mplm_cwk.check_reconstruction()
!! - unit_tests()              *locally defined
!! - destroy()              -> recon1d_mplm_cwk -> recon1d_plm_cw.destroy()
!! - remap_to_sub_grid()    -> recon1d_type.remap_to_sub_grid()
!! - init_parent()          -> init()
!! - reconstruct_parent()   -> recon1d_mplm_cwk.reconstruct()
type, extends (MPLM_CWK) :: EMPLM_CWK

contains
  !> Implementation of the EMPLM_CWK reconstruction
  procedure :: reconstruct => reconstruct
  !> Implementation of unit tests for the EMPLM_CWK reconstruction
  procedure :: unit_tests => unit_tests

end type EMPLM_CWK


  interface
module subroutine reconstruct(this, h, u)
  class(EMPLM_CWK), intent(inout) :: this !< This reconstruction
  real,             intent(in)    :: h(*) !< Grid spacing (thickness) [typically H]
  real,             intent(in)    :: u(*) !< Cell mean values [A]
  ! Local variables
                                    ! differences across the cell [A]

end subroutine reconstruct
logical module function unit_tests(this, verbose, stdout, stderr)
  class(EMPLM_CWK), intent(inout) :: this    !< This reconstruction
  logical,          intent(in)    :: verbose !< True, if verbose
  integer,          intent(in)    :: stdout  !< I/O channel for stdout
  integer,          intent(in)    :: stderr  !< I/O channel for stderr
  ! Local variables

end function unit_tests
  end interface

end module Recon1d_EMPLM_CWK
