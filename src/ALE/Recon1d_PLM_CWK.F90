! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> Piecewise Linear Method 1D reconstruction
!!
!! This implementation of PLM follows Colella and Woodward, 1984, except for assuming
!! uniform cell thicknesses. Cells resort to PCM for extrema including first and last cells in column.
!! The cell-wise reconstructions are limited so that the edge values (which are also the
!! extrema in a cell) are bounded by the neighbor cell means. However, this does not yield
!! monotonic profiles for the whole column.
!!
!! Note that internally the edge values, rather than the PLM slope, are stored to ensure
!! resulting calculations are properly bounded.
module Recon1d_PLM_CWK

use Recon1d_type, only : testing
use Recon1d_PLM_CW, only : PLM_CW

implicit none ; private

public PLM_CWK, testing

!> PLM reconstruction following Colella and Woodward, 1984
!!
!! Implemented by extending recon1d_plm_cw.
!!
!! The source for the methods ultimately used by this class are:
!! - init()                 -> recon1d_plm_cw.init()
!! - reconstruct()             *locally defined
!! - average()              -> recon1d_plm_cw.average()
!! - f()                    -> recon1d_plm_cw.f()
!! - dfdx()                 -> recon1d_plm_cw.dfdx()
!! - x()                    -> recon1d_plm_cw.x()
!! - check_reconstruction() -> recon1d_plm_cw.check_reconstruction()
!! - unit_tests()           -> recon1d_plm_cw.unit_tests()
!! - destroy()              -> recon1d_plm_cw.destroy()
!! - remap_to_sub_grid()    -> recon1d_type.remap_to_sub_grid()
!! - init_parent()          -> init()
!! - reconstruct_parent()   -> reconstruct()
type, extends (PLM_CW) :: PLM_CWK

contains
  !> Implementation of the PLM_CWK reconstruction
  procedure :: reconstruct => reconstruct

end type PLM_CWK


  interface
module subroutine reconstruct(this, h, u)
  class(PLM_CWK), intent(inout) :: this !< This reconstruction
  real,           intent(in)    :: h(*) !< Grid spacing (thickness) [typically H]
  real,           intent(in)    :: u(*) !< Cell mean values [A]
  ! Local variables
                                    ! differences across the cell [A]

end subroutine reconstruct
  end interface

end module Recon1d_PLM_CWK
