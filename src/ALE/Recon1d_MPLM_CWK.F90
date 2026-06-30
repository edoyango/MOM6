! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> Piecewise Linear Method 1D reconstruction in index space
!!
!! This implementation of PLM follows Colella and Woodward, 1984 \cite colella1984, except for assuming
!! uniform resolution so that the method is independent of grid spacing. The cell-wise reconstructions
!! are limited so that the edge values (which are also the extrema in a cell) are bounded by the neighbors.
!! The first and last cells are always limited to PCM.
module Recon1d_MPLM_CWK

use Recon1d_type, only : testing
use Recon1d_PLM_CWK, only : PLM_CWK

implicit none ; private

public MPLM_CWK, testing

!> PLM reconstruction following Colella and Woodward, 1984
!!
!! Implemented by extending recon1d_plm_cwk.
!!
!! The source for the methods ultimately used by this class are:
!! - init()                 -> recon1d_plm_cwk -> recon1d_plm_cw.init()
!! - reconstruct()             *locally defined
!! - average()              -> recon1d_plm_cwk -> recon1d_plm_cw.average()
!! - f()                    -> recon1d_plm_cwk -> recon1d_plm_cw.f()
!! - dfdx()                 -> recon1d_plm_cwk -> recon1d_plm_cw.dfdx()
!! - check_reconstruction()    *locally defined
!! - unit_tests()              *locally defined
!! - destroy()              -> recon1d_plm_cwk -> recon1d_plm_cw.destroy()
!! - remap_to_sub_grid()    -> recon1d_type.remap_to_sub_grid()
!! - init_parent()          -> init()
!! - reconstruct_parent()   -> reconstruct()
type, extends (PLM_CWK) :: MPLM_CWK

contains
  !> Implementation of the MPLM_CWK reconstruction
  procedure :: reconstruct => reconstruct
  !> Implementation of check reconstruction for the MPLM_CWK reconstruction
  procedure :: check_reconstruction => check_reconstruction
  !> Implementation of unit tests for the MPLM_CWK reconstruction
  procedure :: unit_tests => unit_tests

  !> Duplicate interface to reconstruct()
  procedure :: reconstruct_parent => reconstruct
end type MPLM_CWK


  interface
module subroutine reconstruct(this, h, u)
  class(MPLM_CWK), intent(inout) :: this !< This reconstruction
  real,            intent(in)    :: h(*) !< Grid spacing (thickness) [typically H]
  real,            intent(in)    :: u(*) !< Cell mean values [A]
  ! Local variables
                                    ! differences across the cell [A]

end subroutine reconstruct
logical module function check_reconstruction(this, h, u)
  class(MPLM_CWK), intent(in) :: this !< This reconstruction
  real,            intent(in) :: h(*) !< Grid spacing (thickness) [typically H]
  real,            intent(in) :: u(*) !< Cell mean values [A]
  ! Local variables

end function check_reconstruction
logical module function unit_tests(this, verbose, stdout, stderr)
  class(MPLM_CWK), intent(inout) :: this    !< This reconstruction
  logical,         intent(in)    :: verbose !< True, if verbose
  integer,         intent(in)    :: stdout  !< I/O channel for stdout
  integer,         intent(in)    :: stderr  !< I/O channel for stderr
  ! Local variables

end function unit_tests
  end interface

end module Recon1d_MPLM_CWK
