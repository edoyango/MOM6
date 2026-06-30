! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> Module with routines for copying information from a shared dynamic horizontal
!! grid to an ocean-specific horizontal grid and the reverse.
module MOM_transcribe_grid

use MOM_array_transform, only : rotate_array, rotate_array_pair
use MOM_domains,         only : pass_var, pass_vector
use MOM_domains,         only : To_All, SCALAR_PAIR, CGRID_NE, AGRID, BGRID_NE, CORNER
use MOM_dyn_horgrid,     only : dyn_horgrid_type, set_derived_dyn_horgrid
use MOM_dyn_horgrid,     only : rotate_dyngrid=>rotate_dyn_horgrid
use MOM_error_handler,   only : MOM_error, MOM_mesg, FATAL, WARNING
use MOM_grid,            only : ocean_grid_type, set_derived_metrics
use MOM_unit_scaling,    only : unit_scale_type

implicit none ; private

public copy_dyngrid_to_MOM_grid, copy_MOM_grid_to_dyngrid, rotate_dyngrid


  interface
module subroutine copy_dyngrid_to_MOM_grid(dG, oG, US)
  type(dyn_horgrid_type), intent(in)    :: dG  !< Common horizontal grid type
  type(ocean_grid_type),  intent(inout) :: oG  !< Ocean grid type
  type(unit_scale_type),  intent(in)    :: US  !< A dimensional unit scaling type


  ! MOM_grid_init and create_dyn_horgrid are called outside of this routine.
  ! This routine copies over the fields that were set by MOM_initialized_fixed.

  ! Determine the indexing offsets between the grids.
end subroutine copy_dyngrid_to_MOM_grid
module subroutine copy_MOM_grid_to_dyngrid(oG, dG, US)
  type(ocean_grid_type),  intent(in)    :: oG  !< Ocean grid type
  type(dyn_horgrid_type), intent(inout) :: dG  !< Common horizontal grid type
  type(unit_scale_type),  intent(in)    :: US  !< A dimensional unit scaling type


  ! MOM_grid_init and create_dyn_horgrid are called outside of this routine.
  ! This routine copies over the fields that were set by MOM_initialized_fixed.

  ! Determine the indexing offsets between the grids.
end subroutine copy_MOM_grid_to_dyngrid
  end interface

end module MOM_transcribe_grid
