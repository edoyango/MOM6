! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> The "super critical" configuration
module supercritical_initialization

use MOM_dyn_horgrid,    only : dyn_horgrid_type
use MOM_error_handler,  only : MOM_mesg, MOM_error, FATAL, is_root_pe
use MOM_file_parser,    only : get_param, log_version, param_file_type
use MOM_grid,           only : ocean_grid_type
use MOM_open_boundary,  only : ocean_OBC_type, OBC_segment_type, rotate_OBC_segment_direction
use MOM_open_boundary,  only : OBC_DIRECTION_E, OBC_DIRECTION_W
use MOM_time_manager,   only : time_type
use MOM_unit_scaling,   only : unit_scale_type
use MOM_verticalGrid,   only : verticalGrid_type

implicit none ; private

#include <MOM_memory.h>

public supercritical_set_OBC_data


  interface
module subroutine supercritical_set_OBC_data(OBC, G, GV, US, param_file)
  type(ocean_OBC_type),    pointer    :: OBC  !< This open boundary condition type specifies
                                              !! whether, where, and what open boundary
                                              !! conditions are used.
  type(ocean_grid_type),   intent(in) :: G    !< The ocean's grid structure.
  type(verticalGrid_type), intent(in) :: GV   !< The ocean's vertical grid structure
  type(unit_scale_type),   intent(in) :: US   !< A dimensional unit scaling type
  type(param_file_type),   intent(in) :: param_file !< Parameter file structure
  ! Local variables

end subroutine supercritical_set_OBC_data
  end interface

end module supercritical_initialization
