! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> Configures the model for the idealized shelfwave test case.
module shelfwave_initialization

use MOM_domains,        only : sum_across_PEs
use MOM_dyn_horgrid,    only : dyn_horgrid_type
use MOM_error_handler,  only : MOM_mesg, MOM_error, FATAL, WARNING, is_root_pe
use MOM_file_parser,    only : get_param, log_version, param_file_type
use MOM_grid,           only : ocean_grid_type
use MOM_open_boundary,  only : ocean_OBC_type, OBC_NONE, OBC_DIRECTION_W
use MOM_open_boundary,  only : OBC_segment_type, register_OBC
use MOM_open_boundary,  only : OBC_registry_type, rotate_OBC_segment_direction
use MOM_time_manager,   only : time_type, time_to_real
use MOM_unit_scaling,   only : unit_scale_type
use MOM_verticalGrid,   only : verticalGrid_type

implicit none ; private

#include <MOM_memory.h>

character(len=40) :: mdl = "shelfwave_initialization" !< This module's name.

! The following routines are visible to the outside world
public shelfwave_initialize_topography
public shelfwave_set_OBC_data
public register_shelfwave_OBC, shelfwave_OBC_end

!> Control structure for shelfwave open boundaries.
type, public :: shelfwave_OBC_CS ; private
  real :: my_amp        !< Amplitude of the open boundary current inflows [L T-1 ~> m s-1]
  real :: kk            !< Cross-shore wavenumber [km-1] or [m-1]
  real :: ll            !< Longshore wavenumber [km-1] or [m-1]
  real :: alpha         !< Exponential decay rate in the y-direction [km-1] or [m-1]
  real :: omega         !< Frequency of the shelf wave [T-1 ~> s-1]
  logical :: shelfwave_correct_amplitude !< If true, SHELFWAVE_AMPLITUDE gives the actual inflow
                        !! velocity, rather than giving an overall scaling factor for the flow.
end type shelfwave_OBC_CS


  interface
module function register_shelfwave_OBC(param_file, CS, G, US, OBC_Reg)
  type(param_file_type),    intent(in) :: param_file !< parameter file.
  type(shelfwave_OBC_CS),   pointer    :: CS         !< shelfwave control structure.
  type(ocean_grid_type),    intent(in) :: G          !< The ocean's grid structure.
  type(unit_scale_type),    intent(in) :: US         !< A dimensional unit scaling type
  type(OBC_registry_type),  pointer    :: OBC_Reg    !< Open boundary condition registry.
  logical                              :: register_shelfwave_OBC

  ! Local variables

end function register_shelfwave_OBC
module subroutine shelfwave_OBC_end(CS)
  type(shelfwave_OBC_CS), pointer    :: CS         !< shelfwave control structure.

end subroutine shelfwave_OBC_end
module subroutine shelfwave_initialize_topography( D, G, param_file, max_depth, US )
  type(dyn_horgrid_type),          intent(in)  :: G !< The dynamic horizontal grid type
  real, dimension(G%isd:G%ied,G%jsd:G%jed), &
                                   intent(out) :: D !< Ocean bottom depth [Z ~> m]
  type(param_file_type),           intent(in)  :: param_file !< Parameter file structure
  real,                            intent(in)  :: max_depth !< Maximum model depth [Z ~> m]
  type(unit_scale_type),           intent(in)  :: US !< A dimensional unit scaling type

  ! Local variables

end subroutine shelfwave_initialize_topography
module subroutine shelfwave_set_OBC_data(OBC, CS, G, GV, US, h, Time)
  type(ocean_OBC_type),    pointer    :: OBC  !< This open boundary condition type specifies
                                              !! whether, where, and what open boundary
                                              !! conditions are used.
  type(shelfwave_OBC_CS),  pointer    :: CS   !< tidal bay control structure.
  type(ocean_grid_type),   intent(in) :: G    !< The ocean's grid structure.
  type(verticalGrid_type), intent(in) :: GV   !< The ocean's vertical grid structure
  type(unit_scale_type),   intent(in) :: US   !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in) :: h !< layer thickness [H ~> m or kg m-2]
  type(time_type),         intent(in) :: Time !< model time.

  ! The following variables are used to set up the transport in the shelfwave example.
                    ! to compensate for the variable units of the y-coordinate [km axis_unit-1], usually 1 [nondim]
                    ! to account for grid rotation [L T-1 ~> m s-1]

end subroutine shelfwave_set_OBC_data
  end interface

end module shelfwave_initialization
