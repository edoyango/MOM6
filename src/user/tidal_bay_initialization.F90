! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> Configures the model for the "tidal_bay" experiment.
!! tidal_bay = Tidally resonant bay from Zygmunt Kowalik's class on tides.
module tidal_bay_initialization

use MOM_coms,           only : reproducing_sum
use MOM_dyn_horgrid,    only : dyn_horgrid_type
use MOM_error_handler,  only : MOM_mesg, MOM_error, FATAL, WARNING, is_root_pe
use MOM_file_parser,    only : get_param, log_version, param_file_type
use MOM_grid,           only : ocean_grid_type
use MOM_open_boundary,  only : ocean_OBC_type
use MOM_open_boundary,  only : OBC_segment_type, register_OBC
use MOM_open_boundary,  only : OBC_registry_type
use MOM_unit_scaling,   only : unit_scale_type
use MOM_verticalGrid,   only : verticalGrid_type
use MOM_time_manager,   only : time_type, time_to_real

implicit none ; private

#include <MOM_memory.h>

public tidal_bay_set_OBC_data
public register_tidal_bay_OBC

!> Control structure for tidal bay open boundaries.
type, public :: tidal_bay_OBC_CS ; private
  real :: tide_flow = 3.0e6  !< Maximum tidal flux with the tidal bay configuration [L2 Z T-1 ~> m3 s-1]
  real :: tide_period        !< The period associated with the tidal bay configuration [T ~> s]
  real :: tide_ssh_amp       !< The magnitude of the sea surface height anomalies at the inflow
                             !! with the tidal bay configuration [Z ~> m]
end type tidal_bay_OBC_CS


  interface
module function register_tidal_bay_OBC(param_file, CS, US, OBC_Reg)
  type(param_file_type),    intent(in) :: param_file !< parameter file.
  type(tidal_bay_OBC_CS),   intent(inout) :: CS      !< tidal bay control structure.
  type(unit_scale_type),    intent(in) :: US         !< A dimensional unit scaling type
  type(OBC_registry_type),  pointer    :: OBC_Reg    !< OBC registry.
  logical                              :: register_tidal_bay_OBC

end function register_tidal_bay_OBC
module subroutine tidal_bay_set_OBC_data(OBC, CS, G, GV, US, h, Time)
  type(ocean_OBC_type),    pointer    :: OBC  !< This open boundary condition type specifies
                                              !! whether, where, and what open boundary
                                              !! conditions are used.
  type(tidal_bay_OBC_CS),  intent(in) :: CS   !< tidal bay control structure.
  type(ocean_grid_type),   intent(in) :: G    !< The ocean's grid structure.
  type(verticalGrid_type), intent(in) :: GV   !< The ocean's vertical grid structure
  type(unit_scale_type),   intent(in) :: US   !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in) :: h !< layer thickness [H ~> m or kg m-2]
  type(time_type),         intent(in) :: Time !< model time.

  ! The following variables are used to set up the transport in the tidal_bay example.

end subroutine tidal_bay_set_OBC_data
  end interface

end module tidal_bay_initialization
