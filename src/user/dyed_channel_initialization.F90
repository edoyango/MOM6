! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> Initialization for the dyed_channel configuration
module dyed_channel_initialization

use MOM_dyn_horgrid,     only : dyn_horgrid_type
use MOM_error_handler,   only : MOM_mesg, MOM_error, FATAL, WARNING, is_root_pe
use MOM_file_parser,     only : get_param, log_version, param_file_type
use MOM_get_input,       only : directories
use MOM_grid,            only : ocean_grid_type
use MOM_open_boundary,   only : ocean_OBC_type, OBC_NONE
use MOM_open_boundary,   only : OBC_DIRECTION_W, OBC_DIRECTION_N, OBC_DIRECTION_S, OBC_DIRECTION_E
use MOM_open_boundary,   only : OBC_segment_type, register_segment_tracer
use MOM_open_boundary,   only : OBC_registry_type, register_OBC
use MOM_time_manager,    only : time_type, time_to_real
use MOM_tracer_registry, only : tracer_registry_type, tracer_name_lookup
use MOM_tracer_registry, only : tracer_type
use MOM_unit_scaling,    only : unit_scale_type
use MOM_variables,       only : thermo_var_ptrs
use MOM_verticalGrid,    only : verticalGrid_type

implicit none ; private

#include <MOM_memory.h>

public dyed_channel_set_OBC_tracer_data, dyed_channel_OBC_end
public register_dyed_channel_OBC, dyed_channel_update_flow

!> Control structure for dyed-channel open boundaries.
type, public :: dyed_channel_OBC_CS ; private
  real :: zonal_flow = 8.57         !< Mean inflow [L T-1 ~> m s-1]
  real :: tidal_amp = 0.0           !< Sloshing amplitude [L T-1 ~> m s-1]
  real :: frequency  = 0.0          !< Sloshing frequency [T-1 ~> s-1]
  logical :: OBC_transport_bug      !< If true and specified open boundary conditions are being
                                    !! used, use a 1 m (if Boussienesq) or 1 kg m-2 layer thickness
                                    !! instead of the actual thickness.
end type dyed_channel_OBC_CS

integer :: ntr = 0 !< Number of dye tracers
                   !! \todo This is a module variable. Move this variable into the control structure.


  interface
logical module function register_dyed_channel_OBC(param_file, CS, US, OBC_Reg)
  type(param_file_type),     intent(in) :: param_file !< parameter file.
  type(dyed_channel_OBC_CS), pointer    :: CS         !< Dyed channel control structure.
  type(unit_scale_type),     intent(in) :: US         !< A dimensional unit scaling type
  type(OBC_registry_type),   pointer    :: OBC_Reg    !< OBC registry.

  ! Local variables
                          ! recreate the bugs, or if false bugs are only used if actively selected.

end function register_dyed_channel_OBC
module subroutine dyed_channel_OBC_end(CS)
  type(dyed_channel_OBC_CS), pointer :: CS    !< Dyed channel control structure.

end subroutine dyed_channel_OBC_end
module subroutine dyed_channel_set_OBC_tracer_data(OBC, G, GV, param_file, tr_Reg)
  type(ocean_OBC_type),       pointer    :: OBC !< This open boundary condition type specifies
                                                !! whether, where, and what open boundary
                                                !! conditions are used.
  type(ocean_grid_type),      intent(in) :: G   !< The ocean's grid structure.
  type(verticalGrid_type),    intent(in) :: GV  !< The ocean's vertical grid structure.
  type(param_file_type),      intent(in) :: param_file !< A structure indicating the open file
                                                !! to parse for model parameter values.
  type(tracer_registry_type), pointer    :: tr_Reg !< Tracer registry.
  ! Local variables

end subroutine dyed_channel_set_OBC_tracer_data
module subroutine dyed_channel_update_flow(OBC, CS, G, GV, US, h, Time)
  type(ocean_OBC_type),       pointer    :: OBC !< This open boundary condition type specifies
                                                !! whether, where, and what open boundary
                                                !! conditions are used.
  type(dyed_channel_OBC_CS),  pointer    :: CS  !< Dyed channel control structure.
  type(ocean_grid_type),      intent(in) :: G   !< The ocean's grid structure.
  type(verticalGrid_type),    intent(in) :: GV  !< The ocean's vertical grid structure.
  type(unit_scale_type),      intent(in) :: US  !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in) :: h !< layer thickness [H ~> m or kg m-2]
  type(time_type),            intent(in) :: Time !< model time.

  ! Local variables
                    ! reproduce a bug with the older versions of this code [H ~> m or kg m-2]

end subroutine dyed_channel_update_flow
  end interface

end module dyed_channel_initialization
