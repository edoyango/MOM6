! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> Controls where open boundary conditions are applied
module MOM_boundary_update

use MOM_cpu_clock,             only : cpu_clock_id, cpu_clock_begin, cpu_clock_end, CLOCK_ROUTINE
use MOM_diag_mediator,         only : time_type
use MOM_error_handler,         only : MOM_mesg, MOM_error, FATAL, WARNING
use MOM_file_parser,           only : get_param, log_version, param_file_type, log_param
use MOM_grid,                  only : ocean_grid_type
use MOM_dyn_horgrid,           only : dyn_horgrid_type
use MOM_open_boundary,         only : ocean_obc_type, update_OBC_segment_data, chksum_OBC_segments
use MOM_open_boundary,         only : read_OBC_segment_data
use MOM_open_boundary,         only : OBC_registry_type, file_OBC_CS
use MOM_open_boundary,         only : register_file_OBC, file_OBC_end
use MOM_unit_scaling,          only : unit_scale_type
use MOM_tracer_registry,       only : tracer_registry_type
use MOM_variables,             only : thermo_var_ptrs
use MOM_verticalGrid,          only : verticalGrid_type
use DOME_initialization,       only : register_DOME_OBC
use tidal_bay_initialization,  only : tidal_bay_set_OBC_data, register_tidal_bay_OBC
use tidal_bay_initialization,  only : tidal_bay_OBC_CS
use Kelvin_initialization,     only : Kelvin_set_OBC_data, register_Kelvin_OBC
use Kelvin_initialization,     only : Kelvin_OBC_end, Kelvin_OBC_CS
use shelfwave_initialization,  only : shelfwave_set_OBC_data, register_shelfwave_OBC
use shelfwave_initialization,  only : shelfwave_OBC_end, shelfwave_OBC_CS
use dyed_channel_initialization, only : dyed_channel_update_flow, register_dyed_channel_OBC
use dyed_channel_initialization, only : dyed_channel_OBC_end, dyed_channel_OBC_CS

implicit none ; private

#include <MOM_memory.h>

public call_OBC_register, OBC_register_end
public update_OBC_data

!> The control structure for the MOM_boundary_update module
type, public :: update_OBC_CS ; private
  logical :: use_files = .false.        !< If true, use external files for the open boundary.
  logical :: use_Kelvin = .false.       !< If true, use the Kelvin wave open boundary.
  logical :: use_tidal_bay = .false.    !< If true, use the tidal_bay open boundary.
  logical :: use_shelfwave = .false.    !< If true, use the shelfwave open boundary.
  logical :: use_dyed_channel = .false. !< If true, use the dyed channel open boundary.
  logical :: debug_OBCs = .false.       !< If true, write verbose OBC values for debugging purposes.
  logical :: value_update_bug = .true.  !< If true, recover a bug that OBC segment data does not
                                        !! update if all segments use 'value' and none uses 'file'.
  integer :: nk_OBC_debug = 0           !< The number of layers of OBC segment data to write out
                                        !! in full when DEBUG_OBCS is true.
  !>@{ Pointers to the control structures for named OBC specifications
  type(file_OBC_CS), pointer :: file_OBC_CSp => NULL()
  type(Kelvin_OBC_CS), pointer :: Kelvin_OBC_CSp => NULL()
  type(tidal_bay_OBC_CS) :: tidal_bay_OBC
  type(shelfwave_OBC_CS), pointer :: shelfwave_OBC_CSp => NULL()
  type(dyed_channel_OBC_CS), pointer :: dyed_channel_OBC_CSp => NULL()
  !>@}
end type update_OBC_CS

integer :: id_clock_pass !< A CPU time clock ID

! character(len=40)  :: mdl = "MOM_boundary_update" ! This module's name.


  interface
module subroutine call_OBC_register(G, GV, US, param_file, CS, OBC, tr_Reg)
  type(ocean_grid_type),      intent(in) :: G    !< Ocean grid structure
  type(verticalGrid_type),    intent(in) :: GV   !< Ocean vertical grid structure
  type(unit_scale_type),      intent(in) :: US   !< A dimensional unit scaling type
  type(param_file_type),      intent(in) :: param_file !< Parameter file to parse
  type(update_OBC_CS),        pointer    :: CS         !< Control structure for OBCs
  type(ocean_OBC_type),       pointer    :: OBC        !< Open boundary structure
  type(tracer_registry_type), pointer    :: tr_Reg     !< Tracer registry.

  ! Local variables
  ! This include declares and sets the variable "version".
end subroutine call_OBC_register
module subroutine update_OBC_data(OBC, G, GV, US, tv, h, CS, Time)
  type(ocean_grid_type),                     intent(in)    :: G    !< Ocean grid structure
  type(verticalGrid_type),                   intent(in)    :: GV   !< Ocean vertical grid structure
  type(unit_scale_type),                     intent(in)    :: US   !< A dimensional unit scaling type
  type(thermo_var_ptrs),                     intent(in)    :: tv   !< Thermodynamics structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(inout) :: h    !< layer thicknesses [H ~> m or kg m-2]
  type(ocean_OBC_type),                      pointer       :: OBC  !< Open boundary structure
  type(update_OBC_CS),                       pointer       :: CS   !< Control structure for OBCs
  type(time_type),                           intent(in)    :: Time !< Model time

end subroutine update_OBC_data
module subroutine OBC_register_end(CS)
  type(update_OBC_CS),       pointer    :: CS !< Control structure for OBCs

end subroutine OBC_register_end
  end interface

end module MOM_boundary_update
