! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> Initialization for the "sloshing" internal waves configuration.
module sloshing_initialization

use MOM_domains, only : sum_across_PEs
use MOM_dyn_horgrid, only : dyn_horgrid_type
use MOM_error_handler, only : MOM_mesg, MOM_error, FATAL, is_root_pe
use MOM_file_parser, only : get_param, log_version, param_file_type
use MOM_get_input, only : directories
use MOM_grid, only : ocean_grid_type
use MOM_sponge, only : set_up_sponge_field, initialize_sponge, sponge_CS
use MOM_tracer_registry, only : tracer_registry_type
use MOM_unit_scaling, only : unit_scale_type
use MOM_variables, only : thermo_var_ptrs
use MOM_verticalGrid, only : verticalGrid_type

implicit none ; private

#include <MOM_memory.h>

! The following routines are visible to the outside world
public sloshing_initialize_topography
public sloshing_initialize_thickness
public sloshing_initialize_temperature_salinity


  interface
module subroutine sloshing_initialize_topography( D, G, param_file, max_depth )
  type(dyn_horgrid_type),  intent(in)  :: G !< The dynamic horizontal grid type
  real, dimension(G%isd:G%ied,G%jsd:G%jed), &
                           intent(out) :: D !< Ocean bottom depth [Z ~> m]
  type(param_file_type),   intent(in)  :: param_file !< Parameter file structure
  real,                    intent(in)  :: max_depth !< Maximum ocean depth [Z ~> m]

  ! Local variables

end subroutine sloshing_initialize_topography
module subroutine sloshing_initialize_thickness ( h, depth_tot, G, GV, US, param_file, just_read)
  type(ocean_grid_type),   intent(in)  :: G           !< The ocean's grid structure.
  type(verticalGrid_type), intent(in)  :: GV          !< The ocean's vertical grid structure.
  type(unit_scale_type),   intent(in)  :: US          !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(out) :: h           !< The thickness that is being initialized [Z ~> m]
  real, dimension(SZI_(G),SZJ_(G)), &
                           intent(in)  :: depth_tot   !< The nominal total depth of the ocean [Z ~> m]
  type(param_file_type),   intent(in)  :: param_file  !< A structure to parse for model parameter values.
  logical,                 intent(in)  :: just_read   !< If true, this call will only read
                                                      !! parameters without changing h.

  ! Local variables
  ! This include declares and sets the variable "version".

end subroutine sloshing_initialize_thickness
module subroutine sloshing_initialize_temperature_salinity ( T, S, h, G, GV, US, param_file, just_read)
  type(ocean_grid_type),                     intent(in)  :: G  !< Ocean grid structure.
  type(verticalGrid_type),                   intent(in)  :: GV !< The ocean's vertical grid structure.
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(out) :: T  !< Potential temperature [C ~> degC].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(out) :: S  !< Salinity [S ~> ppt].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in)  :: h  !< Layer thickness [Z ~> m].
  type(unit_scale_type),                     intent(in)  :: US !< A dimensional unit scaling type
  type(param_file_type),                     intent(in)  :: param_file !< A structure to parse
                                                               !! for model parameter values.
  logical,                                   intent(in)  :: just_read !< If true, this call will only read
                                                               !! parameters without changing T & S.

  ! Local variables

end subroutine sloshing_initialize_temperature_salinity
  end interface

end module sloshing_initialization
