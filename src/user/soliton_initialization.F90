! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> Initial conditions for the Equatorial Rossby soliton test (Boyd).
module soliton_initialization

use MOM_error_handler, only : MOM_mesg, MOM_error, FATAL, is_root_pe
use MOM_file_parser, only : get_param, log_version, param_file_type
use MOM_get_input, only : directories
use MOM_grid, only : ocean_grid_type
use MOM_unit_scaling, only : unit_scale_type
use MOM_variables, only : thermo_var_ptrs
use MOM_verticalGrid, only : verticalGrid_type
use regrid_consts, only : coordinateMode, DEFAULT_COORDINATE_MODE
use regrid_consts, only : REGRIDDING_LAYER, REGRIDDING_ZSTAR
use regrid_consts, only : REGRIDDING_RHO, REGRIDDING_SIGMA

implicit none ; private

#include <MOM_memory.h>

! Private (module-wise) parameters
character(len=40) :: mdl = "soliton_initialization" !< This module's name.

public soliton_initialize_thickness
public soliton_initialize_velocity


  interface
module subroutine soliton_initialize_thickness(h, depth_tot, G, GV, US, param_file, just_read)
  type(ocean_grid_type),   intent(in)  :: G    !< The ocean's grid structure.
  type(verticalGrid_type), intent(in)  :: GV   !< The ocean's vertical grid structure.
  type(unit_scale_type),   intent(in)  :: US   !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(out) :: h    !< The thickness that is being initialized [Z ~> m]
  real, dimension(SZI_(G),SZJ_(G)), &
                           intent(in)  :: depth_tot !< The nominal total depth of the ocean [Z ~> m]
  type(param_file_type),   intent(in)  :: param_file  !< A structure indicating the open file
                                                      !! to parse for model parameter values.
  logical,                 intent(in)  :: just_read   !< If true, this call will only read
                                                      !! parameters without changing h.

  ! Local variables
  ! This include declares and sets the variable "version".

end subroutine soliton_initialize_thickness
module subroutine soliton_initialize_velocity(u, v, G, GV, US, param_file, just_read)
  type(ocean_grid_type),                      intent(in)  :: G  !< Grid structure
  type(verticalGrid_type),                    intent(in)  :: GV !< The ocean's vertical grid structure
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), intent(out) :: u  !< i-component of velocity [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), intent(out) :: v  !< j-component of velocity [L T-1 ~> m s-1]
  type(unit_scale_type),                      intent(in)  :: US !< A dimensional unit scaling type
  type(param_file_type),   intent(in)  :: param_file  !< A structure indicating the open file
                                                      !! to parse for model parameter values.
  logical,                 intent(in)  :: just_read   !< If true, this call will only read
                                                      !! parameters without changing h.

  ! Local variables

end subroutine soliton_initialize_velocity
  end interface

end module soliton_initialization
