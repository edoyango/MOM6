! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> Configures the model for the geostrophic adjustment test case.
module adjustment_initialization

use MOM_error_handler, only : MOM_mesg, MOM_error, FATAL, is_root_pe
use MOM_file_parser,   only : get_param, log_param, log_version, param_file_type
use MOM_get_input,     only : directories
use MOM_grid,          only : ocean_grid_type
use MOM_unit_scaling,  only : unit_scale_type
use MOM_variables,     only : thermo_var_ptrs
use MOM_verticalGrid,  only : verticalGrid_type
use regrid_consts,     only : coordinateMode, DEFAULT_COORDINATE_MODE
use regrid_consts,     only : REGRIDDING_LAYER, REGRIDDING_ZSTAR
use regrid_consts,     only : REGRIDDING_RHO, REGRIDDING_SIGMA

implicit none ; private

character(len=40) :: mdl = "adjustment_initialization" !< This module's name.

#include <MOM_memory.h>

public adjustment_initialize_thickness
public adjustment_initialize_temperature_salinity

! A note on unit descriptions in comments: MOM6 uses units that can be rescaled for dimensional
! consistency testing. These are noted in comments with units like Z, H, L, and T, along with
! their mks counterparts with notation like "a velocity [Z T-1 ~> m s-1]".  If the units
! vary with the Boussinesq approximation, the Boussinesq variant is given first.


  interface
module subroutine adjustment_initialize_thickness ( h, G, GV, US, param_file, just_read)
  type(ocean_grid_type),   intent(in)  :: G           !< The ocean's grid structure.
  type(verticalGrid_type), intent(in)  :: GV          !< The ocean's vertical grid structure.
  type(unit_scale_type),   intent(in)  :: US          !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(out) :: h           !< The thickness that is being initialized [Z ~> m]
  type(param_file_type),   intent(in)  :: param_file  !< A structure indicating the open file
                                                      !! to parse for model parameter values.
  logical,                 intent(in)  :: just_read   !< If true, this call will only read
                                                      !! parameters without changing h.
  ! Local variables
                            ! negative because it is positive upward.
                            ! positive upward, in depth units [Z ~> m].
                            ! In this subroutine it is hard coded at 1.0 kg m-3 ppt-1.
  ! This include declares and sets the variable "version".

end subroutine adjustment_initialize_thickness
module subroutine adjustment_initialize_temperature_salinity(T, S, h, depth_tot, G, GV, US, param_file, just_read)
  type(ocean_grid_type),   intent(in)  :: G           !< The ocean's grid structure.
  type(verticalGrid_type), intent(in)  :: GV          !< The ocean's vertical grid structure.
  type(unit_scale_type),   intent(in)  :: US          !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(out) :: T           !< The temperature that is being initialized [C ~> degC]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(out) :: S           !< The salinity that is being initialized [S ~> ppt]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in)  :: h           !< The model thicknesses [Z ~> m]
  real, dimension(SZI_(G),SZJ_(G)), &
                           intent(in)  :: depth_tot   !< The nominal total depth of the ocean [Z ~> m]
  type(param_file_type),   intent(in)  :: param_file  !< A structure indicating the open file to
                                                      !! parse for model parameter values.
  logical,                 intent(in)  :: just_read   !< If true, this call will only read
                                                      !! parameters without changing T & S.


end subroutine adjustment_initialize_temperature_salinity
  end interface

end module adjustment_initialization
