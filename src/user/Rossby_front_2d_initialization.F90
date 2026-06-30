! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> Initial conditions for the 2D Rossby front test
module Rossby_front_2d_initialization

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
character(len=40) :: mdl = "Rossby_front_2d_initialization" !< This module's name.
! This include declares and sets the variable "version".
#include "version_variable.h"

public Rossby_front_initialize_thickness
public Rossby_front_initialize_temperature_salinity
public Rossby_front_initialize_velocity

! Parameters defining the initial conditions of this test case
real, parameter :: frontFractionalWidth = 0.5 !< Width of front as fraction of domain [nondim]
real, parameter :: HMLmin = 0.25 !< Shallowest ML as fractional depth of ocean [nondim]
real, parameter :: HMLmax = 0.75 !< Deepest ML as fractional depth of ocean [nondim]


  interface
module subroutine Rossby_front_initialize_thickness(h, G, GV, US, param_file, just_read)
  type(ocean_grid_type),   intent(in)  :: G           !< Grid structure
  type(verticalGrid_type), intent(in)  :: GV          !< Vertical grid structure
  type(unit_scale_type),   intent(in)  :: US          !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(out) :: h           !< The thickness that is being initialized [H ~> m or kg m-2]
  type(param_file_type),   intent(in)  :: param_file  !< A structure indicating the open file
                                                      !! to parse for model parameter values.
  logical,                 intent(in)  :: just_read   !< If true, this call will only read
                                                      !! parameters without changing h.

  ! Local variables

end subroutine Rossby_front_initialize_thickness
module subroutine Rossby_front_initialize_temperature_salinity(T, S, h, G, GV, US, &
                   param_file, just_read)
  type(ocean_grid_type),                     intent(in)  :: G  !< Grid structure
  type(verticalGrid_type),                   intent(in)  :: GV !< The ocean's vertical grid structure.
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(out) :: T  !< Potential temperature [C ~> degC]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(out) :: S  !< Salinity [S ~> ppt]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in)  :: h  !< Thickness [H ~> m or kg m-2]
  type(unit_scale_type),                     intent(in)  :: US !< A dimensional unit scaling type
  type(param_file_type),                     intent(in)  :: param_file   !< Parameter file handle
  logical,                                   intent(in)  :: just_read !< If true, this call will
                                                      !! only read parameters without changing T & S.
  ! Local variables

end subroutine Rossby_front_initialize_temperature_salinity
module subroutine Rossby_front_initialize_velocity(u, v, h, G, GV, US, param_file, just_read)
  type(ocean_grid_type),      intent(in)  :: G  !< Grid structure
  type(verticalGrid_type),    intent(in)  :: GV !< Vertical grid structure
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                              intent(out) :: u  !< i-component of velocity [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                              intent(out) :: v  !< j-component of velocity [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJ_(G), SZK_(GV)), &
                              intent(in)  :: h  !< Thickness [H ~> m or kg m-2]
  type(unit_scale_type),      intent(in)  :: US !< A dimensional unit scaling type
  type(param_file_type),      intent(in)  :: param_file !< A structure indicating the open file
                                                !! to parse for model parameter values.
  logical,                    intent(in)  :: just_read !< If present and true, this call will only
                                                !! read parameters without setting u & v.

                          ! [L2 H-1 T-1 C-1 ~> m s-1 degC-1 or m4 kg-1 s-1 degC-1]

end subroutine Rossby_front_initialize_velocity
real module function yPseudo( G, lat )
  type(ocean_grid_type), intent(in) :: G   !< Grid structure
  real,                  intent(in) :: lat !< Latitude in arbitrary units, often [km]
  ! Local

end function yPseudo
real module function Hml( G, lat, max_depth )
  type(ocean_grid_type), intent(in) :: G   !< Grid structure
  real,                  intent(in) :: lat !< Latitude in arbitrary units, often [km]
  real,                  intent(in) :: max_depth !< The maximum depth of the ocean [Z ~> m] or [H ~> m or kg m-2]
  ! Local

end function Hml
real module function dTdy( G, dT, lat, US )
  type(ocean_grid_type), intent(in) :: G     !< Grid structure
  real,                  intent(in) :: dT    !< Top to bottom temperature difference [C ~> degC]
  real,                  intent(in) :: lat   !< Latitude in the same units as geoLat, often [km]
  type(unit_scale_type), intent(in) :: US    !< A dimensional unit scaling type
  ! Local

end function dTdy
  end interface

end module Rossby_front_2d_initialization
