! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> Initial conditions for an idealized baroclinic zone
module baroclinic_zone_initialization

use MOM_file_parser, only : get_param, log_version, param_file_type
use MOM_file_parser,   only : openParameterBlock, closeParameterBlock
use MOM_grid, only : ocean_grid_type
use MOM_unit_scaling, only : unit_scale_type
use MOM_verticalGrid, only : verticalGrid_type

implicit none ; private

#include <MOM_memory.h>
#include "version_variable.h"

! Private (module-wise) parameters
character(len=40) :: mdl = "baroclinic_zone_initialization" !< This module's name.

public baroclinic_zone_init_temperature_salinity

! A note on unit descriptions in comments: MOM6 uses units that can be rescaled for dimensional
! consistency testing. These are noted in comments with units like Z, H, L, and T, along with
! their mks counterparts with notation like "a velocity [Z T-1 ~> m s-1]".  If the units
! vary with the Boussinesq approximation, the Boussinesq variant is given first.


  interface
module subroutine bcz_params(G, GV, US, param_file, S_ref, dSdz, delta_S, dSdx, T_ref, dTdz, &
                      delta_T, dTdx, L_zone, just_read)
  type(ocean_grid_type),   intent(in)  :: G          !< Grid structure
  type(verticalGrid_type), intent(in)  :: GV         !< The ocean's vertical grid structure.
  type(unit_scale_type),   intent(in)  :: US         !< A dimensional unit scaling type
  type(param_file_type),   intent(in)  :: param_file !< Parameter file handle
  real,                    intent(out) :: S_ref      !< Reference salinity [S ~> ppt]
  real,                    intent(out) :: dSdz       !< Salinity stratification [S Z-1 ~> ppt m-1]
  real,                    intent(out) :: delta_S    !< Salinity difference across baroclinic zone [S ~> ppt]
  real,                    intent(out) :: dSdx       !< Linear salinity gradient, often in [S km-1 ~> ppt km-1]
                                                     !! or [S degrees_E-1 ~> ppt degrees_E-1], depending on
                                                     !! the value of G%x_axis_units
  real,                    intent(out) :: T_ref      !< Reference temperature [C ~> degC]
  real,                    intent(out) :: dTdz       !< Temperature stratification [C Z-1 ~> degC m-1]
  real,                    intent(out) :: delta_T    !< Temperature difference across baroclinic zone [C ~> degC]
  real,                    intent(out) :: dTdx       !< Linear temperature gradient, often in [C km-1 ~> degC km-1]
                                                     !! or [C degrees_E-1 ~> degC degrees_E-1], depending on
                                                     !! the value of G%x_axis_units
  real,                    intent(out) :: L_zone     !< Width of baroclinic zone, often in [km] or [degrees_N],
                                                     !! depending on the value of G%y_axis_units
  logical,                 intent(in)  :: just_read  !< If true, this call will
                                                     !! only read parameters without changing h.

end subroutine bcz_params
module subroutine baroclinic_zone_init_temperature_salinity(T, S, h, depth_tot, G, GV, US, param_file, &
                                                     just_read)
  type(ocean_grid_type),   intent(in)  :: G          !< Grid structure
  type(verticalGrid_type), intent(in)  :: GV         !< The ocean's vertical grid structure.
  type(unit_scale_type),   intent(in)  :: US         !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(out) :: T          !< Potential temperature [C ~> degC]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(out) :: S          !< Salinity [S ~> ppt]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in)  :: h          !< The model thicknesses [Z ~> m]
  real, dimension(SZI_(G),SZJ_(G)), &
                           intent(in)  :: depth_tot  !< The nominal total depth of the ocean [Z ~> m]
  type(param_file_type),   intent(in)  :: param_file !< A structure indicating the open file
                                                     !! to parse for model parameter values.
  logical,                 intent(in)  :: just_read  !< If true, this call will only read
                                                     !! parameters without changing T & S.

                              ! on the value of G%y_axis_units

end subroutine baroclinic_zone_init_temperature_salinity
  end interface

end module baroclinic_zone_initialization
