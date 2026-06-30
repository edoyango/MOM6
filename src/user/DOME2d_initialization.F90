! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> Initialization of the 2D DOME experiment with density water initialized on a coastal shelf.
module DOME2d_initialization

use MOM_ALE_sponge, only : ALE_sponge_CS, set_up_ALE_sponge_field, initialize_ALE_sponge
use MOM_dyn_horgrid, only : dyn_horgrid_type
use MOM_error_handler, only : MOM_mesg, MOM_error, FATAL
use MOM_file_parser, only : get_param, log_version, param_file_type
use MOM_get_input, only : directories
use MOM_grid, only : ocean_grid_type
use MOM_sponge, only : sponge_CS, set_up_sponge_field, initialize_sponge
use MOM_unit_scaling, only : unit_scale_type
use MOM_variables, only : thermo_var_ptrs
use MOM_verticalGrid, only : verticalGrid_type
use regrid_consts, only : coordinateMode, DEFAULT_COORDINATE_MODE
use regrid_consts, only : REGRIDDING_LAYER, REGRIDDING_ZSTAR
use regrid_consts, only : REGRIDDING_RHO, REGRIDDING_SIGMA

implicit none ; private

#include <MOM_memory.h>

! Public functions
public DOME2d_initialize_topography
public DOME2d_initialize_thickness
public DOME2d_initialize_temperature_salinity
public DOME2d_initialize_sponges

! A note on unit descriptions in comments: MOM6 uses units that can be rescaled for dimensional
! consistency testing. These are noted in comments with units like Z, H, L, and T, along with
! their mks counterparts with notation like "a velocity [Z T-1 ~> m s-1]".  If the units
! vary with the Boussinesq approximation, the Boussinesq variant is given first.

character(len=40) :: mdl = "DOME2D_initialization" !< This module's name.


  interface
module subroutine DOME2d_initialize_topography( D, G, param_file, max_depth )
  type(dyn_horgrid_type),  intent(in)  :: G !< The dynamic horizontal grid type
  real, dimension(G%isd:G%ied,G%jsd:G%jed), &
                           intent(out) :: D !< Ocean bottom depth [Z ~> m]
  type(param_file_type),   intent(in)  :: param_file !< Parameter file structure
  real,                    intent(in)  :: max_depth !< Maximum ocean depth [Z ~> m]

  ! Local variables
  ! This include declares and sets the variable "version".

end subroutine DOME2d_initialize_topography
module subroutine DOME2d_initialize_thickness ( h, depth_tot, G, GV, US, param_file, just_read )
  type(ocean_grid_type),   intent(in)  :: G  !< Ocean grid structure
  type(verticalGrid_type), intent(in)  :: GV !< Vertical grid structure
  type(unit_scale_type),   intent(in)  :: US !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(out) :: h           !< The thickness that is being initialized [Z ~> m]
  real, dimension(SZI_(G),SZJ_(G)), &
                           intent(in)  :: depth_tot   !< The nominal total depth of the ocean [Z ~> m]
  type(param_file_type),   intent(in)  :: param_file  !< A structure indicating the open file
                                                      !! to parse for model parameter values.
  logical,                 intent(in)  :: just_read   !< If true, this call will only read
                                                      !! parameters without changing h.

  ! Local variables
                           ! negative because it is positive upward.
                           ! positive upward, in depth units [Z ~> m]

end subroutine DOME2d_initialize_thickness
module subroutine DOME2d_initialize_temperature_salinity ( T, S, h, G, GV, US, param_file, just_read)
  type(ocean_grid_type),                     intent(in)  :: G  !< Ocean grid structure
  type(verticalGrid_type),                   intent(in)  :: GV !< The ocean's vertical grid structure.
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(out) :: T  !< Potential temperature [C ~> degC]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(out) :: S  !< Salinity [S ~> ppt]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in)  :: h  !< Layer thickness [Z ~> m]
  type(unit_scale_type),                     intent(in)  :: US !< A dimensional unit scaling type
  type(param_file_type),                     intent(in)  :: param_file !< Parameter file structure
  logical,                                   intent(in)  :: just_read !< If true, this call will
                                                      !! only read parameters without changing T & S.


end subroutine DOME2d_initialize_temperature_salinity
module subroutine DOME2d_initialize_sponges(G, GV, US, tv, depth_tot, param_file, use_ALE, CSp, ACSp)
  type(ocean_grid_type),   intent(in) :: G  !< Ocean grid structure
  type(verticalGrid_type), intent(in) :: GV !< Vertical grid structure
  type(unit_scale_type),   intent(in) :: US !< A dimensional unit scaling type
  type(thermo_var_ptrs),   intent(in) :: tv !< Thermodynamics structure
  real, dimension(SZI_(G),SZJ_(G)), &
                           intent(in) :: depth_tot  !< The nominal total depth of the ocean [Z ~> m]
  type(param_file_type),   intent(in) :: param_file !< Parameter file structure
  logical,                 intent(in) :: use_ALE !< If true, indicates model is in ALE mode
  type(sponge_CS),         pointer    :: CSp !< Layer-mode sponge structure
  type(ALE_sponge_CS),     pointer    :: ACSp !< ALE-mode sponge structure
  ! Local variables
                                    ! usually negative because it is positive upward.
                                    ! positive upward [Z ~> m].
                                    ! restoring T/S is active [nondim]
                                    ! restoring T/S is active [nondim]

end subroutine DOME2d_initialize_sponges
  end interface

end module DOME2d_initialization
