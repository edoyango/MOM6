! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> Initialization routines for the dense water formation
!! and overflow experiment.
module dense_water_initialization

use MOM_ALE_sponge,    only : ALE_sponge_CS, set_up_ALE_sponge_field, initialize_ALE_sponge
use MOM_dyn_horgrid,   only : dyn_horgrid_type
use MOM_EOS,           only : EOS_type
use MOM_error_handler, only : MOM_error, FATAL
use MOM_file_parser,   only : get_param, param_file_type
use MOM_grid,          only : ocean_grid_type
use MOM_sponge,        only : sponge_CS
use MOM_unit_scaling,  only : unit_scale_type
use MOM_variables,     only : thermo_var_ptrs
use MOM_verticalGrid,  only : verticalGrid_type

implicit none ; private

#include <MOM_memory.h>

public dense_water_initialize_topography
public dense_water_initialize_TS
public dense_water_initialize_sponges

character(len=40) :: mdl = "dense_water_initialization" !< Module name

real, parameter :: default_sill  = 0.2  !< Default depth of the sill [nondim]
real, parameter :: default_shelf = 0.4  !< Default depth of the shelf [nondim]
real, parameter :: default_mld   = 0.25 !< Default depth of the mixed layer [nondim]


  interface
module subroutine dense_water_initialize_topography(D, G, param_file, max_depth)
  type(dyn_horgrid_type),  intent(in)  :: G !< The dynamic horizontal grid type
  real, dimension(G%isd:G%ied,G%jsd:G%jed), &
                           intent(out) :: D !< Ocean bottom depth [Z ~> m]
  type(param_file_type),   intent(in)  :: param_file !< Parameter file structure
  real,                    intent(in)  :: max_depth !< Maximum ocean depth [Z ~> m]

  ! Local variables
                        ! the basin depth [nondim]
                        ! as a fraction the basin depth [nondim]

end subroutine dense_water_initialize_topography
module subroutine dense_water_initialize_TS(G, GV, US, param_file, T, S, h, just_read)
  type(ocean_grid_type),                     intent(in)  :: G !< Horizontal grid control structure
  type(verticalGrid_type),                   intent(in)  :: GV !< Vertical grid control structure
  type(unit_scale_type),                     intent(in)  :: US !< A dimensional unit scaling type
  type(param_file_type),                     intent(in)  :: param_file !< Parameter file structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(out) :: T !< Output temperature [C ~> degC]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(out) :: S !< Output salinity [S ~> ppt]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in)  :: h !< Layer thicknesses [Z ~> m]
  logical,                                   intent(in)  :: just_read !< If true, this call will
                                                      !! only read parameters without changing T & S.
  ! Local variables

end subroutine dense_water_initialize_TS
module subroutine dense_water_initialize_sponges(G, GV, US, tv, depth_tot, param_file, use_ALE, CSp, ACSp)
  type(ocean_grid_type),   intent(in) :: G !< Horizontal grid control structure
  type(verticalGrid_type), intent(in) :: GV !< Vertical grid control structure
  type(unit_scale_type),   intent(in) :: US !< A dimensional unit scaling type
  type(thermo_var_ptrs),   intent(in) :: tv !< Thermodynamic variables
  real, dimension(SZI_(G),SZJ_(G)), &
                           intent(in) :: depth_tot  !< The nominal total depth of the ocean [Z ~> m]
  type(param_file_type),   intent(in) :: param_file !< Parameter file structure
  logical,                 intent(in) :: use_ALE !< ALE flag
  type(sponge_CS),         pointer    :: CSp !< Layered sponge control structure pointer
  type(ALE_sponge_CS),     pointer    :: ACSp !< ALE sponge control structure pointer

  ! Local variables


end subroutine dense_water_initialize_sponges
  end interface

end module dense_water_initialization
