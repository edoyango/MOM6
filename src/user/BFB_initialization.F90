! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> Initialization of the boundary-forced-basing configuration
module BFB_initialization

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

public BFB_set_coord
public BFB_initialize_sponges_southonly

! A note on unit descriptions in comments: MOM6 uses units that can be rescaled for dimensional
! consistency testing. These are noted in comments with units like Z, H, L, and T, along with
! their mks counterparts with notation like "a velocity [Z T-1 ~> m s-1]".  If the units
! vary with the Boussinesq approximation, the Boussinesq variant is given first.


  interface
module subroutine BFB_set_coord(Rlay, g_prime, GV, US, param_file)
  type(verticalGrid_type),  intent(in)  :: GV      !< The ocean's vertical grid structure
  real, dimension(GV%ke),   intent(out) :: Rlay    !< Layer potential density [R ~> kg m-3].
  real, dimension(GV%ke+1), intent(out) :: g_prime !< The reduced gravity at each
                                                   !! interface [L2 Z-1 T-2 ~> m s-2].
  type(unit_scale_type),    intent(in)  :: US      !< A dimensional unit scaling type
  type(param_file_type),    intent(in)  :: param_file !< A structure to parse for run-time parameters

  ! This include declares and sets the variable "version".

end subroutine BFB_set_coord
module subroutine BFB_initialize_sponges_southonly(G, GV, US, use_temperature, tv, depth_tot, param_file, CSp, h)
  type(ocean_grid_type),   intent(in) :: G  !< The ocean's grid structure
  type(verticalGrid_type), intent(in) :: GV !< The ocean's vertical grid structure.
  type(unit_scale_type),   intent(in) :: US !< A dimensional unit scaling type
  logical,                 intent(in) :: use_temperature !< If true, temperature and salinity are used as
                                            !! state variables.
  type(thermo_var_ptrs),   intent(in) :: tv   !< A structure pointing to various thermodynamic variables
  real, dimension(SZI_(G),SZJ_(G)), &
                           intent(in) :: depth_tot !< The nominal total depth of the ocean [Z ~> m]
  type(param_file_type),   intent(in) :: param_file !< A structure to parse for run-time parameters
  type(sponge_CS),         pointer    :: CSp  !< A pointer to the sponge control structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in) :: h    !< Layer thicknesses [H ~> m or kg m-2]

  ! Local variables
  ! This include declares and sets the variable "version".

end subroutine BFB_initialize_sponges_southonly
  end interface

end module BFB_initialization
