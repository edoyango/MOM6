! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> A template of a user to code up customized initial conditions.
module user_initialization

use MOM_error_handler, only : MOM_mesg, MOM_error, FATAL, is_root_pe
use MOM_dyn_horgrid, only : dyn_horgrid_type
use MOM_file_parser, only : get_param, log_version, param_file_type
use MOM_get_input, only : directories
use MOM_grid, only : ocean_grid_type
use MOM_open_boundary, only : ocean_OBC_type, OBC_NONE
use MOM_open_boundary, only : OBC_DIRECTION_E, OBC_DIRECTION_W, OBC_DIRECTION_N
use MOM_open_boundary, only : OBC_DIRECTION_S
use MOM_sponge, only : set_up_sponge_field, initialize_sponge, sponge_CS
use MOM_tracer_registry, only : tracer_registry_type
use MOM_unit_scaling, only : unit_scale_type
use MOM_variables, only : thermo_var_ptrs
use MOM_verticalGrid, only : verticalGrid_type

implicit none ; private

#include <MOM_memory.h>

public USER_set_coord, USER_initialize_topography, USER_initialize_thickness
public USER_initialize_velocity, USER_init_temperature_salinity
public USER_initialize_sponges, USER_set_OBC_data, USER_set_rotation

! A note on unit descriptions in comments: MOM6 uses units that can be rescaled for dimensional
! consistency testing. These are noted in comments with units like Z, H, L, and T, along with
! their mks counterparts with notation like "a velocity [Z T-1 ~> m s-1]".  If the units
! vary with the Boussinesq approximation, the Boussinesq variant is given first.

!> A module variable that should not be used.
!! \todo Move this module variable into a control structure.
logical :: first_call = .true.


  interface
module subroutine USER_set_coord(Rlay, g_prime, GV, US, param_file)
  type(verticalGrid_type),  intent(in)  :: GV      !< The ocean's vertical grid structure
  real, dimension(GV%ke),   intent(out) :: Rlay    !< Layer potential density [R ~> kg m-3].
  real, dimension(GV%ke+1), intent(out) :: g_prime !< The reduced gravity at each
                                                   !! interface [L2 Z-1 T-2 ~> m s-2].
  type(unit_scale_type),    intent(in)  :: US      !< A dimensional unit scaling type
  type(param_file_type),    intent(in)  :: param_file !< A structure indicating the
                                                   !! open file to parse for model
                                                   !! parameter values.

end subroutine USER_set_coord
module subroutine USER_initialize_topography(D, G, param_file, max_depth, US)
  type(dyn_horgrid_type),          intent(in)  :: G !< The dynamic horizontal grid type
  real, dimension(G%isd:G%ied,G%jsd:G%jed), &
                                   intent(out) :: D !< Ocean bottom depth [Z ~> m]
  type(param_file_type),           intent(in)  :: param_file !< Parameter file structure
  real,                            intent(in)  :: max_depth !< Maximum model depth [Z ~> m]
  type(unit_scale_type),           intent(in)  :: US !< A dimensional unit scaling type

end subroutine USER_initialize_topography
module subroutine USER_initialize_thickness(h, G, GV, param_file, just_read)
  type(ocean_grid_type),   intent(in)  :: G  !< The ocean's grid structure.
  type(verticalGrid_type), intent(in)  :: GV !< The ocean's vertical grid structure.
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(out) :: h  !< The thicknesses being initialized [Z ~> m]
  type(param_file_type),   intent(in)  :: param_file !< A structure indicating the open
                                             !! file to parse for model parameter values.
  logical,                 intent(in)  :: just_read !< If true, this call will
                                             !! only read parameters without changing h.

end subroutine USER_initialize_thickness
module subroutine USER_initialize_velocity(u, v, G, GV, US, param_file, just_read)
  type(ocean_grid_type),                       intent(in)  :: G !< Ocean grid structure.
  type(verticalGrid_type),                     intent(in)  :: GV !< The ocean's vertical grid structure.
  real, dimension(SZIB_(G), SZJ_(G),SZK_(GV)), intent(out) :: u !< i-component of velocity [L T-1 ~> m s-1]
  real, dimension(SZI_(G), SZJB_(G),SZK_(GV)), intent(out) :: v !< j-component of velocity [L T-1 ~> m s-1]
  type(unit_scale_type),                       intent(in)  :: US !< A dimensional unit scaling type
  type(param_file_type),                       intent(in)  :: param_file !< A structure indicating the
                                                            !! open file to parse for model
                                                            !! parameter values.
  logical,                                     intent(in)  :: just_read !< If true, this call will
                                                      !! only read parameters without changing u & v.

end subroutine USER_initialize_velocity
module subroutine USER_init_temperature_salinity(T, S, G, GV, param_file, just_read)
  type(ocean_grid_type),                     intent(in)  :: G !< Ocean grid structure.
  type(verticalGrid_type),                   intent(in)  :: GV !< The ocean's vertical grid structure.
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(out) :: T !< Potential temperature [C ~> degC].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(out) :: S !< Salinity [S ~> ppt].
  type(param_file_type),                     intent(in)  :: param_file !< A structure indicating the
                                                            !! open file to parse for model
                                                            !! parameter values.
  logical,                                   intent(in)  :: just_read !< If true, this call will only
                                                           !! read parameters without changing T & S.

end subroutine USER_init_temperature_salinity
module subroutine USER_initialize_sponges(G, GV, use_temp, tv, param_file, CSp, h)
  type(ocean_grid_type),   intent(in) :: G             !< Ocean grid structure.
  type(verticalGrid_type), intent(in) :: GV            !< The ocean's vertical grid structure.
  logical,                 intent(in) :: use_temp      !< If true, temperature and salinity are state variables.
  type(thermo_var_ptrs),   intent(in) :: tv            !< A structure containing pointers
                                                       !! to any available thermodynamic
                                                       !! fields, potential temperature and
                                                       !! salinity or mixed layer density.
                                                       !! Absent fields have NULL ptrs.
  type(param_file_type),   intent(in) :: param_file    !< A structure indicating the
                                                       !! open file to parse for model
                                                       !! parameter values.
  type(sponge_CS),         pointer    :: CSp           !< A pointer to the sponge control structure.
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in) :: h             !< Layer thicknesses [H ~> m or kg m-2].
end subroutine USER_initialize_sponges
module subroutine USER_set_OBC_data(OBC, tv, G, GV, param_file, tr_Reg)
  type(ocean_OBC_type),       pointer    :: OBC   !< This open boundary condition type specifies
                                                  !! whether, where, and what open boundary
                                                  !! conditions are used.
  type(thermo_var_ptrs),      intent(in) :: tv    !< A structure containing pointers to any
                                       !! available thermodynamic fields, including potential
                                       !! temperature and salinity or mixed layer density. Absent
                                       !! fields have NULL ptrs.
  type(ocean_grid_type),      intent(in) :: G     !< The ocean's grid structure.
  type(verticalGrid_type),    intent(in) :: GV    !< The ocean's vertical grid structure.
  type(param_file_type),      intent(in) :: param_file !< A structure indicating the
                                                  !! open file to parse for model
                                                  !! parameter values.
  type(tracer_registry_type), pointer    :: tr_Reg !< Tracer registry.
!  call MOM_error(FATAL, &
!   "USER_initialization.F90, USER_set_OBC_data: " // &
!   "Unmodified user routine called - you must edit the routine to use it")

end subroutine USER_set_OBC_data
module subroutine USER_set_rotation(G, param_file)
  type(ocean_grid_type), intent(inout) :: G    !< The ocean's grid structure
  type(param_file_type), intent(in)    :: param_file !< A structure to parse for run-time parameters
end subroutine USER_set_rotation
module subroutine write_user_log(param_file)
  type(param_file_type), intent(in) :: param_file !< A structure indicating the
                                                  !! open file to parse for model
                                                  !! parameter values.

  ! This include declares and sets the variable "version".

end subroutine write_user_log
  end interface

end module user_initialization
