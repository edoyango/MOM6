! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> Initializes fixed aspects of the related to its vertical coordinate.
module MOM_coord_initialization

use MOM_debugging,        only : chksum
use MOM_EOS,              only : calculate_density, EOS_type
use MOM_error_handler,    only : MOM_mesg, MOM_error, FATAL, WARNING, is_root_pe
use MOM_error_handler,    only : callTree_enter, callTree_leave, callTree_waypoint
use MOM_file_parser,      only : get_param, read_param, log_param, param_file_type, log_version
use MOM_io,               only : create_MOM_file, file_exists
use MOM_io,               only : MOM_netCDF_file, MOM_field
use MOM_io,               only : MOM_read_data, MOM_write_field, vardesc, var_desc, SINGLE_FILE
use MOM_string_functions, only : slasher, uppercase
use MOM_unit_scaling,     only : unit_scale_type
use MOM_variables,        only : thermo_var_ptrs
use MOM_verticalGrid,     only : verticalGrid_type, setVerticalGridAxes
use user_initialization,  only : user_set_coord
use BFB_initialization,   only : BFB_set_coord

implicit none ; private

public MOM_initialize_coord, write_vertgrid_file

! A note on unit descriptions in comments: MOM6 uses units that can be rescaled for dimensional
! consistency testing. These are noted in comments with units like Z, H, L, and T, along with
! their mks counterparts with notation like "a velocity [Z T-1 ~> m s-1]".  If the units
! vary with the Boussinesq approximation, the Boussinesq variant is given first.

character(len=40) :: mdl = "MOM_coord_initialization" !< This module's name.


  interface
module subroutine MOM_initialize_coord(GV, US, PF, tv, max_depth)
  type(verticalGrid_type), intent(inout) :: GV         !< Ocean vertical grid structure.
  type(unit_scale_type),   intent(in)    :: US         !< A dimensional unit scaling type
  type(param_file_type),   intent(in)    :: PF         !< A structure indicating the open file
                                                       !! to parse for model parameter values.
  type(thermo_var_ptrs),   intent(inout) :: tv         !< The thermodynamic variable structure.
  real,                    intent(in)    :: max_depth  !< The ocean's maximum depth [Z ~> m].
  ! Local
  ! This include declares and sets the variable "version".

end subroutine MOM_initialize_coord
module subroutine set_coord_from_gprime(Rlay, g_prime, GV, US, param_file)
  type(verticalGrid_type),  intent(in)  :: GV         !< The ocean's vertical grid structure.
  real, dimension(GV%ke),   intent(out) :: Rlay       !< The layers' target coordinate values
                                                      !! (potential density) [R ~> kg m-3].
  real, dimension(GV%ke+1), intent(out) :: g_prime    !< The reduced gravity across the interfaces
                                                      !! [L2 Z-1 T-2 ~> m s-2].
  type(unit_scale_type),    intent(in)  :: US         !< A dimensional unit scaling type
  type(param_file_type),    intent(in)  :: param_file !< A structure to parse for run-time parameters
  ! Local variables
end subroutine set_coord_from_gprime
module subroutine set_coord_from_layer_density(Rlay, g_prime, GV, US, param_file)
  type(verticalGrid_type),  intent(in)  :: GV         !< The ocean's vertical grid structure.
  real, dimension(GV%ke),   intent(out) :: Rlay       !< The layers' target coordinate values
                                                      !! (potential density) [R ~> kg m-3].
  real, dimension(GV%ke+1), intent(out) :: g_prime    !< The reduced gravity across the interfaces
                                                      !! [L2 Z-1 T-2 ~> m s-2].
  type(unit_scale_type),    intent(in)  :: US         !< A dimensional unit scaling type
  type(param_file_type),    intent(in)  :: param_file !< A structure to parse for run-time parameters

  ! Local variables
end subroutine set_coord_from_layer_density
module subroutine set_coord_from_TS_ref(Rlay, g_prime, GV, US, param_file, eqn_of_state, P_Ref)
  type(verticalGrid_type),  intent(in)  :: GV         !< The ocean's vertical grid structure.
  real, dimension(GV%ke),   intent(out) :: Rlay       !< The layers' target coordinate values
                                                      !! (potential density) [R ~> kg m-3].
  real, dimension(GV%ke+1), intent(out) :: g_prime    !< The reduced gravity across the interfaces
                                                      !! [L2 Z-1 T-2 ~> m s-2].
  type(unit_scale_type),    intent(in)  :: US         !< A dimensional unit scaling type
  type(param_file_type),    intent(in)  :: param_file !< A structure to parse for run-time parameters
  type(EOS_type),           intent(in)  :: eqn_of_state !< Equation of state structure
  real,                     intent(in)  :: P_Ref      !< The coordinate-density reference pressure
                                                      !! [R L2 T-2 ~> Pa].

  ! Local variables
end subroutine set_coord_from_TS_ref
module subroutine set_coord_from_TS_profile(Rlay, g_prime, GV, US, param_file, eqn_of_state, P_Ref)
  type(verticalGrid_type),  intent(in)  :: GV      !< The ocean's vertical grid structure
  real, dimension(GV%ke),   intent(out) :: Rlay    !< Layer potential density [R ~> kg m-3].
  real, dimension(GV%ke+1), intent(out) :: g_prime !< The reduced gravity at each
                                                   !! interface [L2 Z-1 T-2 ~> m s-2].
  type(unit_scale_type),    intent(in)  :: US      !< A dimensional unit scaling type
  type(param_file_type),    intent(in)  :: param_file !< A structure to parse for run-time parameters
  type(EOS_type),           intent(in)  :: eqn_of_state !< Equation of state structure
  real,                     intent(in)  :: P_Ref   !< The coordinate-density reference pressure
                                                   !! [R L2 T-2 ~> Pa].

  ! Local variables

end subroutine set_coord_from_TS_profile
module subroutine set_coord_from_TS_range(Rlay, g_prime, GV, US, param_file, eqn_of_state, P_Ref)
  type(verticalGrid_type),  intent(in)  :: GV      !< The ocean's vertical grid structure
  real, dimension(GV%ke),   intent(out) :: Rlay    !< Layer potential density [R ~> kg m-3].
  real, dimension(GV%ke+1), intent(out) :: g_prime !< The reduced gravity at each
                                                   !! interface [L2 Z-1 T-2 ~> m s-2].
  type(unit_scale_type),    intent(in)  :: US      !< A dimensional unit scaling type
  type(param_file_type),    intent(in)  :: param_file !< A structure to parse for run-time parameters
  type(EOS_type),           intent(in)  :: eqn_of_state !< Equation of state structure
  real,                     intent(in)  :: P_Ref   !< The coordinate-density reference pressure
                                                   !! [R L2 T-2 ~> Pa].

  ! Local variables
                  ! of the range to that in the lighter part of the range.
                  ! Setting this greater than 1 increases the resolution for
                  ! the denser water [nondim].

end subroutine set_coord_from_TS_range
module subroutine set_coord_from_file(Rlay, g_prime, GV, US, param_file)
  type(verticalGrid_type),  intent(in)  :: GV      !< The ocean's vertical grid structure
  real, dimension(GV%ke),   intent(out) :: Rlay    !< Layer potential density [R ~> kg m-3].
  real, dimension(GV%ke+1), intent(out) :: g_prime !< The reduced gravity at each
                                                   !! interface [L2 Z-1 T-2 ~> m s-2].
  type(unit_scale_type),    intent(in)  :: US      !< A dimensional unit scaling type
  type(param_file_type),    intent(in)  :: param_file !< A structure to parse for run-time parameters

  ! Local variables
end subroutine set_coord_from_file
module subroutine set_coord_linear(Rlay, g_prime, GV, US, param_file)
  type(verticalGrid_type),  intent(in)  :: GV      !< The ocean's vertical grid structure
  real, dimension(GV%ke),   intent(out) :: Rlay    !< Layer potential density [R ~> kg m-3].
  real, dimension(GV%ke+1), intent(out) :: g_prime !< The reduced gravity at each
                                                   !! interface [L2 Z-1 T-2 ~> m s-2].
  type(unit_scale_type),    intent(in)  :: US      !< A dimensional unit scaling type
  type(param_file_type),    intent(in)  :: param_file !< A structure to parse for run-time parameters

  ! Local variables
end subroutine set_coord_linear
module subroutine set_coord_to_none(Rlay, g_prime, GV, US, param_file)
  type(verticalGrid_type),  intent(in)  :: GV      !< The ocean's vertical grid structure
  real, dimension(GV%ke),   intent(out) :: Rlay    !< Layer potential density [R ~> kg m-3].
  real, dimension(GV%ke+1), intent(out) :: g_prime !< The reduced gravity at each
                                                   !! interface [L2 Z-1 T-2 ~> m s-2].
  type(unit_scale_type),    intent(in)  :: US      !< A dimensional unit scaling type
  type(param_file_type),    intent(in)  :: param_file !< A structure to parse for run-time parameters
  ! Local variables
end subroutine set_coord_to_none
module subroutine write_vertgrid_file(GV, US, param_file, directory)
  type(verticalGrid_type), intent(in)  :: GV         !< The ocean's vertical grid structure
  type(unit_scale_type),   intent(in)  :: US         !< A dimensional unit scaling type
  type(param_file_type),   intent(in)  :: param_file !< A structure to parse for run-time parameters
  character(len=*),        intent(in)  :: directory  !< The directory into which to place the file.
  ! Local variables

end subroutine write_vertgrid_file
  end interface

end module MOM_coord_initialization
