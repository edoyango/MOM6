! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> This module specifies the initial values and evolving properties of the
!! MOM6 ice shelf, using user-provided code.
module user_shelf_init

use MOM_error_handler, only : MOM_mesg, MOM_error, FATAL, is_root_pe
use MOM_file_parser,   only : get_param, log_version, param_file_type
use MOM_grid,          only : ocean_grid_type
use MOM_time_manager,  only : time_type, set_time, time_type_to_real
use MOM_unit_scaling,  only : unit_scale_type

implicit none ; private

#include <MOM_memory.h>

public USER_initialize_shelf_mass, USER_update_shelf_mass
public USER_init_ice_thickness

! A note on unit descriptions in comments: MOM6 uses units that can be rescaled for dimensional
! consistency testing. These are noted in comments with units like Z, H, L, and T, along with
! their mks counterparts with notation like "a velocity [Z T-1 ~> m s-1]".  If the units
! vary with the Boussinesq approximation, the Boussinesq variant is given first.

!> The control structure for the user_ice_shelf module
type, public :: user_ice_shelf_CS ; private
  real :: Rho_ocean  !< The ocean's typical density [R ~> kg m-3].
  real :: max_draft  !< The maximum ocean draft of the ice shelf [Z ~> m].
  real :: min_draft  !< The minimum ocean draft of the ice shelf [Z ~> m].
  real :: flat_shelf_width !< The range over which the shelf is min_draft thick [km].
  real :: shelf_slope_scale !< The range over which the shelf slopes [km].
  real :: pos_shelf_edge_0 !< The x-position of the shelf edge at time 0 [km].
  real :: shelf_speed  !< The ice shelf speed of translation [km day-1]
  logical :: first_call = .true. !< If true, this module has not been called before.
end type user_ice_shelf_CS


  interface
module subroutine USER_initialize_shelf_mass(mass_shelf, area_shelf_h, h_shelf, hmask, G, US, CS, param_file, new_sim)

  type(ocean_grid_type),   intent(in)  :: G    !< The ocean's grid structure
  real, dimension(SZDI_(G),SZDJ_(G)), &
                           intent(out) :: mass_shelf !< The ice shelf mass per unit area averaged
                                                  !! over the full ocean cell [R Z ~> kg m-2].
  real, dimension(SZDI_(G),SZDJ_(G)), &
                           intent(out) :: h_shelf !< The ice shelf thickness [Z ~> m].
  real, dimension(SZDI_(G),SZDJ_(G)), &
                           intent(out) :: area_shelf_h !< The area per cell covered by the ice shelf [L2 ~> m2].
  real, dimension(SZDI_(G),SZDJ_(G)), &
                           intent(out) :: hmask !< A mask indicating which tracer points are
                                                !! partly or fully covered by an ice-shelf
  type(unit_scale_type),   intent(in)  :: US    !< A structure containing unit conversion factors
  type(user_ice_shelf_CS), pointer     :: CS    !< A pointer to the user ice shelf control structure
  type(param_file_type),   intent(in)  :: param_file !< A structure to parse for run-time parameters
  logical,                 intent(in)  :: new_sim  !< If true, this is a new run; otherwise it is
                                                   !! being started from a restart file.

! This subroutine sets up the initial mass and area covered by the ice shelf.

  ! call MOM_error(FATAL, "USER_shelf_init.F90, USER_set_shelf_mass: " // &
  !  "Unmodified user routine called - you must edit the routine to use it")

end subroutine USER_initialize_shelf_mass
module subroutine USER_init_ice_thickness(h_shelf, area_shelf_h, hmask, G, US, param_file)
  type(ocean_grid_type),   intent(in)  :: G    !< The ocean's grid structure
  real, dimension(SZDI_(G),SZDJ_(G)), &
                           intent(out) :: h_shelf !< The ice shelf thickness [Z ~> m].
  real, dimension(SZDI_(G),SZDJ_(G)), &
                           intent(out) :: area_shelf_h !< The area per cell covered by the ice shelf [L2 ~> m2].
  real, dimension(SZDI_(G),SZDJ_(G)), &
                           intent(out) :: hmask !< A mask indicating which tracer points are
                                                !! partly or fully covered by an ice-shelf [nondim]
  type(unit_scale_type),   intent(in)  :: US    !< A structure containing unit conversion factors
  type(param_file_type),   intent(in)  :: param_file !< A structure to parse for run-time parameters

  ! This subroutine initializes the ice shelf thickness.  Currently it does so
  ! calling USER_initialize_shelf_mass, but this can be revised as needed.
                                                 ! over the full ocean cell [R Z ~> kg m-2].

end subroutine USER_init_ice_thickness
module subroutine USER_update_shelf_mass(mass_shelf, area_shelf_h, h_shelf, hmask, G, CS, Time, new_sim)
  type(ocean_grid_type),   intent(in)    :: G    !< The ocean's grid structure
  real, dimension(SZDI_(G),SZDJ_(G)), &
                           intent(inout) :: mass_shelf !< The ice shelf mass per unit area averaged
                                                  !! over the full ocean cell [R Z ~> kg m-2].
  real, dimension(SZDI_(G),SZDJ_(G)), &
                           intent(inout) :: area_shelf_h !< The area per cell covered by the ice shelf [L2 ~> m2].
  real, dimension(SZDI_(G),SZDJ_(G)), &
                           intent(inout) :: h_shelf !< The ice shelf thickness [Z ~> m].
  real, dimension(SZDI_(G),SZDJ_(G)), &
                           intent(inout) :: hmask !< A mask indicating which tracer points are
                                                  !! partly or fully covered by an ice-shelf [nondim]
  type(user_ice_shelf_CS), pointer       :: CS   !< A pointer to the user ice shelf control structure
  type(time_type),         intent(in)    :: Time !< The current model time
  logical,                 intent(in)    :: new_sim !< If true, this the start of a new run.



end subroutine USER_update_shelf_mass
module subroutine write_user_log(param_file)
  type(param_file_type), intent(in) :: param_file !< A structure to parse for run-time parameters


end subroutine write_user_log
  end interface

end module user_shelf_init
