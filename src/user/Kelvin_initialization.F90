! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> Configures the model for the Kelvin wave experiment.
!!
!! Kelvin = coastally-trapped Kelvin waves from the ROMS examples.
!! Initialize with level surfaces and drive the wave in at the west,
!! radiate out at the east.
module Kelvin_initialization

use MOM_dyn_horgrid,    only : dyn_horgrid_type
use MOM_error_handler,  only : MOM_mesg, MOM_error, FATAL, WARNING, is_root_pe
use MOM_file_parser,    only : get_param, log_version, param_file_type
use MOM_grid,           only : ocean_grid_type
use MOM_open_boundary,  only : ocean_OBC_type, OBC_NONE
use MOM_open_boundary,  only : OBC_segment_type, register_OBC, rotate_OBC_segment_direction
use MOM_open_boundary,  only : OBC_DIRECTION_N, OBC_DIRECTION_E, OBC_DIRECTION_S, OBC_DIRECTION_W
use MOM_open_boundary,  only : OBC_registry_type
use MOM_unit_scaling,   only : unit_scale_type
use MOM_verticalGrid,   only : verticalGrid_type
use MOM_time_manager,   only : time_type, time_to_real

implicit none ; private

#include <MOM_memory.h>

public Kelvin_set_OBC_data, Kelvin_initialize_topography
public register_Kelvin_OBC, Kelvin_OBC_end

! A note on unit descriptions in comments: MOM6 uses units that can be rescaled for dimensional
! consistency testing. These are noted in comments with units like Z, H, L, and T, along with
! their mks counterparts with notation like "a velocity [Z T-1 ~> m s-1]".  If the units
! vary with the Boussinesq approximation, the Boussinesq variant is given first.

!> Control structure for Kelvin wave open boundaries.
type, public :: Kelvin_OBC_CS ; private
  integer :: mode = 0          !< Vertical mode
  real :: coast_angle = 0   !< Angle of coastline [rad]
  real :: coast_offset1 = 0 !< Longshore distance to coastal angle [L ~> m]
  real :: coast_offset2 = 0 !< Offshore distance to coastal angle [L ~> m]
  real :: H0 = 0            !< Bottom depth [Z ~> m]
  real :: F_0               !< Coriolis parameter [T-1 ~> s-1]
  real :: rho_range         !< Density range [R ~> kg m-3]
  real :: rho_0             !< Mean density [R ~> kg m-3]
  real :: wave_period       !< Period of the mode-0 waves [T ~> s]
  real :: ssh_amp           !< Amplitude of the sea surface height forcing for mode-0 waves [Z ~> m]
  real :: inflow_amp        !< Amplitude of the boundary velocity forcing for internal waves [L T-1 ~> m s-1]
  real :: OBC_nudging_time  !< The timescale with which the inflowing open boundary velocities are nudged toward
                            !! their intended values with the Kelvin wave test case [T ~> s], or a negative
                            !! value to retain the value that is set when the OBC segments are initialized.
  logical :: indexing_bugs  !< If true, retain several horizontal indexing bugs that were in the
                            !! original version of Kelvin_set_OBC_data.
end type Kelvin_OBC_CS

! This include declares and sets the variable "version".
#include "version_variable.h"


  interface
logical module function register_Kelvin_OBC(param_file, CS, US, OBC_Reg)
  type(param_file_type),    intent(in) :: param_file !< parameter file.
  type(Kelvin_OBC_CS),      pointer    :: CS         !< Kelvin wave control structure.
  type(unit_scale_type),    intent(in) :: US         !< A dimensional unit scaling type
  type(OBC_registry_type),  pointer    :: OBC_Reg    !< OBC registry.

  ! Local variables
                          ! recreate the bugs, or if false bugs are only used if actively selected.

end function register_Kelvin_OBC
module subroutine Kelvin_OBC_end(CS)
  type(Kelvin_OBC_CS), pointer    :: CS         !< Kelvin wave control structure.

end subroutine Kelvin_OBC_end
module subroutine Kelvin_initialize_topography(D, G, param_file, max_depth, US)
  type(dyn_horgrid_type),          intent(in)  :: G !< The dynamic horizontal grid type
  real, dimension(G%isd:G%ied,G%jsd:G%jed), &
                                   intent(out) :: D !< Ocean bottom depth [Z ~> m]
  type(param_file_type),           intent(in)  :: param_file !< Parameter file structure
  real,                            intent(in)  :: max_depth !< Maximum model depth [Z ~> m]
  type(unit_scale_type),           intent(in)  :: US !< A dimensional unit scaling type

  ! Local variables

end subroutine Kelvin_initialize_topography
module subroutine Kelvin_set_OBC_data(OBC, CS, G, GV, US, h, Time)
  type(ocean_OBC_type),    pointer    :: OBC  !< This open boundary condition type specifies
                                              !! whether, where, and what open boundary
                                              !! conditions are used.
  type(Kelvin_OBC_CS),     pointer    :: CS   !< Kelvin wave control structure.
  type(ocean_grid_type),   intent(in) :: G    !< The ocean's grid structure.
  type(verticalGrid_type), intent(in) :: GV   !< The ocean's vertical grid structure.
  type(unit_scale_type),   intent(in) :: US   !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in) :: h !< layer thickness [H ~> m or kg m-2].
  type(time_type),         intent(in) :: Time !< model time.

  ! The following variables are used to set up the transport in the Kelvin example.

end subroutine Kelvin_set_OBC_data
  end interface

end module Kelvin_initialization
