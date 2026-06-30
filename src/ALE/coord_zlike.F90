! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> Regrid columns for a z-like coordinate (z-star, z-level)
module coord_zlike

use MOM_error_handler, only : MOM_error, FATAL

implicit none ; private

!> Control structure containing required parameters for a z-like coordinate
type, public :: zlike_CS ; private

  !> Number of levels to be generated
  integer :: nk

  !> Minimum thickness allowed for layers, in the same thickness units (perhaps [H ~> m or kg m-2])
  !! that will be used in all subsequent calls to build_zstar_column with this structure.
  real :: min_thickness

  !> Target coordinate resolution, usually in [Z ~> m]
  real, allocatable, dimension(:) :: coordinateResolution
end type zlike_CS

public init_coord_zlike, set_zlike_params, build_zstar_column, end_coord_zlike


  interface
module subroutine init_coord_zlike(CS, nk, coordinateResolution)
  type(zlike_CS),     pointer    :: CS !< Unassociated pointer to hold the control structure
  integer,            intent(in) :: nk !< Number of levels in the grid
  real, dimension(:), intent(in) :: coordinateResolution !< Target coordinate resolution [Z ~> m]

end subroutine init_coord_zlike
module subroutine end_coord_zlike(CS)
  type(zlike_CS), pointer :: CS !< Coordinate control structure

  ! Nothing to do
end subroutine end_coord_zlike
module subroutine set_zlike_params(CS, min_thickness)
  type(zlike_CS), pointer    :: CS !< Coordinate control structure
  real, optional, intent(in) :: min_thickness !< Minimum allowed thickness [H ~> m or kg m-2]

end subroutine set_zlike_params
module subroutine build_zstar_column(CS, depth, total_thickness, zInterface, &
                              z_rigid_top, eta_orig, zScale)
  type(zlike_CS),           intent(in)    :: CS !< Coordinate control structure
  real,                     intent(in)    :: depth !< Depth of ocean bottom (positive downward in the
                                                   !! output units), units may be [Z ~> m] or [H ~> m or kg m-2]
  real,                     intent(in)    :: total_thickness !< Column thickness (positive definite in the same
                                                   !! units as depth) [Z ~> m] or [H ~> m or kg m-2]
  real, dimension(CS%nk+1), intent(inout) :: zInterface !< Absolute positions of interfaces (in the same
                                                   !! units as depth) [Z ~> m] or [H ~> m or kg m-2]
  real, optional,           intent(in)    :: z_rigid_top !< The height of a rigid top (positive upward in the same
                                                   !! units as depth) [Z ~> m] or [H ~> m or kg m-2]
  real, optional,           intent(in)    :: eta_orig !< The actual original height of the top (in the same
                                                   !! units as depth) [Z ~> m] or [H ~> m or kg m-2]
  real, optional,           intent(in)    :: zScale !< Scaling factor from the target coordinate resolution
                                                    !! in Z to desired units for zInterface, perhaps Z_to_H,
                                                    !! often [nondim] or [H Z-1 ~> 1 or kg m-3]
  ! Local variables

end subroutine build_zstar_column
  end interface

end module coord_zlike
