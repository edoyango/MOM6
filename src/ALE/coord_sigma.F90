! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> Regrid columns for the sigma coordinate
module coord_sigma

use MOM_error_handler, only : MOM_error, FATAL

implicit none ; private

!> Control structure containing required parameters for the sigma coordinate
type, public :: sigma_CS ; private

  !> Number of levels
  integer :: nk

  !> Minimum thickness allowed for layers [H ~> m or kg m-2]
  real :: min_thickness

  !> Target coordinate resolution [nondim]
  real, allocatable, dimension(:) :: coordinateResolution
end type sigma_CS

public init_coord_sigma, set_sigma_params, build_sigma_column, end_coord_sigma


  interface
module subroutine init_coord_sigma(CS, nk, coordinateResolution)
  type(sigma_CS),     pointer    :: CS !< Unassociated pointer to hold the control structure
  integer,            intent(in) :: nk !< Number of layers in the grid
  real, dimension(:), intent(in) :: coordinateResolution !< Nominal coordinate resolution [nondim]

end subroutine init_coord_sigma
module subroutine end_coord_sigma(CS)
  type(sigma_CS), pointer :: CS !< Coordinate control structure

  ! nothing to do
end subroutine end_coord_sigma
module subroutine set_sigma_params(CS, min_thickness)
  type(sigma_CS), pointer    :: CS !< Coordinate control structure
  real, optional, intent(in) :: min_thickness !< Minimum allowed thickness [H ~> m or kg m-2]

end subroutine set_sigma_params
module subroutine build_sigma_column(CS, depth, totalThickness, zInterface)
  type(sigma_CS),           intent(in)    :: CS !< Coordinate control structure
  real,                     intent(in)    :: depth !< Depth of ocean bottom (positive [H ~> m or kg m-2])
  real,                     intent(in)    :: totalThickness !< Column thickness (positive [H ~> m or kg m-2])
  real, dimension(CS%nk+1), intent(inout) :: zInterface !< Absolute positions of interfaces [H ~> m or kg m-2]

  ! Local variables

end subroutine build_sigma_column
  end interface

end module coord_sigma
