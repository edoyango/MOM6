! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> Dummy interfaces for writing ODA data
module write_ocean_obs_mod

use ocean_da_types_mod, only : ocean_profile_type
use MOM_time_manager, only : time_type, get_time, set_date

implicit none ; private

public :: open_profile_file, write_profile, close_profile_file, write_ocean_obs_init


  interface
integer module function open_profile_file(name, nvar, grid_lon, grid_lat, thread, fset)
  character(len=*), intent(in) :: name !< File name
  integer, intent(in), optional :: nvar !< Number of variables
  real, dimension(:), optional, intent(in) :: grid_lon !< Longitude [degreeE]
  real, dimension(:), optional, intent(in) :: grid_lat !< Latitude [degreeN]
  integer, optional, intent(in) :: thread !< Thread number
  integer, optional, intent(in) :: fset !< File set

end function open_profile_file
module subroutine write_profile(unit,profile)
  integer, intent(in) :: unit !< File unit
  type(ocean_profile_type), intent(in) :: profile !< Profile to write

end subroutine write_profile
module subroutine close_profile_file(unit)
  integer, intent(in) :: unit !< File unit

end subroutine close_profile_file
module subroutine write_ocean_obs_init()

end subroutine write_ocean_obs_init
  end interface

end module write_ocean_obs_mod
