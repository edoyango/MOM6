! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> This module is used to check the dimensional scaling factors used by the MOM6 ocean model
module MOM_check_scaling

use MOM_error_handler,        only : MOM_error, MOM_mesg, FATAL, WARNING, assert, MOM_get_verbosity
use MOM_unique_scales,        only : check_scaling_uniqueness, scales_to_powers
use MOM_unit_scaling,         only : unit_scale_type
use MOM_verticalGrid,         only : verticalGrid_type

implicit none ; private

! A note on unit descriptions in comments: MOM6 uses units that can be rescaled for dimensional
! consistency testing. These are noted in comments with units like Z, H, L, and T, along with
! their mks counterparts with notation like "a velocity [Z T-1 ~> m s-1]".  If the units
! vary with the Boussinesq approximation, the Boussinesq variant is given first.

public check_MOM6_scaling_factors


  interface
module subroutine check_MOM6_scaling_factors(GV, US)
  type(verticalGrid_type), pointer    :: GV         !< The container for vertical grid data
  type(unit_scale_type),   intent(in) :: US         !< A dimensional unit scaling type

  ! Local variables

  ! If no scaling is being done, simply return.
end subroutine check_MOM6_scaling_factors
module subroutine compose_dimension_list(ns, des, wts)
  integer,                       intent(out)   :: ns     !< The running sum of valid descriptions
  character(len=*), allocatable, intent(inout) :: des(:) !< The unit descriptions that have been converted
  integer,          allocatable, intent(inout) :: wts(:) !< A list of the integer weights for each scaling factor,
                                                         !! perhaps the number of times it occurs in the MOM6 code.

end subroutine compose_dimension_list
module subroutine add_scaling(ns, descs, wts, scaling, weight)
  integer,                       intent(inout) :: ns       !< The running sum of valid descriptions.
  character(len=*), allocatable, intent(inout) :: descs(:) !< The unit descriptions that have been converted
  integer,          allocatable, intent(inout) :: wts(:)   !< A list of the integers for each scaling
  character(len=*),              intent(in)    :: scaling  !< The unit description that will be converted
  integer,             optional, intent(in)    :: weight   !< An optional weight or occurrence count
                                                           !! for this unit description, 1 by default.


end subroutine add_scaling
  end interface

end module MOM_check_scaling
