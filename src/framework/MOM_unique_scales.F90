! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> This module provides tools that can be used to check the uniqueness of the dimensional
!! scaling factors used by the MOM6 ocean model or other models
module MOM_unique_scales

use MOM_error_handler, only : MOM_error, MOM_mesg, FATAL, WARNING, assert, MOM_get_verbosity

implicit none ; private

! A note on unit descriptions in comments: MOM6 uses units that can be rescaled for dimensional
! consistency testing. These are noted in comments with units like Z, H, L, and T, along with
! their mks counterparts with notation like "a velocity [Z T-1 ~> m s-1]".  If the units
! vary with the Boussinesq approximation, the Boussinesq variant is given first.

public check_scaling_uniqueness, scales_to_powers


  interface
module subroutine check_scaling_uniqueness(component, descs, weights, key, scales, max_powers)
  character(len=*),  intent(in) :: component  !< The name of the component (e.g., MOM6) to use in messages
  character(len=*),  intent(in) :: descs(:)   !< The descriptions for each combination of units
  integer,           intent(in) :: weights(:) !< A list of the weights for each described combination
  character(len=*),  intent(in) :: key(:)     !< The key for the unit scaling
  integer,           intent(in) :: scales(:)  !< The powers of 2 that give the scaling for each unit in key
  integer, optional, intent(in) :: max_powers !< The maximum range of powers of 2 to search for
                                              !! suggestions of better scaling factors, or 0 to avoid
                                              !! suggesting improved factors.

  ! Local variables

end subroutine check_scaling_uniqueness
module subroutine encode_dim_powers(scaling, key, dim_powers)

  character(len=*),               intent(in)  :: scaling   !< The unit description that will be converted
  character(len=*), dimension(:), intent(in)  :: key(:)    !< The key for the unit scaling
  integer, dimension(size(key)),  intent(out) :: dim_powers !< The dimensions in scaling of each
                                                           !! element of they key.

  ! Local variables
  ! character(len=128) :: mesg

end subroutine encode_dim_powers
module subroutine scales_to_powers(scale, pow2)
  real,    intent(in)  :: scale(:)  !< The scaling factor for each dimension
  integer, intent(out) :: pow2(:)   !< The exact powers of 2 for each scale, or 0 for non-exact powers of 2.


end subroutine scales_to_powers
integer module function non_unique_scales(scales, list, descs, weights, silent)
  integer,           intent(in) :: scales(:)  !< The power of 2 that gives the scaling factor for each dimension
  integer,           intent(in) :: list(:,:)  !< A list of the integers for each scaling
  character(len=*),  intent(in) :: descs(:)   !< The unit descriptions that have been converted
  integer,           intent(in) :: weights(:) !< A list of the weights for each scaling
  logical, optional, intent(in) :: silent     !< If present and true, do not write any output.

  ! Local variables
                                                ! for the dimensions being tested.

end function non_unique_scales
module function int_array_msg(array)
  integer,  intent(in) :: array(:)  !< The array whose values are to be written.
  character(len=16*size(array)) :: int_array_msg

end function int_array_msg
  end interface

end module MOM_unique_scales
