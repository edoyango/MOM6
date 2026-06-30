! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> Module for supporting the rotation of a field's index map.
!! The implementation of each angle is described below.
!!
!! +90deg: B(i,j) = A(n-j,i)
!!                = transpose, then row reverse
!! 180deg: B(i,j) = A(m-i,n-j)
!!                = row reversal + column reversal
!! -90deg: B(i,j) = A(j,m-i)
!!                = row reverse, then transpose
!!
!! 90 degree rotations change the shape of the field, and are handled
!! separately from 180 degree rotations.
!!
!! It also provides the symmetric_sum functions to do a rotationally invariant
!! sum of the contents of a 1d or 2d array.

module MOM_array_transform

use iso_fortran_env, only : stdout=>output_unit, stderr=>error_unit

implicit none ; private

public rotate_array
public rotate_array_pair
public rotate_vector
public allocate_rotated_array
public symmetric_sum
public symmetric_sum_unit_tests


!> Rotate the elements of an array to the rotated set of indices.
!! Rotation is applied across the first and second axes of the array.
interface rotate_array
  module procedure rotate_array_real_2d
  module procedure rotate_array_real_3d
  module procedure rotate_array_real_4d
  module procedure rotate_array_integer
  module procedure rotate_array_logical
end interface rotate_array


!> Rotate a pair of arrays which map to a rotated set of indices.
!! Rotation is applied across the first and second axes of the array.
!! This rotation should be applied when one field is mapped onto the other.
!! For example, a tracer indexed along u or v face points will map from one
!! to the other after a quarter turn, and back onto itself after a half turn.
interface rotate_array_pair
  module procedure rotate_array_pair_real_2d
  module procedure rotate_array_pair_real_3d
  module procedure rotate_array_pair_integer
end interface rotate_array_pair


!> Rotate an array pair representing the components of a vector.
!! Rotation is applied across the first and second axes of the array.
!! This rotation should be applied when the fields satisfy vector
!! transformation rules.  For example, the u and v components of a velocity
!! will map from one to the other for quarter turns, with a sign change in one
!! component.  A half turn will map elements onto themselves with sign changes
!! in both components.
interface rotate_vector
  module procedure rotate_vector_real_2d
  module procedure rotate_vector_real_3d
  module procedure rotate_vector_real_4d
end interface rotate_vector


!> Allocate an array based on the rotated index map of an unrotated reference array.
interface allocate_rotated_array
  module procedure allocate_rotated_array_real_2d
  module procedure allocate_rotated_array_real_3d
  module procedure allocate_rotated_array_real_4d
  module procedure allocate_rotated_array_integer
end interface allocate_rotated_array


!> Return a rotationally symmetric sum of the elements of an array.
interface symmetric_sum
  module procedure symmetric_sum_1d, symmetric_sum_2d
end interface symmetric_sum



  interface
module subroutine rotate_array_real_2d(A_in, turns, A)
  real, intent(in) :: A_in(:,:) !< Unrotated array [arbitrary]
  integer, intent(in) :: turns  !< Number of quarter turns
  real, intent(out) :: A(:,:)   !< Rotated array [arbitrary]


end subroutine rotate_array_real_2d
module subroutine rotate_array_real_3d(A_in, turns, A)
  real, intent(in) :: A_in(:,:,:) !< Unrotated array [arbitrary]
  integer, intent(in) :: turns    !< Number of quarter turns
  real, intent(out) :: A(:,:,:)   !< Rotated array [arbitrary]


end subroutine rotate_array_real_3d
module subroutine rotate_array_real_4d(A_in, turns, A)
  real, intent(in) :: A_in(:,:,:,:) !< Unrotated array [arbitrary]
  integer, intent(in) :: turns      !< Number of quarter turns
  real, intent(out) :: A(:,:,:,:)   !< Rotated array [arbitrary]


end subroutine rotate_array_real_4d
module subroutine rotate_array_integer(A_in, turns, A)
  integer, intent(in) :: A_in(:,:)  !< Unrotated array
  integer, intent(in) :: turns      !< Number of quarter turns
  integer, intent(out) :: A(:,:)    !< Rotated array


end subroutine rotate_array_integer
module subroutine rotate_array_logical(A_in, turns, A)
  logical, intent(in) :: A_in(:,:)  !< Unrotated array
  integer, intent(in) :: turns      !< Number of quarter turns
  logical, intent(out) :: A(:,:)    !< Rotated array


end subroutine rotate_array_logical
module subroutine rotate_array_pair_real_2d(A_in, B_in, turns, A, B)
  real, intent(in) :: A_in(:,:)   !< Unrotated scalar array pair [arbitrary]
  real, intent(in) :: B_in(:,:)   !< Unrotated scalar array pair [arbitrary]
  integer, intent(in) :: turns    !< Number of quarter turns
  real, intent(out) :: A(:,:)     !< Rotated scalar array pair [arbitrary]
  real, intent(out) :: B(:,:)     !< Rotated scalar array pair [arbitrary]

end subroutine rotate_array_pair_real_2d
module subroutine rotate_array_pair_real_3d(A_in, B_in, turns, A, B)
  real, intent(in) :: A_in(:,:,:)   !< Unrotated scalar array pair [arbitrary]
  real, intent(in) :: B_in(:,:,:)   !< Unrotated scalar array pair [arbitrary]
  integer, intent(in) :: turns      !< Number of quarter turns
  real, intent(out) :: A(:,:,:)     !< Rotated scalar array pair [arbitrary]
  real, intent(out) :: B(:,:,:)     !< Rotated scalar array pair [arbitrary]


end subroutine rotate_array_pair_real_3d
module subroutine rotate_array_pair_integer(A_in, B_in, turns, A, B)
  integer, intent(in) :: A_in(:,:)  !< Unrotated scalar array pair
  integer, intent(in) :: B_in(:,:)  !< Unrotated scalar array pair
  integer, intent(in) :: turns      !< Number of quarter turns
  integer, intent(out) :: A(:,:)    !< Rotated scalar array pair
  integer, intent(out) :: B(:,:)    !< Rotated scalar array pair

end subroutine rotate_array_pair_integer
module subroutine rotate_vector_real_2d(A_in, B_in, turns, A, B)
  real, intent(in) :: A_in(:,:) !< First component of unrotated vector [arbitrary]
  real, intent(in) :: B_in(:,:) !< Second component of unrotated vector [arbitrary]
  integer, intent(in) :: turns  !< Number of quarter turns
  real, intent(out) :: A(:,:)   !< First component of rotated vector [arbitrary]
  real, intent(out) :: B(:,:)   !< Second component of unrotated vector [arbitrary]

end subroutine rotate_vector_real_2d
module subroutine rotate_vector_real_3d(A_in, B_in, turns, A, B)
  real, intent(in) :: A_in(:,:,:) !< First component of unrotated vector [arbitrary]
  real, intent(in) :: B_in(:,:,:) !< Second component of unrotated vector [arbitrary]
  integer, intent(in) :: turns    !< Number of quarter turns
  real, intent(out) :: A(:,:,:)   !< First component of rotated vector [arbitrary]
  real, intent(out) :: B(:,:,:)   !< Second component of unrotated vector [arbitrary]


end subroutine rotate_vector_real_3d
module subroutine rotate_vector_real_4d(A_in, B_in, turns, A, B)
  real, intent(in) :: A_in(:,:,:,:) !< First component of unrotated vector [arbitrary]
  real, intent(in) :: B_in(:,:,:,:) !< Second component of unrotated vector [arbitrary]
  integer, intent(in) :: turns      !< Number of quarter turns
  real, intent(out) :: A(:,:,:,:)   !< First component of rotated vector [arbitrary]
  real, intent(out) :: B(:,:,:,:)   !< Second component of unrotated vector [arbitrary]


end subroutine rotate_vector_real_4d
module subroutine allocate_rotated_array_real_2d(A_in, lb, turns, A)
  ! NOTE: lb must be declared before A_in
  integer, intent(in) :: lb(2)                !< Lower index bounds of A_in
  real, intent(in) :: A_in(lb(1):, lb(2):)    !< Reference array [arbitrary]
  integer, intent(in) :: turns                !< Number of quarter turns
  real, allocatable, intent(inout) :: A(:,:)  !< Array on rotated index [arbitrary]


end subroutine allocate_rotated_array_real_2d
module subroutine allocate_rotated_array_real_3d(A_in, lb, turns, A)
  ! NOTE: lb must be declared before A_in
  integer, intent(in) :: lb(3)                    !< Lower index bounds of A_in
  real, intent(in) :: A_in(lb(1):, lb(2):, lb(3):)  !< Reference array [arbitrary]
  integer, intent(in) :: turns                    !< Number of quarter turns
  real, allocatable, intent(inout) :: A(:,:,:)    !< Array on rotated index [arbitrary]


end subroutine allocate_rotated_array_real_3d
module subroutine allocate_rotated_array_real_4d(A_in, lb, turns, A)
  ! NOTE: lb must be declared before A_in
  integer, intent(in) :: lb(4)                    !< Lower index bounds of A_in
  real, intent(in) :: A_in(lb(1):,lb(2):,lb(3):,lb(4):) !< Reference array [arbitrary]
  integer, intent(in) :: turns                    !< Number of quarter turns
  real, allocatable, intent(inout) :: A(:,:,:,:)  !< Array on rotated index [arbitrary]


end subroutine allocate_rotated_array_real_4d
module subroutine allocate_rotated_array_integer(A_in, lb, turns, A)
  integer, intent(in) :: lb(2)                  !< Lower index bounds of A_in
  integer, intent(in) :: A_in(lb(1):,lb(2):)    !< Reference array
  integer, intent(in) :: turns                  !< Number of quarter turns
  integer, allocatable, intent(inout) :: A(:,:) !< Array on rotated index


end subroutine allocate_rotated_array_integer
module function symmetric_sum_1d(field) result(sum)
  real, dimension(1:), intent(in) :: field !< The field to sum in arbitrary units [A ~> a]
  real :: sum !< The rotationally symmetric sum of the entries in field [A ~> a]

  ! Local variables

end function symmetric_sum_1d
recursive module function symmetric_sum_2d(field) result(sum)
  real, dimension(1:,1:), intent(in) :: field !< The field to sum in arbitrary units [A ~> a]
  real :: sum !< The rotationally symmetric sum of the entries in field [A ~> a]

  ! Local variables

end function symmetric_sum_2d
module function naive_sum_2d(field, abs_val) result(sum)
  real, dimension(1:,1:), intent(in) :: field !< The field to sum in arbitrary units [A ~> a]
  logical, optional,      intent(in) :: abs_val !< If present and true, sum the absolute values
  real :: sum !< The rotation dependent sum of the entries in field [A ~> a]

  ! Local variables

end function naive_sum_2d
logical module function symmetric_sum_unit_tests(verbose)
  ! Arguments
  logical, intent(in) :: verbose !< If true, write results to stdout
  ! Local variables


end function symmetric_sum_unit_tests
  end interface

end module MOM_array_transform
