!> Provides buffers that can dynamically grow as needed. These are primarily intended for the
!! diagnostics which need to store intermediate or partial states of state variables
module MOM_diag_buffers

use MOM_io, only : stdout, stderr

! This file is part of MOM6. See LICENSE.md for the license.

implicit none ; private

public :: diag_buffer_unit_tests_2d, diag_buffer_unit_tests_3d

type, abstract :: buffer_base
end type buffer_base

!> Holds a 2d field
type, extends(buffer_base) :: buffer_2d
  real, dimension(:,:), allocatable :: field !< The actual 2d field to be stored [arbitrary]
end type buffer_2d

!> Holds a 3d field
type, extends(buffer_base) :: buffer_3d
  real, dimension(:,:,:), allocatable :: field !< The actual 3d field to be stored [arbitrary]
end type buffer_3d

!> The base class for the diagnostic buffers in this module
type, abstract :: diag_buffer_base ; private
  integer :: is !< The start slot of the array i-direction
  integer :: js !< The start slot of the array j-direction
  integer :: ie !< The end slot of the array i-direction
  integer :: je !< The end slot of the array j-direction
  real :: fill_value = 0. !< Set the fill value to use when growing the buffer [arbitrary]
  integer, allocatable, dimension(:) :: ids  !< List of diagnostic ids whose slot corresponds to the row in the buffer
  integer :: length = 0 !< The number of slots in the buffer

  contains

  procedure(a_grow), deferred :: grow !< Increase the size of the buffer
  procedure, public :: set_fill_value !< Set the fill value to use when growing the buffer
  procedure, public :: check_capacity_by_id !< Check the size size of the buffer and increase if necessary
  procedure, public :: set_horizontal_extents !< Define the horizontal extents of the arrays
  procedure, public :: mark_available !< Mark that a slot in the buffer can be reused
  procedure, public :: grow_ids !< Increase the size of the vector storing diagnostic ids
  procedure, public :: find_buffer_slot !< Find the slot corresponding to a specific diagnostic id
end type diag_buffer_base

!> Dynamically growing buffer for 2D arrays.
type, extends(diag_buffer_base), public :: diag_buffer_2d ; private
  type(buffer_2d), public, dimension(:), allocatable :: buffer !< The actual 2D buffer which will dynamically grow

  contains

  procedure, public :: grow => grow_2d !< Increase the size of the buffer
  procedure, public :: store => store_2d !< Store a field in the buffer, increasing as necessary
  procedure, public :: set_extents_from_array => set_extents_from_array_2d !< Set extents from array bounds

  final :: finalize_diag_buffer_2d
    !< Finalization stub to improve optimized build time
end type diag_buffer_2d

!> Dynamically growing buffer for 3D arrays.
type, extends(diag_buffer_base), public :: diag_buffer_3d ; private
  type(buffer_3d), public, dimension(:), allocatable :: buffer !< The actual 2D buffer which will dynamically grow
  integer :: ks !< The start slot in the k-dimension
  integer :: ke !< The last slot in the k-dimension

  contains

  procedure, public :: set_vertical_extent !< Set the vertical extents of the buffer
  procedure, public :: grow => grow_3d !< Increase the size of the buffer
  procedure, public :: store => store_3d !< Store a field in the buffer, increasing as necessary
  procedure, public :: set_extents_from_array => set_extents_from_array_3d !< Set extents from array bounds

  final :: finalize_diag_buffer_3d
    !< Finalization stub to improve optimized build time
end type diag_buffer_3d


  interface
module subroutine a_grow(this)
  class(diag_buffer_base), intent(inout) :: this !< The diagnostic buffer
end subroutine a_grow
module subroutine set_fill_value(this, fill_value)
  class(diag_buffer_base), intent(inout) :: this !< The diagnostic buffer
  real,                    intent(in)    :: fill_value !< The fill value to use when growing the buffer [arbitrary]

end subroutine set_fill_value
module subroutine mark_available(this, id)
  class(diag_buffer_base), intent(inout) :: this !< The diagnostic buffer
  integer,                 intent(in)    :: id   !< The diagnostic id

end subroutine mark_available
pure module function find_buffer_slot(this, id) result(slot)
  class(diag_buffer_base), intent(in) :: this !< The diagnostic buffer
  integer, intent(in) :: id !< The diagnostic id

  integer :: slot !< The slot in the buffer corresponding to the diagnostic id

end function find_buffer_slot
module subroutine grow_ids(this)
  class(diag_buffer_base), intent(inout) :: this !< This buffer


end subroutine grow_ids
impure module function check_capacity_by_id(this, id) result(slot)
  class(diag_buffer_base), intent(inout) :: this !< This 2d buffer
  integer,                 intent(in)    :: id   !< The diagnostic id
  integer :: slot

end function check_capacity_by_id
module subroutine set_horizontal_extents(this, is, ie, js, je)
  class(diag_buffer_base), intent(inout) :: this !< The diagnostic buffer
  integer,               intent(in)    :: is !< The start slot of the array i-direction
  integer,               intent(in)    :: ie !< The end slot of the array i-direction
  integer,               intent(in)    :: js !< The start slot of the array j-direction
  integer,               intent(in)    :: je !< The end slot of the array j-direction

end subroutine set_horizontal_extents
module subroutine set_vertical_extent(this, ks, ke)
  class(diag_buffer_3d), intent(inout) :: this !< The diagnostic buffer
  integer,               intent(in)    :: ks !< The start slot of the array k-direction
  integer,               intent(in)    :: ke !< The end slot of the array k-direction

end subroutine set_vertical_extent
module subroutine set_extents_from_array_2d(this, array, fill_value_in)
  class(diag_buffer_2d), intent(inout) :: this !< The diagnostic buffer
  real, dimension(:,:), intent(in)     :: array !< The array whose bounds define the buffer extents
  real, optional,       intent(in)     :: fill_value_in !< Optional fill value

end subroutine set_extents_from_array_2d
module subroutine set_extents_from_array_3d(this, array, fill_value_in)
  class(diag_buffer_3d), intent(inout) :: this !< The diagnostic buffer
  real, dimension(:,:,:), intent(in)   :: array !< The array whose bounds define the buffer extents
  real, optional,         intent(in)   :: fill_value_in !< Optional fill value

end subroutine set_extents_from_array_3d
module subroutine grow_2d(this)
  class(diag_buffer_2d), intent(inout) :: this


  ! Grow the ID array
end subroutine grow_2d
module subroutine store_2d(this, data, id)
  class(diag_buffer_2d), intent(inout) :: this !< This 2d buffer
  real, dimension(:,:),  intent(in)    :: data !< The data to be stored in the buffer [arbitrary]
  integer,               intent(in)    :: id !< The diagnostic id


end subroutine store_2d
module subroutine grow_3d(this)
  class(diag_buffer_3d), intent(inout) :: this


  ! Grow the ID array
end subroutine grow_3d
module subroutine store_3d(this, data, id)
  class(diag_buffer_3d),  intent(inout) :: this !< This 3d buffer
  real, dimension(:,:,:), intent(in)    :: data !< The data to be stored in the buffer [arbitrary]
  integer,                intent(in)    :: id !< The diagnostic id


  ! Find the first slot in the ids array that is 0, i.e. this is a portion of the buffer that can be reused
end subroutine store_3d
module subroutine finalize_diag_buffer_2d(this)
  type(diag_buffer_2d) :: this
    !< Diagnostic buffer
end subroutine finalize_diag_buffer_2d
module subroutine finalize_diag_buffer_3d(this)
  type(diag_buffer_3d) :: this
    !< Diagnostic buffer
end subroutine finalize_diag_buffer_3d
module function diag_buffer_unit_tests_2d(verbose) result(fail)
  logical, intent(in) :: verbose !< If true, write results to stdout
  logical :: fail !< True if any of the unit tests fail

end function diag_buffer_unit_tests_2d
module function diag_buffer_unit_tests_3d(verbose) result(fail)
  logical, intent(in) :: verbose !< If true, write results to stdout
  logical :: fail !< True if any of the unit tests fail

end function diag_buffer_unit_tests_3d
  end interface

end module MOM_diag_buffers
