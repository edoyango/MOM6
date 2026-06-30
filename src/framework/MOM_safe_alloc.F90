! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> Convenience functions for safely allocating memory without
!! accidentally reallocating pointer and causing memory leaks.
module MOM_safe_alloc

implicit none ; private

public safe_alloc_ptr, safe_alloc_alloc

!> Allocate a pointer to a 1-d, 2-d or 3-d array
interface safe_alloc_ptr
  module procedure safe_alloc_ptr_3d_3arg,  safe_alloc_ptr_3d_6arg, safe_alloc_ptr_2d_2arg
  module procedure safe_alloc_ptr_3d, safe_alloc_ptr_2d, safe_alloc_ptr_1d
end interface safe_alloc_ptr

!> Allocate a 2-d or 3-d allocatable array
interface safe_alloc_alloc
  module procedure safe_alloc_allocatable_3d, safe_alloc_allocatable_2d
  module procedure safe_alloc_allocatable_3d_6arg
end interface safe_alloc_alloc

!   This combined interface might work with a later version of Fortran, but
! it fails with the gnu F90 compiler.
!
! interface safe_alloc
!   module procedure safe_alloc_ptr_3d_2arg, safe_alloc_ptr_2d_2arg
!   module procedure safe_alloc_ptr_3d, safe_alloc_ptr_2d, safe_alloc_ptr_1d
!   module procedure safe_alloc_allocatable_3d, safe_alloc_allocatable_2d
! end interface safe_alloc


  interface
module subroutine safe_alloc_ptr_1d(ptr, i1, i2)
  real, dimension(:), pointer :: ptr !< A pointer to allocate
  integer,            intent(in) :: i1 !< The size of the array, or its starting index if i2 is present
  integer, optional,  intent(in) :: i2 !< The ending index of the array
end subroutine safe_alloc_ptr_1d
module subroutine safe_alloc_ptr_2d_2arg(ptr, ni, nj)
  real, dimension(:,:), pointer :: ptr !< A pointer to allocate
  integer, intent(in) :: ni !< The size of the 1st dimension of the array
  integer, intent(in) :: nj !< The size of the 2nd dimension of the array
end subroutine safe_alloc_ptr_2d_2arg
module subroutine safe_alloc_ptr_3d_3arg(ptr, ni, nj, nk)
  real, dimension(:,:,:), pointer :: ptr !< A pointer to allocate
  integer, intent(in) :: ni !< The size of the 1st dimension of the array
  integer, intent(in) :: nj !< The size of the 2nd dimension of the array
  integer, intent(in) :: nk !< The size of the 3rd dimension of the array
end subroutine safe_alloc_ptr_3d_3arg
module subroutine safe_alloc_ptr_2d(ptr, is, ie, js, je)
  real, dimension(:,:), pointer :: ptr !< A pointer to allocate
  integer, intent(in) :: is !< The start index to allocate for the 1st dimension
  integer, intent(in) :: ie !< The end index to allocate for the 1st dimension
  integer, intent(in) :: js !< The start index to allocate for the 2nd dimension
  integer, intent(in) :: je !< The end index to allocate for the 2nd dimension
end subroutine safe_alloc_ptr_2d
module subroutine safe_alloc_ptr_3d(ptr, is, ie, js, je, nk)
  real, dimension(:,:,:), pointer :: ptr !< A pointer to allocate
  integer, intent(in) :: is !< The start index to allocate for the 1st dimension
  integer, intent(in) :: ie !< The end index to allocate for the 1st dimension
  integer, intent(in) :: js !< The start index to allocate for the 2nd dimension
  integer, intent(in) :: je !< The end index to allocate for the 2nd dimension
  integer, intent(in) :: nk !< The size to allocate for the 3rd dimension
end subroutine safe_alloc_ptr_3d
module subroutine safe_alloc_ptr_3d_6arg(ptr, is, ie, js, je, ks, ke)
  real, dimension(:,:,:), pointer :: ptr !< A pointer to allocate
  integer, intent(in) :: is !< The start index to allocate for the 1st dimension
  integer, intent(in) :: ie !< The end index to allocate for the 1st dimension
  integer, intent(in) :: js !< The start index to allocate for the 2nd dimension
  integer, intent(in) :: je !< The end index to allocate for the 2nd dimension
  integer, intent(in) :: ks !< The start index to allocate for the 3rd dimension
  integer, intent(in) :: ke !< The end index to allocate for the 3rd dimension
end subroutine safe_alloc_ptr_3d_6arg
module subroutine safe_alloc_allocatable_2d(ptr, is, ie, js, je)
  real, dimension(:,:), allocatable :: ptr !< An allocatable array to allocate
  integer, intent(in) :: is !< The start index to allocate for the 1st dimension
  integer, intent(in) :: ie !< The end index to allocate for the 1st dimension
  integer, intent(in) :: js !< The start index to allocate for the 2nd dimension
  integer, intent(in) :: je !< The end index to allocate for the 2nd dimension
end subroutine safe_alloc_allocatable_2d
module subroutine safe_alloc_allocatable_3d(ptr, is, ie, js, je, nk)
  real, dimension(:,:,:), allocatable :: ptr !< An allocatable array to allocate
  integer, intent(in) :: is !< The start index to allocate for the 1st dimension
  integer, intent(in) :: ie !< The end index to allocate for the 1st dimension
  integer, intent(in) :: js !< The start index to allocate for the 2nd dimension
  integer, intent(in) :: je !< The end index to allocate for the 2nd dimension
  integer, intent(in) :: nk !< The size to allocate for the 3rd dimension
end subroutine safe_alloc_allocatable_3d
module subroutine safe_alloc_allocatable_3d_6arg(ptr, is, ie, js, je, ks, ke)
  real, dimension(:,:,:), allocatable :: ptr !< An allocatable array to allocate
  integer, intent(in) :: is !< The start index to allocate for the 1st dimension
  integer, intent(in) :: ie !< The end index to allocate for the 1st dimension
  integer, intent(in) :: js !< The start index to allocate for the 2nd dimension
  integer, intent(in) :: je !< The end index to allocate for the 2nd dimension
  integer, intent(in) :: ks !< The start index to allocate for the 3rd dimension
  integer, intent(in) :: ke !< The end index to allocate for the 3rd dimension
end subroutine safe_alloc_allocatable_3d_6arg
  end interface

end module MOM_safe_alloc
