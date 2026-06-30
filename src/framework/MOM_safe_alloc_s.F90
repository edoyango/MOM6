submodule (MOM_safe_alloc) MOM_safe_alloc_s
  implicit none
contains
module procedure safe_alloc_ptr_1d
  if (.not.associated(ptr)) then
    if (present(i2)) then
      allocate(ptr(i1:i2), source=0.0)
    else
      allocate(ptr(i1), source=0.0)
    endif
  endif
end procedure safe_alloc_ptr_1d
module procedure safe_alloc_ptr_2d_2arg
  if (.not.associated(ptr)) then
    allocate(ptr(ni,nj), source=0.0)
  endif
end procedure safe_alloc_ptr_2d_2arg
module procedure safe_alloc_ptr_3d_3arg
  if (.not.associated(ptr)) then
    allocate(ptr(ni,nj,nk), source=0.0)
  endif
end procedure safe_alloc_ptr_3d_3arg
module procedure safe_alloc_ptr_2d
  if (.not.associated(ptr)) then
    allocate(ptr(is:ie,js:je), source=0.0)
  endif
end procedure safe_alloc_ptr_2d
module procedure safe_alloc_ptr_3d
  if (.not.associated(ptr)) then
    allocate(ptr(is:ie,js:je,nk), source=0.0)
  endif
end procedure safe_alloc_ptr_3d
module procedure safe_alloc_ptr_3d_6arg
  if (.not.associated(ptr)) then
    allocate(ptr(is:ie,js:je,ks:ke), source=0.0)
  endif
end procedure safe_alloc_ptr_3d_6arg
module procedure safe_alloc_allocatable_2d
  if (.not.allocated(ptr)) then
    allocate(ptr(is:ie,js:je), source=0.0)
  endif
end procedure safe_alloc_allocatable_2d
module procedure safe_alloc_allocatable_3d
  if (.not.allocated(ptr)) then
    allocate(ptr(is:ie,js:je,nk), source=0.0)
  endif
end procedure safe_alloc_allocatable_3d
module procedure safe_alloc_allocatable_3d_6arg
  if (.not.allocated(ptr)) then
    allocate(ptr(is:ie,js:je,ks:ke), source=0.0)
  endif
end procedure safe_alloc_allocatable_3d_6arg
end submodule MOM_safe_alloc_s
