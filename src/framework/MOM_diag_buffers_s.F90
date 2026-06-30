submodule (MOM_diag_buffers) MOM_diag_buffers_s
  implicit none
contains
module procedure a_grow
end procedure a_grow
module procedure set_fill_value
  this%fill_value = fill_value
end procedure set_fill_value
module procedure mark_available
  integer :: slot
  slot = this%find_buffer_slot(id)
  this%ids(slot) = 0
end procedure mark_available
module procedure find_buffer_slot
  integer, dimension(1) :: temp
  if (allocated(this%ids)) then
    !NOTE: Alternatively could do slot = SUM(findloc(...))
    temp = findloc(this%ids(:), id)
    slot = temp(1)
  else
    slot = 0
  endif

end procedure find_buffer_slot
module procedure grow_ids
  integer, allocatable, dimension(:) :: temp
  integer :: n
  n = this%length

  allocate(temp(n+1))
  if (n>0) temp(1:n) = this%ids(:)
  call move_alloc(temp, this%ids)
end procedure grow_ids
module procedure check_capacity_by_id
  slot = this%find_buffer_slot(id)
  if (slot==0) then
    ! Check to see if there is an open slot
    if (allocated(this%ids)) slot = this%find_buffer_slot(0)
    ! If slot is still 0, then the buffer must grow
    if (slot==0) then
      call this%grow()
      slot = this%length
    endif
    this%ids(slot) = id
  endif
end procedure check_capacity_by_id
module procedure set_horizontal_extents
  this%is = is ; this%ie = ie ; this%js = js ; this%je = je
end procedure set_horizontal_extents
module procedure set_vertical_extent
  this%ks = ks ; this%ke = ke
end procedure set_vertical_extent
module procedure set_extents_from_array_2d
  call this%set_horizontal_extents(lbound(array,1), ubound(array,1), &
                                   lbound(array,2), ubound(array,2))
  if (present(fill_value_in)) call this%set_fill_value(fill_value_in)
end procedure set_extents_from_array_2d
module procedure set_extents_from_array_3d
  call this%set_horizontal_extents(lbound(array,1), ubound(array,1), &
                                   lbound(array,2), ubound(array,2))
  call this%set_vertical_extent(lbound(array,3), ubound(array,3))
  if (present(fill_value_in)) call this%set_fill_value(fill_value_in)
end procedure set_extents_from_array_3d
module procedure grow_2d
  integer :: i, n
  integer :: is, ie, js, je
  type(buffer_2d), dimension(:), allocatable :: new_buffer
  call this%grow_ids()

  is = this%is ; ie = this%ie ; js = this%js ; je = this%je
  n = this%length

  allocate(new_buffer(n+1))
  do i=1,n
    allocate(new_buffer(i)%field(is:ie,js:je))
    new_buffer(i)%field(:,:) = this%buffer(i)%field(:,:)
  enddo
  allocate(new_buffer(n+1)%field(is:ie,js:je), source=this%fill_value)
  call move_alloc(new_buffer, this%buffer)
  this%length = n+1

end procedure grow_2d
module procedure store_2d
  integer :: slot
  slot = this%check_capacity_by_id(id)
  this%buffer(slot)%field(:,:) = data(:,:)
end procedure store_2d
module procedure grow_3d
  integer :: i, n
  integer :: is, ie, js, je, ks, ke
  type(buffer_3d), dimension(:), allocatable :: new_buffer
  call this%grow_ids()

  is = this%is ; ie = this%ie ; js = this%js ; je = this%je ; ks = this%ks ; ke = this%ke
  n = this%length

  allocate(new_buffer(n+1))
  do i=1,n
    allocate(new_buffer(i)%field(is:ie,js:je,ks:ke))
    new_buffer(i)%field(:,:,:) = this%buffer(i)%field(:,:,:)
  enddo
  allocate(new_buffer(n+1)%field(is:ie,js:je,ks:ke), source=this%fill_value)
  call move_alloc(new_buffer, this%buffer)
  this%length = n+1

end procedure grow_3d
module procedure store_3d
  integer :: slot
  slot = this%check_capacity_by_id(id)
  this%buffer(slot)%field(:,:,:) = data(:,:,:)
end procedure store_3d
module procedure finalize_diag_buffer_2d
end procedure finalize_diag_buffer_2d
module procedure finalize_diag_buffer_3d
end procedure finalize_diag_buffer_3d
module procedure diag_buffer_unit_tests_2d
  fail = .false.
  write(stdout,*) '==== MOM_diag_buffers: diag_buffers_unit_tests_2d ==='
  fail = fail .or. new_buffer_2d()
  fail = fail .or. grow_buffer_2d()
  fail = fail .or. fill_value_2d()
  fail = fail .or. store_buffer_2d()
  fail = fail .or. reuse_buffer_2d()

  contains

  !> Ensure properties of a newly initialized buffer
  function new_buffer_2d() result(local_fail)
    type(diag_buffer_2d) :: buffer
    logical :: local_fail !< True if any of the unit tests fail
    local_fail = .false.
    local_fail = local_fail .or. allocated(buffer%buffer)
    if (verbose) write(stdout,*) "new_buffer_2d: ", local_fail
    local_fail = local_fail .or. allocated(buffer%ids)
    if (verbose) write(stdout,*) "new_buffer_2d: ", local_fail
    local_fail = local_fail .or. buffer%length /= 0
    if (verbose) write(stdout,*) "new_buffer_2d: ", local_fail
  end function new_buffer_2d

  !> Test the growing of a buffer
  function grow_buffer_2d() result(local_fail)
    type(diag_buffer_2d) :: buffer
    logical :: local_fail !< True if any of the unit tests fail
    integer, parameter :: is=1, ie=2, js=3, je=6
    integer :: i

    local_fail = .false.

    call buffer%set_horizontal_extents(is=is, ie=ie, js=js, je=je)
    ! Grow the buffer 3 times
    do i=1,3
      call buffer%grow()
      local_fail = local_fail .or. (buffer%length /= i)
      local_fail = local_fail .or. (lbound(buffer%buffer(i)%field, 1) /= is)
      local_fail = local_fail .or. (ubound(buffer%buffer(i)%field, 1) /= ie)
      local_fail = local_fail .or. (lbound(buffer%buffer(i)%field, 2) /= js)
      local_fail = local_fail .or. (ubound(buffer%buffer(i)%field, 2) /= je)
    enddo
    if (verbose) write(stdout,*) "grow_buffer_2d: ", local_fail
  end function grow_buffer_2d

  !> Test that growing new buffer fills the array with a set fill value
  function fill_value_2d() result(local_fail)
    type(diag_buffer_2d) :: buffer
    logical :: local_fail !< True if any of the unit tests fail
    integer, parameter :: is=1, ie=2, js=3, je=6
    real, parameter :: fill_value = -123.456


    local_fail = .false.

    call buffer%set_horizontal_extents(is=is, ie=ie, js=js, je=je)
    call buffer%set_fill_value(fill_value)
    call buffer%grow()
    if (any(buffer%buffer(1)%field(:,:) /= fill_value)) local_fail = .true.
    if (verbose) write(stdout,*) "fill_value_2d: ", local_fail
  end function fill_value_2d

  !> Test storing a buffer based on a unique id
  function store_buffer_2d() result(local_fail)
    type(diag_buffer_2d) :: buffer
    logical :: local_fail !< True if any of the unit tests fail

    integer, parameter :: is=1, ie=2, js=3, je=6, nlen=3
    integer :: i, slot
    real, allocatable, dimension(:,:,:) :: test_2d

    local_fail = .false.

    allocate(test_2d(nlen, is:ie, js:je))
    call random_number(test_2d)
    buffer%is = is
    buffer%ie = ie
    buffer%js = js
    buffer%je = je

    do i=1,nlen
      call buffer%store(test_2d(i,:,:), i*3)
      slot = buffer%find_buffer_slot(i*3)
      local_fail = local_fail .or. ANY(buffer%buffer(slot)%field(:,:) /= test_2d(i,:,:))
    enddo

    if (verbose) write(stdout,*) "store_buffer_2d: ", local_fail
  end function store_buffer_2d

  !> Test the reuse of a buffer. Fill it first like store_buffer_2d. Then,
  !! loop through again, but use the slots of the buffer in the following
  !! order: 2, 1, 3
  function reuse_buffer_2d() result(local_fail)
    type(diag_buffer_2d) :: buffer
    logical :: local_fail !< True if any of the unit tests fail

    integer, parameter :: is=1, ie=2, js=3, je=6, nlen=3
    integer :: i, new_i, id, new_id
    real, dimension(nlen, is:ie, js:je) :: test_2d_first, test_2d_second
    integer, dimension(nlen) :: reorder = [2,1,3]

    local_fail = .false.
    call random_number(test_2d_first)
    call random_number(test_2d_second)

    call buffer%set_horizontal_extents(is=is, ie=ie, js=js, je=je)

    do i=1,nlen
      call buffer%store(test_2d_first(i,:,:), id=i*3)
    enddo

    do i=1,nlen
      new_i = reorder(i)
      ! id and new_id are multiplied by primes to make sure they are unique
      id = reorder(i)*3
      new_id = i*7
      call buffer%mark_available(id=reorder(i)*3)
      call buffer%store(test_2d_second(i,:,:), id=new_id)
      local_fail = local_fail .or. buffer%find_buffer_slot(new_id) /= new_i
      test_2d_first(new_i,:,:) = test_2d_second(i,:,:)
    enddo
    local_fail = local_fail .or. any(buffer%ids /= [14, 7, 21])
    do i=1,nlen
      local_fail = local_fail .or. any(buffer%buffer(i)%field(:,:) /= test_2d_first(i,:,:))
    enddo
    if (verbose) write(stdout,*) "reuse_buffer_2d: ", local_fail
  end function reuse_buffer_2d

end procedure diag_buffer_unit_tests_2d
module procedure diag_buffer_unit_tests_3d
  fail = .false.
  write(stdout,*) '==== MOM_diag_buffers: diag_buffers_unit_tests_3d ==='
  fail = fail .or. new_buffer_3d()
  fail = fail .or. grow_buffer_3d()
  fail = fail .or. fill_value_3d()
  fail = fail .or. store_buffer_3d()
  fail = fail .or. reuse_buffer_3d()

  contains

  !> Ensure properties of a newly initialized buffer
  function new_buffer_3d() result(local_fail)
    type(diag_buffer_3d) :: buffer
    logical :: local_fail !< True if any of the unit tests fail
    local_fail = .false.
    local_fail = local_fail .or. allocated(buffer%buffer)
    local_fail = local_fail .or. allocated(buffer%ids)
    local_fail = local_fail .or. buffer%length /= 0
    if (verbose) write(stdout,*) "new_buffer_3d: ", local_fail
  end function new_buffer_3d

  !> Test the growing of a buffer
  function grow_buffer_3d() result(local_fail)
    type(diag_buffer_3d) :: buffer
    logical :: local_fail !< True if any of the unit tests fail
    integer, parameter :: is=1, ie=2, js=3, je=6, ks=1, ke=10
    integer :: i

    local_fail = .false.

    call buffer%set_horizontal_extents(is=is, ie=ie, js=js, je=je)
    call buffer%set_vertical_extent(ks=ks, ke=ke)
    ! Grow the buffer 3 times
    do i=1,3
      call buffer%grow()
      local_fail = local_fail .or. (buffer%length /= i)
      local_fail = local_fail .or. (lbound(buffer%buffer(i)%field, 1) /= is)
      local_fail = local_fail .or. (ubound(buffer%buffer(i)%field, 1) /= ie)
      local_fail = local_fail .or. (lbound(buffer%buffer(i)%field, 2) /= js)
      local_fail = local_fail .or. (ubound(buffer%buffer(i)%field, 2) /= je)
      local_fail = local_fail .or. (lbound(buffer%buffer(i)%field, 3) /= ks)
      local_fail = local_fail .or. (ubound(buffer%buffer(i)%field, 3) /= ke)
    if (verbose) write(stdout,*) "grow_buffer_3d: ", local_fail
    enddo
    if (verbose) write(stdout,*) "grow_buffer_3d: ", local_fail
  end function grow_buffer_3d

  !> Test that growing new buffer fills the array with a set fill value
  function fill_value_3d() result(local_fail)
    type(diag_buffer_3d) :: buffer
    logical :: local_fail !< True if any of the unit tests fail
    integer, parameter :: is=1, ie=2, js=3, je=6, ks=1, ke=10
    real, parameter :: fill_value = -123.456

    local_fail = .false.

    call buffer%set_horizontal_extents(is=is, ie=ie, js=js, je=je)
    call buffer%set_vertical_extent(ks=ks, ke=ke)
    call buffer%set_fill_value(fill_value)
    call buffer%grow()
    if (any(buffer%buffer(1)%field(:,:,:) /= fill_value)) local_fail = .true.
    if (verbose) write(stdout,*) "fill_value_3d: ", local_fail
  end function fill_value_3d

  !> Test storing a buffer based on a unique id
  function store_buffer_3d() result(local_fail)
    type(diag_buffer_3d) :: buffer
    logical :: local_fail !< True if any of the unit tests fail

    integer, parameter :: is=1, ie=2, js=3, je=6, ks=1, ke=10, nlen=3
    integer :: i, slot
    real, dimension(nlen,is:ie,js:je,ks:ke) :: test_3d

    local_fail = .false.
    call random_number(test_3d)
    buffer%is = is
    buffer%ie = ie
    buffer%js = js
    buffer%je = je
    buffer%ks = ks
    buffer%ke = ke

    do i=1,nlen
      call buffer%store(test_3d(i,:,:,:), i*3)
      slot = buffer%find_buffer_slot(i*3)
      local_fail = local_fail .or. ANY(buffer%buffer(slot)%field(:,:,:) /= test_3d(i,:,:,:))
    enddo

    if (verbose) write(stdout,*) "store_buffer_3d: ", local_fail
  end function store_buffer_3d

  !> Test the reuse of a buffer. Fill it first like store_buffer_3d. Then,
  !! loop through again, but use the slots of the buffer in the following
  !! order: 2, 1, 3
  function reuse_buffer_3d() result(local_fail)
    type(diag_buffer_3d) :: buffer
    logical :: local_fail !< True if any of the unit tests fail

    integer, parameter :: is=1, ie=2, js=3, je=6, ks=1, ke=10, nlen=3
    integer :: i, new_i, id, new_id
    real, dimension(nlen, is:ie, js:je, ks:ke) :: test_3d_first, test_3d_second
    integer, dimension(nlen) :: reorder = [2,1,3]

    local_fail = .false.
    call random_number(test_3d_first)
    call random_number(test_3d_second)

    buffer%is = is
    buffer%ie = ie
    buffer%js = js
    buffer%je = je
    buffer%ks = ks
    buffer%ke = ke

    do i=1,nlen
      call buffer%store(test_3d_first(i,:,:,:), id=i*3)
    enddo

    do i=1,nlen
      new_i = reorder(i)
      ! id and new_id are multiplied by primes to make sure they are unique
      id = reorder(i)*3
      new_id = i*7
      call buffer%mark_available(id=reorder(i)*3)
      call buffer%store(test_3d_second(i,:,:,:), id=new_id)
      local_fail = local_fail .or. buffer%find_buffer_slot(new_id) /= new_i
      test_3d_first(new_i,:,:,:) = test_3d_second(i,:,:,:)
    enddo
    local_fail = local_fail .or. any(buffer%ids /= [14, 7, 21])
    do i=1,nlen
      local_fail = local_fail .or. any(buffer%buffer(i)%field(:,:,:) /= test_3d_first(i,:,:,:))
    enddo
    if (verbose) write(stdout,*) "reuse_buffer_3d: ", local_fail
  end function reuse_buffer_3d

end procedure diag_buffer_unit_tests_3d
end submodule MOM_diag_buffers_s
