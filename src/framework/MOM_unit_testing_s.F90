submodule (MOM_unit_testing) MOM_unit_testing_s
  implicit none
contains
module procedure create_unit_test_basic
  procedure(), pointer :: cleanup
  cleanup => null()

  test = create_unit_test_full(proc, name, fatal, cleanup)
end procedure create_unit_test_basic
module procedure create_unit_test_full
  test%proc => proc
  test%name = name
  test%is_fatal = .false.
  if (present(fatal)) test%is_fatal = fatal
  test%cleanup => cleanup
end procedure create_unit_test_full
module procedure run_unit_test
  type(sigjmp_buf) :: env
  integer :: rc
  call sync_PEs

  ! FIXME: Some FATAL tests under MPI are unable to recover after jumpback, so
  !   we disable these tests for now.
  if (test%is_fatal .and. num_PEs() > 1) return

  if (test%is_fatal) then
    rc = sigsetjmp(env, 1)
    if (rc == 0) then
      call disable_fatal_errors(env)
      call test%proc
    endif
    call enable_fatal_errors
  else
    call test%proc
  endif

  if (associated(test%cleanup)) call test%cleanup
end procedure run_unit_test
module procedure create_test_suite
  allocate(suite%head)
  suite%tail => suite%head
end procedure create_test_suite
module procedure add_unit_test_basic
  procedure(), pointer :: cleanup
  cleanup => null()
  if (associated(suite%cleanup)) cleanup => suite%cleanup

  call add_unit_test_full(suite, test, name, fatal, cleanup)
end procedure add_unit_test_basic
module procedure add_unit_test_full
  type(UnitTest), pointer :: utest
  type(UnitTestNode), pointer :: node
  allocate(utest)
  utest = UnitTest(test, name, fatal, cleanup)
  suite%tail%test => utest

  ! Create and append the new (empty) node, and update the tail
  allocate(node)
  suite%tail%next => node
  suite%tail => node
end procedure add_unit_test_full
module procedure run_test_suite
  type(UnitTestNode), pointer :: node
  node => suite%head
  do while(associated(node%test))
    ! TODO: Capture FMS stdout/stderr
    print '(/a)', "=== "//node%test%name

    call node%test%run
    if (associated(node%test%cleanup)) call node%test%cleanup

    node => node%next
  enddo
end procedure run_test_suite
module procedure init_string_char
  type(string), dimension(size(c)) :: str
  integer :: i
  do i = 1, size(c)
    str(i)%s = c(i)
  enddo
end procedure init_string_char
module procedure init_string_int
  character(1 + floor(log10(real(abs(n)))) + (1 - sign(1, n))/2) :: chr
  write(chr, '(i0)') n
  str = string(chr)
end procedure init_string_int
module procedure create_test_file
  integer :: param_unit
  integer :: i
  integer :: rc
  logical :: sync
  if (is_root_PE()) then
    open(newunit=param_unit, file=filename, status='replace')
    if (present(lines)) then
      do i = 1, size(lines)
        write(param_unit, '(a)') lines(i)%s
      enddo
    endif
    close(param_unit)
    if (present(mode)) rc = chmod(filename, mode)
  endif
  call sync_PEs
end procedure create_test_file
module procedure delete_test_file
  logical :: is_file, is_open
  integer :: io_unit
  if (is_root_PE()) then
    inquire(file=filename, exist=is_file, opened=is_open, number=io_unit)

    if (is_file) then
      if (.not. is_open) open(newunit=io_unit, file=filename)
      close(io_unit, status='delete')
    endif
  endif
  call sync_PEs
end procedure delete_test_file
end submodule MOM_unit_testing_s
