! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

module MOM_unit_testing

use posix, only : chmod
use posix, only : sigsetjmp
use posix, only : sigjmp_buf

use MOM_coms, only : num_PEs, sync_PEs
use MOM_error_handler, only : is_root_pe
use MOM_error_handler, only : disable_fatal_errors
use MOM_error_handler, only : enable_fatal_errors

implicit none ; private

public :: string
public :: create_test_file
public :: delete_test_file
public :: TestSuite


!> String container type
type :: string
  character(len=:), allocatable :: s
    !< Internal character array of string
end type string


!> String constructor
interface string
  module procedure init_string_char
  module procedure init_string_int
end interface string


!> A generalized instance of a unit test function
type :: UnitTest
  private
  procedure(), nopass, pointer :: proc => null()
    !< Unit test function/subroutine
  procedure(), nopass, pointer :: cleanup => null()
    !< Cleanup function to be run after proc
  character(len=:), allocatable :: name
    !< Unit test name (usually set to name of proc)
  logical :: is_fatal
    !< True if proc() is expected to fail
contains
  procedure :: run => run_unit_test
    !< Run the unit test function, proc
end type UnitTest


!> Unit test constructor
interface UnitTest
  module procedure create_unit_test_basic
  module procedure create_unit_test_full
end interface UnitTest


!> Collection of unit tests
type :: TestSuite
  private
  type(UnitTestNode), pointer :: head => null()
    !< Head of the unit test linked list
  type(UnitTestNode), pointer :: tail => null()
    !< Tail of the unit test linked list (pre-allocated and unconfigured)

  ! Public API
  procedure(), nopass, pointer, public :: cleanup => null()
    !< Default cleanup function for unit tests in suite
contains
  private
  procedure :: add_basic => add_unit_test_basic
    !< Add a unit test without a cleanup function
  procedure :: add_full => add_unit_test_full
    !< Add a unit test with an explicit cleanup function
  generic, public :: add => add_basic, add_full
    !< Add a unit test to the test suite
  procedure, public :: run => run_test_suite
    !< Run all unit tests in the suite
end type TestSuite


!> TestSuite constructor
interface TestSuite
  module procedure create_test_suite
end interface TestSuite


!> UnitTest node of TestSuite's linked list
type :: UnitTestNode
  private
  type(UnitTest), pointer :: test => null()
    !< Node contents
  type(UnitTestNode), pointer :: next => null()
    !< Pointer to next node in list
end type UnitTestNode


  interface
module function create_unit_test_basic(proc, name, fatal) result(test)
  procedure() :: proc
    !< Subroutine which defines the unit test
  character(len=*), intent(in) :: name
    !< Name of the unit test
  logical, intent(in), optional :: fatal
    !< True if the test is expected to raise a FATAL error
  type(UnitTest) :: test

end function create_unit_test_basic
module function create_unit_test_full(proc, name, fatal, cleanup) result(test)
  procedure() :: proc
    !< Subroutine which defines the unit test
  character(len=*), intent(in) :: name
    !< Name of the unit test
  logical, optional :: fatal
    !< True if the test is expected to raise a FATAL error
  procedure() :: cleanup
    !< Cleanup subroutine, called after test
  type(UnitTest) :: test

end function create_unit_test_full
module subroutine run_unit_test(test)
  class(UnitTest), intent(in) :: test


end subroutine run_unit_test
module function create_test_suite() result(suite)
  type(TestSuite) :: suite

  ! Setup the head node, but do not populate it
end function create_test_suite
module subroutine add_unit_test_basic(suite, test, name, fatal)
  class(TestSuite), intent(inout) :: suite
  procedure() :: test
  character(len=*), intent(in) :: name
  logical, intent(in), optional :: fatal


end subroutine add_unit_test_basic
module subroutine add_unit_test_full(suite, test, name, fatal, cleanup)
  class(TestSuite), intent(inout) :: suite
  procedure() :: test
  character(len=*), intent(in) :: name
  procedure() :: cleanup
  logical, intent(in), optional :: fatal


  ! Populate the current tail
end subroutine add_unit_test_full
module subroutine run_test_suite(suite)
  class(TestSuite), intent(in) :: suite


end subroutine run_test_suite
module function init_string_char(c) result(str)
  character(len=*), dimension(:), intent(in) :: c
    !< List of character arrays
  type(string), dimension(size(c)) :: str
    !< String output


end function init_string_char
module function init_string_int(n) result(str)
  integer, intent(in) :: n
    !< Integer input
  type(string) :: str
    !< String output

  ! TODO: Estimate this with integer arithmetic

end function init_string_int
module subroutine create_test_file(filename, lines, mode)
  character(len=*), intent(in) :: filename
    !< Name of file to be created
  type(string), intent(in), optional :: lines(:)
    !< list of strings to write to file
  integer, optional, intent(in) :: mode
    !< Permissions of new file


end subroutine create_test_file
module subroutine delete_test_file(filename)
  character(len=*), intent(in) :: filename
    !< Name of file to be deleted


end subroutine delete_test_file
  end interface

end module MOM_unit_testing
