! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> A simple type for keeping track of numerical tests
module numerical_testing_type

implicit none ; private

public testing
public numerical_testing_type_unit_tests

!> Class to assist in unit tests, not to be used outside of Recon1d types
type :: testing
  private
  !> True if any fail has been encountered since this instance of "testing" was created
  logical :: state = .false.
  !> Count of tests checked
  integer :: num_tests_checked = 0
  !> Count of tests failed
  integer :: num_tests_failed = 0
  !> If true, be verbose and write results to stdout. Default True.
  logical :: verbose = .true.
  !> Error channel
  integer, public :: stderr = 0
  !> Standard output channel
  integer, public :: stdout = 6
  !> If true, stop instantly
  logical :: stop_instantly = .false.
  !> If true, ignore fails until ignore_fail=.false.
  logical :: ignore_fail = .false.
  !> Record instances that fail
  integer :: ifailed(100) = 0.
  !> Record label of first instance that failed
  character(len=:), allocatable :: label_first_fail

  contains
    procedure :: test => test           !< Update the testing state
    procedure :: set => set             !< Set attributes
    procedure :: summarize => summarize !< Summarize testing state
    procedure :: real_scalar => real_scalar !< Compare two reals
    procedure :: real_arr => real_arr   !< Compare array of reals
    procedure :: int_arr => int_arr     !< Compare array of integers
end type


  interface
module subroutine test(this, state, label, ignore)
  class(testing),    intent(inout) :: this  !< This testing class
  logical,           intent(in)    :: state !< True to indicate a fail, false otherwise
  character(len=*),  intent(in)    :: label !< Message
  logical, optional, intent(in)    :: ignore !< If present and true, ignore a fail
  ! Local variables

end subroutine test
module subroutine set(this, verbose, stdout, stderr, stop_instantly, ignore_fail)
  class(testing), intent(inout) :: this  !< This testing class
  logical, optional, intent(in) :: verbose !< True or false setting to assign to verbosity
  integer, optional, intent(in) :: stdout !< The stdout channel to use
  integer, optional, intent(in) :: stderr !< The stderr channel to use
  logical, optional, intent(in) :: stop_instantly !< If true, stop immediately on error detection
  logical, optional, intent(in) :: ignore_fail !< If true, ignore fails until this option is set false

end subroutine set
logical module function summarize(this, label)
  class(testing),  intent(inout) :: this  !< This testing class
  character(len=*),   intent(in) :: label !< Message

end function summarize
module subroutine real_scalar(this, u_test, u_true, label, tol, robits, ignore)
  class(testing),  intent(inout) :: this   !< This testing class
  real,               intent(in) :: u_test !< Value to test [A]
  real,               intent(in) :: u_true !< Value to test against (correct answer) [A]
  character(len=*),   intent(in) :: label  !< Message
  real,     optional, intent(in) :: tol    !< The tolerance for differences between u and u_true [A]
  integer,  optional, intent(in) :: robits !< Number of bits of round-off to allow
  logical,  optional, intent(in) :: ignore !< If present and true, ignore a fail
  ! Local variables

end subroutine real_scalar
module subroutine real_arr(this, n, u_test, u_true, label, tol, robits, ignore)
  class(testing),  intent(inout) :: this   !< This testing class
  integer,            intent(in) :: n      !< Number of cells in u
  real, dimension(n), intent(in) :: u_test !< Values to test [A]
  real, dimension(n), intent(in) :: u_true !< Values to test against (correct answer) [A]
  character(len=*),   intent(in) :: label  !< Message
  real,     optional, intent(in) :: tol    !< The tolerance for differences between u and u_true [A]
  integer,  optional, intent(in) :: robits !< Number of bits of round-off to allow
  logical,  optional, intent(in) :: ignore !< If present and true, ignore a fail
  ! Local variables

end subroutine real_arr
module subroutine int_arr(this, n, i_test, i_true, label, ignore)
  class(testing),     intent(inout) :: this   !< This testing class
  integer,               intent(in) :: n      !< Number of cells in u
  integer, dimension(n), intent(in) :: i_test !< Values to test [A]
  integer, dimension(n), intent(in) :: i_true !< Values to test against (correct answer) [A]
  character(len=*),      intent(in) :: label  !< Message
  logical,  optional,    intent(in) :: ignore !< If present and true, ignore a fail
  ! Local variables

end subroutine int_arr
logical module function numerical_testing_type_unit_tests(verbose)
  logical, intent(in) :: verbose !< If true, write results to stdout
  ! Local variables

end function numerical_testing_type_unit_tests
  end interface

end module numerical_testing_type
