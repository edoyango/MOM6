! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> Provides gridded random number capability
module MOM_random

use MOM_hor_index,    only : hor_index_type
use MOM_time_manager, only : time_type, set_date, get_date

use iso_fortran_env,  only : stdout=>output_unit, stderr=>error_unit
use iso_fortran_env, only : int32

implicit none ; private

public :: random_0d_constructor
public :: random_01
public :: random_01_CB
public :: random_norm
public :: random_2d_constructor
public :: random_2d_01
public :: random_2d_norm
public :: random_unit_tests

! Private period parameters for the Mersenne Twister
integer, parameter :: &
    blockSize = 624,          & !< Size of the state vector
    M         = 397,          & !< Pivot element in state vector
    MATRIX_A  = -1727483681,  & !< constant vector a (0x9908b0dfUL)
    UMASK     = ibset(0, 31),  & !< most significant w-r bits (0x80000000UL)
    LMASK     = 2147483647      !< least significant r bits (0x7fffffffUL)

! Private tempering parameters for the Mersenne Twister
integer, parameter :: TMASKB= -1658038656, & !< (0x9d2c5680UL)
                      TMASKC= -272236544     !< (0xefc60000UL)

!> A private type used by the Mersenne Twistor
type randomNumberSequence
  integer                            :: currentElement !< Index into state vector
  integer, dimension(0:blockSize -1) :: state          !< State vector
end type randomNumberSequence

!> Container for pseudo-random number generators
type, public :: PRNG ; private

  !> Scalar random number generator for whole model
  type(randomNumberSequence) :: stream0d

  !> Random number generator for each cell on horizontal grid
  type(randomNumberSequence), dimension(:,:), allocatable :: stream2d

end type PRNG


  interface
real module function random_01(CS)
  type(PRNG), intent(inout) :: CS !< Container for pseudo-random number generators

end function random_01
real module function random_01_CB(ctr, key)
  use iso_fortran_env, only : int64
  integer, intent(in)  :: ctr !< ctr should be incremented each time you call the function
  integer, intent(in)  :: key !< key is like a seed: use a different key for each random stream

end function random_01_CB
real module function random_norm(CS)
  type(PRNG), intent(inout) :: CS !< Container for pseudo-random number generators
  ! Local variables

end function random_norm
module subroutine random_2d_01(CS, HI, rand)
  type(PRNG),           intent(inout) :: CS !< Container for pseudo-random number generators
  type(hor_index_type), intent(in)    :: HI !< Horizontal index structure
  real, dimension(HI%isd:HI%ied,HI%jsd:HI%jed), intent(out) :: rand !< Random numbers between 0 and 1 [nondim]
  ! Local variables

end subroutine random_2d_01
module subroutine random_2d_norm(CS, HI, rand)
  type(PRNG),           intent(inout) :: CS !< Container for pseudo-random number generators
  type(hor_index_type), intent(in)    :: HI !< Horizontal index structure
  real, dimension(HI%isd:HI%ied,HI%jsd:HI%jed), intent(out) :: rand !< Random numbers between 0 and 1 [nondim]
  ! Local variables

end subroutine random_2d_norm
module subroutine random_0d_constructor(CS, Time, seed)
  type(PRNG),      intent(inout) :: CS   !< Container for pseudo-random number generators
  type(time_type), intent(in)    :: Time !< Current model time
  integer,         intent(in)    :: seed !< Seed for PRNG
  ! Local variables

end subroutine random_0d_constructor
module subroutine random_2d_constructor(CS, HI, Time, seed)
  type(PRNG),           intent(inout) :: CS   !< Container for pseudo-random number generators
  type(hor_index_type), intent(in)    :: HI   !< Horizontal index structure
  type(time_type),      intent(in)    :: Time !< Current model time
  integer,              intent(in)    :: seed !< Seed for PRNG
  ! Local variables

end subroutine random_2d_constructor
integer module function seed_from_time(Time)
  type(time_type), intent(in)    :: Time !< Current model time
  ! Local variables

end function seed_from_time
integer module function seed_from_index(HI, i, j)
  type(hor_index_type), intent(in) :: HI !< Horizontal index structure
  integer,              intent(in) :: i !< i-index (of h-cell)
  integer,              intent(in) :: j !< j-index (of h-cell)
  ! Local variables

end function seed_from_index
module subroutine random_destruct(CS)
  type(PRNG), pointer :: CS !< Container for pseudo-random number generators

end subroutine random_destruct
module function new_RandomNumberSequence(seed) result(twister)
  integer, intent(in) :: seed !< Seed to initialize twister
  type(randomNumberSequence) :: twister !< The Mersenne Twister container
  ! Local variables

end function new_RandomNumberSequence
integer module function getRandomInt(twister)
  type(randomNumberSequence), intent(inout) :: twister !< The Mersenne Twister container

end function getRandomInt
double precision module function getRandomReal(twister)
  type(randomNumberSequence), intent(inout) :: twister
  ! Local variables

end function getRandomReal
integer module function mixbits(u, v)
  integer, intent(in) :: u !< An integer
  integer, intent(in) :: v !< An integer

end function mixbits
integer module function twist(u, v)
  integer, intent(in) :: u !< An integer
  integer, intent(in) :: v !< An integer
  ! Local variable

end function twist
module subroutine nextState(twister)
  type(randomNumberSequence), intent(inout) :: twister !< Container for the Mersenne Twister
  ! Local variables

end subroutine nextState
elemental integer module function temper(y)
  integer, intent(in) :: y !< An integer
  ! Local variables

end function temper
logical module function random_unit_tests(verbose)
  logical :: verbose !< True if results should be written to stdout
  ! Local variables
  ! Fake being on a decomposed domain

  ! Fake a decomposed domain
end function random_unit_tests
logical module function test_fn(verbose, good, label, rvalue, ivalue)
  logical,          intent(in) :: verbose !< Verbosity
  logical,          intent(in) :: good !< True if pass, false otherwise
  character(len=*), intent(in) :: label !< Label for messages
  real,             intent(in) :: rvalue !< Result of calculation [nondim]
  integer,          intent(in) :: ivalue !< Result of calculation
  optional :: rvalue, ivalue

end function test_fn
  end interface

end module MOM_random
