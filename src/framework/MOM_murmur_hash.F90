! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> MurmurHash is a non-cryptographic hash function developed by Austin Appleby.
!!
!! This module provides an implementation of the 32-bit MurmurHash3 algorithm.
!! It is used in MOM6 to generate unique hashes of field arrays.  The hash is
!! sensitive to order of elements and can detect changes that would otherwise
!! be missed by the mean/min/max/bitcount tests.
!!
!! Sensitivity to order means that it must be used with care for tests such as
!! processor layout.
!!
!! This implementation assumes data sizes of either 32 or 64 bits.  It cannot
!! be used for smaller types such as strings.
!!
!! https://github.com/aappleby/smhasher
module MOM_murmur_hash

use, intrinsic :: iso_fortran_env, only : int32, int64, real32, real64

implicit none ; private

public :: murmur_hash

!> Return the murmur3 hash of an array.
interface murmur_hash
  procedure murmurhash3_i32
  procedure murmurhash3_i64
  procedure murmurhash3_r32
  procedure murmurhash3_r32_1d
  procedure murmurhash3_r32_2d
  procedure murmurhash3_r32_3d
  procedure murmurhash3_r32_4d
  procedure murmurhash3_r64
  procedure murmurhash3_r64_1d
  procedure murmurhash3_r64_2d
  procedure murmurhash3_r64_3d
  procedure murmurhash3_r64_4d
end interface murmur_hash


  interface
module function murmurhash3_i32(key, seed) result(hash)
  integer(int32), intent(in) :: key(:)
    !< Input array
  integer(int32), intent(in), optional :: seed
    !< Hash seed
  integer(int32) :: hash
    !< Murmur hash of array




end function murmurhash3_i32
module function murmurhash3_i64(key, seed) result(hash)
  integer(int64), intent(in) :: key(:)
    !< Input array
  integer(int32), intent(in), optional :: seed
    !< Hash seed
  integer(int32) :: hash
    !< Murmur hash of array


end function murmurhash3_i64
module function murmurhash3_r32(key, seed) result(hash)
  real(real32), intent(in) :: key
    !< Input array [arbitrary]
  integer(int32), intent(in), optional :: seed
    !< Hash seed
  integer(int32) :: hash
    !< Murmur hash of array


end function murmurhash3_r32
module function murmurhash3_r32_1d(key, seed) result(hash)
  real(real32), intent(in) :: key(:)
    !< Input array [arbitrary]
  integer(int32), intent(in), optional :: seed
    !< Hash seed
  integer(int32) :: hash
    !< Murmur hash of array


end function murmurhash3_r32_1d
module function murmurhash3_r32_2d(key, seed) result(hash)
  real(real32), intent(in) :: key(:,:)
    !< Input array [arbitrary]
  integer(int32), intent(in), optional :: seed
    !< Hash seed
  integer(int32) :: hash
    !< Murmur hash of array


end function murmurhash3_r32_2d
module function murmurhash3_r32_3d(key, seed) result(hash)
  real(real32), intent(in) :: key(:,:,:)
    !< Input array [arbitrary]
  integer(int32), intent(in), optional :: seed
    !< Hash seed
  integer(int32) :: hash
    !< Murmur hash of array


end function murmurhash3_r32_3d
module function murmurhash3_r32_4d(key, seed) result(hash)
  real(real32), intent(in) :: key(:,:,:,:)
    !< Input array [arbitrary]
  integer(int32), intent(in), optional :: seed
    !< Hash seed
  integer(int32) :: hash
    !< Murmur hash of array


end function murmurhash3_r32_4d
module function murmurhash3_r64(key, seed) result(hash)
  real(real64), intent(in) :: key
    !< Input array [arbitrary]
  integer(int32), intent(in), optional :: seed
    !< Hash seed
  integer(int32) :: hash
    !< Murmur hash of array


end function murmurhash3_r64
module function murmurhash3_r64_1d(key, seed) result(hash)
  real(real64), intent(in) :: key(:)
    !< Input array [arbitrary]
  integer(int32), intent(in), optional :: seed
    !< Hash seed
  integer(int32) :: hash
    !< Murmur hash of array


end function murmurhash3_r64_1d
module function murmurhash3_r64_2d(key, seed) result(hash)
  real(real64), intent(in) :: key(:,:)
    !< Input array [arbitrary]
  integer(int32), intent(in), optional :: seed
    !< Hash seed
  integer(int32) :: hash
    !< Murmur hash of array


end function murmurhash3_r64_2d
module function murmurhash3_r64_3d(key, seed) result(hash)
  real(real64), intent(in) :: key(:,:,:)
    !< Input array [arbitrary]
  integer(int32), intent(in), optional :: seed
    !< Hash seed
  integer(int32) :: hash
    !< Murmur hash of array


end function murmurhash3_r64_3d
module function murmurhash3_r64_4d(key, seed) result(hash)
  real(real64), intent(in) :: key(:,:,:,:)
    !< Input array [arbitrary]
  integer(int32), intent(in), optional :: seed
    !< Hash seed
  integer(int32) :: hash
    !< Murmur hash of array


end function murmurhash3_r64_4d
  end interface

end module MOM_murmur_hash
