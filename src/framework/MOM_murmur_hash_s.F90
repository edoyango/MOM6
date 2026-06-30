submodule (MOM_murmur_hash) MOM_murmur_hash_s
  implicit none
contains
module procedure murmurhash3_i32
  integer(int32), parameter :: c1 = int(z'cc9e2d51', kind=int32)
  integer(int32), parameter :: c2 = int(z'1b873593', kind=int32)
  integer(int32), parameter :: c3 = int(z'e6546b64', kind=int32)
  integer(int32), parameter :: c4 = int(z'85ebca6b', kind=int32)
  integer(int32), parameter :: c5 = int(z'c2b2ae35', kind=int32)
  integer :: i
  integer(int32) :: k
  hash = 0
  if (present(seed)) hash = seed

  do i = 1, size(key)
    k = key(i)
    k = k * c1
    k = ishftc(k, 15)
    k = k * c2

    hash = ieor(hash, k)
    hash = ishftc(hash, 13)
    hash = 5 * hash + c3
  enddo

  ! NOTE: This is the point where the algorithm would handle trailing bytes.
  ! Since our arrays are comprised of 4 or 8 byte elements, we skip this part.

  hash = ieor(hash, 4*size(key))

  hash = ieor(hash, ishft(hash, -16))
  hash = hash * c4
  hash = ieor(hash, ishft(hash, -13))
  hash = hash * c5
  hash = ieor(hash, ishft(hash, -16))
end procedure murmurhash3_i32
module procedure murmurhash3_i64
  integer(int32) :: ikey(2*size(key))
  hash = murmur_hash(transfer(key, ikey), seed=seed)
end procedure murmurhash3_i64
module procedure murmurhash3_r32
  integer(int32) :: ikey(1)
  hash = murmur_hash(transfer(key, ikey), seed=seed)
end procedure murmurhash3_r32
module procedure murmurhash3_r32_1d
  integer(int32) :: ikey(size(key))
  hash = murmur_hash(transfer(key, ikey), seed=seed)
end procedure murmurhash3_r32_1d
module procedure murmurhash3_r32_2d
  integer(int32) :: ikey(size(key))
  hash = murmur_hash(transfer(key, ikey), seed=seed)
end procedure murmurhash3_r32_2d
module procedure murmurhash3_r32_3d
  integer(int32) :: ikey(size(key))
  hash = murmur_hash(transfer(key, ikey), seed=seed)
end procedure murmurhash3_r32_3d
module procedure murmurhash3_r32_4d
  integer(int32) :: ikey(size(key))
  hash = murmur_hash(transfer(key, ikey), seed=seed)
end procedure murmurhash3_r32_4d
module procedure murmurhash3_r64
  integer(int32) :: ikey(2)
  hash = murmur_hash(transfer(key, ikey), seed=seed)
end procedure murmurhash3_r64
module procedure murmurhash3_r64_1d
  integer(int32) :: ikey(2*size(key))
  hash = murmur_hash(transfer(key, ikey), seed=seed)
end procedure murmurhash3_r64_1d
module procedure murmurhash3_r64_2d
  integer(int32) :: ikey(2*size(key))
  hash = murmur_hash(transfer(key, ikey), seed=seed)
end procedure murmurhash3_r64_2d
module procedure murmurhash3_r64_3d
  integer(int32) :: ikey(2*size(key))
  hash = murmur_hash(transfer(key, ikey), seed=seed)
end procedure murmurhash3_r64_3d
module procedure murmurhash3_r64_4d
  integer(int32) :: ikey(2*size(key))
  hash = murmur_hash(transfer(key, ikey), seed=seed)
end procedure murmurhash3_r64_4d
end submodule MOM_murmur_hash_s
