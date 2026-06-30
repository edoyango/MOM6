submodule (MOM_array_transform) MOM_array_transform_s
  implicit none
contains
module procedure rotate_array_real_2d
  integer :: m, n
  m = size(A_in, 1)
  n = size(A_in, 2)

  select case (modulo(turns, 4))
    case(0)
      A(:,:) = A_in(:,:)
    case(1)
      A(:,:) = transpose(A_in)
      A(:,:) = A(n:1:-1, :)
    case(2)
      A(:,:) = A_in(m:1:-1, n:1:-1)
    case(3)
      A(:,:) = transpose(A_in(m:1:-1, :))
  end select
end procedure rotate_array_real_2d
module procedure rotate_array_real_3d
  integer :: k
  do k = 1, size(A_in, 3)
    call rotate_array(A_in(:,:,k), turns, A(:,:,k))
  enddo
end procedure rotate_array_real_3d
module procedure rotate_array_real_4d
  integer :: n
  do n = 1, size(A_in, 4)
    call rotate_array(A_in(:,:,:,n), turns, A(:,:,:,n))
  enddo
end procedure rotate_array_real_4d
module procedure rotate_array_integer
  integer :: m, n
  m = size(A_in, 1)
  n = size(A_in, 2)

  select case (modulo(turns, 4))
    case(0)
      A(:,:) = A_in(:,:)
    case(1)
      A(:,:) = transpose(A_in)
      A(:,:) = A(n:1:-1, :)
    case(2)
      A(:,:) = A_in(m:1:-1, n:1:-1)
    case(3)
      A(:,:) = transpose(A_in(m:1:-1, :))
  end select
end procedure rotate_array_integer
module procedure rotate_array_logical
  integer :: m, n
  m = size(A_in, 1)
  n = size(A_in, 2)

  select case (modulo(turns, 4))
    case(0)
      A(:,:) = A_in(:,:)
    case(1)
      A(:,:) = transpose(A_in)
      A(:,:) = A(n:1:-1, :)
    case(2)
      A(:,:) = A_in(m:1:-1, n:1:-1)
    case(3)
      A(:,:) = transpose(A_in(m:1:-1, :))
  end select
end procedure rotate_array_logical
module procedure rotate_array_pair_real_2d
  if (modulo(turns, 2) /= 0) then
    call rotate_array(B_in, turns, A)
    call rotate_array(A_in, turns, B)
  else
    call rotate_array(A_in, turns, A)
    call rotate_array(B_in, turns, B)
  endif
end procedure rotate_array_pair_real_2d
module procedure rotate_array_pair_real_3d
  integer :: k
  do k = 1, size(A_in, 3)
    call rotate_array_pair(A_in(:,:,k), B_in(:,:,k), turns, &
        A(:,:,k), B(:,:,k))
  enddo
end procedure rotate_array_pair_real_3d
module procedure rotate_array_pair_integer
  if (modulo(turns, 2) /= 0) then
    call rotate_array(B_in, turns, A)
    call rotate_array(A_in, turns, B)
  else
    call rotate_array(A_in, turns, A)
    call rotate_array(B_in, turns, B)
  endif
end procedure rotate_array_pair_integer
module procedure rotate_vector_real_2d
  call rotate_array_pair(A_in, B_in, turns, A, B)

  if (modulo(turns, 4) == 1 .or. modulo(turns, 4) == 2) &
    A(:,:) = -A(:,:)

  if (modulo(turns, 4) == 2 .or. modulo(turns, 4) == 3) &
    B(:,:) = -B(:,:)
end procedure rotate_vector_real_2d
module procedure rotate_vector_real_3d
  integer :: k
  do k = 1, size(A_in, 3)
    call rotate_vector(A_in(:,:,k), B_in(:,:,k), turns, A(:,:,k), B(:,:,k))
  enddo
end procedure rotate_vector_real_3d
module procedure rotate_vector_real_4d
  integer :: n
  do n = 1, size(A_in, 4)
    call rotate_vector(A_in(:,:,:,n), B_in(:,:,:,n), turns, &
        A(:,:,:,n), B(:,:,:,n))
  enddo
end procedure rotate_vector_real_4d
module procedure allocate_rotated_array_real_2d
  integer :: ub(2)
  ub(:) = ubound(A_in)

  if (modulo(turns, 2) /= 0) then
    allocate(A(lb(2):ub(2), lb(1):ub(1)))
  else
    allocate(A(lb(1):ub(1), lb(2):ub(2)))
  endif
end procedure allocate_rotated_array_real_2d
module procedure allocate_rotated_array_real_3d
  integer :: ub(3)
  ub(:) = ubound(A_in)

  if (modulo(turns, 2) /= 0) then
    allocate(A(lb(2):ub(2), lb(1):ub(1), lb(3):ub(3)))
  else
    allocate(A(lb(1):ub(1), lb(2):ub(2), lb(3):ub(3)))
  endif
end procedure allocate_rotated_array_real_3d
module procedure allocate_rotated_array_real_4d
  integer:: ub(4)
  ub(:) = ubound(A_in)

  if (modulo(turns, 2) /= 0) then
    allocate(A(lb(2):ub(2), lb(1):ub(1), lb(3):ub(3), lb(4):ub(4)))
  else
    allocate(A(lb(1):ub(1), lb(2):ub(2), lb(3):ub(3), lb(4):ub(4)))
  endif
end procedure allocate_rotated_array_real_4d
module procedure allocate_rotated_array_integer
  integer :: ub(2)
  ub(:) = ubound(A_in)

  if (modulo(turns, 2) /= 0) then
    allocate(A(lb(2):ub(2), lb(1):ub(1)))
  else
    allocate(A(lb(1):ub(1), lb(2):ub(2)))
  endif
end procedure allocate_rotated_array_integer
module procedure symmetric_sum_1d
  integer :: i, szi, szi_2
  szi = size(field, 1)
  szi_2 = szi / 2 ! Note that for an odd number szi_2 is rounded down.
  sum = 0.0
  if (2*szi_2 < szi) sum = field(szi_2+1)
  ! Add pairs of values, working from the inside out.
  do i=szi_2,1,-1
    sum = sum + (field(i) + field(szi+1-i))
  enddo
end procedure symmetric_sum_1d
module procedure symmetric_sum_2d
  real :: quad_sum(2,2) ! The sums in each of the quadrants [A ~> a]
  logical :: odd_i, odd_j
  integer :: ij, szi, szj, szi_2, szj_2, ic, jc
  szi = size(field, 1) ; szj = size(field, 2)
  ! These 5 special cases are equivalent to the general case, but they reduce the use
  ! of complicated logic for common simple cases.
  if ((szi == 1) .and. (szj == 1)) then
    sum = field(1,1)
  elseif ((szi == 2) .and. (szj == 2)) then
    sum = (field(1,1) + field(2,2)) + (field(2,1) + field(1,2))
  elseif ((szi == 3) .and. (szj == 3)) then
    sum = (field(2,2) + ((field(1,2) + field(3,2)) + (field(2,1) + field(2,3)))) + &
          ((field(1,1) + field(3,3)) + (field(3,1) + field(1,3)))
  elseif (szi == 1) then
    sum = symmetric_sum_1d(field(1,:))
  elseif (szj == 1) then
    sum = symmetric_sum_1d(field(:,1))
  else
    ! This is the general case.
    ! Note that for odd numbers szi_2 and szj_2 are rounded down.
    szi_2 = szi / 2
    szj_2 = szj / 2

    odd_i = (2*szi_2 < szi) ! This could be (modulo(szi,2) == 1)
    odd_j = (2*szj_2 < szj)
    ! Start by finding the sums along the central axes if there are an odd number of points.
    if (odd_i .and. odd_j) then
      ic = szi_2+1 ; jc = szj_2+1 ! The index of the central point
      sum = field(ic,jc)
      ! Add pairs of pairs of values, working from the inside out.
      do ij=1,min(szi_2,szj_2)
        sum = sum + ((field(ic-ij,jc) + field(ic+ij,jc)) + (field(ic,jc-ij) + field(ic,jc+ij)))
      enddo
      ! Add extra pairs of values, working from the inside out.
      if (szi_2 > szj_2) then
        do ij=szj_2+1,szi_2
          sum = sum + (field(ic-ij,jc) + field(ic+ij,jc))
        enddo
      elseif (szj_2 > szi_2) then
        do ij=szi_2+1,szj_2
          sum = sum + (field(ic,jc-ij) + field(ic,jc+ij))
        enddo
      endif
    elseif (odd_i) then
      sum = symmetric_sum_1d(field(szi_2+1,1:szj))
    elseif (odd_j) then
      sum = symmetric_sum_1d(field(1:szi,szj_2+1))
    else
      sum = 0.0
    endif

    ! Find the sums in the four quadrants of the array.
    if ((szi_2 > 1) .and. (szj_2 > 1)) then
      ! Use a recursive call to symmetric_sum_2d to determine the sums in the corner quadrants.
      quad_sum(1,1) = symmetric_sum_2d(field(1:szi_2,1:szj_2))
      quad_sum(2,1) = symmetric_sum_2d(field(szi+1-szi_2:szi,1:szj_2))
      quad_sum(1,2) = symmetric_sum_2d(field(1:szi_2,szj+1-szj_2:szj))
      quad_sum(2,2) = symmetric_sum_2d(field(szi+1-szi_2:szi,szj+1-szj_2:szj))
    elseif (szi_2 > 1) then
      quad_sum(1,1) = symmetric_sum_1d(field(1:szi_2,1))
      quad_sum(2,1) = symmetric_sum_1d(field(szi+1-szi_2:szi,1))
      quad_sum(1,2) = symmetric_sum_1d(field(1:szi_2,szj))
      quad_sum(2,2) = symmetric_sum_1d(field(szi+1-szi_2:szi,szj))
    elseif (szj_2 > 1) then
      quad_sum(1,1) = symmetric_sum_1d(field(1,1:szj_2))
      quad_sum(2,1) = symmetric_sum_1d(field(szi,1:szj_2))
      quad_sum(1,2) = symmetric_sum_1d(field(1,szj+1-szj_2:szj))
      quad_sum(2,2) = symmetric_sum_1d(field(szi,szj+1-szj_2:szj))
    else
      quad_sum(1,1) = field(1,1)
      quad_sum(2,1) = field(szi,1)
      quad_sum(1,2) = field(1,szj)
      quad_sum(2,2) = field(szi,szj)
    endif

    sum = sum + ((quad_sum(1,1) + quad_sum(2,2)) + (quad_sum(2,1) + quad_sum(1,2)))
  endif
end procedure symmetric_sum_2d
module procedure naive_sum_2d
  logical :: sum_abs_val
  integer :: i, j, szi, szj
  szi = size(field, 1) ; szj = size(field, 2)
  sum_abs_val = .false. ; if (present(abs_val)) sum_abs_val = abs_val
  sum = 0.0
  if (sum_abs_val) then
    do j=1,szj ; do i=1,szi
      sum = sum + abs(field(i,j))
    enddo ; enddo
  else
    do j=1,szj ; do i=1,szi
      sum = sum + field(i,j)
    enddo ; enddo
  endif
end procedure naive_sum_2d
module procedure symmetric_sum_unit_tests
  character(len=120) :: fail_message !< Blank or a description of the first failed test.
  integer, parameter :: sz=13 ! The maximum size of the test arrays
  real :: array(sz,sz)  ! An array of inexact real values for testing in arbitrary units [A]
  real :: ar_90(sz,sz)  ! Array rotated by 90 degrees in arbitrary units [A]
  real :: ar_180(sz,sz) ! Array rotated by 180 degrees in arbitrary units [A]
  real :: ar_270(sz,sz) ! Array rotated by 270 degrees in arbitrary units [A]
  real :: sum(5)        ! Different versions of sums over a sub-array [A]
  real :: abs_sum       ! The sum of the absolute values of the array [A]
  real :: tol           ! The tolerance for an inexact test [A]
  character(len=120) :: mesg
  integer :: i, j, n, m, r
  logical :: fail
  fail = .false.
  fail_message = ""

  if (verbose) write(stdout,*) '==== MOM_array_transform: symmetric_sum_unit_tests ===='

  ! Fill the array with real numbers that can not be represented exactly.
  do j=1,sz ; do i=1,sz
    array(i,j) = 1.0 / (2.0*(j*sz + i) + 1.0)
    ! Combining positive and negative numbers amplifies differences from the order of arithmetic.
    if (modulo(i+j, 2) == 0) array(i,j) = -array(i,j)
  enddo ; enddo
  call rotate_array_real_2d(array, 1, ar_90)
  call rotate_array_real_2d(array, 2, ar_180)
  call rotate_array_real_2d(array, 3, ar_270)

  do n = 1, sz ; do m = 1, sz
    sum(1) = symmetric_sum(array(1:n,1:m))
    sum(2) = symmetric_sum(ar_90(sz+1-m:sz,1:n))
    sum(3) = symmetric_sum(ar_180(sz+1-n:sz,sz+1-m:sz))
    sum(4) = symmetric_sum(ar_270(1:m,sz+1-n:sz))
    sum(5) = naive_sum_2d(array(1:n,1:m))
    abs_sum = naive_sum_2d(array(1:n,1:m), abs_val=.true.)
    tol = 2.0 * abs_sum * epsilon(abs_sum)
    if (abs(sum(1) - sum(5)) > tol) then
      write(mesg,'(i0," x ",i0," symmetric vs naive sum, sum=",ES13.5," diff=",ES13.5)') &
            n, m, sum(1), sum(5) - sum(1)
      write(stdout,*) "Symmetric_sum_failure: "//trim(mesg)
      write(stderr,*) "Symmetric_sum_failure: "//trim(mesg)
      if (.not.fail) fail_message = mesg ! This is the first failed test.
      fail = .true.
    endif
    do r = 2, 4 ; if (abs(sum(1) - sum(r)) > 0.0) then
      write(mesg,'(i0," x ",i0," with ",i0," degree rotation, sum=",ES13.5," diff=",ES13.5)') &
            n, m, 90*(r-1), sum(1), sum(r) - sum(1)
      write(stdout,*) "Symmetric_sum_failure: "//trim(mesg)
      write(stderr,*) "Symmetric_sum_failure: "//trim(mesg)
      if (.not.fail) fail_message = mesg ! This is the first failed test.
      fail = .true.
    endif ; enddo
  enddo ; enddo

  if (fail) then
    write(stdout,*) "MOM_array_transform: One or more symmetric sum tests has failed."
    write(stderr,*) "MOM_array_transform: One or more symmetric sum tests has failed."
  else
    if (verbose) write(stdout,*) ("MOM_array_transform: All symmetric sum tests have passed.")
  endif
  symmetric_sum_unit_tests = fail

end procedure symmetric_sum_unit_tests
end submodule MOM_array_transform_s
