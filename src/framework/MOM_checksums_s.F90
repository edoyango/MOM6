submodule (MOM_checksums) MOM_checksums_s
  implicit none
contains
module procedure chksum0
  real :: scaling   !< Explicit rescaling factor [a A-1 ~> 1]
  integer :: iounit !< Log IO unit
  real :: rs        !< Rescaled scalar [a]
  integer :: bc     !< Scalar bitcount
  if (checkForNaNs .and. is_NaN(scalar)) &
    call chksum_error(FATAL, 'NaN detected: '//trim(mesg))

  scaling = 1.0
  if (present(unscale)) then ; scaling = unscale
  elseif (present(scale)) then ; scaling = scale ; endif

  iounit = error_unit ; if (present(logunit)) iounit = logunit

  if (calculateStatistics) then
    rs = scaling * scalar
    if (is_root_pe()) &
      call chk_sum_msg(" scalar:", rs, rs, rs, mesg, iounit)
  endif

  if (.not. writeChksums) return

  bc = mod(bitcount(abs(scaling * scalar)), bc_modulus)
  if (is_root_pe()) &
    call chk_sum_msg(" scalar:", bc, mesg, iounit)

  if (writeHash .and. is_root_pe()) &
    write(iounit, '(" scalar: hash=", z8, 1x, a)') &
        murmur_hash(scaling * scalar), mesg
end procedure chksum0
module procedure zchksum
  real, allocatable, dimension(:) :: rescaled_array ! The array with scaling undone [a]
  real :: scaling   ! Explicit rescaling factor [a A-1 ~> 1]
  integer :: iounit !< Log IO unit
  integer :: k
  real :: aMean, aMin, aMax  ! Array mean, global minimum and global maximum [a]
  integer :: bc0
  if (checkForNaNs) then
    if (is_NaN(array(:))) &
      call chksum_error(FATAL, 'NaN detected: '//trim(mesg))
  endif

  scaling = 1.0
  if (present(unscale)) then ; scaling = unscale
  elseif (present(scale)) then ; scaling = scale ; endif

  iounit = error_unit ; if (present(logunit)) iounit = logunit

  if (calculateStatistics) then
    if (present(unscale) .or. present(scale)) then
      allocate(rescaled_array(LBOUND(array,1):UBOUND(array,1)), source=0.0)
      do k=1, size(array, 1)
        rescaled_array(k) = scaling * array(k)
      enddo

      call subStats(rescaled_array, aMean, aMin, aMax)
      deallocate(rescaled_array)
    else
      call subStats(array, aMean, aMin, aMax)
    endif

    if (is_root_pe()) &
      call chk_sum_msg(" column:", aMean, aMin, aMax, mesg, iounit)
  endif

  if (.not. writeChksums) return

  bc0 = subchk(array, scaling)
  if (is_root_pe()) call chk_sum_msg(" column:", bc0, mesg, iounit)

  if (writeHash .and. is_root_pe()) &
    write(iounit, '(" column: hash=", z8, 1x, a)') &
        murmur_hash(scaling * array), mesg

  contains

  integer function subchk(array, unscale)
    real, dimension(:), intent(in) :: array !< The array to be checksummed in
                                            !! arbitrary, possibly rescaled units [A ~> a]
    real, intent(in) :: unscale !< A factor to convert this array back to unscaled units
                                !! for checksums and output [a A-1 ~> 1]
    integer :: k, bc
    subchk = 0
    do k=LBOUND(array, 1), UBOUND(array, 1)
      bc = bitcount(abs(unscale * array(k)))
      subchk = subchk + bc
    enddo
    subchk=mod(subchk, bc_modulus)
  end function subchk

  subroutine subStats(array, aMean, aMin, aMax)
    real, dimension(:), intent(in) :: array !< The array to be checksummed [a]
    real, intent(out) :: aMean !< Array mean [a]
    real, intent(out) :: aMin !< Array minimum [a]
    real, intent(out) :: aMax !< Array maximum [a]

    integer :: k, n

    aMin = array(1)
    aMax = array(1)
    n = 0
    do k=LBOUND(array,1), UBOUND(array,1)
      aMin = min(aMin, array(k))
      aMax = max(aMax, array(k))
      n = n + 1
    enddo
    aMean = sum(array(:)) / real(n)
  end subroutine subStats
end procedure zchksum
module procedure chksum_pair_h_2d
  logical :: vector_pair
  integer :: turns
  type(hor_index_type), pointer :: HI_in
  real, dimension(:,:), pointer :: arrayA_in, arrayB_in ! Rotated arrays [A ~> a]
  vector_pair = .true.
  if (present(scalar_pair)) vector_pair = .not. scalar_pair

  turns = HI%turns
  if (modulo(turns, 4) /= 0) then
    ! Rotate field back to the input grid
    allocate(HI_in)
    call rotate_hor_index(HI, -turns, HI_in)
    allocate(arrayA_in(HI_in%isd:HI_in%ied, HI_in%jsd:HI_in%jed))
    allocate(arrayB_in(HI_in%isd:HI_in%ied, HI_in%jsd:HI_in%jed))

    if (vector_pair) then
      call rotate_vector(arrayA, arrayB, -turns, arrayA_in, arrayB_in)
    else
      call rotate_array_pair(arrayA, arrayB, -turns, arrayA_in, arrayB_in)
    endif
  else
    HI_in => HI
    arrayA_in => arrayA
    arrayB_in => arrayB
  endif

  if (present(haloshift)) then
    call chksum_h_2d(arrayA_in, 'x '//mesg, HI_in, haloshift, omit_corners, &
                     scale=scale, logunit=logunit, unscale=unscale)
    call chksum_h_2d(arrayB_in, 'y '//mesg, HI_in, haloshift, omit_corners, &
                     scale=scale, logunit=logunit, unscale=unscale)
  else
    call chksum_h_2d(arrayA_in, 'x '//mesg, HI_in, scale=scale, logunit=logunit, unscale=unscale)
    call chksum_h_2d(arrayB_in, 'y '//mesg, HI_in, scale=scale, logunit=logunit, unscale=unscale)
  endif
end procedure chksum_pair_h_2d
module procedure chksum_pair_h_3d
  logical :: vector_pair
  integer :: turns
  type(hor_index_type), pointer :: HI_in
  real, dimension(:,:,:), pointer :: arrayA_in, arrayB_in ! Rotated arrays [A ~> a]
  vector_pair = .true.
  if (present(scalar_pair)) vector_pair = .not. scalar_pair

  turns = HI%turns
  if (modulo(turns, 4) /= 0) then
    ! Rotate field back to the input grid
    allocate(HI_in)
    call rotate_hor_index(HI, -turns, HI_in)
    allocate(arrayA_in(HI_in%isd:HI_in%ied, HI_in%jsd:HI_in%jed, size(arrayA, 3)))
    allocate(arrayB_in(HI_in%isd:HI_in%ied, HI_in%jsd:HI_in%jed, size(arrayB, 3)))

    if (vector_pair) then
      call rotate_vector(arrayA, arrayB, -turns, arrayA_in, arrayB_in)
    else
      call rotate_array_pair(arrayA, arrayB, -turns, arrayA_in, arrayB_in)
    endif
  else
    HI_in => HI
    arrayA_in => arrayA
    arrayB_in => arrayB
  endif

  if (present(haloshift)) then
    call chksum_h_3d(arrayA_in, 'x '//mesg, HI_in, haloshift, omit_corners, &
                     scale=scale, logunit=logunit, unscale=unscale)
    call chksum_h_3d(arrayB_in, 'y '//mesg, HI_in, haloshift, omit_corners, &
                     scale=scale, logunit=logunit, unscale=unscale)
  else
    call chksum_h_3d(arrayA_in, 'x '//mesg, HI_in, scale=scale, logunit=logunit, unscale=unscale)
    call chksum_h_3d(arrayB_in, 'y '//mesg, HI_in, scale=scale, logunit=logunit, unscale=unscale)
  endif

  ! NOTE: automatic deallocation of array[AB]_in
end procedure chksum_pair_h_3d
module procedure chksum_h_2d
  real, pointer :: array(:,:)           ! Field array on the input grid [A ~> a]
  real, allocatable, dimension(:,:) :: rescaled_array ! The array with scaling undone [a]
  real, allocatable :: hash_array(:,:)  ! Subarray used to compute hash [a]
  type(hor_index_type), pointer :: HI   ! Horizontal index bounds of the input grid
  real :: scaling   ! Explicit rescaling factor [a A-1 ~> 1]
  integer :: iounit !< Log IO unit
  integer :: i, j
  real :: aMean, aMin, aMax  ! Array mean, global minimum and global maximum [a]
  integer :: bc0, bcSW, bcSE, bcNW, bcNE, hshift
  integer :: bcN, bcS, bcE, bcW
  logical :: do_corners
  integer :: turns                      ! Quarter turns from input to model grid
  turns = HI_m%turns
  if (modulo(turns, 4) /= 0) then
    allocate(HI)
    call rotate_hor_index(HI_m, -turns, HI)
    allocate(array(HI%isd:HI%ied, HI%jsd:HI%jed))
    call rotate_array(array_m, -turns, array)
  else
    HI => HI_m
    array => array_m
  endif

  if (checkForNaNs) then
    if (is_NaN(array(HI%isc:HI%iec,HI%jsc:HI%jec))) &
      call chksum_error(FATAL, 'NaN detected: '//trim(mesg))
!   if (is_NaN(array)) &
!     call chksum_error(FATAL, 'NaN detected in halo: '//trim(mesg))
  endif

  scaling = 1.0
  if (present(unscale)) then ; scaling = unscale
  elseif (present(scale)) then ; scaling = scale ; endif

  iounit = error_unit ; if (present(logunit)) iounit = logunit

  if (calculateStatistics) then
    if (present(unscale) .or. present(scale)) then
      allocate( rescaled_array(LBOUND(array,1):UBOUND(array,1), &
                               LBOUND(array,2):UBOUND(array,2)), source=0.0 )
      do j=HI%jsc,HI%jec ; do i=HI%isc,HI%iec
        rescaled_array(i,j) = scaling*array(i,j)
      enddo ; enddo
      call subStats(HI, rescaled_array, aMean, aMin, aMax)
      deallocate(rescaled_array)
    else
      call subStats(HI, array, aMean, aMin, aMax)
    endif

    if (is_root_pe()) &
      call chk_sum_msg("h-point:", aMean, aMin, aMax, mesg, iounit)
  endif

  if (.not.writeChksums) return

  hshift = default_shift
  if (present(haloshift)) hshift = haloshift
  if (hshift<0) hshift = HI%ied-HI%iec

  if ( HI%isc-hshift<HI%isd .or. HI%iec+hshift>HI%ied .or. &
       HI%jsc-hshift<HI%jsd .or. HI%jec+hshift>HI%jed ) then
    write(0,*) 'chksum_h_2d: haloshift =',hshift
    write(0,*) 'chksum_h_2d: isd,isc,iec,ied=',HI%isd,HI%isc,HI%iec,HI%ied
    write(0,*) 'chksum_h_2d: jsd,jsc,jec,jed=',HI%jsd,HI%jsc,HI%jec,HI%jed
    call chksum_error(FATAL,'Error in chksum_h_2d '//trim(mesg))
  endif

  bc0 = subchk(array, HI, 0, 0, scaling)

  if (hshift==0) then
    if (is_root_pe()) call chk_sum_msg("h-point:", bc0, mesg, iounit)
  else
    do_corners = .true.
    if (present(omit_corners)) do_corners = .not. omit_corners

    if (do_corners) then
      bcSW = subchk(array, HI, -hshift, -hshift, scaling)
      bcSE = subchk(array, HI, hshift, -hshift, scaling)
      bcNW = subchk(array, HI, -hshift, hshift, scaling)
      bcNE = subchk(array, HI, hshift, hshift, scaling)

      if (is_root_pe()) &
        call chk_sum_msg("h-point:", bc0, bcSW, bcSE, bcNW, bcNE, mesg, iounit)
    else
      bcS = subchk(array, HI, 0, -hshift, scaling)
      bcE = subchk(array, HI, hshift, 0, scaling)
      bcW = subchk(array, HI, -hshift, 0, scaling)
      bcN = subchk(array, HI, 0, hshift, scaling)

      if (is_root_pe()) &
        call chk_sum_msg_NSEW("h-point:", bc0, bcN, bcS, bcE, bcW, mesg, iounit)
    endif
  endif

  if (writeHash .and. is_root_pe()) then
    allocate(hash_array(HI%isc:HI%iec, HI%jsc:HI%jec))
    hash_array(:,:) = scaling * array(HI%isc:HI%iec, HI%jsc:HI%jec)

    write(iounit, '("h-point: hash=", z8, 1x, a)') &
        murmur_hash(hash_array), mesg
    deallocate(hash_array)
  endif

  contains
  integer function subchk(array, HI, di, dj, unscale)
    type(hor_index_type), intent(in) ::  HI     !< A horizontal index type
    real, dimension(HI%isd:,HI%jsd:), intent(in) :: array !< The array to be checksummed in
                                 !! arbitrary, possibly rescaled units [A ~> a]
    integer, intent(in) :: di    !< i- direction array shift for this checksum
    integer, intent(in) :: dj    !< j- direction array shift for this checksum
    real, intent(in)    :: unscale !< A factor to convert this array back to unscaled units
                                 !! for checksums and output [a A-1 ~> 1]
    integer :: i, j, bc
    subchk = 0
    do j=HI%jsc+dj,HI%jec+dj ; do i=HI%isc+di,HI%iec+di
      bc = bitcount(abs(unscale*array(i,j)))
      subchk = subchk + bc
    enddo ; enddo
    call sum_across_PEs(subchk)
    subchk=mod(subchk, bc_modulus)
  end function subchk

  subroutine subStats(HI, array, aMean, aMin, aMax)
    type(hor_index_type), intent(in) ::  HI     !< A horizontal index type
    real, dimension(HI%isd:,HI%jsd:), intent(in) :: array !< The array to be checksummed [a]
    real, intent(out) :: aMean !< Array mean [a]
    real, intent(out) :: aMin !< Array minimum [a]
    real, intent(out) :: aMax !< Array maximum [a]

    integer :: i, j, n

    aMin = array(HI%isc,HI%jsc)
    aMax = array(HI%isc,HI%jsc)
    n = 0
    do j=HI%jsc,HI%jec ; do i=HI%isc,HI%iec
      aMin = min(aMin, array(i,j))
      aMax = max(aMax, array(i,j))
      n = n + 1
    enddo ; enddo
    aMean = reproducing_sum(array(HI%isc:HI%iec,HI%jsc:HI%jec))
    call sum_across_PEs(n)
    call min_across_PEs(aMin)
    call max_across_PEs(aMax)
    aMean = aMean / real(n)
  end subroutine subStats

end procedure chksum_h_2d
module procedure chksum_pair_B_2d
  logical :: sym
  logical :: vector_pair
  integer :: turns
  type(hor_index_type), pointer :: HI_in
  real, dimension(:,:), pointer :: arrayA_in, arrayB_in ! Rotated arrays [A ~> a]
  vector_pair = .true.
  if (present(scalar_pair)) vector_pair = .not. scalar_pair

  turns = HI%turns
  if (modulo(turns, 4) /= 0) then
    ! Rotate field back to the input grid
    allocate(HI_in)
    call rotate_hor_index(HI, -turns, HI_in)
    allocate(arrayA_in(HI_in%IsdB:HI_in%IedB, HI_in%JsdB:HI_in%JedB))
    allocate(arrayB_in(HI_in%IsdB:HI_in%IedB, HI_in%JsdB:HI_in%JedB))

    if (vector_pair) then
      call rotate_vector(arrayA, arrayB, -turns, arrayA_in, arrayB_in)
    else
      call rotate_array_pair(arrayA, arrayB, -turns, arrayA_in, arrayB_in)
    endif
  else
    HI_in => HI
    arrayA_in => arrayA
    arrayB_in => arrayB
  endif

  sym = .false. ; if (present(symmetric)) sym = symmetric

  if (present(haloshift)) then
    call chksum_B_2d(arrayA_in, 'x '//mesg, HI_in, haloshift, symmetric=sym, &
                     omit_corners=omit_corners, scale=scale, logunit=logunit, unscale=unscale)
    call chksum_B_2d(arrayB_in, 'y '//mesg, HI_in, haloshift, symmetric=sym, &
                     omit_corners=omit_corners, scale=scale, logunit=logunit, unscale=unscale)
  else
    call chksum_B_2d(arrayA_in, 'x '//mesg, HI_in, symmetric=sym, &
                     scale=scale, logunit=logunit, unscale=unscale)
    call chksum_B_2d(arrayB_in, 'y '//mesg, HI_in, symmetric=sym, &
                     scale=scale, logunit=logunit, unscale=unscale)
  endif

end procedure chksum_pair_B_2d
module procedure chksum_pair_B_3d
  logical :: vector_pair
  integer :: turns
  type(hor_index_type), pointer :: HI_in
  real, dimension(:,:,:), pointer :: arrayA_in, arrayB_in ! Rotated arrays [A ~> a]
  vector_pair = .true.
  if (present(scalar_pair)) vector_pair = .not. scalar_pair

  turns = HI%turns
  if (modulo(turns, 4) /= 0) then
    ! Rotate field back to the input grid
    allocate(HI_in)
    call rotate_hor_index(HI, -turns, HI_in)
    allocate(arrayA_in(HI_in%IsdB:HI_in%IedB, HI_in%JsdB:HI_in%JedB, size(arrayA, 3)))
    allocate(arrayB_in(HI_in%IsdB:HI_in%IedB, HI_in%JsdB:HI_in%JedB, size(arrayB, 3)))

    if (vector_pair) then
      call rotate_vector(arrayA, arrayB, -turns, arrayA_in, arrayB_in)
    else
      call rotate_array_pair(arrayA, arrayB, -turns, arrayA_in, arrayB_in)
    endif
  else
    HI_in => HI
    arrayA_in => arrayA
    arrayB_in => arrayB
  endif

  if (present(haloshift)) then
    call chksum_B_3d(arrayA_in, 'x '//mesg, HI_in, haloshift, symmetric, &
                     omit_corners, scale=scale, logunit=logunit, unscale=unscale)
    call chksum_B_3d(arrayB_in, 'y '//mesg, HI_in, haloshift, symmetric, &
                     omit_corners, scale=scale, logunit=logunit, unscale=unscale)
  else
    call chksum_B_3d(arrayA_in, 'x '//mesg, HI_in, symmetric=symmetric, &
                     scale=scale, logunit=logunit, unscale=unscale)
    call chksum_B_3d(arrayB_in, 'y '//mesg, HI_in, symmetric=symmetric, &
                     scale=scale, logunit=logunit, unscale=unscale)
  endif
end procedure chksum_pair_B_3d
module procedure chksum_B_2d
  real, pointer :: array(:,:)           ! Field array on the input grid [A ~> a]
  real, allocatable, dimension(:,:) :: rescaled_array ! The array with scaling undone [a]
  real, allocatable :: hash_array(:,:)  ! Subarray used to compute hash [a]
  type(hor_index_type), pointer :: HI   ! Horizontal index bounds of the input grid
  real :: scaling   ! Explicit rescaling factor [a A-1 ~> 1]
  integer :: iounit !< Log IO unit
  integer :: i, j, Is, Js
  real :: aMean, aMin, aMax  ! Array mean, global minimum and global maximum [a]
  integer :: bc0, bcSW, bcSE, bcNW, bcNE, hshift
  integer :: bcN, bcS, bcE, bcW
  logical :: do_corners, sym, sym_stats
  integer :: turns                      ! Quarter turns from input to model grid
  turns = HI_m%turns
  if (modulo(turns, 4) /= 0) then
    allocate(HI)
    call rotate_hor_index(HI_m, -turns, HI)
    allocate(array(HI%IsdB:HI%IedB, HI%JsdB:HI%JedB))
    call rotate_array(array_m, -turns, array)
  else
    HI => HI_m
    array => array_m
  endif

  if (checkForNaNs) then
    if (is_NaN(array(HI%IscB:HI%IecB,HI%JscB:HI%JecB))) &
      call chksum_error(FATAL, 'NaN detected: '//trim(mesg))
!   if (is_NaN(array)) &
!     call chksum_error(FATAL, 'NaN detected in halo: '//trim(mesg))
  endif

  scaling = 1.0
  if (present(unscale)) then ; scaling = unscale
  elseif (present(scale)) then ; scaling = scale ; endif

  iounit = error_unit ; if (present(logunit)) iounit = logunit
  sym_stats = .false. ; if (present(symmetric)) sym_stats = symmetric
  if (present(haloshift)) then ; if (haloshift > 0) sym_stats = .true. ; endif

  if (calculateStatistics) then
    if (present(unscale) .or. present(scale)) then
      allocate( rescaled_array(LBOUND(array,1):UBOUND(array,1), &
                               LBOUND(array,2):UBOUND(array,2)), source=0.0 )
      Is = HI%isc ; if (sym_stats) Is = HI%isc-1
      Js = HI%jsc ; if (sym_stats) Js = HI%jsc-1
      do J=Js,HI%JecB ; do I=Is,HI%IecB
        rescaled_array(I,J) = scaling*array(I,J)
      enddo ; enddo
      call subStats(HI, rescaled_array, sym_stats, aMean, aMin, aMax)
      deallocate(rescaled_array)
    else
      call subStats(HI, array, sym_stats, aMean, aMin, aMax)
    endif
    if (is_root_pe()) &
      call chk_sum_msg("B-point:", aMean, aMin, aMax, mesg, iounit)
  endif

  if (.not.writeChksums) return

  hshift = default_shift
  if (present(haloshift)) hshift = haloshift
  if (hshift<0) hshift = HI%ied-HI%iec

  if ( HI%iscB-hshift<HI%isdB .or. HI%iecB+hshift>HI%iedB .or. &
       HI%jscB-hshift<HI%jsdB .or. HI%jecB+hshift>HI%jedB ) then
    write(0,*) 'chksum_B_2d: haloshift =',hshift
    write(0,*) 'chksum_B_2d: isd,isc,iec,ied=',HI%isdB,HI%iscB,HI%iecB,HI%iedB
    write(0,*) 'chksum_B_2d: jsd,jsc,jec,jed=',HI%jsdB,HI%jscB,HI%jecB,HI%jedB
    call chksum_error(FATAL,'Error in chksum_B_2d '//trim(mesg))
  endif

  bc0 = subchk(array, HI, 0, 0, scaling)

  sym = .false. ; if (present(symmetric)) sym = symmetric

  if ((hshift==0) .and. .not.sym) then
    if (is_root_pe()) call chk_sum_msg("B-point:", bc0, mesg, iounit)
  else
    do_corners = .true.
    if (present(omit_corners)) do_corners = .not. omit_corners

    if (do_corners) then
      if (sym) then
        bcSW = subchk(array, HI, -hshift-1, -hshift-1, scaling)
        bcSE = subchk(array, HI, hshift, -hshift-1, scaling)
        bcNW = subchk(array, HI, -hshift-1, hshift, scaling)
      else
        bcSW = subchk(array, HI, -hshift, -hshift, scaling)
        bcSE = subchk(array, HI, hshift, -hshift, scaling)
        bcNW = subchk(array, HI, -hshift, hshift, scaling)
      endif
      bcNE = subchk(array, HI, hshift, hshift, scaling)

      if (is_root_pe()) &
        call chk_sum_msg("B-point:", bc0, bcSW, bcSE, bcNW, bcNE, mesg, iounit)
    else
      bcS = subchk(array, HI, 0, -hshift, scaling)
      bcE = subchk(array, HI, hshift, 0, scaling)
      bcW = subchk(array, HI, -hshift, 0, scaling)
      bcN = subchk(array, HI, 0, hshift, scaling)

      if (is_root_pe()) &
        call chk_sum_msg_NSEW("B-point:", bc0, bcN, bcS, bcE, bcW, mesg, iounit)
    endif
  endif

  if (writeHash .and. is_root_pe()) then
    allocate(hash_array(HI%isc:HI%iec, HI%jsc:HI%jec))
    hash_array(:,:) = scaling * array(HI%isc:HI%iec, HI%jsc:HI%jec)

    write(iounit, '("B-point: hash=", z8, 1x, a)') &
        murmur_hash(hash_array), mesg
    deallocate(hash_array)
  endif

  contains

  integer function subchk(array, HI, di, dj, unscale)
    type(hor_index_type), intent(in) ::  HI     !< A horizontal index type
    real, dimension(HI%IsdB:,HI%JsdB:), intent(in) :: array !< The array to be checksummed in
                                 !! arbitrary, possibly rescaled units [A ~> a]
    integer, intent(in) :: di    !< i- direction array shift for this checksum
    integer, intent(in) :: dj    !< j- direction array shift for this checksum
    real, intent(in)    :: unscale !< A factor to convert this array back to unscaled units
                                 !! for checksums and output [a A-1 ~> 1]
    integer :: i, j, bc
    subchk = 0
    ! This line deliberately uses the h-point computational domain.
    do J=HI%jsc+dj,HI%jec+dj ; do I=HI%isc+di,HI%iec+di
      bc = bitcount(abs(unscale*array(I,J)))
      subchk = subchk + bc
    enddo ; enddo
    call sum_across_PEs(subchk)
    subchk=mod(subchk, bc_modulus)
  end function subchk

  subroutine subStats(HI, array, sym_stats, aMean, aMin, aMax)
    type(hor_index_type), intent(in) ::  HI     !< A horizontal index type
    real, dimension(HI%IsdB:,HI%JsdB:), intent(in) :: array !< The array to be checksummed [a]
    logical,          intent(in) :: sym_stats !< If true, evaluate the statistics on the
                                              !! full symmetric computational domain.
    real, intent(out) :: aMean !< Array mean [a]
    real, intent(out) :: aMin !< Array minimum [a]
    real, intent(out) :: aMax !< Array maximum [a]

    integer :: i, j, n, IsB, JsB

    IsB = HI%isc ; if (sym_stats) IsB = HI%isc-1
    JsB = HI%jsc ; if (sym_stats) JsB = HI%jsc-1

    aMin = array(HI%isc,HI%jsc) ; aMax = aMin
    do J=JsB,HI%JecB ; do I=IsB,HI%IecB
      aMin = min(aMin, array(I,J))
      aMax = max(aMax, array(I,J))
    enddo ; enddo
    ! This line deliberately uses the h-point computational domain.
    aMean = reproducing_sum(array(HI%isc:HI%iec,HI%jsc:HI%jec))
    n = (1 + HI%jec - HI%jsc) * (1 + HI%iec - HI%isc)
    call sum_across_PEs(n)
    call min_across_PEs(aMin)
    call max_across_PEs(aMax)
    aMean = aMean / real(n)
  end subroutine subStats

end procedure chksum_B_2d
module procedure chksum_uv_2d
  logical :: vector_pair
  integer :: turns
  type(hor_index_type), pointer :: HI_in
  real, dimension(:,:), pointer :: arrayU_in, arrayV_in ! Rotated arrays [A ~> a]
  vector_pair = .true.
  if (present(scalar_pair)) vector_pair = .not. scalar_pair

  turns = HI%turns
  if (modulo(turns, 4) /= 0) then
    ! Rotate field back to the input grid
    allocate(HI_in)
    call rotate_hor_index(HI, -turns, HI_in)
    allocate(arrayU_in(HI_in%IsdB:HI_in%IedB, HI_in%jsd:HI_in%jed))
    allocate(arrayV_in(HI_in%isd:HI_in%ied, HI_in%JsdB:HI_in%JedB))

    if (vector_pair) then
      call rotate_vector(arrayU, arrayV, -turns, arrayU_in, arrayV_in)
    else
      call rotate_array_pair(arrayU, arrayV, -turns, arrayU_in, arrayV_in)
    endif
  else
    HI_in => HI
    arrayU_in => arrayU
    arrayV_in => arrayV
  endif

  if (present(haloshift)) then
    call chksum_u_2d(arrayU_in, 'u '//mesg, HI_in, haloshift, symmetric, &
                     omit_corners, scale=scale, logunit=logunit, unscale=unscale)
    call chksum_v_2d(arrayV_in, 'v '//mesg, HI_in, haloshift, symmetric, &
                     omit_corners, scale=scale, logunit=logunit, unscale=unscale)
  else
    call chksum_u_2d(arrayU_in, 'u '//mesg, HI_in, symmetric=symmetric, &
                     scale=scale, logunit=logunit, unscale=unscale)
    call chksum_v_2d(arrayV_in, 'v '//mesg, HI_in, symmetric=symmetric, &
                     scale=scale, logunit=logunit, unscale=unscale)
  endif
end procedure chksum_uv_2d
module procedure chksum_uv_3d
  logical :: vector_pair
  integer :: turns
  type(hor_index_type), pointer :: HI_in
  real, dimension(:,:,:), pointer :: arrayU_in, arrayV_in ! Rotated arrays [A ~> a]
  vector_pair = .true.
  if (present(scalar_pair)) vector_pair = .not. scalar_pair

  turns = HI%turns
  if (modulo(turns, 4) /= 0) then
    ! Rotate field back to the input grid
    allocate(HI_in)
    call rotate_hor_index(HI, -turns, HI_in)
    allocate(arrayU_in(HI_in%IsdB:HI_in%IedB, HI_in%jsd:HI_in%jed, size(arrayU, 3)))
    allocate(arrayV_in(HI_in%isd:HI_in%ied, HI_in%JsdB:HI_in%JedB, size(arrayV, 3)))

    if (vector_pair) then
      call rotate_vector(arrayU, arrayV, -turns, arrayU_in, arrayV_in)
    else
      call rotate_array_pair(arrayU, arrayV, -turns, arrayU_in, arrayV_in)
    endif
  else
    HI_in => HI
    arrayU_in => arrayU
    arrayV_in => arrayV
  endif

  if (present(haloshift)) then
    call chksum_u_3d(arrayU_in, 'u '//mesg, HI_in, haloshift, symmetric, &
                     omit_corners, scale=scale, logunit=logunit, unscale=unscale)
    call chksum_v_3d(arrayV_in, 'v '//mesg, HI_in, haloshift, symmetric, &
                     omit_corners, scale=scale, logunit=logunit, unscale=unscale)
  else
    call chksum_u_3d(arrayU_in, 'u '//mesg, HI_in, symmetric=symmetric, &
                     scale=scale, logunit=logunit, unscale=unscale)
    call chksum_v_3d(arrayV_in, 'v '//mesg, HI_in, symmetric=symmetric, &
                     scale=scale, logunit=logunit, unscale=unscale)
  endif
end procedure chksum_uv_3d
module procedure chksum_u_2d
  real, pointer :: array(:,:)           ! Field array on the input grid [A ~> a]
  real, allocatable, dimension(:,:) :: rescaled_array ! The array with scaling undone [a]
  real, allocatable :: hash_array(:,:)  ! Subarray used to compute hash [a]
  type(hor_index_type), pointer :: HI   ! Horizontal index bounds of the input grid
  real :: scaling   ! Explicit rescaling factor [a A-1 ~> 1]
  integer :: iounit !< Log IO unit
  integer :: i, j, Is
  real :: aMean, aMin, aMax  ! Array mean, global minimum and global maximum [a]
  integer :: bc0, bcSW, bcSE, bcNW, bcNE, hshift
  integer :: bcN, bcS, bcE, bcW
  logical :: do_corners, sym, sym_stats
  integer :: turns                      ! Quarter turns from input to model grid
  turns = HI_m%turns
  if (modulo(turns, 4) /= 0) then
    allocate(HI)
    call rotate_hor_index(HI_m, -turns, HI)
    if (modulo(turns, 2) /= 0) then
      ! Arrays originating from v-points must be handled by vchksum
      allocate(array(HI%isd:HI%ied, HI%JsdB:HI%JedB))
      call rotate_array(array_m, -turns, array)
      call vchksum(array, mesg, HI, haloshift, symmetric, omit_corners, &
                   scale=scale, logunit=logunit, unscale=unscale)
      return
    else
      allocate(array(HI%IsdB:HI%IedB, HI%jsd:HI%jed))
      call rotate_array(array_m, -turns, array)
    endif
  else
    HI => HI_m
    array => array_m
  endif

  if (checkForNaNs) then
    if (is_NaN(array(HI%IscB:HI%IecB,HI%jsc:HI%jec))) &
      call chksum_error(FATAL, 'NaN detected: '//trim(mesg))
!   if (is_NaN(array)) &
!     call chksum_error(FATAL, 'NaN detected in halo: '//trim(mesg))
  endif

  scaling = 1.0
  if (present(unscale)) then ; scaling = unscale
  elseif (present(scale)) then ; scaling = scale ; endif

  iounit = error_unit ; if (present(logunit)) iounit = logunit
  sym_stats = .false. ; if (present(symmetric)) sym_stats = symmetric
  if (present(haloshift)) then ; if (haloshift > 0) sym_stats = .true. ; endif

  if (calculateStatistics) then
    if (present(unscale) .or. present(scale)) then
      allocate( rescaled_array(LBOUND(array,1):UBOUND(array,1), &
                               LBOUND(array,2):UBOUND(array,2)), source=0.0 )
      Is = HI%isc ; if (sym_stats) Is = HI%isc-1
      do j=HI%jsc,HI%jec ; do I=Is,HI%IecB
        rescaled_array(I,j) = scaling*array(I,j)
      enddo ; enddo
      call subStats(HI, rescaled_array, sym_stats, aMean, aMin, aMax)
      deallocate(rescaled_array)
    else
      call subStats(HI, array, sym_stats, aMean, aMin, aMax)
    endif

    if (is_root_pe()) &
      call chk_sum_msg("u-point:", aMean, aMin, aMax, mesg, iounit)
  endif

  if (.not.writeChksums) return

  hshift = default_shift
  if (present(haloshift)) hshift = haloshift
  if (hshift<0) hshift = HI%iedB-HI%iecB

  if ( HI%iscB-hshift<HI%isdB .or. HI%iecB+hshift>HI%iedB .or. &
       HI%jsc-hshift<HI%jsd .or. HI%jec+hshift>HI%jed ) then
    write(0,*) 'chksum_u_2d: haloshift =',hshift
    write(0,*) 'chksum_u_2d: isd,isc,iec,ied=',HI%isdB,HI%iscB,HI%iecB,HI%iedB
    write(0,*) 'chksum_u_2d: jsd,jsc,jec,jed=',HI%jsd,HI%jsc,HI%jec,HI%jed
    call chksum_error(FATAL,'Error in chksum_u_2d '//trim(mesg))
  endif

  bc0 = subchk(array, HI, 0, 0, scaling)

  sym = .false. ; if (present(symmetric)) sym = symmetric

  if ((hshift==0) .and. .not.sym) then
    if (is_root_pe()) call chk_sum_msg("u-point:", bc0, mesg, iounit)
  else
    do_corners = .true.
    if (present(omit_corners)) do_corners = .not. omit_corners

    if (hshift==0) then
      bcW = subchk(array, HI, -hshift-1, 0, scaling)
      if (is_root_pe()) call chk_sum_msg_W("u-point:", bc0, bcW, mesg, iounit)
    elseif (do_corners) then
      if (sym) then
        bcSW = subchk(array, HI, -hshift-1, -hshift, scaling)
        bcNW = subchk(array, HI, -hshift-1, hshift, scaling)
      else
        bcSW = subchk(array, HI, -hshift, -hshift, scaling)
        bcNW = subchk(array, HI, -hshift, hshift, scaling)
      endif
      bcSE = subchk(array, HI, hshift, -hshift, scaling)
      bcNE = subchk(array, HI, hshift, hshift, scaling)

      if (is_root_pe()) &
        call chk_sum_msg("u-point:", bc0, bcSW, bcSE, bcNW, bcNE, mesg, iounit)
    else
      bcS = subchk(array, HI, 0, -hshift, scaling)
      bcE = subchk(array, HI, hshift, 0, scaling)
      if (sym) then
        bcW = subchk(array, HI, -hshift-1, 0, scaling)
      else
        bcW = subchk(array, HI, -hshift, 0, scaling)
      endif
      bcN = subchk(array, HI, 0, hshift, scaling)

      if (is_root_pe()) &
        call chk_sum_msg_NSEW("u-point:", bc0, bcN, bcS, bcE, bcW, mesg, iounit)
    endif
  endif

  if (writeHash .and. is_root_pe()) then
    allocate(hash_array(HI%isc:HI%iec, HI%jsc:HI%jec))
    hash_array(:,:) = scaling * array(HI%isc:HI%iec, HI%jsc:HI%jec)

    write(iounit, '("u-point: hash=", z8, 1x, a)') &
        murmur_hash(hash_array), mesg
    deallocate(hash_array)
  endif

  contains

  integer function subchk(array, HI, di, dj, unscale)
    type(hor_index_type), intent(in) ::  HI     !< A horizontal index type
    real, dimension(HI%IsdB:,HI%jsd:), intent(in) :: array !< The array to be checksummed in
                                 !! arbitrary, possibly rescaled units [A ~> a]
    integer, intent(in) :: di    !< i- direction array shift for this checksum
    integer, intent(in) :: dj    !< j- direction array shift for this checksum
    real, intent(in)    :: unscale !< A factor to convert this array back to unscaled units
                                 !! for checksums and output [a A-1 ~> 1]
    integer :: i, j, bc
    subchk = 0
    ! This line deliberately uses the h-point computational domain.
    do j=HI%jsc+dj,HI%jec+dj ; do I=HI%isc+di,HI%iec+di
      bc = bitcount(abs(unscale*array(I,j)))
      subchk = subchk + bc
    enddo ; enddo
    call sum_across_PEs(subchk)
    subchk=mod(subchk, bc_modulus)
  end function subchk

  subroutine subStats(HI, array, sym_stats, aMean, aMin, aMax)
    type(hor_index_type), intent(in) ::  HI     !< A horizontal index type
    real, dimension(HI%IsdB:,HI%jsd:), intent(in) :: array !< The array to be checksummed [a]
    logical,          intent(in) :: sym_stats !< If true, evaluate the statistics on the
                                              !! full symmetric computational domain.
    real, intent(out) :: aMean !< Array mean [a]
    real, intent(out) :: aMin !< Array minimum [a]
    real, intent(out) :: aMax !< Array maximum [a]

    integer :: i, j, n, IsB

    IsB = HI%isc ; if (sym_stats) IsB = HI%isc-1

    aMin = array(HI%isc,HI%jsc) ; aMax = aMin
    do j=HI%jsc,HI%jec ; do I=IsB,HI%IecB
      aMin = min(aMin, array(I,j))
      aMax = max(aMax, array(I,j))
    enddo ; enddo
    ! This line deliberately uses the h-point computational domain.
    aMean = reproducing_sum(array(HI%isc:HI%iec,HI%jsc:HI%jec))
    n = (1 + HI%jec - HI%jsc) * (1 + HI%iec - HI%isc)
    call sum_across_PEs(n)
    call min_across_PEs(aMin)
    call max_across_PEs(aMax)
    aMean = aMean / real(n)
  end subroutine subStats

end procedure chksum_u_2d
module procedure chksum_v_2d
  real, pointer :: array(:,:)           ! Field array on the input grid [A ~> a]
  real, allocatable, dimension(:,:) :: rescaled_array ! The array with scaling undone [a]
  real, allocatable :: hash_array(:,:)  ! Subarray used to compute hash [a]
  type(hor_index_type), pointer :: HI   ! Horizontal index bounds of the input grid
  real :: scaling   ! Explicit rescaling factor [a A-1 ~> 1]
  integer :: iounit !< Log IO unit
  integer :: i, j, Js
  real :: aMean, aMin, aMax  ! Array mean, global minimum and global maximum [a]
  integer :: bc0, bcSW, bcSE, bcNW, bcNE, hshift
  integer :: bcN, bcS, bcE, bcW
  logical :: do_corners, sym, sym_stats
  integer :: turns                      ! Quarter turns from input to model grid
  turns = HI_m%turns
  if (modulo(turns, 4) /= 0) then
    allocate(HI)
    call rotate_hor_index(HI_m, -turns, HI)
    if (modulo(turns, 2) /= 0) then
      ! Arrays originating from u-points must be handled by uchksum
      allocate(array(HI%IsdB:HI%IedB, HI%jsd:HI%jed))
      call rotate_array(array_m, -turns, array)
      call uchksum(array, mesg, HI, haloshift, symmetric, omit_corners, &
                   scale=scale, logunit=logunit, unscale=unscale)
      return
    else
      allocate(array(HI%isd:HI%ied, HI%JsdB:HI%JedB))
      call rotate_array(array_m, -turns, array)
    endif
  else
    HI => HI_m
    array => array_m
  endif

  if (checkForNaNs) then
    if (is_NaN(array(HI%isc:HI%iec,HI%JscB:HI%JecB))) &
      call chksum_error(FATAL, 'NaN detected: '//trim(mesg))
!   if (is_NaN(array)) &
!     call chksum_error(FATAL, 'NaN detected in halo: '//trim(mesg))
  endif

  scaling = 1.0
  if (present(unscale)) then ; scaling = unscale
  elseif (present(scale)) then ; scaling = scale ; endif

  iounit = error_unit ; if (present(logunit)) iounit = logunit
  sym_stats = .false. ; if (present(symmetric)) sym_stats = symmetric
  if (present(haloshift)) then ; if (haloshift > 0) sym_stats = .true. ; endif

  if (calculateStatistics) then
    if (present(unscale) .or. present(scale)) then
      allocate( rescaled_array(LBOUND(array,1):UBOUND(array,1), &
                               LBOUND(array,2):UBOUND(array,2)), source=0.0 )
      Js = HI%jsc ; if (sym_stats) Js = HI%jsc-1
      do J=Js,HI%JecB ; do i=HI%isc,HI%iec
        rescaled_array(i,J) = scaling*array(i,J)
      enddo ; enddo
      call subStats(HI, rescaled_array, sym_stats, aMean, aMin, aMax)
      deallocate(rescaled_array)
    else
      call subStats(HI, array, sym_stats, aMean, aMin, aMax)
    endif

    if (is_root_pe()) &
      call chk_sum_msg("v-point:", aMean, aMin, aMax, mesg, iounit)
  endif

  if (.not.writeChksums) return

  hshift = default_shift
  if (present(haloshift)) hshift = haloshift
  if (hshift<0) hshift = HI%ied-HI%iec

  if ( HI%isc-hshift<HI%isd .or. HI%iec+hshift>HI%ied .or. &
       HI%jscB-hshift<HI%jsdB .or. HI%jecB+hshift>HI%jedB ) then
    write(0,*) 'chksum_v_2d: haloshift =',hshift
    write(0,*) 'chksum_v_2d: isd,isc,iec,ied=',HI%isd,HI%isc,HI%iec,HI%ied
    write(0,*) 'chksum_v_2d: jsd,jsc,jec,jed=',HI%jsdB,HI%jscB,HI%jecB,HI%jedB
    call chksum_error(FATAL,'Error in chksum_v_2d '//trim(mesg))
  endif

  bc0 = subchk(array, HI, 0, 0, scaling)

  sym = .false. ; if (present(symmetric)) sym = symmetric

  if ((hshift==0) .and. .not.sym) then
    if (is_root_pe()) call chk_sum_msg("v-point:", bc0, mesg, iounit)
  else
    do_corners = .true.
    if (present(omit_corners)) do_corners = .not. omit_corners

    if (hshift==0) then
      bcS = subchk(array, HI, 0, -hshift-1, scaling)
      if (is_root_pe()) call chk_sum_msg_S("v-point:", bc0, bcS, mesg, iounit)
    elseif (do_corners) then
      if (sym) then
        bcSW = subchk(array, HI, -hshift, -hshift-1, scaling)
        bcSE = subchk(array, HI, hshift, -hshift-1, scaling)
      else
        bcSW = subchk(array, HI, -hshift, -hshift, scaling)
        bcSE = subchk(array, HI, hshift, -hshift, scaling)
      endif
      bcNW = subchk(array, HI, -hshift, hshift, scaling)
      bcNE = subchk(array, HI, hshift, hshift, scaling)

      if (is_root_pe()) &
        call chk_sum_msg("v-point:", bc0, bcSW, bcSE, bcNW, bcNE, mesg, iounit)
    else
      if (sym) then
        bcS = subchk(array, HI, 0, -hshift-1, scaling)
      else
        bcS = subchk(array, HI, 0, -hshift, scaling)
      endif
      bcE = subchk(array, HI, hshift, 0, scaling)
      bcW = subchk(array, HI, -hshift, 0, scaling)
      bcN = subchk(array, HI, 0, hshift, scaling)

      if (is_root_pe()) &
        call chk_sum_msg_NSEW("v-point:", bc0, bcN, bcS, bcE, bcW, mesg, iounit)
    endif
  endif

  if (writeHash .and. is_root_pe()) then
    allocate(hash_array(HI%isc:HI%iec, HI%jsc:HI%jec))
    hash_array(:,:) = scaling * array(HI%isc:HI%iec, HI%jsc:HI%jec)

    write(iounit, '("v-point: hash=", z8, 1x, a)') &
        murmur_hash(hash_array), mesg
    deallocate(hash_array)
  endif

  contains

  integer function subchk(array, HI, di, dj, unscale)
    type(hor_index_type), intent(in) ::  HI     !< A horizontal index type
    real, dimension(HI%isd:,HI%JsdB:), intent(in) :: array !< The array to be checksummed in
                                 !! arbitrary, possibly rescaled units [A ~> a]
    integer, intent(in) :: di    !< i- direction array shift for this checksum
    integer, intent(in) :: dj    !< j- direction array shift for this checksum
    real, intent(in)    :: unscale !< A factor to convert this array back to unscaled units
                                 !! for checksums and output [a A-1 ~> 1]
    integer :: i, j, bc
    subchk = 0
    ! This line deliberately uses the h-point computational domain.
    do J=HI%jsc+dj,HI%jec+dj ; do i=HI%isc+di,HI%iec+di
      bc = bitcount(abs(unscale*array(i,J)))
      subchk = subchk + bc
    enddo ; enddo
    call sum_across_PEs(subchk)
    subchk=mod(subchk, bc_modulus)
  end function subchk

  subroutine subStats(HI, array, sym_stats, aMean, aMin, aMax)
    type(hor_index_type), intent(in) ::  HI     !< A horizontal index type
    real, dimension(HI%isd:,HI%JsdB:), intent(in) :: array !< The array to be checksummed [a]
    logical,          intent(in) :: sym_stats !< If true, evaluate the statistics on the
                                              !! full symmetric computational domain.
    real, intent(out) :: aMean !< Array mean [a]
    real, intent(out) :: aMin !< Array minimum [a]
    real, intent(out) :: aMax !< Array maximum [a]

    integer :: i, j, n, JsB

    JsB = HI%jsc ; if (sym_stats) JsB = HI%jsc-1

    aMin = array(HI%isc,HI%jsc) ; aMax = aMin
    do J=JsB,HI%JecB ; do i=HI%isc,HI%iec
      aMin = min(aMin, array(i,J))
      aMax = max(aMax, array(i,J))
    enddo ; enddo
    ! This line deliberately uses the h-computational domain.
    aMean = reproducing_sum(array(HI%isc:HI%iec,HI%jsc:HI%jec))
    n = (1 + HI%jec - HI%jsc) * (1 + HI%iec - HI%isc)
    call sum_across_PEs(n)
    call min_across_PEs(aMin)
    call max_across_PEs(aMax)
    aMean = aMean / real(n)
  end subroutine subStats

end procedure chksum_v_2d
module procedure chksum_h_3d
  real, pointer :: array(:,:,:)         ! Field array on the input grid [A ~> a]
  real, allocatable, dimension(:,:,:) :: rescaled_array ! The array with scaling undone [a]
  real, allocatable :: hash_array(:,:,:)  ! Subarray used to compute hash [a]
  type(hor_index_type), pointer :: HI   ! Horizontal index bounds of the input grid
  real :: scaling   ! Explicit rescaling factor [a A-1 ~> 1]
  integer :: iounit !< Log IO unit
  integer :: i, j, k
  real :: aMean, aMin, aMax  ! Array mean, global minimum and global maximum [a]
  integer :: bc0, bcSW, bcSE, bcNW, bcNE, hshift
  integer :: bcN, bcS, bcE, bcW
  logical :: do_corners
  integer :: turns                      ! Quarter turns from input to model grid
  turns = HI_m%turns
  if (modulo(turns, 4) /= 0) then
    allocate(HI)
    call rotate_hor_index(HI_m, -turns, HI)
    allocate(array(HI%isd:HI%ied, HI%jsd:HI%jed, size(array_m, 3)))
    call rotate_array(array_m, -turns, array)
  else
    HI => HI_m
    array => array_m
  endif

  if (checkForNaNs) then
    if (is_NaN(array(HI%isc:HI%iec,HI%jsc:HI%jec,:))) &
      call chksum_error(FATAL, 'NaN detected: '//trim(mesg))
!   if (is_NaN(array)) &
!     call chksum_error(FATAL, 'NaN detected in halo: '//trim(mesg))
  endif

  scaling = 1.0
  if (present(unscale)) then ; scaling = unscale
  elseif (present(scale)) then ; scaling = scale ; endif

  iounit = error_unit ; if (present(logunit)) iounit = logunit

  if (calculateStatistics) then
    if (present(unscale) .or. present(scale)) then
      allocate( rescaled_array(LBOUND(array,1):UBOUND(array,1), &
                               LBOUND(array,2):UBOUND(array,2), &
                               LBOUND(array,3):UBOUND(array,3)), source=0.0 )
      do k=1,size(array,3) ; do j=HI%jsc,HI%jec ; do i=HI%isc,HI%iec
        rescaled_array(i,j,k) = scaling*array(i,j,k)
      enddo ; enddo ; enddo

      call subStats(HI, rescaled_array, aMean, aMin, aMax)
      deallocate(rescaled_array)
    else
      call subStats(HI, array, aMean, aMin, aMax)
    endif

    if (is_root_pe()) &
      call chk_sum_msg("h-point:", aMean, aMin, aMax, mesg, iounit)
  endif

  if (.not.writeChksums) return

  hshift = default_shift
  if (present(haloshift)) hshift = haloshift
  if (hshift<0) hshift = HI%ied-HI%iec

  if ( HI%isc-hshift<HI%isd .or. HI%iec+hshift>HI%ied .or. &
       HI%jsc-hshift<HI%jsd .or. HI%jec+hshift>HI%jed ) then
    write(0,*) 'chksum_h_3d: haloshift =',hshift
    write(0,*) 'chksum_h_3d: isd,isc,iec,ied=',HI%isd,HI%isc,HI%iec,HI%ied
    write(0,*) 'chksum_h_3d: jsd,jsc,jec,jed=',HI%jsd,HI%jsc,HI%jec,HI%jed
    call chksum_error(FATAL,'Error in chksum_h_3d '//trim(mesg))
  endif

  bc0 = subchk(array, HI, 0, 0, scaling)

  if (hshift==0) then
    if (is_root_pe()) call chk_sum_msg("h-point:", bc0, mesg, iounit)
  else
    do_corners = .true.
    if (present(omit_corners)) do_corners = .not. omit_corners

    if (do_corners) then
      bcSW = subchk(array, HI, -hshift, -hshift, scaling)
      bcSE = subchk(array, HI, hshift, -hshift, scaling)
      bcNW = subchk(array, HI, -hshift, hshift, scaling)
      bcNE = subchk(array, HI, hshift, hshift, scaling)

      if (is_root_pe()) &
        call chk_sum_msg("h-point:", bc0, bcSW, bcSE, bcNW, bcNE, mesg, iounit)
    else
      bcS = subchk(array, HI, 0, -hshift, scaling)
      bcE = subchk(array, HI, hshift, 0, scaling)
      bcW = subchk(array, HI, -hshift, 0, scaling)
      bcN = subchk(array, HI, 0, hshift, scaling)

      if (is_root_pe()) &
        call chk_sum_msg_NSEW("h-point:", bc0, bcN, bcS, bcE, bcW, mesg, iounit)
    endif
  endif

  if (writeHash .and. is_root_pe()) then
    allocate(hash_array(HI%isc:HI%iec, HI%jsc:HI%jec, size(array, 3)))
    hash_array(:,:,:) = scaling * array(HI%isc:HI%iec, HI%jsc:HI%jec, :)

    write(iounit, '("h-point: hash=", z8, 1x, a)') &
        murmur_hash(hash_array), mesg
    deallocate(hash_array)
  endif

  contains

  integer function subchk(array, HI, di, dj, unscale)
    type(hor_index_type), intent(in) ::  HI     !< A horizontal index type
    real, dimension(HI%isd:,HI%jsd:,:), intent(in) :: array !< The array to be checksummed in
                                 !! arbitrary, possibly rescaled units [A ~> a]
    integer, intent(in) :: di    !< i- direction array shift for this checksum
    integer, intent(in) :: dj    !< j- direction array shift for this checksum
    real, intent(in)    :: unscale !< A factor to convert this array back to unscaled units
                                 !! for checksums and output [a A-1 ~> 1]
    integer :: i, j, k, bc
    subchk = 0
    do k=LBOUND(array,3),UBOUND(array,3) ; do j=HI%jsc+dj,HI%jec+dj ; do i=HI%isc+di,HI%iec+di
      bc = bitcount(abs(unscale*array(i,j,k)))
      subchk = subchk + bc
    enddo ; enddo ; enddo
    call sum_across_PEs(subchk)
    subchk=mod(subchk, bc_modulus)
  end function subchk

  subroutine subStats(HI, array, aMean, aMin, aMax)
    type(hor_index_type), intent(in) ::  HI     !< A horizontal index type
    real, dimension(HI%isd:,HI%jsd:,:), intent(in) :: array !< The array to be checksummed [a]
    real, intent(out) :: aMean !<  Array mean [a]
    real, intent(out) :: aMin !< Array minimum [a]
    real, intent(out) :: aMax !< Array maximum [a]

    integer :: i, j, k, n

    aMin = array(HI%isc,HI%jsc,1)
    aMax = array(HI%isc,HI%jsc,1)
    n = 0
    do k=LBOUND(array,3),UBOUND(array,3) ; do j=HI%jsc,HI%jec ; do i=HI%isc,HI%iec
      aMin = min(aMin, array(i,j,k))
      aMax = max(aMax, array(i,j,k))
      n = n + 1
    enddo ; enddo ; enddo
    aMean = reproducing_sum(array(HI%isc:HI%iec,HI%jsc:HI%jec,:))
    call sum_across_PEs(n)
    call min_across_PEs(aMin)
    call max_across_PEs(aMax)
    aMean = aMean / real(n)
  end subroutine subStats

end procedure chksum_h_3d
module procedure chksum_B_3d
  real, pointer :: array(:,:,:)         ! Field array on the input grid [A ~> a]
  real, allocatable, dimension(:,:,:) :: rescaled_array ! The array with scaling undone [a]
  real, allocatable :: hash_array(:,:,:)  ! Subarray used to compute hash [a]
  type(hor_index_type), pointer :: HI   ! Horizontal index bounds of the input grid
  real :: scaling   ! Explicit rescaling factor [a A-1 ~> 1]
  integer :: iounit !< Log IO unit
  integer :: i, j, k, Is, Js
  real :: aMean, aMin, aMax  ! Array mean, global minimum and global maximum [a]
  integer :: bc0, bcSW, bcSE, bcNW, bcNE, hshift
  integer :: bcN, bcS, bcE, bcW
  logical :: do_corners, sym, sym_stats
  integer :: turns                      ! Quarter turns from input to model grid
  turns = HI_m%turns
  if (modulo(turns, 4) /= 0) then
    allocate(HI)
    call rotate_hor_index(HI_m, -turns, HI)
    allocate(array(HI%IsdB:HI%IedB, HI%JsdB:HI%JedB, size(array_m, 3)))
    call rotate_array(array_m, -turns, array)
  else
    HI => HI_m
    array => array_m
  endif

  if (checkForNaNs) then
    if (is_NaN(array(HI%IscB:HI%IecB,HI%JscB:HI%JecB,:))) &
      call chksum_error(FATAL, 'NaN detected: '//trim(mesg))
!   if (is_NaN(array)) &
!     call chksum_error(FATAL, 'NaN detected in halo: '//trim(mesg))
  endif

  scaling = 1.0
  if (present(unscale)) then ; scaling = unscale
  elseif (present(scale)) then ; scaling = scale ; endif

  iounit = error_unit ; if (present(logunit)) iounit = logunit
  sym_stats = .false. ; if (present(symmetric)) sym_stats = symmetric
  if (present(haloshift)) then ; if (haloshift > 0) sym_stats = .true. ; endif

  if (calculateStatistics) then
    if (present(unscale) .or. present(scale)) then
      allocate( rescaled_array(LBOUND(array,1):UBOUND(array,1), &
                               LBOUND(array,2):UBOUND(array,2), &
                               LBOUND(array,3):UBOUND(array,3)), source=0.0 )
      Is = HI%isc ; if (sym_stats) Is = HI%isc-1
      Js = HI%jsc ; if (sym_stats) Js = HI%jsc-1
      do k=1,size(array,3) ; do J=Js,HI%JecB ; do I=Is,HI%IecB
        rescaled_array(I,J,k) = scaling*array(I,J,k)
      enddo ; enddo ; enddo
      call subStats(HI, rescaled_array, sym_stats, aMean, aMin, aMax)
      deallocate(rescaled_array)
    else
      call subStats(HI, array, sym_stats, aMean, aMin, aMax)
    endif

    if (is_root_pe()) &
      call chk_sum_msg("B-point:", aMean, aMin, aMax, mesg, iounit)
  endif

  if (.not.writeChksums) return

  hshift = default_shift
  if (present(haloshift)) hshift = haloshift
  if (hshift<0) hshift = HI%ied-HI%iec

  if ( HI%isc-hshift<HI%isd .or. HI%iec+hshift>HI%ied .or. &
       HI%jsc-hshift<HI%jsd .or. HI%jec+hshift>HI%jed ) then
    write(0,*) 'chksum_B_3d: haloshift =',hshift
    write(0,*) 'chksum_B_3d: isd,isc,iec,ied=',HI%isd,HI%isc,HI%iec,HI%ied
    write(0,*) 'chksum_B_3d: jsd,jsc,jec,jed=',HI%jsd,HI%jsc,HI%jec,HI%jed
    call chksum_error(FATAL,'Error in chksum_B_3d '//trim(mesg))
  endif

  bc0 = subchk(array, HI, 0, 0, scaling)

  sym = .false. ; if (present(symmetric)) sym = symmetric

  if ((hshift==0) .and. .not.sym) then
    if (is_root_pe()) call chk_sum_msg("B-point:", bc0, mesg, iounit)
  else
  do_corners = .true.
    if (present(omit_corners)) do_corners = .not. omit_corners

    if (do_corners) then
      if (sym) then
        bcSW = subchk(array, HI, -hshift-1, -hshift-1, scaling)
        bcSE = subchk(array, HI, hshift, -hshift-1, scaling)
        bcNW = subchk(array, HI, -hshift-1, hshift, scaling)
      else
        bcSW = subchk(array, HI, -hshift-1, -hshift-1, scaling)
        bcSE = subchk(array, HI, hshift, -hshift-1, scaling)
        bcNW = subchk(array, HI, -hshift-1, hshift, scaling)
      endif
      bcNE = subchk(array, HI, hshift, hshift, scaling)

      if (is_root_pe()) &
        call chk_sum_msg("B-point:", bc0, bcSW, bcSE, bcNW, bcNE, mesg, iounit)
    else
      if (sym) then
        bcS = subchk(array, HI, 0, -hshift-1, scaling)
        bcW = subchk(array, HI, -hshift-1, 0, scaling)
      else
        bcS = subchk(array, HI, 0, -hshift, scaling)
        bcW = subchk(array, HI, -hshift, 0, scaling)
      endif
      bcE = subchk(array, HI, hshift, 0, scaling)
      bcN = subchk(array, HI, 0, hshift, scaling)

      if (is_root_pe()) &
        call chk_sum_msg_NSEW("B-point:", bc0, bcN, bcS, bcE, bcW, mesg, iounit)
    endif
  endif

  if (writeHash .and. is_root_pe()) then
    allocate(hash_array(HI%isc:HI%iec, HI%jsc:HI%jec, size(array, 3)))
    hash_array(:,:,:) = scaling * array(HI%isc:HI%iec, HI%jsc:HI%jec, :)

    write(iounit, '("B-point: hash=", z8, 1x, a)') &
        murmur_hash(hash_array), mesg
    deallocate(hash_array)
  endif

  contains

  integer function subchk(array, HI, di, dj, unscale)
    type(hor_index_type), intent(in) ::  HI     !< A horizontal index type
    real, dimension(HI%IsdB:,HI%JsdB:,:), intent(in) :: array !< The array to be checksummed in
                                 !! arbitrary, possibly rescaled units [A ~> a]
    integer, intent(in) :: di    !< i- direction array shift for this checksum
    integer, intent(in) :: dj    !< j- direction array shift for this checksum
    real, intent(in)    :: unscale !< A factor to convert this array back to unscaled units
                                 !! for checksums and output [a A-1 ~> 1]
    integer :: i, j, k, bc
    subchk = 0
    ! This line deliberately uses the h-point computational domain.
    do k=LBOUND(array,3),UBOUND(array,3) ; do J=HI%jsc+dj,HI%jec+dj ; do I=HI%isc+di,HI%iec+di
      bc = bitcount(abs(unscale*array(I,J,k)))
      subchk = subchk + bc
    enddo ; enddo ; enddo
    call sum_across_PEs(subchk)
    subchk=mod(subchk, bc_modulus)
  end function subchk

  subroutine subStats(HI, array, sym_stats, aMean, aMin, aMax)
    type(hor_index_type), intent(in) ::  HI     !< A horizontal index type
    real, dimension(HI%IsdB:,HI%JsdB:,:), intent(in) :: array !< The array to be checksummed [a]
    logical,          intent(in) :: sym_stats !< If true, evaluate the statistics on the
                                              !! full symmetric computational domain.
    real, intent(out) :: aMean !< Array mean [a]
    real, intent(out) :: aMin !< Array minimum [a]
    real, intent(out) :: aMax !< Array maximum [a]

    integer :: i, j, k, n, IsB, JsB

    IsB = HI%isc ; if (sym_stats) IsB = HI%isc-1
    JsB = HI%jsc ; if (sym_stats) JsB = HI%jsc-1

    aMin = array(HI%isc,HI%jsc,1) ; aMax = aMin
    do k=LBOUND(array,3),UBOUND(array,3) ; do J=JsB,HI%JecB ; do I=IsB,HI%IecB
      aMin = min(aMin, array(I,J,k))
      aMax = max(aMax, array(I,J,k))
    enddo ; enddo ; enddo
    aMean = reproducing_sum(array(HI%isc:HI%iec,HI%jsc:HI%jec,:))
    n = (1 + HI%jec - HI%jsc) * (1 + HI%iec - HI%isc) * size(array,3)
    call sum_across_PEs(n)
    call min_across_PEs(aMin)
    call max_across_PEs(aMax)
    aMean = aMean / real(n)
  end subroutine subStats

end procedure chksum_B_3d
module procedure chksum_u_3d
  real, pointer :: array(:,:,:)         ! Field array on the input grid [A ~> a]
  real, allocatable, dimension(:,:,:) :: rescaled_array ! The array with scaling undone [a]
  real, allocatable :: hash_array(:,:,:)  ! Subarray used to compute hash [a]
  type(hor_index_type), pointer :: HI   ! Horizontal index bounds of the input grid
  real :: scaling   ! Explicit rescaling factor [a A-1 ~> 1]
  integer :: iounit !< Log IO unit
  integer :: i, j, k, Is
  real :: aMean, aMin, aMax  ! Array mean, global minimum and global maximum [a]
  integer :: bc0, bcSW, bcSE, bcNW, bcNE, hshift
  integer :: bcN, bcS, bcE, bcW
  logical :: do_corners, sym, sym_stats
  integer :: turns                      ! Quarter turns from input to model grid
  turns = HI_m%turns
  if (modulo(turns, 4) /= 0) then
    allocate(HI)
    call rotate_hor_index(HI_m, -turns, HI)
    if (modulo(turns, 2) /= 0) then
      ! Arrays originating from v-points must be handled by vchksum
      allocate(array(HI%isd:HI%ied, HI%JsdB:HI%JedB, size(array_m, 3)))
      call rotate_array(array_m, -turns, array)
      call vchksum(array, mesg, HI, haloshift, symmetric, omit_corners, &
                   scale=scale, logunit=logunit, unscale=unscale)
      return
    else
      allocate(array(HI%IsdB:HI%IedB, HI%jsd:HI%jed, size(array_m, 3)))
      call rotate_array(array_m, -turns, array)
    endif
  else
    HI => HI_m
    array => array_m
  endif

  if (checkForNaNs) then
    if (is_NaN(array(HI%IscB:HI%IecB,HI%jsc:HI%jec,:))) &
      call chksum_error(FATAL, 'NaN detected: '//trim(mesg))
!   if (is_NaN(array)) &
!     call chksum_error(FATAL, 'NaN detected in halo: '//trim(mesg))
  endif

  scaling = 1.0
  if (present(unscale)) then ; scaling = unscale
  elseif (present(scale)) then ; scaling = scale ; endif

  iounit = error_unit ; if (present(logunit)) iounit = logunit
  sym_stats = .false. ; if (present(symmetric)) sym_stats = symmetric
  if (present(haloshift)) then ; if (haloshift > 0) sym_stats = .true. ; endif

  if (calculateStatistics) then
    if (present(unscale) .or. present(scale)) then
      allocate( rescaled_array(LBOUND(array,1):UBOUND(array,1), &
                               LBOUND(array,2):UBOUND(array,2), &
                               LBOUND(array,3):UBOUND(array,3)), source=0.0 )
      Is = HI%isc ; if (sym_stats) Is = HI%isc-1
      do k=1,size(array,3) ; do j=HI%jsc,HI%jec ; do I=Is,HI%IecB
        rescaled_array(I,j,k) = scaling*array(I,j,k)
      enddo ; enddo ; enddo
      call subStats(HI, rescaled_array, sym_stats, aMean, aMin, aMax)
      deallocate(rescaled_array)
    else
      call subStats(HI, array, sym_stats, aMean, aMin, aMax)
    endif
    if (is_root_pe()) &
      call chk_sum_msg("u-point:", aMean, aMin, aMax, mesg, iounit)
  endif

  if (.not.writeChksums) return

  hshift = default_shift
  if (present(haloshift)) hshift = haloshift
  if (hshift<0) hshift = HI%ied-HI%iec

  if ( HI%isc-hshift<HI%isd .or. HI%iec+hshift>HI%ied .or. &
       HI%jsc-hshift<HI%jsd .or. HI%jec+hshift>HI%jed ) then
    write(0,*) 'chksum_u_3d: haloshift =',hshift
    write(0,*) 'chksum_u_3d: isd,isc,iec,ied=',HI%isd,HI%isc,HI%iec,HI%ied
    write(0,*) 'chksum_u_3d: jsd,jsc,jec,jed=',HI%jsd,HI%jsc,HI%jec,HI%jed
    call chksum_error(FATAL,'Error in chksum_u_3d '//trim(mesg))
  endif

  bc0 = subchk(array, HI, 0, 0, scaling)

  sym = .false. ; if (present(symmetric)) sym = symmetric

  if ((hshift==0) .and. .not.sym) then
    if (is_root_pe()) call chk_sum_msg("u-point:", bc0, mesg, iounit)
  else
    do_corners = .true.
    if (present(omit_corners)) do_corners = .not. omit_corners

    if (hshift==0) then
      bcW = subchk(array, HI, -hshift-1, 0, scaling)
      if (is_root_pe()) call chk_sum_msg_W("u-point:", bc0, bcW, mesg, iounit)
    elseif (do_corners) then
      if (sym) then
        bcSW = subchk(array, HI, -hshift-1, -hshift, scaling)
        bcNW = subchk(array, HI, -hshift-1, hshift, scaling)
      else
        bcSW = subchk(array, HI, -hshift, -hshift, scaling)
        bcNW = subchk(array, HI, -hshift, hshift, scaling)
      endif
      bcSE = subchk(array, HI, hshift, -hshift, scaling)
      bcNE = subchk(array, HI, hshift, hshift, scaling)

      if (is_root_pe()) &
        call chk_sum_msg("u-point:", bc0, bcSW, bcSE, bcNW, bcNE, mesg, iounit)
    else
      bcS = subchk(array, HI, 0, -hshift, scaling)
      bcE = subchk(array, HI, hshift, 0, scaling)
      if (sym) then
        bcW = subchk(array, HI, -hshift-1, 0, scaling)
      else
        bcW = subchk(array, HI, -hshift, 0, scaling)
      endif
      bcN = subchk(array, HI, 0, hshift, scaling)

      if (is_root_pe()) &
        call chk_sum_msg_NSEW("u-point:", bc0, bcN, bcS, bcE, bcW, mesg, iounit)
    endif
  endif

  if (writeHash .and. is_root_pe()) then
    allocate(hash_array(HI%isc:HI%iec, HI%jsc:HI%jec, size(array, 3)))
    hash_array(:,:,:) = scaling * array(HI%isc:HI%iec, HI%jsc:HI%jec, :)

    write(iounit, '("u-point: hash=", z8, 1x, a)') &
        murmur_hash(hash_array), mesg
    deallocate(hash_array)
  endif

  contains

  integer function subchk(array, HI, di, dj, unscale)
    type(hor_index_type), intent(in) ::  HI     !< A horizontal index type
    real, dimension(HI%IsdB:,HI%jsd:,:), intent(in) :: array !< The array to be checksummed in
                                 !! arbitrary, possibly rescaled units [A ~> a]
    integer, intent(in) :: di    !< i- direction array shift for this checksum
    integer, intent(in) :: dj    !< j- direction array shift for this checksum
    real, intent(in)    :: unscale !< A factor to convert this array back to unscaled units
                                 !! for checksums and output [a A-1 ~> 1]
    integer :: i, j, k, bc
    subchk = 0
    ! This line deliberately uses the h-point computational domain.
    do k=LBOUND(array,3),UBOUND(array,3) ; do j=HI%jsc+dj,HI%jec+dj ; do I=HI%isc+di,HI%iec+di
      bc = bitcount(abs(unscale*array(I,j,k)))
      subchk = subchk + bc
    enddo ; enddo ; enddo
    call sum_across_PEs(subchk)
    subchk=mod(subchk, bc_modulus)
  end function subchk

  subroutine subStats(HI, array, sym_stats, aMean, aMin, aMax)
    type(hor_index_type), intent(in) ::  HI     !< A horizontal index type
    real, dimension(HI%IsdB:,HI%jsd:,:), intent(in) :: array !< The array to be checksummed [a]
    logical,          intent(in) :: sym_stats !< If true, evaluate the statistics on the
                                              !! full symmetric computational domain.
    real, intent(out) :: aMean !< Array mean [a]
    real, intent(out) :: aMin !< Array minimum [a]
    real, intent(out) :: aMax !< Array maximum [a]

    integer :: i, j, k, n, IsB

    IsB = HI%isc ; if (sym_stats) IsB = HI%isc-1

    aMin = array(HI%isc,HI%jsc,1) ; aMax = aMin
    do k=LBOUND(array,3),UBOUND(array,3) ; do j=HI%jsc,HI%jec ; do I=IsB,HI%IecB
      aMin = min(aMin, array(I,j,k))
      aMax = max(aMax, array(I,j,k))
    enddo ; enddo ; enddo
    ! This line deliberately uses the h-point computational domain.
    aMean = reproducing_sum(array(HI%isc:HI%iec,HI%jsc:HI%jec,:))
    n = (1 + HI%jec - HI%jsc) * (1 + HI%iec - HI%isc) * size(array,3)
    call sum_across_PEs(n)
    call min_across_PEs(aMin)
    call max_across_PEs(aMax)
    aMean = aMean / real(n)
  end subroutine subStats

end procedure chksum_u_3d
module procedure chksum_v_3d
  real, pointer :: array(:,:,:)         ! Field array on the input grid [A ~> a]
  real, allocatable, dimension(:,:,:) :: rescaled_array ! The array with scaling undone [a]
  real, allocatable :: hash_array(:,:,:)  ! Subarray used to compute hash [a]
  type(hor_index_type), pointer :: HI   ! Horizontal index bounds of the input grid
  real :: scaling   ! Explicit rescaling factor [a A-1 ~> 1]
  integer :: iounit !< Log IO unit
  integer :: i, j, k, Js
  integer :: bc0, bcSW, bcSE, bcNW, bcNE, hshift
  integer :: bcN, bcS, bcE, bcW
  real :: aMean, aMin, aMax  ! Array mean, global minimum and global maximum [a]
  logical :: do_corners, sym, sym_stats
  integer :: turns                      ! Quarter turns from input to model grid
  turns = HI_m%turns
  if (modulo(turns, 4) /= 0) then
    allocate(HI)
    call rotate_hor_index(HI_m, -turns, HI)
    if (modulo(turns, 2) /= 0) then
      ! Arrays originating from u-points must be handled by uchksum
      allocate(array(HI%IsdB:HI%IedB, HI%jsd:HI%jed, size(array_m, 3)))
      call rotate_array(array_m, -turns, array)
      call uchksum(array, mesg, HI, haloshift, symmetric, omit_corners, &
                   scale=scale, logunit=logunit, unscale=unscale)
      return
    else
      allocate(array(HI%isd:HI%ied, HI%JsdB:HI%JedB, size(array_m, 3)))
      call rotate_array(array_m, -turns, array)
    endif
  else
    HI => HI_m
    array => array_m
  endif

  if (checkForNaNs) then
    if (is_NaN(array(HI%isc:HI%iec,HI%JscB:HI%JecB,:))) &
      call chksum_error(FATAL, 'NaN detected: '//trim(mesg))
!   if (is_NaN(array)) &
!     call chksum_error(FATAL, 'NaN detected in halo: '//trim(mesg))
  endif

  scaling = 1.0
  if (present(unscale)) then ; scaling = unscale
  elseif (present(scale)) then ; scaling = scale ; endif

  iounit = error_unit ; if (present(logunit)) iounit = logunit
  sym_stats = .false. ; if (present(symmetric)) sym_stats = symmetric
  if (present(haloshift)) then ; if (haloshift > 0) sym_stats = .true. ; endif

  if (calculateStatistics) then
    if (present(unscale) .or. present(scale)) then
      allocate( rescaled_array(LBOUND(array,1):UBOUND(array,1), &
                               LBOUND(array,2):UBOUND(array,2), &
                               LBOUND(array,3):UBOUND(array,3)), source=0.0 )
      Js = HI%jsc ; if (sym_stats) Js = HI%jsc-1
      do k=1,size(array,3) ; do J=Js,HI%JecB ; do i=HI%isc,HI%iec
        rescaled_array(i,J,k) = scaling*array(i,J,k)
      enddo ; enddo ; enddo
      call subStats(HI, rescaled_array, sym_stats, aMean, aMin, aMax)
      deallocate(rescaled_array)
    else
      call subStats(HI, array, sym_stats, aMean, aMin, aMax)
    endif
    if (is_root_pe()) &
      call chk_sum_msg("v-point:", aMean, aMin, aMax, mesg, iounit)
  endif

  if (.not.writeChksums) return

  hshift = default_shift
  if (present(haloshift)) hshift = haloshift
  if (hshift<0) hshift = HI%ied-HI%iec

  if ( HI%isc-hshift<HI%isd .or. HI%iec+hshift>HI%ied .or. &
       HI%jsc-hshift<HI%jsd .or. HI%jec+hshift>HI%jed ) then
    write(0,*) 'chksum_v_3d: haloshift =',hshift
    write(0,*) 'chksum_v_3d: isd,isc,iec,ied=',HI%isd,HI%isc,HI%iec,HI%ied
    write(0,*) 'chksum_v_3d: jsd,jsc,jec,jed=',HI%jsd,HI%jsc,HI%jec,HI%jed
    call chksum_error(FATAL,'Error in chksum_v_3d '//trim(mesg))
  endif

  bc0 = subchk(array, HI, 0, 0, scaling)

  sym = .false. ; if (present(symmetric)) sym = symmetric

  if ((hshift==0) .and. .not.sym) then
    if (is_root_pe()) call chk_sum_msg("v-point:", bc0, mesg, iounit)
  else
    do_corners = .true.
    if (present(omit_corners)) do_corners = .not. omit_corners

    if (hshift==0) then
      bcS = subchk(array, HI, 0, -hshift-1, scaling)
      if (is_root_pe()) call chk_sum_msg_S("v-point:", bc0, bcS, mesg, iounit)
    elseif (do_corners) then
      if (sym) then
        bcSW = subchk(array, HI, -hshift, -hshift-1, scaling)
        bcSE = subchk(array, HI, hshift, -hshift-1, scaling)
      else
        bcSW = subchk(array, HI, -hshift, -hshift, scaling)
        bcSE = subchk(array, HI, hshift, -hshift, scaling)
      endif
      bcNW = subchk(array, HI, -hshift, hshift, scaling)
      bcNE = subchk(array, HI, hshift, hshift, scaling)

      if (is_root_pe()) &
        call chk_sum_msg("v-point:", bc0, bcSW, bcSE, bcNW, bcNE, mesg, iounit)
    else
      if (sym) then
        bcS = subchk(array, HI, 0, -hshift-1, scaling)
      else
        bcS = subchk(array, HI, 0, -hshift, scaling)
      endif
      bcE = subchk(array, HI, hshift, 0, scaling)
      bcW = subchk(array, HI, -hshift, 0, scaling)
      bcN = subchk(array, HI, 0, hshift, scaling)

      if (is_root_pe()) &
        call chk_sum_msg_NSEW("v-point:", bc0, bcN, bcS, bcE, bcW, mesg, iounit)
    endif
  endif

  if (writeHash .and. is_root_pe()) then
    allocate(hash_array(HI%isc:HI%iec, HI%jsc:HI%jec, size(array, 3)))
    hash_array(:,:,:) = scaling * array(HI%isc:HI%iec, HI%jsc:HI%jec, :)

    write(iounit, '("v-point: hash=", z8, 1x, a)') &
        murmur_hash(hash_array), mesg
    deallocate(hash_array)
  endif

  contains

  integer function subchk(array, HI, di, dj, unscale)
    type(hor_index_type), intent(in) ::  HI     !< A horizontal index type
    real, dimension(HI%isd:,HI%JsdB:,:), intent(in) :: array !< The array to be checksummed in
                                 !! arbitrary, possibly rescaled units [A ~> a]
    integer, intent(in) :: di    !< i- direction array shift for this checksum
    integer, intent(in) :: dj    !< j- direction array shift for this checksum
    real, intent(in)    :: unscale !< A factor to convert this array back to unscaled units
                                 !! for checksums and output [a A-1 ~> 1]
    integer :: i, j, k, bc
    subchk = 0
    ! This line deliberately uses the h-point computational domain.
    do k=LBOUND(array,3),UBOUND(array,3) ; do J=HI%jsc+dj,HI%jec+dj ; do i=HI%isc+di,HI%iec+di
      bc = bitcount(abs(unscale*array(i,J,k)))
      subchk = subchk + bc
    enddo ; enddo ; enddo
    call sum_across_PEs(subchk)
    subchk=mod(subchk, bc_modulus)
  end function subchk

  !subroutine subStats(HI, array, mesg, sym_stats)
  subroutine subStats(HI, array, sym_stats, aMean, aMin, aMax)
    type(hor_index_type), intent(in) ::  HI     !< A horizontal index type
    real, dimension(HI%isd:,HI%JsdB:,:), intent(in) :: array !< The array to be checksummed [a]
    logical,          intent(in) :: sym_stats !< If true, evaluate the statistics on the
                                              !! full symmetric computational domain.
    real, intent(out) :: aMean   !< Mean of array over domain [a]
    real, intent(out) :: aMin    !< Minimum of array over domain [a]
    real, intent(out) :: aMax    !< Maximum of array over domain [a]

    integer :: i, j, k, n, JsB

    JsB = HI%jsc ; if (sym_stats) JsB = HI%jsc-1

    aMin = array(HI%isc,HI%jsc,1) ; aMax = aMin
    do k=LBOUND(array,3),UBOUND(array,3) ; do J=JsB,HI%JecB ; do i=HI%isc,HI%iec
      aMin = min(aMin, array(i,J,k))
      aMax = max(aMax, array(i,J,k))
    enddo ; enddo ; enddo
    ! This line deliberately uses the h-point computational domain.
    aMean = reproducing_sum(array(HI%isc:HI%iec,HI%jsc:HI%jec,:))
    n = (1 + HI%jec - HI%jsc) * (1 + HI%iec - HI%isc) * size(array,3)
    call sum_across_PEs(n)
    call min_across_PEs(aMin)
    call max_across_PEs(aMax)
    aMean = aMean / real(n)
  end subroutine subStats

end procedure chksum_v_3d
module procedure chksum1d
  integer :: is, ie, i, bc, sum1, sum_bc, ioUnit
  real :: sum  ! The global sum of the array [A]
  real, allocatable :: sum_here(:) ! The sum on each PE [A]
  logical :: compare
  integer :: pe_num   ! pe number of the data
  integer :: nPEs     ! Total number of processsors
  is = LBOUND(array,1) ; ie = UBOUND(array,1)
  if (present(start_i)) is = start_i
  if (present(end_i)) ie = end_i
  compare = .true. ; if (present(compare_PEs)) compare = compare_PEs
  iounit = error_unit ; if (present(logunit)) iounit = logunit

  sum = 0.0 ; sum_bc = 0
  do i=is,ie
    sum = sum + array(i)
    bc = bitcount(ABS(array(i)))
    sum_bc = sum_bc + bc
  enddo

  pe_num = pe_here() + 1 - root_pe() ; nPEs = num_pes()
  allocate(sum_here(nPEs), source=0.0) ; sum_here(pe_num) = sum
  call sum_across_PEs(sum_here,nPEs)

  sum1 = sum_bc
  call sum_across_PEs(sum1)

  if (.not.compare) then
    sum = 0.0
    do i=1,nPEs ; sum = sum + sum_here(i) ; enddo
    sum_bc = sum1
  elseif (is_root_pe()) then
    if (sum1 /= nPEs*sum_bc) &
      write(iounit, '(A40," bitcounts do not match across PEs: ",I12,1X,I12)') &
            mesg, sum1, nPEs*sum_bc
    do i=1,nPEs ; if (sum /= sum_here(i)) then
      write(iounit, '(A40," PE ",I0," sum mismatches root_PE: ",3(ES22.13,1X))') &
            mesg, i, sum_here(i), sum, sum_here(i)-sum
    endif ; enddo
  endif
  deallocate(sum_here)

  if (is_root_pe()) &
    write(iounit,'(A50,1X,ES25.16,1X,I12)') mesg, sum, sum_bc

end procedure chksum1d
module procedure chksum2d
  integer :: xs, xe, ys, ye, i, j, sum1, bc, iounit
  real :: sum  ! The global sum of the array [A]
  iounit = error_unit ; if (present(logunit)) iounit = logunit

  xs = LBOUND(array,1) ; xe = UBOUND(array,1)
  ys = LBOUND(array,2) ; ye = UBOUND(array,2)

  sum = 0.0 ; sum1 = 0
  do i=xs,xe ; do j=ys,ye
    bc = bitcount(abs(array(i,j)))
    sum1 = sum1 + bc
  enddo ; enddo
  call sum_across_PEs(sum1)

  sum = reproducing_sum(array(:,:))

  if (is_root_pe()) &
    write(iounit,'(A50,1X,ES25.16,1X,I12)') mesg, sum, sum1
!    write(iounit,'(A40,1X,Z16.16,1X,Z16.16,1X,ES25.16,1X,I12)') &
!      mesg, sum, sum1, sum, sum1

end procedure chksum2d
module procedure chksum3d
  integer :: xs, xe, ys, ye, zs, ze, i, j, k, bc, sum1, iounit
  real :: sum  ! The global sum of the array [A]
  iounit = error_unit ; if (present(logunit)) iounit = logunit

  xs = LBOUND(array,1) ; xe = UBOUND(array,1)
  ys = LBOUND(array,2) ; ye = UBOUND(array,2)
  zs = LBOUND(array,3) ; ze = UBOUND(array,3)

  sum = 0.0 ; sum1 = 0
  do i=xs,xe ; do j=ys,ye ; do k=zs,ze
    bc = bitcount(ABS(array(i,j,k)))
    sum1 = sum1 + bc
  enddo ; enddo ; enddo

  call sum_across_PEs(sum1)
  sum = reproducing_sum(array(:,:,:))

  if (is_root_pe()) &
    write(iounit, '(A50,1X,ES25.16,1X,I12)') mesg, sum, sum1
!    write(iounit, '(A40,1X,Z16.16,1X,Z16.16,1X,ES25.16,1X,I12)') &
!      mesg, sum, sum1, sum, sum1

end procedure chksum3d
module procedure is_NaN_0d
  if (((x < 0.0) .and. (x >= 0.0)) .or. &
            (.not.(x < 0.0) .and. .not.(x >= 0.0))) then
    is_NaN_0d = .true.
  else
    is_NaN_0d = .false.
  endif

end procedure is_NaN_0d
module procedure is_NaN_1d
  integer :: i, n
  logical :: global_check
  n = 0
  do i = LBOUND(x,1), UBOUND(x,1)
    if (is_NaN_0d(x(i))) n = n + 1
  enddo
  global_check = .true.
  if (present(skip_mpp)) global_check = .not.skip_mpp

  if (global_check) call sum_across_PEs(n)
  is_NaN_1d = .false.
  if (n>0) is_NaN_1d = .true.

end procedure is_NaN_1d
module procedure is_NaN_2d
  integer :: i, j, n
  n = 0
  do j = LBOUND(x,2), UBOUND(x,2) ; do i = LBOUND(x,1), UBOUND(x,1)
    if (is_NaN_0d(x(i,j))) n = n + 1
  enddo ; enddo
  call sum_across_PEs(n)
  is_NaN_2d = .false.
  if (n>0) is_NaN_2d = .true.

end procedure is_NaN_2d
module procedure is_NaN_3d
  integer :: i, j, k, n
  n = 0
  do k = LBOUND(x,3), UBOUND(x,3)
    do j = LBOUND(x,2), UBOUND(x,2) ; do i = LBOUND(x,1), UBOUND(x,1)
      if (is_NaN_0d(x(i,j,k))) n = n + 1
    enddo ; enddo
  enddo
  call sum_across_PEs(n)
  is_NaN_3d = .false.
  if (n>0) is_NaN_3d = .true.

end procedure is_NaN_3d
module procedure field_checksum_real_0d
  real :: scale_fac  ! A local copy of unscale if it is present [a A-1 ~> 1] or 1 otherwise
  if (present(turns)) call MOM_error(FATAL, "Rotation not supported for 0d fields.")

  scale_fac = 1.0 ; if (present(unscale)) scale_fac = unscale

  chksum = field_chksum(scale_fac*field, pelist=pelist, mask_val=mask_val)
end procedure field_checksum_real_0d
module procedure field_checksum_real_1d
  real :: scale_fac  ! A local copy of unscale if it is present [a A-1 ~> 1] or 1 otherwise
  if (present(turns)) call MOM_error(FATAL, "Rotation not supported for 1d fields.")

  scale_fac = 1.0 ; if (present(unscale)) scale_fac = unscale

  chksum = field_chksum(scale_fac*field(:), pelist=pelist, mask_val=mask_val)
end procedure field_checksum_real_1d
module procedure field_checksum_real_2d
  real, allocatable :: field_rot(:,:)  ! A rotated version of field, with the same units [A ~> a]
  integer :: qturns ! The number of quarter turns through which to rotate field
  logical :: do_unscale ! If true, unscale the variable before it is checksummed
  qturns = 0
  if (present(turns)) &
    qturns = modulo(turns, 4)

  do_unscale = .false. ; if (present(unscale)) do_unscale = (unscale /= 1.0)

  if (qturns == 0) then
    if (do_unscale) then
      chksum = field_chksum(unscale*field(:,:), pelist=pelist, mask_val=mask_val)
    else
      chksum = field_chksum(field, pelist=pelist, mask_val=mask_val)
    endif
  else
    call allocate_rotated_array(field, [1,1], qturns, field_rot)
    call rotate_array(field, qturns, field_rot)
    if (do_unscale) field_rot(:,:) = unscale*field_rot(:,:)
    chksum = field_chksum(field_rot, pelist=pelist, mask_val=mask_val)
    deallocate(field_rot)
  endif
end procedure field_checksum_real_2d
module procedure field_checksum_real_3d
  real, allocatable :: field_rot(:,:,:)  ! A rotated version of field, with the same units [A ~> a]
  integer :: qturns ! The number of quarter turns through which to rotate field
  logical :: do_unscale ! If true, unscale the variable before it is checksummed
  qturns = 0
  if (present(turns)) &
    qturns = modulo(turns, 4)

  do_unscale = .false. ; if (present(unscale)) do_unscale = (unscale /= 1.0)

  if (qturns == 0) then
    if (do_unscale) then
      chksum = field_chksum(unscale*field(:,:,:), pelist=pelist, mask_val=mask_val)
    else
      chksum = field_chksum(field, pelist=pelist, mask_val=mask_val)
    endif
  else
    call allocate_rotated_array(field, [1,1,1], qturns, field_rot)
    call rotate_array(field, qturns, field_rot)
    if (do_unscale) field_rot(:,:,:) = unscale*field_rot(:,:,:)
    chksum = field_chksum(field_rot, pelist=pelist, mask_val=mask_val)
    deallocate(field_rot)
  endif
end procedure field_checksum_real_3d
module procedure field_checksum_real_4d
  real, allocatable :: field_rot(:,:,:,:)  ! A rotated version of field, with the same units [A ~> a]
  integer :: qturns ! The number of quarter turns through which to rotate field
  logical :: do_unscale ! If true, unscale the variable before it is checksummed
  qturns = 0
  if (present(turns)) &
    qturns = modulo(turns, 4)

  do_unscale = .false. ; if (present(unscale)) do_unscale = (unscale /= 1.0)

  if (qturns == 0) then
    if (do_unscale) then
      chksum = field_chksum(unscale*field(:,:,:,:), pelist=pelist, mask_val=mask_val)
    else
      chksum = field_chksum(field, pelist=pelist, mask_val=mask_val)
    endif
  else
    call allocate_rotated_array(field, [1,1,1,1], qturns, field_rot)
    call rotate_array(field, qturns, field_rot)
    if (do_unscale) field_rot(:,:,:,:) = unscale*field_rot(:,:,:,:)
    chksum = field_chksum(field_rot, pelist=pelist, mask_val=mask_val)
    deallocate(field_rot)
  endif
end procedure field_checksum_real_4d
module procedure chk_sum_msg1
  if (is_root_pe()) &
    write(iounit, '(a,1(a,i10,1x),a)') fmsg, " c=", bc0, trim(mesg)
end procedure chk_sum_msg1
module procedure chk_sum_msg5
  if (is_root_pe()) write(iounit, '(A,5(A,I10,1X),A)') &
    fmsg, " c=", bc0, "sw=", bcSW, "se=", bcSE, "nw=", bcNW, "ne=", bcNE, trim(mesg)
end procedure chk_sum_msg5
module procedure chk_sum_msg_NSEW
  if (is_root_pe()) write(iounit, '(A,5(A,I10,1X),A)') &
    fmsg, " c=", bc0, "N=", bcN, "S=", bcS, "E=", bcE, "W=", bcW, trim(mesg)
end procedure chk_sum_msg_NSEW
module procedure chk_sum_msg_S
  if (is_root_pe()) write(iounit, '(A,2(A,I10,1X),A)') &
    fmsg, " c=", bc0, "S=", bcS, trim(mesg)
end procedure chk_sum_msg_S
module procedure chk_sum_msg_W
  if (is_root_pe()) write(iounit, '(A,2(A,I10,1X),A)') &
    fmsg, " c=", bc0, "W=", bcW, trim(mesg)
end procedure chk_sum_msg_W
module procedure chk_sum_msg2
  if (is_root_pe()) write(iounit, '(A,2(A,I9,1X),A)') &
    fmsg, " c=", bc0, "s/w=", bcSW, trim(mesg)
end procedure chk_sum_msg2
module procedure chk_sum_msg3
  if (is_root_pe()) write(iounit, '(A,3(A,ES25.16,1X),A)') &
    fmsg, " mean=", aMean, "min=", (0. + aMin), "max=", (0. + aMax), trim(mesg)
end procedure chk_sum_msg3
module procedure MOM_checksums_init
# include "version_variable.h"
  character(len=40)  :: mdl = "MOM_checksums" ! This module's name.
  call log_version(param_file, mdl, version)

end procedure MOM_checksums_init
module procedure chksum_error
  call MOM_error(signal, message)
end procedure chksum_error
module procedure bitcount
  integer, parameter :: xk = kind(x)  !< Kind type of x
  bitcount = popcnt(transfer(x, 1_xk))
end procedure bitcount
end submodule MOM_checksums_s
