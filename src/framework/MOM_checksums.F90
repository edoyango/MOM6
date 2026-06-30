! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> Routines to calculate checksums of various array and vector types
module MOM_checksums

use MOM_array_transform, only : rotate_array, rotate_array_pair, rotate_vector
use MOM_array_transform, only : allocate_rotated_array
use MOM_coms,            only : PE_here, root_PE, num_PEs, sum_across_PEs
use MOM_coms,            only : min_across_PEs, max_across_PEs
use MOM_coms,            only : reproducing_sum, field_chksum
use MOM_error_handler,   only : MOM_error, FATAL, is_root_pe
use MOM_file_parser,     only : log_version, param_file_type
use MOM_hor_index,       only : hor_index_type, rotate_hor_index
use MOM_murmur_hash,     only : murmur_hash

use iso_fortran_env,     only : error_unit, int32, int64

implicit none ; private

public :: chksum0, zchksum, rotated_field_chksum, field_checksum
public :: hchksum, Bchksum, uchksum, vchksum, qchksum, is_NaN, chksum
public :: hchksum_pair, uvchksum, Bchksum_pair
public :: MOM_checksums_init

! A note on unit descriptions in comments: MOM6 uses units that can be rescaled for dimensional
! consistency testing. These are noted in comments with units like Z, H, L, and T, along with
! their mks counterparts with notation like "a velocity [Z T-1 ~> m s-1]".  If the units
! vary with the Boussinesq approximation, the Boussinesq variant is given first.
! The functions in this module work with variables with arbitrary units, in which case the
! arbitrary rescaled units are indicated with [A ~> a], while the unscaled units are just [a].

!> Checksums a pair of arrays (2d or 3d) staggered at tracer points
interface hchksum_pair
  module procedure chksum_pair_h_2d, chksum_pair_h_3d
end interface

!> Checksums a pair velocity arrays (2d or 3d) staggered at C-grid locations
interface uvchksum
  module procedure chksum_uv_2d, chksum_uv_3d
end interface

!> Checksums an array (2d or 3d) staggered at C-grid u points.
interface uchksum
  module procedure chksum_u_2d, chksum_u_3d
end interface

!> Checksums an array (2d or 3d) staggered at C-grid v points.
interface vchksum
  module procedure chksum_v_2d, chksum_v_3d
end interface

!> Checksums a pair of arrays (2d or 3d) staggered at corner points
interface Bchksum_pair
  module procedure chksum_pair_B_2d, chksum_pair_B_3d
end interface

!> Checksums an array (2d or 3d) staggered at tracer points.
interface hchksum
  module procedure chksum_h_2d, chksum_h_3d
end interface

!> Checksums an array (2d or 3d) staggered at corner points.
interface Bchksum
  module procedure chksum_B_2d, chksum_B_3d
end interface

!> This is an older interface that has been renamed Bchksum
interface qchksum
  module procedure chksum_B_2d, chksum_B_3d
end interface

!> This is an older interface for 1-, 2-, or 3-D checksums
interface chksum
  module procedure chksum1d, chksum2d, chksum3d
end interface

!> Write a message with either checksums or numerical statistics of arrays
interface chk_sum_msg
  module procedure chk_sum_msg1, chk_sum_msg2, chk_sum_msg3, chk_sum_msg5
end interface

!> Returns .true. if any element of x is a NaN, and .false. otherwise.
interface is_NaN
  module procedure is_NaN_0d, is_NaN_1d, is_NaN_2d, is_NaN_3d
end interface

!> Compute the checksum on all elements of a field that may need to be rotated or unscaled.
!! This interface uses the field_chksum function that is used to verify file contents, which
!! may differ from the bitcount function used for other checksums in this module.
interface rotated_field_chksum
  module procedure field_checksum_real_0d
  module procedure field_checksum_real_1d
  module procedure field_checksum_real_2d
  module procedure field_checksum_real_3d
  module procedure field_checksum_real_4d
end interface rotated_field_chksum


!> Compute the checksum on all elements of a field that may need to be rotated or unscaled.
!! This interface uses the field_chksum function that is used to verify file contents, which
!! may differ from the bitcount function used for other checksums in this module.
interface field_checksum
  module procedure field_checksum_real_0d
  module procedure field_checksum_real_1d
  module procedure field_checksum_real_2d
  module procedure field_checksum_real_3d
  module procedure field_checksum_real_4d
end interface field_checksum

integer, parameter :: bc_modulus = 1000000000 !< Modulus of checksum bitcount
integer, parameter :: default_shift=0 !< The default array shift
logical :: calculateStatistics=.true. !< If true, report min, max and mean.
logical :: writeChksums=.true. !< If true, report the bitcount checksum
logical :: checkForNaNs=.true. !< If true, checks array for NaNs and cause
                               !! FATAL error if any are found
logical :: writeHash = .false. !< If true, report the murmur hash
  !! NOTE: writeHash is currently disabled due to non-compliant diagnostics.


  interface
module subroutine chksum0(scalar, mesg, scale, logunit, unscale)
  real,              intent(in) :: scalar  !< The array to be checksummed in
                                           !! arbitrary, possibly rescaled units [A ~> a]
  character(len=*),  intent(in) :: mesg    !< An identifying message
  real,    optional, intent(in) :: scale   !< A factor to convert this array back to unscaled units
                                           !! for checksums and output [a A-1 ~> 1]
  integer, optional, intent(in) :: logunit !< IO unit for checksum logging
  real,    optional, intent(in) :: unscale !< A factor to convert this array back to unscaled units
                                           !! for checksums and output [a A-1 ~> 1].
                                           !! Here scale and unscale are synonymous, but unscale
                                           !! takes precedence if both are present.

  ! Local variables
  ! In the following comments, [A] is used to indicate the arbitrary, possibly rescaled units
  ! of the input array while [a] indicates the unscaled (e.g., mks) units that should be used
  ! for checksums and output

end subroutine chksum0
module subroutine zchksum(array, mesg, scale, logunit, unscale)
  real, dimension(:), intent(in) :: array   !< The array to be checksummed in
                                            !! arbitrary, possibly rescaled units [A ~> a]
  character(len=*),   intent(in) :: mesg    !< An identifying message
  real,     optional, intent(in) :: scale   !< A factor to convert this array back to unscaled units
                                            !! for checksums and output [a A-1 ~> 1]
  integer,  optional, intent(in) :: logunit !< IO unit for checksum logging
  real,     optional, intent(in) :: unscale !< A factor to convert this array back to unscaled units
                                            !! for checksums and output [a A-1 ~> 1].
                                            !! Here scale and unscale are synonymous, but unscale
                                            !! takes precedence if both are present.

  ! Local variables
  ! In the following comments, [A] is used to indicate the arbitrary, possibly rescaled units
  ! of the input array while [a] indicates the unscaled (e.g., mks) units that should be used
  ! for checksums and output

end subroutine zchksum
module subroutine chksum_pair_h_2d(mesg, arrayA, arrayB, HI, haloshift, omit_corners, &
                            scale, logunit, scalar_pair, unscale)
  character(len=*),                 intent(in) :: mesg !< Identifying messages
  type(hor_index_type),   target,   intent(in) :: HI     !< A horizontal index type
  real, dimension(HI%isd:,HI%jsd:), target, intent(in) :: arrayA !< The first array to be checksummed in
                                                            !! arbitrary, possibly rescaled units [A ~> a]
  real, dimension(HI%isd:,HI%jsd:), target, intent(in) :: arrayB !< The second array to be checksummed in
                                                            !! arbitrary, possibly rescaled units [A ~> a]
  integer,                optional, intent(in) :: haloshift !< The width of halos to check (default 0)
  logical,                optional, intent(in) :: omit_corners !< If true, avoid checking diagonal shifts
  real,                   optional, intent(in) :: scale     !< A factor to convert these arrays back to unscaled
                                                            !! units for checksums and output [a A-1 ~> 1]
  integer,                optional, intent(in) :: logunit   !< IO unit for checksum logging
  logical,                optional, intent(in) :: scalar_pair !< If true, then the arrays describe
                                                            !! a scalar, rather than vector
  real,                   optional, intent(in) :: unscale   !< A factor to convert this array back to unscaled units
                                                            !! for checksums and output [a A-1 ~> 1].
                                                            !! Here scale and unscale are synonymous, but unscale
                                                            !! takes precedence if both are present.

end subroutine chksum_pair_h_2d
module subroutine chksum_pair_h_3d(mesg, arrayA, arrayB, HI, haloshift, omit_corners, &
                            scale, logunit, scalar_pair, unscale)
  character(len=*),                    intent(in) :: mesg !< Identifying messages
  type(hor_index_type),      target,   intent(in) :: HI   !< A horizontal index type
  real, dimension(HI%isd:,HI%jsd:, :), target, intent(in) :: arrayA !< The first array to be checksummed in
                                                             !! arbitrary, possibly rescaled units [A ~> a]
  real, dimension(HI%isd:,HI%jsd:, :), target, intent(in) :: arrayB !< The second array to be checksummed in
                                                             !! arbitrary, possibly rescaled units [A ~> a]
  integer,                   optional, intent(in) :: haloshift !< The width of halos to check (default 0)
  logical,                   optional, intent(in) :: omit_corners !< If true, avoid checking diagonal shifts
  real,                      optional, intent(in) :: scale     !< A factor to convert this array back to unscaled units
                                                               !! for checksums and output [a A-1 ~> 1]
  integer,                   optional, intent(in) :: logunit   !< IO unit for checksum logging
  logical,                   optional, intent(in) :: scalar_pair !< If true, then the arrays describe
                                                               !! a scalar, rather than vector
  real,                      optional, intent(in) :: unscale   !< A factor to convert this array back to unscaled units
                                                               !! for checksums and output [a A-1 ~> 1].
                                                               !! Here scale and unscale are synonymous, but unscale
                                                               !! takes precedence if both are present.
  ! Local variables

end subroutine chksum_pair_h_3d
module subroutine chksum_h_2d(array_m, mesg, HI_m, haloshift, omit_corners, scale, logunit, unscale)
  type(hor_index_type), target, intent(in) :: HI_m         !< Horizontal index bounds of the model grid
  real, dimension(HI_m%isd:,HI_m%jsd:), target, intent(in) :: array_m !< Field array on the model grid in
                                                           !! arbitrary, possibly rescaled units [A ~> a]
  character(len=*),                intent(in) :: mesg      !< An identifying message
  integer,               optional, intent(in) :: haloshift !< The width of halos to check (default 0)
  logical,               optional, intent(in) :: omit_corners !< If true, avoid checking diagonal shifts
  real,                  optional, intent(in) :: scale     !< A factor to convert this array back to unscaled units
                                                           !! for checksums and output [a A-1 ~> 1]
  integer,               optional, intent(in) :: logunit   !< IO unit for checksum logging
  real,                  optional, intent(in) :: unscale   !< A factor to convert this array back to unscaled units
                                                           !! for checksums and output [a A-1 ~> 1].
                                                           !! Here scale and unscale are synonymous, but unscale
                                                           !! takes precedence if both are present.

  ! Local variables
  ! In the following comments, [A] is used to indicate the arbitrary, possibly rescaled units
  ! of the input array while [a] indicates the unscaled (e.g., mks) units that should be used
  ! for checksums and output


  ! Rotate array to the input grid
end subroutine chksum_h_2d
module subroutine chksum_pair_B_2d(mesg, arrayA, arrayB, HI, haloshift, symmetric, &
                            omit_corners, scale, logunit, scalar_pair, unscale)
  character(len=*),                 intent(in) :: mesg   !< Identifying messages
  type(hor_index_type),   target,   intent(in) :: HI     !< A horizontal index type
  real, dimension(HI%isd:,HI%jsd:), target, intent(in) :: arrayA !< The first array to be checksummed in
                                                            !! arbitrary, possibly rescaled units [A ~> a]
  real, dimension(HI%isd:,HI%jsd:), target, intent(in) :: arrayB !< The second array to be checksummed in
                                                            !! arbitrary, possibly rescaled units [A ~> a]
  logical,                optional, intent(in) :: symmetric !< If true, do the checksums on the full
                                                            !! symmetric computational domain.
  integer,                optional, intent(in) :: haloshift !< The width of halos to check (default 0)
  logical,                optional, intent(in) :: omit_corners !< If true, avoid checking diagonal shifts
  real,                   optional, intent(in) :: scale     !< A factor to convert this array back to unscaled units
                                                            !! for checksums and output [a A-1 ~> 1]
  integer,                optional, intent(in) :: logunit   !< IO unit for checksum logging
  logical,                optional, intent(in) :: scalar_pair !< If true, then the arrays describe
                                                            !! a scalar, rather than vector
  real,                   optional, intent(in) :: unscale   !< A factor to convert this array back to unscaled units
                                                            !! for checksums and output [a A-1 ~> 1].
                                                            !! Here scale and unscale are synonymous, but unscale
                                                            !! takes precedence if both are present.


end subroutine chksum_pair_B_2d
module subroutine chksum_pair_B_3d(mesg, arrayA, arrayB, HI, haloshift, symmetric, &
                            omit_corners, scale, logunit, scalar_pair, unscale)
  character(len=*),                    intent(in) :: mesg !< Identifying messages
  type(hor_index_type),      target,   intent(in) :: HI     !< A horizontal index type
  real, dimension(HI%IsdB:,HI%JsdB:, :), target, intent(in) :: arrayA !< The first array to be checksummed in
                                                               !! arbitrary, possibly rescaled units [A ~> a]
  real, dimension(HI%IsdB:,HI%JsdB:, :), target, intent(in) :: arrayB !< The second array to be checksummed in
                                                               !! arbitrary, possibly rescaled units [A ~> a]
  integer,                   optional, intent(in) :: haloshift !< The width of halos to check (default 0)
  logical,                   optional, intent(in) :: symmetric !< If true, do the checksums on the full
                                                               !! symmetric computational domain.
  logical,                   optional, intent(in) :: omit_corners !< If true, avoid checking diagonal shifts
  real,                      optional, intent(in) :: scale     !< A factor to convert this array back to unscaled units
                                                               !! for checksums and output [a A-1 ~> 1]
  integer,                   optional, intent(in) :: logunit   !< IO unit for checksum logging
  logical,                   optional, intent(in) :: scalar_pair !< If true, then the arrays describe
                                                               !! a scalar, rather than vector
  real,                      optional, intent(in) :: unscale   !< A factor to convert this array back to unscaled units
                                                               !! for checksums and output [a A-1 ~> 1].
                                                               !! Here scale and unscale are synonymous, but unscale
                                                               !! takes precedence if both are present.
  ! Local variables

end subroutine chksum_pair_B_3d
module subroutine chksum_B_2d(array_m, mesg, HI_m, haloshift, symmetric, omit_corners, &
                       scale, logunit, unscale)
  type(hor_index_type), target, intent(in) :: HI_m     !< A horizontal index type
  real, dimension(HI_m%IsdB:,HI_m%JsdB:), &
                        target, intent(in) :: array_m !< The array to be checksummed in
                                                !! arbitrary, possibly rescaled units [A ~> a]
  character(len=*),     intent(in) :: mesg      !< An identifying message
  integer,    optional, intent(in) :: haloshift !< The width of halos to check (default 0)
  logical,    optional, intent(in) :: symmetric !< If true, do the checksums on the
                                                !! full symmetric computational domain.
  logical,    optional, intent(in) :: omit_corners !< If true, avoid checking diagonal shifts
  real,       optional, intent(in) :: scale     !< A factor to convert this array back to unscaled units
                                                !! for checksums and output [a A-1 ~> 1]
  integer,    optional, intent(in) :: logunit   !< IO unit for checksum logging
  real,       optional, intent(in) :: unscale   !< A factor to convert this array back to unscaled units
                                                !! for checksums and output [a A-1 ~> 1].
                                                !! Here scale and unscale are synonymous, but unscale
                                                !! takes precedence if both are present.

  ! Local variables
  ! In the following comments, [A] is used to indicate the arbitrary, possibly rescaled units
  ! of the input array while [a] indicates the unscaled (e.g., mks) units that should be used
  ! for checksums and output

  ! Rotate array to the input grid
end subroutine chksum_B_2d
module subroutine chksum_uv_2d(mesg, arrayU, arrayV, HI, haloshift, symmetric, &
                        omit_corners, scale, logunit, scalar_pair, unscale)
  character(len=*),                  intent(in) :: mesg   !< Identifying messages
  type(hor_index_type),    target,   intent(in) :: HI     !< A horizontal index type
  real, dimension(HI%IsdB:,HI%jsd:), target, intent(in) :: arrayU !< The u-component array to be checksummed in
                                                             !! arbitrary, possibly rescaled units [A ~> a]
  real, dimension(HI%isd:,HI%JsdB:), target, intent(in) :: arrayV !< The v-component array to be checksummed in
                                                             !! arbitrary, possibly rescaled units [A ~> a]
  integer,                 optional, intent(in) :: haloshift !< The width of halos to check (default 0)
  logical,                 optional, intent(in) :: symmetric !< If true, do the checksums on the full
                                                             !! symmetric computational domain.
  logical,                 optional, intent(in) :: omit_corners !< If true, avoid checking diagonal shifts
  real,                    optional, intent(in) :: scale     !< A factor to convert these arrays back to unscaled
                                                             !! units for checksums and output [a A-1 ~> 1]
  integer,                 optional, intent(in) :: logunit   !< IO unit for checksum logging
  logical,                 optional, intent(in) :: scalar_pair !< If true, then the arrays describe a
                                                             !! a scalar, rather than vector
  real,                    optional, intent(in) :: unscale   !< A factor to convert this array back to unscaled units
                                                             !! for checksums and output [a A-1 ~> 1].
                                                             !! Here scale and unscale are synonymous, but unscale
                                                             !! takes precedence if both are present.
  ! Local variables

end subroutine chksum_uv_2d
module subroutine chksum_uv_3d(mesg, arrayU, arrayV, HI, haloshift, symmetric, &
                        omit_corners, scale, logunit, scalar_pair, unscale)
  character(len=*),                    intent(in) :: mesg   !< Identifying messages
  type(hor_index_type),      target,   intent(in) :: HI     !< A horizontal index type
  real, dimension(HI%IsdB:,HI%jsd:,:), target, intent(in) :: arrayU !< The u-component array to be checksummed in
                                                               !! arbitrary, possibly rescaled units [A ~> a]
  real, dimension(HI%isd:,HI%JsdB:,:), target, intent(in) :: arrayV !< The v-component array to be checksummed in
                                                               !! arbitrary, possibly rescaled units [A ~> a]
  integer,                   optional, intent(in) :: haloshift !< The width of halos to check (default 0)
  logical,                   optional, intent(in) :: symmetric !< If true, do the checksums on the full
                                                               !! symmetric computational domain.
  logical,                   optional, intent(in) :: omit_corners !< If true, avoid checking diagonal shifts
  real,                      optional, intent(in) :: scale     !< A factor to convert these arrays back to unscaled
                                                               !! units for checksums and output [a A-1 ~> 1]
  integer,                   optional, intent(in) :: logunit   !< IO unit for checksum logging
  logical,                 optional, intent(in) :: scalar_pair !< If true, then the arrays describe a
                                                               !! a scalar, rather than vector
  real,                      optional, intent(in) :: unscale   !< A factor to convert this array back to unscaled units
                                                               !! for checksums and output [a A-1 ~> 1].
                                                               !! Here scale and unscale are synonymous, but unscale
                                                               !! takes precedence if both are present.
  ! Local variables

end subroutine chksum_uv_3d
module subroutine chksum_u_2d(array_m, mesg, HI_m, haloshift, symmetric, omit_corners, &
                       scale, logunit, unscale)
  type(hor_index_type),  target,   intent(in) :: HI_m      !< A horizontal index type
  real, dimension(HI_m%IsdB:,HI_m%jsd:), target, intent(in) :: array_m !< The array to be checksummed in
                                                           !! arbitrary, possibly rescaled units [A ~> a]
  character(len=*),                intent(in) :: mesg      !< An identifying message
  integer,               optional, intent(in) :: haloshift !< The width of halos to check (default 0)
  logical,               optional, intent(in) :: symmetric !< If true, do the checksums on the full
                                                           !! symmetric computational domain.
  logical,               optional, intent(in) :: omit_corners !< If true, avoid checking diagonal shifts
  real,                  optional, intent(in) :: scale     !< A factor to convert this array back to unscaled units
                                                           !! for checksums and output [a A-1 ~> 1]
  integer,               optional, intent(in) :: logunit   !< IO unit for checksum logging
  real,                  optional, intent(in) :: unscale   !< A factor to convert this array back to unscaled units
                                                           !! for checksums and output [a A-1 ~> 1].
                                                           !! Here scale and unscale are synonymous, but unscale
                                                           !! takes precedence if both are present.

  ! Local variables
  ! In the following comments, [A] is used to indicate the arbitrary, possibly rescaled units
  ! of the input array while [a] indicates the unscaled (e.g., mks) units that should be used
  ! for checksums and output

  ! Rotate array to the input grid
end subroutine chksum_u_2d
module subroutine chksum_v_2d(array_m, mesg, HI_m, haloshift, symmetric, omit_corners, &
                       scale, logunit, unscale)
  type(hor_index_type),  target,   intent(in) :: HI_m      !< A horizontal index type
  real, dimension(HI_m%isd:,HI_m%JsdB:), target, intent(in) :: array_m !< The array to be checksummed in
                                                           !! arbitrary, possibly rescaled units [A ~> a]
  character(len=*),                intent(in) :: mesg      !< An identifying message
  integer,               optional, intent(in) :: haloshift !< The width of halos to check (default 0)
  logical,               optional, intent(in) :: symmetric !< If true, do the checksums on the full
                                                           !! symmetric computational domain.
  logical,               optional, intent(in) :: omit_corners !< If true, avoid checking diagonal shifts
  real,                  optional, intent(in) :: scale     !< A factor to convert this array back to unscaled units
                                                           !! for checksums and output [a A-1 ~> 1]
  integer,               optional, intent(in) :: logunit   !< IO unit for checksum logging
  real,                  optional, intent(in) :: unscale   !< A factor to convert this array back to unscaled units
                                                           !! for checksums and output [a A-1 ~> 1].
                                                           !! Here scale and unscale are synonymous, but unscale
                                                           !! takes precedence if both are present.

  ! Local variables
  ! In the following comments, [A] is used to indicate the arbitrary, possibly rescaled units
  ! of the input array while [a] indicates the unscaled (e.g., mks) units that should be used
  ! for checksums and output

  ! Rotate array to the input grid
end subroutine chksum_v_2d
module subroutine chksum_h_3d(array_m, mesg, HI_m, haloshift, omit_corners, scale, logunit, unscale)
  type(hor_index_type),    target,   intent(in) :: HI_m !< A horizontal index type
  real, dimension(HI_m%isd:,HI_m%jsd:,:), target, intent(in) :: array_m !< The array to be checksummed in
                                                             !! arbitrary, possibly rescaled units [A ~> a]
  character(len=*),                  intent(in) :: mesg      !< An identifying message
  integer,                 optional, intent(in) :: haloshift !< The width of halos to check (default 0)
  logical,                 optional, intent(in) :: omit_corners !< If true, avoid checking diagonal shifts
  real,                    optional, intent(in) :: scale     !< A factor to convert this array back to unscaled units
                                                             !! for checksums and output [a A-1 ~> 1]
  integer,                 optional, intent(in) :: logunit   !< IO unit for checksum logging
  real,                    optional, intent(in) :: unscale   !< A factor to convert this array back to unscaled units
                                                             !! for checksums and output [a A-1 ~> 1].
                                                             !! Here scale and unscale are synonymous, but unscale
                                                             !! takes precedence if both are present.

  ! Local variables
  ! In the following comments, [A] is used to indicate the arbitrary, possibly rescaled units
  ! of the input array while [a] indicates the unscaled (e.g., mks) units that should be used
  ! for checksums and output

  ! Rotate array to the input grid
end subroutine chksum_h_3d
module subroutine chksum_B_3d(array_m, mesg, HI_m, haloshift, symmetric, omit_corners, &
                       scale, logunit, unscale)
  type(hor_index_type),     target,   intent(in) :: HI_m !< A horizontal index type
  real, dimension(HI_m%IsdB:,HI_m%JsdB:,:), target, intent(in) :: array_m !< The array to be checksummed in
                                                              !! arbitrary, possibly rescaled units [A ~> a]
  character(len=*),                   intent(in) :: mesg      !< An identifying message
  integer,                  optional, intent(in) :: haloshift !< The width of halos to check (default 0)
  logical,                  optional, intent(in) :: symmetric !< If true, do the checksums on the full
                                                              !! symmetric computational domain.
  logical,                  optional, intent(in) :: omit_corners !< If true, avoid checking diagonal shifts
  real,                     optional, intent(in) :: scale     !< A factor to convert this array back to unscaled units
                                                              !! for checksums and output [a A-1 ~> 1]
  integer,                  optional, intent(in) :: logunit   !< IO unit for checksum logging
  real,                     optional, intent(in) :: unscale   !< A factor to convert this array back to unscaled units
                                                              !! for checksums and output [a A-1 ~> 1].
                                                              !! Here scale and unscale are synonymous, but unscale
                                                              !! takes precedence if both are present.

  ! Local variables
  ! In the following comments, [A] is used to indicate the arbitrary, possibly rescaled units
  ! of the input array while [a] indicates the unscaled (e.g., mks) units that should be used
  ! for checksums and output

  ! Rotate array to the input grid
end subroutine chksum_B_3d
module subroutine chksum_u_3d(array_m, mesg, HI_m, haloshift, symmetric, omit_corners, &
                       scale, logunit, unscale)
  type(hor_index_type),    target,   intent(in) :: HI_m !< A horizontal index type
  real, dimension(HI_m%isdB:,HI_m%Jsd:,:), target, intent(in) :: array_m !< The array to be checksummed in
                                                             !! arbitrary, possibly rescaled units [A ~> a]
  character(len=*),                  intent(in) :: mesg      !< An identifying message
  integer,                 optional, intent(in) :: haloshift !< The width of halos to check (default 0)
  logical,                 optional, intent(in) :: symmetric !< If true, do the checksums on the full
                                                             !! symmetric computational domain.
  logical,                 optional, intent(in) :: omit_corners !< If true, avoid checking diagonal shifts
  real,                    optional, intent(in) :: scale     !< A factor to convert this array back to unscaled units
                                                             !! for checksums and output [a A-1 ~> 1]
  integer,                 optional, intent(in) :: logunit   !< IO unit for checksum logging
  real,                    optional, intent(in) :: unscale   !< A factor to convert this array back to unscaled units
                                                             !! for checksums and output [a A-1 ~> 1].
                                                             !! Here scale and unscale are synonymous, but unscale
                                                             !! takes precedence if both are present.

  ! Local variables
  ! In the following comments, [A] is used to indicate the arbitrary, possibly rescaled units
  ! of the input array while [a] indicates the unscaled (e.g., mks) units that should be used
  ! for checksums and output

  ! Rotate array to the input grid
end subroutine chksum_u_3d
module subroutine chksum_v_3d(array_m, mesg, HI_m, haloshift, symmetric, omit_corners, &
                       scale, logunit, unscale)
  type(hor_index_type),    target,   intent(in) :: HI_m      !< A horizontal index type
  real, dimension(HI_m%isd:,HI_m%JsdB:,:), target, intent(in) :: array_m !< The array to be checksummed in
                                                             !! arbitrary, possibly rescaled units [A ~> a]
  character(len=*),                  intent(in) :: mesg      !< An identifying message
  integer,                 optional, intent(in) :: haloshift !< The width of halos to check (default 0)
  logical,                 optional, intent(in) :: symmetric !< If true, do the checksums on the full
                                                             !! symmetric computational domain.
  logical,                 optional, intent(in) :: omit_corners !< If true, avoid checking diagonal shifts
  real,                    optional, intent(in) :: scale     !< A factor to convert this array back to unscaled units
                                                             !! for checksums and output [a A-1 ~> 1]
  integer,                 optional, intent(in) :: logunit   !< IO unit for checksum logging
  real,                    optional, intent(in) :: unscale   !< A factor to convert this array back to unscaled units
                                                             !! for checksums and output [a A-1 ~> 1].
                                                             !! Here scale and unscale are synonymous, but unscale
                                                             !! takes precedence if both are present.

  ! Local variables
  ! In the following comments, [A] is used to indicate the arbitrary, possibly rescaled units
  ! of the input array while [a] indicates the unscaled (e.g., mks) units that should be used
  ! for checksums and output

  ! Rotate array to the input grid
end subroutine chksum_v_3d
module subroutine chksum1d(array, mesg, start_i, end_i, compare_PEs, logunit)
  real, dimension(:), intent(in) :: array   !< The array to be summed (index starts at 1) in arbitrary units [A].
  character(len=*),   intent(in) :: mesg    !< An identifying message.
  integer, optional,  intent(in) :: start_i !< The starting index for the sum (default 1)
  integer, optional,  intent(in) :: end_i   !< The ending index for the sum (default all)
  logical, optional,  intent(in) :: compare_PEs !< If true, compare across PEs instead of summing
                                                !! and list the root_PE value (default true)
  integer, optional,  intent(in) :: logunit !< IO unit for checksum logging


end subroutine chksum1d
module subroutine chksum2d(array, mesg, logunit)

  real, dimension(:,:), intent(in) :: array !< The array to be checksummed in arbitrary units [A]
  character(len=*),     intent(in) :: mesg  !< An identifying message
  integer,    optional, intent(in) :: logunit !< IO unit for checksum logging


end subroutine chksum2d
module subroutine chksum3d(array, mesg, logunit)

  real, dimension(:,:,:), intent(in) :: array !< The array to be checksummed in arbitrary units [A]
  character(len=*),       intent(in) :: mesg  !< An identifying message
  integer,      optional, intent(in) :: logunit !< IO unit for checksum logging


end subroutine chksum3d
module function is_NaN_0d(x)
  real, intent(in) :: x !< The value to be checked for NaNs in arbitrary units [A]
  logical :: is_NaN_0d

 !is_NaN_0d = (((x < 0.0) .and. (x >= 0.0)) .or. &
 !          (.not.(x < 0.0) .and. .not.(x >= 0.0)))
end function is_NaN_0d
module function is_NaN_1d(x, skip_mpp)
  real, dimension(:), intent(in) :: x !< The array to be checked for NaNs in arbitrary units [A]
  logical,  optional, intent(in) :: skip_mpp  !< If true, only check this array only
                                              !! on the local PE (default false).
  logical :: is_NaN_1d


end function is_NaN_1d
module function is_NaN_2d(x)
  real, dimension(:,:), intent(in) :: x !< The array to be checked for NaNs in arbitrary units [A]
  logical :: is_NaN_2d


end function is_NaN_2d
module function is_NaN_3d(x)
  real, dimension(:,:,:), intent(in) :: x !< The array to be checked for NaNs in arbitrary units [A]
  logical :: is_NaN_3d


end function is_NaN_3d
module function field_checksum_real_0d(field, pelist, mask_val, turns, unscale) &
    result(chksum)
  real,              intent(in) :: field      !< Input scalar to be checksummed in arbitrary,
                                              !! possibly rescaled units [A ~> a]
  integer, optional, intent(in) :: pelist(:)  !< PE list of ranks to checksum
  real,    optional, intent(in) :: mask_val   !< FMS mask value [nondim]
  integer, optional, intent(in) :: turns      !< Number of quarter turns
  real,    optional, intent(in) :: unscale    !< A factor to convert this array back to
                                              !! unscaled units for checksums [a A-1 ~> 1]
  integer(kind=int64) :: chksum               !< checksum of scalar


end function field_checksum_real_0d
module function field_checksum_real_1d(field, pelist, mask_val, turns, unscale) &
    result(chksum)
  real, dimension(:), intent(in) :: field     !< Input array to be checksummed in arbitrary,
                                              !! possibly rescaled units [A ~> a]
  integer,  optional, intent(in) :: pelist(:) !< PE list of ranks to checksum
  real,     optional, intent(in) :: mask_val  !< FMS mask value [nondim]
  integer,  optional, intent(in) :: turns     !< Number of quarter turns
  real,     optional, intent(in) :: unscale   !< A factor to convert this array back to
                                              !! unscaled units for checksums [a A-1 ~> 1]
  integer(kind=int64) :: chksum               !< checksum of array


end function field_checksum_real_1d
module function field_checksum_real_2d(field, pelist, mask_val, turns, unscale) &
    result(chksum)
  real, dimension(:,:),     intent(in) :: field     !< Unrotated input field to be checksummed in
                                                    !! arbitrary, possibly rescaled units [A ~> a]
  integer,        optional, intent(in) :: pelist(:) !< PE list of ranks to checksum
  real,           optional, intent(in) :: mask_val  !< FMS mask value [nondim]
  integer,        optional, intent(in) :: turns     !< Number of quarter turns
  real,           optional, intent(in) :: unscale   !< A factor to convert this array back to
                                                    !! unscaled units for checksums [a A-1 ~> 1]
  integer(kind=int64) :: chksum                     !< checksum of array

  ! Local variables

end function field_checksum_real_2d
module function field_checksum_real_3d(field, pelist, mask_val, turns, unscale) &
    result(chksum)
  real, dimension(:,:,:),   intent(in) :: field     !< Unrotated input field to be checksummed in
                                                    !! arbitrary, possibly rescaled units [A ~> a]
  integer,        optional, intent(in) :: pelist(:) !< PE list of ranks to checksum
  real,           optional, intent(in) :: mask_val  !< FMS mask value [nondim]
  integer,        optional, intent(in) :: turns     !< Number of quarter turns
  real,           optional, intent(in) :: unscale   !< A factor to convert this array back to
                                                    !! unscaled units for checksums [a A-1 ~> 1]
  integer(kind=int64) :: chksum                     !< checksum of array

  ! Local variables

end function field_checksum_real_3d
module function field_checksum_real_4d(field, pelist, mask_val, turns, unscale) &
    result(chksum)
  real, dimension(:,:,:,:), intent(in) :: field     !< Unrotated input field to be checksummed in
                                                    !! arbitrary, possibly rescaled units [A ~> a]
  integer,        optional, intent(in) :: pelist(:) !< PE list of ranks to checksum
  real,           optional, intent(in) :: mask_val  !< FMS mask value [nondim]
  integer,        optional, intent(in) :: turns     !< Number of quarter turns
  real,           optional, intent(in) :: unscale   !< A factor to convert this array back to
                                                    !! unscaled units for checksums [a A-1 ~> 1]
  integer(kind=int64) :: chksum                     !< checksum of array

  ! Local variables

end function field_checksum_real_4d
module subroutine chk_sum_msg1(fmsg, bc0, mesg, iounit)
  character(len=*), intent(in) :: fmsg !< A checksum code-location specific preamble
  character(len=*), intent(in) :: mesg !< An identifying message supplied by top-level caller
  integer,          intent(in) :: bc0  !< The bitcount of the non-shifted array
  integer,          intent(in) :: iounit !< Checksum logger IO unit

end subroutine chk_sum_msg1
module subroutine chk_sum_msg5(fmsg, bc0, bcSW, bcSE, bcNW, bcNE, mesg, iounit)
  character(len=*), intent(in) :: fmsg !< A checksum code-location specific preamble
  character(len=*), intent(in) :: mesg !< An identifying message supplied by top-level caller
  integer,          intent(in) :: bc0  !< The bitcount of the non-shifted array
  integer,          intent(in) :: bcSW !< The bitcount for SW shifted array
  integer,          intent(in) :: bcSE !< The bitcount for SE shifted array
  integer,          intent(in) :: bcNW !< The bitcount for NW shifted array
  integer,          intent(in) :: bcNE !< The bitcount for NE shifted array
  integer,          intent(in) :: iounit !< Checksum logger IO unit

end subroutine chk_sum_msg5
module subroutine chk_sum_msg_NSEW(fmsg, bc0, bcN, bcS, bcE, bcW, mesg, iounit)
  character(len=*), intent(in) :: fmsg !< A checksum code-location specific preamble
  character(len=*), intent(in) :: mesg !< An identifying message supplied by top-level caller
  integer,          intent(in) :: bc0  !< The bitcount of the non-shifted array
  integer,          intent(in) :: bcN !< The bitcount for N shifted array
  integer,          intent(in) :: bcS !< The bitcount for S shifted array
  integer,          intent(in) :: bcE !< The bitcount for E shifted array
  integer,          intent(in) :: bcW !< The bitcount for W shifted array
  integer,          intent(in) :: iounit !< Checksum logger IO unit

end subroutine chk_sum_msg_NSEW
module subroutine chk_sum_msg_S(fmsg, bc0, bcS, mesg, iounit)
  character(len=*), intent(in) :: fmsg !< A checksum code-location specific preamble
  character(len=*), intent(in) :: mesg !< An identifying message supplied by top-level caller
  integer,          intent(in) :: bc0  !< The bitcount of the non-shifted array
  integer,          intent(in) :: bcS  !< The bitcount of the south-shifted array
  integer,          intent(in) :: iounit !< Checksum logger IO unit

end subroutine chk_sum_msg_S
module subroutine chk_sum_msg_W(fmsg, bc0, bcW, mesg, iounit)
  character(len=*), intent(in) :: fmsg !< A checksum code-location specific preamble
  character(len=*), intent(in) :: mesg !< An identifying message supplied by top-level caller
  integer,          intent(in) :: bc0  !< The bitcount of the non-shifted array
  integer,          intent(in) :: bcW  !< The bitcount of the west-shifted array
  integer,          intent(in) :: iounit !< Checksum logger IO unit

end subroutine chk_sum_msg_W
module subroutine chk_sum_msg2(fmsg, bc0, bcSW, mesg, iounit)
  character(len=*), intent(in) :: fmsg !< A checksum code-location specific preamble
  character(len=*), intent(in) :: mesg !< An identifying message supplied by top-level caller
  integer,          intent(in) :: bc0  !< The bitcount of the non-shifted array
  integer,          intent(in) :: bcSW !< The bitcount of the southwest-shifted array
  integer,          intent(in) :: iounit !< Checksum logger IO unit

end subroutine chk_sum_msg2
module subroutine chk_sum_msg3(fmsg, aMean, aMin, aMax, mesg, iounit)
  character(len=*), intent(in) :: fmsg !< A checksum code-location specific preamble
  character(len=*), intent(in) :: mesg !< An identifying message supplied by top-level caller
  real,             intent(in) :: aMean !< The mean value of the array in arbitrary units [A]
  real,             intent(in) :: aMin !< The minimum value of the array [A]
  real,             intent(in) :: aMax !< The maximum value of the array [A]
  integer,          intent(in) :: iounit !< Checksum logger IO unit

  ! NOTE: We add zero to aMin and aMax to remove any negative zeros.
  ! This is due to inconsistencies of signed zero in local vs MPI calculations.

end subroutine chk_sum_msg3
module subroutine MOM_checksums_init(param_file)
  type(param_file_type),   intent(in)    :: param_file !< A structure to parse for run-time parameters
  ! This include declares and sets the variable "version".

end subroutine MOM_checksums_init
module subroutine chksum_error(signal, message)
  ! Wrapper for MOM_error to help place specific break points in debuggers
  integer, intent(in) :: signal !< An error severity level, such as FATAL or WARNING
  character(len=*), intent(in) :: message !< An error message
end subroutine chksum_error
integer module function bitcount(x)
  real, intent(in) :: x !< Number to be bitcount in arbitrary units [A]


  ! NOTE: Assumes that reals and integers of kind=xk are the same size
end function bitcount
  end interface

end module MOM_checksums
