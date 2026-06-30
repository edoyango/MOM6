! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> Thin interfaces to non-domain-oriented mpp communication subroutines
module MOM_coms_infra

use iso_fortran_env, only : int32, int64

use mpp_mod, only : mpp_pe, mpp_root_pe, mpp_npes, mpp_set_root_pe
use mpp_mod, only : mpp_set_current_pelist, mpp_get_current_pelist
use mpp_mod, only : mpp_broadcast, mpp_sync, mpp_sync_self, mpp_chksum
use mpp_mod, only : mpp_sum, mpp_max, mpp_min
use memutils_mod, only : print_memuse_stats
use fms_mod, only : fms_end, fms_init

implicit none ; private

public :: PE_here, root_PE, num_PEs, set_rootPE, Set_PElist, Get_PElist, sync_PEs
public :: broadcast, sum_across_PEs, min_across_PEs, max_across_PEs
public :: any_across_PEs, all_across_PEs
public :: field_chksum, MOM_infra_init, MOM_infra_end

! This module provides interfaces to the non-domain-oriented communication
! subroutines.

!> Communicate an array, string or scalar from one PE to others
interface broadcast
  module procedure broadcast_char, broadcast_int32_0D, broadcast_int64_0D, broadcast_int1D
  module procedure broadcast_real0D, broadcast_real1D, broadcast_real2D, broadcast_real3D
end interface broadcast

!> Compute a checksum for a field distributed over a PE list.  If no PE list is
!! provided, then the current active PE list is used.
interface field_chksum
  module procedure field_chksum_real_0d
  module procedure field_chksum_real_1d
  module procedure field_chksum_real_2d
  module procedure field_chksum_real_3d
  module procedure field_chksum_real_4d
end interface field_chksum

!> Find the sum of field across PEs, and update PEs with the sums.
interface sum_across_PEs
  module procedure sum_across_PEs_int4_0d
  module procedure sum_across_PEs_int4_1d
  module procedure sum_across_PEs_int4_2d
  module procedure sum_across_PEs_int8_0d
  module procedure sum_across_PEs_int8_1d
  module procedure sum_across_PEs_int8_2d
  module procedure sum_across_PEs_real_0d
  module procedure sum_across_PEs_real_1d
  module procedure sum_across_PEs_real_2d
end interface sum_across_PEs

!> Find the maximum value of field across PEs, and update PEs with the values.
interface max_across_PEs
  module procedure max_across_PEs_int_0d
  module procedure max_across_PEs_real_0d
  module procedure max_across_PEs_real_1d
end interface max_across_PEs

!> Find the minimum value of field across PEs, and update PEs with the values.
interface min_across_PEs
  module procedure min_across_PEs_int_0d
  module procedure min_across_PEs_real_0d
  module procedure min_across_PEs_real_1d
end interface min_across_PEs


  interface
module function PE_here() result(pe)
  integer :: pe   !< PE ID of the current process
end function PE_here
module function root_PE() result(pe)
  integer :: pe   !< root PE ID
end function root_PE
module function num_PEs() result(npes)
  integer :: npes   !< Number of PEs
end function num_PEs
module subroutine set_rootPE(pe)
  integer, intent(in) :: pe   !< ID of the PE to be assigned as root
end subroutine set_rootPE
module subroutine Set_PEList(pelist, no_sync)
  integer, optional, intent(in) :: pelist(:)  !< List of PEs to set for communication
  logical, optional, intent(in) :: no_sync    !< Do not sync after list update.
end subroutine Set_PEList
module subroutine Get_PEList(pelist, name, commID)
  integer,                    intent(out) :: pelist(:) !< List of PE IDs of the current PE list
  character(len=*), optional, intent(out) :: name   !< Name of PE list
  integer,          optional, intent(out) :: commID !< Communicator ID of PE list

end subroutine Get_PEList
module subroutine sync_PEs(pelist)
  integer, optional, intent(in) :: pelist(:)  !< The list of PEs to be synced

end subroutine sync_PEs
module subroutine broadcast_char(dat, length, from_PE, PElist, blocking)
  character(len=*),  intent(inout) :: dat(:)    !< The data to communicate and destination
  integer,           intent(in)    :: length    !< The length of each string
  integer, optional, intent(in)    :: from_PE   !< The source PE, by default the root PE
  integer, optional, intent(in)    :: PElist(:) !< The list of participating PEs, by default the
                                                !! active PE set as previously set via Set_PElist.
  logical, optional, intent(in)    :: blocking  !< If true, barriers are added around the call


end subroutine broadcast_char
module subroutine broadcast_int64_0D(dat, from_PE, PElist, blocking)
  integer(kind=int64),   intent(inout) :: dat       !< The data to communicate and destination
  integer,     optional, intent(in)    :: from_PE   !< The source PE, by default the root PE
  integer,     optional, intent(in)    :: PElist(:) !< The list of participating PEs, by default the
                                                    !! active PE set as previously set via Set_PElist.
  logical,     optional, intent(in)    :: blocking  !< If true, barriers are added around the call


end subroutine broadcast_int64_0D
module subroutine broadcast_int32_0D(dat, from_PE, PElist, blocking)
  integer(kind=int32),   intent(inout) :: dat       !< The data to communicate and destination
  integer,     optional, intent(in)    :: from_PE   !< The source PE, by default the root PE
  integer,     optional, intent(in)    :: PElist(:) !< The list of participating PEs, by default the
                                                    !! active PE set as previously set via Set_PElist.
  logical,     optional, intent(in)    :: blocking  !< If true, barriers are added around the call


end subroutine broadcast_int32_0D
module subroutine broadcast_int1D(dat, length, from_PE, PElist, blocking)
  integer, dimension(:), intent(inout) :: dat       !< The data to communicate and destination
  integer,               intent(in)    :: length    !< The number of data elements
  integer,     optional, intent(in)    :: from_PE   !< The source PE, by default the root PE
  integer,     optional, intent(in)    :: PElist(:) !< The list of participating PEs, by default the
                                                    !! active PE set as previously set via Set_PElist.
  logical,     optional, intent(in)    :: blocking  !< If true, barriers are added around the call


end subroutine broadcast_int1D
module subroutine broadcast_real0D(dat, from_PE, PElist, blocking)
  real,                 intent(inout) :: dat       !< The data to communicate and destination
  integer,    optional, intent(in)    :: from_PE   !< The source PE, by default the root PE
  integer,    optional, intent(in)    :: PElist(:) !< The list of participating PEs, by default the
                                                   !! active PE set as previously set via Set_PElist.
  logical,    optional, intent(in)    :: blocking  !< If true, barriers are added around the call


end subroutine broadcast_real0D
module subroutine broadcast_real1D(dat, length, from_PE, PElist, blocking)
  real, dimension(:),   intent(inout) :: dat       !< The data to communicate and destination
  integer,              intent(in)    :: length    !< The number of data elements
  integer,    optional, intent(in)    :: from_PE   !< The source PE, by default the root PE
  integer,    optional, intent(in)    :: PElist(:) !< The list of participating PEs, by default the
                                                   !! active PE set as previously set via Set_PElist.
  logical,    optional, intent(in)    :: blocking  !< If true, barriers are added around the call


end subroutine broadcast_real1D
module subroutine broadcast_real2D(dat, length, from_PE, PElist, blocking)
  real, dimension(:,:), intent(inout) :: dat       !< The data to communicate and destination
  integer,              intent(in)    :: length    !< The total number of data elements
  integer,    optional, intent(in)    :: from_PE   !< The source PE, by default the root PE
  integer,    optional, intent(in)    :: PElist(:) !< The list of participating PEs, by default the
                                                   !! active PE set as previously set via Set_PElist.
  logical,    optional, intent(in)    :: blocking  !< If true, barriers are added around the call


end subroutine broadcast_real2D
module subroutine broadcast_real3D(dat, length, from_PE, PElist, blocking)
  real, dimension(:,:,:), intent(inout) :: dat       !< The data to communicate and destination
  integer,              intent(in)    :: length    !< The total number of data elements
  integer,    optional, intent(in)    :: from_PE   !< The source PE, by default the root PE
  integer,    optional, intent(in)    :: PElist(:) !< The list of participating PEs, by default the
                                                   !! active PE set as previously set via Set_PElist.
  logical,    optional, intent(in)    :: blocking  !< If true, barriers are added around the call


end subroutine broadcast_real3D
module function field_chksum_real_0d(field, pelist, mask_val) result(chksum)
  real,              intent(in) :: field      !< Input scalar
  integer, optional, intent(in) :: pelist(:)  !< PE list of ranks to checksum
  real,    optional, intent(in) :: mask_val   !< FMS mask value
  integer(kind=int64) :: chksum               !< checksum of array

end function field_chksum_real_0d
module function field_chksum_real_1d(field, pelist, mask_val) result(chksum)
  real, dimension(:), intent(in) :: field     !< Input array
  integer,  optional, intent(in) :: pelist(:) !< PE list of ranks to checksum
  real,     optional, intent(in) :: mask_val  !< FMS mask value
  integer(kind=int64) :: chksum               !< checksum of array

end function field_chksum_real_1d
module function field_chksum_real_2d(field, pelist, mask_val) result(chksum)
  real, dimension(:,:), intent(in) :: field     !< Unrotated input field
  integer,    optional, intent(in) :: pelist(:) !< PE list of ranks to checksum
  real,       optional, intent(in) :: mask_val  !< FMS mask value
  integer(kind=int64) :: chksum                 !< checksum of array

end function field_chksum_real_2d
module function field_chksum_real_3d(field, pelist, mask_val) result(chksum)
  real, dimension(:,:,:), intent(in) :: field     !< Unrotated input field
  integer,      optional, intent(in) :: pelist(:) !< PE list of ranks to checksum
  real,         optional, intent(in) :: mask_val  !< FMS mask value
  integer(kind=int64) :: chksum               !< checksum of array

end function field_chksum_real_3d
module function field_chksum_real_4d(field, pelist, mask_val) result(chksum)
  real, dimension(:,:,:,:), intent(in) :: field     !< Unrotated input field
  integer,        optional, intent(in) :: pelist(:) !< PE list of ranks to checksum
  real,           optional, intent(in) :: mask_val  !< FMS mask value
  integer(kind=int64) :: chksum               !< checksum of array

end function field_chksum_real_4d
module subroutine sum_across_PEs_int4_0d(field, pelist)
  integer(kind=int32), intent(inout) :: field     !< Value on this PE, and the sum across PEs upon return
  integer,   optional, intent(in)    :: pelist(:) !< List of PEs to work with

end subroutine sum_across_PEs_int4_0d
module subroutine sum_across_PEs_int4_1d(field, length, pelist)
  integer(kind=int32), dimension(:), intent(inout) :: field     !< The values to add, the sums upon return
  integer,                           intent(in)    :: length    !< Number of elements in field to add
  integer,                 optional, intent(in)    :: pelist(:) !< List of PEs to work with

end subroutine sum_across_PEs_int4_1d
module subroutine sum_across_PEs_int4_2d(field, length, pelist)
  integer(kind=int32), dimension(:,:), intent(inout) :: field     !< The values to add, the sums upon return
  integer,                             intent(in)    :: length    !< Number of elements in field to add
  integer,                 optional,   intent(in)    :: pelist(:) !< List of PEs to work with

end subroutine sum_across_PEs_int4_2d
module subroutine sum_across_PEs_int8_0d(field, pelist)
  integer(kind=int64), intent(inout) :: field     !< Value on this PE, and the sum across PEs upon return
  integer,   optional, intent(in)    :: pelist(:) !< List of PEs to work with

end subroutine sum_across_PEs_int8_0d
module subroutine sum_across_PEs_int8_1d(field, length, pelist)
  integer(kind=int64), dimension(:), intent(inout) :: field     !< The values to add, the sums upon return
  integer,                           intent(in)    :: length    !< Number of elements in field to add
  integer,                 optional, intent(in)    :: pelist(:) !< List of PEs to work with

end subroutine sum_across_PEs_int8_1d
module subroutine sum_across_PEs_int8_2d(field, length, pelist)
  integer(kind=int64), &
           dimension(:,:), intent(inout) :: field     !< The values to add, the sums upon return
  integer,                 intent(in)    :: length    !< The total number of positions to sum, usually
                                                      !! the product of the array sizes.
  integer,       optional, intent(in)    :: pelist(:) !< List of PEs to work with

end subroutine sum_across_PEs_int8_2d
module subroutine sum_across_PEs_real_0d(field, pelist)
  real,              intent(inout) :: field     !< Value on this PE, and the sum across PEs upon return
  integer, optional, intent(in)    :: pelist(:) !< List of PEs to work with

end subroutine sum_across_PEs_real_0d
module subroutine sum_across_PEs_real_1d(field, length, pelist)
  real, dimension(:), intent(inout) :: field     !< The values to add, the sums upon return
  integer,            intent(in)    :: length    !< Number of elements in field to add
  integer,  optional, intent(in)    :: pelist(:) !< List of PEs to work with

end subroutine sum_across_PEs_real_1d
module subroutine sum_across_PEs_real_2d(field, length, pelist)
  real, dimension(:,:), intent(inout) :: field     !< The values to add, the sums upon return
  integer,              intent(in)    :: length    !< The total number of positions to sum, usually
                                                   !! the product of the array sizes.
  integer,    optional, intent(in)    :: pelist(:) !< List of PEs to work with

end subroutine sum_across_PEs_real_2d
module subroutine max_across_PEs_int_0d(field, pelist)
  integer,           intent(inout) :: field     !< The values to compare, the maximum upon return
  integer, optional, intent(in)    :: pelist(:) !< List of PEs to work with

end subroutine max_across_PEs_int_0d
module subroutine max_across_PEs_real_0d(field, pelist)
  real,              intent(inout) :: field     !< The values to compare, the maximum upon return
  integer, optional, intent(in)    :: pelist(:) !< List of PEs to work with

end subroutine max_across_PEs_real_0d
module subroutine max_across_PEs_real_1d(field, length, pelist)
  real, dimension(:), intent(inout) :: field     !< The list of values being compared, with the
                                                 !! maxima in each position upon return
  integer,            intent(in)    :: length    !< Number of elements in field to compare
  integer,  optional, intent(in)    :: pelist(:) !< List of PEs to work with

end subroutine max_across_PEs_real_1d
module subroutine min_across_PEs_int_0d(field, pelist)
  integer,           intent(inout) :: field     !< The values to compare, the minimum upon return
  integer, optional, intent(in)    :: pelist(:) !< List of PEs to work with

end subroutine min_across_PEs_int_0d
module subroutine min_across_PEs_real_0d(field, pelist)
  real,              intent(inout) :: field     !< The values to compare, the minimum upon return
  integer, optional, intent(in)    :: pelist(:) !< List of PEs to work with
end subroutine min_across_PEs_real_0d
module subroutine min_across_PEs_real_1d(field, length, pelist)
  real, dimension(:), intent(inout) :: field     !< The list of values being compared, with the
                                                 !! minima in each position upon return
  integer,            intent(in)    :: length    !< Number of elements in field to compare
  integer,  optional, intent(in)    :: pelist(:) !< List of PEs to work with

end subroutine min_across_PEs_real_1d
module function any_across_PEs(field, pelist)
  logical, intent(in)           :: field      !< Local PE value
  integer, optional, intent(in) :: pelist(:)  !< List of PEs to work with
  logical :: any_across_PEs


  ! FMS1 does not support logical collectives, so integer flags are used.
end function any_across_PEs
module function all_across_PEs(field, pelist)
  logical, intent(in)           :: field      !< Local PE value
  integer, optional, intent(in) :: pelist(:)  !< List of PEs to work with
  logical :: all_across_PEs


  ! FMS1 does not support logical collectives, so integer flags are used.
end function all_across_PEs
module subroutine MOM_infra_init(localcomm)
  integer, optional, intent(in) :: localcomm  !< Communicator ID to initialize
end subroutine MOM_infra_init
module subroutine MOM_infra_end
end subroutine MOM_infra_end
  end interface

end module MOM_coms_infra
