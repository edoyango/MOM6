! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> Describes the decomposed MOM domain and has routines for communications across PEs
module MOM_domain_infra

use MOM_coms_infra,  only : PE_here, root_PE, num_PEs
use MOM_cpu_clock_infra, only : cpu_clock_begin, cpu_clock_end
use MOM_error_infra, only : MOM_error=>MOM_err, NOTE, WARNING, FATAL

use mpp_domains_mod, only : domain2D, domain1D
use mpp_domains_mod, only : mpp_define_io_domain, mpp_define_domains, mpp_deallocate_domain
use mpp_domains_mod, only : mpp_get_domain_components, mpp_get_domain_extents, mpp_get_layout
use mpp_domains_mod, only : mpp_get_compute_domain, mpp_get_data_domain, mpp_get_global_domain
use mpp_domains_mod, only : mpp_get_boundary, mpp_update_domains
use mpp_domains_mod, only : mpp_start_update_domains, mpp_complete_update_domains
use mpp_domains_mod, only : mpp_create_group_update, mpp_do_group_update
use mpp_domains_mod, only : mpp_reset_group_update_field, mpp_group_update_initialized
use mpp_domains_mod, only : mpp_start_group_update, mpp_complete_group_update
use mpp_domains_mod, only : mpp_compute_block_extent, mpp_compute_extent
use mpp_domains_mod, only : mpp_broadcast_domain, mpp_redistribute, mpp_global_field
use mpp_domains_mod, only : AGRID, BGRID_NE, CGRID_NE, SCALAR_PAIR, BITWISE_EXACT_SUM
use mpp_domains_mod, only : CYCLIC_GLOBAL_DOMAIN
use mpp_domains_mod, only : FOLD_NORTH_EDGE, FOLD_SOUTH_EDGE, FOLD_EAST_EDGE, FOLD_WEST_EDGE
use mpp_domains_mod, only : To_East => WUPDATE, To_West => EUPDATE, Omit_Corners => EDGEUPDATE
use mpp_domains_mod, only : To_North => SUPDATE, To_South => NUPDATE
use mpp_domains_mod, only : CENTER, CORNER, NORTH_FACE => NORTH, EAST_FACE => EAST
use fms_io_utils_mod, only : file_exists, parse_mask_table
use fms_affinity_mod, only : fms_affinity_init, fms_affinity_set, fms_affinity_get

! This subroutine is not in MOM6/src but may be required by legacy drivers
! use mpp_domains_mod, only : global_field_sum => mpp_global_sum

! The `group_pass_type` fields are never accessed, so we keep it as an FMS type
use mpp_domains_mod, only : group_pass_type => mpp_group_update_type

implicit none ; private

! These types are inherited from mpp, but are treated as opaque here.
public :: domain2D, domain1D, group_pass_type
! These interfaces are actually implemented or have explicit interfaces in this file.
public :: create_MOM_domain, clone_MOM_domain, get_domain_components, get_domain_extent
public :: deallocate_MOM_domain, get_global_shape, compute_block_extent, compute_extent
public :: pass_var, pass_vector, fill_symmetric_edges, rescale_comp_data
public :: pass_var_start, pass_var_complete, pass_vector_start, pass_vector_complete
public :: create_group_pass, do_group_pass, start_group_pass, complete_group_pass
public :: redistribute_array, broadcast_domain, same_domain, global_field
public :: get_simple_array_i_ind, get_simple_array_j_ind
public :: MOM_thread_affinity_set, set_MOM_thread_affinity
! These are encoding constant parmeters with self-explanatory names.
public :: To_East, To_West, To_North, To_South, To_All, Omit_Corners
public :: AGRID, BGRID_NE, CGRID_NE, SCALAR_PAIR
public :: CORNER, CENTER, NORTH_FACE, EAST_FACE
public :: set_domain, nullify_domain
! These are no longer used by MOM6 because the reproducing sum works so well, but they are
! still referenced by some of the non-GFDL couplers.
! public :: global_field_sum, BITWISE_EXACT_SUM

!> Do a halo update on an array
interface pass_var
  module procedure pass_var_3d, pass_var_2d
end interface pass_var

!> Do a halo update on a pair of arrays representing the two components of a vector
interface pass_vector
  module procedure pass_vector_3d, pass_vector_2d
end interface pass_vector

!> Initiate a non-blocking halo update on an array
interface pass_var_start
  module procedure pass_var_start_3d, pass_var_start_2d
end interface pass_var_start

!> Complete a non-blocking halo update on an array
interface pass_var_complete
  module procedure pass_var_complete_3d, pass_var_complete_2d
end interface pass_var_complete

!> Initiate a halo update on a pair of arrays representing the two components of a vector
interface pass_vector_start
  module procedure pass_vector_start_3d, pass_vector_start_2d
end interface pass_vector_start

!> Complete a halo update on a pair of arrays representing the two components of a vector
interface pass_vector_complete
  module procedure pass_vector_complete_3d, pass_vector_complete_2d
end interface pass_vector_complete

!> Set up a group of halo updates
interface create_group_pass
  module procedure create_var_group_pass_2d
  module procedure create_var_group_pass_3d
  module procedure create_vector_group_pass_2d
  module procedure create_vector_group_pass_3d
end interface create_group_pass

!> Do a set of halo updates that fill in the values at the duplicated edges
!! of a staggered symmetric memory domain
interface fill_symmetric_edges
  module procedure fill_vector_symmetric_edges_2d !, fill_vector_symmetric_edges_3d
!   module procedure fill_scalar_symmetric_edges_2d, fill_scalar_symmetric_edges_3d
end interface fill_symmetric_edges

!> Rescale the values of an array in its computational domain by a constant factor
interface rescale_comp_data
  module procedure rescale_comp_data_4d, rescale_comp_data_3d, rescale_comp_data_2d
end interface rescale_comp_data

!> Pass an array from one MOM domain to another
interface redistribute_array
  module procedure redistribute_array_2d, redistribute_array_3d, redistribute_array_4d
end interface redistribute_array

!> Copy one MOM_domain_type into another
interface clone_MOM_domain
  module procedure clone_MD_to_MD, clone_MD_to_d2D
end interface clone_MOM_domain

!> Extract the 1-d domain components from a MOM_domain or domain2d
interface get_domain_components
  module procedure get_domain_components_MD, get_domain_components_d2D
end interface get_domain_components

!> Returns the index ranges that have been stored in a MOM_domain_type
interface get_domain_extent
  module procedure get_domain_extent_MD, get_domain_extent_d2D
end interface get_domain_extent


!> The MOM_domain_type contains information about the domain decomposition.
type, public :: MOM_domain_type
  character(len=64) :: name     !< The name of this domain
  type(domain2D), pointer :: mpp_domain => NULL() !< The FMS domain with halos
                                !! on this processor, centered at h points.
  type(domain2D), pointer :: mpp_domain_d(:) => NULL() !< A coarse FMS domain with halos
                                !! on this processor, centered at h points.
  integer :: niglobal           !< The total horizontal i-domain size.
  integer :: njglobal           !< The total horizontal j-domain size.
  integer :: nihalo             !< The i-halo size in memory.
  integer :: njhalo             !< The j-halo size in memory.
  logical :: symmetric          !< True if symmetric memory is used with this domain.
  logical :: nonblocking_updates  !< If true, non-blocking halo updates are
                                !! allowed.  The default is .false. (for now).
  logical :: thin_halo_updates  !< If true, optional arguments may be used to
                                !! specify the width of the halos that are
                                !! updated with each call.
  integer :: layout(2)          !< This domain's processor layout.  This is
                                !! saved to enable the construction of related
                                !! new domains with different resolutions or
                                !! other properties.
  integer :: io_layout(2)       !< The IO-layout used with this domain.
  integer :: X_FLAGS            !< Flag that specifies the properties of the
                                !! domain in the i-direction in a define_domain call.
  integer :: Y_FLAGS            !< Flag that specifies the properties of the
                                !! domain in the j-direction in a define_domain call.
  logical, pointer :: maskmap(:,:) => NULL() !< A pointer to an array indicating
                                !! which logical processors are actually used for
                                !! the ocean code. The other logical processors
                                !! would be contain only land points and are not
                                !! assigned to actual processors. This need not be
                                !! assigned if all logical processors are used.
  integer :: turns              !< Number of quarter-turns from input to this grid.
  type(MOM_domain_type), pointer :: domain_in => NULL()
                                !< Reference to unrotated domain (if turned)
end type MOM_domain_type

integer, parameter :: To_All = To_East + To_West + To_North + To_South !< A flag for passing in all directions


  interface
module subroutine pass_var_3d(array, MOM_dom, sideflag, complete, position, halo, &
                       clock)
  real, dimension(:,:,:), intent(inout) :: array    !< The array which is having its halos points
                                                    !! exchanged.
  type(MOM_domain_type),  intent(inout) :: MOM_dom  !< The MOM_domain_type containing the mpp_domain
                                                    !! needed to determine where data should be
                                                    !! sent.
  integer,      optional, intent(in)    :: sideflag !< An optional integer indicating which
      !! directions the data should be sent.  It is TO_ALL or the sum of any of TO_EAST, TO_WEST,
      !! TO_NORTH, and TO_SOUTH.  For example, TO_EAST sends the data to the processor to the east,
      !! sothe halos on the western side are filled.  TO_ALL is the default if sideflag is omitted.
  logical,      optional, intent(in)    :: complete !< An optional argument indicating whether the
                                                    !! halo updates should be completed before
                                                    !! progress resumes. Omitting complete is the
                                                    !! same as setting complete to .true.
  integer,      optional, intent(in)    :: position !< An optional argument indicating the position.
                                                    !! This is CENTER by default and is often CORNER,
                                                    !! but could also be EAST_FACE or NORTH_FACE.
  integer,      optional, intent(in)    :: halo     !< The size of the halo to update - the full
                                                    !! halo by default.
  integer,      optional, intent(in)    :: clock    !< The handle for a cpu time clock that should be
                                                    !! started then stopped to time this routine.


end subroutine pass_var_3d
module subroutine pass_var_2d(array, MOM_dom, sideflag, complete, position, halo, inner_halo, clock)
  real, dimension(:,:),  intent(inout) :: array    !< The array which is having its halos points
                                                   !! exchanged.
  type(MOM_domain_type), intent(inout) :: MOM_dom  !< The MOM_domain_type containing the mpp_domain
                                                   !! needed to determine where data should be sent.
  integer,     optional, intent(in)    :: sideflag !< An optional integer indicating which
      !! directions the data should be sent. It is TO_ALL or the sum of any of TO_EAST, TO_WEST,
      !! TO_NORTH, and TO_SOUTH.  For example, TO_EAST sends the data to the processor to the east,
      !! so the halos on the western side are filled.  TO_ALL is the default if sideflag is omitted.
  logical,     optional, intent(in)    :: complete !< An optional argument indicating whether the
                                                   !! halo updates should be completed before
                                                   !! progress resumes.  Omitting complete is the
                                                   !! same as setting complete to .true.
  integer,     optional, intent(in)    :: position !< An optional argument indicating the position.
                                                   !! This is CENTER by default and is often CORNER,
                                                   !! but could also be EAST_FACE or NORTH_FACE.
  integer,     optional, intent(in)    :: halo     !< The size of the halo to update - the full halo
                                                   !! by default.
  integer,     optional, intent(in)    :: inner_halo !< The size of an inner halo to avoid updating,
                                                   !! or 0 to avoid updating symmetric memory
                                                   !! computational domain points.  Setting this >=0
                                                   !! also enforces that complete=.true.
  integer,     optional, intent(in)    :: clock    !< The handle for a cpu time clock that should be
                                                   !! started then stopped to time this routine.

  ! Local variables

end subroutine pass_var_2d
module function pass_var_start_2d(array, MOM_dom, sideflag, position, complete, halo, &
                           clock)
  real, dimension(:,:),   intent(inout) :: array    !< The array which is having its halos points
                                                    !! exchanged.
  type(MOM_domain_type),  intent(inout) :: MOM_dom  !< The MOM_domain_type containing the mpp_domain
                                                    !! needed to determine where data should be
                                                    !! sent.
  integer,      optional, intent(in)    :: sideflag !< An optional integer indicating which
      !! directions the data should be sent. It is TO_ALL or the sum of any of TO_EAST, TO_WEST,
      !! TO_NORTH, and TO_SOUTH.  For example, TO_EAST sends the data to the processor to the east,
      !! so the halos on the western side are filled.  TO_ALL is the default if sideflag is omitted.
  integer,      optional, intent(in)    :: position !< An optional argument indicating the position.
                                                    !! This is CENTER by default and is often CORNER,
                                                    !! but could also be EAST_FACE or NORTH_FACE.
  logical,      optional, intent(in)    :: complete !< An optional argument indicating whether the
                                                    !! halo updates should be completed before
                                                    !! progress resumes.  Omitting complete is the
                                                    !! same as setting complete to .true.
  integer,      optional, intent(in)    :: halo     !< The size of the halo to update - the full
                                                    !! halo by default.
  integer,      optional, intent(in)    :: clock    !< The handle for a cpu time clock that should be
                                                    !! started then stopped to time this routine.
  integer                               :: pass_var_start_2d  !<The integer index for this update.


end function pass_var_start_2d
module function pass_var_start_3d(array, MOM_dom, sideflag, position, complete, halo, &
                           clock)
  real, dimension(:,:,:), intent(inout) :: array    !< The array which is having its halos points
                                                    !! exchanged.
  type(MOM_domain_type),  intent(inout) :: MOM_dom  !< The MOM_domain_type containing the mpp_domain
                                                    !! needed to determine where data should be
                                                    !! sent.
  integer,      optional, intent(in)    :: sideflag !< An optional integer indicating which
      !! directions the data should be sent. It is TO_ALL or the sum of any of TO_EAST, TO_WEST,
      !! TO_NORTH, and TO_SOUTH.  For example, TO_EAST sends the data to the processor to the east,
      !! so the halos on the western side are filled.  TO_ALL is the default if sideflag is omitted.
  integer,      optional, intent(in)    :: position !< An optional argument indicating the position.
                                                    !! This is CENTER by default and is often CORNER,
                                                    !! but could also be EAST_FACE or NORTH_FACE.
  logical,      optional, intent(in)    :: complete !< An optional argument indicating whether the
                                                    !! halo updates should be completed before
                                                    !! progress resumes.  Omitting complete is the
                                                    !! same as setting complete to .true.
  integer,      optional, intent(in)    :: halo     !< The size of the halo to update - the full
                                                    !! halo by default.
  integer,      optional, intent(in)    :: clock    !< The handle for a cpu time clock that should be
                                                    !! started then stopped to time this routine.
  integer                               :: pass_var_start_3d  !< The integer index for this update.


end function pass_var_start_3d
module subroutine pass_var_complete_2d(id_update, array, MOM_dom, sideflag, position, halo, &
                                clock)
  integer,                intent(in)    :: id_update !< The integer id of this update which has
                                                    !! been returned from a previous call to
                                                    !! pass_var_start.
  real, dimension(:,:),   intent(inout) :: array    !< The array which is having its halos points
                                                    !! exchanged.
  type(MOM_domain_type),  intent(inout) :: MOM_dom  !< The MOM_domain_type containing the mpp_domain
                                                    !! needed to determine where data should be
                                                    !! sent.
  integer,      optional, intent(in)    :: sideflag !< An optional integer indicating which
      !! directions the data should be sent. It is TO_ALL or the sum of any of TO_EAST, TO_WEST,
      !! TO_NORTH, and TO_SOUTH.  For example, TO_EAST sends the data to the processor to the east,
      !! so the halos on the western side are filled.  TO_ALL is the default if sideflag is omitted.
  integer,      optional, intent(in)    :: position !< An optional argument indicating the position.
                                                    !! This is CENTER by default and is often CORNER,
                                                    !! but could also be EAST_FACE or NORTH_FACE.
  integer,      optional, intent(in)    :: halo     !< The size of the halo to update - the full
                                                    !! halo by default.
  integer,      optional, intent(in)    :: clock    !< The handle for a cpu time clock that should be
                                                    !! started then stopped to time this routine.


end subroutine pass_var_complete_2d
module subroutine pass_var_complete_3d(id_update, array, MOM_dom, sideflag, position, halo, &
                                clock)
  integer,                intent(in)    :: id_update !< The integer id of this update which has
                                                    !! been returned from a previous call to
                                                    !! pass_var_start.
  real, dimension(:,:,:), intent(inout) :: array    !< The array which is having its halos points
                                                    !! exchanged.
  type(MOM_domain_type),  intent(inout) :: MOM_dom  !< The MOM_domain_type containing the mpp_domain
                                                    !! needed to determine where data should be
                                                    !! sent.
  integer,      optional, intent(in)    :: sideflag !< An optional integer indicating which
      !! directions the data should be sent. It is TO_ALL or the sum of any of TO_EAST, TO_WEST,
      !! TO_NORTH, and TO_SOUTH.  For example, TO_EAST sends the data to the processor to the east,
      !! so the halos on the western side are filled.  TO_ALL is the default if sideflag is omitted.
  integer,      optional, intent(in)    :: position !< An optional argument indicating the position.
                                                    !! This is CENTER by default and is often CORNER,
                                                    !! but could also be EAST_FACE or NORTH_FACE.
  integer,      optional, intent(in)    :: halo     !< The size of the halo to update - the full
                                                    !! halo by default.
  integer,      optional, intent(in)    :: clock    !< The handle for a cpu time clock that should be
                                                    !! started then stopped to time this routine.


end subroutine pass_var_complete_3d
module subroutine pass_vector_2d(u_cmpt, v_cmpt, MOM_dom, direction, stagger, complete, halo, &
                          clock)
  real, dimension(:,:),  intent(inout) :: u_cmpt    !< The nominal zonal (u) component of the vector
                                                    !! pair which is having its halos points
                                                    !! exchanged.
  real, dimension(:,:),  intent(inout) :: v_cmpt    !< The nominal meridional (v) component of the
                                                    !! vector pair which is having its halos points
                                                    !! exchanged.
  type(MOM_domain_type), intent(inout) :: MOM_dom   !< The MOM_domain_type containing the mpp_domain
                                                    !! needed to determine where data should be
                                                    !! sent.
  integer,     optional, intent(in)    :: direction !< An optional integer indicating which
      !! directions the data should be sent.  It is TO_ALL or the sum of any of TO_EAST, TO_WEST,
      !! TO_NORTH, and TO_SOUTH, possibly plus SCALAR_PAIR if these are paired non-directional
      !! scalars discretized at the typical vector component locations.  For example, TO_EAST sends
      !! the data to the processor to the east, so the halos on the western side are filled. TO_ALL
      !! is the default if omitted.
  integer,     optional, intent(in)    :: stagger   !< An optional flag, which may be one of A_GRID,
                     !! BGRID_NE, or CGRID_NE, indicating where the two components of the vector are
                     !! discretized. Omitting stagger is the same as setting it to CGRID_NE.
  logical,     optional, intent(in)    :: complete  !< An optional argument indicating whether the
                                     !! halo updates should be completed before progress resumes.
                                     !! Omitting complete is the same as setting complete to .true.
  integer,     optional, intent(in)    :: halo      !< The size of the halo to update - the full
                                                    !! halo by default.
  integer,     optional, intent(in)    :: clock     !< The handle for a cpu time clock that should be
                                                    !! started then stopped to time this routine.

  ! Local variables

end subroutine pass_vector_2d
module subroutine fill_vector_symmetric_edges_2d(u_cmpt, v_cmpt, MOM_dom, stagger, scalar, &
                                          clock)
  real, dimension(:,:),  intent(inout) :: u_cmpt  !< The nominal zonal (u) component of the vector
                                                  !! pair which is having its halos points
                                                  !! exchanged.
  real, dimension(:,:),  intent(inout) :: v_cmpt  !< The nominal meridional (v) component of the
                                                  !! vector pair which is having its halos points
                                                  !! exchanged.
  type(MOM_domain_type), intent(inout) :: MOM_dom !< The MOM_domain_type containing the mpp_domain
                                                  !! needed to determine where data should be
                                                  !! sent.
  integer,     optional, intent(in)    :: stagger !< An optional flag, which may be one of A_GRID,
                     !! BGRID_NE, or CGRID_NE, indicating where the two components of the vector are
                     !! discretized. Omitting stagger is the same as setting it to CGRID_NE.
  logical,     optional, intent(in)    :: scalar  !< An optional argument indicating whether.
  integer,     optional, intent(in)    :: clock    !< The handle for a cpu time clock that should be
                                                   !! started then stopped to time this routine.

  ! Local variables

end subroutine fill_vector_symmetric_edges_2d
module subroutine pass_vector_3d(u_cmpt, v_cmpt, MOM_dom, direction, stagger, complete, halo, &
                          clock)
  real, dimension(:,:,:), intent(inout) :: u_cmpt   !< The nominal zonal (u) component of the vector
                                                    !! pair which is having its halos points
                                                    !! exchanged.
  real, dimension(:,:,:), intent(inout) :: v_cmpt   !< The nominal meridional (v) component of the
                                                    !! vector pair which is having its halos points
                                                    !! exchanged.
  type(MOM_domain_type),  intent(inout) :: MOM_dom  !< The MOM_domain_type containing the mpp_domain
                                                    !! needed to determine where data should be
                                                    !! sent.
  integer,      optional, intent(in)    :: direction !< An optional integer indicating which
      !! directions the data should be sent.  It is TO_ALL or the sum of any of TO_EAST, TO_WEST,
      !! TO_NORTH, and TO_SOUTH, possibly plus SCALAR_PAIR if these are paired non-directional
      !! scalars discretized at the typical vector component locations.  For example, TO_EAST sends
      !! the data to the processor to the east, so the halos on the western side are filled. TO_ALL
      !! is the default if omitted.
  integer,      optional, intent(in)    :: stagger  !< An optional flag, which may be one of A_GRID,
                     !! BGRID_NE, or CGRID_NE, indicating where the two components of the vector are
                     !! discretized. Omitting stagger is the same as setting it to CGRID_NE.
  logical,      optional, intent(in)    :: complete !< An optional argument indicating whether the
                                     !! halo updates should be completed before progress resumes.
                                     !! Omitting complete is the same as setting complete to .true.
  integer,      optional, intent(in)    :: halo     !< The size of the halo to update - the full
                                                    !! halo by default.
  integer,      optional, intent(in)    :: clock    !< The handle for a cpu time clock that should be
                                                    !! started then stopped to time this routine.

  ! Local variables

end subroutine pass_vector_3d
module function pass_vector_start_2d(u_cmpt, v_cmpt, MOM_dom, direction, stagger, complete, halo, &
                              clock)
  real, dimension(:,:),   intent(inout) :: u_cmpt   !< The nominal zonal (u) component of the vector
                                                    !! pair which is having its halos points
                                                    !! exchanged.
  real, dimension(:,:),   intent(inout) :: v_cmpt   !< The nominal meridional (v) component of the
                                                    !! vector pair which is having its halos points
                                                    !! exchanged.
  type(MOM_domain_type),  intent(inout) :: MOM_dom  !< The MOM_domain_type containing the mpp_domain
                                                    !! needed to determine where data should be
                                                    !! sent.
  integer,      optional, intent(in)    :: direction !< An optional integer indicating which
      !! directions the data should be sent.  It is TO_ALL or the sum of any of TO_EAST, TO_WEST,
      !! TO_NORTH, and TO_SOUTH, possibly plus SCALAR_PAIR if these are paired non-directional
      !! scalars discretized at the typical vector component locations.  For example, TO_EAST sends
      !! the data to the processor to the east, so the halos on the western side are filled. TO_ALL
      !! is the default if omitted.
  integer,      optional, intent(in)    :: stagger  !< An optional flag, which may be one of A_GRID,
                     !! BGRID_NE, or CGRID_NE, indicating where the two components of the vector are
                     !! discretized. Omitting stagger is the same as setting it to CGRID_NE.
  logical,      optional, intent(in)    :: complete !< An optional argument indicating whether the
                                     !! halo updates should be completed before progress resumes.
                                     !! Omitting complete is the same as setting complete to .true.
  integer,      optional, intent(in)    :: halo     !< The size of the halo to update - the full
                                                    !! halo by default.
  integer,      optional, intent(in)    :: clock    !< The handle for a cpu time clock that should be
                                                    !! started then stopped to time this routine.
  integer                               :: pass_vector_start_2d !< The integer index for this
                                                                !! update.

  ! Local variables

end function pass_vector_start_2d
module function pass_vector_start_3d(u_cmpt, v_cmpt, MOM_dom, direction, stagger, complete, halo, &
                              clock)
  real, dimension(:,:,:), intent(inout) :: u_cmpt   !< The nominal zonal (u) component of the vector
                                                    !! pair which is having its halos points
                                                    !! exchanged.
  real, dimension(:,:,:), intent(inout) :: v_cmpt   !< The nominal meridional (v) component of the
                                                    !! vector pair which is having its halos points
                                                    !! exchanged.
  type(MOM_domain_type),  intent(inout) :: MOM_dom  !< The MOM_domain_type containing the mpp_domain
                                                    !! needed to determine where data should be
                                                    !! sent.
  integer,      optional, intent(in)    :: direction !< An optional integer indicating which
      !! directions the data should be sent.  It is TO_ALL or the sum of any of TO_EAST, TO_WEST,
      !! TO_NORTH, and TO_SOUTH, possibly plus SCALAR_PAIR if these are paired non-directional
      !! scalars discretized at the typical vector component locations.  For example, TO_EAST sends
      !! the data to the processor to the east, so the halos on the western side are filled. TO_ALL
      !! is the default if omitted.
  integer,      optional, intent(in)    :: stagger  !< An optional flag, which may be one of A_GRID,
                     !! BGRID_NE, or CGRID_NE, indicating where the two components of the vector are
                     !! discretized. Omitting stagger is the same as setting it to CGRID_NE.
  logical,      optional, intent(in)    :: complete !< An optional argument indicating whether the
                                     !! halo updates should be completed before progress resumes.
                                     !! Omitting complete is the same as setting complete to .true.
  integer,      optional, intent(in)    :: halo     !< The size of the halo to update - the full
                                                    !! halo by default.
  integer,      optional, intent(in)    :: clock    !< The handle for a cpu time clock that should be
                                                    !! started then stopped to time this routine.
  integer                               :: pass_vector_start_3d !< The integer index for this
                                                                !! update.
  ! Local variables

end function pass_vector_start_3d
module subroutine pass_vector_complete_2d(id_update, u_cmpt, v_cmpt, MOM_dom, direction, stagger, halo, &
                                   clock)
  integer,                intent(in)    :: id_update !< The integer id of this update which has been
                                                    !! returned from a previous call to
                                                    !! pass_var_start.
  real, dimension(:,:),   intent(inout) :: u_cmpt   !< The nominal zonal (u) component of the vector
                                                    !! pair which is having its halos points
                                                    !! exchanged.
  real, dimension(:,:),   intent(inout) :: v_cmpt   !< The nominal meridional (v) component of the
                                                    !! vector pair which is having its halos points
                                                    !! exchanged.
  type(MOM_domain_type),  intent(inout) :: MOM_dom  !< The MOM_domain_type containing the mpp_domain
                                                    !! needed to determine where data should be
                                                    !! sent.
  integer,      optional, intent(in)    :: direction !< An optional integer indicating which
      !! directions the data should be sent.  It is TO_ALL or the sum of any of TO_EAST, TO_WEST,
      !! TO_NORTH, and TO_SOUTH, possibly plus SCALAR_PAIR if these are paired non-directional
      !! scalars discretized at the typical vector component locations.  For example, TO_EAST sends
      !! the data to the processor to the east, so the halos on the western side are filled. TO_ALL
      !! is the default if omitted.
  integer,      optional, intent(in)    :: stagger  !< An optional flag, which may be one of A_GRID,
                     !! BGRID_NE, or CGRID_NE, indicating where the two components of the vector are
                     !! discretized. Omitting stagger is the same as setting it to CGRID_NE.
  integer,      optional, intent(in)    :: halo     !< The size of the halo to update - the full
                                                    !! halo by default.
  integer,      optional, intent(in)    :: clock    !< The handle for a cpu time clock that should be
                                                    !! started then stopped to time this routine.
  ! Local variables

end subroutine pass_vector_complete_2d
module subroutine pass_vector_complete_3d(id_update, u_cmpt, v_cmpt, MOM_dom, direction, stagger, halo, &
                                   clock)
  integer,                intent(in)    :: id_update !< The integer id of this update which has been
                                                    !! returned from a previous call to
                                                    !! pass_var_start.
  real, dimension(:,:,:), intent(inout) :: u_cmpt   !< The nominal zonal (u) component of the vector
                                                    !! pair which is having its halos points
                                                    !! exchanged.
  real, dimension(:,:,:), intent(inout) :: v_cmpt   !< The nominal meridional (v) component of the
                                                    !! vector pair which is having its halos points
                                                    !! exchanged.
  type(MOM_domain_type),  intent(inout) :: MOM_dom  !< The MOM_domain_type containing the mpp_domain
                                                    !! needed to determine where data should be
                                                    !! sent.
  integer,      optional, intent(in)    :: direction !< An optional integer indicating which
      !! directions the data should be sent.  It is TO_ALL or the sum of any of TO_EAST, TO_WEST,
      !! TO_NORTH, and TO_SOUTH, possibly plus SCALAR_PAIR if these are paired non-directional
      !! scalars discretized at the typical vector component locations.  For example, TO_EAST sends
      !! the data to the processor to the east, so the halos on the western side are filled. TO_ALL
      !! is the default if omitted.
  integer,      optional, intent(in)    :: stagger  !< An optional flag, which may be one of A_GRID,
                     !! BGRID_NE, or CGRID_NE, indicating where the two components of the vector are
                     !! discretized. Omitting stagger is the same as setting it to CGRID_NE.
  integer,      optional, intent(in)    :: halo     !< The size of the halo to update - the full
                                                    !! halo by default.
  integer,      optional, intent(in)    :: clock    !< The handle for a cpu time clock that should be
                                                    !! started then stopped to time this routine.
  ! Local variables

end subroutine pass_vector_complete_3d
module subroutine create_var_group_pass_2d(group, array, MOM_dom, sideflag, position, &
                                    halo, clock)
  type(group_pass_type),  intent(inout) :: group    !< The data type that store information for
                                                    !! group update. This data will be used in
                                                    !! do_group_pass.
  real, dimension(:,:),   intent(inout) :: array    !< The array which is having its halos points
                                                    !! exchanged.
  type(MOM_domain_type),  intent(inout) :: MOM_dom  !< The MOM_domain_type containing the mpp_domain
                                                    !! needed to determine where data should be
                                                    !! sent.
  integer,      optional, intent(in)    :: sideflag !< An optional integer indicating which
      !! directions the data should be sent. It is TO_ALL or the sum of any of TO_EAST, TO_WEST,
      !! TO_NORTH, and TO_SOUTH.  For example, TO_EAST sends the data to the processor to the east,
      !! so the halos on the western side are filled.  TO_ALL is the default if sideflag is omitted.
  integer,      optional, intent(in)    :: position !< An optional argument indicating the position.
                                                    !! This is CENTER by default and is often CORNER,
                                                    !! but could also be EAST_FACE or NORTH_FACE.
  integer,      optional, intent(in)    :: halo     !< The size of the halo to update - the full
                                                    !! halo by default.
  integer,      optional, intent(in)    :: clock    !< The handle for a cpu time clock that should be
                                                    !! started then stopped to time this routine.
  ! Local variables

end subroutine create_var_group_pass_2d
module subroutine create_var_group_pass_3d(group, array, MOM_dom, sideflag, position, halo, &
                                    clock)
  type(group_pass_type),  intent(inout) :: group    !< The data type that store information for
                                                    !! group update. This data will be used in
                                                    !! do_group_pass.
  real, dimension(:,:,:), intent(inout) :: array    !< The array which is having its halos points
                                                    !! exchanged.
  type(MOM_domain_type),  intent(inout) :: MOM_dom  !< The MOM_domain_type containing the mpp_domain
                                                    !! needed to determine where data should be
                                                    !! sent.
  integer,      optional, intent(in)    :: sideflag !< An optional integer indicating which
      !! directions the data should be sent. It is TO_ALL or the sum of any of TO_EAST, TO_WEST,
      !! TO_NORTH, and TO_SOUTH.  For example, TO_EAST sends the data to the processor to the east,
      !! so the halos on the western side are filled.  TO_ALL is the default if sideflag is omitted.
  integer,      optional, intent(in)    :: position !< An optional argument indicating the position.
                                                    !! This is CENTER by default and is often CORNER,
                                                    !! but could also be EAST_FACE or NORTH_FACE.
  integer,      optional, intent(in)    :: halo     !< The size of the halo to update - the full
                                                    !! halo by default.
  integer,      optional, intent(in)    :: clock    !< The handle for a cpu time clock that should be
                                                    !! started then stopped to time this routine.
  ! Local variables

end subroutine create_var_group_pass_3d
module subroutine create_vector_group_pass_2d(group, u_cmpt, v_cmpt, MOM_dom, direction, stagger, halo, &
                                       clock)
  type(group_pass_type),  intent(inout) :: group    !< The data type that store information for
                                                    !! group update. This data will be used in
                                                    !! do_group_pass.
  real, dimension(:,:),   intent(inout) :: u_cmpt   !< The nominal zonal (u) component of the vector
                                                    !! pair which is having its halos points
                                                    !! exchanged.
  real, dimension(:,:),   intent(inout) :: v_cmpt   !< The nominal meridional (v) component of the
                                                    !! vector pair which is having its halos points
                                                    !! exchanged.

  type(MOM_domain_type),  intent(inout) :: MOM_dom  !< The MOM_domain_type containing the mpp_domain
                                                    !! needed to determine where data should be
                                                    !! sent
  integer,      optional, intent(in)    :: direction !< An optional integer indicating which
      !! directions the data should be sent.  It is TO_ALL or the sum of any of TO_EAST, TO_WEST,
      !! TO_NORTH, and TO_SOUTH, possibly plus SCALAR_PAIR if these are paired non-directional
      !! scalars discretized at the typical vector component locations.  For example, TO_EAST sends
      !! the data to the processor to the east, so the halos on the western side are filled. TO_ALL
      !! is the default if omitted.
  integer,      optional, intent(in)    :: stagger  !< An optional flag, which may be one of A_GRID,
                     !! BGRID_NE, or CGRID_NE, indicating where the two components of the vector are
                     !! discretized. Omitting stagger is the same as setting it to CGRID_NE.
  integer,      optional, intent(in)    :: halo     !< The size of the halo to update - the full
                                                    !! halo by default.
  integer,      optional, intent(in)    :: clock    !< The handle for a cpu time clock that should be
                                                    !! started then stopped to time this routine.
  ! Local variables

end subroutine create_vector_group_pass_2d
module subroutine create_vector_group_pass_3d(group, u_cmpt, v_cmpt, MOM_dom, direction, stagger, halo, &
                                       clock)
  type(group_pass_type),  intent(inout) :: group    !< The data type that store information for
                                                    !! group update. This data will be used in
                                                    !! do_group_pass.
  real, dimension(:,:,:), intent(inout) :: u_cmpt   !< The nominal zonal (u) component of the vector
                                                    !! pair which is having its halos points
                                                    !! exchanged.
  real, dimension(:,:,:), intent(inout) :: v_cmpt   !< The nominal meridional (v) component of the
                                                    !! vector pair which is having its halos points
                                                    !! exchanged.

  type(MOM_domain_type),  intent(inout) :: MOM_dom  !< The MOM_domain_type containing the mpp_domain
                                                    !! needed to determine where data should be
                                                    !! sent.
  integer,      optional, intent(in)    :: direction !< An optional integer indicating which
      !! directions the data should be sent.  It is TO_ALL or the sum of any of TO_EAST, TO_WEST,
      !! TO_NORTH, and TO_SOUTH, possibly plus SCALAR_PAIR if these are paired non-directional
      !! scalars discretized at the typical vector component locations.  For example, TO_EAST sends
      !! the data to the processor to the east, so the halos on the western side are filled. TO_ALL
      !! is the default if omitted.
  integer,      optional, intent(in)    :: stagger  !< An optional flag, which may be one of A_GRID,
                     !! BGRID_NE, or CGRID_NE, indicating where the two components of the vector are
                     !! discretized. Omitting stagger is the same as setting it to CGRID_NE.
  integer,      optional, intent(in)    :: halo     !< The size of the halo to update - the full
                                                    !! halo by default.
  integer,      optional, intent(in)    :: clock    !< The handle for a cpu time clock that should be
                                                    !! started then stopped to time this routine.

  ! Local variables

end subroutine create_vector_group_pass_3d
module subroutine do_group_pass(group, MOM_dom, clock)
  type(group_pass_type), intent(inout) :: group     !< The data type that store information for
                                                    !! group update. This data will be used in
                                                    !! do_group_pass.
  type(MOM_domain_type), intent(inout) :: MOM_dom   !< The MOM_domain_type containing the mpp_domain
                                                    !! needed to determine where data should be
                                                    !! sent.
  integer,     optional, intent(in)    :: clock     !< The handle for a cpu time clock that should be
                                                    !! started then stopped to time this routine.

end subroutine do_group_pass
module subroutine start_group_pass(group, MOM_dom, clock)
  type(group_pass_type), intent(inout) :: group    !< The data type that store information for
                                                   !! group update. This data will be used in
                                                   !! do_group_pass.
  type(MOM_domain_type), intent(inout) :: MOM_dom  !< The MOM_domain_type containing the mpp_domain
                                                   !! needed to determine where data should be
                                                   !! sent.
  integer,     optional, intent(in)    :: clock    !< The handle for a cpu time clock that should be
                                                   !! started then stopped to time this routine.


end subroutine start_group_pass
module subroutine complete_group_pass(group, MOM_dom, clock)
  type(group_pass_type), intent(inout) :: group    !< The data type that store information for
                                                   !! group update. This data will be used in
                                                   !! do_group_pass.
  type(MOM_domain_type), intent(inout) :: MOM_dom  !< The MOM_domain_type containing the mpp_domain
                                                   !! needed to determine where data should be
                                                   !! sent.
  integer,     optional, intent(in)    :: clock    !< The handle for a cpu time clock that should be
                                                   !! started then stopped to time this routine.

end subroutine complete_group_pass
module subroutine redistribute_array_2d(Domain1, array1, Domain2, array2, complete)
  type(domain2d), &
           intent(in)  :: Domain1 !< The MOM domain from which to extract information.
  real, dimension(:,:), intent(in) :: array1 !< The array from which to extract information.
  type(domain2d), &
           intent(in)  :: Domain2 !< The MOM domain receiving information.
  real, dimension(:,:), intent(out) :: array2 !< The array receiving information.
  logical, optional, intent(in) :: complete  !< If true, finish communication before proceeding.

  ! Local variables

end subroutine redistribute_array_2d
module subroutine redistribute_array_3d(Domain1, array1, Domain2, array2, complete)
  type(domain2d), &
           intent(in)  :: Domain1 !< The MOM domain from which to extract information.
  real, dimension(:,:,:), intent(in) :: array1 !< The array from which to extract information.
  type(domain2d), &
           intent(in)  :: Domain2 !< The MOM domain receiving information.
  real, dimension(:,:,:), intent(out) :: array2 !< The array receiving information.
  logical, optional, intent(in) :: complete  !< If true, finish communication before proceeding.

  ! Local variables

end subroutine redistribute_array_3d
module subroutine redistribute_array_4d(Domain1, array1, Domain2, array2, complete)
  type(domain2d), &
           intent(in)  :: Domain1 !< The MOM domain from which to extract information.
  real, dimension(:,:,:,:), intent(in) :: array1 !< The array from which to extract information.
  type(domain2d), &
           intent(in)  :: Domain2 !< The MOM domain receiving information.
  real, dimension(:,:,:,:), intent(out) :: array2 !< The array receiving information.
  logical, optional, intent(in) :: complete  !< If true, finish communication before proceeding.

  ! Local variables

end subroutine redistribute_array_4d
module subroutine rescale_comp_data_4d(domain, array, scale, zero_zeros)
  type(MOM_domain_type),    intent(in)    :: domain !< MOM domain from which to extract information
  real, dimension(:,:,:,:), intent(inout) :: array  !< The array which is having the data in its
                                                    !! computational domain rescaled
  real,                     intent(in)    :: scale  !< A scaling factor by which to multiply the
                                                    !! values in the computational domain of array
  logical,        optional, intent(in)    :: zero_zeros !< If present and true, convert negative zeros
                                                    !! into ordinary signless zeros.

end subroutine rescale_comp_data_4d
module subroutine rescale_comp_data_3d(domain, array, scale, zero_zeros)
  type(MOM_domain_type),  intent(in)    :: domain !< MOM domain from which to extract information
  real, dimension(:,:,:), intent(inout) :: array  !< The array which is having the data in its
                                                  !! computational domain rescaled
  real,                   intent(in)    :: scale  !< A scaling factor by which to multiply the
                                                  !! values in the computational domain of array
  logical,      optional, intent(in)    :: zero_zeros !< If present and true, convert negative zeros
                                                  !! into ordinary signless zeros.

end subroutine rescale_comp_data_3d
module subroutine rescale_comp_data_2d(domain, array, scale, zero_zeros)
  type(MOM_domain_type), intent(in)    :: domain !< MOM domain from which to extract information
  real, dimension(:,:),  intent(inout) :: array  !< The array which is having the data in its
                                                 !! computational domain rescaled
  real,                  intent(in)    :: scale  !< A scaling factor by which to multiply the
                                                 !! values in the computational domain of array
  logical,      optional, intent(in)   :: zero_zeros !< If present and true, convert negative zeros
                                                  !! into ordinary signless zeros.

end subroutine rescale_comp_data_2d
module subroutine create_MOM_domain(MOM_dom, n_global, n_halo, reentrant, tripolar_N, layout, io_layout, &
                             domain_name, mask_table, symmetric, thin_halos, nonblocking)
  type(MOM_domain_type),      pointer    :: MOM_dom   !< A pointer to the MOM_domain_type being defined here.
  integer, dimension(2),      intent(in) :: n_global  !< The number of points on the global grid in
                                                      !! the i- and j-directions
  integer, dimension(2),      intent(in) :: n_halo    !< The number of halo points on each processor
  logical, dimension(2),      intent(in) :: reentrant !< If true the grid is periodic in the i- and j- directions
  logical,                    intent(in) :: tripolar_N !< If true the grid uses northern tripolar connectivity
  integer, dimension(2),      intent(in) :: layout    !< The layout of logical PEs in the i- and j-directions.
  integer, dimension(2), optional, intent(in) :: io_layout !< The layout for parallel input and output.
  character(len=*), optional, intent(in) :: domain_name !< A name for this domain, "MOM" if missing.
  character(len=*), optional, intent(in) :: mask_table !< The full relative or absolute path to the mask table.
  logical,          optional, intent(in) :: symmetric !< If present, this specifies whether this domain
                                                      !! uses symmetric memory, or true if missing.
  logical,          optional, intent(in) :: thin_halos !< If present, this specifies whether to permit the use of
                                                      !! thin halo updates, or true if missing.
  logical,          optional, intent(in) :: nonblocking !< If present, this specifies whether to permit the use of
                                                      !! nonblocking halo updates, or false if missing.

  ! local variables

end subroutine create_MOM_domain
module subroutine deallocate_MOM_domain(MOM_domain, cursory)
  type(MOM_domain_type), pointer :: MOM_domain !< A pointer to the MOM_domain_type being deallocated
  logical,  optional, intent(in) :: cursory    !< If true do not deallocate fields associated
                                               !! with the underlying infrastructure

end subroutine deallocate_MOM_domain
module function MOM_thread_affinity_set()
  ! Local variables
  !$ integer :: ocean_nthreads       ! Number of openMP threads
  !$ integer :: omp_get_num_threads  ! An openMP function that returns the number of threads
  logical :: MOM_thread_affinity_set

end function MOM_thread_affinity_set
module subroutine set_MOM_thread_affinity(ocean_nthreads, ocean_hyper_thread)
  integer, intent(in) :: ocean_nthreads     !< Number of openMP threads to use for the ocean model
  logical, intent(in) :: ocean_hyper_thread !< If true, use hyper threading

  ! Local variables
  !$ integer :: omp_get_thread_num, omp_get_num_threads !< These are the results of openMP functions

  !$ call fms_affinity_init()  ! fms_affinity_init can be safely called more than once.
  !$ call fms_affinity_set('OCEAN', ocean_hyper_thread, ocean_nthreads)
  !$ call omp_set_num_threads(ocean_nthreads)
  !$OMP PARALLEL
  !$ write(6,*) "MOM_domains_mod OMPthreading ", fms_affinity_get(), omp_get_thread_num(), omp_get_num_threads()
  !$ flush(6)
  !$OMP END PARALLEL
end subroutine set_MOM_thread_affinity
module subroutine get_domain_components_MD(MOM_dom, x_domain, y_domain)
  type(MOM_domain_type),    intent(in)    :: MOM_dom  !< The MOM_domain whose contents are being extracted
  type(domain1D), optional, intent(inout) :: x_domain !< The 1-d logical x-domain
  type(domain1D), optional, intent(inout) :: y_domain !< The 1-d logical y-domain

end subroutine get_domain_components_MD
module subroutine get_domain_components_d2D(domain, x_domain, y_domain)
  type(domain2D),           intent(in)    :: domain  !< The 2D domain whose contents are being extracted
  type(domain1D), optional, intent(inout) :: x_domain !< The 1-d logical x-domain
  type(domain1D), optional, intent(inout) :: y_domain !< The 1-d logical y-domain

end subroutine get_domain_components_d2D
module subroutine clone_MD_to_MD(MD_in, MOM_dom, min_halo, halo_size, symmetric, domain_name, &
                          turns, refine, extra_halo, io_layout)
  type(MOM_domain_type), target, intent(in) :: MD_in  !< An existing MOM_domain
  type(MOM_domain_type), pointer :: MOM_dom
                                  !< A pointer to a MOM_domain that will be
                                  !! allocated if it is unassociated, and will have data
                                  !! copied from MD_in
  integer, dimension(2), &
               optional, intent(inout) :: min_halo !< If present, this sets the
                                  !! minimum halo size for this domain in the i- and j-
                                  !! directions, and returns the actual halo size used.
  integer,     optional, intent(in)    :: halo_size !< If present, this sets the halo
                                  !! size for the domain in the i- and j-directions.
                                  !! min_halo and halo_size can not both be present.
  logical,     optional, intent(in)    :: symmetric !< If present, this specifies
                                  !! whether the new domain is symmetric, regardless of
                                  !! whether the macro SYMMETRIC_MEMORY_ is defined.
  character(len=*), &
               optional, intent(in)    :: domain_name !< A name for the new domain, copied
                                  !! from MD_in if missing.
  integer, optional, intent(in) :: turns   !< Number of quarter turns
  integer, optional, intent(in) :: refine  !< A factor by which to enhance the grid resolution.
  integer, optional, intent(in) :: extra_halo !< An extra number of points in the halos
                                  !! compared with MD_in
  integer, optional, intent(in) :: io_layout(2)
    !< A user-defined IO layout to replace the domain's IO layout


                                             ! The sum of exni must equal MOM_dom%niglobal.
                                             ! The sum of exni must equal MOM_dom%niglobal.

end subroutine clone_MD_to_MD
module subroutine clone_MD_to_d2D(MD_in, mpp_domain, min_halo, halo_size, symmetric, &
                           domain_name, turns, xextent, yextent, coarsen)
  type(MOM_domain_type), intent(in)    :: MD_in !< An existing MOM_domain to be cloned
  type(domain2d),        intent(inout) :: mpp_domain !< The new mpp_domain to be set up
  integer, dimension(2), &
               optional, intent(inout) :: min_halo !< If present, this sets the
                                  !! minimum halo size for this domain in the i- and j-
                                  !! directions, and returns the actual halo size used.
  integer,     optional, intent(in)    :: halo_size !< If present, this sets the halo
                                  !! size for the domain in the i- and j-directions.
                                  !! min_halo and halo_size can not both be present.
  logical,     optional, intent(in)    :: symmetric !< If present, this specifies
                                  !! whether the new domain is symmetric, regardless of
                                  !! whether the macro SYMMETRIC_MEMORY_ is defined or
                                  !! whether MD_in is symmetric.
  character(len=*), &
               optional, intent(in)    :: domain_name !< A name for the new domain, "MOM"
                                  !! if missing.
  integer, optional, intent(in) :: turns  !< Number of quarter turns - not implemented here.
  integer, optional, intent(in) :: coarsen !< A factor by which to coarsen this grid.
                                  !! The default of 1 is for no coarsening.
  integer, dimension(:), optional, intent(in) :: xextent !< The number of grid points in the
                                  !! tracer computational domain for division of the x-layout.
  integer, dimension(:), optional, intent(in) :: yextent !< The number of grid points in the
                                  !! tracer computational domain for division of the y-layout.


end subroutine clone_MD_to_d2D
module subroutine get_domain_extent_MD(Domain, isc, iec, jsc, jec, isd, ied, jsd, jed, &
                                isg, ieg, jsg, jeg, idg_offset, jdg_offset, &
                                symmetric, local_indexing, index_offset, coarsen)
  type(MOM_domain_type), &
                     intent(in)  :: Domain !< The MOM domain from which to extract information
  integer,           intent(out) :: isc    !< The start i-index of the computational domain
  integer,           intent(out) :: iec    !< The end i-index of the computational domain
  integer,           intent(out) :: jsc    !< The start j-index of the computational domain
  integer,           intent(out) :: jec    !< The end j-index of the computational domain
  integer,           intent(out) :: isd    !< The start i-index of the data domain
  integer,           intent(out) :: ied    !< The end i-index of the data domain
  integer,           intent(out) :: jsd    !< The start j-index of the data domain
  integer,           intent(out) :: jed    !< The end j-index of the data domain
  integer, optional, intent(out) :: isg    !< The start i-index of the global domain
  integer, optional, intent(out) :: ieg    !< The end i-index of the global domain
  integer, optional, intent(out) :: jsg    !< The start j-index of the global domain
  integer, optional, intent(out) :: jeg    !< The end j-index of the global domain
  integer, optional, intent(out) :: idg_offset !< The offset between the corresponding global and
                                           !! data i-index spaces.
  integer, optional, intent(out) :: jdg_offset !< The offset between the corresponding global and
                                           !! data j-index spaces.
  logical, optional, intent(out) :: symmetric  !< True if symmetric memory is used.
  logical, optional, intent(in)  :: local_indexing !< If true, local tracer array indices start at 1,
                                           !! as in most MOM6 code.  The default is true.
  integer, optional, intent(in)  :: index_offset   !< A fixed additional offset to all indices. This
                                           !! can be useful for some types of debugging with
                                           !! dynamic memory allocation.  The default is 0.
  integer, optional, intent(in)  :: coarsen !< The index of the factor by which the grid is coarsened.
                                           !!  The default is 0, for no coarsening.

  ! Local variables

end subroutine get_domain_extent_MD
module subroutine get_domain_extent_d2D(Domain, isc, iec, jsc, jec, isd, ied, jsd, jed)
  type(domain2d),    intent(in)  :: Domain !< The MOM domain from which to extract information
  integer,           intent(out) :: isc    !< The start i-index of the computational domain
  integer,           intent(out) :: iec    !< The end i-index of the computational domain
  integer,           intent(out) :: jsc    !< The start j-index of the computational domain
  integer,           intent(out) :: jec    !< The end j-index of the computational domain
  integer, optional, intent(out) :: isd    !< The start i-index of the data domain
  integer, optional, intent(out) :: ied    !< The end i-index of the data domain
  integer, optional, intent(out) :: jsd    !< The start j-index of the data domain
  integer, optional, intent(out) :: jed    !< The end j-index of the data domain

  ! Local variables

end subroutine get_domain_extent_d2D
module subroutine get_simple_array_i_ind(domain, size, is, ie, symmetric)
  type(MOM_domain_type), intent(in)  :: domain !< MOM domain from which to extract information
  integer,               intent(in)  :: size   !< The i-array size
  integer,               intent(out) :: is     !< The computational domain starting i-index.
  integer,               intent(out) :: ie     !< The computational domain ending i-index.
  logical,     optional, intent(in)  :: symmetric !< If present, indicates whether symmetric sizes
                                               !! can be considered.
  ! Local variables

end subroutine get_simple_array_i_ind
module subroutine get_simple_array_j_ind(domain, size, js, je, symmetric)
  type(MOM_domain_type), intent(in)  :: domain !< MOM domain from which to extract information
  integer,               intent(in)  :: size   !< The j-array size
  integer,               intent(out) :: js     !< The computational domain starting j-index.
  integer,               intent(out) :: je     !< The computational domain ending j-index.
  logical,     optional, intent(in)  :: symmetric !< If present, indicates whether symmetric sizes
                                               !! can be considered.
  ! Local variables

end subroutine get_simple_array_j_ind
module subroutine invert(array)
  integer, dimension(:), intent(inout) :: array !< The 1-d array to invert
end subroutine invert
module subroutine get_global_shape(domain, niglobal, njglobal)
  type(MOM_domain_type), intent(in)  :: domain   !< MOM domain from which to extract information
  integer,               intent(out) :: niglobal !< i-index global size of h-point arrays
  integer,               intent(out) :: njglobal !< j-index global size of h-point arrays

end subroutine get_global_shape
module subroutine compute_block_extent(isg, ieg, ndivs, ibegin, iend)
  integer,               intent(in)  :: isg    !< The starting index of the global index space
  integer,               intent(in)  :: ieg    !< The ending index of the global index space
  integer,               intent(in)  :: ndivs  !< The number of divisions
  integer, dimension(:), intent(out) :: ibegin !< The starting index of each division
  integer, dimension(:), intent(out) :: iend   !< The ending index of each division

end subroutine compute_block_extent
module subroutine compute_extent(isg, ieg, ndivs, ibegin, iend)
  integer,               intent(in)  :: isg    !< The starting index of the global index space
  integer,               intent(in)  :: ieg    !< The ending index of the global index space
  integer,               intent(in)  :: ndivs  !< The number of divisions
  integer, dimension(:), intent(out) :: ibegin !< The starting index of each division
  integer, dimension(:), intent(out) :: iend   !< The ending index of each division

end subroutine compute_extent
module subroutine broadcast_domain(domain)
  type(domain2d),  intent(inout) :: domain !< The domain2d type that will be shared across PEs.

end subroutine broadcast_domain
module subroutine global_field(domain, local, global)
  type(domain2d),       intent(inout) :: domain !< The domain2d type that describes the decomposition
  real, dimension(:,:), intent(in)    :: local  !< The portion of the array on the local PE
  real, dimension(:,:), intent(out)   :: global !< The whole global array

end subroutine global_field
logical module function same_domain(domain_a, domain_b)
  type(domain2D), intent(in) :: domain_a !< The first domain in the comparison
  type(domain2D), intent(in) :: domain_b !< The second domain in the comparison

  ! Local variables

  ! This routine currently does a few checks for consistent domains; more could be added.
end function same_domain
module subroutine get_layout_extents(Domain, extent_i, extent_j)
  type(MOM_domain_type), intent(in)  :: domain !< MOM domain from which to extract information
  integer, dimension(:), allocatable, intent(inout) :: extent_i  !< The number of points in the
                                               !! i-direction in each i-row of the layout
  integer, dimension(:), allocatable, intent(inout) :: extent_j  !< The number of points in the
                                               !! j-direction in each j-row of the layout

end subroutine get_layout_extents
module subroutine set_domain(Domain)
  type(MOM_domain_type), intent(in) :: Domain
    !< MOM domain to be designated as the internal FMS I/O domain

  ! FMS2 does not have domain-based internal FMS I/O operations, so this
  ! function does nothing.
end subroutine set_domain
module subroutine nullify_domain
  ! No internal FMS I/O domain can be assigned, so this function does nothing.
end subroutine nullify_domain
  end interface

end module MOM_domain_infra
