! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> Describes the decomposed MOM domain and has routines for communications across PEs
module MOM_domains

use MOM_coms_infra,       only : MOM_infra_init, MOM_infra_end
use MOM_coms_infra,       only : PE_here, root_PE, num_PEs, broadcast
use MOM_coms_infra,       only : sum_across_PEs, min_across_PEs, max_across_PEs
use MOM_domain_infra,     only : MOM_domain_type, domain2D, domain1D, group_pass_type
use MOM_domain_infra,     only : create_MOM_domain, clone_MOM_domain, deallocate_MOM_domain
use MOM_domain_infra,     only : get_domain_extent, get_domain_components, same_domain
use MOM_domain_infra,     only : compute_block_extent, get_global_shape
use MOM_domain_infra,     only : pass_var, pass_vector, fill_symmetric_edges
use MOM_domain_infra,     only : pass_var_start, pass_var_complete
use MOM_domain_infra,     only : pass_vector_start, pass_vector_complete
use MOM_domain_infra,     only : create_group_pass, do_group_pass
use MOM_domain_infra,     only : start_group_pass, complete_group_pass
use MOM_domain_infra,     only : rescale_comp_data, global_field, redistribute_array, broadcast_domain
use MOM_domain_infra,     only : MOM_thread_affinity_set, set_MOM_thread_affinity
use MOM_domain_infra,     only : AGRID, BGRID_NE, CGRID_NE, SCALAR_PAIR
use MOM_domain_infra,     only : CORNER, CENTER, NORTH_FACE, EAST_FACE
use MOM_domain_infra,     only : To_East, To_West, To_North, To_South, To_All, Omit_Corners
use MOM_domain_infra,     only : compute_extent
use MOM_error_handler,    only : MOM_error, MOM_mesg, NOTE, WARNING, FATAL, is_root_pe
use MOM_file_parser,      only : get_param, log_param, log_version, param_file_type
use MOM_io_infra,         only : file_exists, read_field, open_ASCII_file, close_file, WRITEONLY_FILE
use MOM_string_functions, only : slasher
use MOM_cpu_clock,        only : cpu_clock_id, cpu_clock_begin, cpu_clock_end, CLOCK_ROUTINE
use MOM_unit_scaling,     only : unit_scale_type

implicit none ; private

public :: MOM_infra_init, MOM_infra_end
!  Domain types and creation and destruction routines
public :: MOM_domain_type, domain2D, domain1D
public :: MOM_domains_init, create_MOM_domain, clone_MOM_domain, deallocate_MOM_domain
public :: MOM_thread_affinity_set, set_MOM_thread_affinity
public :: MOM_define_layout
!  Domain query routines
public :: get_domain_extent, get_domain_components, get_global_shape, same_domain
public :: PE_here, root_PE, num_PEs
!  Blocks are not actively used in MOM6, so this routine could be deprecated.
public :: compute_block_extent
!  Single call communication routines
public :: pass_var, pass_vector, fill_symmetric_edges, broadcast
!  Non-blocking communication routines
public :: pass_var_start, pass_var_complete, pass_vector_start, pass_vector_complete
!  Multi-variable group communication routines and type
public :: create_group_pass, do_group_pass, group_pass_type, start_group_pass, complete_group_pass
!  Global reduction routines
public :: sum_across_PEs, min_across_PEs, max_across_PEs
public :: global_field, redistribute_array, broadcast_domain
!  Simple index-convention-invariant array manipulation routine
public :: rescale_comp_data
!> These encoding constants are used to indicate the staggering of scalars and vectors
public :: AGRID, BGRID_NE, CGRID_NE, SCALAR_PAIR
!> These encoding constants are used to indicate the discretization position of a variable
public :: CORNER, CENTER, NORTH_FACE, EAST_FACE
!> These encoding constants indicate communication patterns.  In practice they can be added.
public :: To_East, To_West, To_North, To_South, To_All, Omit_Corners


  interface
module subroutine MOM_domains_init(MOM_dom, param_file, symmetric, static_memory, &
                            NIHALO, NJHALO, NIGLOBAL, NJGLOBAL, NIPROC, NJPROC, &
                            min_halo, domain_name, include_name, param_suffix, US, MOM_dom_unmasked)
  type(MOM_domain_type),           pointer       :: MOM_dom      !< A pointer to the MOM_domain_type
                                                                 !! being defined here.
  type(param_file_type),           intent(in)    :: param_file   !< A structure to parse for
                                                                 !! run-time parameters
  logical, optional,               intent(in)    :: symmetric    !< If present, this specifies
                                                  !! whether this domain is symmetric, regardless of
                                                  !! whether the macro SYMMETRIC_MEMORY_ is defined.
  logical, optional,               intent(in)    :: static_memory !< If present and true, this
                                                  !! domain type is set up for static memory and
                                                  !! error checking of various input values is
                                                  !! performed against those in the input file.
  integer, optional,               intent(in)    :: NIHALO       !< Default halo sizes, required
                                                                 !! with static memory.
  integer, optional,               intent(in)    :: NJHALO       !< Default halo sizes, required
                                                                 !! with static memory.
  integer, optional,               intent(in)    :: NIGLOBAL     !< Total domain sizes, required
                                                                 !! with static memory.
  integer, optional,               intent(in)    :: NJGLOBAL     !< Total domain sizes, required
                                                                 !! with static memory.
  integer, optional,               intent(in)    :: NIPROC       !< Processor counts, required with
                                                                 !! static memory.
  integer, optional,               intent(in)    :: NJPROC       !< Processor counts, required with
                                                                 !! static memory.
  integer, dimension(2), optional, intent(inout) :: min_halo     !< If present, this sets the
                                            !! minimum halo size for this domain in the i- and j-
                                            !! directions, and returns the actual halo size used.
  character(len=*),      optional, intent(in)    :: domain_name  !< A name for this domain, "MOM"
                                                                 !! if missing.
  character(len=*),      optional, intent(in)    :: include_name !< A name for model's include file,
                                                                 !! "MOM_memory.h" if missing.
  character(len=*),      optional, intent(in)    :: param_suffix !< A suffix to apply to
                                                                 !! layout-specific parameters.
  type(unit_scale_type), optional, pointer       :: US           !< A dimensional unit scaling type
  type(MOM_domain_type), optional, pointer       :: MOM_dom_unmasked !< Unmasked MOM domain instance.
                                                                 !! Set to null if masking is not enabled.

  ! Local variables
  !$ integer :: ocean_nthreads       ! Number of openMP threads
  !$ logical :: ocean_omp_hyper_thread ! If true use openMP hyper-threads
                            ! width of the halos that are updated with each call.

  ! This include declares and sets the variable "version".

end subroutine MOM_domains_init
module subroutine MOM_define_layout(n_global, ndivs, layout)
  integer, dimension(2), intent(in)  :: n_global !< The total number of gridpoints in 2 directions
  integer,               intent(in)  :: ndivs    !< The total number of (logical) PEs
  integer, dimension(2), intent(out) :: layout   !< The generated layout of PEs

  ! Local variables

  ! At present, this algorithm is a copy of mpp_define_layout, but it could perhaps be improved?

end subroutine MOM_define_layout
module subroutine gen_auto_mask_table(n_global, reentrant, tripolar_N, npes, param_file, inputdir, filename, layout, US)
  integer, dimension(2), intent(in)         :: n_global   !< The total number of gridpoints in 2 directions
  logical, dimension(2), intent(in)         :: reentrant  !< True if the x- and y- directions are periodic.
  logical,               intent(in)         :: tripolar_N !< A flag indicating whether there is n. tripolar connectivity
  integer,               intent(in)         :: npes       !< The desired number of active PEs.
  type(param_file_type), intent(in)         :: param_file !< A structure to parse for run-time parameters
  character(len=128),    intent(in)         :: inputdir   !< INPUTDIR parameter
  character(len=:), allocatable, intent(in) :: filename   !< Mask table file path (to be auto-generated.)
  integer, dimension(2), intent(out)        :: layout     !< The generated layout of PEs (incl. masked blocks)
  type(unit_scale_type), optional, pointer  :: US         !< A dimensional unit scaling type

  ! Local variables

end subroutine gen_auto_mask_table
module subroutine determine_land_blocks(mask, nx, ny, idiv, jdiv, ibuf, jbuf, num_masked_blocks, mask_table)
  integer, dimension(:,:), intent(in)   :: mask     !< cell masks based on depth and MINIMUM_DEPTH
  integer, intent(in)                   :: nx       !< Total number of gridpoints in x-dir (global)
  integer, intent(in)                   :: ny       !< Total number of gridpoints in y-dir (global)
  integer, intent(in)                   :: idiv     !< number of divisions along x-dir
  integer, intent(in)                   :: jdiv     !< number of divisions along y-dir
  integer, intent(in)                   :: ibuf     !< number of buffer cells in x-dir.
                                                    !! (not necessarily the same as NIHALO)
  integer, intent(in)                   :: jbuf     !< number of buffer cells in y-dir.
                                                    !! (not necessarily the same as NJHALO)
  integer, intent(out)                  :: num_masked_blocks !< the final number of masked blocks
  integer, intent(out), optional        :: mask_table(:,:) !< the resulting array of mask_table
  ! integer

end subroutine determine_land_blocks
module subroutine write_auto_mask_file(mask_table, layout, npes, filename)
  integer, intent(in) :: mask_table(:,:)      !> mask table array to be written out.
  integer, dimension(2), intent(in) :: layout !> PE layout
  integer, intent(in) :: npes                 !> Number of divisions (incl. eliminated ones)
  character(len=:), allocatable, intent(in) :: filename !> file name for the mask_table to be written
  ! local

  ! Eliminate only enough blocks to ensure that the number of active blocks precisely matches the target npes.
end subroutine write_auto_mask_file
  end interface

end module MOM_domains
