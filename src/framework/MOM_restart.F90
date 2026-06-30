! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> The MOM6 facility for reading and writing restart files, and querying what has been read.
module MOM_restart

use, intrinsic :: iso_fortran_env, only : int64
use MOM_array_transform, only : rotate_array, rotate_vector, rotate_array_pair
use MOM_checksums, only : chksum => field_checksum
use MOM_domains, only : PE_here, num_PEs, AGRID, BGRID_NE, CGRID_NE
use MOM_error_handler, only : MOM_error, MOM_mesg, FATAL, WARNING, NOTE, is_root_pe, MOM_get_verbosity
use MOM_file_parser, only : get_param, log_param, log_version, param_file_type
use MOM_grid, only : ocean_grid_type
use MOM_io, only : create_MOM_file, file_exists
use MOM_io, only : MOM_infra_file, MOM_field
use MOM_io, only : MOM_read_data, read_data, MOM_write_field, field_exists
use MOM_io, only : vardesc, var_desc, query_vardesc, modify_vardesc, get_filename_appendix
use MOM_io, only : MULTIPLE, READONLY_FILE, SINGLE_FILE
use MOM_io, only : CENTER, CORNER, NORTH_FACE, EAST_FACE
use MOM_io, only : axis_info, get_axis_info
use MOM_string_functions, only : lowercase
use MOM_time_manager,  only : time_type, time_type_to_real, real_to_time
use MOM_time_manager,  only : days_in_month, get_date, set_date
use MOM_verticalGrid,  only : verticalGrid_type

implicit none ; private

public restart_init, restart_end, restore_state, register_restart_field
public copy_restart_var, copy_restart_vector
public save_restart, query_initialized, set_initialized, only_read_from_restarts
public restart_registry_lock, restart_init_end, vardesc
public restart_files_exist, determine_is_new_run, is_new_run
public register_restart_field_as_obsolete, register_restart_pair
public lock_check

! A note on unit descriptions in comments: MOM6 uses units that can be rescaled for dimensional
! consistency testing. These are noted in comments with units like Z, H, L, and T, along with
! their mks counterparts with notation like "a velocity [Z T-1 ~> m s-1]".  If the units
! vary with the Boussinesq approximation, the Boussinesq variant is given first.
! The functions in this module work with variables with arbitrary units, in which case the
! arbitrary rescaled units are indicated with [A ~> a], while the unscaled units are just [a].

!> A type for making arrays of pointers to 4-d arrays
type p4d
  real, dimension(:,:,:,:), pointer :: p => NULL() !< A pointer to a 4d array in arbitrary rescaled units [A ~> a]
end type p4d

!> A type for making arrays of pointers to 3-d arrays
type p3d
  real, dimension(:,:,:), pointer :: p => NULL() !< A pointer to a 3d array in arbitrary rescaled units [A ~> a]
end type p3d

!> A type for making arrays of pointers to 2-d arrays
type p2d
  real, dimension(:,:), pointer :: p => NULL() !< A pointer to a 2d array in arbitrary rescaled units [A ~> a]
end type p2d

!> A type for making arrays of pointers to 1-d arrays
type p1d
  real, dimension(:), pointer :: p => NULL() !< A pointer to a 1d array in arbitrary rescaled units [A ~> a]
end type p1d

!> A type for making arrays of pointers to scalars
type p0d
  real, pointer :: p => NULL() !< A pointer to a scalar in arbitrary rescaled units [A ~> a]
end type p0d

!> A structure with information about a single restart field
type field_restart
  type(vardesc) :: vars         !< Description of a field that is to be read from or written
                                !! to the restart file.
  logical :: mand_var           !< If .true. the run will abort if this field is not successfully
                                !! read from the restart file.
  logical :: initialized        !< .true. if this field has been read from the restart file.
  character(len=32) :: var_name !< A name by which a variable may be queried.
  real    :: conv = 1.0         !< A factor by which a restart field should be multiplied before it
                                !! is written to a restart file, usually to convert it to MKS or
                                !! other standard units [a A-1 ~> 1].  When read, the restart field
                                !! is multiplied by the reciprocal of this factor.
end type field_restart

!> A structure to store information about restart fields that are no longer used
type obsolete_restart
  character(len=32) :: field_name       !< Name of restart field that is no longer in use
  character(len=32) :: replacement_name !< Name of replacement restart field, if applicable
end type obsolete_restart

!> A restart registry and the control structure for restarts
type, public :: MOM_restart_CS ; private
  logical :: initialized = .false. !< True if this control structure has been initialized.
  logical :: restart    !< restart is set to .true. if the run has been started from a full restart
                        !! file.  Otherwise some fields must be initialized approximately.
  integer :: novars = 0 !< The number of restart fields that have been registered.
  integer :: num_obsolete_vars = 0  !< The number of obsolete restart fields that have been registered.
  logical :: parallel_restartfiles  !< If true, the IO layout is used to group processors that write
                                    !! to the same restart file or each processor writes its own
                                    !! (numbered) restart file.  If false, a single restart file is
                                    !! generated after internally combining output from all PEs.
  logical :: new_run                !< If true, the input filenames and restart file existence will
                                    !! result in a new run that is not initialized from restart files.
  logical :: new_run_set = .false.  !< If true, new_run has been determined for this restart_CS.
  logical :: checksum_required      !< If true, require the restart checksums to match and error out otherwise.
                                    !! Users may want to avoid this comparison if for example the restarts are
                                    !! made from a run with a different mask_table than the current run,
                                    !! in which case the checksums will not match and cause crash.
  logical :: symmetric_checksums    !< If true, do the restart checksums on all the edge points for
                                    !! a non-reentrant grid.  Setting this to true requires that
                                    !! SYMMETRIC_MEMORY_ is defined at compile time.
  logical :: unsigned_zeros         !< If true, convert any negative zeros that would be written to
                                    !! the restart file into ordinary unsigned zeros.  This does not
                                    !! change answers, but it can be helpful in comparing restart
                                    !! files after grid rotation, for example.
  logical :: reentrant_x            !< If true, the domain is reentrant in the x-direction.  This is only
                                    !! used here to determine the extent of the restart checksums.
  logical :: reentrant_y            !< If true, the domain is reentrant in the y-direction.  This is only
                                    !! used here to determine the extent of the restart checksums.
  character(len=240) :: restartfile !< The name or name root for MOM restart files.
  integer :: turns                  !< Number of quarter turns from input to model domain
  logical :: locked = .false.       !< If true this registry has been locked and no further restart
                                    !! fields can be added without explicitly unlocking the registry.

  !> An array of descriptions of the registered fields
  type(field_restart), pointer :: restart_field(:) => NULL()

  !> An array of obsolete restart fields
  type(obsolete_restart), pointer :: restart_obsolete(:) => NULL()

  !>@{ Pointers to the fields that have been registered for restarts
  type(p0d), pointer :: var_ptr0d(:) => NULL()
  type(p1d), pointer :: var_ptr1d(:) => NULL()
  type(p2d), pointer :: var_ptr2d(:) => NULL()
  type(p3d), pointer :: var_ptr3d(:) => NULL()
  type(p4d), pointer :: var_ptr4d(:) => NULL()
  !>@}
  integer :: max_fields !< The maximum number of restart fields
end type MOM_restart_CS

!> Register fields for restarts
interface register_restart_field
  module procedure register_restart_field_ptr4d, register_restart_field_4d
  module procedure register_restart_field_ptr3d, register_restart_field_3d
  module procedure register_restart_field_ptr2d, register_restart_field_2d
  module procedure register_restart_field_ptr1d, register_restart_field_1d
  module procedure register_restart_field_ptr0d, register_restart_field_0d
end interface

!> Register a pair of restart fields whose rotations map onto each other
interface register_restart_pair
  module procedure register_restart_pair_ptr2d
  module procedure register_restart_pair_ptr3d
  module procedure register_restart_pair_ptr4d
end interface register_restart_pair

!> Indicate whether a field has been read from a restart file
interface query_initialized
  module procedure query_initialized_name
  module procedure query_initialized_0d, query_initialized_0d_name
  module procedure query_initialized_1d, query_initialized_1d_name
  module procedure query_initialized_2d, query_initialized_2d_name
  module procedure query_initialized_3d, query_initialized_3d_name
  module procedure query_initialized_4d, query_initialized_4d_name
end interface

!> Specify that a field has been initialized, even if it was not read from a restart file
interface set_initialized
  module procedure set_initialized_name, set_initialized_0d_name
  module procedure set_initialized_1d_name, set_initialized_2d_name
  module procedure set_initialized_3d_name, set_initialized_4d_name
end interface

!> Copy the restart variable with the specified name into an array, perhaps after rotation
interface copy_restart_var
  module procedure copy_restart_var_3d
end interface copy_restart_var

!> Copy the restart vector component variables with the specified names into a pair of arrays,
!! perhaps after rotation
interface copy_restart_vector
  module procedure copy_restart_vector_3d
end interface copy_restart_vector

!> Read optional variables from restart files.
interface only_read_from_restarts
  module procedure only_read_restart_field_4d
  module procedure only_read_restart_field_3d
  module procedure only_read_restart_field_2d
!  module procedure only_read_restart_field_1d
!  module procedure only_read_restart_field_0d
  module procedure only_read_restart_pair_3d
end interface


  interface
module subroutine register_restart_field_as_obsolete(field_name, replacement_name, CS)
  character(*), intent(in) :: field_name       !< Name of restart field that is no longer in use
  character(*), intent(in) :: replacement_name !< Name of replacement restart field, if applicable
  type(MOM_restart_CS), intent(inout) :: CS    !< MOM restart control struct

end subroutine register_restart_field_as_obsolete
module subroutine register_restart_field_ptr3d(f_ptr, var_desc, mandatory, CS, conversion)
  real, dimension(:,:,:), &
                      target, intent(in) :: f_ptr     !< A pointer to the field to be read or written
                                                      !! in arbitrary rescaled units [A ~> a]
  type(vardesc),              intent(in) :: var_desc  !< A structure with metadata about this variable
  logical,                    intent(in) :: mandatory !< If true, the run will abort if this field is not
                                                      !! successfully read from the restart file.
  type(MOM_restart_CS),       intent(inout) :: CS     !< MOM restart control struct
  real,             optional, intent(in) :: conversion !< A factor to multiply a restart field by
                                                      !! before it is written [a A-1 ~> 1], 1 by default.

end subroutine register_restart_field_ptr3d
module subroutine register_restart_field_ptr4d(f_ptr, var_desc, mandatory, CS, conversion)
  real, dimension(:,:,:,:), &
                      target, intent(in) :: f_ptr     !< A pointer to the field to be read or written
                                                      !! in arbitrary rescaled units [A ~> a]
  type(vardesc),              intent(in) :: var_desc  !< A structure with metadata about this variable
  logical,                    intent(in) :: mandatory !< If true, the run will abort if this field is not
                                                      !! successfully read from the restart file.
  type(MOM_restart_CS),       intent(inout) :: CS     !< MOM restart control struct
  real,             optional, intent(in) :: conversion !< A factor to multiply a restart field by
                                                      !! before it is written [a A-1 ~> 1], 1 by default.

end subroutine register_restart_field_ptr4d
module subroutine register_restart_field_ptr2d(f_ptr, var_desc, mandatory, CS, conversion)
  real, dimension(:,:), &
                      target, intent(in) :: f_ptr     !< A pointer to the field to be read or written
                                                      !! in arbitrary rescaled units [A ~> a]
  type(vardesc),              intent(in) :: var_desc  !< A structure with metadata about this variable
  logical,                    intent(in) :: mandatory !< If true, the run will abort if this field is not
                                                      !! successfully read from the restart file.
  type(MOM_restart_CS),       intent(inout) :: CS     !< MOM restart control struct
  real,             optional, intent(in) :: conversion !< A factor to multiply a restart field by
                                                      !! before it is written [a A-1 ~> 1], 1 by default.

end subroutine register_restart_field_ptr2d
module subroutine register_restart_field_ptr1d(f_ptr, var_desc, mandatory, CS, conversion)
  real, dimension(:), target, intent(in) :: f_ptr     !< A pointer to the field to be read or written
                                                      !! in arbitrary rescaled units [A ~> a]
  type(vardesc),              intent(in) :: var_desc  !< A structure with metadata about this variable
  logical,                    intent(in) :: mandatory !< If true, the run will abort if this field is not
                                                      !! successfully read from the restart file.
  type(MOM_restart_CS),       intent(inout) :: CS     !< MOM restart control struct
  real,             optional, intent(in) :: conversion !< A factor to multiply a restart field by
                                                      !! before it is written [a A-1 ~> 1], 1 by default.

end subroutine register_restart_field_ptr1d
module subroutine register_restart_field_ptr0d(f_ptr, var_desc, mandatory, CS, conversion)
  real,               target, intent(in) :: f_ptr     !< A pointer to the field to be read or written
                                                      !! in arbitrary rescaled units [A ~> a]
  type(vardesc),              intent(in) :: var_desc  !< A structure with metadata about this variable
  logical,                    intent(in) :: mandatory !< If true, the run will abort if this field is not
                                                      !! successfully read from the restart file.
  type(MOM_restart_CS),       intent(inout) :: CS     !< MOM restart control struct
  real,             optional, intent(in) :: conversion !< A factor to multiply a restart field by
                                                      !! before it is written [a A-1 ~> 1], 1 by default.

end subroutine register_restart_field_ptr0d
module subroutine register_restart_pair_ptr2d(a_ptr, b_ptr, a_desc, b_desc, &
                mandatory, CS, conversion, scalar_pair)
  real, dimension(:,:), target, intent(in) :: a_ptr   !< First field pointer
                                                      !! in arbitrary rescaled units [A ~> a]
  real, dimension(:,:), target, intent(in) :: b_ptr   !< Second field pointer
                                                      !! in arbitrary rescaled units [A ~> a]
  type(vardesc),                intent(in) :: a_desc  !< First field descriptor
  type(vardesc),                intent(in) :: b_desc  !< Second field descriptor
  logical,                      intent(in) :: mandatory !< If true, abort if field is missing
  type(MOM_restart_CS),      intent(inout) :: CS      !< MOM restart control structure
  real,               optional, intent(in) :: conversion !< A factor to multiply a restart field by
                                                      !! before it is written [a A-1 ~> 1], 1 by default.
  logical,            optional, intent(in) :: scalar_pair !< If true, the arrays describe a pair of
                                                      !! scalars, instead of vector components
                                                      !! whose signs change when rotated

  ! Local variables
                          ! including sign changes to account for grid rotation [a A-1 ~> 1]

end subroutine register_restart_pair_ptr2d
module subroutine register_restart_pair_ptr3d(a_ptr, b_ptr, a_desc, b_desc, &
                mandatory, CS, conversion, scalar_pair)
  real, dimension(:,:,:), target, intent(in) :: a_ptr !< First field pointer
                                                      !! in arbitrary rescaled units [A ~> a]
  real, dimension(:,:,:), target, intent(in) :: b_ptr !< Second field pointer
                                                      !! in arbitrary rescaled units [A ~> a]
  type(vardesc),                intent(in) :: a_desc  !< First field descriptor
  type(vardesc),                intent(in) :: b_desc  !< Second field descriptor
  logical,                      intent(in) :: mandatory !< If true, abort if field is missing
  type(MOM_restart_CS),      intent(inout) :: CS      !< MOM restart control structure
  real,               optional, intent(in) :: conversion !< A factor to multiply a restart field by
                                                      !! before it is written [a A-1 ~> 1], 1 by default.
  logical,            optional, intent(in) :: scalar_pair !< If true, the arrays describe a pair of
                                                      !! scalars, instead of vector components
                                                      !! whose signs change when rotated

  ! Local variables
                          ! including sign changes to account for grid rotation [a A-1 ~> 1]

end subroutine register_restart_pair_ptr3d
module subroutine register_restart_pair_ptr4d(a_ptr, b_ptr, a_desc, b_desc, &
                mandatory, CS, conversion, scalar_pair)
  real, dimension(:,:,:,:), target, intent(in) :: a_ptr !< First field pointer
                                                      !! in arbitrary rescaled units [A ~> a]
  real, dimension(:,:,:,:), target, intent(in) :: b_ptr !< Second field pointer
                                                      !! in arbitrary rescaled units [A ~> a]
  type(vardesc),                intent(in) :: a_desc  !< First field descriptor
  type(vardesc),                intent(in) :: b_desc  !< Second field descriptor
  logical,                      intent(in) :: mandatory !< If true, abort if field is missing
  type(MOM_restart_CS),      intent(inout) :: CS      !< MOM restart control structure
  real,               optional, intent(in) :: conversion !< A factor to multiply a restart field by
                                                      !! before it is written [a A-1 ~> 1], 1 by default.
  logical,            optional, intent(in) :: scalar_pair !< If true, the arrays describe a pair of
                                                      !! scalars, instead of vector components
                                                      !! whose signs change when rotated

  ! Local variables
                          ! including sign changes to account for grid rotation [a A-1 ~> 1]

end subroutine register_restart_pair_ptr4d
module subroutine set_conversion_pair(u_conv, v_conv, turns, conversion, scalar_pair)
  real,   intent(out) :: u_conv !< A factor to multiply the u-component of a vector by before it is
                                !! written, including sign changes due to grid rotation [a A-1 ~> 1]
  real,   intent(out) :: v_conv !< A factor to multiply the u-component of a vector by before it is
                                !! written, including sign changes due to grid rotation [a A-1 ~> 1]
  integer, intent(in) :: turns  !< Number of quarter turns from input to model domain
  real,    optional, intent(in) :: conversion !< A factor to multiply a restart field by
                                !! before it is written [a A-1 ~> 1], 1 by default.
  logical, optional, intent(in) :: scalar_pair !< If true, the arrays describe a pair of scalars,
                                 !! instead of vector components whose signs change when rotated

  ! Local variables

end subroutine set_conversion_pair
module subroutine register_restart_field_4d(f_ptr, name, mandatory, CS, longname, units, conversion, &
                                     hor_grid, z_grid, t_grid, extra_axes)
  real, dimension(:,:,:,:), &
                      target, intent(in) :: f_ptr     !< A pointer to the field to be read or written
                                                      !! in arbitrary rescaled units [A ~> a]
  character(len=*),           intent(in) :: name      !< variable name to be used in the restart file
  logical,                    intent(in) :: mandatory !< If true, the run will abort if this field is not
                                                      !! successfully read from the restart file.
  type(MOM_restart_CS),       intent(inout) :: CS     !< MOM restart control struct
  character(len=*), optional, intent(in) :: longname  !< variable long name
  character(len=*), optional, intent(in) :: units     !< variable units
  real,             optional, intent(in) :: conversion !< A factor to multiply a restart field by
                                                      !! before it is written [a A-1 ~> 1], 1 by default.
  character(len=*), optional, intent(in) :: hor_grid  !< variable horizontal staggering, 'h' if absent
  character(len=*), optional, intent(in) :: z_grid    !< variable vertical staggering, 'L' if absent
  character(len=*), optional, intent(in) :: t_grid    !< time description: s, p, or 1, 's' if absent
  type(axis_info),  dimension(:), &
                    optional, intent(in) :: extra_axes !< dimensions other than space-time


  ! first 2 dimensions in dim_names are reserved for i,j
  ! so extra_dimensions are shifted to index 3.
  ! this is designed not to break the behavior in SIS2
  ! (see register_restart_field_4d in SIS_restart.F90)
end subroutine register_restart_field_4d
module subroutine register_restart_field_3d(f_ptr, name, mandatory, CS, longname, units, conversion, &
                                     hor_grid, z_grid, t_grid, extra_axes)
  real, dimension(:,:,:), &
                      target, intent(in) :: f_ptr     !< A pointer to the field to be read or written
                                                      !! in arbitrary rescaled units [A ~> a]
  character(len=*),           intent(in) :: name      !< variable name to be used in the restart file
  logical,                    intent(in) :: mandatory !< If true, the run will abort if this field is not
                                                      !! successfully read from the restart file.
  type(MOM_restart_CS),       intent(inout) :: CS     !< MOM restart control struct
  character(len=*), optional, intent(in) :: longname  !< variable long name
  character(len=*), optional, intent(in) :: units     !< variable units
  real,             optional, intent(in) :: conversion !< A factor to multiply a restart field by
                                                      !! before it is written [a A-1 ~> 1], 1 by default.
  character(len=*), optional, intent(in) :: hor_grid  !< variable horizontal staggering, 'h' if absent
  character(len=*), optional, intent(in) :: z_grid    !< variable vertical staggering, 'L' if absent
  character(len=*), optional, intent(in) :: t_grid    !< time description: s, p, or 1, 's' if absent
  type(axis_info),  dimension(:), &
                    optional, intent(in) :: extra_axes !< dimensions other than space-time


  ! first 2 dimensions in dim_names are reserved for i,j
  ! so extra_dimensions are shifted to index 3.
  ! this is designed not to break the behavior in SIS2
  ! (see register_restart_field_4d in SIS_restart.F90)
end subroutine register_restart_field_3d
module subroutine register_restart_field_2d(f_ptr, name, mandatory, CS, longname, units, conversion, &
                                     hor_grid, z_grid, t_grid)
  real, dimension(:,:), &
                      target, intent(in) :: f_ptr     !< A pointer to the field to be read or written
                                                      !! in arbitrary rescaled units [A ~> a]
  character(len=*),           intent(in) :: name      !< variable name to be used in the restart file
  logical,                    intent(in) :: mandatory !< If true, the run will abort if this field is not
                                                      !! successfully read from the restart file.
  type(MOM_restart_CS),       intent(inout) :: CS     !< MOM restart control struct
  character(len=*), optional, intent(in) :: longname  !< variable long name
  character(len=*), optional, intent(in) :: units     !< variable units
  real,             optional, intent(in) :: conversion !< A factor to multiply a restart field by
                                                      !! before it is written [a A-1 ~> 1], 1 by default.
  character(len=*), optional, intent(in) :: hor_grid  !< variable horizontal staggering, 'h' if absent
  character(len=*), optional, intent(in) :: z_grid    !< variable vertical staggering, '1' if absent
  character(len=*), optional, intent(in) :: t_grid    !< time description: s, p, or 1, 's' if absent


end subroutine register_restart_field_2d
module subroutine register_restart_field_1d(f_ptr, name, mandatory, CS, longname, units, conversion, &
                                     hor_grid, z_grid, t_grid)
  real, dimension(:), target, intent(in) :: f_ptr     !< A pointer to the field to be read or written
                                                      !! in arbitrary rescaled units [A ~> a]
  character(len=*),           intent(in) :: name      !< variable name to be used in the restart file
  logical,                    intent(in) :: mandatory !< If true, the run will abort if this field is not
                                                      !! successfully read from the restart file.
  type(MOM_restart_CS),       intent(inout) :: CS     !< MOM restart control struct
  character(len=*), optional, intent(in) :: longname  !< variable long name
  character(len=*), optional, intent(in) :: units     !< variable units
  real,             optional, intent(in) :: conversion !< A factor to multiply a restart field by
                                                      !! before it is written [a A-1 ~> 1], 1 by default.
  character(len=*), optional, intent(in) :: hor_grid  !< variable horizontal staggering, '1' if absent
  character(len=*), optional, intent(in) :: z_grid    !< variable vertical staggering, 'L' if absent
  character(len=*), optional, intent(in) :: t_grid    !< time description: s, p, or 1, 's' if absent


end subroutine register_restart_field_1d
module subroutine register_restart_field_0d(f_ptr, name, mandatory, CS, longname, units, conversion, &
                                     t_grid)
  real,               target, intent(in) :: f_ptr     !< A pointer to the field to be read or written
                                                      !! in arbitrary rescaled units [A ~> a]
  character(len=*),           intent(in) :: name      !< variable name to be used in the restart file
  logical,                    intent(in) :: mandatory !< If true, the run will abort if this field is not
                                                      !! successfully read from the restart file.
  type(MOM_restart_CS),       intent(inout) :: CS     !< MOM restart control struct
  character(len=*), optional, intent(in) :: longname  !< variable long name
  character(len=*), optional, intent(in) :: units     !< variable units
  real,             optional, intent(in) :: conversion !< A factor to multiply a restart field by
                                                      !! before it is written [a A-1 ~> 1], 1 by default.
  character(len=*), optional, intent(in) :: t_grid    !< time description: s, p, or 1, 's' if absent


end subroutine register_restart_field_0d
module function query_initialized_name(name, CS) result(query_initialized)
  character(len=*),     intent(in) :: name  !< The name of the field that is being queried
  type(MOM_restart_CS), intent(in) :: CS    !< MOM restart control struct
  logical :: query_initialized


end function query_initialized_name
module function query_initialized_0d(f_ptr, CS) result(query_initialized)
  real,         target, intent(in) :: f_ptr !< A pointer to the field that is being queried [arbitrary]
  type(MOM_restart_CS), intent(in) :: CS    !< MOM restart control struct
  logical :: query_initialized


end function query_initialized_0d
module function query_initialized_1d(f_ptr, CS) result(query_initialized)
  real, dimension(:), target, intent(in) :: f_ptr !< A pointer to the field that is being queried [arbitrary]
  type(MOM_restart_CS),       intent(in) :: CS    !< MOM restart control struct
  logical :: query_initialized


end function query_initialized_1d
module function query_initialized_2d(f_ptr, CS) result(query_initialized)
  real, dimension(:,:), &
                target, intent(in) :: f_ptr !< A pointer to the field that is being queried [arbitrary]
  type(MOM_restart_CS), intent(in) :: CS    !< MOM restart control struct
  logical :: query_initialized


end function query_initialized_2d
module function query_initialized_3d(f_ptr, CS) result(query_initialized)
  real, dimension(:,:,:), &
                target, intent(in) :: f_ptr !< A pointer to the field that is being queried [arbitrary]
  type(MOM_restart_CS), intent(in) :: CS    !< MOM restart control struct
  logical :: query_initialized


end function query_initialized_3d
module function query_initialized_4d(f_ptr, CS) result(query_initialized)
  real, dimension(:,:,:,:),  &
                target, intent(in) :: f_ptr !< A pointer to the field that is being queried [arbitrary]
  type(MOM_restart_CS), intent(in) :: CS    !< MOM restart control struct
  logical :: query_initialized


end function query_initialized_4d
module function query_initialized_0d_name(f_ptr, name, CS) result(query_initialized)
  real,         target, intent(in) :: f_ptr !< The field that is being queried [arbitrary]
  character(len=*),     intent(in) :: name  !< The name of the field that is being queried
  type(MOM_restart_CS), intent(in) :: CS    !< MOM restart control struct
  logical :: query_initialized


end function query_initialized_0d_name
module function query_initialized_1d_name(f_ptr, name, CS) result(query_initialized)
  real, dimension(:),  &
                target, intent(in) :: f_ptr !< The field that is being queried [arbitrary]
  character(len=*),     intent(in) :: name  !< The name of the field that is being queried
  type(MOM_restart_CS), intent(in) :: CS    !< MOM restart control struct
  logical :: query_initialized


end function query_initialized_1d_name
module function query_initialized_2d_name(f_ptr, name, CS) result(query_initialized)
  real, dimension(:,:),  &
                target, intent(in) :: f_ptr !< The field that is being queried [arbitrary]
  character(len=*),     intent(in) :: name  !< The name of the field that is being queried
  type(MOM_restart_CS), intent(in) :: CS    !< MOM restart control struct
  logical :: query_initialized


end function query_initialized_2d_name
module function query_initialized_3d_name(f_ptr, name, CS) result(query_initialized)
  real, dimension(:,:,:),  &
                target, intent(in) :: f_ptr !< The field that is being queried [arbitrary]
  character(len=*),     intent(in) :: name  !< The name of the field that is being queried
  type(MOM_restart_CS), intent(in) :: CS    !< MOM restart control struct
  logical :: query_initialized


end function query_initialized_3d_name
module function query_initialized_4d_name(f_ptr, name, CS) result(query_initialized)
  real, dimension(:,:,:,:),  &
                target, intent(in) :: f_ptr !< The field that is being queried [arbitrary]
  character(len=*),     intent(in) :: name  !< The name of the field that is being queried
  type(MOM_restart_CS), intent(in) :: CS    !< MOM restart control struct
  logical :: query_initialized


end function query_initialized_4d_name
module subroutine set_initialized_name(name, CS)
  character(len=*),     intent(in)    :: name  !< The name of the field that is being set
  type(MOM_restart_CS), intent(inout) :: CS    !< MOM restart control struct


end subroutine set_initialized_name
module subroutine set_initialized_0d_name(f_ptr, name, CS)
  real,         target, intent(in)    :: f_ptr !< The variable that has been initialized [arbitrary]
  character(len=*),     intent(in)    :: name  !< The name of the field that has been initialized
  type(MOM_restart_CS), intent(inout) :: CS    !< MOM restart control struct


end subroutine set_initialized_0d_name
module subroutine set_initialized_1d_name(f_ptr, name, CS)
  real, dimension(:),  &
                target, intent(in)    :: f_ptr !< The array that has been initialized [arbitrary]
  character(len=*),     intent(in)    :: name  !< The name of the field that has been initialized
  type(MOM_restart_CS), intent(inout) :: CS    !< MOM restart control struct


end subroutine set_initialized_1d_name
module subroutine set_initialized_2d_name(f_ptr, name, CS)
  real, dimension(:,:),  &
                target, intent(in)    :: f_ptr !< The array that has been initialized [arbitrary]
  character(len=*),     intent(in)    :: name  !< The name of the field that has been initialized
  type(MOM_restart_CS), intent(inout) :: CS    !< MOM restart control struct


end subroutine set_initialized_2d_name
module subroutine set_initialized_3d_name(f_ptr, name, CS)
  real, dimension(:,:,:),  &
                target, intent(in)    :: f_ptr !< The array that has been initialized [arbitrary]
  character(len=*),     intent(in)    :: name  !< The name of the field that has been initialized
  type(MOM_restart_CS), intent(inout) :: CS    !< MOM restart control struct


end subroutine set_initialized_3d_name
module subroutine set_initialized_4d_name(f_ptr, name, CS)
  real, dimension(:,:,:,:),  &
                target, intent(in)    :: f_ptr !< The array that has been initialized [arbitrary]
  character(len=*),     intent(in)    :: name  !< The name of the field that has been initialized
  type(MOM_restart_CS), intent(inout) :: CS    !< MOM restart control struct


end subroutine set_initialized_4d_name
module subroutine only_read_restart_field_4d(varname, f_ptr, G, CS, position, filename, directory, success, scale)
  character(len=*),                intent(in)    :: varname   !< The variable name to be used in the restart file
  real, dimension(:,:,:,:),        intent(inout) :: f_ptr     !< The array for the field to be read
                                                              !! in arbitrary rescaled units [A ~> a]
  type(ocean_grid_type),           intent(in)    :: G         !< The ocean's grid structure
  type(MOM_restart_CS),            intent(in)    :: CS        !< MOM restart control struct
  integer,               optional, intent(in)    :: position  !< A coded integer indicating the horizontal
                                                              !! position of this variable
  character(len=*),      optional, intent(in)    :: filename  !< The list of restart file names or a single
                                                              !! character 'r' to read automatically named files
  character(len=*),      optional, intent(in)    :: directory !< The directory in which to seek restart files.
  logical,               optional, intent(out)   :: success   !< True if the field was read successfully
  real,                  optional, intent(in)    :: scale     !< A factor by which the field will be scaled
                                                              !! [A a-1 ~> 1] to convert from the units in
                                                              !! the file to the internal units of this field

  ! Local variables

end subroutine only_read_restart_field_4d
module subroutine only_read_restart_field_3d(varname, f_ptr, G, CS, position, filename, directory, success, scale)
  character(len=*),                intent(in)    :: varname   !< The variable name to be used in the restart file
  real, dimension(:,:,:),          intent(inout) :: f_ptr     !< The array for the field to be read
                                                              !! in arbitrary rescaled units [A ~> a]
  type(ocean_grid_type),           intent(in)    :: G         !< The ocean's grid structure
  type(MOM_restart_CS),            intent(in)    :: CS        !< MOM restart control struct
  integer,               optional, intent(in)    :: position  !< A coded integer indicating the horizontal
                                                              !! position of this variable
  character(len=*),      optional, intent(in)    :: filename  !< The list of restart file names or a single
                                                              !! character 'r' to read automatically named files
  character(len=*),      optional, intent(in)    :: directory !< The directory in which to seek restart files.
  logical,               optional, intent(out)   :: success   !< True if the field was read successfully
  real,                  optional, intent(in)    :: scale     !< A factor by which the field will be scaled
                                                              !! [A a-1 ~> 1] to convert from the units in
                                                              !! the file to the internal units of this field

  ! Local variables

end subroutine only_read_restart_field_3d
module subroutine only_read_restart_field_2d(varname, f_ptr, G, CS, position, filename, directory, success, scale)
  character(len=*),                intent(in)    :: varname   !< The variable name to be used in the restart file
  real, dimension(:,:),            intent(inout) :: f_ptr     !< The array for the field to be read
                                                              !! in arbitrary rescaled units [A ~> a]
  type(ocean_grid_type),           intent(in)    :: G         !< The ocean's grid structure
  type(MOM_restart_CS),            intent(in)    :: CS        !< MOM restart control struct
  integer,               optional, intent(in)    :: position  !< A coded integer indicating the horizontal
                                                              !! position of this variable
  character(len=*),      optional, intent(in)    :: filename  !< The list of restart file names or a single
                                                              !! character 'r' to read automatically named files
  character(len=*),      optional, intent(in)    :: directory !< The directory in which to seek restart files.
  logical,               optional, intent(out)   :: success   !< True if the field was read successfully
  real,                  optional, intent(in)    :: scale     !< A factor by which the field will be scaled
                                                              !! [A a-1 ~> 1] to convert from the units in
                                                              !! the file to the internal units of this field

  ! Local variables

end subroutine only_read_restart_field_2d
module subroutine only_read_restart_pair_3d(a_ptr, b_ptr, a_name, b_name, G, CS, &
                                     stagger, filename, directory, success, scale)
  real, dimension(:,:,:),          intent(inout) :: a_ptr     !< The array for the first field to be read
                                                              !! in arbitrary rescaled units [A ~> a]
  real, dimension(:,:,:),          intent(inout) :: b_ptr     !< The array for the second field to be read
                                                              !! in arbitrary rescaled units [A ~> a]
  character(len=*),                intent(in)    :: a_name    !< The first variable name to be used in the restart file
  character(len=*),                intent(in)    :: b_name    !< The second variable name to be used in the restart file
  type(ocean_grid_type),           intent(in)    :: G         !< The ocean's grid structure
  type(MOM_restart_CS),            intent(in)    :: CS        !< MOM restart control struct
  integer,               optional, intent(in)    :: stagger   !< A coded integer indicating the horizontal
                                                              !! position of this pair of variables
  character(len=*),      optional, intent(in)    :: filename  !< The list of restart file names or a single
                                                              !! character 'r' to read automatically named files
  character(len=*),      optional, intent(in)    :: directory !< The directory in which to seek restart files.
  logical,               optional, intent(out)   :: success   !< True if the field was read successfully
  real,                  optional, intent(in)    :: scale     !< A factor by which the fields will be scaled
                                                              !! [A a-1 ~> 1] to convert from the units in
                                                              !! the file to the internal units of this field

  ! Local variables

end subroutine only_read_restart_pair_3d
module function find_var_in_restart_files(varname, G, CS, file_path, filename, directory, is_global) result (found)
  character(len=*),                intent(in)    :: varname   !< The variable name to be used in the restart file
  type(ocean_grid_type),           intent(in)    :: G         !< The ocean's grid structure
  type(MOM_restart_CS),            intent(in)    :: CS        !< MOM restart control struct
  character(len=:),   allocatable, intent(out)   :: file_path !< The full path to the file in which the
                                                              !! variable is found
  character(len=*),      optional, intent(in)    :: filename  !< The list of restart file names or a single
                                                              !! character 'r' to read automatically named files
  character(len=*),      optional, intent(in)    :: directory !< The directory in which to seek restart files.
  logical,               optional, intent(out)   :: is_global !< True if the file is global.
  logical :: found !< True if the named variable was found in the restart files.

  ! Local variables

end function find_var_in_restart_files
module subroutine copy_restart_var_3d(var, name, CS, unrotate)
  real, dimension(:,:,:), intent(inout) :: var   !< The field that is being copied [arbitrary]
  character(len=*),       intent(in)    :: name  !< The name of the field that is being copied
  type(MOM_restart_CS),   intent(in)    :: CS    !< MOM restart control struct
  logical, optional,      intent(in)    :: unrotate !< If present and true, the output is on an unrotated grid.


end subroutine copy_restart_var_3d
module subroutine copy_restart_vector_3d(u_var, v_var, u_name, v_name, CS, unrotate, scalar_pair)
  real, dimension(:,:,:), intent(inout) :: u_var !< The u-component of the field that is being copied [arbitrary]
  real, dimension(:,:,:), intent(inout) :: v_var !< The u-component of the field that is being copied [arbitrary]
  character(len=*),       intent(in)    :: u_name !< The name of the u-component of the field that is being copied
  character(len=*),       intent(in)    :: v_name !< The name of the v-component of the field that is being copied
  type(MOM_restart_CS),   intent(in)    :: CS    !< MOM restart control struct
  logical, optional,      intent(in)    :: unrotate !< If present and true, the output is on an unrotated grid.
  logical,      optional, intent(in)    :: scalar_pair !< If true, the arrays describe a pair of
                                                 !! scalars, instead of vector components
                                                 !! whose signs change when rotated


end subroutine copy_restart_vector_3d
logical module function size_mismatch_3d(var_a, var_b, turns, size_msg)
  real,    intent(in) :: var_a(:,:,:)   !< The first field being compared
  real,    intent(in) :: var_b(:,:,:)   !< The second field being compared
  integer, intent(in) :: turns          !< Number of quarter turns from input to model domain
  character(len=256), intent(out) :: size_msg  !< The array sizes

end function size_mismatch_3d
module subroutine save_restart(directory, time, G, CS, time_stamped, filename, GV, num_rest_files, write_IC)
  character(len=*),        intent(in)    :: directory !< The directory where the restart files
                                                  !! are to be written
  type(time_type),         intent(in)    :: time  !< The current model time
  type(ocean_grid_type),   intent(inout) :: G     !< The ocean's grid structure as seen from the driver.
  type(MOM_restart_CS),    intent(inout) :: CS    !< MOM restart control struct
  logical,       optional, intent(in)    :: time_stamped !< If present and true, add time-stamp
                                                  !! to the restart file names
  character(len=*), optional, intent(in) :: filename !< A filename that overrides the name in CS%restartfile
  type(verticalGrid_type), &
                 optional, intent(in)    :: GV    !< The ocean's vertical grid structure
  integer,       optional, intent(out)   :: num_rest_files !< number of restart files written
  logical,       optional, intent(in)    :: write_IC !< If present and true, initial conditions
                                                  !! are being written

  ! Local variables
                                        ! are to be read from the restart file.
                                        ! each variable that will be written.
                                        ! to the name of files after the first.
                                        ! and the variables already in a file.
                                        ! starting position of each variable in a file's record,
                                        ! based on the use of NetCDF 3.6 or later.  For earlier
                                        ! versions of NetCDF, the value was 2147483647_int64.
                                        ! current and next files.

end subroutine save_restart
module subroutine restore_state(filename, directory, day, G, CS)
  character(len=*),      intent(in)  :: filename  !< The list of restart file names or a single
                                                  !! character 'r' to read automatically named files
  character(len=*),      intent(in)  :: directory !< The directory in which to find restart files
  type(time_type),       intent(out) :: day       !< The time of the restarted run
  type(ocean_grid_type), intent(in)  :: G         !< The ocean's grid structure
  type(MOM_restart_CS),  intent(inout) :: CS      !< MOM restart control struct

  ! Local variables
                 ! from the units in the file to the internal units of this field
                             ! explicitly in filename) that are open.



end subroutine restore_state
module function restart_files_exist(filename, directory, G, CS)
  character(len=*),      intent(in)  :: filename  !< The list of restart file names or a single
                                                  !! character 'r' to read automatically named files
  character(len=*),      intent(in)  :: directory !< The directory in which to find restart files
  type(ocean_grid_type), intent(in)  :: G         !< The ocean's grid structure
  type(MOM_restart_CS),  intent(in)  :: CS        !< MOM restart control struct
  logical :: restart_files_exist                  !< The function result, which indicates whether
                                                  !! any of the explicitly or automatically named
                                                  !! restart files exist in directory

end function restart_files_exist
module function determine_is_new_run(filename, directory, G, CS) result(is_new_run)
  character(len=*),      intent(in)  :: filename  !< The list of restart file names or a single
                                                  !! character 'r' to read automatically named files
  character(len=*),      intent(in)  :: directory !< The directory in which to find restart files
  type(ocean_grid_type), intent(in)  :: G         !< The ocean's grid structure
  type(MOM_restart_CS),  intent(inout) :: CS      !< MOM restart control struct
  logical :: is_new_run                           !< The function result, which indicates whether
                                                  !! this is a new run, based on the value of
                                                  !! filename and whether restart files exist

end function determine_is_new_run
module function is_new_run(CS)
  type(MOM_restart_CS), intent(in) :: CS  !< MOM restart control struct

  logical :: is_new_run                !< The function result, which had been stored in CS during
                                       !! a previous call to determine_is_new_run

end function is_new_run
module function open_restart_units(filename, directory, G, CS, IO_handles, file_paths, &
                            global_files) result(num_files)
  character(len=*),      intent(in)  :: filename  !< The list of restart file names or a single
                                                  !! character 'r' to read automatically named files
  character(len=*),      intent(in)  :: directory !< The directory in which to find restart files
  type(ocean_grid_type), intent(in)  :: G         !< The ocean's grid structure
  type(MOM_restart_CS),  intent(in)  :: CS        !< MOM restart control struct

  type(MOM_infra_file), dimension(:), &
               optional, intent(out) :: IO_handles !< The I/O handles of all opened files
  character(len=*), dimension(:), &
               optional, intent(out) :: file_paths   !< The full paths to open files
  logical, dimension(:), &
               optional, intent(out) :: global_files !< True if a file is global

  integer :: num_files  !< The number of files (both automatically named restart
                        !! files and others explicitly in filename) that have been opened.

  ! Local variables
                                  ! additional restart files.
                             ! been opened using their numbered suffix.
                             ! current file name.

end function open_restart_units
module function get_num_restart_files(filenames, directory, G, CS, file_paths) result(num_files)
  character(len=*),      intent(in)  :: filenames !< The list of restart file names or a single
                                                  !! character 'r' to read automatically named files
  character(len=*),      intent(in)  :: directory !< The directory in which to find restart files
  type(ocean_grid_type), intent(in)  :: G         !< The ocean's grid structure
  type(MOM_restart_CS),  intent(in)  :: CS        !< MOM restart control struct
  character(len=*), dimension(:), &
               optional, intent(out) :: file_paths !< The full paths to the restart files.

  integer :: num_files  !< The function result, the number of files (both automatically named
                        !! restart files and others explicitly in filename) that have been opened

end function get_num_restart_files
module subroutine restart_init(param_file, CS, restart_root)
  type(param_file_type), intent(in) :: param_file !< A structure to parse for run-time parameters
  type(MOM_restart_CS),  pointer    :: CS !< A pointer to a MOM_restart_CS object that is allocated here
  character(len=*), optional, &
                         intent(in) :: restart_root !< A filename root that overrides the value
                                          !! set by RESTARTFILE to enable the use of this module by
                                          !! other components than MOM.


  ! This include declares and sets the variable "version".

end subroutine restart_init
module subroutine lock_check(CS, var_desc, name)
  type(MOM_restart_CS),       intent(in) :: CS        !< A MOM_restart_CS object (intent in)
  type(vardesc),    optional, intent(in) :: var_desc  !< A structure with metadata about this variable
  character(len=*), optional, intent(in) :: name      !< variable name to be used in the restart file


end subroutine lock_check
module subroutine restart_registry_lock(CS, unlocked)
  type(MOM_restart_CS), intent(inout) :: CS        !< A MOM_restart_CS object (intent inout)
  logical, optional,    intent(in)    :: unlocked  !< If present and true, unlock the registry

end subroutine restart_registry_lock
module subroutine restart_init_end(CS)
  type(MOM_restart_CS),  pointer    :: CS !< A pointer to a MOM_restart_CS object

end subroutine restart_init_end
module subroutine restart_end(CS)
  type(MOM_restart_CS),  pointer    :: CS !< A pointer to a MOM_restart_CS object

end subroutine restart_end
module subroutine restart_error(CS)
  type(MOM_restart_CS),  intent(in) :: CS   !< MOM restart control struct


end subroutine restart_error
module subroutine get_checksum_loop_ranges(G, CS, pos, isL, ieL, jsL, jeL)
  type(ocean_grid_type), intent(in)  :: G   !< The ocean's grid structure
  type(MOM_restart_CS),  intent(in)  :: CS  !< MOM restart control structure
  integer,               intent(in)  :: pos !< A coded integer indicating the horizontal staggering
                                            !! of a variable
  integer,               intent(out) :: isL !< i-start for checksum
  integer,               intent(out) :: ieL !< i-end for checksum
  integer,               intent(out) :: jsL !< j-start for checksum
  integer,               intent(out) :: jeL !< j-end for checksum

  ! Regular non-symmetric compute domain
end subroutine get_checksum_loop_ranges
module function get_variable_byte_size(pos, z_grid, t_grid, G, num_z) result(var_sz)
  integer,               intent(in) :: pos      !< An integer indicating the horizontal staggering position
  character(len=8),      intent(in) :: z_grid   !< The vertical grid string to interpret
  character(len=8),      intent(in) :: t_grid   !< A time string to interpret
  type(ocean_grid_type), intent(in) :: G        !< The ocean's grid structure
  integer,               intent(in) :: num_z    !< The number of vertical layers in the grid
  integer(kind=int64) :: var_sz !< The function result, the size in bytes of a variable

  ! Local variables

end function get_variable_byte_size
  end interface

end module MOM_restart
