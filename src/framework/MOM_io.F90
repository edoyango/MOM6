! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> This module contains I/O framework code
module MOM_io

use MOM_array_transform,  only : allocate_rotated_array, rotate_array
use MOM_array_transform,  only : rotate_array_pair, rotate_vector
use MOM_domains,          only : MOM_domain_type, domain1D, broadcast, get_domain_components
use MOM_domains,          only : rescale_comp_data, num_PEs, AGRID, BGRID_NE, CGRID_NE
use MOM_dyn_horgrid,      only : dyn_horgrid_type
use MOM_ensemble_manager, only : get_ensemble_id
use MOM_error_handler,    only : MOM_error, NOTE, FATAL, WARNING, is_root_PE
use MOM_file_parser,      only : log_version, param_file_type
use MOM_grid,             only : ocean_grid_type
use MOM_io_infra,         only : read_field, read_vector
use MOM_io_infra,         only : read_data => read_field ! Deprecated
use MOM_io_infra,         only : read_field_chksum
use MOM_io_infra,         only : file_exists
use MOM_io_infra,         only : open_ASCII_file, close_file, file_is_open
use MOM_io_infra,         only : get_field_size, field_exists, get_field_atts
use MOM_io_infra,         only : get_axis_data, get_filename_suffix
use MOM_io_infra,         only : write_version
use MOM_io_infra,         only : MOM_namelist_file, check_namelist_error, io_infra_init, io_infra_end
use MOM_io_infra,         only : APPEND_FILE, ASCII_FILE, MULTIPLE, NETCDF_FILE, OVERWRITE_FILE
use MOM_io_infra,         only : READONLY_FILE, SINGLE_FILE, WRITEONLY_FILE
use MOM_io_infra,         only : CENTER, CORNER, NORTH_FACE, EAST_FACE
use MOM_io_file,          only : MOM_file, MOM_infra_file, MOM_netcdf_file
use MOM_io_file,          only : MOM_axis, MOM_field
use MOM_string_functions, only : lowercase, slasher
use MOM_verticalGrid,     only : verticalGrid_type

use iso_fortran_env,      only : int32, int64, stdout_iso=>output_unit, stderr_iso=>error_unit
use netcdf,               only : NF90_open, NF90_inq_varid, NF90_inq_varids, NF90_inquire, NF90_close
use netcdf,               only : NF90_inquire_variable, NF90_get_var, NF90_get_att, NF90_inquire_attribute
use netcdf,               only : NF90_strerror, NF90_inquire_dimension
use netcdf,               only : NF90_NOWRITE, NF90_NOERR, NF90_GLOBAL, NF90_ENOTATT, NF90_CHAR

! The following are not used in MOM6, but may be used by externals (e.g. SIS2).
use MOM_io_infra, only : axistype   ! still used but soon to be nuked
use MOM_io_infra, only : fieldtype
use MOM_io_infra, only : file_type
use MOM_io_infra, only : get_file_info
use MOM_io_infra, only : get_file_fields
use MOM_io_infra, only : get_file_times
use MOM_io_infra, only : open_file
use MOM_io_infra, only : write_field

implicit none ; private

! These interfaces are actually implemented in this file.
public :: create_MOM_file, reopen_MOM_file, cmor_long_std, ensembler, MOM_io_init
public :: MOM_field
public :: MOM_write_field, var_desc, modify_vardesc, query_vardesc, position_from_horgrid
public :: open_namelist_file, check_namelist_error, check_nml_error
public :: get_var_sizes, verify_variable_units, num_timelevels, read_variable, read_attribute
public :: open_file_to_read, close_file_to_read
! The following are simple pass throughs of routines from MOM_io_infra or other modules.
public :: file_exists, open_ASCII_file, close_file
public :: MOM_file, MOM_infra_file, MOM_netcdf_file
public :: field_exists, get_filename_appendix
public :: fieldtype, field_size, get_field_atts
public :: axistype, get_axis_data
public :: MOM_read_data, MOM_read_vector, read_field_chksum
public :: read_netCDF_data
public :: slasher, write_version_number
public :: io_infra_init, io_infra_end
public :: stdout_if_root
public :: get_var_axes_info
public :: get_axis_info
! This is used to set up information descibing non-domain-decomposed axes.
public :: axis_info, set_axis_info, delete_axis_info
! This is used to set up global file attributes
public :: attribute_info, set_attribute_info, delete_attribute_info
! This API is here just to support potential use by non-FMS drivers, and should not persist.
public :: read_data
!> These encoding constants are used to indicate the file format
public :: ASCII_FILE, NETCDF_FILE
!> These encoding constants are used to indicate whether the file is domain decomposed
public :: MULTIPLE, SINGLE_FILE
!> These encoding constants are used to indicate the access mode for a file
public :: APPEND_FILE, OVERWRITE_FILE, READONLY_FILE, WRITEONLY_FILE
!> These encoding constants are used to indicate the discretization position of a variable
public :: CENTER, CORNER, NORTH_FACE, EAST_FACE

! The following are not used in MOM6, but may be used by externals (e.g. SIS2).
public :: create_file
public :: reopen_file
public :: file_type
public :: open_file
public :: get_file_info
public :: get_file_fields
public :: get_file_times

!> Read a field from file using the infrastructure I/O.
interface MOM_read_data
  module procedure MOM_read_data_0d
  module procedure MOM_read_data_0d_int
  module procedure MOM_read_data_1d
  module procedure MOM_read_data_1d_int
  module procedure MOM_read_data_2d
  module procedure MOM_read_data_2d_region
  module procedure MOM_read_data_3d
  module procedure MOM_read_data_3d_region
  module procedure MOM_read_data_4d
end interface MOM_read_data

!> Read a vector from file using the infrastructure I/O.
interface MOM_read_vector
  module procedure MOM_read_vector_2d
  module procedure MOM_read_vector_3d
end interface MOM_read_vector

!> Read a field using native netCDF I/O
!!
!! This function is primarily used for unstructured data which may contain
!! content that cannot be parsed by infrastructure I/O.
interface read_netCDF_data
  ! NOTE: Only 2D I/O is currently used; this should be expanded as needed.
  module procedure read_netCDF_data_2d
end interface read_netCDF_data

!> Write a registered field to an output file, potentially with rotation
interface MOM_write_field
  module procedure MOM_write_field_legacy_4d
  module procedure MOM_write_field_legacy_3d
  module procedure MOM_write_field_legacy_2d
  module procedure MOM_write_field_legacy_1d
  module procedure MOM_write_field_legacy_0d
  module procedure MOM_write_field_4d
  module procedure MOM_write_field_3d
  module procedure MOM_write_field_2d
  module procedure MOM_write_field_1d
  module procedure MOM_write_field_0d
end interface MOM_write_field

!> Read an entire named variable from a named netCDF file using netCDF calls directly, rather
!! than any infrastructure routines and broadcast it from the root PE to the other PEs.
interface read_variable
  module procedure read_variable_0d, read_variable_0d_int
  module procedure read_variable_1d, read_variable_1d_int
  module procedure read_variable_2d, read_variable_3d
end interface read_variable

!> Read a global or variable attribute from a named netCDF file using netCDF calls
!! directly, in some cases reading from the root PE before broadcasting to the other PEs.
interface read_attribute
  module procedure read_attribute_str, read_attribute_real
  module procedure read_attribute_int32, read_attribute_int64
end interface read_attribute

!> Type that stores information that can be used to create a non-decomposed axis.
type :: axis_info
  character(len=32)  :: name = ""       !< The name of this axis for use in files
  character(len=256) :: longname = ""   !< A longer name describing this axis
  character(len=48)  :: units = ""      !< The units of the axis labels
  character(len=8)   :: cartesian = "N" !< A variable indicating which direction
                                        !! this axis corresponds with. Valid values
                                        !! include 'X', 'Y', 'Z', 'T', and 'N' for none.
  integer            :: sense = 0       !< This is 1 for axes whose values increase upward, or -1
                                        !! if they increase downward.  The default, 0, is ignored.
  integer            :: ax_size = 0     !< The number of elements in this axis
  real, allocatable, dimension(:) :: ax_data !< The values of the data on the axis [arbitrary]
end type axis_info

!> Type for describing a 3-d variable for output
type, public :: vardesc
  character(len=64)  :: name               !< Variable name in a NetCDF file
  character(len=48)  :: units              !< Physical dimensions of the variable
  character(len=240) :: longname           !< Long name of the variable
  character(len=8)   :: hor_grid           !< Horizontal grid:  u, v, h, q, Cu, Cv, T, Bu, or 1
  character(len=8)   :: z_grid             !< Vertical grid:  L, i, or 1
  character(len=8)   :: t_grid             !< Time description: s, p, or 1
  character(len=64)  :: cmor_field_name    !< CMOR name
  character(len=64)  :: cmor_units         !< CMOR physical dimensions of the variable
  character(len=240) :: cmor_longname      !< CMOR long name of the variable
  real               :: conversion         !< for unit conversions, such as needed to convert
                                           !! from intensive to extensive [various] or [a A-1 ~> 1]
                                           !! to undo internal dimensional rescaling
  character(len=32)  :: dim_names(5)       !< The names in the file of the axes for this variable
  integer            :: position = -1      !< An integer encoding the horizontal position, it may
                                           !! CENTER, CORNER, EAST_FACE, NORTH_FACE, or 0.
  type(axis_info) :: extra_axes(5)         !< dimensions other than space-time
end type vardesc

!> Type that stores for a global file attribute
type :: attribute_info ; private
  character(len=:), allocatable :: name    !< The name of this attribute
  character(len=:), allocatable :: att_val !< The values of this attribute
end type attribute_info

integer, public :: stdout = stdout_iso  !< standard output unit
integer, public :: stderr = stderr_iso  !< standard output unit

! A note on unit descriptions in comments: MOM6 uses units that can be rescaled for dimensional
! consistency testing. These are noted in comments with units like Z, H, L, and T, along with
! their mks counterparts with notation like "a velocity [Z T-1 ~> m s-1]".  If the units
! vary with the Boussinesq approximation, the Boussinesq variant is given first.
! The functions in this module work with variables with arbitrary units, in which case the
! arbitrary rescaled units are indicated with [A ~> a], while the unscaled units are just [a].


  interface
module subroutine create_file(IO_handle, filename, vars, novars, fields, threading, &
    timeunit, G, dG, GV, checksums, extra_axes, global_atts)
  type(file_type), intent(inout) :: IO_handle
    !< Handle for a files or fileset that is to be opened or reopened for
    !! writing
  character(len=*), intent(in) :: filename
    !< full path to the file to create
  type(vardesc), intent(in) :: vars(:)
    !< structures describing fields written to filename
  integer, intent(in) :: novars
    !< number of fields written to filename
  type(fieldtype), intent(inout) :: fields(:)
    !< array of fieldtypes for each variable
  integer, optional, intent(in) :: threading
    !< SINGLE_FILE or MULTIPLE
  real, optional, intent(in) :: timeunit
    !< length of the units for time [s]. The default value is 86400.0, for 1
    !! day.
  type(ocean_grid_type), optional, intent(in) :: G
    !< ocean horizontal grid structure; G or dG is required if the new file
    !! uses any horizontal grid axes.
  type(dyn_horgrid_type), optional, intent(in) :: dG
    !< dynamic horizontal grid structure; G or dG is required if the new file
    !! uses any horizontal grid axes.
  type(verticalGrid_type), optional, intent(in) :: GV
    !< ocean vertical grid structure, which is ! required if the new file uses
    !! any vertical grid axes.
  integer(kind=int64), optional, intent(in) :: checksums(:,:)
    !< checksums of vars
  type(axis_info), optional, intent(in) :: extra_axes(:)
    !< Types with information about some axes that might be used in this file
  type(attribute_info), optional, intent(in) :: global_atts(:)
    !< Global attributes to write to this file


end subroutine create_file
module subroutine create_MOM_file(IO_handle, filename, vars, novars, fields, &
    threading, timeunit, G, dG, GV, checksums, extra_axes, global_atts)
  class(MOM_file),       intent(inout) :: IO_handle  !< Handle for a files or fileset that is to be
                                                     !! opened or reopened for writing
  character(len=*),      intent(in)    :: filename   !< full path to the file to create
  type(vardesc),         intent(in)    :: vars(:)    !< structures describing fields written to filename
  integer,               intent(in)    :: novars     !< number of fields written to filename
  type(MOM_field),       intent(inout) :: fields(:)  !< array of fieldtypes for each variable
  integer, optional,     intent(in)    :: threading  !< SINGLE_FILE or MULTIPLE
  real, optional,        intent(in)    :: timeunit   !< length of the units for time [s]. The
                                                     !! default value is 86400.0, for 1 day.
  type(ocean_grid_type),   optional, intent(in) :: G !< ocean horizontal grid structure; G or dG
                                                     !! is required if the new file uses any
                                                     !! horizontal grid axes.
  type(dyn_horgrid_type),  optional, intent(in) :: dG !< dynamic horizontal grid structure; G or dG
                                                     !! is required if the new file uses any
                                                     !! horizontal grid axes.
  type(verticalGrid_type), optional, intent(in) :: GV !< ocean vertical grid structure, which is
                                                     !! required if the new file uses any
                                                     !! vertical grid axes.
  integer(kind=int64),     optional, intent(in) :: checksums(:,:) !< checksums of vars
  type(axis_info),         dimension(:), &
                           optional, intent(in) :: extra_axes !< Types with information about
                                                     !! some axes that might be used in this file
  type(attribute_info),    optional, intent(in) :: global_atts(:) !< Global attributes to
                                                     !! write to this file


end subroutine create_MOM_file
module subroutine reopen_file(IO_handle, filename, vars, novars, fields, threading, &
    timeunit, G, dG, GV, extra_axes, global_atts)
  type(file_type), intent(inout) :: IO_handle
    !< Handle for a file or fileset that is to be opened or reopened for
    !! writing
  character(len=*), intent(in) :: filename
    !< full path to the file to create
  type(vardesc), intent(in) :: vars(:)
    !< structures describing fields written to filename
  integer, intent(in) :: novars
    !< number of fields written to filename
  type(fieldtype), intent(inout) :: fields(:)
    !< array of fieldtypes for each variable
  integer, optional, intent(in) :: threading
    !< SINGLE_FILE or MULTIPLE
  real, optional, intent(in) :: timeunit
    !< length of the units for time [s]. The default value is 86400.0, for 1
    !! day.
  type(ocean_grid_type), optional, intent(in) :: G
    !< ocean horizontal grid structure; G or dG is required if a new file uses
    !! any horizontal grid axes.
  type(dyn_horgrid_type), optional, intent(in) :: dG
    !< dynamic horizontal grid structure; G or dG is required if a new file
    !! uses any horizontal grid axes.
  type(verticalGrid_type), optional, intent(in) :: GV
    !< ocean vertical grid structure, which is required if a new file uses any
    !! vertical grid axes.
  type(axis_info), optional, intent(in) :: extra_axes(:)
    !< Types with information about some axes that might be used in this file
  type(attribute_info), optional, intent(in) :: global_atts(:)
    !< Global attributes to write to this file

    !< Wrapper to MOM file
    !< Wrapper to MOM fields

end subroutine reopen_file
module subroutine reopen_MOM_file(IO_handle, filename, vars, novars, fields, &
    threading, timeunit, G, dG, GV, extra_axes, global_atts)
  class(MOM_file),       intent(inout) :: IO_handle  !< Handle for a file or fileset that is to be
                                                     !! opened or reopened for writing
  character(len=*),      intent(in)    :: filename   !< full path to the file to create
  type(vardesc),         intent(in)    :: vars(:)    !< structures describing fields written to filename
  integer,               intent(in)    :: novars     !< number of fields written to filename
  type(MOM_field),       intent(inout) :: fields(:)  !< array of fieldtypes for each variable
  integer, optional,     intent(in)    :: threading  !< SINGLE_FILE or MULTIPLE
  real, optional,        intent(in)    :: timeunit   !< length of the units for time [s]. The
                                                     !! default value is 86400.0, for 1 day.
  type(ocean_grid_type),   optional, intent(in) :: G !< ocean horizontal grid structure; G or dG
                                                     !! is required if a new file uses any
                                                     !! horizontal grid axes.
  type(dyn_horgrid_type),  optional, intent(in) :: dG !< dynamic horizontal grid structure; G or dG
                                                     !! is required if a new file uses any
                                                     !! horizontal grid axes.
  type(verticalGrid_type), optional, intent(in) :: GV !< ocean vertical grid structure, which is
                                                     !! required if a new file uses any
                                                     !! vertical grid axes.
  type(axis_info),         optional, intent(in) :: extra_axes(:)  !< Types with information about
                                                     !! some axes that might be used in this file
  type(attribute_info),    optional, intent(in) :: global_atts(:) !< Global attributes to
                                                     !! write to this file


end subroutine reopen_MOM_file
integer module function stdout_if_root()
end function stdout_if_root
module function num_timelevels(filename, varname, min_dims) result(n_time)
  character(len=*),  intent(in) :: filename   !< name of the file to read
  character(len=*),  intent(in) :: varname    !< variable whose number of time levels
                                              !! are to be returned
  integer, optional, intent(in) :: min_dims   !< The minimum number of dimensions a variable must have
                                              !! if it has a time dimension.  If the variable has 1 less
                                              !! dimension than this, then 0 is returned.
  integer :: n_time                           !< number of time levels varname has in filename


end function num_timelevels
module subroutine get_var_sizes(filename, varname, ndims, sizes, match_case, caller, all_read, dim_names, ncid_in)
  character(len=*),      intent(in)  :: filename   !< Name of the file to read, used here in messages
  character(len=*),      intent(in)  :: varname    !< The variable name, used here for messages
  integer,               intent(out) :: ndims      !< The number of dimensions to the variable
  integer, dimension(:), intent(out) :: sizes      !< The dimension sizes, or 0 for extra values
  logical,     optional, intent(in)  :: match_case !< If false, allow for variables name matches to be
                                                   !! case insensitive, but take a perfect match if
                                                   !! found.  The default is true.
  character(len=*), optional, intent(in) :: caller !< The name of a calling routine for use in error messages
  logical,     optional, intent(in)  :: all_read   !< If present and true, all PEs that call this
                                                   !! routine actually do the read, otherwise only
                                                   !! root PE reads and then it broadcasts the results.
  character(len=*), dimension(:), &
               optional, intent(out) :: dim_names  !< The names of the dimensions for this variable
  integer,     optional, intent(in)  :: ncid_in    !< The netCDF ID of an open file.  If absent, the
                                                   !! file is opened and closed within this routine.


end subroutine get_var_sizes
module subroutine read_var_sizes(filename, varname, ndims, sizes, match_case, caller, dim_names, ncid_in)
  character(len=*),      intent(in)  :: filename   !< Name of the file to read, used here in messages
  character(len=*),      intent(in)  :: varname    !< The variable name, used here for messages
  integer,               intent(out) :: ndims      !< The number of dimensions to the variable
  integer, dimension(:), intent(out) :: sizes      !< The dimension sizes, or 0 for extra values
  logical,     optional, intent(in)  :: match_case !< If false, allow for variables name matches to be
                                                   !! case insensitive, but take a perfect match if
                                                   !! found.  The default is true.
  character(len=*), &
               optional, intent(in)  :: caller     !< The name of a calling routine for use in error messages
  character(len=*), dimension(:), &
               optional, intent(out) :: dim_names  !< The names of the dimensions for this variable
  integer,     optional, intent(in)  :: ncid_in    !< The netCDF ID of an open file.  If absent, the
                                                   !! file is opened and closed within this routine.

end subroutine read_var_sizes
module subroutine read_variable_0d(filename, varname, var, ncid_in, scale)
  character(len=*),  intent(in)    :: filename !< The name of the file to read
  character(len=*),  intent(in)    :: varname  !< The variable name of the data in the file
  real,              intent(inout) :: var      !< The scalar into which to read the data in arbitrary units [A ~> a]
  integer, optional, intent(in)    :: ncid_in  !< The netCDF ID of an open file.  If absent, the
                                               !! file is opened and closed within this routine
  real,    optional, intent(in)    :: scale    !< A scaling factor that the variable is multiplied by
                                               !! before it is returned to convert from the units in the file
                                               !! to the internal units for this variable [A a-1 ~> 1]

end subroutine read_variable_0d
module subroutine read_variable_1d(filename, varname, var, ncid_in, scale)
  character(len=*),   intent(in)    :: filename !< The name of the file to read
  character(len=*),   intent(in)    :: varname  !< The variable name of the data in the file
  real, dimension(:), intent(inout) :: var      !< The 1-d array into which to read the data in arbitrary units [A ~> a]
  integer,  optional, intent(in)    :: ncid_in  !< The netCDF ID of an open file.  If absent, the
                                                !! file is opened and closed within this routine
  real,     optional, intent(in)    :: scale    !< A scaling factor that the variable is multiplied by
                                                !! before it is returned to convert from the units in the file
                                                !! to the internal units for this variable [A a-1 ~> 1]

end subroutine read_variable_1d
module subroutine read_variable_0d_int(filename, varname, var, ncid_in)
  character(len=*),  intent(in)    :: filename !< The name of the file to read
  character(len=*),  intent(in)    :: varname  !< The variable name of the data in the file
  integer,           intent(inout) :: var      !< The scalar into which to read the data
  integer, optional, intent(in)    :: ncid_in  !< The netCDF ID of an open file.  If absent, the
                                               !! file is opened and closed within this routine.

end subroutine read_variable_0d_int
module subroutine read_variable_1d_int(filename, varname, var, ncid_in)
  character(len=*),      intent(in)    :: filename !< The name of the file to read
  character(len=*),      intent(in)    :: varname  !< The variable name of the data in the file
  integer, dimension(:), intent(inout) :: var      !< The 1-d array into which to read the data
  integer, optional,     intent(in)    :: ncid_in  !< The netCDF ID of an open file.  If absent, the
                                                   !! file is opened and closed within this routine.

end subroutine read_variable_1d_int
module subroutine read_variable_2d(filename, varname, var, start, nread, ncid_in)
  character(len=*), intent(in) :: filename  !< Name of file to be read
  character(len=*), intent(in) :: varname   !< Name of variable to be read
  real, intent(out)            :: var(:,:)  !< Output array of variable [arbitrary]
  integer, optional, intent(in) :: start(:) !< Starting index on each axis.
  integer, optional, intent(in) :: nread(:) !< Number of values to be read along each axis
  integer, optional, intent(in) :: ncid_in  !< netCDF ID of an opened file.
              !! If absent, the file is opened and closed within this routine.


  ! Validate shape of start and nread
end subroutine read_variable_2d
module subroutine read_variable_3d(filename, varname, var, start, nread, ncid_in)
  character(len=*), intent(in) :: filename  !< Name of file to be read
  character(len=*), intent(in) :: varname   !< Name of variable to be read
  real, intent(out)            :: var(:,:,:)  !< Output array of variable [arbitrary]
  integer, optional, intent(in) :: start(:) !< Starting index on each axis.
  integer, optional, intent(in) :: nread(:) !< Number of values to be read along each axis
  integer, optional, intent(in) :: ncid_in  !< netCDF ID of an opened file.
              !! If absent, the file is opened and closed within this routine.


  ! Validate shape of start and nread
end subroutine read_variable_3d
module subroutine read_attribute_str(filename, attname, att_val, varname, found, all_read, ncid_in)
  character(len=*),           intent(in)  :: filename !< Name of the file to read
  character(len=*),           intent(in)  :: attname  !< Name of the attribute to read
  character(:), allocatable,  intent(out) :: att_val  !< The value of the attribute
  character(len=*), optional, intent(in)  :: varname  !< The name of the variable whose attribute will
                                                      !! be read. If missing, read a global attribute.
  logical,          optional, intent(out) :: found    !< Returns true if the attribute is found
  logical,          optional, intent(in)  :: all_read !< If present and true, all PEs that call this
                                                      !! routine actually do the read, otherwise only
                                                      !! root PE reads and then broadcasts the results.
  integer,          optional, intent(in)  :: ncid_in  !< The netCDF ID of an open file.  If absent, the
                                                      !! file is opened and closed within this routine.

end subroutine read_attribute_str
module subroutine read_attribute_int32(filename, attname, att_val, varname, found, all_read, ncid_in)
  character(len=*),           intent(in)  :: filename !< Name of the file to read
  character(len=*),           intent(in)  :: attname  !< Name of the attribute to read
  integer(kind=int32),        intent(out) :: att_val  !< The value of the attribute
  character(len=*), optional, intent(in)  :: varname  !< The name of the variable whose attribute will
                                                      !! be read. If missing, read a global attribute.
  logical,          optional, intent(out) :: found    !< Returns true if the attribute is found
  logical,          optional, intent(in)  :: all_read !< If present and true, all PEs that call this
                                                      !! routine actually do the read, otherwise only
                                                      !! root PE reads and then broadcasts the results.
  integer,        optional, intent(in)    :: ncid_in  !< The netCDF ID of an open file.  If absent, the
                                                      !! file is opened and closed within this routine.

end subroutine read_attribute_int32
module subroutine read_attribute_int64(filename, attname, att_val, varname, found, all_read, ncid_in)
  character(len=*),           intent(in)  :: filename !< Name of the file to read
  character(len=*),           intent(in)  :: attname  !< Name of the attribute to read
  integer(kind=int64),        intent(out) :: att_val  !< The value of the attribute
  character(len=*), optional, intent(in)  :: varname  !< The name of the variable whose attribute will
                                                      !! be read. If missing, read a global attribute.
  logical,          optional, intent(out) :: found    !< Returns true if the attribute is found
  logical,          optional, intent(in)  :: all_read !< If present and true, all PEs that call this
                                                      !! routine actually do the read, otherwise only
                                                      !! root PE reads and then broadcasts the results.
  integer,        optional, intent(in)    :: ncid_in  !< The netCDF ID of an open file.  If absent, the
                                                      !! file is opened and closed within this routine.

end subroutine read_attribute_int64
module subroutine read_attribute_real(filename, attname, att_val, varname, found, all_read, ncid_in)
  character(len=*),           intent(in)  :: filename !< Name of the file to read
  character(len=*),           intent(in)  :: attname  !< Name of the attribute to read
  real,                       intent(out) :: att_val  !< The value of the attribute [arbitrary]
  character(len=*), optional, intent(in)  :: varname  !< The name of the variable whose attribute will
                                                      !! be read. If missing, read a global attribute.
  logical,          optional, intent(out) :: found    !< Returns true if the attribute is found
  logical,          optional, intent(in)  :: all_read !< If present and true, all PEs that call this
                                                      !! routine actually do the read, otherwise only
                                                      !! root PE reads and then broadcasts the results.
  integer,          optional, intent(in)  :: ncid_in  !< The netCDF ID of an open file.  If absent, the
                                                      !! file is opened and closed within this routine.

end subroutine read_attribute_real
module subroutine open_file_to_read(filename, ncid, success)
  character(len=*),  intent(in)  :: filename   !< path and name of the file to open for reading
  integer,           intent(out) :: ncid       !< The netcdf handle for the file
  logical, optional, intent(out) :: success    !< Returns true if the file was opened, or if this
                                               !! argument is not present, failure is fatal error.
  ! Local variables

end subroutine open_file_to_read
module subroutine close_file_to_read(ncid, filename)
  integer,                    intent(inout) :: ncid       !< The netcdf handle for the file to close
  character(len=*), optional, intent(in)    :: filename   !< path and name of the file to close
end subroutine close_file_to_read
module subroutine get_varid(varname, ncid, filename, varid, match_case, found)
  character(len=*),  intent(in)  :: varname    !< The name of the variable that is being sought
  integer,           intent(in)  :: ncid       !< The open netcdf handle for the file
  character(len=*),  intent(in)  :: filename   !< name of the file to read, used here in messages
  integer,           intent(out) :: varid      !< The netcdf handle for the variable
  logical, optional, intent(in)  :: match_case !< If false, allow for variables name matches to be
                                               !! case insensitive, but take a perfect match if
                                               !! found.  The default is true.
  logical, optional, intent(out) :: found      !< Returns true if the attribute is found


end subroutine get_varid
module subroutine verify_variable_units(filename, varname, expected_units, msg, ierr, alt_units)
  character(len=*),           intent(in)    :: filename  !< File name
  character(len=*),           intent(in)    :: varname   !< Variable name
  character(len=*),           intent(in)    :: expected_units !< Expected units of variable
  character(len=*),           intent(inout) :: msg       !< Message to use for errors
  logical,                    intent(out)   :: ierr      !< True if an error occurs
  character(len=*), optional, intent(in)    :: alt_units !< Alterate acceptable units of variable

  ! Local variables

end subroutine verify_variable_units
module function var_desc(name, units, longname, hor_grid, z_grid, t_grid, cmor_field_name, &
                  cmor_units, cmor_longname, conversion, caller, position, dim_names, &
                  extra_axes, fixed) result(vd)
  character(len=*),           intent(in) :: name            !< variable name
  character(len=*), optional, intent(in) :: units           !< variable units
  character(len=*), optional, intent(in) :: longname        !< variable long name
  character(len=*), optional, intent(in) :: hor_grid        !< A character string indicating the horizontal
                                                            !! position of this variable
  character(len=*), optional, intent(in) :: z_grid          !< variable vertical staggering
  character(len=*), optional, intent(in) :: t_grid          !< time description: s, p, or 1
  character(len=*), optional, intent(in) :: cmor_field_name !< CMOR name
  character(len=*), optional, intent(in) :: cmor_units      !< CMOR physical dimensions of variable
  character(len=*), optional, intent(in) :: cmor_longname   !< CMOR long name
  real            , optional, intent(in) :: conversion      !< for unit conversions, such as needed to
                                                            !! convert from intensive to extensive
                                                            !! [various] or [a A-1 ~> 1]
  character(len=*), optional, intent(in) :: caller          !< The calling routine for error messages
  integer,          optional, intent(in) :: position        !< A coded integer indicating the horizontal position
                                                            !! of this variable if it has such dimensions.
                                                            !! Valid values include CORNER, CENTER, EAST_FACE
                                                            !! NORTH_FACE, and 0 for no horizontal dimensions.
  character(len=*), dimension(:), &
                    optional, intent(in) :: dim_names       !< The names of the dimensions of this variable
  type(axis_info),  dimension(:), &
                    optional, intent(in) :: extra_axes      !< dimensions other than space-time
  logical,          optional, intent(in) :: fixed           !< If true, this does not evolve with time
  type(vardesc)                          :: vd              !< vardesc type that is created

end function var_desc
module subroutine modify_vardesc(vd, name, units, longname, hor_grid, z_grid, t_grid, &
                 cmor_field_name, cmor_units, cmor_longname, conversion, caller, position, dim_names, &
                 extra_axes)
  type(vardesc),              intent(inout) :: vd              !< vardesc type that is modified
  character(len=*), optional, intent(in)    :: name            !< name of variable
  character(len=*), optional, intent(in)    :: units           !< units of variable
  character(len=*), optional, intent(in)    :: longname        !< long name of variable
  character(len=*), optional, intent(in)    :: hor_grid        !< horizontal staggering of variable
  character(len=*), optional, intent(in)    :: z_grid          !< vertical staggering of variable
  character(len=*), optional, intent(in)    :: t_grid          !< time description: s, p, or 1
  character(len=*), optional, intent(in)    :: cmor_field_name !< CMOR name
  character(len=*), optional, intent(in)    :: cmor_units      !< CMOR physical dimensions of variable
  character(len=*), optional, intent(in)    :: cmor_longname   !< CMOR long name
  real            , optional, intent(in)    :: conversion      !< A multiplicative factor for unit conversions,
                                                               !! such as needed to convert from intensive to
                                                               !! extensive or dimensional consistency testing
                                                               !! [various] or [a A-1 ~> 1]
  character(len=*), optional, intent(in)    :: caller          !< The calling routine for error messages
  integer,          optional, intent(in)    :: position        !< A coded integer indicating the horizontal position
                                                               !! of this variable if it has such dimensions.
                                                               !! Valid values include CORNER, CENTER, EAST_FACE
                                                               !! NORTH_FACE, and 0 for no horizontal dimensions.
  character(len=*), dimension(:), &
                    optional, intent(in)    :: dim_names       !< The names of the dimensions of this variable
  type(axis_info),  dimension(:), &
                    optional, intent(in)    :: extra_axes      !< dimensions other than space-time


end subroutine modify_vardesc
integer module function position_from_horgrid(hor_grid)
  character(len=*), intent(in)    :: hor_grid        !< horizontal staggering of variable

end function position_from_horgrid
module subroutine set_axis_info(axis, name, units, longname, ax_size, ax_data, cartesian, sense)
  type(axis_info),              intent(inout) :: axis  !< A type with information about a named axis
  character(len=*),             intent(in)    :: name  !< The name of this axis for use in files
  character(len=*),   optional, intent(in)    :: units !< The units of the axis labels
  character(len=*),   optional, intent(in)    :: longname  !< Long name of the axis variable
  integer,            optional, intent(in)    :: ax_size !< The number of elements in this axis
  real, dimension(:), optional, intent(in)    :: ax_data !< The values of the data on the axis [arbitrary]
  character(len=*),   optional, intent(in)    :: cartesian !< A variable indicating which direction this axis
                                                       !! axis corresponds with. Valid values
                                                       !! include 'X', 'Y', 'Z', 'T', and 'N' (the default) for none.
  integer,            optional, intent(in)    :: sense !< This is 1 for axes whose values increase upward, or -1
                                                       !! if they increase downward.  The default, 0, is ignored.

end subroutine set_axis_info
module subroutine delete_axis_info(axes)
  type(axis_info), dimension(:), intent(inout) :: axes  !< An array with information about named axes

end subroutine delete_axis_info
module subroutine get_axis_info(axis,name,longname,units,cartesian,ax_size,ax_data)
  type(axis_info), intent(in) :: axis                               !< An axis type
  character(len=*), intent(out), optional    :: name                !< The axis name.
  character(len=*), intent(out), optional    :: longname            !< The axis longname.
  character(len=*), intent(out), optional    :: units               !< The axis units.
  character(len=*), intent(out), optional    :: cartesian           !< The cartesian attribute
                                                                    !! of the axis [X,Y,Z,T].
  integer,          intent(out), optional   :: ax_size              !< The size of the axis.
  real, optional, allocatable, dimension(:), intent(out) :: ax_data !< The axis label data [arbitrary]

end subroutine get_axis_info
module subroutine set_attribute_info(attribute, name, str_value)
  type(attribute_info), intent(inout) :: attribute !< A type with information about a named attribute
  character(len=*),     intent(in)    :: name      !< The name of this attribute for use in files
  character(len=*),     intent(in)    :: str_value !< The value of this attribute

end subroutine set_attribute_info
module subroutine delete_attribute_info(atts)
  type(attribute_info), dimension(:), intent(inout) :: atts  !< An array of global attributes

end subroutine delete_attribute_info
module function cmor_long_std(longname) result(std_name)
  character(len=*), intent(in) :: longname  !< The CMOR longname being converted
  character(len=len(longname)) :: std_name  !< The CMOR standard name generated from longname


end function cmor_long_std
module subroutine query_vardesc(vd, name, units, longname, hor_grid, z_grid, t_grid, &
                         cmor_field_name, cmor_units, cmor_longname, conversion, caller, &
                         extra_axes, position, dim_names)
  type(vardesc),              intent(in)  :: vd                 !< vardesc type that is queried
  character(len=*), optional, intent(out) :: name               !< name of variable
  character(len=*), optional, intent(out) :: units              !< units of variable
  character(len=*), optional, intent(out) :: longname           !< long name of variable
  character(len=*), optional, intent(out) :: hor_grid           !< horizontal staggering of variable
  character(len=*), optional, intent(out) :: z_grid             !< verticle staggering of variable
  character(len=*), optional, intent(out) :: t_grid             !< time description: s, p, or 1
  character(len=*), optional, intent(out) :: cmor_field_name    !< CMOR name
  character(len=*), optional, intent(out) :: cmor_units         !< CMOR physical dimensions of variable
  character(len=*), optional, intent(out) :: cmor_longname      !< CMOR long name
  real            , optional, intent(out) :: conversion         !< for unit conversions, such as needed to
                                                                !! convert from intensive to extensive
                                                                !! [various] or [a A-1 ~> 1]
  character(len=*), optional, intent(in)  :: caller             !< calling routine?
  type(axis_info),  dimension(5), &
                    optional, intent(out) :: extra_axes      !< dimensions other than space-time
  integer,          optional, intent(out) :: position        !< A coded integer indicating the horizontal position
                                                            !! of this variable if it has such dimensions.
                                                            !! Valid values include CORNER, CENTER, EAST_FACE
                                                            !! NORTH_FACE, and 0 for no horizontal dimensions.
  character(len=*), dimension(:), &
                    optional, intent(out) :: dim_names       !< The names of the dimensions of this variable

end subroutine query_vardesc
module subroutine MOM_read_data_0d(filename, fieldname, data, timelevel, scale, MOM_Domain, &
                            global_file, file_may_be_4d)
  character(len=*), intent(in)  :: filename     !< Input filename
  character(len=*), intent(in)  :: fieldname    !< Field variable name
  real, intent(inout)           :: data         !< Field value in arbitrary units [A ~> a]
  integer, optional, intent(in) :: timelevel    !< Time level to read in file
  real, optional, intent(in)    :: scale        !< A scaling factor that the variable is multiplied by
                                                !! before it is returned to convert from the units in the file
                                                !! to the internal units for this variable [A a-1 ~> 1]
  type(MOM_domain_type), optional, intent(in) :: MOM_Domain !< Model domain decomposition
  logical, optional, intent(in) :: global_file    !< If true, read from a single file
  logical, optional, intent(in) :: file_may_be_4d !< If true, fields may be stored
                                                  !! as 4d arrays in the file.

end subroutine MOM_read_data_0d
module subroutine MOM_read_data_0d_int(filename, fieldname, data, timelevel)
  character(len=*), intent(in) :: filename    !< Input filename
  character(len=*), intent(in) :: fieldname   !< Field variable name
  integer, intent(inout) :: data              !< Field value
  integer, optional, intent(in) :: timelevel  !< Time level to read in file

end subroutine MOM_read_data_0d_int
module subroutine MOM_read_data_1d(filename, fieldname, data, timelevel, scale, MOM_Domain, &
                            global_file, file_may_be_4d)
  character(len=*), intent(in)  :: filename   !< Input filename
  character(len=*), intent(in)  :: fieldname  !< Field variable name
  real, dimension(:), intent(inout) :: data   !< Field value in arbitrary units [A ~> a]
  integer, optional, intent(in) :: timelevel  !< Time level to read in file
  real, optional, intent(in)    :: scale      !< A scaling factor that the variable is multiplied by
                                              !! before it is returned to convert from the units in the file
                                              !! to the internal units for this variable [A a-1 ~> 1]
  type(MOM_domain_type), optional, intent(in) :: MOM_Domain !< Model domain decomposition
  logical, optional, intent(in) :: global_file    !< If true, read from a single file
  logical, optional, intent(in) :: file_may_be_4d !< If true, fields may be stored
                                                  !! as 4d arrays in the file.

end subroutine MOM_read_data_1d
module subroutine MOM_read_data_1d_int(filename, fieldname, data, timelevel)
  character(len=*), intent(in) :: filename    !< Input filename
  character(len=*), intent(in) :: fieldname   !< Field variable name
  integer, dimension(:), intent(inout) :: data  !< Field value
  integer, optional, intent(in) :: timelevel  !< Time level to read in file

end subroutine MOM_read_data_1d_int
module subroutine MOM_read_data_2d(filename, fieldname, data, MOM_Domain, timelevel, position, &
                            scale, global_file, file_may_be_4d, turns)
  character(len=*), intent(in)  :: filename  !< Input filename
  character(len=*), intent(in)  :: fieldname !< Field variable name
  real, dimension(:,:), intent(inout) :: data   !< Field value in arbitrary units [A ~> a]
  type(MOM_domain_type), target, &
                     intent(in) :: MOM_Domain !< Model domain decomposition
  integer, optional, intent(in) :: timelevel !< Time level to read in file
  integer, optional, intent(in) :: position  !< Grid positioning flag
  real, optional, intent(in)    :: scale     !< A scaling factor that the variable is multiplied by
                                             !! before it is returned to convert from the units in the file
                                             !! to the internal units for this variable [A a-1 ~> 1]
  logical, optional, intent(in) :: global_file    !< If true, read from a single file
  logical, optional, intent(in) :: file_may_be_4d !< If true, fields may be stored
                                                  !! as 4d arrays in the file.
  integer, optional, intent(in) :: turns        !< Number of quarter-turns to rotate the data.  If absent
                                                !! the number of turns is taken from MOM_Domain.

  ! Local variables

end subroutine MOM_read_data_2d
module subroutine read_netCDF_data_2d(filename, fieldname, values, MOM_Domain, &
                            timelevel, position, rescale, turns)
  character(len=*), intent(in) :: filename
    !< Input filename
  character(len=*), intent(in)  :: fieldname
    !< Field variable name
  real, intent(inout) :: values(:,:)
    !< Field values read from the file.  It would be intent(out) but for the
    !! need to preserve any initialized values in the halo regions.
  type(MOM_domain_type), intent(in) :: MOM_Domain
    !< Model domain decomposition
  integer, optional, intent(in) :: timelevel
    !< Time level to read in file
  integer, optional, intent(in) :: position
    !< Grid positioning flag
  real, optional, intent(in) :: rescale
    !< Rescale factor, omitting this is the same as setting it to 1.
  integer, optional, intent(in) :: turns
    !< Number of quarter-turns to rotate the data.  If absent the number of turns is taken
    !! from MOM_Domain.

    ! Number of quarter-turns from input to model grid
    ! Field array on the unrotated input grid
    ! netCDF file handle

  ! General-purpose IO will require the following arguments, but they are not
  ! yet implemented, so we raise an error if they are present.

  ! Fields are currently assumed on cell centers, and position is unsupported
end subroutine read_netCDF_data_2d
module subroutine MOM_read_data_2d_region(filename, fieldname, data, start, nread, MOM_domain, &
                                   no_domain, scale, turns)
  character(len=*), intent(in)  :: filename   !< Input filename
  character(len=*), intent(in)  :: fieldname  !< Field variable name
  real, dimension(:,:), intent(inout) :: data !< Field value in arbitrary units [A ~> a]
  integer, dimension(:), intent(in) :: start  !< Starting index for each axis.
                                              !! In 2d, start(3:4) must be 1.
  integer, dimension(:), intent(in) :: nread  !< Number of values to read along each axis.
                                              !! In 2d, nread(3:4) must be 1.
  type(MOM_domain_type), optional, intent(in) :: MOM_Domain !< Model domain decomposition
  logical, optional, intent(in) :: no_domain  !< If true, field does not use
                                              !! domain decomposion.
  real, optional, intent(in)    :: scale      !< A scaling factor that the variable is multiplied by
                                              !! before it is returned to convert from the units in the file
                                              !! to the internal units for this variable [A a-1 ~> 1]
  integer, optional, intent(in) :: turns      !< Number of quarter turns from
                                              !! input to model grid


end subroutine MOM_read_data_2d_region
module subroutine MOM_read_data_3d(filename, fieldname, data, MOM_Domain, timelevel, position, &
                            scale, global_file, file_may_be_4d, turns)
  character(len=*), intent(in)  :: filename     !< Input filename
  character(len=*), intent(in)  :: fieldname    !< Field variable name
  real, dimension(:,:,:), intent(inout) :: data !< Field value in arbitrary units [A ~> a]
  type(MOM_domain_type), target, &
                     intent(in) :: MOM_Domain   !< Model domain decomposition
  integer, optional, intent(in) :: timelevel    !< Time level to read in file
  integer, optional, intent(in) :: position     !< Grid positioning flag
  real, optional, intent(in)    :: scale        !< A scaling factor that the variable is multiplied by
                                                !! before it is returned to convert from the units in the file
                                                !! to the internal units for this variable [A a-1 ~> 1]
  logical, optional, intent(in) :: global_file  !< If true, read from a single file
  logical, optional, intent(in) :: file_may_be_4d !< If true, fields may be stored
                                                  !! as 4d arrays in the file.
  integer, optional, intent(in) :: turns        !< Number of quarter-turns to rotate the data.  If absent
                                                !! the number of turns is taken from MOM_Domain.

  ! Local variables

end subroutine MOM_read_data_3d
module subroutine MOM_read_data_3d_region(filename, fieldname, data, start, nread, MOM_domain, &
                                   no_domain, scale, turns)
  character(len=*), intent(in)  :: filename   !< Input filename
  character(len=*), intent(in)  :: fieldname  !< Field variable name
  real, dimension(:,:,:), intent(inout) :: data !< Field value in arbitrary units [A ~> a]
  integer, dimension(:), intent(in) :: start  !< Starting index for each axis.
  integer, dimension(:), intent(in) :: nread  !< Number of values to read along each axis.
  type(MOM_domain_type), optional, intent(in) :: MOM_Domain !< Model domain decomposition
  logical, optional, intent(in) :: no_domain  !< If true, field does not use
                                              !! domain decomposion.
  real, optional, intent(in)    :: scale      !< A scaling factor that the variable is multiplied by
                                              !! before it is returned to convert from the units in the file
                                              !! to the internal units for this variable [A a-1 ~> 1]
  integer, optional, intent(in) :: turns      !< Number of quarter turns from
                                              !! input to model grid


end subroutine MOM_read_data_3d_region
module subroutine MOM_read_data_4d(filename, fieldname, data, MOM_Domain, &
                            timelevel, position, scale, global_file, turns)
  character(len=*), intent(in) :: filename      !< Input filename
  character(len=*), intent(in) :: fieldname     !< Field variable name
  real, dimension(:,:,:,:), intent(inout) :: data !< Field value in arbitrary units [A ~> a]
  type(MOM_domain_type), target, &
                     intent(in) :: MOM_Domain   !< Model domain decomposition
  integer, optional, intent(in) :: timelevel    !< Time level to read in file
  integer, optional, intent(in) :: position     !< Grid positioning flag
  real,    optional, intent(in) :: scale        !< A scaling factor that the variable is multiplied by
                                                !! before it is returned to convert from the units in the file
                                                !! to the internal units for this variable [A a-1 ~> 1]
  logical, optional, intent(in) :: global_file  !< If true, read from a single file
  integer, optional, intent(in) :: turns        !< Number of quarter-turns to rotate the data.  If absent
                                                !! the number of turns is taken from MOM_Domain.

  ! Local variables

end subroutine MOM_read_data_4d
module subroutine MOM_read_vector_2d(filename, u_fieldname, v_fieldname, u_data, v_data, MOM_Domain, &
                              timelevel, stagger, scalar_pair, scale, turns)
  character(len=*), intent(in) :: filename      !< Input filename
  character(len=*), intent(in) :: u_fieldname   !< Field variable name in u
  character(len=*), intent(in) :: v_fieldname   !< Field variable name in v
  real, dimension(:,:), intent(inout) :: u_data !< Field value at u points in arbitrary units [A ~> a]
  real, dimension(:,:), intent(inout) :: v_data !< Field value at v points in arbitrary units  [A ~> a]
  type(MOM_domain_type), target, &
                     intent(in) :: MOM_Domain   !< Model domain decomposition
  integer, optional, intent(in) :: timelevel    !< Time level to read in file
  integer, optional, intent(in) :: stagger      !< Grid staggering flag
  logical, optional, intent(in) :: scalar_pair  !< True if tuple is not a vector
  real,    optional, intent(in) :: scale        !< A scaling factor that the vector is multiplied by
                                                !! before it is returned to convert from the units in the file
                                                !! to the internal units for this variable [A a-1 ~> 1]
  integer, optional, intent(in) :: turns        !< Number of quarter-turns to rotate the data.  If absent
                                                !! the number of turns is taken from MOM_Domain.

  ! Local variables

end subroutine MOM_read_vector_2d
module subroutine MOM_read_vector_3d(filename, u_fieldname, v_fieldname, u_data, v_data, MOM_Domain, &
                              timelevel, stagger, scalar_pair, scale, turns)
  character(len=*), intent(in) :: filename      !< Input filename
  character(len=*), intent(in) :: u_fieldname   !< Field variable name in u
  character(len=*), intent(in) :: v_fieldname   !< Field variable name in v
  real, dimension(:,:,:), intent(inout) :: u_data !< Field value in u in arbitrary units [A ~> a]
  real, dimension(:,:,:), intent(inout) :: v_data !< Field value in v in arbitrary units [A ~> a]
  type(MOM_domain_type), target, &
                     intent(in) :: MOM_Domain   !< Model domain decomposition
  integer, optional, intent(in) :: timelevel    !< Time level to read in file
  integer, optional, intent(in) :: stagger      !< Grid staggering flag
  logical, optional, intent(in) :: scalar_pair  !< True if tuple is not a vector
  real,    optional, intent(in) :: scale        !< A scaling factor that the vector is multiplied by
                                                !! before it is returned to convert from the units in the file
                                                !! to the internal units for this variable [A a-1 ~> 1]
  integer, optional, intent(in) :: turns        !< Number of quarter-turns to rotate the data.  If absent
                                                !! the number of turns is taken from MOM_Domain.

  ! Local variables

end subroutine MOM_read_vector_3d
module subroutine MOM_write_field_legacy_4d(IO_handle, field_md, MOM_domain, field, tstamp, tile_count, &
                              fill_value, turns, scale, unscale, zero_zeros)
  type(file_type),          intent(inout) :: IO_handle  !< Handle for a file that is open for writing
  type(fieldtype),          intent(in)    :: field_md   !< Field type with metadata
  type(MOM_domain_type),    intent(in)    :: MOM_domain !< The MOM_Domain that describes the decomposition
  real, dimension(:,:,:,:), intent(inout) :: field      !< Unrotated field to write in arbitrary units [A ~> a]
  real,           optional, intent(in)    :: tstamp     !< Model timestamp, often in [days]
  integer,        optional, intent(in)    :: tile_count !< PEs per tile (default: 1)
  real,           optional, intent(in)    :: fill_value !< Missing data fill value in the units used in the file [a]
  integer,        optional, intent(in)    :: turns      !< Number of quarter-turns to rotate the data
  real,           optional, intent(in)    :: scale      !< A scaling factor that the field is multiplied by before
                                                        !! it is written [a A-1 ~> 1], for example to convert it
                                                        !! from its internal units to the desired units for output
  real,           optional, intent(in)    :: unscale    !< A scaling factor that the field is multiplied by before
                                                        !! it is written [a A-1 ~> 1], for example to convert it
                                                        !! from its internal units to the desired units for output.
                                                        !! Here scale and unscale are synonymous, but unscale
                                                        !! takes precedence if both are present.
  logical,        optional, intent(in)    :: zero_zeros !< If present and true, convert negative zeros
                                                        !! into ordinary signless zeros.

  ! Local variables
                                           ! rescaled [A ~> a] then [a]

end subroutine MOM_write_field_legacy_4d
module subroutine MOM_write_field_legacy_3d(IO_handle, field_md, MOM_domain, field, tstamp, tile_count, &
                              fill_value, turns, scale, unscale, zero_zeros)
  type(file_type),        intent(inout) :: IO_handle  !< Handle for a file that is open for writing
  type(fieldtype),        intent(in)    :: field_md   !< Field type with metadata
  type(MOM_domain_type),  intent(in)    :: MOM_domain !< The MOM_Domain that describes the decomposition
  real, dimension(:,:,:), intent(inout) :: field      !< Unrotated field to write in arbitrary units [A ~> a]
  real,         optional, intent(in)    :: tstamp     !< Model timestamp, often in [days]
  integer,      optional, intent(in)    :: tile_count !< PEs per tile (default: 1)
  real,         optional, intent(in)    :: fill_value !< Missing data fill value in the units used in the file [a]
  integer,      optional, intent(in)    :: turns      !< Number of quarter-turns to rotate the data
  real,         optional, intent(in)    :: scale      !< A scaling factor that the field is multiplied by before
                                                      !! it is written [a A-1 ~> 1], for example to convert it
                                                      !! from its internal units to the desired units for output
  real,         optional, intent(in)    :: unscale    !< A scaling factor that the field is multiplied by before
                                                      !! it is written [a A-1 ~> 1], for example to convert it
                                                      !! from its internal units to the desired units for output.
                                                      !! Here scale and unscale are synonymous, but unscale
                                                      !! takes precedence if both are present.
  logical,      optional, intent(in)    :: zero_zeros !< If present and true, convert negative zeros
                                                      !! into ordinary signless zeros.

  ! Local variables
                                         ! rescaled [A ~> a] then [a]

end subroutine MOM_write_field_legacy_3d
module subroutine MOM_write_field_legacy_2d(IO_handle, field_md, MOM_domain, field, tstamp, tile_count, &
                              fill_value, turns, scale, unscale, zero_zeros)
  type(file_type),        intent(inout) :: IO_handle  !< Handle for a file that is open for writing
  type(fieldtype),        intent(in)    :: field_md   !< Field type with metadata
  type(MOM_domain_type),  intent(in)    :: MOM_domain !< The MOM_Domain that describes the decomposition
  real, dimension(:,:),   intent(inout) :: field      !< Unrotated field to write in arbitrary units [A ~> a]
  real,         optional, intent(in)    :: tstamp     !< Model timestamp, often in [days]
  integer,      optional, intent(in)    :: tile_count !< PEs per tile (default: 1)
  real,         optional, intent(in)    :: fill_value !< Missing data fill value
  integer,      optional, intent(in)    :: turns      !< Number of quarter-turns to rotate the data
  real,         optional, intent(in)    :: scale      !< A scaling factor that the field is multiplied by before
                                                      !! it is written [a A-1 ~> 1], for example to convert it
                                                      !! from its internal units to the desired units for output
  real,         optional, intent(in)    :: unscale    !< A scaling factor that the field is multiplied by before
                                                      !! it is written [a A-1 ~> 1], for example to convert it
                                                      !! from its internal units to the desired units for output.
                                                      !! Here scale and unscale are synonymous, but unscale
                                                      !! takes precedence if both are present.
  logical,      optional, intent(in)    :: zero_zeros !< If present and true, convert negative zeros
                                                      !! into ordinary signless zeros.

  ! Local variables
                                       ! rescaled [A ~> a] then [a]

end subroutine MOM_write_field_legacy_2d
module subroutine MOM_write_field_legacy_1d(IO_handle, field_md, field, tstamp, fill_value, scale, unscale, zero_zeros)
  type(file_type),        intent(inout) :: IO_handle  !< Handle for a file that is open for writing
  type(fieldtype),        intent(in)    :: field_md   !< Field type with metadata
  real, dimension(:),     intent(in)    :: field      !< Field to write in arbitrary units [A ~> a]
  real,         optional, intent(in)    :: tstamp     !< Model timestamp, often in [days]
  real,         optional, intent(in)    :: fill_value !< Missing data fill value [a]
  real,         optional, intent(in)    :: scale      !< A scaling factor that the field is multiplied by before
                                                      !! it is written [a A-1 ~> 1], for example to convert it
                                                      !! from its internal units to the desired units for output
  real,         optional, intent(in)    :: unscale    !< A scaling factor that the field is multiplied by before
                                                      !! it is written [a A-1 ~> 1], for example to convert it
                                                      !! from its internal units to the desired units for output.
                                                      !! Here scale and unscale are synonymous, but unscale
                                                      !! takes precedence if both are present.
  logical,      optional, intent(in)    :: zero_zeros !< If present and true, convert negative zeros
                                                      !! into ordinary signless zeros.

  ! Local variables

end subroutine MOM_write_field_legacy_1d
module subroutine MOM_write_field_legacy_0d(IO_handle, field_md, field, tstamp, fill_value, scale, unscale, zero_zeros)
  type(file_type),        intent(inout) :: IO_handle  !< Handle for a file that is open for writing
  type(fieldtype),        intent(in)    :: field_md   !< Field type with metadata
  real,                   intent(in)    :: field      !< Field to write in arbitrary units [A ~> a]
  real,         optional, intent(in)    :: tstamp     !< Model timestamp, often in [days]
  real,         optional, intent(in)    :: fill_value !< Missing data fill value [a]
  real,         optional, intent(in)    :: scale      !< A scaling factor that the field is multiplied by before
                                                      !! it is written [a A-1 ~> 1], for example to convert it
                                                      !! from its internal units to the desired units for output
  real,         optional, intent(in)    :: unscale    !< A scaling factor that the field is multiplied by before
                                                      !! it is written [a A-1 ~> 1], for example to convert it
                                                      !! from its internal units to the desired units for output.
                                                      !! Here scale and unscale are synonymous, but unscale
                                                      !! takes precedence if both are present.
  logical,      optional, intent(in)    :: zero_zeros !< If present and true, convert negative zeros
                                                      !! into ordinary signless zeros.

  ! Local variables

end subroutine MOM_write_field_legacy_0d
module subroutine MOM_write_field_4d(IO_handle, field_md, MOM_domain, field, tstamp, tile_count, &
                              fill_value, turns, scale, unscale, zero_zeros)
  class(MOM_file),          intent(inout) :: IO_handle  !< Handle for a file that is open for writing
  type(MOM_field),          intent(in)    :: field_md   !< Field type with metadata
  type(MOM_domain_type),    intent(in)    :: MOM_domain !< The MOM_Domain that describes the decomposition
  real, dimension(:,:,:,:), intent(inout) :: field      !< Unrotated field to write in arbitrary units [A ~> a]
  real,           optional, intent(in)    :: tstamp     !< Model timestamp, often in [days]
  integer,        optional, intent(in)    :: tile_count !< PEs per tile (default: 1)
  real,           optional, intent(in)    :: fill_value !< Missing data fill value [a]
  integer,        optional, intent(in)    :: turns      !< Number of quarter-turns to rotate the data
  real,           optional, intent(in)    :: scale      !< A scaling factor that the field is multiplied by before
                                                        !! it is written [a A-1 ~> 1], for example to convert it
                                                        !! from its internal units to the desired units for output
  real,           optional, intent(in)    :: unscale    !< A scaling factor that the field is multiplied by before
                                                        !! it is written [a A-1 ~> 1], for example to convert it
                                                        !! from its internal units to the desired units for output.
                                                        !! Here scale and unscale are synonymous, but unscale
                                                        !! takes precedence if both are present.
  logical,        optional, intent(in)    :: zero_zeros !< If present and true, convert negative zeros
                                                        !! into ordinary signless zeros.

  ! Local variables

end subroutine MOM_write_field_4d
module subroutine MOM_write_field_3d(IO_handle, field_md, MOM_domain, field, tstamp, tile_count, &
                              fill_value, turns, scale, unscale, zero_zeros)
  class(MOM_file),        intent(inout) :: IO_handle  !< Handle for a file that is open for writing
  type(MOM_field),        intent(in)    :: field_md   !< Field type with metadata
  type(MOM_domain_type),  intent(in)    :: MOM_domain !< The MOM_Domain that describes the decomposition
  real, dimension(:,:,:), intent(inout) :: field      !< Unrotated field to write in arbitrary units [A ~> a]
  real,         optional, intent(in)    :: tstamp     !< Model timestamp, often in [days]
  integer,      optional, intent(in)    :: tile_count !< PEs per tile (default: 1)
  real,         optional, intent(in)    :: fill_value !< Missing data fill value [a]
  integer,      optional, intent(in)    :: turns      !< Number of quarter-turns to rotate the data
  real,         optional, intent(in)    :: scale      !< A scaling factor that the field is multiplied by before
                                                      !! it is written [a A-1 ~> 1], for example to convert it
                                                      !! from its internal units to the desired units for output
  real,         optional, intent(in)    :: unscale    !< A scaling factor that the field is multiplied by before
                                                      !! it is written [a A-1 ~> 1], for example to convert it
                                                      !! from its internal units to the desired units for output.
                                                      !! Here scale and unscale are synonymous, but unscale
                                                      !! takes precedence if both are present.
  logical,      optional, intent(in)    :: zero_zeros !< If present and true, convert negative zeros
                                                      !! into ordinary signless zeros.

  ! Local variables

end subroutine MOM_write_field_3d
module subroutine MOM_write_field_2d(IO_handle, field_md, MOM_domain, field, tstamp, tile_count, &
                              fill_value, turns, scale, unscale, zero_zeros)
  class(MOM_file),        intent(inout) :: IO_handle  !< Handle for a file that is open for writing
  type(MOM_field),        intent(in)    :: field_md   !< Field type with metadata
  type(MOM_domain_type),  intent(in)    :: MOM_domain !< The MOM_Domain that describes the decomposition
  real, dimension(:,:),   intent(inout) :: field      !< Unrotated field to write in arbitrary units [A ~> a]
  real,         optional, intent(in)    :: tstamp     !< Model timestamp, often in [days]
  integer,      optional, intent(in)    :: tile_count !< PEs per tile (default: 1)
  real,         optional, intent(in)    :: fill_value !< Missing data fill value [a]
  integer,      optional, intent(in)    :: turns      !< Number of quarter-turns to rotate the data
  real,         optional, intent(in)    :: scale      !< A scaling factor that the field is multiplied by before
                                                      !! it is written [a A-1 ~> 1], for example to convert it
                                                      !! from its internal units to the desired units for output
  real,         optional, intent(in)    :: unscale    !< A scaling factor that the field is multiplied by before
                                                      !! it is written [a A-1 ~> 1], for example to convert it
                                                      !! from its internal units to the desired units for output.
                                                      !! Here scale and unscale are synonymous, but unscale
                                                      !! takes precedence if both are present.
  logical,      optional, intent(in)    :: zero_zeros !< If present and true, convert negative zeros
                                                      !! into ordinary signless zeros.

  ! Local variables

end subroutine MOM_write_field_2d
module subroutine MOM_write_field_1d(IO_handle, field_md, field, tstamp, fill_value, scale, unscale, zero_zeros)
  class(MOM_file),        intent(inout) :: IO_handle  !< Handle for a file that is open for writing
  type(MOM_field),        intent(in)    :: field_md   !< Field type with metadata
  real, dimension(:),     intent(in)    :: field      !< Field to write in arbitrary units [A ~> a]
  real,         optional, intent(in)    :: tstamp     !< Model timestamp, often in [days]
  real,         optional, intent(in)    :: fill_value !< Missing data fill value [a]
  real,         optional, intent(in)    :: scale      !< A scaling factor that the field is multiplied by before
                                                      !! it is written [a A-1 ~> 1], for example to convert it
                                                      !! from its internal units to the desired units for output
  real,         optional, intent(in)    :: unscale    !< A scaling factor that the field is multiplied by before
                                                      !! it is written [a A-1 ~> 1], for example to convert it
                                                      !! from its internal units to the desired units for output.
                                                      !! Here scale and unscale are synonymous, but unscale
                                                      !! takes precedence if both are present.
  logical,      optional, intent(in)    :: zero_zeros !< If present and true, convert negative zeros
                                                      !! into ordinary signless zeros.

  ! Local variables

end subroutine MOM_write_field_1d
module subroutine MOM_write_field_0d(IO_handle, field_md, field, tstamp, fill_value, scale, unscale, zero_zeros)
  class(MOM_file),        intent(inout) :: IO_handle  !< Handle for a file that is open for writing
  type(MOM_field),        intent(in)    :: field_md   !< Field type with metadata
  real,                   intent(in)    :: field      !< Field to write in arbitrary units [A ~> a]
  real,         optional, intent(in)    :: tstamp     !< Model timestamp, often in [days]
  real,         optional, intent(in)    :: fill_value !< Missing data fill value [a]
  real,         optional, intent(in)    :: scale      !< A scaling factor that the field is multiplied by before
                                                      !! it is written [a A-1 ~> 1], for example to convert it
                                                      !! from its internal units to the desired units for output
  real,         optional, intent(in)    :: unscale    !< A scaling factor that the field is multiplied by before
                                                      !! it is written [a A-1 ~> 1], for example to convert it
                                                      !! from its internal units to the desired units for output.
                                                      !! Here scale and unscale are synonymous, but unscale
                                                      !! takes precedence if both are present.
  logical,      optional, intent(in)    :: zero_zeros !< If present and true, convert negative zeros
                                                      !! into ordinary signless zeros.

  ! Local variables

end subroutine MOM_write_field_0d
module subroutine field_size(filename, fieldname, sizes, field_found, no_domain, ndims, ncid_in)
  character(len=*),      intent(in)    :: filename  !< The name of the file to read
  character(len=*),      intent(in)    :: fieldname !< The name of the variable whose sizes are returned
  integer, dimension(:), intent(inout) :: sizes     !< The sizes of the variable in each dimension
  logical,     optional, intent(out)   :: field_found !< This indicates whether the field was found in
                                                    !! the input file.  Without this argument, there
                                                    !! is a fatal error if the field is not found.
  logical,     optional, intent(in)    :: no_domain !< If present and true, do not check for file
                                                    !! names with an appended tile number.  If
                                                    !! ndims is present, the default changes to true.
  integer,     optional, intent(out)   :: ndims     !< The number of dimensions to the variable
  integer,     optional, intent(in)    :: ncid_in   !< The netCDF ID of an open file.  If absent, the
                                                    !! file is opened and closed within this routine.

end subroutine field_size
module subroutine safe_string_copy(str1, str2, fieldnm, caller)
  character(len=*),           intent(in)  :: str1    !< The string being copied
  character(len=*),           intent(out) :: str2    !< The string being copied into
  character(len=*), optional, intent(in)  :: fieldnm !< The name of the field for error messages
  character(len=*), optional, intent(in)  :: caller  !< The calling routine for error messages

end subroutine safe_string_copy
module function ensembler(name, ens_no_in) result(en_nm)
  character(len=*),  intent(in) :: name       !< The name to be modified
  integer, optional, intent(in) :: ens_no_in  !< The number of the current ensemble member
  character(len=len(name)) :: en_nm  !< The name encoded with the ensemble number

  ! This function replaces "%#E" or "%E" with the ensemble number anywhere it
  ! occurs in name, with %E using 4 or 6 digits (depending on the ensemble size)
  ! and %#E using # digits, where # is a number from 1 to 9.


end function ensembler
module subroutine get_filename_appendix(suffix)
  character(len=*), intent(out) :: suffix !< A string to append to filenames

end subroutine get_filename_appendix
module subroutine write_version_number(version, tag, unit)
  character(len=*),           intent(in) :: version !< A string that contains the routine name and version
  character(len=*), optional, intent(in) :: tag  !< A tag name to add to the message
  integer,          optional, intent(in) :: unit !< An alternate unit number for output

end subroutine write_version_number
module function open_namelist_file(file) result(unit)
  character(len=*), optional, intent(in) :: file !< The file to open, by default "input.nml"
  integer                                :: unit !< The opened unit number of the namelist file
end function open_namelist_file
module function check_nml_error(IOstat, nml_name) result(ierr)
  integer,          intent(in) :: IOstat   !< An I/O status field from a namelist read call
  character(len=*), intent(in) :: nml_name !< The name of the namelist
  integer :: ierr    !< A copy of IOstat that is returned to preserve legacy function behavior
end function check_nml_error
module subroutine MOM_io_init(param_file)
  type(param_file_type), intent(in) :: param_file  !< structure indicating the open file to
                                                   !! parse for model parameter values.

  ! This include declares and sets the variable "version".

end subroutine MOM_io_init
module subroutine get_var_axes_info(filename, fieldname, axes_info)
  character(len=*), intent(in) ::            filename  !< A filename from which to read
  character(len=*), intent(in) ::            fieldname !< The name of the field to read
  type(axis_info), dimension(4), intent(inout) :: axes_info !< A returned array of field axis information

  !! local variables
  !! cartesian axis data


end subroutine get_var_axes_info
  end interface

end module MOM_io
