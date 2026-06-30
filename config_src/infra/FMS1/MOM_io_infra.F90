! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> This module contains a thin inteface to mpp and fms I/O code
module MOM_io_infra

use MOM_domain_infra,     only : MOM_domain_type, rescale_comp_data, AGRID, BGRID_NE, CGRID_NE
use MOM_domain_infra,     only : domain2d, domain1d, CENTER, CORNER, NORTH_FACE, EAST_FACE
use MOM_error_infra,      only : MOM_error=>MOM_err, NOTE, FATAL, WARNING

use fms_mod,              only : write_version_number, open_namelist_file, check_nml_error
use fms_io_mod,           only : file_exist, field_exist, field_size, read_data
use fms_io_mod,           only : fms_io_exit, get_filename_appendix
use mpp_io_mod,           only : mpp_open, mpp_close, mpp_flush
use mpp_io_mod,           only : mpp_write_meta, mpp_write, mpp_read
use mpp_io_mod,           only : mpp_get_atts, mpp_attribute_exist
use mpp_io_mod,           only : mpp_get_axes, axistype, mpp_get_axis_data
use mpp_io_mod,           only : mpp_get_axis_length
use mpp_io_mod,           only : mpp_get_fields, fieldtype
use mpp_io_mod,           only : mpp_get_info, mpp_get_times
use mpp_io_mod,           only : mpp_io_init
use mpp_mod,              only : stdout_if_root=>stdout
! These are encoding constants.
use mpp_io_mod,           only : APPEND_FILE=>MPP_APPEND, WRITEONLY_FILE=>MPP_WRONLY
use mpp_io_mod,           only : OVERWRITE_FILE=>MPP_OVERWR, READONLY_FILE=>MPP_RDONLY
use mpp_io_mod,           only : NETCDF_FILE=>MPP_NETCDF, ASCII_FILE=>MPP_ASCII
use mpp_io_mod,           only : MULTIPLE=>MPP_MULTI, SINGLE_FILE=>MPP_SINGLE
use mpp_mod,              only : lowercase
use iso_fortran_env,      only : int64

implicit none ; private

! These interfaces are actually implemented or have explicit interfaces in this file.
public :: open_file, open_ASCII_file, file_is_open, close_file, flush_file, file_exists
public :: get_file_info, get_file_fields, get_file_times, get_filename_suffix
public :: read_field, read_vector, write_metadata, write_field
public :: field_exists, get_field_atts, get_field_size, read_field_chksum
public :: get_axis_size, get_axis_data, set_axis_data
public :: io_infra_init, io_infra_end, MOM_namelist_file, check_namelist_error, write_version
public :: stdout_if_root
! These types are inherited from underlying infrastructure code, to act as containers for
! information about fields and axes, respectively, and are opaque to this module.
public :: fieldtype, axistype
! These are encoding constant parmeters.
public :: ASCII_FILE, NETCDF_FILE, SINGLE_FILE, MULTIPLE
public :: APPEND_FILE, READONLY_FILE, OVERWRITE_FILE, WRITEONLY_FILE
public :: CENTER, CORNER, NORTH_FACE, EAST_FACE

!> Indicate whether a file exists, perhaps with domain decomposition
interface file_exists
  module procedure FMS_file_exists
  module procedure MOM_file_exists
end interface

!> Open a file (or fileset) for parallel or single-file I/).
interface open_file
  module procedure open_file_type, open_file_unit
end interface open_file

!> Read a data field from a file
interface read_field
  module procedure read_field_4d
  module procedure read_field_3d, read_field_3d_region
  module procedure read_field_2d, read_field_2d_region
  module procedure read_field_1d, read_field_1d_int
  module procedure read_field_0d, read_field_0d_int
end interface read_field

!> Write a registered field to an output file
interface write_field
  module procedure write_field_4d
  module procedure write_field_3d
  module procedure write_field_2d
  module procedure write_field_1d
  module procedure write_field_0d
  module procedure MOM_write_axis
end interface write_field

!> Read a pair of data fields representing the two components of a vector from a file
interface read_vector
  module procedure MOM_read_vector_3d
  module procedure MOM_read_vector_2d
end interface read_vector

!> Write metadata about a variable or axis to a file and store it for later reuse
interface write_metadata
  module procedure write_metadata_axis, write_metadata_field, write_metadata_global
end interface write_metadata

!> Close a file (or fileset).  If the file handle does not point to an open file,
!! close_file simply returns without doing anything.
interface close_file
  module procedure close_file_type, close_file_unit
end interface close_file

!> Ensure that the output stream associated with a file handle is fully sent to disk
interface flush_file
  module procedure flush_file_type, flush_file_unit
end interface flush_file

!> Type for holding a handle to an open file and related information
type, public :: file_type ; private
  integer :: unit = -1 !< The framework identfier or netCDF unit number of an output file
  character(len=:), allocatable :: filename !< The path to this file, if it is open
  logical :: open_to_read  = .false. !< If true, this file or fileset can be read
  logical :: open_to_write = .false. !< If true, this file or fileset can be written to
end type file_type


  interface
module subroutine read_field_chksum(field, chksum, valid_chksum)
  type(fieldtype),     intent(in)  :: field !< The field whose checksum attribute is to be read.
  integer(kind=int64), intent(out) :: chksum !< The checksum for the field.
  logical,             intent(out) :: valid_chksum  !< If true, chksum has been successfully read.
  ! Local variables

end subroutine read_field_chksum
logical module function MOM_file_exists(filename, MOM_Domain)
  character(len=*),       intent(in) :: filename   !< The name of the file being inquired about
  type(MOM_domain_type),  intent(in) :: MOM_Domain !< The MOM_Domain that describes the decomposition

! This function uses the fms_io function file_exist to determine whether
! a named file (or its decomposed variant) exists.

end function MOM_file_exists
logical module function FMS_file_exists(filename)
  character(len=*),         intent(in) :: filename  !< The name of the file being inquired about
  ! This function uses the fms_io function file_exist to determine whether
  ! a named file (or its decomposed variant) exists.

end function FMS_file_exists
logical module function file_is_open(IO_handle)
  type(file_type), intent(in) :: IO_handle !< Handle to a file to inquire about

end function file_is_open
module subroutine close_file_type(IO_handle)
  type(file_type), intent(inout) :: IO_handle   !< The I/O handle for the file to be closed

end subroutine close_file_type
module subroutine close_file_unit(unit)
  integer, intent(inout) :: unit   !< The I/O unit for the file to be closed

end subroutine close_file_unit
module subroutine flush_file_type(file)
  type(file_type), intent(in) :: file    !< The I/O handle for the file to flush

end subroutine flush_file_type
module subroutine flush_file_unit(unit)
  integer, intent(in) :: unit    !< The I/O unit for the file to flush

end subroutine flush_file_unit
module subroutine io_infra_init(maxunits)
  integer,   optional, intent(in) :: maxunits !< An optional maximum number of file
                                              !! unit numbers that can be used.
end subroutine io_infra_init
module subroutine io_infra_end()
end subroutine io_infra_end
module function MOM_namelist_file(file) result(unit)
  character(len=*), optional, intent(in) :: file !< The file to open, by default "input.nml".
  integer                                :: unit !< The opened unit number of the namelist file
end function MOM_namelist_file
module subroutine check_namelist_error(IOstat, nml_name)
  integer,          intent(in) :: IOstat   !< An I/O status field from a namelist read call
  character(len=*), intent(in) :: nml_name !< The name of the namelist
end subroutine check_namelist_error
module subroutine write_version(version, tag, unit)
  character(len=*),           intent(in) :: version !< A string that contains the routine name and version
  character(len=*), optional, intent(in) :: tag  !< A tag name to add to the message
  integer,          optional, intent(in) :: unit !< An alternate unit number for output

end subroutine write_version
module subroutine open_file_unit(unit, filename, action, form, threading, fileset, nohdrs, domain, MOM_domain)
  integer,                  intent(out) :: unit   !< The I/O unit for the opened file
  character(len=*),         intent(in)  :: filename !< The name of the file being opened
  integer,        optional, intent(in)  :: action !< A flag indicating whether the file can be read
                                                  !! or written to and how to handle existing files.
  integer,        optional, intent(in)  :: form   !< A flag indicating the format of a new file.  The
                                                  !! default is ASCII_FILE, but NETCDF_FILE is also common.
  integer,        optional, intent(in)  :: threading !< A flag indicating whether one (SINGLE_FILE)
                                                  !! or multiple PEs (MULTIPLE) participate in I/O.
                                                  !! With the default, the root PE does I/O.
  integer,        optional, intent(in)  :: fileset !< A flag indicating whether multiple PEs doing I/O due
                                                  !! to threading=MULTIPLE write to the same file (SINGLE_FILE)
                                                  !! or to one file per PE (MULTIPLE, the default).
  logical,        optional, intent(in)  :: nohdrs !< If nohdrs is .TRUE., headers are not written to
                                                  !! ASCII files.  The default is .false.
  type(domain2d), optional, intent(in)  :: domain !< A domain2d type that describes the decomposition
  type(MOM_domain_type), optional, intent(in) :: MOM_Domain !< A MOM_Domain that describes the decomposition

end subroutine open_file_unit
module subroutine open_file_type(IO_handle, filename, action, MOM_domain, threading, fileset)
  type(file_type),          intent(inout) :: IO_handle !< The handle for the opened file
  character(len=*),         intent(in)    :: filename !< The path name of the file being opened
  integer,        optional, intent(in)    :: action !< A flag indicating whether the file can be read
                                                    !! or written to and how to handle existing files.
                                                    !! The default is WRITE_ONLY.
  type(MOM_domain_type), &
                  optional, intent(in)    :: MOM_Domain !< A MOM_Domain that describes the decomposition
  integer,        optional, intent(in)    :: threading !< A flag indicating whether one (SINGLE_FILE)
                                                    !! or multiple PEs (MULTIPLE) participate in I/O.
                                                    !! With the default, the root PE does I/O.
  integer,        optional, intent(in)    :: fileset !< A flag indicating whether multiple PEs doing I/O due
                                                    !! to threading=MULTIPLE write to the same file (SINGLE_FILE)
                                                    !! or to one file per PE (MULTIPLE, the default).

end subroutine open_file_type
module subroutine open_ASCII_file(unit, file, action, threading, fileset)
  integer,                  intent(out) :: unit   !< The I/O unit for the opened file
  character(len=*),         intent(in)  :: file   !< The name of the file being opened
  integer,        optional, intent(in)  :: action !< A flag indicating whether the file can be read
                                                  !! or written to and how to handle existing files.
  integer,        optional, intent(in)  :: threading !< A flag indicating whether one (SINGLE_FILE)
                                                  !! or multiple PEs (MULTIPLE) participate in I/O.
                                                  !! With the default, the root PE does I/O.
  integer,        optional, intent(in)  :: fileset !< A flag indicating whether multiple PEs doing I/O due
                                                  !! to threading=MULTIPLE write to the same file (SINGLE_FILE)
                                                  !! or to one file per PE (MULTIPLE, the default).

end subroutine open_ASCII_file
module subroutine get_filename_suffix(suffix)
  character(len=*), intent(out) :: suffix !< A string to append to filenames

end subroutine get_filename_suffix
module subroutine get_file_info(IO_handle, ndim, nvar, ntime)
  type(file_type),    intent(in)  :: IO_handle !< Handle for a file that is open for I/O
  integer,  optional, intent(out) :: ndim  !< The number of dimensions in the file
  integer,  optional, intent(out) :: nvar  !< The number of variables in the file
  integer,  optional, intent(out) :: ntime !< The number of time levels in the file

  ! Local variables

end subroutine get_file_info
module subroutine get_file_times(IO_handle, time_values, ntime)
  type(file_type),                 intent(in)    :: IO_handle !< Handle for a file that is open for I/O
  real, allocatable, dimension(:), intent(inout) :: time_values !< The real times for the records in file.
  integer,               optional, intent(out)   :: ntime !< The number of time levels in the file


end subroutine get_file_times
module subroutine get_file_fields(IO_handle, fields)
  type(file_type),               intent(in)    :: IO_handle !< Handle for a file that is open for I/O
  type(fieldtype), dimension(:), intent(inout) :: fields !< Field-type descriptions of all of
                                                         !! the variables in a file.
end subroutine get_file_fields
module subroutine get_field_atts(field, name, units, longname, checksum)
  type(fieldtype),            intent(in)  :: field !< The field to extract information from
  character(len=*), optional, intent(out) :: name  !< The variable name
  character(len=*), optional, intent(out) :: units !< The units of the variable
  character(len=*), optional, intent(out) :: longname  !< The long name of the variable
  integer(kind=int64),  dimension(:), &
                    optional, intent(out) :: checksum !< The checksums of the variable in a file
end subroutine get_field_atts
module function field_exists(filename, field_name, domain, no_domain, MOM_domain)
  character(len=*),                 intent(in) :: filename   !< The name of the file being inquired about
  character(len=*),                 intent(in) :: field_name !< The name of the field being sought
  type(domain2d), target, optional, intent(in) :: domain     !< A domain2d type that describes the decomposition
  logical,                optional, intent(in) :: no_domain  !< This file does not use domain decomposition
  type(MOM_domain_type),  optional, intent(in) :: MOM_Domain !< A MOM_Domain that describes the decomposition
  logical                                      :: field_exists !< True if filename exists and field_name is in filename

end function field_exists
module subroutine get_field_size(filename, fieldname, sizes, field_found, no_domain)
  character(len=*),      intent(in)    :: filename  !< The name of the file to read
  character(len=*),      intent(in)    :: fieldname !< The name of the variable whose sizes are returned
  integer, dimension(:), intent(inout) :: sizes     !< The sizes of the variable in each dimension
  logical,     optional, intent(out)   :: field_found !< This indicates whether the field was found in
                                                    !! the input file.  Without this argument, there
                                                    !! is a fatal error if the field is not found.
  logical,     optional, intent(in)    :: no_domain !< If present and true, do not check for file
                                                    !! names with an appended tile number

end subroutine get_field_size
module function get_axis_size(axis) result(axis_size)
  type(axistype), intent(in) :: axis
    !< Infra axis
  integer :: axis_size
    !< Axis size

end function get_axis_size
module subroutine get_axis_data(axis, axis_name, axis_data)
  type(axistype), intent(in) :: axis
    !< Infra axis
  character(len=256), intent(out) :: axis_name
    !< Axis name
  real, dimension(:), intent(out) :: axis_data
    !< Axis points

end subroutine get_axis_data
module subroutine set_axis_data(axis, axis_name, axis_data)
  type(axistype), intent(inout) :: axis
    !< Target axis
  character(len=256), intent(in) :: axis_name
    !< Target axis name
  real, intent(in) :: axis_data(:)
    !< Target axis values

end subroutine set_axis_data
module subroutine read_field_0d(filename, fieldname, data, timelevel, scale, MOM_Domain, &
                         global_file, file_may_be_4d)
  character(len=*),       intent(in)    :: filename  !< The name of the file to read
  character(len=*),       intent(in)    :: fieldname !< The variable name of the data in the file
  real,                   intent(inout) :: data      !< The 1-dimensional array into which the data
  integer,      optional, intent(in)    :: timelevel !< The time level in the file to read
  real,         optional, intent(in)    :: scale     !< A scaling factor that the field is multiplied
                                                     !! by before it is returned.
  type(MOM_domain_type), &
                optional, intent(in)    :: MOM_Domain !< The MOM_Domain that describes the decomposition
  logical,      optional, intent(in)    :: global_file !< If true, read from a single global file
  logical,      optional, intent(in)    :: file_may_be_4d !< If true, this file may have 4-d arrays,
                                                     !! in which case a more elaborate set of calls
                                                     !! is needed to read it due to FMS limitations.

  ! Local variables

end subroutine read_field_0d
module subroutine read_field_1d(filename, fieldname, data, timelevel, scale, MOM_Domain, &
                            global_file, file_may_be_4d)
  character(len=*),       intent(in)    :: filename  !< The name of the file to read
  character(len=*),       intent(in)    :: fieldname !< The variable name of the data in the file
  real, dimension(:),     intent(inout) :: data      !< The 1-dimensional array into which the data
  integer,      optional, intent(in)    :: timelevel !< The time level in the file to read
  real,         optional, intent(in)    :: scale     !< A scaling factor that the field is multiplied
                                                     !! by before they are returned.
  type(MOM_domain_type), &
                optional, intent(in)    :: MOM_Domain !< The MOM_Domain that describes the decomposition
  logical,      optional, intent(in)    :: global_file !< If true, read from a single global file
  logical,      optional, intent(in)    :: file_may_be_4d !< If true, this file may have 4-d arrays,
                                                     !! in which case a more elaborate set of calls
                                                     !! is needed to read it due to FMS limitations.

  ! Local variables

end subroutine read_field_1d
module subroutine read_field_2d(filename, fieldname, data, MOM_Domain, &
                         timelevel, position, scale, global_file, file_may_be_4d)
  character(len=*),       intent(in)    :: filename  !< The name of the file to read
  character(len=*),       intent(in)    :: fieldname !< The variable name of the data in the file
  real, dimension(:,:),   intent(inout) :: data      !< The 2-dimensional array into which the data
                                                     !! should be read
  type(MOM_domain_type),  intent(in)    :: MOM_Domain !< The MOM_Domain that describes the decomposition
  integer,      optional, intent(in)    :: timelevel !< The time level in the file to read
  integer,      optional, intent(in)    :: position  !< A flag indicating where this data is located
  real,         optional, intent(in)    :: scale     !< A scaling factor that the field is multiplied
                                                     !! by before it is returned.
  logical,      optional, intent(in)    :: global_file !< If true, read from a single global file
  logical,      optional, intent(in)    :: file_may_be_4d !< If true, this file may have 4-d arrays,
                                                     !! in which case a more elaborate set of calls
                                                     !! is needed to read it due to FMS limitations.

  ! Local variables

end subroutine read_field_2d
module subroutine read_field_2d_region(filename, fieldname, data, start, nread, MOM_domain, &
                                no_domain, scale)
  character(len=*),       intent(in)    :: filename  !< The name of the file to read
  character(len=*),       intent(in)    :: fieldname !< The variable name of the data in the file
  real, dimension(:,:),   intent(inout) :: data      !< The 2-dimensional array into which the data
                                                     !! should be read
  integer, dimension(:),  intent(in)    :: start     !< The starting index to read in each of 4
                                                     !! dimensions.  For this 2-d read, the 3rd
                                                     !! and 4th values are always 1.
  integer, dimension(:),  intent(in)    :: nread     !< The number of points to read in each of 4
                                                     !! dimensions.  For this 2-d read, the 3rd
                                                     !! and 4th values are always 1.
  type(MOM_domain_type), &
                optional, intent(in)    :: MOM_Domain !< The MOM_Domain that describes the decomposition
  logical,      optional, intent(in)    :: no_domain !< If present and true, this variable does not
                                                     !! use domain decomposion.
  real,         optional, intent(in)    :: scale     !< A scaling factor that the field is multiplied
                                                     !! by before it is returned.

end subroutine read_field_2d_region
module subroutine read_field_3d(filename, fieldname, data, MOM_Domain, &
                            timelevel, position, scale, global_file, file_may_be_4d)
  character(len=*),       intent(in)    :: filename  !< The name of the file to read
  character(len=*),       intent(in)    :: fieldname !< The variable name of the data in the file
  real, dimension(:,:,:), intent(inout) :: data      !< The 3-dimensional array into which the data
                                                     !! should be read
  type(MOM_domain_type),  intent(in)    :: MOM_Domain !< The MOM_Domain that describes the decomposition
  integer,      optional, intent(in)    :: timelevel !< The time level in the file to read
  integer,      optional, intent(in)    :: position  !< A flag indicating where this data is located
  real,         optional, intent(in)    :: scale     !< A scaling factor that the field is multiplied
                                                     !! by before it is returned.
  logical,      optional, intent(in)    :: global_file !< If true, read from a single global file
  logical,      optional, intent(in)    :: file_may_be_4d !< If true, this file may have 4-d arrays,
                                                     !! in which case a more elaborate set of calls
                                                     !! is needed to read it due to FMS limitations.

  ! Local variables

end subroutine read_field_3d
module subroutine read_field_3d_region(filename, fieldname, data, start, nread, MOM_domain, &
                                no_domain, scale)
  character(len=*),       intent(in)    :: filename  !< The name of the file to read
  character(len=*),       intent(in)    :: fieldname !< The variable name of the data in the file
  real, dimension(:,:,:),   intent(inout) :: data    !< The 3-dimensional array into which the data
                                                     !! should be read
  integer, dimension(:),  intent(in)    :: start     !< The starting index to read in each of 4
                                                     !! dimensions.  For this 3-d read, the
                                                     !! 4th values are always 1.
  integer, dimension(:),  intent(in)    :: nread     !< The number of points to read in each of 4
                                                     !! dimensions.  For this 3-d read, the
                                                     !! 4th values are always 1.
  type(MOM_domain_type), &
                optional, intent(in)    :: MOM_Domain !< The MOM_Domain that describes the decomposition
  logical,      optional, intent(in)    :: no_domain !< If present and true, this variable does not
                                                     !! use domain decomposion.
  real,         optional, intent(in)    :: scale     !< A scaling factor that the field is multiplied
                                                     !! by before it is returned.

end subroutine read_field_3d_region
module subroutine read_field_4d(filename, fieldname, data, MOM_Domain, &
                            timelevel, position, scale, global_file)
  character(len=*),       intent(in)    :: filename  !< The name of the file to read
  character(len=*),       intent(in)    :: fieldname !< The variable name of the data in the file
  real, dimension(:,:,:,:), intent(inout) :: data    !< The 4-dimensional array into which the data
                                                     !! should be read
  type(MOM_domain_type),  intent(in)    :: MOM_Domain !< The MOM_Domain that describes the decomposition
  integer,      optional, intent(in)    :: timelevel !< The time level in the file to read
  integer,      optional, intent(in)    :: position  !< A flag indicating where this data is located
  real,         optional, intent(in)    :: scale     !< A scaling factor that the field is multiplied
                                                     !! by before it is returned.
  logical,      optional, intent(in)    :: global_file !< If true, read from a single global file

  ! Local variables

  ! This single call does not work for a 4-d array due to FMS limitations, so multiple calls are
  ! needed.
  ! call read_data(filename, fieldname, data, MOM_Domain%mpp_domain, &
  !                timelevel=timelevel, position=position)

end subroutine read_field_4d
module subroutine read_field_0d_int(filename, fieldname, data, timelevel)
  character(len=*),       intent(in)    :: filename  !< The name of the file to read
  character(len=*),       intent(in)    :: fieldname !< The variable name of the data in the file
  integer,                intent(inout) :: data      !< The 1-dimensional array into which the data
  integer,      optional, intent(in)    :: timelevel !< The time level in the file to read

end subroutine read_field_0d_int
module subroutine read_field_1d_int(filename, fieldname, data, timelevel)
  character(len=*),       intent(in)    :: filename  !< The name of the file to read
  character(len=*),       intent(in)    :: fieldname !< The variable name of the data in the file
  integer, dimension(:),  intent(inout) :: data      !< The 1-dimensional array into which the data
  integer,      optional, intent(in)    :: timelevel !< The time level in the file to read

end subroutine read_field_1d_int
module subroutine MOM_read_vector_2d(filename, u_fieldname, v_fieldname, u_data, v_data, MOM_Domain, &
                              timelevel, stagger, scalar_pair, scale)
  character(len=*),       intent(in)    :: filename  !< The name of the file to read
  character(len=*),       intent(in)    :: u_fieldname !< The variable name of the u data in the file
  character(len=*),       intent(in)    :: v_fieldname !< The variable name of the v data in the file
  real, dimension(:,:),   intent(inout) :: u_data    !< The 2 dimensional array into which the
                                                     !! u-component of the data should be read
  real, dimension(:,:),   intent(inout) :: v_data    !< The 2 dimensional array into which the
                                                     !! v-component of the data should be read
  type(MOM_domain_type),  intent(in)    :: MOM_Domain !< The MOM_Domain that describes the decomposition
  integer,      optional, intent(in)    :: timelevel !< The time level in the file to read
  integer,      optional, intent(in)    :: stagger   !< A flag indicating where this vector is discretized
  logical,      optional, intent(in)    :: scalar_pair !< If true, a pair of scalars are to be read
  real,         optional, intent(in)    :: scale     !< A scaling factor that the fields are multiplied
                                                     !! by before they are returned.

end subroutine MOM_read_vector_2d
module subroutine MOM_read_vector_3d(filename, u_fieldname, v_fieldname, u_data, v_data, MOM_Domain, &
                              timelevel, stagger, scalar_pair, scale)
  character(len=*),       intent(in)    :: filename  !< The name of the file to read
  character(len=*),       intent(in)    :: u_fieldname !< The variable name of the u data in the file
  character(len=*),       intent(in)    :: v_fieldname !< The variable name of the v data in the file
  real, dimension(:,:,:), intent(inout) :: u_data    !< The 3 dimensional array into which the
                                                     !! u-component of the data should be read
  real, dimension(:,:,:), intent(inout) :: v_data    !< The 3 dimensional array into which the
                                                     !! v-component of the data should be read
  type(MOM_domain_type),  intent(in)    :: MOM_Domain !< The MOM_Domain that describes the decomposition
  integer,      optional, intent(in)    :: timelevel !< The time level in the file to read
  integer,      optional, intent(in)    :: stagger   !< A flag indicating where this vector is discretized
  logical,      optional, intent(in)    :: scalar_pair !< If true, a pair of scalars are to be read.cretized
  real,         optional, intent(in)    :: scale     !< A scaling factor that the fields are multiplied
                                                     !! by before they are returned.


end subroutine MOM_read_vector_3d
module subroutine write_field_4d(IO_handle, field_md, MOM_domain, field, tstamp, tile_count, fill_value)
  type(file_type),          intent(in)    :: IO_handle  !< Handle for a file that is open for writing
  type(fieldtype),          intent(in)    :: field_md   !< Field type with metadata
  type(MOM_domain_type),    intent(in)    :: MOM_domain !< The MOM_Domain that describes the decomposition
  real, dimension(:,:,:,:), intent(inout) :: field      !< Field to write
  real,           optional, intent(in)    :: tstamp     !< Model time of this field
  integer,        optional, intent(in)    :: tile_count !< PEs per tile (default: 1)
  real,           optional, intent(in)    :: fill_value !< Missing data fill value

end subroutine write_field_4d
module subroutine write_field_3d(IO_handle, field_md, MOM_domain, field, tstamp, tile_count, fill_value)
  type(file_type),        intent(in)    :: IO_handle  !< Handle for a file that is open for writing
  type(fieldtype),        intent(in)    :: field_md   !< Field type with metadata
  type(MOM_domain_type),  intent(in)    :: MOM_domain !< The MOM_Domain that describes the decomposition
  real, dimension(:,:,:), intent(inout) :: field      !< Field to write
  real,         optional, intent(in)    :: tstamp     !< Model time of this field
  integer,      optional, intent(in)    :: tile_count !< PEs per tile (default: 1)
  real,         optional, intent(in)    :: fill_value !< Missing data fill value

end subroutine write_field_3d
module subroutine write_field_2d(IO_handle, field_md, MOM_domain, field, tstamp, tile_count, fill_value)
  type(file_type),        intent(in)    :: IO_handle  !< Handle for a file that is open for writing
  type(fieldtype),        intent(in)    :: field_md   !< Field type with metadata
  type(MOM_domain_type),  intent(in)    :: MOM_domain !< The MOM_Domain that describes the decomposition
  real, dimension(:,:),   intent(inout) :: field      !< Field to write
  real,         optional, intent(in)    :: tstamp     !< Model time of this field
  integer,      optional, intent(in)    :: tile_count !< PEs per tile (default: 1)
  real,         optional, intent(in)    :: fill_value !< Missing data fill value

end subroutine write_field_2d
module subroutine write_field_1d(IO_handle, field_md, field, tstamp)
  type(file_type),        intent(in)    :: IO_handle  !< Handle for a file that is open for writing
  type(fieldtype),        intent(in)    :: field_md   !< Field type with metadata
  real, dimension(:),     intent(in)    :: field      !< Field to write
  real,         optional, intent(in)    :: tstamp     !< Model time of this field

end subroutine write_field_1d
module subroutine write_field_0d(IO_handle, field_md, field, tstamp)
  type(file_type),        intent(in)    :: IO_handle  !< Handle for a file that is open for writing
  type(fieldtype),        intent(in)    :: field_md   !< Field type with metadata
  real,                   intent(in)    :: field      !< Field to write
  real,         optional, intent(in)    :: tstamp     !< Model time of this field

end subroutine write_field_0d
module subroutine MOM_write_axis(IO_handle, axis)
  type(file_type), intent(in) :: IO_handle  !< Handle for a file that is open for writing
  type(axistype),  intent(in) :: axis       !< An axis type variable with information to write

end subroutine MOM_write_axis
module subroutine write_metadata_axis(IO_handle, axis, name, units, longname, cartesian, sense, domain, &
                               data, edge_axis, calendar)
  type(file_type),            intent(in)    :: IO_handle  !< Handle for a file that is open for writing
  type(axistype),             intent(inout) :: axis  !< The axistype where this information is stored.
  character(len=*),           intent(in)    :: name  !< The name in the file of this axis
  character(len=*),           intent(in)    :: units !< The units of this axis
  character(len=*),           intent(in)    :: longname !< The long description of this axis
  character(len=*), optional, intent(in)    :: cartesian !< A variable indicating which direction
                                                     !! this axis corresponds with. Valid values
                                                     !! include 'X', 'Y', 'Z', 'T', and 'N' for none.
  integer,          optional, intent(in)    :: sense !< This is 1 for axes whose values increase upward, or
                                                     !! -1 if they increase downward.
  type(domain1D),   optional, intent(in)    :: domain !< The domain decomposion for this axis
  real, dimension(:), optional, intent(in)  :: data   !< The coordinate values of the points on this axis
  logical,          optional, intent(in)    :: edge_axis !< If true, this axis marks an edge of the tracer cells
  character(len=*), optional, intent(in)    :: calendar !< The name of the calendar used with a time axis

end subroutine write_metadata_axis
module subroutine write_metadata_field(IO_handle, field, axes, name, units, longname, &
                                pack, standard_name, checksum)
  type(file_type),            intent(in)    :: IO_handle  !< Handle for a file that is open for writing
  type(fieldtype),            intent(inout) :: field !< The fieldtype where this information is stored
  type(axistype), dimension(:), intent(in)  :: axes  !< Handles for the axis used for this variable
  character(len=*),           intent(in)    :: name  !< The name in the file of this variable
  character(len=*),           intent(in)    :: units !< The units of this variable
  character(len=*),           intent(in)    :: longname !< The long description of this variable
  integer,          optional, intent(in)    :: pack  !< A precision reduction factor with which the
                                                     !! variable.  The default, 1, has no reduction,
                                                     !! but 2 is not uncommon.
  character(len=*), optional, intent(in)    :: standard_name !< The standard (e.g., CMOR) name for this variable
  integer(kind=int64), dimension(:), &
                    optional, intent(in)    :: checksum !< Checksum values that can be used to verify reads.


end subroutine write_metadata_field
module subroutine write_metadata_global(IO_handle, name, attribute)
  type(file_type),            intent(in)    :: IO_handle !< Handle for a file that is open for writing
  character(len=*),           intent(in)    :: name      !< The name in the file of this global attribute
  character(len=*),           intent(in)    :: attribute !< The value of this attribute

end subroutine write_metadata_global
  end interface

end module MOM_io_infra
