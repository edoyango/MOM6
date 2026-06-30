! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> MOM6 interface to netCDF operations
module MOM_netcdf

use, intrinsic :: iso_fortran_env, only : real32, real64

use netcdf, only : nf90_create, nf90_open, nf90_close
use netcdf, only : nf90_sync
use netcdf, only : NF90_CLOBBER, NF90_NOCLOBBER, NF90_WRITE, NF90_NOWRITE
use netcdf, only : nf90_enddef
use netcdf, only : nf90_def_dim, nf90_def_var
use netcdf, only : NF90_UNLIMITED
use netcdf, only : nf90_get_var
use netcdf, only : nf90_put_var, nf90_put_att
use netcdf, only : NF90_FLOAT, NF90_DOUBLE
use netcdf, only : nf90_strerror, NF90_NOERR
use netcdf, only : NF90_GLOBAL
use netcdf, only : nf90_inquire, nf90_inquire_dimension, nf90_inquire_variable
use netcdf, only : nf90_inq_dimids, nf90_inq_varids
use netcdf, only : NF90_MAX_NAME

use MOM_error_handler, only : MOM_error, FATAL
use MOM_io_infra, only : READONLY_FILE, WRITEONLY_FILE
use MOM_io_infra, only : APPEND_FILE, OVERWRITE_FILE

implicit none ; private

public :: netcdf_file_type
public :: netcdf_axis
public :: netcdf_field
public :: open_netcdf_file
public :: close_netcdf_file
public :: flush_netcdf_file
public :: register_netcdf_axis
public :: register_netcdf_field
public :: write_netcdf_field
public :: write_netcdf_axis
public :: write_netcdf_attribute
public :: get_netcdf_size
public :: get_netcdf_fields
public :: get_netcdf_filename
public :: read_netcdf_field


!> Internal time value used to indicate an uninitialized time
real, parameter :: NULLTIME = -1
! NOTE: For now, we use the FMS-compatible value, but may change in the future.


!> netCDF file abstraction
type :: netcdf_file_type
  private
  integer :: ncid
    !< netCDF file ID
  character(len=:), allocatable :: filename
    !< netCDF filename
  logical :: define_mode
    !< True if file is in define mode.
  integer :: time_id
    !< Time axis variable ID
  real :: time
    !< Current model time
  integer :: time_level
    !< Current time level for output
end type netcdf_file_type


!> Dimension axis for a netCDF file
type :: netcdf_axis
  private
  character(len=:), allocatable, public :: label
    !< Axis label name
  real, allocatable :: points(:)
    !< Grid points along the axis
  integer :: dimid
    !< netCDF dimension ID associated with axis
  integer :: varid
    !< netCDF variable ID associated with axis
end type netcdf_axis


!> Field variable for a netCDF file
type netcdf_field
  private
  character(len=:), allocatable, public :: label
    !< Variable name
  integer :: varid
    !< netCDF variable ID for field
end type netcdf_field


!> Write values to a field of a netCDF file
interface write_netcdf_field
  module procedure write_netcdf_field_4d
  module procedure write_netcdf_field_3d
  module procedure write_netcdf_field_2d
  module procedure write_netcdf_field_1d
  module procedure write_netcdf_field_0d
end interface write_netcdf_field


  interface
module subroutine open_netcdf_file(handle, filename, mode)
  type(netcdf_file_type), intent(inout) :: handle
    !< netCDF file handle
  character(len=*), intent(in) :: filename
    !< netCDF filename
  integer, intent(in), optional :: mode
    !< Input MOM I/O mode

    ! MOM I/O mode
    ! netCDF creation mode
    ! nf90_create return code
    ! netCDF error message buffer

  ! I/O configuration
end subroutine open_netcdf_file
module subroutine close_netcdf_file(handle)
  type(netcdf_file_type), intent(in) :: handle


end subroutine close_netcdf_file
module subroutine flush_netcdf_file(handle)
  type(netcdf_file_type), intent(in) :: handle


end subroutine flush_netcdf_file
module subroutine enable_netcdf_write(handle)
  type(netcdf_file_type), intent(inout) :: handle


end subroutine enable_netcdf_write
module function register_netcdf_field(handle, label, axes, longname, units) &
    result(field)
  type(netcdf_file_type), intent(in) :: handle
    !< netCDF file handle
  character(len=*), intent(in) :: label
    !< netCDF field name in the file
  type(netcdf_axis), intent(in) :: axes(:)
    !< Axes along which field is defined
  character(len=*), intent(in) :: longname
    !< Long name of the netCDF field
  character(len=*), intent(in) :: units
    !< Field units of measurement
  type(netcdf_field) :: field
    !< netCDF field

    ! netCDF function return code
    ! Loop index
    ! netCDF dimension IDs of axes
    ! netCDF data type

  ! Gather the axis netCDF dimension IDs
end function register_netcdf_field
module function register_netcdf_axis(handle, label, units, longname, points, &
    cartesian, sense) result(axis)
  type(netcdf_file_type), intent(inout) :: handle
    !< netCDF file handle
  character(len=*), intent(in) :: label
    !< netCDF axis name in the file
  character(len=*), intent(in), optional :: units
    !< Axis units of measurement
  character(len=*), intent(in), optional :: longname
    !< Long name of the axis
  real, intent(in), optional :: points(:)
    !< Values of axis points (for fixed axes)
  character(len=*), intent(in), optional :: cartesian
    !< Character denoting axis direction: X, Y, Z, T, or N for none
  integer, intent(in), optional :: sense
    !< Axis direction; +1 if axis increases upward or -1 if downward

  type(netcdf_axis) :: axis
    !< netCDF coordinate axis

    ! netCDF external data type
    ! netCDF function return code
    ! True if the axis is unlimited in size (e.g. time)
    ! Either the number of points in the axis, or unlimited flag
    ! Axis direction; +1 if axis increases upward or -1 if downward
    ! CF-compiant value of sense attribute (as 'positive')

  ! Create the axis dimension
end function register_netcdf_axis
module subroutine write_netcdf_field_4d(handle, field, values, time)
  type(netcdf_file_type), intent(inout) :: handle
    !< netCDF file handle
  type(netcdf_field), intent(in) :: field
    !< Field metadata
  real, intent(in) :: values(:,:,:,:)
    !< Field values
  real, intent(in), optional :: time
    !< Timestep index to write data

    ! netCDF return code
    ! Start indices, if timestep is included

  ! Verify write mode
end subroutine write_netcdf_field_4d
module subroutine write_netcdf_field_3d(handle, field, values, time)
  type(netcdf_file_type), intent(inout) :: handle
    !< netCDF file handle
  type(netcdf_field), intent(in) :: field
    !< Field metadata
  real, intent(in) :: values(:,:,:)
    !< Field values
  real, intent(in), optional :: time
    !< Timestep index to write data

    ! netCDF return code
    ! Start indices, if timestep is included

  ! Verify write mode
end subroutine write_netcdf_field_3d
module subroutine write_netcdf_field_2d(handle, field, values, time)
  type(netcdf_file_type), intent(inout) :: handle
    !< netCDF file handle
  type(netcdf_field), intent(in) :: field
    !< Field metadata
  real, intent(in) :: values(:,:)
    !< Field values
  real, intent(in), optional :: time
    !< Timestep index to write data

    ! netCDF return code
    ! Start indices, if timestep is included

  ! Verify write mode
end subroutine write_netcdf_field_2d
module subroutine write_netcdf_field_1d(handle, field, values, time)
  type(netcdf_file_type), intent(inout) :: handle
    !< netCDF file handle
  type(netcdf_field), intent(in) :: field
    !< Field metadata
  real, intent(in) :: values(:)
    !< Field values
  real, intent(in), optional :: time
    !< Timestep index to write data

    ! netCDF return code
    ! Start indices, if timestep is included

  ! Verify write mode
end subroutine write_netcdf_field_1d
module subroutine write_netcdf_field_0d(handle, field, scalar, time)
  type(netcdf_file_type), intent(inout) :: handle
    !< netCDF file handle
  type(netcdf_field), intent(in) :: field
    !< Field metadata
  real, intent(in) :: scalar
    !< Field values
  real, intent(in), optional :: time
    !< Timestep index to write data

    ! netCDF return code
    ! Start indices, if timestep is included

  ! Verify write mode
end subroutine write_netcdf_field_0d
module subroutine write_netcdf_axis(handle, axis)
  type(netcdf_file_type), intent(inout) :: handle
    !< netCDF file handle
  type(netcdf_axis), intent(in) :: axis
    !< field variable

    ! netCDF return code

  ! Verify write mode
end subroutine write_netcdf_axis
module subroutine write_netcdf_attribute(handle, label, attribute)
  type(netcdf_file_type), intent(in) :: handle
    !< netCDF file handle
  character(len=*), intent(in) :: label
    !< File attribute
  character(len=*), intent(in) :: attribute
    !< File attribute value

    ! netCDF return code

end subroutine write_netcdf_attribute
module subroutine get_netcdf_size(handle, ndims, nvars, nsteps)
  type(netcdf_file_type), intent(in) :: handle
    !< netCDF input file
  integer, intent(out), optional :: ndims
    !< number of dimensions in the file
  integer, intent(out), optional :: nvars
    !< number of variables in the file
  integer, intent(out), optional :: nsteps
    !< number of values in the file's unlimited axis

    ! netCDF return code
    ! netCDF dimension ID for unlimited time axis

end subroutine get_netcdf_size
module subroutine get_netcdf_fields(handle, axes, fields)
  type(netcdf_file_type), intent(inout) :: handle
    !< netCDF file handle
  type(netcdf_axis), intent(inout), allocatable :: axes(:)
    !< netCDF file axes
  type(netcdf_field), intent(inout), allocatable :: fields(:)
    !< netCDF file fields

    ! Number of netCDF dimensions
    ! Number of netCDF dimensions
    ! netCDF variables in handle
    ! Number of fields in the file (i.e. non-axis variables)
    ! netCDF dimension IDs of file
    ! netCDF variable IDs of file
    ! netCDF dimension ID for the unlimited axis variable, if present
    ! Index of the unlimited axis in axes(:), if present
    ! Current dimension or variable label
    ! Current dimension length
    ! netCDF return code
    ! Group-based counts for nf90_inq_* (unused)
    ! True if the current variable is an axis

    ! Flag indicating exclusion of parent groups in netCDF file
    ! NOTE: This must be passed as a variable, and cannot be declared as a
    !   parameter.

end subroutine get_netcdf_fields
module function get_netcdf_filename(handle)
  type(netcdf_file_type), intent(in) :: handle !< A netCDF file handle
  character(len=:), allocatable :: get_netcdf_filename !< The name of the file that this handle refers to.

end function get_netcdf_filename
module subroutine read_netcdf_field(handle, field, values, bounds)
  type(netcdf_file_type), intent(in) :: handle
  type(netcdf_field), intent(in) :: field
  real, intent(out) :: values(:,:)
  integer, optional, intent(in) :: bounds(2,2)

    ! netCDF return code
    ! Axis start index
    ! Axis index count

end subroutine read_netcdf_field
module subroutine update_netcdf_timestep(handle, time)
  type(netcdf_file_type), intent(inout) :: handle
    !< netCDF file handle
  real, intent(in) :: time
    !< New model time

    !< Time axis start index array
    !< netCDF return code

end subroutine update_netcdf_timestep
module subroutine check_netcdf_call(ncerr, header, message)
  integer, intent(in) :: ncerr
    !< netCDF error code
  character(len=*), intent(in) :: header
    !< Message header (usually calling subroutine)
  character(len=*), intent(in) :: message
    !< Error message (usually action which instigated the error)

    ! Full error message, including netCDF message

end subroutine check_netcdf_call
  end interface

end module MOM_netcdf
