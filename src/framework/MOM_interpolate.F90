! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> This module provides added functionality to the FMS temporal and spatial interpolation routines
module MOM_interpolate

use MOM_array_transform, only : allocate_rotated_array, rotate_array
use MOM_error_handler,   only : MOM_error, FATAL
use MOM_interp_infra,    only : time_interp_extern, init_external_field=>init_extern_field
use MOM_interp_infra,    only : time_interp_external_init=>time_interp_extern_init
use MOM_interp_infra,    only : horiz_interp_type
use MOM_interp_infra,    only : get_external_field_info_infra => get_external_field_info
use MOM_interp_infra,    only : run_horiz_interp, build_horiz_interp_weights
use MOM_interp_infra,    only : external_field
use MOM_io_infra,        only : axistype
use MOM_io_infra,        only : get_axis_size, get_axis_data
use MOM_io,              only : axis_info, set_axis_info
use MOM_time_manager, only : time_type, set_date, operator(+), operator(<), operator(>)

implicit none ; private

!> Data type used to store information about forcing datasets that are time series
!! E.g. how do we align the data in the model with the time axis in the file?
type, public :: forcing_timeseries_dataset
    character(len=200) :: file_name  !< name of file containing river flux forcing
    logical :: l_time_varying        !< .true. => forcing is dependent on model time, .false. => static forcing
    ! logical :: l_FMS_modulo        !< .true. => let FMS handle determining time level to read (e.g. for climatologies)
    type(time_type) :: data_forcing  !< convert data_forcing_year to time type
    type(time_type) :: data_start    !< convert data_start_year to time type
    type(time_type) :: data_end      !< convert data_end_year to time type
    type(time_type) :: m2d_offset    !< add to model time to get data time
end type forcing_timeseries_dataset

public :: time_interp_external, init_external_field, time_interp_external_init
public :: get_external_field_info
public :: horiz_interp_type, run_horiz_interp, build_horiz_interp_weights
public :: external_field
public :: forcing_timeseries_set_time_type_vars
public :: map_model_time_to_forcing_time

!> Read a field based on model time, and rotate to the model domain.
interface time_interp_external
  module procedure time_interp_external_0d
  module procedure time_interp_external_2d
  module procedure time_interp_external_3d
end interface time_interp_external


  interface
module subroutine time_interp_external_0d(field, time, data_in, verbose, scale)
  type(external_field), intent(in) :: field    !< Handle for time interpolated field
  type(time_type),   intent(in)    :: time     !< The target time for the data
  real,              intent(inout) :: data_in  !< The interpolated value in arbitrary units [A ~> a]
  logical, optional, intent(in)    :: verbose  !< If true, write verbose output for debugging
  real,    optional, intent(in)    :: scale    !< A scaling factor that new values of data_in are
                                               !! multiplied by before it is returned [A a-1 ~> 1]

  ! Store the input value in case the scaling factor is perfectly invertable.
end subroutine time_interp_external_0d
module subroutine time_interp_external_2d(field, time, data_in, interp, &
                                   verbose, horz_interp, mask_out, turns, scale)
  type(external_field), intent(in)    :: field    !< Handle for time interpolated field
  type(time_type),      intent(in)    :: time     !< The target time for the data
  real, dimension(:,:), intent(inout) :: data_in  !< The array in which to store the interpolated
                                                  !! values in arbitrary units [A ~> a]
  integer,    optional, intent(in)    :: interp   !< A flag indicating the temporal interpolation method
  logical,    optional, intent(in)    :: verbose  !< If true, write verbose output for debugging
  type(horiz_interp_type), &
              optional, intent(in)    :: horz_interp !< A structure to control horizontal interpolation
  logical, dimension(:,:), &
              optional, intent(out)   :: mask_out !< An array that is true where there is valid data
  integer,    optional, intent(in)    :: turns    !< Number of quarter turns to rotate the data
  real,       optional, intent(in)    :: scale    !< A scaling factor that new values of data_in are
                                                  !! multiplied by before it is returned [A a-1 ~> 1]


  ! TODO: Mask rotation requires logical array rotation support
end subroutine time_interp_external_2d
module subroutine time_interp_external_3d(field, time, data_in, interp, &
                                   verbose, horz_interp, mask_out, turns, scale)
  type(external_field), intent(in)      :: field    !< Handle for time interpolated field
  type(time_type),        intent(in)    :: time     !< The target time for the data
  real, dimension(:,:,:), intent(inout) :: data_in  !< The array in which to store the interpolated
                                                    !! values in arbitrary units [A ~> a]
  integer,      optional, intent(in)    :: interp   !< A flag indicating the temporal interpolation method
  logical,      optional, intent(in)    :: verbose  !< If true, write verbose output for debugging
  type(horiz_interp_type), &
                optional, intent(in)    :: horz_interp !< A structure to control horizontal interpolation
  logical, dimension(:,:,:), &
                optional, intent(out)   :: mask_out !< An array that is true where there is valid data
  integer,      optional, intent(in)    :: turns    !< Number of quarter turns to rotate the data
  real,         optional, intent(in)    :: scale    !< A scaling factor that new values of data_in are
                                                    !! multiplied by before it is returned [A a-1 ~> 1]


  ! TODO: Mask rotation requires logical array rotation support
end subroutine time_interp_external_3d
module subroutine forcing_timeseries_set_time_type_vars(data_start_year, data_end_year, data_ref_year, &
  model_ref_year, data_forcing_year, forcing_dataset)

  integer,                          intent(in)    :: data_start_year    !< first year of data to read
                                                                        !! (this is ignored for static forcing)
  integer,                          intent(in)    :: data_end_year      !< last year of data to read
                                                                        !! (this is ignored for static forcing)
  integer,                          intent(in)    :: data_ref_year      !< for time-varying forcing, align
                                                                        !! data_ref_year in file with
                                                                        !! model_ref_year in model
  integer,                          intent(in)    :: model_ref_year     !< for time-varying forcing, align
                                                                        !! data_ref_year in file with
                                                                        !! model_ref_year in model
  integer,                          intent(in)    :: data_forcing_year  !< for static forcing, read file at this
                                                                        !! date (this is ignored for time-varying
                                                                        !! forcing)
  type(forcing_timeseries_dataset), intent(inout) :: forcing_dataset    !< information about forcing file

end subroutine forcing_timeseries_set_time_type_vars
module function map_model_time_to_forcing_time(Time, forcing_dataset)

  type(time_type),                  intent(in)  :: Time             !< Model time
  type(forcing_timeseries_dataset), intent(in)  :: forcing_dataset  !< information about forcing file
  type(time_type) :: map_model_time_to_forcing_time                 !< time to read forcing file

end function map_model_time_to_forcing_time
module subroutine get_external_field_info(field, size, axes, missing)
  type(external_field), intent(in) :: field
    !< Handle for time interpolated external field returned from a previous
    !! call to init_external_field()
  integer, optional, intent(inout) :: size(4)
    !< Dimension sizes for the input data
  type(axis_info), optional, intent(inout) :: axes(4)
    !< Axis types for the input data
  real, optional, intent(inout) :: missing
    !< Missing value for the input data

    ! Axis as represented in the infra
    ! Axis name
    ! Axis points

    ! Axis index
    ! Axis size

end subroutine get_external_field_info
  end interface

end module MOM_interpolate
