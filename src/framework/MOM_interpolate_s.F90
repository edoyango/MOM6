submodule (MOM_interpolate) MOM_interpolate_s
  implicit none
contains
module procedure time_interp_external_0d
  real :: data_in_pre_scale ! The input data before rescaling [a]
  real :: I_scale ! The inverse of scale [a A-1 ~> 1]
  data_in_pre_scale = data_in
  I_scale = 1.0
  if (present(scale)) then ; if ((scale /= 1.0) .and. (scale /= 0.0)) then
    ! Because time_interp_extern has the ability to only set some values, but no clear
    ! mechanism to determine which values have been set, the input data has to
    ! be unscaled so that it will have the right values when it is returned.
    I_scale = 1.0 / scale
    data_in = data_in * I_scale
  endif ; endif

  call time_interp_extern(field, time, data_in, verbose=verbose)

  if (present(scale)) then ; if (scale /= 1.0) then
    ! Rescale data that has been newly set and restore the scaling of unset data.
    if (data_in == I_scale * data_in_pre_scale) then
      data_in = data_in_pre_scale
    else
      data_in = scale * data_in
    endif
  endif ; endif

end procedure time_interp_external_0d
module procedure time_interp_external_2d
  real, allocatable :: data_in_pre_scale(:,:) ! The input data before rescaling [a]
  real, allocatable :: data_pre_rot(:,:)      ! The unscaled input data before rotation [a]
  real    :: I_scale ! The inverse of scale [a A-1 ~> 1]
  integer :: qturns ! The number of quarter turns to rotate the data
  integer :: i, j
  if (present(mask_out)) &
    call MOM_error(FATAL, "Rotation of masked output not yet support")

  if (present(scale)) then ; if ((scale /= 1.0) .and. (scale /= 0.0)) then
    ! Because time_interp_extern has the ability to only set some values, but no clear mechanism
    ! to determine which values have been set, the input data has to be unscaled so that it will
    ! have the right values when it is returned.  It may be a problem for some compiler settings
    ! if there are NaNs in data_in, but they will not spread.
    if (abs(fraction(scale)) /= 1.0) then
      ! This scaling factor may not be perfectly invertable, so store the input value
      allocate(data_in_pre_scale, source=data_in)
    endif
    I_scale = 1.0 / scale
    data_in(:,:) = I_scale * data_in(:,:)
  endif ; endif

  qturns = 0 ; if (present(turns)) qturns = modulo(turns, 4)

  if (qturns == 0) then
    call time_interp_extern(field, time, data_in, interp=interp, &
                            verbose=verbose, horz_interp=horz_interp)
  else
    call allocate_rotated_array(data_in, [1,1], -qturns, data_pre_rot)
    call time_interp_extern(field, time, data_pre_rot, interp=interp, &
                            verbose=verbose, horz_interp=horz_interp)
    call rotate_array(data_pre_rot, turns, data_in)
    deallocate(data_pre_rot)
  endif

  if (present(scale)) then ; if (scale /= 1.0) then
    ! Rescale data that has been newly set and restore the scaling of unset data.
    if ((abs(fraction(scale)) /= 1.0) .and. (scale /= 0.0)) then
      do j=LBOUND(data_in,2),UBOUND(data_in,2) ; do i=LBOUND(data_in,1),UBOUND(data_in,1)
        ! This handles the case where scale is not exactly invertable for data
        ! values that have not been modified by time_interp_extern.
        if (data_in(i,j) == I_scale * data_in_pre_scale(i,j)) then
          data_in(i,j) = data_in_pre_scale(i,j)
        else
          data_in(i,j) = scale * data_in(i,j)
        endif
      enddo ; enddo
    else
      data_in(:,:) = scale * data_in(:,:)
    endif
  endif ; endif

end procedure time_interp_external_2d
module procedure time_interp_external_3d
  real, allocatable :: data_in_pre_scale(:,:,:) ! The input data before rescaling [a]
  real, allocatable :: data_pre_rot(:,:,:)      ! The unscaled input data before rotation [a]
  real    :: I_scale ! The inverse of scale [a A-1 ~> 1]
  integer :: qturns  ! The number of quarter turns to rotate the data
  integer :: i, j, k
  if (present(mask_out)) &
    call MOM_error(FATAL, "Rotation of masked output not yet support")

  if (present(scale)) then ; if ((scale /= 1.0) .and. (scale /= 0.0)) then
    ! Because time_interp_extern has the ability to only set some values, but no clear mechanism
    ! to determine which values have been set, the input data has to be unscaled so that it will
    ! have the right values when it is returned.  It may be a problem for some compiler settings
    ! if there are NaNs in data_in, but they will not spread.
    if (abs(fraction(scale)) /= 1.0) then
      ! This scaling factor may not be perfectly invertable, so store the input value
      allocate(data_in_pre_scale, source=data_in)
    endif
    I_scale = 1.0 / scale
    data_in(:,:,:) = I_scale * data_in(:,:,:)
  endif ; endif

  qturns = 0 ; if (present(turns)) qturns = modulo(turns, 4)

  if (qturns == 0) then
    call time_interp_extern(field, time, data_in, interp=interp, &
                            verbose=verbose, horz_interp=horz_interp)
  else
    call allocate_rotated_array(data_in, [1,1,1], -qturns, data_pre_rot)
    call time_interp_extern(field, time, data_pre_rot, interp=interp, &
                            verbose=verbose, horz_interp=horz_interp)
    call rotate_array(data_pre_rot, turns, data_in)
    deallocate(data_pre_rot)
  endif

  if (present(scale)) then ; if (scale /= 1.0) then
    ! Rescale data that has been newly set and restore the scaling of unset data.
    if ((abs(fraction(scale)) /= 1.0) .and. (scale /= 0.0)) then
      do k=LBOUND(data_in,3),UBOUND(data_in,3)
        do j=LBOUND(data_in,2),UBOUND(data_in,2)
          do i=LBOUND(data_in,1),UBOUND(data_in,1)
            ! This handles the case where scale is not exactly invertable for data
            ! values that have not been modified by time_interp_extern.
            if (data_in(i,j,k) == I_scale * data_in_pre_scale(i,j,k)) then
              data_in(i,j,k) = data_in_pre_scale(i,j,k)
            else
              data_in(i,j,k) = scale * data_in(i,j,k)
            endif
          enddo
        enddo
      enddo
    else
      data_in(:,:,:) = scale * data_in(:,:,:)
    endif
  endif ; endif

end procedure time_interp_external_3d
module procedure forcing_timeseries_set_time_type_vars
  if (forcing_dataset%l_time_varying) then
    forcing_dataset%data_start = set_date(data_start_year, 1, 1)
    forcing_dataset%data_end = set_date(data_end_year, 1, 1)
    forcing_dataset%m2d_offset = set_date(data_ref_year - model_ref_year, 1, 1)
  else
    forcing_dataset%data_forcing = set_date(data_forcing_year, 1, 1)
  endif

end procedure forcing_timeseries_set_time_type_vars
module procedure map_model_time_to_forcing_time
  if (forcing_dataset%l_time_varying) then
    map_model_time_to_forcing_time = Time + forcing_dataset%m2d_offset
    ! If Time + offset is not between data_start and data_end, use whichever of those values is closer
    if (map_model_time_to_forcing_time < forcing_dataset%data_start) &
      map_model_time_to_forcing_time = forcing_dataset%data_start
    if (map_model_time_to_forcing_time > forcing_dataset%data_end) &
      map_model_time_to_forcing_time = forcing_dataset%data_end
  else
    map_model_time_to_forcing_time = forcing_dataset%data_forcing
  endif

end procedure map_model_time_to_forcing_time
module procedure get_external_field_info
  type(axistype) :: axes_infra(4)
  character(len=256) :: axis_name
  real, allocatable :: ax_data(:)
  integer :: n
  integer :: ax_size
  if (present(axes)) then
    call get_external_field_info_infra(field, size=size, axes=axes_infra, &
        missing=missing)
    ! TODO: Most of these methods were written to expect four dimensions.
    do n=1,4
      ! Convert axistype to axis_info
      ax_size = get_axis_size(axes_infra(n))
      allocate(ax_data(ax_size))
      call get_axis_data(axes_infra(n), axis_name, ax_data)
      call set_axis_info(axes(n), trim(axis_name), ax_data=ax_data)
      deallocate(ax_data)
    enddo
  else
    call get_external_field_info_infra(field, size=size, missing=missing)
  endif
end procedure get_external_field_info
end submodule MOM_interpolate_s
