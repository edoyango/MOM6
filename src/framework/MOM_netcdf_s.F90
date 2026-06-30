submodule (MOM_netcdf) MOM_netcdf_s
  implicit none
contains
module procedure open_netcdf_file
  integer :: io_mode
  integer :: cmode
  integer :: rc
  character(len=:), allocatable :: msg
  io_mode = WRITEONLY_FILE
  if (present(mode)) io_mode = mode

  ! Translate the MOM I/O config to the netCDF mode
  select case(io_mode)
    case (WRITEONLY_FILE)
      rc = nf90_create(filename, nf90_noclobber, handle%ncid)
      handle%define_mode = .true.
    case (OVERWRITE_FILE)
      rc = nf90_create(filename, nf90_clobber, handle%ncid)
      handle%define_mode = .true.
    case (APPEND_FILE)
      rc = nf90_open(filename, nf90_write, handle%ncid)
      handle%define_mode = .false.
    case (READONLY_FILE)
      rc = nf90_open(filename, nf90_nowrite, handle%ncid)
      handle%define_mode = .false.
    case default
      call MOM_error(FATAL, &
          'open_netcdf_file: File ' // filename // ': Unknown mode.')
  end select
  call check_netcdf_call(rc, 'open_netcdf_file', 'File ' // filename)

  handle%filename = filename

  ! FMS writes the filename as an attribute
  if (any(io_mode == [WRITEONLY_FILE, OVERWRITE_FILE])) &
    call write_netcdf_attribute(handle, 'filename', filename)
end procedure open_netcdf_file
module procedure close_netcdf_file
  integer :: rc
  rc = nf90_close(handle%ncid)
  call check_netcdf_call(rc, 'close_netcdf_file', &
      'File "' // handle%filename // '"')
end procedure close_netcdf_file
module procedure flush_netcdf_file
  integer :: rc
  rc = nf90_sync(handle%ncid)
  call check_netcdf_call(rc, 'flush_netcdf_file', &
    'File "' // handle%filename // '"')
end procedure flush_netcdf_file
module procedure enable_netcdf_write
  integer :: rc
  if (handle%define_mode) then
    rc = nf90_enddef(handle%ncid)
    call check_netcdf_call(rc, 'enable_netcdf_write', &
        'File "' // handle%filename // '"')
    handle%define_mode = .false.
  endif
end procedure enable_netcdf_write
module procedure register_netcdf_field
  integer :: rc
  integer :: i
  integer, allocatable :: dimids(:)
  integer :: xtype
  allocate(dimids(size(axes)))
  dimids(:) = [(axes(i)%dimid, i = 1, size(axes))]

  field%label = label

  ! Determine the corresponding netCDF data type
  ! TODO: Support a `pack`-like argument
  select case (kind(1.0))
    case (real32)
      xtype = NF90_FLOAT
    case (real64)
      xtype = NF90_DOUBLE
    case default
      call MOM_error(FATAL, "register_netcdf_field: Unknown kind(real).")
  end select

  ! Register the field variable
  rc = nf90_def_var(handle%ncid, label, xtype, dimids, field%varid)
  call check_netcdf_call(rc, 'register_netcdf_field', &
      'File "' // handle%filename // '", Field "' // label // '"')

  ! Assign attributes

  rc = nf90_put_att(handle%ncid, field%varid, 'long_name', longname)
  call check_netcdf_call(rc, 'register_netcdf_field', &
    'Attribute "long_name" of variable "' // label // '" in file "' &
    // handle%filename // '"')

  rc = nf90_put_att(handle%ncid, field%varid, 'units', units)
  call check_netcdf_call(rc, 'register_netcdf_field', &
    'Attribute "units" of variable "' // label // '" in file "' &
    // handle%filename // '"')
end procedure register_netcdf_field
module procedure register_netcdf_axis
  integer :: xtype
  integer :: rc
  logical :: unlimited
  integer :: axis_size
  integer :: axis_sense
  character(len=:), allocatable :: sense_attr
  unlimited = .false.
  if (present(cartesian)) then
    if (cartesian == 'T') unlimited = .true.
  endif

  ! Either the axis is explicitly set with data or is declared as unlimited
  if (present(points) .eqv. unlimited) then
    call MOM_error(FATAL, &
        "Axis must either have explicit points or be a time axis ('T').")
  endif

  axis%label = label

  if (present(points)) then
    axis_size = size(points)
    allocate(axis%points(axis_size))
    axis%points(:) = points(:)
  else
    axis_size = NF90_UNLIMITED
  endif

  rc = nf90_def_dim(handle%ncid, label, axis_size, axis%dimid)
  call check_netcdf_call(rc, 'register_netcdf_axis', &
      'Dimension "' // label // '" in file "' // handle%filename // '"')

  ! Determine the corresponding netCDF data type
  ! TODO: Support a `pack`-like argument
  select case (kind(1.0))
    case (real32)
      xtype = NF90_FLOAT
    case (real64)
      xtype = NF90_DOUBLE
    case default
      call MOM_error(FATAL, "register_netcdf_axis: Unknown kind(real).")
  end select

  ! Create a variable corresponding to the axis
  rc = nf90_def_var(handle%ncid, label, xtype, axis%dimid, axis%varid)
  call check_netcdf_call(rc, 'register_netcdf_axis', &
      'Variable ' // label // ' in file ' // handle%filename)

  ! Define the time axis, if available
  if (unlimited) then
    handle%time_id = axis%varid
    handle%time_level = 0
    handle%time = NULLTIME
  endif

  ! Assign attributes if present
  if (present(longname)) then
    rc = nf90_put_att(handle%ncid, axis%varid, 'long_name', longname)
    call check_netcdf_call(rc, 'register_netcdf_axis', &
      'Attribute ''long_name'' of variable ' // label // ' in file ' &
      // handle%filename)
  endif

  if (present(units)) then
    rc = nf90_put_att(handle%ncid, axis%varid, 'units', units)
    call check_netcdf_call(rc, 'register_netcdf_axis', &
      'Attribute ''units'' of variable ' // label // ' in file ' &
      // handle%filename)
  endif

  if (present(cartesian)) then
    rc = nf90_put_att(handle%ncid, axis%varid, 'cartesian_axis', cartesian)
    call check_netcdf_call(rc, 'register_netcdf_axis', &
      'Attribute ''cartesian_axis'' of variable ' // label // ' in file ' &
      // handle%filename)
  endif

  axis_sense = 0
  if (present(sense)) axis_sense = sense

  if (axis_sense /= 0) then
    select case (axis_sense)
      case (1)
        sense_attr = 'up'
      case (-1)
        sense_attr = 'down'
      case default
        call MOM_error(FATAL, 'register_netcdf_axis: sense must be either ' &
          // '0, 1, or -1.')
    end select
    rc = nf90_put_att(handle%ncid, axis%varid, 'positive', sense_attr)
    call check_netcdf_call(rc, 'register_netcdf_axis', &
      'Attribute "positive" of variable "' // label // '" in file "' &
      // handle%filename // '"')
  endif
end procedure register_netcdf_axis
module procedure write_netcdf_field_4d
  integer :: rc
  integer :: start(5)
  if (handle%define_mode) &
    call enable_netcdf_write(handle)

  if (present(time)) then
    call update_netcdf_timestep(handle, time)
    start(:4) = 1
    start(5) = handle%time_level
    rc = nf90_put_var(handle%ncid, field%varid, values, start)
  else
    rc = nf90_put_var(handle%ncid, field%varid, values)
  endif
  call check_netcdf_call(rc, 'write_netcdf_file', &
      'File "' // handle%filename // '", Field "' // field%label // '"')
end procedure write_netcdf_field_4d
module procedure write_netcdf_field_3d
  integer :: rc
  integer :: start(4)
  if (handle%define_mode) &
    call enable_netcdf_write(handle)

  if (present(time)) then
    call update_netcdf_timestep(handle, time)
    start(:3) = 1
    start(4) = handle%time_level
    rc = nf90_put_var(handle%ncid, field%varid, values, start)
  else
    rc = nf90_put_var(handle%ncid, field%varid, values)
  endif
  call check_netcdf_call(rc, 'write_netcdf_file', &
      'File "' // handle%filename // '", Field "' // field%label // '"')
end procedure write_netcdf_field_3d
module procedure write_netcdf_field_2d
  integer :: rc
  integer :: start(3)
  if (handle%define_mode) &
    call enable_netcdf_write(handle)

  if (present(time)) then
    call update_netcdf_timestep(handle, time)
    start(:2) = 1
    start(3) = handle%time_level
    rc = nf90_put_var(handle%ncid, field%varid, values, start)
  else
    rc = nf90_put_var(handle%ncid, field%varid, values)
  endif
  call check_netcdf_call(rc, 'write_netcdf_file', &
      'File "' // handle%filename // '", Field "' // field%label // '"')
end procedure write_netcdf_field_2d
module procedure write_netcdf_field_1d
  integer :: rc
  integer :: start(2)
  if (handle%define_mode) &
    call enable_netcdf_write(handle)

  if (present(time)) then
    call update_netcdf_timestep(handle, time)
    start(1) = 1
    start(2) = handle%time_level
    rc = nf90_put_var(handle%ncid, field%varid, values, start)
  else
    rc = nf90_put_var(handle%ncid, field%varid, values)
  endif
  call check_netcdf_call(rc, 'write_netcdf_file', &
      'File "' // handle%filename // '", Field "' // field%label // '"')
end procedure write_netcdf_field_1d
module procedure write_netcdf_field_0d
  integer :: rc
  integer :: start(1)
  if (handle%define_mode) &
    call enable_netcdf_write(handle)

  if (present(time)) then
    call update_netcdf_timestep(handle, time)
    start(1) = handle%time_level
    rc = nf90_put_var(handle%ncid, field%varid, scalar, start)
  else
    rc = nf90_put_var(handle%ncid, field%varid, scalar)
  endif
  call check_netcdf_call(rc, 'write_netcdf_file', &
      'File "' // handle%filename // '", Field "' // field%label // '"')
end procedure write_netcdf_field_0d
module procedure write_netcdf_axis
  integer :: rc
  if (handle%define_mode) &
    call enable_netcdf_write(handle)

  rc = nf90_put_var(handle%ncid, axis%varid, axis%points)
  call check_netcdf_call(rc, 'write_netcdf_axis', &
      'File "' // handle%filename // '", Axis "' // axis%label // '"')
end procedure write_netcdf_axis
module procedure write_netcdf_attribute
  integer :: rc
  rc = nf90_put_att(handle%ncid, NF90_GLOBAL, label, attribute)
  call check_netcdf_call(rc, 'write_netcdf_attribute', &
      'File "' // handle%filename // '", Attribute "' // label // '"')
end procedure write_netcdf_attribute
module procedure get_netcdf_size
  integer :: rc
  integer :: unlimited_dimid
  rc = nf90_inquire(handle%ncid, &
      nDimensions=ndims, &
      nVariables=nvars, &
      unlimitedDimId=unlimited_dimid &
  )
  call check_netcdf_call(rc, 'get_netcdf_size', &
      'File "' // handle%filename // '"')

  rc = nf90_inquire_dimension(handle%ncid, unlimited_dimid, len=nsteps)
  call check_netcdf_call(rc, 'get_netcdf_size', &
      'File "' // handle%filename // '"')
end procedure get_netcdf_size
module procedure get_netcdf_fields
  integer :: ndims
  integer :: nvars
  type(netcdf_field), allocatable :: vars(:)
  integer :: nfields
  integer, allocatable :: dimids(:)
  integer, allocatable :: varids(:)
  integer :: unlim_dimid
  integer :: unlim_index
  character(len=NF90_MAX_NAME) :: label
  integer :: len
  integer :: rc
  integer :: grp_ndims, grp_nvars
  logical :: is_axis
  integer :: i, j, n
  integer, save :: no_parent_groups = 0
  rc = nf90_inquire(handle%ncid, &
      nDimensions=ndims, &
      nVariables=nvars, &
      unlimitedDimId=unlim_dimid &
  )
  call check_netcdf_call(rc, 'get_netcdf_fields', &
      'File "' // handle%filename // '"')

  allocate(dimids(ndims))
  rc = nf90_inq_dimids(handle%ncid, grp_ndims, dimids, no_parent_groups)
  call check_netcdf_call(rc, 'get_netcdf_fields', &
      'File "' // handle%filename // '"')

  allocate(varids(nvars))
  rc = nf90_inq_varids(handle%ncid, grp_nvars, varids)
  call check_netcdf_call(rc, 'get_netcdf_fields', &
      'File "' // trim(handle%filename) // '"')

  ! Initialize unlim_index with an unreachable value (outside [1,ndims])
  unlim_index = -1

  allocate(axes(ndims))
  do i = 1, ndims
    rc = nf90_inquire_dimension(handle%ncid, dimids(i), name=label, len=len)
    call check_netcdf_call(rc, 'get_netcdf_fields', &
        'File "' // trim(handle%filename) // '"')

    ! Check for the unlimited axis
    if (dimids(i) == unlim_dimid) unlim_index = i

    axes(i)%dimid = dimids(i)
    axes(i)%label = trim(label)
    allocate(axes(i)%points(len))
  enddo

  ! We cannot know if every axis also has a variable representation, so we
  ! over-allocate vars(:) and fill as fields are identified.
  allocate(vars(nvars))

  nfields = 0
  do i = 1, nvars
    rc = nf90_inquire_variable(handle%ncid, varids(i), name=label)
    call check_netcdf_call(rc, 'get_netcdf_fields', &
        'File "' // trim(handle%filename) // '"')

    ! Check if variable is an axis
    is_axis = .false.
    do j = 1, ndims
      if (label == axes(j)%label) then
        rc = nf90_get_var(handle%ncid, varids(i), axes(j)%points)
        call check_netcdf_call(rc, 'get_netcdf_fields', &
            'File "' // trim(handle%filename) // '"')
        axes(j)%varid = varids(i)

        if (j == unlim_index) then
          handle%time_id = varids(i)
          handle%time_level = size(axes(j)%points)
          handle%time = NULLTIME
        endif

        is_axis = .true.
        exit
      endif
    enddo
    if (is_axis) cycle

    nfields = nfields + 1
    vars(nfields)%label = trim(label)
    vars(nfields)%varid = varids(i)
  enddo

  allocate(fields(nfields))
  fields(:) = vars(:nfields)
end procedure get_netcdf_fields
module procedure get_netcdf_filename
  get_netcdf_filename = handle%filename

end procedure get_netcdf_filename
module procedure read_netcdf_field
  integer :: rc
  integer :: istart(2)
  integer :: icount(2)
  if (present(bounds)) then
    istart(:) = bounds(1,:)
    icount(:) = bounds(2,:) - bounds(1,:) + 1
    rc = nf90_get_var(handle%ncid, field%varid, values, start=istart, count=icount)
  else
    rc = nf90_get_var(handle%ncid, field%varid, values)
  endif
  call check_netcdf_call(rc, 'read_netcdf_field', &
      'File "' // trim(handle%filename) // '", Field "' // trim(field%label) // '"')
end procedure read_netcdf_field
module procedure update_netcdf_timestep
  integer :: start(1)
  integer :: rc
  if (time > handle%time + epsilon(time)) then
    handle%time = time
    handle%time_level = handle%time_level + 1

    ! Write new value to time axis
    start = [handle%time_level]
    rc = nf90_put_var(handle%ncid, handle%time_id, time, start=start)
    call check_netcdf_call(rc, 'update_netcdf_timestep', &
        'File "' // handle%filename // '"')
  endif
end procedure update_netcdf_timestep
module procedure check_netcdf_call
  character(len=:), allocatable :: errmsg
  if (ncerr /= NF90_NOERR) then
    errmsg = trim(header) // ": " // trim(message) // new_line('/') &
      // trim(nf90_strerror(ncerr))
    call MOM_error(FATAL, errmsg)
  endif
end procedure check_netcdf_call
end submodule MOM_netcdf_s
