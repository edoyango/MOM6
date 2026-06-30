submodule (MOM_io_infra) MOM_io_infra_s
  implicit none
contains
module procedure read_field_chksum
  integer(kind=int64), dimension(3) :: checksum_file
  checksum_file(:) = -1
  valid_chksum = mpp_attribute_exist(field, "checksum")
  if (valid_chksum) then
    call get_field_atts(field, checksum=checksum_file)
    chksum = checksum_file(1)
  else
    chksum = -1
  endif
end procedure read_field_chksum
module procedure MOM_file_exists
  MOM_file_exists = file_exist(filename, MOM_Domain%mpp_domain)

end procedure MOM_file_exists
module procedure FMS_file_exists
  FMS_file_exists = file_exist(filename)
end procedure FMS_file_exists
module procedure file_is_open
  file_is_open = (IO_handle%unit >= 0)
end procedure file_is_open
module procedure close_file_type
  call mpp_close(IO_handle%unit)
  if (allocated(IO_handle%filename)) deallocate(IO_handle%filename)
  IO_handle%open_to_read = .false. ; IO_handle%open_to_write = .false.
end procedure close_file_type
module procedure close_file_unit
  call mpp_close(unit)
end procedure close_file_unit
module procedure flush_file_type
  call mpp_flush(file%unit)
end procedure flush_file_type
module procedure flush_file_unit
  call mpp_flush(unit)
end procedure flush_file_unit
module procedure io_infra_init
  call mpp_io_init(maxunit=maxunits)
end procedure io_infra_init
module procedure io_infra_end
  call fms_io_exit()
end procedure io_infra_end
module procedure MOM_namelist_file
  unit = open_namelist_file(file)
end procedure MOM_namelist_file
module procedure check_namelist_error
  integer :: ierr
  ierr = check_nml_error(IOstat, nml_name)
end procedure check_namelist_error
module procedure write_version
  call write_version_number(version, tag, unit)
end procedure write_version
module procedure open_file_unit
  if (present(MOM_Domain)) then
    call mpp_open(unit, filename, action=action, form=form, threading=threading, fileset=fileset, &
                  nohdrs=nohdrs, domain=MOM_Domain%mpp_domain)
  else
    call mpp_open(unit, filename, action=action, form=form, threading=threading, fileset=fileset, &
                  nohdrs=nohdrs, domain=domain)
  endif
end procedure open_file_unit
module procedure open_file_type
  if (present(MOM_Domain)) then
    call mpp_open(IO_handle%unit, filename, action=action, form=NETCDF_FILE, threading=threading, &
                  fileset=fileset, domain=MOM_Domain%mpp_domain)
  else
    call mpp_open(IO_handle%unit, filename, action=action, form=NETCDF_FILE, threading=threading, &
                  fileset=fileset)
  endif
  IO_handle%filename = trim(filename)
  if (present(action)) then
    if (action == READONLY_FILE) then
      IO_handle%open_to_read = .true. ; IO_handle%open_to_write = .false.
    else
      IO_handle%open_to_read = .false. ; IO_handle%open_to_write = .true.
    endif
  else
    IO_handle%open_to_read = .false. ; IO_handle%open_to_write = .true.
  endif

end procedure open_file_type
module procedure open_ASCII_file
  call mpp_open(unit, file, action=action, form=ASCII_FILE, threading=threading, fileset=fileset, &
                  nohdrs=.true.)

end procedure open_ASCII_file
module procedure get_filename_suffix
  call get_filename_appendix(suffix)
end procedure get_filename_suffix
module procedure get_file_info
  integer :: ndims, nvars, natts, ntimes
  call mpp_get_info(IO_handle%unit, ndims, nvars, natts, ntimes )

  if (present(ndim)) ndim = ndims
  if (present(nvar)) nvar = nvars
  if (present(ntime)) ntime = ntimes

end procedure get_file_info
module procedure get_file_times
  integer :: ntimes
  if (allocated(time_values)) deallocate(time_values)
  call get_file_info(IO_handle, ntime=ntimes)
  if (present(ntime)) ntime = ntimes
  if (ntimes > 0) then
    allocate(time_values(ntimes))
    call mpp_get_times(IO_handle%unit, time_values)
  endif
end procedure get_file_times
module procedure get_file_fields
  call mpp_get_fields(IO_handle%unit, fields)
end procedure get_file_fields
module procedure get_field_atts
  call mpp_get_atts(field, name=name, units=units, longname=longname, checksum=checksum)
end procedure get_field_atts
module procedure field_exists
  if (present(MOM_domain)) then
    field_exists = field_exist(filename, field_name, domain=MOM_domain%mpp_domain, no_domain=no_domain)
  else
    field_exists = field_exist(filename, field_name, domain=domain, no_domain=no_domain)
  endif

end procedure field_exists
module procedure get_field_size
  call field_size(filename, fieldname, sizes, field_found=field_found, no_domain=no_domain)

end procedure get_field_size
module procedure get_axis_size
  axis_size = mpp_get_axis_length(axis)
end procedure get_axis_size
module procedure get_axis_data
  call mpp_get_atts(axis, name=axis_name)
  call mpp_get_axis_data(axis, axis_data)
end procedure get_axis_data
module procedure set_axis_data
  call MOM_error(FATAL, "set_axis_data in FMS1 is not yet implemented.")
end procedure set_axis_data
module procedure read_field_0d
  character(len=80)  :: varname             ! The name of a variable in the file
  type(fieldtype), allocatable :: fields(:) ! An array of types describing all the variables in the file
  logical :: use_fms_read_data, file_is_global
  integer :: n, unit, ndim, nvar, natt, ntime
  use_fms_read_data = .true. ; if (present(file_may_be_4d)) use_fms_read_data = .not.file_may_be_4d
  file_is_global = .true. ; if (present(global_file)) file_is_global = global_file

  if (.not.use_fms_read_data) then
    if (file_is_global) then
      call mpp_open(unit, trim(filename), form=NETCDF_FILE, action=READONLY_FILE, &
                    threading=MULTIPLE, fileset=SINGLE_FILE) !, domain=MOM_Domain%mpp_domain )
    else
      call mpp_open(unit, trim(filename), form=NETCDF_FILE, action=READONLY_FILE, &
                    threading=MULTIPLE, fileset=MULTIPLE, domain=MOM_Domain%mpp_domain )
    endif
    call mpp_get_info(unit, ndim, nvar, natt, ntime)
    allocate(fields(nvar))
    call mpp_get_fields(unit, fields(1:nvar))
    do n=1, nvar
      call mpp_get_atts(fields(n), name=varname)
      if (lowercase(trim(varname)) == lowercase(trim(fieldname))) then
        ! Maybe something should be done depending on the value of ntime.
        call mpp_read(unit, fields(n), data, timelevel)
        exit
      endif
    enddo

    deallocate(fields)
    call mpp_close(unit)
  elseif (present(MOM_Domain)) then
    call read_data(filename, fieldname, data, MOM_Domain%mpp_domain, timelevel=timelevel)
  else
    call read_data(filename, fieldname, data, timelevel=timelevel, no_domain=.true.)
  endif

  if (present(scale)) then ; if (scale /= 1.0) then
    data = scale*data
  endif ; endif
end procedure read_field_0d
module procedure read_field_1d
  character(len=80)  :: varname             ! The name of a variable in the file
  type(fieldtype), allocatable :: fields(:) ! An array of types describing all the variables in the file
  logical :: use_fms_read_data, file_is_global
  integer :: n, unit, ndim, nvar, natt, ntime
  use_fms_read_data = .true. ; if (present(file_may_be_4d)) use_fms_read_data = .not.file_may_be_4d
  file_is_global = .true. ; if (present(global_file)) file_is_global = global_file

  if (.not.use_fms_read_data) then
    if (file_is_global) then
      call mpp_open(unit, trim(filename), form=NETCDF_FILE, action=READONLY_FILE, &
                    threading=MULTIPLE, fileset=SINGLE_FILE) !, domain=MOM_Domain%mpp_domain )
    else
      call mpp_open(unit, trim(filename), form=NETCDF_FILE, action=READONLY_FILE, &
                    threading=MULTIPLE, fileset=MULTIPLE, domain=MOM_Domain%mpp_domain )
    endif
    call mpp_get_info(unit, ndim, nvar, natt, ntime)
    allocate(fields(nvar))
    call mpp_get_fields(unit, fields(1:nvar))
    do n=1, nvar
      call mpp_get_atts(fields(n), name=varname)
      if (lowercase(trim(varname)) == lowercase(trim(fieldname))) then
        call MOM_error(NOTE, "Reading 1-d variable "//trim(fieldname)//" from file "//trim(filename))
        ! Maybe something should be done depending on the value of ntime.
        call mpp_read(unit, fields(n), data, timelevel)
        exit
      endif
    enddo
    if ((n == nvar+1) .or. (nvar < 1)) call MOM_error(WARNING, &
      "read_field apparently did not find 1-d variable "//trim(fieldname)//" in file "//trim(filename))

    deallocate(fields)
    call mpp_close(unit)
  elseif (present(MOM_Domain)) then
    call read_data(filename, fieldname, data, MOM_Domain%mpp_domain, timelevel=timelevel)
  else
    call read_data(filename, fieldname, data, timelevel=timelevel, no_domain=.true.)
  endif

  if (present(scale)) then ; if (scale /= 1.0) then
    data(:) = scale*data(:)
  endif ; endif
end procedure read_field_1d
module procedure read_field_2d
  character(len=80)  :: varname             ! The name of a variable in the file
  type(fieldtype), allocatable :: fields(:) ! An array of types describing all the variables in the file
  logical :: use_fms_read_data, file_is_global
  integer :: n, unit, ndim, nvar, natt, ntime
  use_fms_read_data = .true. ; if (present(file_may_be_4d)) use_fms_read_data = .not.file_may_be_4d
  file_is_global = .true. ; if (present(global_file)) file_is_global = global_file

  if (use_fms_read_data) then
    call read_data(filename, fieldname, data, MOM_Domain%mpp_domain, &
                   timelevel=timelevel, position=position)
  else
    if (file_is_global) then
      call mpp_open(unit, trim(filename), form=NETCDF_FILE, action=READONLY_FILE, &
                    threading=MULTIPLE, fileset=SINGLE_FILE) !, domain=MOM_Domain%mpp_domain )
    else
      call mpp_open(unit, trim(filename), form=NETCDF_FILE, action=READONLY_FILE, &
                    threading=MULTIPLE, fileset=MULTIPLE, domain=MOM_Domain%mpp_domain )
    endif
    call mpp_get_info(unit, ndim, nvar, natt, ntime)
    allocate(fields(nvar))
    call mpp_get_fields(unit, fields(1:nvar))
    do n=1, nvar
      call mpp_get_atts(fields(n), name=varname)
      if (lowercase(trim(varname)) == lowercase(trim(fieldname))) then
        call MOM_error(NOTE, "Reading 2-d variable "//trim(fieldname)//" from file "//trim(filename))
        ! Maybe something should be done depending on the value of ntime.
        call mpp_read(unit, fields(n), MOM_Domain%mpp_domain, data, timelevel)
        exit
      endif
    enddo
    if ((n == nvar+1) .or. (nvar < 1)) call MOM_error(WARNING, &
      "read_field apparently did not find 2-d variable "//trim(fieldname)//" in file "//trim(filename))

    deallocate(fields)
    call mpp_close(unit)
  endif

  if (present(scale)) then ; if (scale /= 1.0) then
    call rescale_comp_data(MOM_Domain, data, scale)
  endif ; endif
end procedure read_field_2d
module procedure read_field_2d_region
  if (present(MOM_Domain)) then
    call read_data(filename, fieldname, data, start, nread, domain=MOM_Domain%mpp_domain, &
                   no_domain=no_domain)
  else
    call read_data(filename, fieldname, data, start, nread, no_domain=no_domain)
  endif

  if (present(scale)) then ; if (scale /= 1.0) then
    if (present(MOM_Domain)) then
      call rescale_comp_data(MOM_Domain, data, scale)
    else
      ! Dangerously rescale the whole array
      data(:,:) = scale*data(:,:)
    endif
  endif ; endif
end procedure read_field_2d_region
module procedure read_field_3d
  character(len=80)  :: varname             ! The name of a variable in the file
  type(fieldtype), allocatable :: fields(:) ! An array of types describing all the variables in the file
  logical :: use_fms_read_data, file_is_global
  integer :: n, unit, ndim, nvar, natt, ntime
  use_fms_read_data = .true. ; if (present(file_may_be_4d)) use_fms_read_data = .not.file_may_be_4d
  file_is_global = .true. ; if (present(global_file)) file_is_global = global_file

  if (use_fms_read_data) then
    call read_data(filename, fieldname, data, MOM_Domain%mpp_domain, &
                   timelevel=timelevel, position=position)
  else
    if (file_is_global) then
      call mpp_open(unit, trim(filename), form=NETCDF_FILE, action=READONLY_FILE, &
                    threading=MULTIPLE, fileset=SINGLE_FILE) !, domain=MOM_Domain%mpp_domain )
    else
      call mpp_open(unit, trim(filename), form=NETCDF_FILE, action=READONLY_FILE, &
                    threading=MULTIPLE, fileset=MULTIPLE, domain=MOM_Domain%mpp_domain )
    endif
    call mpp_get_info(unit, ndim, nvar, natt, ntime)
    allocate(fields(nvar))
    call mpp_get_fields(unit, fields(1:nvar))
    do n=1, nvar
      call mpp_get_atts(fields(n), name=varname)
      if (lowercase(trim(varname)) == lowercase(trim(fieldname))) then
        call MOM_error(NOTE, "Reading 3-d variable "//trim(fieldname)//" from file "//trim(filename))
        ! Maybe something should be done depending on the value of ntime.
        call mpp_read(unit, fields(n), MOM_Domain%mpp_domain, data, timelevel)
        exit
      endif
    enddo
    if ((n == nvar+1) .or. (nvar < 1)) call MOM_error(WARNING, &
      "read_field apparently did not find 3-d variable "//trim(fieldname)//" in file "//trim(filename))

    deallocate(fields)
    call mpp_close(unit)
  endif

  if (present(scale)) then ; if (scale /= 1.0) then
    call rescale_comp_data(MOM_Domain, data, scale)
  endif ; endif
end procedure read_field_3d
module procedure read_field_3d_region
  if (present(MOM_Domain)) then
    call read_data(filename, fieldname, data, start, nread, domain=MOM_Domain%mpp_domain, &
                   no_domain=no_domain)
  else
    call read_data(filename, fieldname, data, start, nread, no_domain=no_domain)
  endif

  if (present(scale)) then ; if (scale /= 1.0) then
    if (present(MOM_Domain)) then
      call rescale_comp_data(MOM_Domain, data, scale)
    else
      ! Dangerously rescale the whole array
      data(:,:,:) = scale*data(:,:,:)
    endif
  endif ; endif
end procedure read_field_3d_region
module procedure read_field_4d
  character(len=80)  :: varname             ! The name of a variable in the file
  type(fieldtype), allocatable :: fields(:) ! An array of types describing all the variables in the file
  logical :: file_is_global
  integer :: n, unit, ndim, nvar, natt, ntime
  file_is_global = .true. ; if (present(global_file)) file_is_global = global_file

  if (file_is_global) then
    call mpp_open(unit, trim(filename), form=NETCDF_FILE, action=READONLY_FILE, &
                  threading=MULTIPLE, fileset=SINGLE_FILE) !, domain=MOM_Domain%mpp_domain )
  else
    call mpp_open(unit, trim(filename), form=NETCDF_FILE, action=READONLY_FILE, &
                  threading=MULTIPLE, fileset=MULTIPLE, domain=MOM_Domain%mpp_domain )
  endif
  call mpp_get_info(unit, ndim, nvar, natt, ntime)
  allocate(fields(nvar))
  call mpp_get_fields(unit, fields(1:nvar))
  do n=1, nvar
    call mpp_get_atts(fields(n), name=varname)
    if (lowercase(trim(varname)) == lowercase(trim(fieldname))) then
        call MOM_error(NOTE, "Reading 4-d variable "//trim(fieldname)//" from file "//trim(filename))
      ! Maybe something should be done depending on the value of ntime.
      call mpp_read(unit, fields(n), MOM_Domain%mpp_domain, data, timelevel)
      exit
    endif
  enddo
  if ((n == nvar+1) .or. (nvar < 1)) call MOM_error(WARNING, &
    "read_field apparently did not find 4-d variable "//trim(fieldname)//" in file "//trim(filename))

  deallocate(fields)
  call mpp_close(unit)

  if (present(scale)) then ; if (scale /= 1.0) then
    call rescale_comp_data(MOM_Domain, data, scale)
  endif ; endif
end procedure read_field_4d
module procedure read_field_0d_int
  call read_data(filename, fieldname, data, timelevel=timelevel, no_domain=.true.)
end procedure read_field_0d_int
module procedure read_field_1d_int
  call read_data(filename, fieldname, data, timelevel=timelevel, no_domain=.true.)
end procedure read_field_1d_int
module procedure MOM_read_vector_2d
  integer :: u_pos, v_pos
  u_pos = EAST_FACE ; v_pos = NORTH_FACE
  if (present(stagger)) then
    if (stagger == CGRID_NE) then ; u_pos = EAST_FACE ; v_pos = NORTH_FACE
    elseif (stagger == BGRID_NE) then ; u_pos = CORNER ; v_pos = CORNER
    elseif (stagger == AGRID) then ; u_pos = CENTER ; v_pos = CENTER ; endif
  endif

  call read_data(filename, u_fieldname, u_data, MOM_Domain%mpp_domain, &
                 timelevel=timelevel, position=u_pos)
  call read_data(filename, v_fieldname, v_data, MOM_Domain%mpp_domain, &
                 timelevel=timelevel, position=v_pos)

  if (present(scale)) then ; if (scale /= 1.0) then
    call rescale_comp_data(MOM_Domain, u_data, scale)
    call rescale_comp_data(MOM_Domain, v_data, scale)
  endif ; endif

end procedure MOM_read_vector_2d
module procedure MOM_read_vector_3d
  integer :: u_pos, v_pos
  u_pos = EAST_FACE ; v_pos = NORTH_FACE
  if (present(stagger)) then
    if (stagger == CGRID_NE) then ; u_pos = EAST_FACE ; v_pos = NORTH_FACE
    elseif (stagger == BGRID_NE) then ; u_pos = CORNER ; v_pos = CORNER
    elseif (stagger == AGRID) then ; u_pos = CENTER ; v_pos = CENTER ; endif
  endif

  call read_data(filename, u_fieldname, u_data, MOM_Domain%mpp_domain, &
                 timelevel=timelevel, position=u_pos)
  call read_data(filename, v_fieldname, v_data, MOM_Domain%mpp_domain, &
                 timelevel=timelevel, position=v_pos)

  if (present(scale)) then ; if (scale /= 1.0) then
    call rescale_comp_data(MOM_Domain, u_data, scale)
    call rescale_comp_data(MOM_Domain, v_data, scale)
  endif ; endif

end procedure MOM_read_vector_3d
module procedure write_field_4d
  call mpp_write(IO_handle%unit, field_md, MOM_domain%mpp_domain, field, tstamp=tstamp, &
                 tile_count=tile_count, default_data=fill_value)
end procedure write_field_4d
module procedure write_field_3d
  call mpp_write(IO_handle%unit, field_md, MOM_domain%mpp_domain, field, tstamp=tstamp, &
                   tile_count=tile_count, default_data=fill_value)
end procedure write_field_3d
module procedure write_field_2d
  call mpp_write(IO_handle%unit, field_md, MOM_domain%mpp_domain, field, tstamp=tstamp, &
                   tile_count=tile_count, default_data=fill_value)
end procedure write_field_2d
module procedure write_field_1d
  call mpp_write(IO_handle%unit, field_md, field, tstamp=tstamp)
end procedure write_field_1d
module procedure write_field_0d
  call mpp_write(IO_handle%unit, field_md, field, tstamp=tstamp)
end procedure write_field_0d
module procedure MOM_write_axis
  call mpp_write(IO_handle%unit, axis)

end procedure MOM_write_axis
module procedure write_metadata_axis
  call mpp_write_meta(IO_handle%unit, axis, name, units, longname, cartesian=cartesian, sense=sense, &
                      domain=domain, data=data, calendar=calendar)
end procedure write_metadata_axis
module procedure write_metadata_field
  call mpp_write_meta(IO_handle%unit, field, axes, name, units, longname, &
                      pack=pack, standard_name=standard_name, checksum=checksum)
  ! unused opt. args: min=min, max=max, fill=fill, scale=scale, add=add, &

end procedure write_metadata_field
module procedure write_metadata_global
  call mpp_write_meta(IO_handle%unit, name, cval=attribute)
end procedure write_metadata_global
end submodule MOM_io_infra_s
