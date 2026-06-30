submodule (MOM_io_infra) MOM_io_infra_s
  implicit none
contains
module procedure read_field_chksum
  chksum = -1
  valid_chksum = field%valid_chksum
  if (valid_chksum) chksum = field%chksum_read

end procedure read_field_chksum
module procedure MOM_file_exists
  type(FmsNetcdfDomainFile_t) :: fileobj
  MOM_file_exists = fms2_open_file(fileobj, filename, "read", MOM_Domain%mpp_domain)
  if (MOM_file_exists) call fms2_close_file(fileobj)
end procedure MOM_file_exists
module procedure FMS_file_exists
  FMS_file_exists = fms2_file_exist(filename)
end procedure FMS_file_exists
module procedure file_is_open
  file_is_open = ((IO_handle%unit >= 0) .or. associated(IO_handle%fileobj))
end procedure file_is_open
module procedure close_file_type
  if (associated(IO_handle%fileobj)) then
    call fms2_close_file(IO_handle%fileobj)
    deallocate(IO_handle%fileobj)
  endif
  if (allocated(IO_handle%filename)) deallocate(IO_handle%filename)
  IO_handle%open_to_read = .false. ; IO_handle%open_to_write = .false.
  IO_handle%num_times = 0 ; IO_handle%file_time = 0.0
end procedure close_file_type
module procedure close_file_unit
  logical :: unit_is_open
  inquire(iounit, opened=unit_is_open)
  if (unit_is_open) close(iounit)
end procedure close_file_unit
module procedure flush_file
  if (associated(IO_handle%fileobj)) then
    call fms2_flush_file(IO_handle%fileobj)
  endif
end procedure flush_file
module procedure io_infra_init
end procedure io_infra_init
module procedure io_infra_end
end procedure io_infra_end
module procedure MOM_namelist_file
  character(len=:), allocatable :: nmlpath
  character(len=:), allocatable :: nmlpath_pe
  if (present(filepath)) then
    nmlpath = trim(filepath)
  else
    ! FMS1 first checks for a namelist unique to the PE list, `input_{}.nml`.
    ! If not found, it defaults to `input.nml`.
    nmlpath_pe = 'input_' // trim(mpp_get_current_pelist_name()) // '.nml'
    if (file_exists(nmlpath_pe)) then
      nmlpath = nmlpath_pe
    else
      nmlpath = 'input.nml'
    endif
  endif
  call open_ASCII_file(iounit, nmlpath, action=READONLY_FILE)
end procedure MOM_namelist_file
module procedure check_namelist_error
  integer :: ierr
  ierr = check_nml_error(IOstat, nml_name)
end procedure check_namelist_error
module procedure write_version
  call write_version_number(version, tag, unit)
end procedure write_version
module procedure open_file
  type(FmsNetcdfDomainFile_t) :: fileobj_read ! A handle to a domain-decomposed file for obtaining information
  logical :: success         ! If true, the file was opened successfully
  integer :: file_mode       ! An integer that encodes whether the file is to be opened for
  character(len=40)  :: mode ! A character string that encodes whether the file is to be opened for
  character(len=:), allocatable :: filename_tmp  ! A copy of filename with .nc appended if necessary.
  character(len=256) :: dim_unlim_name ! name of the unlimited dimension in the file
  integer :: index_nc
  if (IO_handle%open_to_write) then
    call MOM_error(WARNING, "open_file called for file "//trim(filename)//&
        " with an IO_handle that is already open to to write.")
    return
  endif
  if (IO_handle%open_to_read) then
    call MOM_error(FATAL, "open_file called for file "//trim(filename)//&
        " with an IO_handle that is already open to to read.")
  endif

  file_mode = WRITEONLY_FILE ; if (present(action)) file_mode = action

  ! Domains are currently required to use FMS I/O.
  ! NOTE: We restrict FMS2 IO usage to domain-based files due to issues with
  ! string-based attributes in certain compilers.
  ! But we may relax this requirement in the future.
  if (.not. present(MOM_Domain)) &
    call MOM_error(FATAL, 'open_file: FMS I/O requires a domain input.')

  if (.not.associated(IO_handle%fileobj)) allocate (IO_handle%fileobj)

  ! The FMS1 interface automatically appends .nc if necessary, but FMS2 interface does not.
  index_nc = index(trim(filename), ".nc")
  if (index_nc > 0) then
    filename_tmp = trim(filename)
  else
    filename_tmp = trim(filename)//".nc"
    if (is_root_PE()) call MOM_error(WARNING, "Open_file is appending .nc to the filename "//trim(filename))
  endif

  if (file_mode == WRITEONLY_FILE) then ; mode = "write"
  elseif (file_mode == APPEND_FILE) then ; mode = "append"
  elseif (file_mode == OVERWRITE_FILE) then ; mode = "overwrite"
  elseif (file_mode == READONLY_FILE) then ; mode = "read"
  else
    call MOM_error(FATAL, "open_file called with unrecognized action.")
  endif

  IO_handle%num_times = 0
  IO_handle%file_time = 0.0
  if ((file_mode == APPEND_FILE) .and. file_exists(filename_tmp, MOM_Domain)) then
    ! Determine the latest file time and number of records so far.
    success = fms2_open_file(fileObj_read, trim(filename_tmp), "read", MOM_domain%mpp_domain)
    dim_unlim_name = find_unlimited_dimension_name(fileObj_read)
    if (len_trim(dim_unlim_name) > 0) &
      call get_dimension_size(fileObj_read, trim(dim_unlim_name), IO_handle%num_times)
    if (IO_handle%num_times > 0) &
      call fms2_read_data(fileObj_read, trim(dim_unlim_name), IO_handle%file_time, &
                          unlim_dim_level=IO_handle%num_times)
    call fms2_close_file(fileObj_read)
  endif

  success = fms2_open_file(IO_handle%fileobj, trim(filename_tmp), trim(mode), MOM_domain%mpp_domain)
  if (.not.success) call MOM_error(FATAL, "Unable to open file "//trim(filename_tmp))
  IO_handle%filename = trim(filename)

  if (file_mode == READONLY_FILE) then
    IO_handle%open_to_read = .true. ; IO_handle%open_to_write = .false.
  else
    IO_handle%open_to_read = .false. ; IO_handle%open_to_write = .true.
  endif

end procedure open_file
module procedure open_ASCII_file
  integer :: action_flag
  integer :: threading_flag
  integer :: fileset_flag
  logical :: exists
  logical :: is_open
  character(len=6) :: action_arg, position_arg
  character(len=:), allocatable :: filename
  action_flag = WRITEONLY_FILE
  if (present(action)) action_flag = action

  action_arg = 'write'
  if (action_flag == READONLY_FILE) action_arg = 'read'

  position_arg = 'rewind'
  if (action_flag == APPEND_FILE) position_arg = 'append'

  ! Threading configuration

  threading_flag = SINGLE_FILE
  if (present(threading)) threading_flag = threading

  fileset_flag = MULTIPLE
  if (present(fileset)) fileset_flag = fileset

  ! Force fileset to be consistent with threading (as in FMS1)
  if (threading_flag == SINGLE_FILE) fileset_flag = SINGLE_FILE

  ! Construct the distributed filename, if needed
  filename = file
  if (fileset_flag == MULTIPLE) then
    if (mpp_npes() > 10000) then
      write(filename, '(a,".",i6.6)') trim(filename), mpp_pe() - mpp_root_pe()
    else
      write(filename, '(a,".",i4.4)') trim(filename), mpp_pe() - mpp_root_pe()
    endif
  endif

  inquire(file=filename, exist=exists)
  if (exists .and. action_flag == WRITEONLY_FILE) &
    call MOM_error(WARNING, 'open_ASCII_file: File ' // trim(filename) // &
                            ' opened WRITEONLY already exists!')

  open(newunit=unit, file=filename, action=trim(action_arg), &
       position=trim(position_arg))

  ! This checks if open() failed but did not raise a runtime error.
  inquire(unit, opened=is_open)
  if (.not. is_open) &
    call MOM_error(FATAL, &
        'open_ASCII_file: File "' // trim(filename) // '" failed to open.')

  ! NOTE: There are two possible mpp_write_meta functions in FMS1:
  ! - call mpp_write_meta( unit, 'filename', cval=mpp_file(unit)%name)
  ! - call mpp_write_meta( unit, 'NumFilesInSet', ival=nfiles)
  ! I'm not convinced we actually want these, but note them here in case.
end procedure open_ASCII_file
module procedure get_filename_suffix
  call get_filename_appendix(suffix)
end procedure get_filename_suffix
module procedure get_file_info
  character(len=256) :: dim_unlim_name ! name of the unlimited dimension in the file
  integer :: ndims, nvars, natts, ntimes
  if (present(ndim)) ndim = get_num_dimensions(IO_handle%fileobj)
  if (present(nvar)) nvar = get_num_variables(IO_handle%fileobj)
  if (present(ntime)) then
    ntime = 0
    dim_unlim_name = find_unlimited_dimension_name(IO_handle%fileobj)
    if (len_trim(dim_unlim_name) > 0) &
      call get_dimension_size(IO_handle%fileobj, trim(dim_unlim_name), ntime)
  endif
end procedure get_file_info
module procedure get_file_times
  character(len=256) :: dim_unlim_name ! name of the unlimited dimension in the file
  integer :: ntimes  ! The number of time levels in the file
  if (allocated(time_values)) deallocate(time_values)
  call get_file_info(IO_handle, ntime=ntimes)
  if (present(ntime)) ntime = ntimes
  if (ntimes > 0) then
    allocate(time_values(ntimes))
    dim_unlim_name = find_unlimited_dimension_name(IO_handle%fileobj)
    if (len_trim(dim_unlim_name) > 0) &
      call fms2_read_data(IO_handle%fileobj, trim(dim_unlim_name), time_values)
  endif
end procedure get_file_times
module procedure get_file_fields
  character(len=256),  dimension(size(fields)) :: var_names ! The names of all variables
  character(len=256)  :: units    ! The units of a variable as recorded in the file
  character(len=2048) :: longname ! The long-name of a variable as recorded in the file
  character(len=64)   :: checksum_char ! The hexadecimal checksum read from the file
  integer(kind=int64), dimension(3) :: checksum_file ! The checksums for a variable in the file
  integer :: nvar  ! The number of variables in the file
  integer :: i
  nvar = size(fields)
  ! Local variables
  call get_variable_names(IO_handle%fileobj, var_names)
  do i=1,nvar
    fields(i)%name = trim(var_names(i))
    longname = ""
    if (variable_att_exists(IO_handle%fileobj, var_names(i), "long_name")) &
      call get_variable_attribute(IO_handle%fileobj, var_names(i), "long_name", longname)
    fields(i)%longname = trim(longname)
    units = ""
    if (variable_att_exists(IO_handle%fileobj, var_names(i), "units")) &
      call get_variable_attribute(IO_handle%fileobj, var_names(i), "units", units)
    fields(i)%units = trim(units)

    fields(i)%valid_chksum = variable_att_exists(IO_handle%fileobj, var_names(i), "checksum")
    if (fields(i)%valid_chksum) then
      call get_variable_attribute(IO_handle%fileobj, var_names(i), 'checksum', checksum_char)
      ! If there are problems, there might need to be code added to handle commas.
      read (checksum_char(1:16), '(Z16)') fields(i)%chksum_read
    endif
  enddo
end procedure get_file_fields
module procedure get_field_atts
  if (present(name)) name = trim(field%name)
  if (present(units)) units = trim(field%units)
  if (present(longname)) longname = trim(field%longname)
  if (present(checksum)) checksum = field%chksum_read

end procedure get_field_atts
module procedure field_exists
  type(FmsNetcdfDomainFile_t) :: fileObj_dd ! A handle to a domain-decomposed file for obtaining information
  type(FmsNetcdfFile_t) :: fileObj_simple   ! A handle to a non-domain-decomposed file for obtaining information
  logical :: success         ! If true, the file was opened successfully
  logical :: domainless      ! If true, this file does not use a domain-decomposed file.
  domainless = .not.(present(MOM_domain) .or. present(domain))
  if (present(no_domain)) then
    if (domainless .and. .not.no_domain) call MOM_error(FATAL, &
        "field_exists: When no_domain is present and false, a domain must be supplied in query about "//&
        trim(field_name)//" in file "//trim(filename))
    domainless = no_domain
  endif

  field_exists = .false.
  if (file_exists(filename)) then
    if (domainless) then
      success = fms2_open_file(fileObj_simple, trim(filename), "read")
      if (success) then
        field_exists = variable_exists(fileObj_simple, field_name)
        call fms2_close_file(fileObj_simple)
      endif
    else
      if (present(MOM_domain)) then
        success = fms2_open_file(fileObj_dd, trim(filename), "read", MOM_domain%mpp_domain)
      else
        success = fms2_open_file(fileObj_dd, trim(filename), "read", domain)
      endif
      if (success) then
        field_exists = variable_exists(fileobj_dd, field_name)
        call fms2_close_file(fileObj_dd)
      endif
    endif
  endif
end procedure field_exists
module procedure get_field_size
  type(FmsNetcdfFile_t) :: fileobj_read ! A handle to a non-domain-decomposed file for obtaining information
  logical :: success         ! If true, the file was opened successfully
  logical :: field_exists    ! True if filename exists and field_name is in filename
  integer :: i, ndims
  character(len=512), allocatable :: dimnames(:)  ! Field dimension names
  logical, allocatable :: is_x(:), is_y(:), is_t(:)     ! True if index matches axis type
  integer :: size_indices(4)        ! Mapping of size index to FMS1 convention
  integer :: idx, swap
  field_exists = .false.
  if (file_exists(filename)) then
    success = fms2_open_file(fileObj_read, trim(filename), "read")
    if (success) then
      field_exists = variable_exists(fileobj_read, fieldname)
      if (field_exists) then
        ndims = get_variable_num_dimensions(fileobj_read, fieldname)
        if (ndims > size(sizes)) call MOM_error(FATAL, &
          "get_field_size called with too few sizes for "//trim(fieldname)//" in "//trim(filename))
        call get_variable_size(fileobj_read, fieldname, sizes(1:ndims))

        do i=ndims+1,size(sizes) ; sizes(i) = 0 ; enddo

        ! If sizes exceeds ndims, then we fallback to the FMS1 convention
        ! where sizes has at least 4 dimension, and try to position values.
        if (size(sizes) > ndims)  then
          ! Assume FMS1 positioning rules: (nx, ny, nz, nt, ...)
          if (size(sizes) < 4) &
            call MOM_error(FATAL, "If sizes(:) exceeds field dimensions, "&
                &"then its length must be at least 4.")

          ! Fall back to the FMS1 default values of 1 (from mpp field%size)
          sizes(ndims+1:) = 1

          ! Gather the field dimension names
          allocate(dimnames(ndims))
          dimnames(:) = ""
          call get_variable_dimension_names(fileObj_read, trim(fieldname), &
                                            dimnames)

          ! Test the dimensions against standard (x,y,t) names and attributes
          allocate(is_x(ndims), is_y(ndims), is_t(ndims))
          is_x(:) = .false.
          is_y(:) = .false.
          is_t(:) = .false.
          call categorize_axes(fileObj_read, filename, ndims, dimnames, &
                               is_x, is_y, is_t)

          ! Currently no z-test is supported, so disable assignment with 0
          size_indices = [ &
              find_index(is_x), &
              find_index(is_y), &
              0, &
              find_index(is_t) &
          ]

          do i = 1, size(size_indices)
            idx = size_indices(i)
            if (idx > 0) then
              swap = sizes(i)
              sizes(i) = sizes(idx)
              sizes(idx) = swap
            endif
          enddo

          deallocate(is_x, is_y, is_t)
          deallocate(dimnames)
        endif
      endif
    endif
  endif
  if (present(field_found)) field_found = field_exists
end procedure get_field_size
module procedure find_index
  integer :: i
  loc = 0
  do i = 1, size(vec)
    if (vec(i)) then
      loc = i
      exit
    endif
  enddo
end procedure find_index
module procedure get_axis_size
  axis_size = size(axis%ax_data)
end procedure get_axis_size
module procedure get_axis_data
  integer :: i
  if (allocated(axis%ax_data)) then
    if (size(axis%ax_data) > size(axis_data)) &
      call MOM_error(FATAL, "get_axis_data called with too small of an " &
          // "output data array for " // trim(axis%name) // ".")
    do i=1,size(axis%ax_data)
      axis_data(i) = axis%ax_data(i)
    enddo
  endif

  axis_name = axis%name
end procedure get_axis_data
module procedure set_axis_data
  axis%name = axis_name

  if (allocated(axis%ax_data)) deallocate(axis%ax_data)
  allocate(axis%ax_data(size(axis_data)))

  axis%ax_data(:) = axis_data(:)

  ! NOTE: We do not yet consider domain-decomposed axes.
  axis%domain_decomposed = .false.
end procedure set_axis_data
module procedure read_field_0d
  type(FmsNetcdfFile_t)       :: fileObj ! A handle to a non-domain-decomposed file
  type(FmsNetcdfDomainFile_t) :: fileobj_DD ! A handle to a domain-decomposed file object
  character(len=96) :: var_to_read ! Name of variable to read from the netcdf file
  logical :: has_time_dim          ! True if the variable has an unlimited time axis.
  logical :: success               ! True if the file was successfully opened
  if (present(MOM_Domain)) then
    ! Open the FMS2 file-set.
    success = fms2_open_file(fileobj_DD, filename, "read", MOM_domain%mpp_domain)
    if (.not.success) call MOM_error(FATAL, "Failed to open "//trim(filename))

    ! Find the matching case-insensitive variable name in the file and prepare to read it.
    call prepare_to_read_var(fileobj_DD, fieldname, "read_field_0d: ", filename, &
                             var_to_read, has_time_dim, timelevel)

    ! Read the data.
    if (present(timelevel) .and. has_time_dim) then
      call fms2_read_data(fileobj_DD, var_to_read, data, unlim_dim_level=timelevel)
    else
      call fms2_read_data(fileobj_DD, var_to_read, data)
    endif

    ! Close the file-set.
    if (check_if_open(fileobj_DD)) call fms2_close_file(fileobj_DD)
  else
    ! Open the FMS2 file-set.
    success = fms2_open_file(fileObj, trim(filename), "read")
    if (.not.success) call MOM_error(FATAL, "Failed to open "//trim(filename))

    ! Find the matching case-insensitive variable name in the file, and determine whether it
    ! has a time dimension.
    call find_varname_in_file(fileObj, fieldname, "read_field_0d: ", filename, &
                              var_to_read, has_time_dim, timelevel)

    ! Read the data.
    if (present(timelevel) .and. has_time_dim) then
      call fms2_read_data(fileobj, var_to_read, data, unlim_dim_level=timelevel)
    else
      call fms2_read_data(fileobj, var_to_read, data)
    endif

    ! Close the file-set.
    if (check_if_open(fileobj)) call fms2_close_file(fileobj)
  endif

  if (present(scale)) then ; if (scale /= 1.0) then
    data = scale*data
  endif ; endif

end procedure read_field_0d
module procedure read_field_1d
  type(FmsNetcdfFile_t)       :: fileObj ! A handle to a non-domain-decomposed file
  type(FmsNetcdfDomainFile_t) :: fileobj_DD ! A handle to a domain-decomposed file object
  character(len=96) :: var_to_read ! Name of variable to read from the netcdf file
  logical :: has_time_dim          ! True if the variable has an unlimited time axis.
  logical :: success               ! True if the file was successfully opened
  if (present(MOM_Domain)) then
    ! Open the FMS2 file-set.
    success = fms2_open_file(fileobj_DD, filename, "read", MOM_domain%mpp_domain)
    if (.not.success) call MOM_error(FATAL, "Failed to open "//trim(filename))

    ! Find the matching case-insensitive variable name in the file and prepare to read it.
    call prepare_to_read_var(fileobj_DD, fieldname, "read_field_1d: ", filename, &
                             var_to_read, has_time_dim, timelevel)

    ! Read the data.
    if (present(timelevel) .and. has_time_dim) then
      call fms2_read_data(fileobj_DD, var_to_read, data, unlim_dim_level=timelevel)
    else
      call fms2_read_data(fileobj_DD, var_to_read, data)
    endif

    ! Close the file-set.
    if (check_if_open(fileobj_DD)) call fms2_close_file(fileobj_DD)
  else
    ! Open the FMS2 file-set.
    success = fms2_open_file(fileObj, trim(filename), "read")
    if (.not.success) call MOM_error(FATAL, "Failed to open "//trim(filename))

    ! Find the matching case-insensitive variable name in the file, and determine whether it
    ! has a time dimension.
    call find_varname_in_file(fileObj, fieldname, "read_field_1d: ", filename, &
                              var_to_read, has_time_dim, timelevel)

    ! Read the data.
    if (present(timelevel) .and. has_time_dim) then
      call fms2_read_data(fileobj, var_to_read, data, unlim_dim_level=timelevel)
    else
      call fms2_read_data(fileobj, var_to_read, data)
    endif

    ! Close the file-set.
    if (check_if_open(fileobj)) call fms2_close_file(fileobj)
  endif

  if (present(scale)) then ; if (scale /= 1.0) then
    data(:) = scale*data(:)
  endif ; endif

end procedure read_field_1d
module procedure read_field_2d
  type(FmsNetcdfDomainFile_t) :: fileobj ! A handle to a domain-decomposed file object
  character(len=96) :: var_to_read ! Name of variable to read from the netcdf file
  logical :: has_time_dim          ! True if the variable has an unlimited time axis.
  logical :: success               ! True if the file was successfully opened
  success = fms2_open_file(fileobj, filename, "read", MOM_domain%mpp_domain)
  if (.not.success) call MOM_error(FATAL, "Failed to open "//trim(filename))

  ! Find the matching case-insensitive variable name in the file and prepare to read it.
  call prepare_to_read_var(fileobj, fieldname, "read_field_2d: ", filename, &
                           var_to_read, has_time_dim, timelevel, position)

  ! Read the data.
  if (present(timelevel) .and. has_time_dim) then
    call fms2_read_data(fileobj, var_to_read, data, unlim_dim_level=timelevel)
  else
    call fms2_read_data(fileobj, var_to_read, data)
  endif

  ! Close the file-set.
  if (check_if_open(fileobj)) call fms2_close_file(fileobj)

  if (present(scale)) then ; if (scale /= 1.0) then
    call rescale_comp_data(MOM_Domain, data, scale)
  endif ; endif

end procedure read_field_2d
module procedure read_field_2d_region
  type(FmsNetcdfFile_t)       :: fileObj ! A handle to a non-domain-decomposed file
  type(FmsNetcdfDomainFile_t) :: fileobj_DD ! A handle to a domain-decomposed file object
  character(len=96) :: var_to_read ! Name of variable to read from the netcdf file
  logical :: success               ! True if the file was successfully opened
  if (present(MOM_Domain)) then
    ! Open the FMS2 file-set.
    success = fms2_open_file(fileobj_DD, filename, "read", MOM_domain%mpp_domain)
    if (.not.success) call MOM_error(FATAL, "Failed to open "//trim(filename))

    ! Find the matching case-insensitive variable name in the file and prepare to read it.
    call prepare_to_read_var(fileobj_DD, fieldname, "read_field_2d_region: ", &
                             filename, var_to_read)

    ! Read the data.
    call fms2_read_data(fileobj_DD, var_to_read, data, corner=start(1:2), edge_lengths=nread(1:2))

    ! Close the file-set.
    if (check_if_open(fileobj_DD)) call fms2_close_file(fileobj_DD)
  else
    ! Open the FMS2 file-set.
    success = fms2_open_file(fileObj, trim(filename), "read")
    if (.not.success) call MOM_error(FATAL, "Failed to open "//trim(filename))

    ! Find the matching case-insensitive variable name in the file, and determine whether it
    ! has a time dimension.
    call find_varname_in_file(fileObj, fieldname, "read_field_2d_region: ", filename, var_to_read)

    ! Read the data.
    call fms2_read_data(fileobj, var_to_read, data, corner=start(1:2), edge_lengths=nread(1:2))

    ! Close the file-set.
    if (check_if_open(fileobj)) call fms2_close_file(fileobj)
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
  type(FmsNetcdfDomainFile_t) :: fileobj ! A handle to a domain-decomposed file object
  character(len=96) :: var_to_read ! Name of variable to read from the netcdf file
  logical :: has_time_dim          ! True if the variable has an unlimited time axis.
  logical :: success               ! True if the file was successfully opened
  success = fms2_open_file(fileobj, filename, "read", MOM_domain%mpp_domain)
  if (.not.success) call MOM_error(FATAL, "Failed to open "//trim(filename))

  ! Find the matching case-insensitive variable name in the file and prepare to read it.
  call prepare_to_read_var(fileobj, fieldname, "read_field_3d: ", filename, &
                           var_to_read, has_time_dim, timelevel, position)

  ! Read the data.
  if (present(timelevel) .and. has_time_dim) then
    call fms2_read_data(fileobj, var_to_read, data, unlim_dim_level=timelevel)
  else
    call fms2_read_data(fileobj, var_to_read, data)
  endif

  ! Close the file-set.
  if (check_if_open(fileobj)) call fms2_close_file(fileobj)

  if (present(scale)) then ; if (scale /= 1.0) then
    call rescale_comp_data(MOM_Domain, data, scale)
  endif ; endif

end procedure read_field_3d
module procedure read_field_3d_region
  type(FmsNetcdfFile_t)       :: fileObj ! A handle to a non-domain-decomposed file
  type(FmsNetcdfDomainFile_t) :: fileobj_DD ! A handle to a domain-decomposed file object
  character(len=96) :: var_to_read ! Name of variable to read from the netcdf file
  logical :: success               ! True if the file was successfully opened
  if (present(MOM_Domain)) then
    ! Open the FMS2 file-set.
    success = fms2_open_file(fileobj_DD, filename, "read", MOM_domain%mpp_domain)
    if (.not.success) call MOM_error(FATAL, "Failed to open "//trim(filename))

    ! Find the matching case-insensitive variable name in the file and prepare to read it.
    call prepare_to_read_var(fileobj_DD, fieldname, "read_field_2d_region: ", &
                             filename, var_to_read)

    ! Read the data.
    call fms2_read_data(fileobj_DD, var_to_read, data, corner=start(1:3), edge_lengths=nread(1:3))

    ! Close the file-set.
    if (check_if_open(fileobj_DD)) call fms2_close_file(fileobj_DD)
  else
    ! Open the FMS2 file-set.
    success = fms2_open_file(fileObj, trim(filename), "read")
    if (.not.success) call MOM_error(FATAL, "Failed to open "//trim(filename))

    ! Find the matching case-insensitive variable name in the file, and determine whether it
    ! has a time dimension.
    call find_varname_in_file(fileObj, fieldname, "read_field_2d_region: ", filename, var_to_read)

    ! Read the data.
    call fms2_read_data(fileobj, var_to_read, data, corner=start(1:3), edge_lengths=nread(1:3))

    ! Close the file-set.
    if (check_if_open(fileobj)) call fms2_close_file(fileobj)
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
  type(FmsNetcdfDomainFile_t) :: fileobj ! A handle to a domain-decomposed file object
  logical :: has_time_dim          ! True if the variable has an unlimited time axis.
  character(len=96) :: var_to_read ! Name of variable to read from the netcdf file
  logical :: success  ! True if the file was successfully opened
  success = fms2_open_file(fileobj, filename, "read", MOM_domain%mpp_domain)
  if (.not.success) call MOM_error(FATAL, "Failed to open "//trim(filename))

  ! Find the matching case-insensitive variable name in the file and prepare to read it.
  call prepare_to_read_var(fileobj, fieldname, "read_field_4d: ", filename, &
                           var_to_read, has_time_dim, timelevel, position)

  ! Read the data.
  if (present(timelevel) .and. has_time_dim) then
    call fms2_read_data(fileobj, var_to_read, data, unlim_dim_level=timelevel)
  else
    call fms2_read_data(fileobj, var_to_read, data)
  endif

  ! Close the file-set.
  if (check_if_open(fileobj)) call fms2_close_file(fileobj)

  if (present(scale)) then ; if (scale /= 1.0) then
    call rescale_comp_data(MOM_Domain, data, scale)
  endif ; endif

end procedure read_field_4d
module procedure read_field_0d_int
  type(FmsNetcdfFile_t) :: fileObj ! A handle to a non-domain-decomposed file
  logical :: has_time_dim          ! True if the variable has an unlimited time axis.
  character(len=96) :: var_to_read ! Name of variable to read from the netcdf file
  logical :: success               ! If true, the file was opened successfully
  success = fms2_open_file(fileObj, trim(filename), "read")
  if (.not.success) call MOM_error(FATAL, "Failed to open "//trim(filename))

  ! Find the matching case-insensitive variable name in the file, and determine whether it
  ! has a time dimension.
  call find_varname_in_file(fileObj, fieldname, "read_field_0d_int: ", filename, &
                            var_to_read, has_time_dim, timelevel)

  ! Read the data.
  if (present(timelevel) .and. has_time_dim) then
    call fms2_read_data(fileobj, var_to_read, data, unlim_dim_level=timelevel)
  else
    call fms2_read_data(fileobj, var_to_read, data)
  endif

  ! Close the file-set.
  if (check_if_open(fileobj)) call fms2_close_file(fileobj)
end procedure read_field_0d_int
module procedure read_field_1d_int
  type(FmsNetcdfFile_t) :: fileObj ! A handle to a non-domain-decomposed file for obtaining information
  logical :: has_time_dim          ! True if the variable has an unlimited time axis.
  character(len=96) :: var_to_read ! Name of variable to read from the netcdf file
  logical :: success               ! If true, the file was opened successfully
  success = fms2_open_file(fileObj, trim(filename), "read")
  if (.not.success) call MOM_error(FATAL, "Failed to open "//trim(filename))

  ! Find the matching case-insensitive variable name in the file, and determine whether it
  ! has a time dimension.
  call find_varname_in_file(fileObj, fieldname, "read_field_1d_int: ", filename, &
                            var_to_read, has_time_dim, timelevel)

  ! Read the data.
  if (present(timelevel) .and. has_time_dim) then
    call fms2_read_data(fileobj, var_to_read, data, unlim_dim_level=timelevel)
  else
    call fms2_read_data(fileobj, var_to_read, data)
  endif

  ! Close the file-set.
  if (check_if_open(fileobj)) call fms2_close_file(fileobj)
end procedure read_field_1d_int
module procedure read_vector_2d
  type(FmsNetcdfDomainFile_t) :: fileobj ! A handle to a domain-decomposed file object
  logical :: has_time_dim           ! True if the variables have an unlimited time axis.
  character(len=96) :: u_var, v_var ! Name of u and v variables to read from the netcdf file
  logical :: success                ! True if the file was successfully opened
  integer :: u_pos, v_pos           ! Flags indicating the positions of the u- and v- components.
  u_pos = EAST_FACE ; v_pos = NORTH_FACE
  if (present(stagger)) then
    if (stagger == CGRID_NE) then ; u_pos = EAST_FACE ; v_pos = NORTH_FACE
    elseif (stagger == BGRID_NE) then ; u_pos = CORNER ; v_pos = CORNER
    elseif (stagger == AGRID) then ; u_pos = CENTER ; v_pos = CENTER ; endif
  endif

  ! Open the FMS2 file-set.
  success = fms2_open_file(fileobj, filename, "read", MOM_domain%mpp_domain)
  if (.not.success) call MOM_error(FATAL, "Failed to open "//trim(filename))

  ! Find the matching case-insensitive u- and v-variable names in the file and prepare to read them.
  call prepare_to_read_var(fileobj, u_fieldname, "read_vector_2d: ", filename, &
                           u_var, has_time_dim, timelevel, position=u_pos)
  call prepare_to_read_var(fileobj, v_fieldname, "read_vector_2d: ", filename, &
                           v_var, has_time_dim, timelevel, position=v_pos)

  ! Read the u-data and v-data. There would already been an error message for one
  ! of the variables if they are inconsistent in having an unlimited dimension.
  if (present(timelevel) .and. has_time_dim) then
    call fms2_read_data(fileobj, u_var, u_data, unlim_dim_level=timelevel)
    call fms2_read_data(fileobj, v_var, v_data, unlim_dim_level=timelevel)
  else
    call fms2_read_data(fileobj, u_var, u_data)
    call fms2_read_data(fileobj, v_var, v_data)
  endif

  ! Close the file-set.
  if (check_if_open(fileobj)) call fms2_close_file(fileobj)

  if (present(scale)) then ; if (scale /= 1.0) then
    call rescale_comp_data(MOM_Domain, u_data, scale)
    call rescale_comp_data(MOM_Domain, v_data, scale)
  endif ; endif

end procedure read_vector_2d
module procedure read_vector_3d
  type(FmsNetcdfDomainFile_t) :: fileobj ! A handle to a domain-decomposed file object
  logical :: has_time_dim           ! True if the variables have an unlimited time axis.
  character(len=96) :: u_var, v_var ! Name of u and v variables to read from the netcdf file
  logical :: success                ! True if the file was successfully opened
  integer :: u_pos, v_pos           ! Flags indicating the positions of the u- and v- components.
  u_pos = EAST_FACE ; v_pos = NORTH_FACE
  if (present(stagger)) then
    if (stagger == CGRID_NE) then ; u_pos = EAST_FACE ; v_pos = NORTH_FACE
    elseif (stagger == BGRID_NE) then ; u_pos = CORNER ; v_pos = CORNER
    elseif (stagger == AGRID) then ; u_pos = CENTER ; v_pos = CENTER ; endif
  endif

  ! Open the FMS2 file-set.
  success = fms2_open_file(fileobj, filename, "read", MOM_domain%mpp_domain)
  if (.not.success) call MOM_error(FATAL, "Failed to open "//trim(filename))

  ! Find the matching case-insensitive u- and v-variable names in the file and prepare to read them.
  call prepare_to_read_var(fileobj, u_fieldname, "read_vector_3d: ", filename, &
                           u_var, has_time_dim, timelevel, position=u_pos)
  call prepare_to_read_var(fileobj, v_fieldname, "read_vector_3d: ", filename, &
                           v_var, has_time_dim, timelevel, position=v_pos)

  ! Read the u-data and v-data, dangerously assuming either both or neither have time dimensions.
  ! There would already been an error message for one of the variables if they are inconsistent.
  if (present(timelevel) .and. has_time_dim) then
    call fms2_read_data(fileobj, u_var, u_data, unlim_dim_level=timelevel)
    call fms2_read_data(fileobj, v_var, v_data, unlim_dim_level=timelevel)
  else
    call fms2_read_data(fileobj, u_var, u_data)
    call fms2_read_data(fileobj, v_var, v_data)
  endif

  ! Close the file-set.
  if (check_if_open(fileobj)) call fms2_close_file(fileobj)

  if (present(scale)) then ; if (scale /= 1.0) then
    call rescale_comp_data(MOM_Domain, u_data, scale)
    call rescale_comp_data(MOM_Domain, v_data, scale)
  endif ; endif

end procedure read_vector_3d
module procedure find_varname_in_file
  logical :: variable_found ! Is a case-insensitive version of the variable found in the netCDF file?
  character(len=256), allocatable, dimension(:) :: var_names ! The names of all the variables in the netCDF file
  character(len=256), allocatable :: dim_names(:) ! The names of a variable's dimensions
  integer :: nvars          ! The number of variables in the file
  integer :: dim_unlim_size ! The current size of the unlimited (time) dimension in the file.
  integer :: num_var_dims   ! The number of dimensions a variable has in the file.
  integer :: time_dim       ! The position of the unlimited (time) dimension for a variable, or -1
  integer :: i
  if (.not.check_if_open(fileobj))  &
    call MOM_error(FATAL, trim(err_header)//trim(filename)//" was not open in call to find_varname_in_file.")

  ! Search for the variable in the file, looking for the case-sensitive name first.
  if (variable_exists(fileobj, trim(fieldname))) then
    var_to_read = trim(fieldname)
  else ! Look for case-insensitive variable name matches.
    nvars = get_num_variables(fileobj)
    if (nvars < 1) call MOM_error(FATAL, "nvars is less than 1 for file "//trim(filename))
    allocate(var_names(nvars))
    call get_variable_names(fileobj, var_names)

    ! search for the variable in the file
    variable_found = .false.
    do i=1,nvars
      if (lowercase(trim(var_names(i))) == lowercase(trim(fieldname))) then
        variable_found = .true.
        var_to_read = trim(var_names(i))
        exit
      endif
    enddo
    if (.not.(variable_found)) &
      call MOM_error(FATAL, trim(err_header)//trim(fieldname)//" not found in "//trim(filename))
    deallocate(var_names)
  endif

  ! FMS2 can not handle a timelevel argument if the variable does not have one in the file,
  ! so some error checking and logic are required.
  if (present(has_time_dim) .or. present(timelevel)) then
    time_dim = -1

    num_var_dims = get_variable_num_dimensions(fileobj, trim(var_to_read))
    allocate(dim_names(num_var_dims)) ; dim_names(:) = ""
    call get_variable_dimension_names(fileobj, trim(var_to_read), dim_names)

    do i=1,num_var_dims
      if (is_dimension_unlimited(fileobj, dim_names(i))) then
        time_dim = i
        if (present(timelevel)) then
          call get_dimension_size(fileobj, dim_names(i), dim_unlim_size)
          if ((timelevel > dim_unlim_size) .and. is_root_PE()) call MOM_error(FATAL, &
                trim(err_header)//"Attempting to read a time level of "//trim(var_to_read)//&
                " that exceeds the size of the time dimension in "//trim(filename))
        endif
        exit
      endif
    enddo
    deallocate(dim_names)

    if (present(timelevel) .and. (time_dim < 0) .and. is_root_PE()) &
      call MOM_error(WARNING, trim(err_header)//"time level specified, but the variable "//&
                   trim(var_to_read)//" does not have an unlimited dimension in "//trim(filename))
    if ((.not.present(timelevel)) .and. (time_dim > 0) .and. is_root_PE()) &
      call MOM_error(WARNING, trim(err_header)//"The variable "//trim(var_to_read)//&
                    " has an unlimited dimension in "//trim(filename)//" but no time level is specified.")
    if (present(has_time_dim)) has_time_dim = (time_dim > 0)
  endif

end procedure find_varname_in_file
module procedure prepare_to_read_var
  logical :: variable_found ! Is a case-insensitive version of the variable found in the netCDF file?
  character(len=256), allocatable, dimension(:) :: var_names ! The names of all the variables in the netCDF file
  character(len=256), allocatable :: dim_names(:) ! The names of a variable's dimensions
  integer :: nvars          ! The number of variables in the file.
  integer :: dim_unlim_size ! The current size of the unlimited (time) dimension in the file.
  integer :: num_var_dims   ! The number of dimensions a variable has in the file.
  integer :: time_dim       ! The position of the unlimited (time) dimension for a variable, or -1
  integer :: i
  if (.not.check_if_open(fileobj))  &
    call MOM_error(FATAL, trim(err_header)//trim(filename)//" was not open in call to prepare_to_read_var.")

  ! Search for the variable in the file, looking for the case-sensitive name first.
  if (variable_exists(fileobj, trim(fieldname))) then
    var_to_read = trim(fieldname)
  else  ! Look for case-insensitive variable name matches.
    nvars = get_num_variables(fileobj)
    if (nvars < 1) call MOM_error(FATAL, "nvars is less than 1 for file "//trim(filename))
    allocate(var_names(nvars))
    call get_variable_names(fileobj, var_names)

    variable_found = .false.
    do i=1,nvars
      if (lowercase(trim(var_names(i))) == lowercase(trim(fieldname))) then
        variable_found = .true.
        var_to_read = trim(var_names(i))
        exit
      endif
    enddo
    if (.not.(variable_found)) &
      call MOM_error(FATAL, trim(err_header)//trim(fieldname)//" not found in "//trim(filename))
    deallocate(var_names)
  endif

  ! FMS2 can not handle a timelevel argument if the variable does not have one in the file,
  ! so some error checking and logic are required.
  if (present(has_time_dim) .or. present(timelevel)) then
    time_dim = -1

    num_var_dims = get_variable_num_dimensions(fileobj, trim(var_to_read))
    allocate(dim_names(num_var_dims)) ; dim_names(:) = ""
    call get_variable_dimension_names(fileobj, trim(var_to_read), dim_names)

    do i=1,num_var_dims
      if (is_dimension_unlimited(fileobj, dim_names(i))) then
        time_dim = i
        if (present(timelevel)) then
          call get_dimension_size(fileobj, dim_names(i), dim_unlim_size)
          if ((timelevel > dim_unlim_size) .and. is_root_PE()) call MOM_error(FATAL, &
                trim(err_header)//"Attempting to read a time level of "//trim(var_to_read)//&
                " that exceeds the size of the time dimension in "//trim(filename))
        endif
        exit
      endif
    enddo
    deallocate(dim_names)

    if (present(timelevel) .and. (time_dim < 0) .and. is_root_PE()) &
      call MOM_error(WARNING, trim(err_header)//"time level specified, but the variable "//&
                   trim(var_to_read)//" does not have an unlimited dimension in "//trim(filename))
    if ((.not.present(timelevel)) .and. (time_dim > 0) .and. is_root_PE()) &
      call MOM_error(WARNING, trim(err_header)//"The variable "//trim(var_to_read)//&
                    " has an unlimited dimension in "//trim(filename)//" but no time level is specified.")
    if (present(has_time_dim)) has_time_dim = (time_dim > 0)
  endif

  ! Registering the variable axes essentially just specifies the discrete position of this variable.
  call MOM_register_variable_axes(fileobj, var_to_read, filename, position)

end procedure prepare_to_read_var
module procedure MOM_register_variable_axes
  character(len=256), allocatable, dimension(:) :: dim_names ! variable dimension names
  integer, allocatable, dimension(:) :: dimSizes ! variable dimension sizes
  logical, allocatable, dimension(:) :: is_x ! Is this a (likely domain-decomposed) x-axis
  logical, allocatable, dimension(:) :: is_y ! Is this a (likely domain-decomposed) y-axis
  logical, allocatable, dimension(:) :: is_t ! Is this a time axis or another unlimited axis
  integer :: ndims ! number of dimensions
  integer :: xPos, yPos ! Discrete positions for x and y axes. Default is CENTER
  integer :: i
  xPos = CENTER ; yPos = CENTER
  if (present(position)) then
    if ((position == CORNER) .or. (position == EAST_FACE)) xPos = EAST_FACE
    if ((position == CORNER) .or. (position == NORTH_FACE)) yPos = NORTH_FACE
  endif

  ! get variable dimension names and lengths
  ndims = get_variable_num_dimensions(fileObj, trim(variableName))
  allocate(dimSizes(ndims))
  allocate(dim_names(ndims))
  allocate(is_x(ndims)) ; is_x(:) = .false.
  allocate(is_y(ndims)) ; is_y(:) = .false.
  allocate(is_t(ndims)) ; is_t(:) = .false.
  call get_variable_size(fileObj, trim(variableName), dimSizes)
  call get_variable_dimension_names(fileObj, trim(variableName), dim_names)
  call categorize_axes(fileObj, filename, ndims, dim_names, is_x, is_y, is_t)

  ! register the axes
  do i=1,ndims
    if ( .not.is_dimension_registered(fileobj, trim(dim_names(i))) ) then
      if (is_x(i)) then
        call register_axis(fileObj, trim(dim_names(i)), "x", domain_position=xPos)
      elseif (is_y(i)) then
        call register_axis(fileObj, trim(dim_names(i)), "y", domain_position=yPos)
      else
        call register_axis(fileObj, trim(dim_names(i)), dimSizes(i))
      endif
    endif
  enddo

  deallocate(dimSizes, dim_names, is_x, is_y, is_t)
end procedure MOM_register_variable_axes
module procedure categorize_axes
  character(len=128) :: cartesian ! A flag indicating a Cartesian direction - usually a single character.
  character(len=512) :: dim_list  ! A concatenated list of dimension names.
  character(len=128) :: units ! units corresponding to a specific variable dimension
  logical :: x_found, y_found ! Indicate whether an x- or y- dimension have been found.
  integer :: i
  x_found = .false. ; y_found = .false.
  is_x(:) = .false. ; is_y(:) = .false.
  do i=1,ndims
    is_t(i) = is_dimension_unlimited(fileObj, trim(dim_names(i)))
    ! First look for indicative variable attributes
    if (.not.is_t(i)) then
      if (variable_exists(fileobj, trim(dim_names(i)))) then
        cartesian = ""
        if (variable_att_exists(fileobj, trim(dim_names(i)), "cartesian_axis")) then
          call get_variable_attribute(fileobj, trim(dim_names(i)), "cartesian_axis", cartesian(1:1))
        elseif (variable_att_exists(fileobj, trim(dim_names(i)), "axis")) then
          call get_variable_attribute(fileobj, trim(dim_names(i)), "axis", cartesian(1:1))
        endif
        cartesian = adjustl(cartesian)
        if ((index(cartesian, "X") == 1) .or. (index(cartesian, "x") == 1)) is_x(i) = .true.
        if ((index(cartesian, "Y") == 1) .or. (index(cartesian, "y") == 1)) is_y(i) = .true.
        if ((index(cartesian, "T") == 1) .or. (index(cartesian, "t") == 1)) is_t(i) = .true.
      endif
    endif
    if (is_x(i)) x_found = .true.
    if (is_y(i)) y_found = .true.
  enddo

  if (.not.(x_found .and. y_found)) then
    ! Next look for hints from axis names for uncharacterized axes
    do i=1,ndims ; if (.not.(is_x(i) .or. is_y(i) .or. is_t(i))) then
      call categorize_axis_from_name(dim_names(i), is_x(i), is_y(i))
      if (is_x(i)) x_found = .true.
      if (is_y(i)) y_found = .true.
    endif ; enddo
  endif

  if (.not.(x_found .and. y_found)) then
    ! Look for hints from CF-compliant axis units for uncharacterized axes
    do i=1,ndims ; if (.not.(is_x(i) .or. is_y(i) .or. is_t(i))) then
      if (variable_exists(fileobj, trim(dim_names(i)))) then
        call get_variable_units(fileobj, trim(dim_names(i)), units)
        call categorize_axis_from_units(units, is_x(i), is_y(i))
      endif
      if (is_x(i)) x_found = .true.
      if (is_y(i)) y_found = .true.
    endif ; enddo
  endif

  if (.not.(x_found .and. y_found) .and. ((ndims>2) .or. ((ndims==2) .and. .not.is_t(ndims)))) then
    ! This is a case where one would expect to find x-and y-dimensions, but none have been found.
    if (is_root_pe()) then
      dim_list = trim(dim_names(1))//", "//trim(dim_names(2))
      do i=3,ndims ; dim_list = trim(dim_list)//", "//trim(dim_names(i)) ; enddo
      call MOM_error(WARNING, "categorize_axes: Failed to identify x- and y- axes in the axis list ("//&
                     trim(dim_list)//") of a variable being read from "//trim(filename))
    endif
  endif

end procedure categorize_axes
module procedure categorize_axis_from_units
  is_x = .false. ; is_y = .false.
  select case (lowercase(trim(unit_string)))
    case ("degrees_north"); is_y = .true.
    case ("degree_north") ; is_y = .true.
    case ("degrees_n")    ; is_y = .true.
    case ("degree_n")     ; is_y = .true.
    case ("degreen")      ; is_y = .true.
    case ("degreesn")     ; is_y = .true.
    case ("degrees_east") ; is_x = .true.
    case ("degree_east")  ; is_x = .true.
    case ("degreese")     ; is_x = .true.
    case ("degreee")      ; is_x = .true.
    case ("degree_e")     ; is_x = .true.
    case ("degrees_e")    ; is_x = .true.
    case default ; is_x = .false. ; is_y = .false.
  end select

end procedure categorize_axis_from_units
module procedure categorize_axis_from_name
  is_x = .false. ; is_y = .false.
  select case(trim(lowercase(dimname)))
    case ("grid_x_t")  ; is_x = .true.
    case ("nx")        ; is_x = .true.
    case ("nxp")       ; is_x = .true.
    case ("longitude") ; is_x = .true.
    case ("long")      ; is_x = .true.
    case ("lon")       ; is_x = .true.
    case ("lonh")      ; is_x = .true.
    case ("lonq")      ; is_x = .true.
    case ("xh")        ; is_x = .true.
    case ("xq")        ; is_x = .true.
    case ("i")         ; is_x = .true.

    case ("grid_y_t")  ; is_y = .true.
    case ("ny")        ; is_y = .true.
    case ("nyp")       ; is_y = .true.
    case ("latitude")  ; is_y = .true.
    case ("lat")       ; is_y = .true.
    case ("lath")      ; is_y = .true.
    case ("latq")      ; is_y = .true.
    case ("yh")        ; is_y = .true.
    case ("yq")        ; is_y = .true.
    case ("j")         ; is_y = .true.

    case default ; is_x = .false. ; is_y = .false.
  end select

end procedure categorize_axis_from_name
module procedure write_field_4d
  integer :: time_index
  if (present(tstamp)) then
    time_index = write_time_if_later(IO_handle, tstamp)
    call write_data(IO_handle%fileobj, trim(field_md%name), field, unlim_dim_level=time_index)
  else
    call write_data(IO_handle%fileobj, trim(field_md%name), field)
  endif
end procedure write_field_4d
module procedure write_field_3d
  integer :: time_index
  if (present(tstamp)) then
    time_index = write_time_if_later(IO_handle, tstamp)
    call write_data(IO_handle%fileobj, trim(field_md%name), field, unlim_dim_level=time_index)
  else
    call write_data(IO_handle%fileobj, trim(field_md%name), field)
  endif
end procedure write_field_3d
module procedure write_field_2d
  integer :: time_index
  if (present(tstamp)) then
    time_index = write_time_if_later(IO_handle, tstamp)
    call write_data(IO_handle%fileobj, trim(field_md%name), field, unlim_dim_level=time_index)
  else
    call write_data(IO_handle%fileobj, trim(field_md%name), field)
  endif
end procedure write_field_2d
module procedure write_field_1d
  integer :: time_index
  if (present(tstamp)) then
    time_index = write_time_if_later(IO_handle, tstamp)
    call write_data(IO_handle%fileobj, trim(field_md%name), field, unlim_dim_level=time_index)
  else
    call write_data(IO_handle%fileobj, trim(field_md%name), field)
  endif
end procedure write_field_1d
module procedure write_field_0d
  integer :: time_index
  if (present(tstamp)) then
    time_index = write_time_if_later(IO_handle, tstamp)
    call write_data(IO_handle%fileobj, trim(field_md%name), field, unlim_dim_level=time_index)
  else
    call write_data(IO_handle%fileobj, trim(field_md%name), field)
  endif
end procedure write_field_0d
module procedure write_time_if_later
  character(len=256) :: dim_unlim_name ! name of the unlimited dimension in the file
  if ((field_time > IO_handle%file_time) .or. (IO_handle%num_times == 0)) then
    IO_handle%file_time = field_time
    IO_handle%num_times = IO_handle%num_times + 1
    dim_unlim_name = find_unlimited_dimension_name(IO_handle%fileobj)
    if (len_trim(dim_unlim_name) > 0) &
      call write_data(IO_handle%fileobj, trim(dim_unlim_name), [field_time], &
                      corner=[IO_handle%num_times], edge_lengths=[1])
  endif

  write_time_if_later = IO_handle%num_times
end procedure write_time_if_later
module procedure MOM_write_axis
  integer :: is, ie
  if (axis%domain_decomposed) then
    ! FMS2 does not domain-decompose 1d arrays, so we explicitly slice it
    call get_global_io_domain_indices(IO_handle%fileobj, trim(axis%name), is, ie)
    call write_data(IO_handle%fileobj, trim(axis%name), axis%ax_data(is:ie))
  else
    call write_data(IO_handle%fileobj, trim(axis%name), axis%ax_data)
  endif
end procedure MOM_write_axis
module procedure write_metadata_axis
  character(len=:), allocatable :: cart ! A left-adjusted and trimmed copy of cartesian
  logical :: is_x, is_y, is_t  ! If true, this is a domain-decomposed axis in one of the directions.
  integer :: position    ! A flag indicating the axis staggering position.
  integer :: i, isc, iec, global_size
  if (is_dimension_registered(IO_handle%fileobj, trim(name))) then
    call MOM_error(FATAL, "write_metadata_axis was called more than once for axis "//trim(name)//&
                          " in file "//trim(IO_handle%filename))
    return
  endif

  axis%name = trim(name)
  if (present(data) .and. allocated(axis%ax_data)) call MOM_error(FATAL, &
        "Data is already allocated in a call to write_metadata_axis for axis "//&
        trim(name)//" in file "//trim(IO_handle%filename))

  is_x = .false. ; is_y = .false. ; is_t = .false.
  position = CENTER
  if (present(cartesian)) then
    cart = trim(adjustl(cartesian))
    if ((index(cart, "X") == 1) .or. (index(cart, "x") == 1)) is_x = .true.
    if ((index(cart, "Y") == 1) .or. (index(cart, "y") == 1)) is_y = .true.
    if ((index(cart, "T") == 1) .or. (index(cart, "t") == 1)) is_t = .true.
  endif

  ! For now, we assume that all horizontal axes are domain-decomposed.
  if (is_x .or. is_y) &
    axis%domain_decomposed = .true.

  if (is_x) then
    if (present(edge_axis)) then ; if (edge_axis) position = EAST_FACE ; endif
    call register_axis(IO_handle%fileobj, trim(name), 'x', domain_position=position)
  elseif (is_y) then
    if (present(edge_axis)) then ; if (edge_axis) position = NORTH_FACE ; endif
    call register_axis(IO_handle%fileobj, trim(name), 'y', domain_position=position)
  elseif (is_t .and. .not.present(data)) then
    ! This is the unlimited (time) dimension.
    call register_axis(IO_handle%fileobj, trim(name), unlimited)
  else
    if (.not.present(data)) call MOM_error(FATAL,"MOM_io:register_diagnostic_axis: "//&
                      "An axis_length argument is required to register the axis "//trim(name))
    call register_axis(IO_handle%fileobj, trim(name), size(data))
  endif

  if (present(data)) then
    ! With FMS2, the data for the axis labels has to match the computational domain on this PE.
    if (present(domain)) then
      ! The commented-out code on the next ~11 lines runs but there is missing data in the output file
      ! call mpp_get_compute_domain(domain, isc, iec)
      ! call mpp_get_global_domain(domain, size=global_size)
      ! if (size(data) == global_size) then
      !   allocate(axis%ax_data(iec+1-isc)) ; axis%ax_data(:) = data(isc:iec)
      !   ! A simpler set of labels: do i=1,iec-isc ; axis%ax_data(i) = real(isc + i) - 1.0 ; enddo
      ! elseif (size(data) == global_size+1) then
      !   ! This is an edge axis.  Note the effective SW indexing convention here.
      !   allocate(axis%ax_data(iec+2-isc)) ; axis%ax_data(:) = data(isc:iec+1)
      !   ! A simpler set of labels: do i=1,iec+1-isc ; axis%ax_data(i) = real(isc + i) - 1.5 ; enddo
      ! else
      !   call MOM_error(FATAL, "Unexpected size of data for "//trim(name)//" in write_metadata_axis.")
      ! endif

      ! This works for a simple 1x1 IO layout, but gives errors for nontrivial IO layouts
      allocate(axis%ax_data(size(data))) ; axis%ax_data(:) = data(:)

    else  ! Store the entire array of axis labels.
      allocate(axis%ax_data(size(data))) ; axis%ax_data(:) = data(:)
    endif
  endif


  ! Now create the variable that describes this axis.
  call register_field(IO_handle%fileobj, trim(name), "double", dimensions=(/name/))
  if (len_trim(longname) > 0) &
    call register_variable_attribute(IO_handle%fileobj, trim(name), 'long_name', &
                                     trim(longname), len_trim(longname))
  if (len_trim(units) > 0) &
    call register_variable_attribute(IO_handle%fileobj, trim(name), 'units', &
                                     trim(units), len_trim(units))
  if (present(cartesian)) &
    call register_variable_attribute(IO_handle%fileobj, trim(name), 'cartesian_axis', &
                                     trim(cartesian), len_trim(cartesian))
  if (present(sense)) &
    call register_variable_attribute(IO_handle%fileobj, trim(name), 'sense', sense)
end procedure write_metadata_axis
module procedure write_metadata_field
  character(len=256), dimension(size(axes)) :: dim_names ! The names of the dimensions
  character(len=16) :: prec_string     ! A string specifying the precision with which to save this variable
  character(len=64) :: checksum_string ! checksum character array created from checksum argument
  integer :: i, ndims
  ndims = size(axes)
  do i=1,ndims ; dim_names(i) = trim(axes(i)%name) ; enddo
  prec_string = "double" ; if (present(pack)) then ; if (pack > 1) prec_string = "float" ; endif
  call register_field(IO_handle%fileobj, trim(name), trim(prec_string), dimensions=dim_names)
  if (len_trim(longname) > 0) &
    call register_variable_attribute(IO_handle%fileobj, trim(name), 'long_name', &
                                     trim(longname), len_trim(longname))
  if (len_trim(units) > 0) &
    call register_variable_attribute(IO_handle%fileobj, trim(name), 'units', &
                                     trim(units), len_trim(units))
  if (present(standard_name)) &
    call register_variable_attribute(IO_handle%fileobj, trim(name), 'standard_name', &
                                     trim(standard_name), len_trim(standard_name))
  if (present(checksum)) then
    write (checksum_string,'(Z16)') checksum(1) ! Z16 is the hexadecimal format code
    call register_variable_attribute(IO_handle%fileobj, trim(name), "checksum", &
                                     trim(checksum_string), len_trim(checksum_string))
  endif

  ! Store information in the field-type, regardless of which interfaces are used.
  field%name = trim(name)
  field%longname = trim(longname)
  field%units = trim(units)
  field%chksum_read = -1
  field%valid_chksum = .false.

end procedure write_metadata_field
module procedure write_metadata_global
  call register_global_attribute(IO_handle%fileobj, name, attribute, len_trim(attribute))
end procedure write_metadata_global
module procedure find_unlimited_dimension_name
  integer :: ndims
  character(len=256), allocatable :: dim_names(:)
  integer :: i
  ndims = get_num_dimensions(fileobj)
  allocate(dim_names(ndims))
  call get_dimension_names(fileobj, dim_names)

  do i = 1, ndims
    if (is_dimension_unlimited(fileobj, dim_names(i))) then
      label = trim(dim_names(i))
      exit
    endif
  enddo
  deallocate(dim_names)

  if (.not. allocated(label)) &
    label = ''
end procedure find_unlimited_dimension_name
module procedure lowercase
  integer, parameter :: co=iachar('a')-iachar('A') ! case offset
  integer :: k
  lowercase = input_string
  do k=1, len_trim(input_string)
    if (lowercase(k:k) >= 'A' .and. lowercase(k:k) <= 'Z') &
        lowercase(k:k) = achar(ichar(lowercase(k:k))+co)
  enddo
end procedure lowercase
end submodule MOM_io_infra_s
