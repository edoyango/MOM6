submodule (MOM_io) MOM_io_s
  implicit none
contains
module procedure create_file
  type(MOM_infra_file) :: new_file
  type(MOM_field) :: new_fields(novars)
  new_file%handle_infra = IO_handle

  call create_MOM_file(new_file, filename, vars, novars, new_fields, &
      threading=threading, timeunit=timeunit, G=G, dG=dG, GV=GV, &
      checksums=checksums, extra_axes=extra_axes, global_atts=global_atts)

  IO_handle = new_file%handle_infra
  call new_file%get_file_fieldtypes(fields(:novars))
end procedure create_file
module procedure create_MOM_file
  logical        :: use_lath, use_lonh, use_latq, use_lonq, use_time
  logical        :: use_layer, use_int, use_periodic
  logical        :: one_file, domain_set, dim_found
  logical, dimension(:), allocatable :: use_extra_axis
  type(MOM_axis) :: axis_lath, axis_latq, axis_lonh, axis_lonq
  type(MOM_axis) :: axis_layer, axis_int, axis_time, axis_periodic
  type(MOM_axis), dimension(:), allocatable :: more_axes ! Axes generated from extra_axes
  type(MOM_axis) :: axes(5)     ! The axes of a variable
  type(MOM_domain_type), pointer :: Domain => NULL()
  type(domain1d) :: x_domain, y_domain
  integer        :: position, numaxes, pack, thread, k, n, m
  integer        :: num_extra_dims ! The number of extra possible dimensions from extra_axes
  integer        :: isg, ieg, jsg, jeg, IsgB, IegB, JsgB, JegB
  integer        :: var_periods, num_periods=0
  real, dimension(:), allocatable :: axis_val ! Axis label values [various]
  real, pointer, dimension(:) :: &
    gridLatT => NULL(), & ! The latitude of T or B points for the purpose of labeling
    gridLatB => NULL(), & ! the output axes, often in units of [degrees_N] or [km] or [m].
    gridLonT => NULL(), & ! The longitude of T or B points for the purpose of labeling
    gridLonB => NULL()    ! the output axes, often in units of [degrees_E] or [km] or [m].
  character(len=40) :: time_units, x_axis_units, y_axis_units
  character(len=8)  :: t_grid, t_grid_read
  character(len=64) :: ax_name(5)  ! The axis names of a variable
  use_lath  = .false. ; use_lonh     = .false.
  use_latq  = .false. ; use_lonq     = .false.
  use_time  = .false. ; use_periodic = .false.
  use_layer = .false. ; use_int      = .false.
  num_extra_dims = 0
  if (present(extra_axes)) then
    num_extra_dims = size(extra_axes)
    if (num_extra_dims > 0) then
      allocate(use_extra_axis(num_extra_dims)) ; use_extra_axis = .false.
      allocate(more_axes(num_extra_dims))
    endif
  endif

  thread = SINGLE_FILE
  if (PRESENT(threading)) thread = threading

  domain_set = .false.
  if (present(G)) then
    domain_set = .true. ; Domain => G%Domain
    gridLatT => G%gridLatT ; gridLatB => G%gridLatB
    gridLonT => G%gridLonT ; gridLonB => G%gridLonB
    x_axis_units = G%x_axis_units ; y_axis_units = G%y_axis_units
    isg = G%isg ; ieg = G%ieg ; jsg = G%jsg ; jeg = G%jeg
    IsgB = G%IsgB ; IegB = G%IegB ; JsgB = G%JsgB ; JegB = G%JegB
  elseif (present(dG)) then
    domain_set = .true. ; Domain => dG%Domain
    gridLatT => dG%gridLatT ; gridLatB => dG%gridLatB
    gridLonT => dG%gridLonT ; gridLonB => dG%gridLonB
    x_axis_units = dG%x_axis_units ; y_axis_units = dG%y_axis_units
    isg = dG%isg ; ieg = dG%ieg ; jsg = dG%jsg ; jeg = dG%jeg
    IsgB = dG%IsgB ; IegB = dG%IegB ; JsgB = dG%JsgB ; JegB = dG%JegB
  endif

  one_file = .true.
  if (domain_set) one_file = (thread == SINGLE_FILE)

  if (one_file) then
    if (domain_set) then
      call IO_handle%open(filename, action=OVERWRITE_FILE, &
          MOM_domain=domain, threading=thread, fileset=SINGLE_FILE)
    else
      call IO_handle%open(filename, action=OVERWRITE_FILE, threading=thread, &
          fileset=SINGLE_FILE)
    endif
  else
    call IO_handle%open(filename, action=OVERWRITE_FILE, MOM_domain=Domain, &
        threading=thread, fileset=thread)
  endif

! Define the coordinates.
  do k=1,novars
    position = vars(k)%position
    if (position == -1) position = position_from_horgrid(vars(k)%hor_grid)
    select case (position)
      case (CENTER) ; use_lath = .true. ; use_lonh = .true.
      case (CORNER) ; use_latq = .true. ; use_lonq = .true.
      case (EAST_FACE) ; use_lath = .true. ; use_lonq = .true.
      case (NORTH_FACE) ; use_latq = .true. ; use_lonh = .true.
      case (0) ! Do nothing.
      case default
        call MOM_error(WARNING, "MOM_io create_file: "//trim(vars(k)%name)//" has an unrecognized value of postion")
    end select
    select case (vars(k)%z_grid)
      case ('L') ; use_layer = .true.
      case ('i') ; use_int = .true.
      case ('1') ! Do nothing.
      case default
        call MOM_error(FATAL, "MOM_io create_file: "//trim(vars(k)%name)//&
                        " has unrecognized z_grid "//trim(vars(k)%z_grid))
    end select
    t_grid = adjustl(vars(k)%t_grid)
    select case (t_grid(1:1))
      case ('s', 'a', 'm') ; use_time = .true.
      case ('p') ; use_periodic = .true.
        if (len_trim(t_grid(2:8)) <= 0) call MOM_error(FATAL, &
          "MOM_io create_file: No periodic axis length was specified in "//&
          trim(vars(k)%t_grid) // " in the periodic axes of variable "//&
          trim(vars(k)%name)//" in file "//trim(filename))
        var_periods = -9999999
        t_grid_read = adjustl(t_grid(2:8))
        read(t_grid_read,*) var_periods
        if (var_periods == -9999999) call MOM_error(FATAL, &
          "MOM_io create_file: Failed to read the number of periods from "//&
          trim(vars(k)%t_grid) // " in the periodic axes of variable "//&
          trim(vars(k)%name)//" in file "//trim(filename))
        if (var_periods < 1) call MOM_error(FATAL, "MOM_io create_file: "//&
           "variable "//trim(vars(k)%name)//" in file "//trim(filename)//&
           " uses a periodic time axis, and must have a positive "//&
           "value for the number of periods in "//vars(k)%t_grid )
        if ((num_periods > 0) .and. (var_periods /= num_periods)) &
          call MOM_error(FATAL, "MOM_io create_file: "//&
            "Only one value of the number of periods can be used in the "//&
            "create_file call for file "//trim(filename)//".  The second is "//&
            "variable "//trim(vars(k)%name)//" with t_grid "//vars(k)%t_grid )

        num_periods = var_periods
      case ('1') ! Do nothing.
      case default
        call MOM_error(WARNING, "MOM_io create_file: "//trim(vars(k)%name)//&
                        " has unrecognized t_grid "//trim(vars(k)%t_grid))
    end select

    do n=1,5 ; if (len_trim(vars(k)%dim_names(n)) > 0) then
      dim_found = .false.
      do m=1,num_extra_dims
        if (lowercase(trim(vars(k)%dim_names(n))) == lowercase(trim(extra_axes(m)%name))) then
          use_extra_axis(m) = .true.
          dim_found = .true.
          exit
        endif
      enddo
      if (.not.dim_found) call MOM_error(FATAL, "Unable to find a match for dimension "//&
          trim(vars(k)%dim_names(n))//" for variable "//trim(vars(k)%name)//" in file "//trim(filename))
    endif ; enddo
  enddo

  if ((use_lath .or. use_lonh .or. use_latq .or. use_lonq)) then
    if (.not.domain_set) call MOM_error(FATAL, "create_file: "//&
      "An ocean_grid_type or dyn_horgrid_type is required to create a file with a horizontal coordinate.")

    call get_domain_components(Domain, x_domain, y_domain)
  endif
  if ((use_layer .or. use_int) .and. .not.present(GV)) call MOM_error(FATAL, &
    "create_file: A vertical grid type is required to create a file with a vertical coordinate.")

  if (use_lath) &
    axis_lath = IO_handle%register_axis("lath", units=y_axis_units, longname="Latitude", &
                        cartesian='Y', domain=y_domain, data=gridLatT(jsg:jeg))
  if (use_lonh) &
    axis_lonh = IO_handle%register_axis("lonh", units=x_axis_units, longname="Longitude", &
                        cartesian='X', domain=x_domain, data=gridLonT(isg:ieg))
  if (use_latq) &
    axis_latq = IO_handle%register_axis("latq", units=y_axis_units, longname="Latitude", &
                        cartesian='Y', domain=y_domain, data=gridLatB(JsgB:JegB), edge_axis=.true.)
  if (use_lonq) &
    axis_lonq = IO_handle%register_axis("lonq", units=x_axis_units, longname="Longitude", &
                        cartesian='X', domain=x_domain, data=gridLonB(IsgB:IegB), edge_axis=.true.)
  if (use_layer) &
    axis_layer = IO_handle%register_axis("Layer", units=trim(GV%zAxisUnits), &
                        longname="Layer "//trim(GV%zAxisLongName), cartesian='Z', &
                        sense=1, data=GV%sLayer(1:GV%ke))
  if (use_int) &
    axis_int = IO_handle%register_axis("Interface", units=trim(GV%zAxisUnits), &
                        longname="Interface "//trim(GV%zAxisLongName), cartesian='Z', &
                        sense=1, data=GV%sInterface(1:GV%ke+1))

  if (use_time) then ; if (present(timeunit)) then
    ! Set appropriate units, depending on the value.
    if (timeunit < 0.0) then
      time_units = "days" ! The default value.
    elseif ((timeunit >= 0.99) .and. (timeunit < 1.01)) then
      time_units = "seconds"
    elseif ((timeunit >= 3599.0) .and. (timeunit < 3601.0)) then
      time_units = "hours"
    elseif ((timeunit >= 86399.0) .and. (timeunit < 86401.0)) then
      time_units = "days"
    elseif ((timeunit >= 3.0e7) .and. (timeunit < 3.2e7)) then
      time_units = "years"
    else
      write(time_units,'(es8.2," s")') timeunit
    endif

    axis_time = IO_handle%register_axis("Time", units=time_units, longname="Time", cartesian='T')
  else
    axis_time = IO_handle%register_axis("Time", units="days", longname="Time", cartesian='T')
  endif ; endif

  if (use_periodic) then
    if (num_periods <= 1) call MOM_error(FATAL, "MOM_io create_file: "//&
      "num_periods for file "//trim(filename)//" must be at least 1.")
    ! Define a periodic axis with unit labels.
    allocate(axis_val(num_periods))
    do k=1,num_periods ; axis_val(k) = real(k) ; enddo
    axis_periodic = IO_handle%register_axis("Period", units="nondimensional", &
        longname="Periods for cyclical variables", cartesian='T', data=axis_val)
    deallocate(axis_val)
  endif

  do m=1,num_extra_dims ; if (use_extra_axis(m)) then
    if (allocated(extra_axes(m)%ax_data)) then
      more_axes(m) = IO_handle%register_axis(extra_axes(m)%name, units=extra_axes(m)%units, &
                          longname=extra_axes(m)%longname, cartesian=extra_axes(m)%cartesian, &
                          sense=extra_axes(m)%sense, data=extra_axes(m)%ax_data)
    elseif (trim(extra_axes(m)%cartesian) == "T") then
      more_axes(m) = IO_handle%register_axis(extra_axes(m)%name, units=extra_axes(m)%units, &
                          longname=extra_axes(m)%longname, cartesian=extra_axes(m)%cartesian)
    else
      ! FMS requires that non-time axes have variables that label their values, even if they are trivial.
      allocate (axis_val(extra_axes(m)%ax_size))
      do k=1,extra_axes(m)%ax_size ; axis_val(k) = real(k) ; enddo
      more_axes(m) = IO_handle%register_axis(extra_axes(m)%name, units=extra_axes(m)%units, &
                          longname=extra_axes(m)%longname, cartesian=extra_axes(m)%cartesian, &
                          sense=extra_axes(m)%sense, data=axis_val)
      deallocate(axis_val)
    endif
  endif ; enddo

  do k=1,novars
    numaxes = 0
    position = vars(k)%position
    if (position == -1) position = position_from_horgrid(vars(k)%hor_grid)
    select case (position)
      case (CENTER)
        numaxes = 2 ; axes(1) = axis_lonh ; axes(2) = axis_lath ; ax_name(1) = "lonh" ; ax_name(2) = "lath"
      case (CORNER)
        numaxes = 2 ; axes(1) = axis_lonq ; axes(2) = axis_latq ; ax_name(1) = "lonq" ; ax_name(2) = "latq"
      case (EAST_FACE)
        numaxes = 2 ; axes(1) = axis_lonq ; axes(2) = axis_lath ; ax_name(1) = "lonq" ; ax_name(2) = "lath"
      case (NORTH_FACE)
        numaxes = 2 ; axes(1) = axis_lonh ; axes(2) = axis_latq ; ax_name(1) = "lonh" ; ax_name(2) = "latq"
      case (0) ! Do nothing.
      case default
        call MOM_error(WARNING, "MOM_io create_file: "//trim(vars(k)%name)//&
                        " has unrecognized position, hor_grid = "//trim(vars(k)%hor_grid))
    end select
    select case (vars(k)%z_grid)
      case ('L') ; numaxes = numaxes+1 ; axes(numaxes) = axis_layer ; ax_name(numaxes) = "Layer"
      case ('i') ; numaxes = numaxes+1 ; axes(numaxes) = axis_int ; ax_name(numaxes) = "Interface"
      case ('1') ! Do nothing.
      case default
        call MOM_error(FATAL, "MOM_io create_file: "//trim(vars(k)%name)//&
                        " has unrecognized z_grid "//trim(vars(k)%z_grid))
    end select

    do n=1,numaxes
      if ( (len_trim(vars(k)%dim_names(n)) > 0) .and. (trim(ax_name(n)) /= trim(vars(k)%dim_names(n))) ) &
        call MOM_error(WARNING, "MOM_io create_file: dimension "//trim(ax_name(n))//&
               " of variable "//trim(vars(k)%name)//" in "//trim(filename)//&
               " is being set inconsistently as "//trim(vars(k)%dim_names(n)))
    enddo
    do n=numaxes+1,5 ; if (len_trim(vars(k)%dim_names(n)) > 0) then
      dim_found = .false.
      do m=1,num_extra_dims
        if (lowercase(trim(vars(k)%dim_names(n))) == lowercase(trim(extra_axes(m)%name))) then
          numaxes = numaxes+1 ; axes(numaxes) = more_axes(m)
          exit
        endif
      enddo
    endif ; enddo

    t_grid = adjustl(vars(k)%t_grid)
    select case (t_grid(1:1))
      case ('s', 'a', 'm') ; numaxes = numaxes+1 ; axes(numaxes) = axis_time
      case ('p')           ; numaxes = numaxes+1 ; axes(numaxes) = axis_periodic
      case ('1') ! Do nothing.
      case default
        call MOM_error(WARNING, "MOM_io create_file: "//trim(vars(k)%name)//&
                        " has unrecognized t_grid "//trim(vars(k)%t_grid))
    end select

    pack = 1
    if (present(checksums)) then
      fields(k) = IO_handle%register_field(axes(1:numaxes), vars(k)%name, vars(k)%units, &
                          vars(k)%longname, pack=pack, checksum=checksums(k,:), conversion=vars(k)%conversion)
    else
      fields(k) = IO_handle%register_field(axes(1:numaxes), vars(k)%name, vars(k)%units, &
                          vars(k)%longname, pack=pack, conversion=vars(k)%conversion)
    endif
  enddo

  if (present(global_atts)) then
    do n=1,size(global_atts)
      if (allocated(global_atts(n)%name) .and. allocated(global_atts(n)%att_val)) &
        call IO_handle%write_attribute(global_atts(n)%name, global_atts(n)%att_val)
    enddo
  endif

  ! Now write the variables with the axis label values
  if (use_lath) call IO_handle%write_field(axis_lath)
  if (use_latq) call IO_handle%write_field(axis_latq)
  if (use_lonh) call IO_handle%write_field(axis_lonh)
  if (use_lonq) call IO_handle%write_field(axis_lonq)
  if (use_layer) call IO_handle%write_field(axis_layer)
  if (use_int) call IO_handle%write_field(axis_int)
  if (use_periodic) call IO_handle%write_field(axis_periodic)
  do m=1,num_extra_dims ; if (use_extra_axis(m)) then
    call IO_handle%write_field(more_axes(m))
  endif ; enddo

  if (num_extra_dims > 0) then
    deallocate(use_extra_axis, more_axes)
  endif
end procedure create_MOM_file
module procedure reopen_file
  type(MOM_infra_file) :: mfile
  type(MOM_field), allocatable :: mfields(:)
  integer :: i
  mfile%handle_infra = IO_handle
  allocate(mfields(size(fields)))

  call reopen_MOM_file(mfile, filename, vars, novars, mfields, &
      threading=threading, timeunit=timeunit, G=G, dG=dG, GV=GV, &
      extra_axes=extra_axes, global_atts=global_atts)

  IO_handle = mfile%handle_infra
  call get_file_fields(IO_handle, fields)
end procedure reopen_file
module procedure reopen_MOM_file
  type(MOM_domain_type), pointer :: Domain => NULL()
  character(len=200) :: check_name, mesg
  integer :: length, nvar, thread
  logical :: exists, one_file, domain_set
  thread = SINGLE_FILE
  if (PRESENT(threading)) thread = threading

  ! For single-file IO, only the root PE is required to set up the fields.
  ! This permits calls by either the root PE or all PEs
  if (.not. is_root_PE() .and. thread == SINGLE_FILE) return

  ! For multiple IO domains, we would need additional functionality:
  ! * Identify ranks as IO PEs
  ! * Determine the filename of
  ! Neither of these tasks should be handed by MOM6, so we cannot safely use
  ! this function.  A framework-specific `inquire()` function is needed.
  ! Until it exists, we will disable this function.
  if (thread == MULTIPLE) &
    call MOM_error(FATAL, 'reopen_MOM_file does not yet support files with ' &
      // 'multiple I/O domains.')

  check_name = filename
  length = len(trim(check_name))
  if (check_name(length-2:length) /= ".nc") check_name = trim(check_name)//".nc"
  if (thread /= SINGLE_FILE) check_name = trim(check_name)//".0000"

  inquire(file=check_name,EXIST=exists)

  if (.not.exists) then
    call create_MOM_file(IO_handle, filename, vars, novars, fields, &
        threading, timeunit, G=G, dG=dG, GV=GV, extra_axes=extra_axes, &
        global_atts=global_atts)
  else

    domain_set = .false.
    if (present(G)) then
      domain_set = .true. ; Domain => G%Domain
    elseif (present(dG)) then
      domain_set = .true. ; Domain => dG%Domain
    endif

    one_file = .true.
    if (domain_set) one_file = (thread == SINGLE_FILE)

    if (one_file) then
      call IO_handle%open(filename, APPEND_FILE, threading=thread)
    else
      call IO_handle%open(filename, APPEND_FILE, MOM_domain=Domain)
    endif
    if (.not. IO_handle%file_is_open()) return

    call IO_handle%get_file_info(nvar=nvar)

    if (nvar == -1) then
      write (mesg,*) "Reopening file ",trim(filename)," apparently had ",nvar,&
                     " variables. Clobbering and creating file with ",novars," instead."
      call MOM_error(WARNING,"MOM_io: "//mesg)
      call create_MOM_file(IO_handle, filename, vars, novars, fields, &
          threading, timeunit, G=G, dG=dG, GV=GV, extra_axes=extra_axes, &
          global_atts=global_atts)
    elseif (nvar /= novars) then
      write (mesg,*) "Reopening file ",trim(filename)," with ",novars,&
                     " variables instead of ",nvar,"."
      call MOM_error(FATAL,"MOM_io: "//mesg)
    endif

    if (nvar > 0) call IO_handle%get_file_fields(fields(1:nvar))
  endif
end procedure reopen_MOM_file
module procedure stdout_if_root
  stdout_if_root = 0
  if (is_root_PE()) stdout_if_root = stdout
end procedure stdout_if_root
module procedure num_timelevels
  character(len=256) :: msg
  integer :: ndims
  integer :: sizes(8)
  n_time = -1

  ! To do almost the same via MOM_io_infra calls, we could do the following:
  !   found = field_exists(filename, varname)
  !   if (found) then
  !     call open_file(ncid, filename, action=READONLY_FILE, form=NETCDF_FILE, threading=MULTIPLE)
  !     call get_file_info(ncid, ntime=n_time)
  !   endif
  ! However, this does not handle the case where the time axis for the variable is not the record
  ! axis and min_dims is not used.

  call get_var_sizes(filename, varname, ndims, sizes, match_case=.false., caller="num_timelevels")

  if (ndims > 0) n_time = sizes(ndims)

  if (present(min_dims)) then
    if (ndims < min_dims-1) then
      write(msg, '(I0)') min_dims
      call MOM_error(WARNING, "num_timelevels: variable "//trim(varname)//" in file "//&
          trim(filename)//" has fewer than min_dims = "//trim(msg)//" dimensions.")
      n_time = -1
    elseif (ndims == min_dims - 1) then
      n_time = 0
    endif
  endif

end procedure num_timelevels
module procedure get_var_sizes
  logical :: do_read, do_broadcast
  integer, allocatable :: size_msg(:)  ! An array combining the number of dimensions and the sizes.
  integer :: n, nval
  do_read = is_root_pe()
  if (present(all_read)) do_read = all_read .or. do_read
  do_broadcast = .true. ; if (present(all_read)) do_broadcast = .not.all_read

  if (do_read) call read_var_sizes(filename, varname, ndims, sizes, match_case, caller, dim_names, ncid_in)

  if (do_broadcast) then
    ! Distribute the sizes from the root PE.
    nval = size(sizes) + 1

    allocate(size_msg(nval))
    size_msg(1) = ndims
    do n=2,nval ; size_msg(n) = sizes(n-1) ; enddo

    call broadcast(size_msg, nval, blocking=.true.)

    ndims = size_msg(1)
    do n=2,nval ;  sizes(n-1) = size_msg(n) ; enddo
    deallocate(size_msg)

    if (present(dim_names) .and. (ndims > 0)) then
      nval = min(ndims, size(dim_names))
      call broadcast(dim_names(1:nval), len(dim_names(1)), blocking=.true.)
    endif
  endif

end procedure get_var_sizes
module procedure read_var_sizes
  character(len=256) :: hdr, dimname
  integer, allocatable :: dimids(:)
  integer :: varid, ncid, n, status
  logical :: success, found
  hdr = "get_var_size: " ; if (present(caller)) hdr = trim(hdr)//": "
  sizes(:) = 0 ; ndims = -1

  if (present(ncid_in)) then
    ncid = ncid_in
  else
    call open_file_to_read(filename, ncid, success=success)
    if (.not.success) then
      call MOM_error(WARNING, "Unsuccessfully attempted to open file "//trim(filename))
      return
    endif
  endif

  ! Get the dimension sizes of the variable varname.
  call get_varid(varname, ncid, filename, varid, match_case=match_case, found=found)
  if (.not.found) then
    call MOM_error(WARNING, "Could not find variable "//trim(varname)//" in file "//trim(filename))
    return
  endif

  status = NF90_inquire_variable(ncid, varid, ndims=ndims)
  if (status /= NF90_NOERR) then
    call MOM_error(WARNING, trim(hdr) // trim(NF90_STRERROR(status)) //&
      " Getting number of dimensions of "//trim(varname)//" in "//trim(filename))
    return
  endif
  if (ndims < 1) return

  allocate(dimids(ndims))
  status = NF90_inquire_variable(ncid, varid, dimids=dimids(1:ndims))
  if (status /= NF90_NOERR) then
    call MOM_error(WARNING, trim(hdr) // trim(NF90_STRERROR(status)) //&
      " Getting dimension IDs for "//trim(varname)//" in "//trim(filename))
    deallocate(dimids) ; return
  endif

  do n = 1, min(ndims,size(sizes))
    status = NF90_Inquire_Dimension(ncid, dimids(n), name=dimname, len=sizes(n))
    if (status /= NF90_NOERR) call MOM_error(WARNING, trim(hdr) // trim(NF90_STRERROR(status)) //&
        " Getting dimension length for "//trim(varname)//" in "//trim(filename))
    if (present(dim_names)) then
      if (n <= size(dim_names)) dim_names(n) = trim(dimname)
    endif
  enddo
  deallocate(dimids)

  if (.not.present(ncid_in)) call close_file_to_read(ncid, filename)

end procedure read_var_sizes
module procedure read_variable_0d
  integer :: varid, ncid, rc
  character(len=256) :: hdr
  hdr = "read_variable_0d"

  if (is_root_pe()) then
    if (present(ncid_in)) then
      ncid = ncid_in
    else
      call open_file_to_read(filename, ncid)
    endif

    call get_varid(varname, ncid, filename, varid, match_case=.false.)
    if (varid < 0) call MOM_error(FATAL, "Unable to get netCDF varid for "//trim(varname)//&
                                         " in "//trim(filename))
    rc = NF90_get_var(ncid, varid, var)
    if (rc /= NF90_NOERR) call MOM_error(FATAL, trim(hdr) // trim(NF90_STRERROR(rc)) //&
          " Difficulties reading "//trim(varname)//" from "//trim(filename))

    if (.not.present(ncid_in)) call close_file_to_read(ncid, filename)

    if (present(scale)) var = scale * var
  endif

  call broadcast(var, blocking=.true.)
end procedure read_variable_0d
module procedure read_variable_1d
  integer :: varid, ncid, rc
  character(len=256) :: hdr
  hdr = "read_variable_1d"

  if (is_root_pe()) then
    if (present(ncid_in)) then
      ncid = ncid_in
    else
      call open_file_to_read(filename, ncid)
    endif

    call get_varid(varname, ncid, filename, varid, match_case=.false.)
    if (varid < 0) call MOM_error(FATAL, "Unable to get netCDF varid for "//trim(varname)//&
                                         " in "//trim(filename))
    rc = NF90_get_var(ncid, varid, var)
    if (rc /= NF90_NOERR) call MOM_error(FATAL, trim(hdr) // trim(NF90_STRERROR(rc)) //&
          " Difficulties reading "//trim(varname)//" from "//trim(filename))

    if (.not.present(ncid_in)) call close_file_to_read(ncid, filename)

    if (present(scale)) then ; if (scale /= 1.0) then
      var(:) = scale * var(:)
    endif ; endif
  endif

  call broadcast(var, size(var), blocking=.true.)
end procedure read_variable_1d
module procedure read_variable_0d_int
  integer :: varid, ncid, rc
  character(len=256) :: hdr
  hdr = "read_variable_0d_int"

  if (is_root_pe()) then
    if (present(ncid_in)) then
      ncid = ncid_in
    else
      call open_file_to_read(filename, ncid)
    endif

    call get_varid(varname, ncid, filename, varid, match_case=.false.)
    if (varid < 0) call MOM_error(FATAL, "Unable to get netCDF varid for "//trim(varname)//&
                                         " in "//trim(filename))
    rc = NF90_get_var(ncid, varid, var)
    if (rc /= NF90_NOERR) call MOM_error(FATAL, trim(hdr) // trim(NF90_STRERROR(rc)) //&
          " Difficulties reading "//trim(varname)//" from "//trim(filename))

    if (.not.present(ncid_in)) call close_file_to_read(ncid, filename)
  endif

  call broadcast(var, blocking=.true.)
end procedure read_variable_0d_int
module procedure read_variable_1d_int
  integer :: varid, ncid, rc
  character(len=256) :: hdr
  hdr = "read_variable_1d_int"

  if (is_root_pe()) then
    if (present(ncid_in)) then
      ncid = ncid_in
    else
      call open_file_to_read(filename, ncid)
    endif

    call get_varid(varname, ncid, filename, varid, match_case=.false.)
    if (varid < 0) call MOM_error(FATAL, "Unable to get netCDF varid for "//trim(varname)//&
                                         " in "//trim(filename))
    rc = NF90_get_var(ncid, varid, var)
    if (rc /= NF90_NOERR) call MOM_error(FATAL, trim(hdr) // trim(NF90_STRERROR(rc)) //&
          " Difficulties reading "//trim(varname)//" from "//trim(filename))

    if (.not.present(ncid_in)) call close_file_to_read(ncid, filename)
  endif

  call broadcast(var, size(var), blocking=.true.)
end procedure read_variable_1d_int
module procedure read_variable_2d
  integer :: ncid, varid
  integer :: field_ndims, dim_len
  integer, allocatable :: field_dimids(:), field_shape(:)
  integer, allocatable :: field_start(:), field_nread(:)
  integer :: i, rc
  character(len=*), parameter :: hdr = "read_variable_2d: "
  if (present(start)) then
    if (size(start) < 2) &
      call MOM_error(FATAL, hdr // trim(nf90_strerror(rc)) &
        // " start must have at least two dimensions.")
  endif

  if (present(nread)) then
    if (size(nread) < 2) &
      call MOM_error(FATAL, hdr // trim(nf90_strerror(rc)) &
        // " nread must have at least two dimensions.")

    if (any(nread(3:) > 1)) &
      call MOM_error(FATAL, hdr // trim(nf90_strerror(rc)) &
        // " nread may only read a single level in higher dimensions.")
  endif

  ! Since start and nread may be reshaped, we cannot rely on netCDF to ensure
  ! that their lengths are equivalent, and must do it here.
  if (present(start) .and. present(nread)) then
    if (size(start) /= size(nread)) &
      call MOM_error(FATAL, hdr // trim(nf90_strerror(rc)) &
        // " start and nread must have the same length.")
  endif

  ! Open and read `varname` from `filename`
  if (is_root_pe()) then
    if (present(ncid_in)) then
      ncid = ncid_in
    else
      call open_file_to_Read(filename, ncid)
    endif

    call get_varid(varname, ncid, filename, varid, match_case=.false.)
    if (varid < 0) call MOM_error(FATAL, "Unable to get netCDF varid for "//trim(varname)//&
                                         " in "//trim(filename))

    ! Query for the dimensionality of the input field
    rc = nf90_inquire_variable(ncid, varid, ndims=field_ndims)
    if (rc /= NF90_NOERR) call MOM_error(FATAL, hdr // trim(nf90_strerror(rc)) //&
          ": Difficulties reading "//trim(varname)//" from "//trim(filename))

    ! Confirm that field is at least 2d
    if (field_ndims < 2) &
      call MOM_error(FATAL, hdr // trim(nf90_strerror(rc)) // " " // &
          trim(varname) // " from " // trim(filename) // " is not a 2D field.")

    ! If start and nread are present, then reshape them to match field dims
    if (present(start) .or. present(nread)) then
      allocate(field_shape(field_ndims))
      allocate(field_dimids(field_ndims))

      rc = nf90_inquire_variable(ncid, varid, dimids=field_dimids)
      if (rc /= NF90_NOERR) call MOM_error(FATAL, hdr // trim(nf90_strerror(rc)) //&
            ": Difficulties reading "//trim(varname)//" from "//trim(filename))

      do i = 1, field_ndims
        rc = nf90_inquire_dimension(ncid, field_dimids(i), len=dim_len)
        if (rc /= NF90_NOERR) &
          call MOM_error(FATAL, hdr // trim(nf90_strerror(rc)) &
              // ": Difficulties reading dimensions from " // trim(filename))
        field_shape(i) = dim_len
      enddo

      ! Reshape start(:) and nreads(:) in case ranks differ
      allocate(field_start(field_ndims))
      field_start(:) = 1
      if (present(start)) then
        dim_len = min(size(start), size(field_start))
        field_start(:dim_len) = start(:dim_len)
      endif

      allocate(field_nread(field_ndims))
      field_nread(:2) = field_shape(:2)
      field_nread(3:) = 1
      if (present(nread)) field_nread(:2) = nread(:2)

      rc = nf90_get_var(ncid, varid, var, field_start, field_nread)

      deallocate(field_start)
      deallocate(field_nread)
      deallocate(field_shape)
      deallocate(field_dimids)
    else
      rc = nf90_get_var(ncid, varid, var)
    endif

    if (rc /= NF90_NOERR) call MOM_error(FATAL, hdr // trim(nf90_strerror(rc)) //&
          " Difficulties reading "//trim(varname)//" from "//trim(filename))

    if (.not.present(ncid_in)) call close_file_to_read(ncid, filename)
  endif

  call broadcast(var, size(var), blocking=.true.)
end procedure read_variable_2d
module procedure read_variable_3d
  integer :: ncid, varid
  integer :: field_ndims, dim_len
  integer, allocatable :: field_dimids(:), field_shape(:)
  integer, allocatable :: field_start(:), field_nread(:)
  integer :: i, rc
  character(len=*), parameter :: hdr = "read_variable_3d: "
  if (present(start)) then
    if (size(start) < 2) &
      call MOM_error(FATAL, hdr // trim(nf90_strerror(rc)) &
        // " start must have at least two dimensions.")
  endif

  if (present(nread)) then
    if (size(nread) < 2) &
      call MOM_error(FATAL, hdr // trim(nf90_strerror(rc)) &
        // " nread must have at least two dimensions.")

    if (any(nread(3:) > 1)) &
      call MOM_error(FATAL, hdr // trim(nf90_strerror(rc)) &
        // " nread may only read a single level in higher dimensions.")
  endif

  ! Since start and nread may be reshaped, we cannot rely on netCDF to ensure
  ! that their lengths are equivalent, and must do it here.
  if (present(start) .and. present(nread)) then
    if (size(start) /= size(nread)) &
      call MOM_error(FATAL, hdr // trim(nf90_strerror(rc)) &
        // " start and nread must have the same length.")
  endif

  ! Open and read `varname` from `filename`
  if (is_root_pe()) then
    if (present(ncid_in)) then
      ncid = ncid_in
    else
      call open_file_to_Read(filename, ncid)
    endif

    call get_varid(varname, ncid, filename, varid, match_case=.false.)
    if (varid < 0) call MOM_error(FATAL, "Unable to get netCDF varid for "//trim(varname)//&
                                         " in "//trim(filename))

    ! Query for the dimensionality of the input field
    rc = nf90_inquire_variable(ncid, varid, ndims=field_ndims)
    if (rc /= NF90_NOERR) call MOM_error(FATAL, hdr // trim(nf90_strerror(rc)) //&
          ": Difficulties reading "//trim(varname)//" from "//trim(filename))

    ! Confirm that field is at least 2d
    if (field_ndims < 2) &
      call MOM_error(FATAL, hdr // trim(nf90_strerror(rc)) // " " // &
          trim(varname) // " from " // trim(filename) // " is not a 2D field.")

    ! If start and nread are present, then reshape them to match field dims
    if (present(start) .or. present(nread)) then
      allocate(field_shape(field_ndims))
      allocate(field_dimids(field_ndims))

      rc = nf90_inquire_variable(ncid, varid, dimids=field_dimids)
      if (rc /= NF90_NOERR) call MOM_error(FATAL, hdr // trim(nf90_strerror(rc)) //&
            ": Difficulties reading "//trim(varname)//" from "//trim(filename))

      do i = 1, field_ndims
        rc = nf90_inquire_dimension(ncid, field_dimids(i), len=dim_len)
        if (rc /= NF90_NOERR) &
          call MOM_error(FATAL, hdr // trim(nf90_strerror(rc)) &
              // ": Difficulties reading dimensions from " // trim(filename))
        field_shape(i) = dim_len
      enddo

      ! Reshape start(:) and nreads(:) in case ranks differ
      allocate(field_start(field_ndims))
      field_start(:) = 1
      if (present(start)) then
        dim_len = min(size(start), size(field_start))
        field_start(:dim_len) = start(:dim_len)
      endif

      allocate(field_nread(field_ndims))
      field_nread(:3) = field_shape(:3)
      !field_nread(3:) = 1
      if (present(nread)) field_nread(:3) = nread(:3)

      rc = nf90_get_var(ncid, varid, var, field_start, field_nread)

      deallocate(field_start)
      deallocate(field_nread)
      deallocate(field_shape)
      deallocate(field_dimids)
    else
      rc = nf90_get_var(ncid, varid, var)
    endif

    if (rc /= NF90_NOERR) call MOM_error(FATAL, hdr // trim(nf90_strerror(rc)) //&
          " Difficulties reading "//trim(varname)//" from "//trim(filename))

    if (.not.present(ncid_in)) call close_file_to_read(ncid, filename)
  endif

  call broadcast(var, size(var), blocking=.true.)
end procedure read_variable_3d
module procedure read_attribute_str
  logical :: do_read, do_broadcast
  integer :: rc, ncid, varid, att_type, att_len, info(2)
  character(len=256) :: hdr, att_str
  character(len=:), dimension(:), allocatable :: tmp_str
  hdr = "read_attribute_str"
  att_len = 0

  do_read = is_root_pe() ; if (present(all_read)) do_read = all_read .or. do_read
  do_broadcast = .true. ; if (present(all_read)) do_broadcast = .not.all_read

  if (do_read) then
    if (present(ncid_in)) then
      ncid = ncid_in
    else
      call open_file_to_read(filename, ncid, success=found)
      if (present(found)) then ; if (.not.found) do_read = .false. ; endif
    endif
  endif

  if (do_read) then
    rc = NF90_ENOTATT ; att_len = 0
    if (present(varname)) then  ! Read a variable attribute
      call get_varid(varname, ncid, filename, varid, match_case=.false., found=found)
      att_str = "att "//trim(attname)//" for "//trim(varname)//" from "//trim(filename)
    else   ! Read a global attribute
      varid = NF90_GLOBAL
      att_str = "global att "//trim(attname)//" from "//trim(filename)
    endif
    if ((varid > 0) .or. (varid == NF90_GLOBAL)) then ! The named variable does exist, and found would be true.
      rc = NF90_inquire_attribute(ncid, varid, attname, xtype=att_type, len=att_len)
      if ((.not. present(found)) .or. (rc /= NF90_ENOTATT)) then
        if ((rc /= NF90_NOERR) .and. (rc /= NF90_ENOTATT)) &
          call MOM_error(FATAL, trim(hdr) // trim(NF90_STRERROR(rc)) //" Error getting info for "//trim(att_str))
        if (att_type /= NF90_CHAR) &
          call MOM_error(FATAL, trim(hdr)//": Attribute data type is not a char for "//trim(att_str))
  !      if (att_len > len(att_val)) &
  !        call MOM_error(FATAL, trim(hdr)//": Insufficiently long string passed in to read "//trim(att_str))
        allocate(character(att_len) :: att_val)

        if (rc == NF90_NOERR) then
          rc = NF90_get_att(ncid, varid, attname, att_val)
          if ((rc /= NF90_NOERR) .and. (rc /= NF90_ENOTATT)) &
            call MOM_error(FATAL, trim(hdr) // trim(NF90_STRERROR(rc)) //" Difficulties reading "//trim(att_str))
        endif
      endif
    endif
    if (present(found)) found = (rc == NF90_NOERR)

    if (.not.present(ncid_in)) call close_file_to_read(ncid, filename)
  endif

  if (do_broadcast) then
    ! Communicate the string length
    info(1) = att_len ; info(2) = 0 ; if (do_read .and. found) info(2) = 1
    call broadcast(info, 2, blocking=.true.)
    if (present(found)) then
      found = (info(2) /= 0)
      if (.not. found) return
    endif
    att_len = info(1)

    if (att_len > 0) then
      ! These extra copies are here because broadcast only supports arrays of strings.
      allocate(character(att_len) :: tmp_str(1))
      if (.not.do_read) allocate(character(att_len) :: att_val)
      if (do_read) tmp_str(1) = att_val
      call broadcast(tmp_str, att_len, blocking=.true.)
      att_val = tmp_str(1)
    elseif (.not.allocated(att_val)) then
      allocate(character(4) :: att_val) ; att_val = ''
    endif
  elseif (.not.allocated(att_val)) then
    allocate(character(4) :: att_val) ; att_val = ''
  endif
end procedure read_attribute_str
module procedure read_attribute_int32
  logical :: do_read, do_broadcast
  integer :: rc, ncid, varid, is_found
  character(len=256) :: hdr
  hdr = "read_attribute_int32"
  att_val = 0

  do_read = is_root_pe() ; if (present(all_read)) do_read = all_read .or. do_read
  do_broadcast = .true. ; if (present(all_read)) do_broadcast = .not.all_read

  if (do_read) then
    if (present(ncid_in)) then
      ncid = ncid_in
    else
      call open_file_to_read(filename, ncid, success=found)
      if (present(found)) then ; if (.not.found) do_read = .false. ; endif
    endif
  endif

  if (do_read) then
    rc = NF90_ENOTATT
    if (present(varname)) then  ! Read a variable attribute
      call get_varid(varname, ncid, filename, varid, match_case=.false., found=found)
      if (varid >= 0) then ! The named variable does exist, and found would be true.
        rc = NF90_get_att(ncid, varid, attname, att_val)
        if ((rc /= NF90_NOERR) .and. (rc /= NF90_ENOTATT)) &
          call MOM_error(FATAL, trim(hdr) // trim(NF90_STRERROR(rc)) //" Difficulties reading att "//&
                trim(attname)//" for "//trim(varname)//" from "//trim(filename))
      endif
    else  ! Read a global attribute
      rc = NF90_get_att(ncid, NF90_GLOBAL, attname, att_val)
      if ((rc /= NF90_NOERR) .and. (rc /= NF90_ENOTATT)) &
        call MOM_error(FATAL, trim(hdr) // trim(NF90_STRERROR(rc)) //&
                " Difficulties reading global att "//trim(attname)//" from "//trim(filename))
    endif
    if (present(found)) found = (rc == NF90_NOERR)

    if (.not.present(ncid_in)) call close_file_to_read(ncid, filename)
  endif

  if (do_broadcast) then
    if (present(found)) then
      is_found = 0 ; if (is_root_pe() .and. found) is_found = 1
      call broadcast(is_found, blocking=.false.)
    endif
    call broadcast(att_val, blocking=.true.)
    if (present(found)) found = (is_found /= 0)
  endif

end procedure read_attribute_int32
module procedure read_attribute_int64
  logical :: do_read, do_broadcast
  integer :: rc, ncid, varid, is_found
  character(len=256) :: hdr
  hdr = "read_attribute_int64"
  att_val = 0

  do_read = is_root_pe() ; if (present(all_read)) do_read = all_read .or. do_read
  do_broadcast = .true. ; if (present(all_read)) do_broadcast = .not.all_read

  if (do_read) then
    if (present(ncid_in)) then
      ncid = ncid_in
    else
      call open_file_to_read(filename, ncid, success=found)
      if (present(found)) then ; if (.not.found) do_read = .false. ; endif
    endif
  endif

  if (do_read) then
    rc = NF90_ENOTATT
    if (present(varname)) then  ! Read a variable attribute
      call get_varid(varname, ncid, filename, varid, match_case=.false., found=found)
      if (varid >= 0) then ! The named variable does exist, and found would be true.
        rc = NF90_get_att(ncid, varid, attname, att_val)
        if ((rc /= NF90_NOERR) .and. (rc /= NF90_ENOTATT)) &
          call MOM_error(FATAL, trim(hdr) // trim(NF90_STRERROR(rc)) //" Difficulties reading att "//&
                trim(attname)//" for "//trim(varname)//" from "//trim(filename))
      endif
    else  ! Read a global attribute
      rc = NF90_get_att(ncid, NF90_GLOBAL, attname, att_val)
      if ((rc /= NF90_NOERR) .and. (rc /= NF90_ENOTATT)) &
        call MOM_error(FATAL, trim(hdr) // trim(NF90_STRERROR(rc)) //&
                " Difficulties reading global att "//trim(attname)//" from "//trim(filename))
    endif
    if (present(found)) found = (rc == NF90_NOERR)

    rc = NF90_close(ncid)
  endif

  if (do_broadcast) then
    if (present(found)) then
      is_found = 0 ; if (is_root_pe() .and. found) is_found = 1
      call broadcast(is_found, blocking=.false.)
    endif
    call broadcast(att_val, blocking=.true.)
    if (present(found)) found = (is_found /= 0)
  endif

end procedure read_attribute_int64
module procedure read_attribute_real
  logical :: do_read, do_broadcast
  integer :: rc, ncid, varid, is_found
  character(len=256) :: hdr
  hdr = "read_attribute_real"
  att_val = 0.0

  do_read = is_root_pe() ; if (present(all_read)) do_read = all_read .or. do_read
  do_broadcast = .true. ; if (present(all_read)) do_broadcast = .not.all_read

  if (do_read) then
    if (present(ncid_in)) then
      ncid = ncid_in
    else
      call open_file_to_read(filename, ncid, success=found)
      if (present(found)) then ; if (.not.found) do_read = .false. ; endif
    endif
  endif

  if (do_read) then
    rc = NF90_ENOTATT
    if (present(varname)) then  ! Read a variable attribute
      call get_varid(varname, ncid, filename, varid, match_case=.false., found=found)
      if (varid >= 0) then ! The named variable does exist, and found would be true.
        rc = NF90_get_att(ncid, varid, attname, att_val)
        if ((rc /= NF90_NOERR) .and. (rc /= NF90_ENOTATT)) &
          call MOM_error(FATAL, trim(hdr) // trim(NF90_STRERROR(rc)) //" Difficulties reading att "//&
                trim(attname)//" for "//trim(varname)//" from "//trim(filename))
      endif
    else  ! Read a global attribute
      rc = NF90_get_att(ncid, NF90_GLOBAL, attname, att_val)
      if ((rc /= NF90_NOERR) .and. (rc /= NF90_ENOTATT)) &
        call MOM_error(FATAL, trim(hdr) // trim(NF90_STRERROR(rc)) //&
                " Difficulties reading global att "//trim(attname)//" from "//trim(filename))
    endif
    if (present(found)) found = (rc == NF90_NOERR)

    if (.not.present(ncid_in)) call close_file_to_read(ncid, filename)
  endif

  if (do_broadcast) then
    if (present(found)) then
      is_found = 0 ; if (is_root_pe() .and. found) is_found = 1
      call broadcast(is_found, blocking=.false.)
    endif
    call broadcast(att_val, blocking=.true.)
    if (present(found)) found = (is_found /= 0)
  endif

end procedure read_attribute_real
module procedure open_file_to_read
  integer rc
  rc = NF90_open(trim(filename), NF90_NOWRITE, ncid)
  if (present(success)) then
    success = (rc == NF90_NOERR)
  elseif (rc /= NF90_NOERR) then
    call MOM_error(FATAL, "Difficulties opening "//trim(filename)//" - "//trim(NF90_STRERROR(rc)) )
  endif

end procedure open_file_to_read
module procedure close_file_to_read
  integer :: rc
  if (ncid >= 0) then
    rc = NF90_close(ncid)
    if (present(filename) .and. (rc /= NF90_NOERR)) then
      call MOM_error(WARNING, "Difficulties closing "//trim(filename)//": "//trim(NF90_STRERROR(rc)))
    elseif (rc /= NF90_NOERR) then
      call MOM_error(WARNING, "Difficulties closing file: "//trim(NF90_STRERROR(rc)))
    endif
  endif
  ncid = -1
end procedure close_file_to_read
module procedure get_varid
  logical :: var_found, insensitive
  character(len=256) :: name
  integer, allocatable :: varids(:)
  integer :: nvars, status, n
  varid = -1
  var_found = .false.
  insensitive = .false. ; if (present(match_case)) insensitive = .not.match_case

  if (insensitive) then
    ! This code ounddoes a case-insensitive search for a variable in the file.
    status = NF90_inquire(ncid, nVariables=nvars)
    if (present(found) .and. ((status /= NF90_NOERR) .or. (nvars < 1))) then
      found = .false. ; return
    elseif (status /= NF90_NOERR) then
      call MOM_error(FATAL, "get_varid:  Difficulties getting the number of variables in file "//&
          trim(filename)//" - "//trim(NF90_STRERROR(status)))
    elseif (nvars < 1) then
      call MOM_error(FATAL, "get_varid: There appear not to be any variables in "//trim(filename))
    endif

    allocate(varids(nvars))

    status = nf90_inq_varids(ncid, nvars, varids)
    if (status /= NF90_NOERR) then
      call MOM_error(WARNING, "get_varid: Difficulties getting the variable IDs in file "//&
          trim(filename)//" - "//trim(NF90_STRERROR(status)))
      nvars = -1  ! Full error handling will occur after the do-loop.
    endif

    do n = 1,nvars
      status = nf90_inquire_variable(ncid, varids(n), name=name)
      if (status /= NF90_NOERR) then
        call MOM_error(WARNING, "get_varid:  Difficulties getting a variable name in file "//&
            trim(filename)//" - "//trim(NF90_STRERROR(status)))
      endif

      if (trim(lowercase(name)) == trim(lowercase(varname))) then
        if (var_found) then
          call MOM_error(WARNING, "get_varid: Two variables match the case-insensitive name "//&
                  trim(varname)//" in file "//trim(filename))
          ! Replace the first variable if the second one is a case-sensitive match
          if (trim(name) == trim(varname)) varid = varids(n)
        else
          varid = varids(n) ; var_found = .true.
        endif
      endif
    enddo
    if (present(found)) found = var_found
    if ((.not.var_found) .and. .not.present(found)) call MOM_error(FATAL, &
        "get_varid: variable "//trim(varname)//" was not found in file "//trim(filename))

    deallocate(varids)
  else
    status = NF90_INQ_VARID(ncid, trim(varname), varid)
    if (present(found)) found = (status == NF90_NOERR)
    if ((status /= NF90_NOERR) .and. .not.present(found)) then
      call MOM_error(FATAL, "get_varid: Difficulties getting a variable id for "//&
          trim(varname)//" in file "//trim(filename)//" - "//trim(NF90_STRERROR(status)))
    endif
  endif

end procedure get_varid
module procedure verify_variable_units
  character (len=200) :: units
  logical :: units_correct, success
  integer :: i, ncid, status, vid
  if (.not.is_root_pe()) then ! Only the root PE should do the verification.
    ierr = .false. ; msg = '' ; return
  endif

  ierr = .true.
  call open_file_to_read(filename, ncid, success)
  if (.not.success) then
    msg = 'File not found: '//trim(filename)
    return
  endif

  status = NF90_INQ_VARID(ncid, trim(varname), vid)
  if (status /= NF90_NOERR) then
    msg = 'Var not found: '//trim(varname)
  else
    status = NF90_GET_ATT(ncid, vid, "units", units)
    if (status /= NF90_NOERR) then
      msg = 'Attribute not found: units'
    else
      ! NF90_GET_ATT can return attributes with null characters, which TRIM will not truncate.
      ! This loop replaces any null characters with a space so that the subsequent check
      ! between the read units and the expected units will pass
      do i=1,LEN_TRIM(units)
        if (units(i:i) == CHAR(0)) units(i:i) = " "
      enddo

      units_correct = (trim(units) == trim(expected_units))
      if (present(alt_units)) then
        units_correct = units_correct .or. (trim(units) == trim(alt_units))
      endif
      if (units_correct) then
        ierr = .false.
        msg = ''
      else
        msg = 'Units incorrect: '//trim(units)//' /= '//trim(expected_units)
      endif
    endif
  endif

  status = NF90_close(ncid)

end procedure verify_variable_units
module procedure var_desc
  character(len=120) :: cllr
  cllr = "var_desc"
  if (present(caller)) cllr = trim(caller)

  call safe_string_copy(name, vd%name, "vd%name", cllr)

  vd%longname = "" ; vd%units = ""
  vd%hor_grid = 'h' ; vd%position = CENTER ; vd%z_grid = 'L' ; vd%t_grid = 's'
  if (present(dim_names)) vd%z_grid = '1' ! In this case the names are used to set the non-horizontal axes
  if (present(fixed)) then ; if (fixed) vd%t_grid = '1' ; endif

  vd%cmor_field_name  =  ""
  vd%cmor_units       =  ""
  vd%cmor_longname    =  ""
  vd%conversion       =  1.0
  vd%dim_names(:)     =  ""

  call modify_vardesc(vd, units=units, longname=longname, hor_grid=hor_grid, &
                      z_grid=z_grid, t_grid=t_grid, position=position, dim_names=dim_names, &
                      cmor_field_name=cmor_field_name, cmor_units=cmor_units, &
                      cmor_longname=cmor_longname, conversion=conversion, caller=cllr, &
                      extra_axes=extra_axes)

end procedure var_desc
module procedure modify_vardesc
  character(len=120) :: cllr
  integer :: n
  cllr = "mod_vardesc" ; if (present(caller)) cllr = trim(caller)

  if (present(name))      call safe_string_copy(name, vd%name, "vd%name", cllr)

  if (present(longname))  call safe_string_copy(longname, vd%longname, &
                               "vd%longname of "//trim(vd%name), cllr)
  if (present(units))     call safe_string_copy(units, vd%units,       &
                               "vd%units of "//trim(vd%name), cllr)
  if (present(position)) then
    vd%position = position
    select case (position)
      case (CENTER)     ; vd%hor_grid = 'T'
      case (CORNER)     ; vd%hor_grid = 'Bu'
      case (EAST_FACE)  ; vd%hor_grid = 'Cu'
      case (NORTH_FACE) ; vd%hor_grid = 'Cv'
      case (0)          ; vd%hor_grid = '1'
      case default
        call MOM_error(FATAL, "modify_vardesc: "//trim(vd%name)//" has unrecognized position argument")
    end select
  endif
  if (present(hor_grid)) then
    call safe_string_copy(hor_grid, vd%hor_grid, "vd%hor_grid of "//trim(vd%name), cllr)
    vd%position = position_from_horgrid(vd%hor_grid)
    if (present(caller) .and. (vd%position == -1)) then
      call MOM_error(FATAL, "modify_vardesc called by "//trim(caller)//": "//trim(vd%name)//&
                   " has an unrecognized hor_grid argument "//trim(vd%hor_grid))
    elseif (vd%position == -1) then
      call MOM_error(FATAL, "modify_vardesc called with bad hor_grid argument "//trim(vd%hor_grid))
    endif
  endif
  if (present(z_grid))    call safe_string_copy(z_grid, vd%z_grid,     &
                               "vd%z_grid of "//trim(vd%name), cllr)
  if (present(t_grid))    call safe_string_copy(t_grid, vd%t_grid,     &
                               "vd%t_grid of "//trim(vd%name), cllr)

  if (present(cmor_field_name)) call safe_string_copy(cmor_field_name, vd%cmor_field_name, &
                                     "vd%cmor_field_name of "//trim(vd%name), cllr)
  if (present(cmor_units))      call safe_string_copy(cmor_units, vd%cmor_units, &
                                     "vd%cmor_units of "//trim(vd%name), cllr)
  if (present(cmor_longname))   call safe_string_copy(cmor_longname, vd%cmor_longname, &
                                     "vd%cmor_longname of "//trim(vd%name), cllr)

  if (present(conversion)) vd%conversion = conversion

  if (present(dim_names)) then
    do n=1,min(5,size(dim_names)) ; if (len_trim(dim_names(n)) > 0) then
      call safe_string_copy(dim_names(n), vd%dim_names(n), "vd%dim_names of "//trim(vd%name), cllr)
    endif ; enddo
  endif

  if (present(extra_axes)) then
    do n=1,size(extra_axes) ; if (len_trim(extra_axes(n)%name) > 0) then
      vd%extra_axes(n) = extra_axes(n)
    endif ; enddo
  endif

end procedure modify_vardesc
module procedure position_from_horgrid
  select case (trim(hor_grid))
    case ('h')   ; position_from_horgrid = CENTER
    case ('q')   ; position_from_horgrid = CORNER
    case ('u')   ; position_from_horgrid = EAST_FACE
    case ('v')   ; position_from_horgrid = NORTH_FACE
    case ('T')   ; position_from_horgrid = CENTER
    case ('Bu')  ; position_from_horgrid = CORNER
    case ('Cu')  ; position_from_horgrid = EAST_FACE
    case ('Cv')  ; position_from_horgrid = NORTH_FACE
    case ('1')   ; position_from_horgrid = 0
    case default ; position_from_horgrid = -1 ! This is a bad-value flag.
  end select
end procedure position_from_horgrid
module procedure set_axis_info
  call safe_string_copy(name, axis%name, "axis%name of "//trim(name), "set_axis_info")
  ! Set the default values.
  axis%longname = trim(axis%name) ; axis%units = "" ; axis%cartesian = "N" ; axis%sense = 0

  if (present(longname))  call safe_string_copy(longname, axis%longname, &
                                                "axis%longname of "//trim(name), "set_axis_info")
  if (present(units))     call safe_string_copy(units, axis%units, &
                                                "axis%units of "//trim(name), "set_axis_info")
  if (present(cartesian)) call safe_string_copy(cartesian, axis%cartesian, &
                                                "axis%cartesian of "//trim(name), "set_axis_info")
  if (present(sense)) axis%sense = sense

  if (.not.(present(ax_size) .or. present(ax_data)) ) then
    call MOM_error(FATAL, "set_axis_info called for "//trim(name)//&
      "without either an ax_size or an ax_data argument.")
  elseif (present(ax_size) .and. present(ax_data)) then
    if (size(ax_data) /= ax_size) call MOM_error(FATAL, "set_axis_info called for "//trim(name)//&
      "with an inconsistent value of ax_size and size of ax_data.")
  endif

  if (present(ax_size)) then
    axis%ax_size = ax_size
  else
    axis%ax_size = size(ax_data)
  endif
  if (present(ax_data)) then
    allocate(axis%ax_data(axis%ax_size)) ; axis%ax_data(:) = ax_data(:)
  endif

end procedure set_axis_info
module procedure delete_axis_info
  integer :: n
  do n=1,size(axes)
    axes(n)%name = "" ; axes(n)%longname = "" ; axes(n)%units = "" ; axes(n)%cartesian = "N"
    axes(n)%sense = 0 ; axes(n)%ax_size = 0
    if (allocated(axes(n)%ax_data)) deallocate(axes(n)%ax_data)
  enddo
end procedure delete_axis_info
module procedure get_axis_info
  if (present(ax_data)) then
    if (allocated(ax_data)) deallocate(ax_data)
    allocate(ax_data(axis%ax_size))
    ax_data(:) = axis%ax_data
  endif

  if (present(name)) name = axis%name
  if (present(longname)) longname = axis%longname
  if (present(units)) units = axis%units
  if (present(cartesian)) cartesian = axis%cartesian
  if (present(ax_size)) ax_size = axis%ax_size

end procedure get_axis_info
module procedure set_attribute_info
  attribute%name = trim(name)
  attribute%att_val = trim(str_value)
end procedure set_attribute_info
module procedure delete_attribute_info
  integer :: n
  do n=1,size(atts)
    if (allocated(atts(n)%name)) deallocate(atts(n)%name)
    if (allocated(atts(n)%att_val)) deallocate(atts(n)%att_val)
  enddo
end procedure delete_attribute_info
module procedure cmor_long_std
  integer :: k
  std_name = lowercase(longname)

  do k=1, len_trim(std_name)
    if (std_name(k:k) == ' ') std_name(k:k) = '_'
  enddo

end procedure cmor_long_std
module procedure query_vardesc
  integer :: n
  integer, parameter :: nmax_extraaxes = 5
  character(len=120) :: cllr, varname
  cllr = "mod_vardesc"
  if (present(caller)) cllr = trim(caller)

  if (present(name))      call safe_string_copy(vd%name, name,         &
                               "vd%name of "//trim(vd%name), cllr)
  if (present(longname))  call safe_string_copy(vd%longname, longname, &
                               "vd%longname of "//trim(vd%name), cllr)
  if (present(units))     call safe_string_copy(vd%units, units,       &
                               "vd%units of "//trim(vd%name), cllr)
  if (present(hor_grid))  call safe_string_copy(vd%hor_grid, hor_grid, &
                               "vd%hor_grid of "//trim(vd%name), cllr)
  if (present(z_grid))    call safe_string_copy(vd%z_grid, z_grid,     &
                               "vd%z_grid of "//trim(vd%name), cllr)
  if (present(t_grid))    call safe_string_copy(vd%t_grid, t_grid,     &
                               "vd%t_grid of "//trim(vd%name), cllr)

  if (present(cmor_field_name)) call safe_string_copy(vd%cmor_field_name, cmor_field_name, &
                                     "vd%cmor_field_name of "//trim(vd%name), cllr)
  if (present(cmor_units))      call safe_string_copy(vd%cmor_units, cmor_units,          &
                                     "vd%cmor_units of "//trim(vd%name), cllr)
  if (present(cmor_longname))   call safe_string_copy(vd%cmor_longname, cmor_longname, &
                                     "vd%cmor_longname of "//trim(vd%name), cllr)

  if (present(conversion)) conversion = vd%conversion

  if (present(position)) then
    position = vd%position
    if (position == -1) position = position_from_horgrid(vd%hor_grid)
  endif
  if (present(dim_names)) then
    do n=1,min(5,size(dim_names))
      call safe_string_copy(vd%dim_names(n), dim_names(n), "vd%dim_names of "//trim(vd%name), cllr)
    enddo
  endif

  if (present(extra_axes)) then
    ! save_restart expects 5 extra axes (can be empty)
    do n=1, nmax_extraaxes
      if (vd%extra_axes(n)%ax_size>=1) then
        extra_axes(n) = vd%extra_axes(n)
      else
        ! return an empty axis
        write(varname,"('dummy',i1.1)") n
        call set_axis_info(extra_axes(n), name=trim(varname), ax_size=1)
      endif
    enddo
  endif

end procedure query_vardesc
module procedure MOM_read_data_0d
  call read_field(filename, fieldname, data, &
                  timelevel=timelevel, scale=scale, MOM_Domain=MOM_Domain, &
                  global_file=global_file, file_may_be_4d=file_may_be_4d)
end procedure MOM_read_data_0d
module procedure MOM_read_data_0d_int
  call read_field(filename, fieldname, data, timelevel=timelevel)
end procedure MOM_read_data_0d_int
module procedure MOM_read_data_1d
  call read_field(filename, fieldname, data, &
                  timelevel=timelevel, scale=scale, MOM_Domain=MOM_Domain, &
                  global_file=global_file, file_may_be_4d=file_may_be_4d)

end procedure MOM_read_data_1d
module procedure MOM_read_data_1d_int
  call read_field(filename, fieldname, data, timelevel=timelevel)
end procedure MOM_read_data_1d_int
module procedure MOM_read_data_2d
  integer :: qturns   ! Number of quarter-turns from input to model grid
  real, allocatable :: data_in(:,:)  ! Field array on the input grid in arbitrary units [A ~> a]
  type(MOM_domain_type), pointer :: domain_ptr => NULL()  ! Pointer to the unrotated domain for reading
  qturns = MOM_domain%turns ; if (present(turns)) qturns = modulo(turns, 4)

  domain_ptr => MOM_Domain
  if (associated(MOM_Domain%domain_in) .and. (qturns /= 0)) domain_ptr => MOM_Domain%domain_in

  if (qturns == 0) then
    call read_field(filename, fieldname, data, MOM_Domain, &
                    timelevel=timelevel, position=position, scale=scale, &
                    global_file=global_file, file_may_be_4d=file_may_be_4d)
  else
    call allocate_rotated_array(data, [1,1], -qturns, data_in)
    call rotate_array(data, -qturns, data_in)
    call read_field(filename, fieldname, data_in, domain_ptr, &
                    timelevel=timelevel, position=position, scale=scale, &
                    global_file=global_file, file_may_be_4d=file_may_be_4d)
    call rotate_array(data_in, qturns, data)
    deallocate(data_in)
  endif

end procedure MOM_read_data_2d
module procedure read_netCDF_data_2d
  integer :: qturns
  real, allocatable :: values_in(:,:)
  type(MOM_netcdf_file) :: handle
  if (present(position)) &
    call MOM_error(FATAL, 'read_netCDF_data: position is not yet supported.')

  ! Timelevels are not yet supported
  if (present(timelevel)) &
    call MOM_error(FATAL, 'read_netCDF_data: timelevel is not yet supported.')

  call handle%open(filename, action=READONLY_FILE, MOM_domain=MOM_domain)
  call handle%update()

  qturns = MOM_domain%turns ; if (present(turns)) qturns = modulo(turns, 4)

  if (qturns == 0) then
    call handle%read(fieldname, values, rescale=rescale)
  else
    call allocate_rotated_array(values, [1,1], -qturns, values_in)
    call rotate_array(values, -qturns, values_in)
    call handle%read(fieldname, values_in, rescale=rescale)
    call rotate_array(values_in, qturns, values)
    deallocate(values_in)
  endif

  call handle%close()
end procedure read_netCDF_data_2d
module procedure MOM_read_data_2d_region
  integer :: qturns                   ! Number of quarter turns
  real, allocatable :: data_in(:,:)   ! Field array on the input grid in arbitrary units [A ~> a]
  qturns = 0
  if (present(turns)) qturns = modulo(turns, 4)

  if (qturns == 0) then
    call read_field(filename, fieldname, data, start, nread, &
                    MOM_Domain=MOM_Domain, no_domain=no_domain, scale=scale)
  else
    call allocate_rotated_array(data, [1,1], -qturns, data_in)
    call rotate_array(data, -qturns, data_in)
    if (associated(MOM_Domain%domain_in)) then
      call read_field(filename, fieldname, data_in, start, nread, &
                      MOM_Domain=MOM_Domain%domain_in, no_domain=no_domain, scale=scale)
    else
      call read_field(filename, fieldname, data_in, start, nread, &
                      MOM_Domain=MOM_Domain, no_domain=no_domain, scale=scale)
    endif
    call rotate_array(data_in, qturns, data)
    deallocate(data_in)
  endif
end procedure MOM_read_data_2d_region
module procedure MOM_read_data_3d
  integer :: qturns   ! Number of quarter-turns from input to model grid
  real, allocatable :: data_in(:,:,:)  ! Field array on the input grid in arbitrary units [A ~> a]
  type(MOM_domain_type), pointer :: domain_ptr => NULL()  ! Pointer to the unrotated domain for reading
  domain_ptr => MOM_Domain
  if (associated(MOM_Domain%domain_in) .and. (qturns /= 0)) domain_ptr => MOM_Domain%domain_in

  qturns = MOM_domain%turns ; if (present(turns)) qturns = modulo(turns, 4)
  if (qturns == 0) then
    call read_field(filename, fieldname, data, MOM_Domain, &
                    timelevel=timelevel, position=position, scale=scale, &
                    global_file=global_file, file_may_be_4d=file_may_be_4d)
  else
    call allocate_rotated_array(data, [1,1,1], -qturns, data_in)
    call rotate_array(data, -qturns, data_in)
    call read_field(filename, fieldname, data_in, domain_ptr, &
                    timelevel=timelevel, position=position, scale=scale, &
                    global_file=global_file, file_may_be_4d=file_may_be_4d)
    call rotate_array(data_in, qturns, data)
    deallocate(data_in)
  endif

end procedure MOM_read_data_3d
module procedure MOM_read_data_3d_region
  integer :: qturns                   ! Number of quarter turns
  real, allocatable :: data_in(:,:,:)   ! Field array on the input grid in arbitrary units [A ~> a]
  qturns = 0
  if (present(turns)) qturns = modulo(turns, 4)

  if (qturns == 0) then
    call read_field(filename, fieldname, data, start, nread, &
                    MOM_Domain=MOM_Domain, no_domain=no_domain, scale=scale)
  else
    call allocate_rotated_array(data, [1,1,1], -qturns, data_in)
    call rotate_array(data, -qturns, data_in)
    if (associated(MOM_Domain%domain_in)) then
      call read_field(filename, fieldname, data_in, start, nread, &
                      MOM_Domain=MOM_Domain%domain_in, no_domain=no_domain, scale=scale)
    else
      call read_field(filename, fieldname, data_in, start, nread, &
                      MOM_Domain=MOM_Domain, no_domain=no_domain, scale=scale)
    endif
    call rotate_array(data_in, qturns, data)
    deallocate(data_in)
  endif
end procedure MOM_read_data_3d_region
module procedure MOM_read_data_4d
  integer :: qturns   ! Number of quarter-turns from input to model grid
  real, allocatable :: data_in(:,:,:,:)  ! Field array on the input grid in arbitrary units [A ~> a]
  type(MOM_domain_type), pointer :: domain_ptr => NULL()  ! Pointer to the unrotated domain for reading
  qturns = MOM_domain%turns ; if (present(turns)) qturns = modulo(turns, 4)

  domain_ptr => MOM_Domain
  if (associated(MOM_Domain%domain_in) .and. (qturns /= 0)) domain_ptr => MOM_Domain%domain_in

  if (qturns == 0) then
    call read_field(filename, fieldname, data, MOM_Domain, &
        timelevel=timelevel, position=position, scale=scale, &
        global_file=global_file)
  else
    ! Read field along the input grid and rotate to the model grid
    call allocate_rotated_array(data, [1,1,1,1], -qturns, data_in)
    call rotate_array(data, -qturns, data_in)
    call read_field(filename, fieldname, data_in, domain_ptr, timelevel=timelevel, &
                    position=position, scale=scale, global_file=global_file)
    call rotate_array(data_in, qturns, data)
    deallocate(data_in)
  endif

end procedure MOM_read_data_4d
module procedure MOM_read_vector_2d
  integer :: qturns ! Number of quarter-turns from input to model grid
  real, allocatable :: u_data_in(:,:), v_data_in(:,:)   ! [uv] on the input grid in arbitrary units [A ~> a]
  type(MOM_domain_type), pointer :: domain_ptr => NULL()  ! Pointer to the unrotated domain for reading
  qturns = MOM_domain%turns ; if (present(turns)) qturns = modulo(turns, 4)

  domain_ptr => MOM_Domain
  if (associated(MOM_Domain%domain_in) .and. (qturns /= 0)) domain_ptr => MOM_Domain%domain_in

  if (qturns == 0) then
    call read_vector(filename, u_fieldname, v_fieldname, &
                     u_data, v_data, MOM_domain, timelevel=timelevel, stagger=stagger, &
                     scalar_pair=scalar_pair, scale=scale)
  else
    call allocate_rotated_array(u_data, [1,1], -qturns, u_data_in)
    call allocate_rotated_array(v_data, [1,1], -qturns, v_data_in)
    if (scalar_pair) then
      call rotate_array_pair(u_data, v_data, -qturns, u_data_in, v_data_in)
    else
      call rotate_vector(u_data, v_data, -qturns, u_data_in, v_data_in)
    endif
    call read_vector(filename, u_fieldname, v_fieldname, u_data_in, v_data_in, &
                     domain_ptr, timelevel=timelevel, &
                     stagger=stagger, scalar_pair=scalar_pair, scale=scale)
    if (scalar_pair) then
      call rotate_array_pair(u_data_in, v_data_in, qturns, u_data, v_data)
    else
      call rotate_vector(u_data_in, v_data_in, qturns, u_data, v_data)
    endif
    deallocate(v_data_in)
    deallocate(u_data_in)
  endif

end procedure MOM_read_vector_2d
module procedure MOM_read_vector_3d
  integer :: qturns ! Number of quarter-turns from input to model grid
  real, allocatable :: u_data_in(:,:,:), v_data_in(:,:,:) ! [uv] on the input grid in arbitrary units [A ~> a]
  type(MOM_domain_type), pointer :: domain_ptr => NULL()  ! Pointer to the unrotated domain for reading
  qturns = MOM_domain%turns ; if (present(turns)) qturns = modulo(turns, 4)

  domain_ptr => MOM_Domain
  if (associated(MOM_Domain%domain_in) .and. (qturns /= 0)) domain_ptr => MOM_Domain%domain_in

  if (qturns == 0) then
    call read_vector(filename, u_fieldname, v_fieldname, &
                     u_data, v_data, MOM_domain, timelevel=timelevel, stagger=stagger, &
                     scalar_pair=scalar_pair, scale=scale)
  else
    call allocate_rotated_array(u_data, [1,1,1], -qturns, u_data_in)
    call allocate_rotated_array(v_data, [1,1,1], -qturns, v_data_in)
    if (scalar_pair) then
      call rotate_array_pair(u_data, v_data, -qturns, u_data_in, v_data_in)
    else
      call rotate_vector(u_data, v_data, -qturns, u_data_in, v_data_in)
    endif
    call read_vector(filename, u_fieldname, v_fieldname, u_data_in, v_data_in, &
                     domain_ptr, timelevel=timelevel, &
                     stagger=stagger, scalar_pair=scalar_pair, scale=scale)
    if (scalar_pair) then
      call rotate_array_pair(u_data_in, v_data_in, qturns, u_data, v_data)
    else
      call rotate_vector(u_data_in, v_data_in, qturns, u_data, v_data)
    endif
    deallocate(v_data_in)
    deallocate(u_data_in)
  endif

end procedure MOM_read_vector_3d
module procedure MOM_write_field_legacy_4d
  real, allocatable :: field_rot(:,:,:,:)  ! A rotated version of field, with the same units [a] or
  real :: scale_fac ! A scaling factor to use before writing the array [a A-1 ~> 1]
  integer :: qturns ! The number of quarter turns through which to rotate field
  qturns = 0 ; if (present(turns)) qturns = modulo(turns, 4)
  scale_fac = 1.0 ; if (present(scale)) scale_fac = scale
  if (present(unscale)) scale_fac = unscale

  if ((qturns == 0) .and. (scale_fac == 1.0) .and. .not.present(zero_zeros)) then
    call write_field(IO_handle, field_md, MOM_domain, field, tstamp=tstamp, &
                         tile_count=tile_count, fill_value=fill_value)
  else
    call allocate_rotated_array(field, [1,1,1,1], qturns, field_rot)
    call rotate_array(field, qturns, field_rot)
    call rescale_comp_data(MOM_Domain, field_rot, scale_fac, zero_zeros)
    call write_field(IO_handle, field_md, MOM_domain, field_rot, tstamp=tstamp, &
                         tile_count=tile_count, fill_value=fill_value)
    deallocate(field_rot)
  endif
end procedure MOM_write_field_legacy_4d
module procedure MOM_write_field_legacy_3d
  real, allocatable :: field_rot(:,:,:)  ! A rotated version of field, with the same units [a] or
  real :: scale_fac ! A scaling factor to use before writing the array [a A-1 ~> 1]
  integer :: qturns ! The number of quarter turns through which to rotate field
  qturns = 0 ; if (present(turns)) qturns = modulo(turns, 4)
  scale_fac = 1.0 ; if (present(scale)) scale_fac = scale
  if (present(unscale)) scale_fac = unscale

  if ((qturns == 0) .and. (scale_fac == 1.0) .and. .not.present(zero_zeros)) then
    call write_field(IO_handle, field_md, MOM_domain, field, tstamp=tstamp, &
                         tile_count=tile_count, fill_value=fill_value)
  else
    call allocate_rotated_array(field, [1,1,1], qturns, field_rot)
    call rotate_array(field, qturns, field_rot)
    call rescale_comp_data(MOM_Domain, field_rot, scale_fac, zero_zeros)
    call write_field(IO_handle, field_md, MOM_domain, field_rot, tstamp=tstamp, &
                         tile_count=tile_count, fill_value=fill_value)
    deallocate(field_rot)
  endif
end procedure MOM_write_field_legacy_3d
module procedure MOM_write_field_legacy_2d
  real, allocatable :: field_rot(:,:)  ! A rotated version of field, with the same units [a] or
  real :: scale_fac ! A scaling factor to use before writing the array [a A-1 ~> 1]
  integer :: qturns ! The number of quarter turns through which to rotate field
  qturns = 0 ; if (present(turns)) qturns = modulo(turns, 4)
  scale_fac = 1.0 ; if (present(scale)) scale_fac = scale
  if (present(unscale)) scale_fac = unscale

  if ((qturns == 0) .and. (scale_fac == 1.0) .and. .not.present(zero_zeros)) then
    call write_field(IO_handle, field_md, MOM_domain, field, tstamp=tstamp, &
                         tile_count=tile_count, fill_value=fill_value)
  else
    call allocate_rotated_array(field, [1,1], qturns, field_rot)
    call rotate_array(field, qturns, field_rot)
    call rescale_comp_data(MOM_Domain, field_rot, scale_fac, zero_zeros)
    call write_field(IO_handle, field_md, MOM_domain, field_rot, tstamp=tstamp, &
                         tile_count=tile_count, fill_value=fill_value)
    deallocate(field_rot)
  endif
end procedure MOM_write_field_legacy_2d
module procedure MOM_write_field_legacy_1d
  real, dimension(:), allocatable :: array ! A rescaled copy of field [a]
  real :: scale_fac ! A scaling factor to use before writing the array [a A-1 ~> 1]
  logical :: design_zeros ! If true, convert negative zeros into ordinary signless zeros.
  integer :: i
  scale_fac = 1.0 ; if (present(scale)) scale_fac = scale
  if (present(unscale)) scale_fac = unscale

  design_zeros = .false. ; if (present(zero_zeros)) design_zeros = zero_zeros

  if ((scale_fac == 1.0) .and. (.not.design_zeros)) then
    call write_field(IO_handle, field_md, field, tstamp=tstamp)
  else
    allocate(array(size(field)))
    array(:) = scale_fac * field(:)
    if (present(fill_value)) then
      do i=1,size(field) ; if (field(i) == fill_value) array(i) = fill_value ; enddo
    endif
    if (design_zeros) then ! Convert negative zeros into zeros
      do i=1,size(field) ; if (array(i) == 0.0) array(i) = 0.0 ; enddo
    endif
    call write_field(IO_handle, field_md, array, tstamp=tstamp)
    deallocate(array)
  endif
end procedure MOM_write_field_legacy_1d
module procedure MOM_write_field_legacy_0d
  real :: scale_fac  ! A scaling factor to use before writing the field [a A-1 ~> 1]
  real :: scaled_val ! A rescaled copy of field [a]
  scale_fac = 1.0 ; if (present(scale)) scale_fac = scale
  if (present(unscale)) scale_fac = unscale

  scaled_val = field * scale_fac

  if (present(fill_value)) then ; if (field == fill_value) scaled_val = fill_value ; endif
  if (present(zero_zeros)) then ; if (zero_zeros .and. (scaled_val == 0.0)) scaled_val = 0.0 ; endif

  call write_field(IO_handle, field_md, scaled_val, tstamp=tstamp)
end procedure MOM_write_field_legacy_0d
module procedure MOM_write_field_4d
  real, allocatable :: field_rot(:,:,:,:)  ! A rotated version of field, with the same units or rescaled [a]
  real :: scale_fac ! A scaling factor to use before writing the array [a A-1 ~> 1]
  integer :: qturns ! The number of quarter turns through which to rotate field
  qturns = 0 ; if (present(turns)) qturns = modulo(turns, 4)
  scale_fac = 1.0 ; if (present(scale)) scale_fac = scale
  if (present(unscale)) scale_fac = unscale

  if ((qturns == 0) .and. (scale_fac == 1.0) .and. .not.present(zero_zeros)) then
    call IO_handle%write_field(field_md, MOM_domain, field, tstamp=tstamp, &
                         tile_count=tile_count, fill_value=fill_value)
  else
    call allocate_rotated_array(field, [1,1,1,1], qturns, field_rot)
    call rotate_array(field, qturns, field_rot)
    call rescale_comp_data(MOM_Domain, field_rot, scale_fac, zero_zeros)
    call IO_handle%write_field(field_md, MOM_domain, field_rot, tstamp=tstamp, &
                         tile_count=tile_count, fill_value=fill_value)
    deallocate(field_rot)
  endif
end procedure MOM_write_field_4d
module procedure MOM_write_field_3d
  real, allocatable :: field_rot(:,:,:)  ! A rotated version of field, with the same units or rescaled [a]
  real :: scale_fac ! A scaling factor to use before writing the array [a A-1 ~> 1]
  integer :: qturns ! The number of quarter turns through which to rotate field
  qturns = 0 ; if (present(turns)) qturns = modulo(turns, 4)
  scale_fac = 1.0 ; if (present(scale)) scale_fac = scale
  if (present(unscale)) scale_fac = unscale

  if ((qturns == 0) .and. (scale_fac == 1.0) .and. .not.present(zero_zeros)) then
    call IO_handle%write_field(field_md, MOM_domain, field, tstamp=tstamp, &
                         tile_count=tile_count, fill_value=fill_value)
  else
    call allocate_rotated_array(field, [1,1,1], qturns, field_rot)
    call rotate_array(field, qturns, field_rot)
    call rescale_comp_data(MOM_Domain, field_rot, scale_fac, zero_zeros)
    call IO_handle%write_field(field_md, MOM_domain, field_rot, tstamp=tstamp, &
                         tile_count=tile_count, fill_value=fill_value)
    deallocate(field_rot)
  endif
end procedure MOM_write_field_3d
module procedure MOM_write_field_2d
  real, allocatable :: field_rot(:,:)  ! A rotated version of field, with the same units or rescaled [a]
  real :: scale_fac ! A scaling factor to use before writing the array [a A-1 ~> 1]
  integer :: qturns ! The number of quarter turns through which to rotate field
  qturns = 0 ; if (present(turns)) qturns = modulo(turns, 4)
  scale_fac = 1.0 ; if (present(scale)) scale_fac = scale
  if (present(unscale)) scale_fac = unscale

  if ((qturns == 0) .and. (scale_fac == 1.0) .and. .not.present(zero_zeros)) then
    call IO_handle%write_field(field_md, MOM_domain, field, tstamp=tstamp, &
                         tile_count=tile_count, fill_value=fill_value)
  else
    call allocate_rotated_array(field, [1,1], qturns, field_rot)
    call rotate_array(field, qturns, field_rot)
    call rescale_comp_data(MOM_Domain, field_rot, scale_fac, zero_zeros)
    call IO_handle%write_field(field_md, MOM_domain, field_rot, tstamp=tstamp, &
                         tile_count=tile_count, fill_value=fill_value)
    deallocate(field_rot)
  endif
end procedure MOM_write_field_2d
module procedure MOM_write_field_1d
  real, dimension(:), allocatable :: array ! A rescaled copy of field in arbtrary unscaled units [a]
  real :: scale_fac ! A scaling factor to use before writing the array [a A-1 ~> 1]
  logical :: design_zeros ! If true, convert negative zeros into ordinary signless zeros.
  integer :: i
  scale_fac = 1.0 ; if (present(scale)) scale_fac = scale
  if (present(unscale)) scale_fac = unscale

  design_zeros = .false. ; if (present(zero_zeros)) design_zeros = zero_zeros

  if ((scale_fac == 1.0) .and. (.not.design_zeros)) then
    call IO_handle%write_field(field_md, field, tstamp=tstamp)
  else
    allocate(array(size(field)))
    array(:) = scale_fac * field(:)
    if (present(fill_value)) then
      do i=1,size(field) ; if (field(i) == fill_value) array(i) = fill_value ; enddo
    endif
    if (design_zeros) then ! Convert negative zeros into zeros
      do i=1,size(field) ; if (array(i) == 0.0) array(i) = 0.0 ; enddo
    endif
    call IO_handle%write_field(field_md, array, tstamp=tstamp)
    deallocate(array)
  endif
end procedure MOM_write_field_1d
module procedure MOM_write_field_0d
  real :: scale_fac  ! A scaling factor to use before writing the field [a A-1 ~> 1]
  real :: scaled_val ! A rescaled copy of field in arbtrary unscaled units [a]
  scale_fac = 1.0 ; if (present(scale)) scale_fac = scale
  if (present(unscale)) scale_fac = unscale

  scaled_val = field * scale_fac

  if (present(fill_value)) then ; if (field == fill_value) scaled_val = fill_value ; endif
  if (present(zero_zeros)) then ; if (zero_zeros .and. (scaled_val == 0.0)) scaled_val = 0.0 ; endif

  call IO_handle%write_field(field_md, scaled_val, tstamp=tstamp)
end procedure MOM_write_field_0d
module procedure field_size
  if (present(ndims)) then
    if (present(no_domain)) then ; if (.not.no_domain) call MOM_error(FATAL, &
          "field_size does not support the ndims argument when no_domain is present and false.")
    endif
    call get_var_sizes(filename, fieldname, ndims, sizes, match_case=.false., ncid_in=ncid_in)
    if (present(field_found)) field_found = (ndims >= 0)
    if ((ndims < 0) .and. .not.present(field_found)) then
      call MOM_error(FATAL, "Variable "//trim(fieldname)//" not found in "//trim(filename) )
    endif
  else
    call get_field_size(filename, fieldname, sizes, field_found=field_found, no_domain=no_domain)
  endif

end procedure field_size
module procedure safe_string_copy
  if (len(trim(str1)) > len(str2)) then
    if (present(fieldnm) .and. present(caller)) then
      call MOM_error(FATAL, trim(caller)//" attempted to copy the overly long string "//&
                     trim(str1)//" into "//trim(fieldnm))
    else
      call MOM_error(FATAL, "safe_string_copy: The string "//trim(str1)//&
                     " is longer than its intended target.")
    endif
  endif
  str2 = trim(str1)
end procedure safe_string_copy
module procedure ensembler
  character(len=len(name)) :: tmp
  character(10) :: ens_num_char
  character(3)  :: code_str
  integer :: ens_no
  integer :: n, is
  en_nm = trim(name)
  if (index(name,"%") == 0) return

  if (present(ens_no_in)) then
    ens_no = ens_no_in
  else
    ens_no = get_ensemble_id()
  endif

  write(ens_num_char, '(I0)') ens_no
  do
    is = index(en_nm,"%E")
    if (is == 0) exit
    if (len(en_nm) < len(trim(en_nm)) + len(trim(ens_num_char)) - 2) &
      call MOM_error(FATAL, "MOM_io ensembler: name "//trim(name)// &
      " is not long enough for %E expansion for ens_no "//trim(ens_num_char))
    tmp = en_nm(1:is-1)//trim(ens_num_char)//trim(en_nm(is+2:))
    en_nm = tmp
  enddo

  if (index(name,"%") == 0) return

  write(ens_num_char, '(I10.10)') ens_no
  do n=1,9 ; do
    write(code_str, '("%",I1,"E")') n

    is = index(en_nm,code_str)
    if (is == 0) exit
    if (ens_no < 10**n) then
      if (len(en_nm) < len(trim(en_nm)) + n-3) call MOM_error(FATAL, &
        "MOM_io ensembler: name "//trim(name)//" is not long enough for %E expansion.")
      tmp = en_nm(1:is-1)//trim(ens_num_char(11-n:10))//trim(en_nm(is+3:))
    else
      call MOM_error(FATAL, "MOM_io ensembler: Ensemble number is too large "//&
          "to be encoded with "//code_str//" in "//trim(name))
    endif
    en_nm = tmp
  enddo ; enddo

end procedure ensembler
module procedure get_filename_appendix
  call get_filename_suffix(suffix)
end procedure get_filename_appendix
module procedure write_version_number
  call write_version(version, tag, unit)
end procedure write_version_number
module procedure open_namelist_file
  unit = MOM_namelist_file(file)
end procedure open_namelist_file
module procedure check_nml_error
  call check_namelist_error(IOstat, nml_name)
  ierr = IOstat
end procedure check_nml_error
module procedure MOM_io_init
# include "version_variable.h"
  character(len=40)  :: mdl = "MOM_io" ! This module's name.
  call log_version(param_file, mdl, version)

end procedure MOM_io_init
module procedure get_var_axes_info
  integer ::  rcode
  logical :: success
  integer ::  ncid, varid, ndims
  integer :: id, jd, kd
  integer, dimension(4) :: dims, dim_id
  character(len=128)  :: dim_name(4)
  integer, dimension(1) :: start, count
  real, allocatable, dimension(:) :: x ! x-axis labels, often [degrees_E] or [km] or [m]
  real, allocatable, dimension(:) :: y ! y-axis labels, often [degrees_N] or [km] or [m]
  real, allocatable, dimension(:) :: z ! vertical axis labels [various], often [m] or [kg m-3]
  call open_file_to_read(filename, ncid, success=success)

  rcode = NF90_INQ_VARID(ncid, trim(fieldname), varid)
  if (rcode /= 0) call MOM_error(FATAL,"error finding variable "//trim(fieldname)//&
                                 " in file "//trim(filename)//" in hinterp_extrap")

  rcode = NF90_INQUIRE_VARIABLE(ncid, varid, ndims=ndims, dimids=dims)
  if (rcode /= 0) call MOM_error(FATAL, "Error inquiring about the dimensions of "//trim(fieldname)//&
                                 " in file "//trim(filename)//" in hinterp_extrap")
  if (ndims < 3) call MOM_error(FATAL,"Variable "//trim(fieldname)//" in file "//trim(filename)// &
                                " has too few dimensions to be read as a 3-d array.")
  rcode = NF90_INQUIRE_DIMENSION(ncid, dims(1), dim_name(1), len=id)
  if (rcode /= 0) call MOM_error(FATAL,"error reading dimension 1 data for "// &
                trim(fieldname)//" in file "// trim(filename)//" in hinterp_extrap")
  rcode = NF90_INQ_VARID(ncid, dim_name(1), dim_id(1))
  if (rcode /= 0) call MOM_error(FATAL,"error finding variable "//trim(dim_name(1))//&
                                 " in file "//trim(filename)//" in hinterp_extrap")
  rcode = NF90_INQUIRE_DIMENSION(ncid, dims(2), dim_name(2), len=jd)
  if (rcode /= 0) call MOM_error(FATAL,"error reading dimension 2 data for "// &
                trim(fieldname)//" in file "// trim(filename)//" in hinterp_extrap")
  rcode = NF90_INQ_VARID(ncid, dim_name(2), dim_id(2))
  if (rcode /= 0) call MOM_error(FATAL,"error finding variable "//trim(dim_name(2))//&
                                 " in file "//trim(filename)//" in hinterp_extrap")
  rcode = NF90_INQUIRE_DIMENSION(ncid, dims(3), dim_name(3), len=kd)
  if (rcode /= 0) call MOM_error(FATAL,"error reading dimension 3 data for "// &
                trim(fieldname)//" in file "// trim(filename)//" in hinterp_extrap")
  rcode = NF90_INQ_VARID(ncid, dim_name(3), dim_id(3))
  if (rcode /= 0) call MOM_error(FATAL,"error finding variable "//trim(dim_name(3))//&
                                 " in file "//trim(filename)//" in hinterp_extrap")
  allocate(x(id), y(jd), z(kd))

  start = 1 ; count = 1 ; count(1) = id
  rcode = NF90_GET_VAR(ncid, dim_id(1), x, start, count)
  if (rcode /= 0) call MOM_error(FATAL,"error reading dimension 1 values for var_name "// &
                trim(fieldname)//",dim_name "//trim(dim_name(1))//" in file "// trim(filename)//" in hinterp_extrap")
  start = 1 ; count = 1 ; count(1) = jd
  rcode = NF90_GET_VAR(ncid, dim_id(2), y, start, count)
  if (rcode /= 0) call MOM_error(FATAL,"error reading dimension 2 values for var_name "// &
                trim(fieldname)//",dim_name "//trim(dim_name(2))//" in file "// trim(filename)//" in  hinterp_extrap")
  start = 1 ; count = 1 ; count(1) = kd
  rcode = NF90_GET_VAR(ncid, dim_id(3), z, start, count)
  if (rcode /= 0) call MOM_error(FATAL,"error reading dimension 3 values for var_name "// &
                trim(fieldname//",dim_name "//trim(dim_name(3)))//" in file "// trim(filename)//" in  hinterp_extrap")

  call set_axis_info(axes_info(1), name=trim(dim_name(1)), ax_size=id, ax_data=x,cartesian='X')
  call set_axis_info(axes_info(2), name=trim(dim_name(2)), ax_size=jd, ax_data=y,cartesian='Y')
  call set_axis_info(axes_info(3), name=trim(dim_name(3)), ax_size=kd, ax_data=z,cartesian='Z')

  call close_file_to_read(ncid, filename)

  deallocate(x,y,z)

end procedure get_var_axes_info
end submodule MOM_io_s
