submodule (MOM_io_file) MOM_io_file_s
  implicit none
contains
module procedure initialize_axis_list_infra
  allocate(list%head)
  list%tail => list%head
end procedure initialize_axis_list_infra
module procedure append_axis_list_infra
  type(axis_node_infra), pointer :: empty_node
  list%tail%label = label
  list%tail%axis = axis

  ! Extend list to next empty node
  allocate(empty_node)
  list%tail%next => empty_node
  list%tail => empty_node
end procedure append_axis_list_infra
module procedure get_axis_list_infra
  type(axis_node_infra), pointer :: node
  node => list%head
  do while(associated(node%next))
    if (node%label == label) exit
    node => node%next
  enddo
  if (.not. associated(node)) &
    call MOM_error(FATAL, "axis associated with " // label // " not found.")

  axis = node%axis
end procedure get_axis_list_infra
module procedure finalize_axis_list_infra
  type(axis_node_infra), pointer :: node, next_node
  node => list%head
  do while(associated(node))
    next_node => node
    node => node%next
    deallocate(next_node)
  enddo
end procedure finalize_axis_list_infra
module procedure initialize_axis_list_nc
  allocate(list%head)
  list%tail => list%head
end procedure initialize_axis_list_nc
module procedure append_axis_list_nc
  type(axis_node_nc), pointer :: empty_node
  list%tail%label = label
  list%tail%axis = axis

  ! Extend list to next empty node
  allocate(empty_node)
  list%tail%next => empty_node
  list%tail => empty_node
end procedure append_axis_list_nc
module procedure get_axis_list_nc
  type(axis_node_nc), pointer :: node
  node => list%head
  do while(associated(node%next))
    if (node%label == label) exit
    node => node%next
  enddo
  if (.not. associated(node)) &
    call MOM_error(FATAL, "axis associated with " // label // " not found.")

  axis = node%axis
end procedure get_axis_list_nc
module procedure finalize_axis_list_nc
  type(axis_node_nc), pointer :: node, next_node
  node => list%head
  do while(associated(node))
    next_node => node
    node => node%next
    deallocate(next_node)
  enddo
end procedure finalize_axis_list_nc
module procedure initialize_field_list_infra
  allocate(list%head)
  list%tail => list%head
end procedure initialize_field_list_infra
module procedure append_field_list_infra
  type(field_node_infra), pointer :: empty_node
  list%tail%label = label
  list%tail%field = field

  ! Extend list to next empty node
  allocate(empty_node)
  list%tail%next => empty_node
  list%tail => empty_node
end procedure append_field_list_infra
module procedure get_field_list_infra
  type(field_node_infra), pointer :: node
  node => list%head
  do while(associated(node%next))
    if (node%label == label) exit
    node => node%next
  enddo
  if (.not. associated(node)) &
    call MOM_error(FATAL, "field associated with " // label // " not found.")

  field = node%field
end procedure get_field_list_infra
module procedure finalize_field_list_infra
  type(field_node_infra), pointer :: node, next_node
  node => list%head
  do while(associated(node))
    next_node => node
    node => node%next
    deallocate(next_node)
  enddo
end procedure finalize_field_list_infra
module procedure initialize_field_list_nc
  allocate(list%head)
  list%tail => list%head
end procedure initialize_field_list_nc
module procedure append_field_list_nc
  type(field_node_nc), pointer :: empty_node
  list%tail%label = label
  list%tail%field = field

  ! Extend list to next empty node
  allocate(empty_node)
  list%tail%next => empty_node
  list%tail => empty_node
end procedure append_field_list_nc
module procedure get_field_list_nc
  type(field_node_nc), pointer :: node
  node => list%head
  do while(associated(node%next))
    if (node%label == label) exit
    node => node%next
  enddo
  if (.not. associated(node)) &
    call MOM_error(FATAL, "field associated with " // label // " not found.")

  field = node%field
end procedure get_field_list_nc
module procedure finalize_field_list_nc
  type(field_node_nc), pointer :: node, next_node
  node => list%head
  do while(associated(node))
    next_node => node
    node => node%next
    deallocate(next_node)
  enddo
end procedure finalize_field_list_nc
module procedure open_file_infra
  logical :: use_single_file_domain
  use_single_file_domain = .false.
  if (present(MOM_domain) .and. present(fileset)) then
    if (fileset == SINGLE_FILE) &
      use_single_file_domain = .true.
  endif

  if (use_single_file_domain) then
    call clone_MOM_domain(MOM_domain, handle%domain, io_layout=[1,1])
    call open_file(handle%handle_infra, filename, action=action, &
        MOM_domain=handle%domain, threading=threading, fileset=fileset)
  else
    call open_file(handle%handle_infra, filename, action=action, &
        MOM_domain=MOM_domain, threading=threading, fileset=fileset)
  endif

  call handle%axes%init()
  call handle%fields%init()
end procedure open_file_infra
module procedure close_file_infra
  if (associated(handle%domain)) &
    call deallocate_MOM_domain(handle%domain)

  call close_file(handle%handle_infra)
  call handle%axes%finalize()
  call handle%fields%finalize()
end procedure close_file_infra
module procedure flush_file_infra
  call flush_file(handle%handle_infra)
end procedure flush_file_infra
module procedure register_axis_infra
  type(axistype) :: ax_infra
  call write_metadata(handle%handle_infra, ax_infra, label, units, longname, &
      cartesian=cartesian, sense=sense, domain=domain, data=data, &
      edge_axis=edge_axis, calendar=calendar)

  call handle%axes%append(ax_infra, label)
  axis%label = label
end procedure register_axis_infra
module procedure register_field_infra
  type(fieldtype) :: field_infra
  type(axistype), allocatable :: field_axes(:)
  integer :: i
  allocate(field_axes(size(axes)))
  do i = 1, size(axes)
    field_axes(i) = handle%axes%get(axes(i)%label)
  enddo

  call write_metadata(handle%handle_infra, field_infra, field_axes, label, &
      units, longname, pack=pack, standard_name=standard_name, checksum=checksum)

  call handle%fields%append(field_infra, label)
  field%label = label
  field%conversion = 1.0 ; if (present(conversion)) field%conversion = conversion
end procedure register_field_infra
module procedure write_field_4d_infra
  type(fieldtype) :: field_infra
  real, allocatable :: unscaled_field(:,:,:,:) ! An unscaled version of field for output [a]
  field_infra = handle%fields%get(field_md%label)
  if (field_md%conversion == 1.0) then
    call write_field(handle%handle_infra, field_infra, MOM_domain, field, &
        tstamp=tstamp, tile_count=tile_count, fill_value=fill_value)
  else
    allocate(unscaled_field, source=field)
    unscaled_field(:,:,:,:) = field_md%conversion * field(:,:,:,:)
    call write_field(handle%handle_infra, field_infra, MOM_domain, unscaled_field, &
        tstamp=tstamp, tile_count=tile_count, fill_value=fill_value)
    deallocate(unscaled_field)
  endif
end procedure write_field_4d_infra
module procedure write_field_3d_infra
  type(fieldtype) :: field_infra
  real, allocatable :: unscaled_field(:,:,:) ! An unscaled version of field for output [a]
  field_infra = handle%fields%get(field_md%label)
  if (field_md%conversion == 1.0) then
    call write_field(handle%handle_infra, field_infra, MOM_domain, field, &
        tstamp=tstamp, tile_count=tile_count, fill_value=fill_value)
  else
    allocate(unscaled_field, source=field)
    unscaled_field(:,:,:) = field_md%conversion * field(:,:,:)
    call write_field(handle%handle_infra, field_infra, MOM_domain, unscaled_field, &
        tstamp=tstamp, tile_count=tile_count, fill_value=fill_value)
    deallocate(unscaled_field)
  endif

end procedure write_field_3d_infra
module procedure write_field_2d_infra
  type(fieldtype) :: field_infra
  real, allocatable :: unscaled_field(:,:) ! An unscaled version of field for output [a]
  field_infra = handle%fields%get(field_md%label)
  if (field_md%conversion == 1.0) then
    call write_field(handle%handle_infra, field_infra, MOM_domain, field, &
        tstamp=tstamp, tile_count=tile_count, fill_value=fill_value)
  else
    allocate(unscaled_field, source=field)
    unscaled_field(:,:) = field_md%conversion * field(:,:)
    call write_field(handle%handle_infra, field_infra, MOM_domain, unscaled_field, &
        tstamp=tstamp, tile_count=tile_count, fill_value=fill_value)
    deallocate(unscaled_field)
  endif
end procedure write_field_2d_infra
module procedure write_field_1d_infra
  type(fieldtype) :: field_infra
  real, allocatable :: unscaled_field(:) ! An unscaled version of field for output [a]
  field_infra = handle%fields%get(field_md%label)
  if (field_md%conversion == 1.0) then
    call write_field(handle%handle_infra, field_infra, field, tstamp=tstamp)
  else
    allocate(unscaled_field, source=field)
    unscaled_field(:) = field_md%conversion * field(:)
    call write_field(handle%handle_infra, field_infra, unscaled_field, tstamp=tstamp)
    deallocate(unscaled_field)
  endif
end procedure write_field_1d_infra
module procedure write_field_0d_infra
  type(fieldtype) :: field_infra
  real :: unscaled_field ! An unscaled version of field for output [a]
  field_infra = handle%fields%get(field_md%label)
  unscaled_field = field_md%conversion*field
  call write_field(handle%handle_infra, field_infra, unscaled_field, tstamp=tstamp)
end procedure write_field_0d_infra
module procedure write_field_axis_infra
  type(axistype) :: axis_infra
  axis_infra = handle%axes%get(axis%label)
  call write_field(handle%handle_infra, axis_infra)
end procedure write_field_axis_infra
module procedure write_attribute_infra
  call write_metadata(handle%handle_infra, name, attribute)
end procedure write_attribute_infra
module procedure file_is_open_infra
  file_is_open_infra = fms2_file_is_open(handle%handle_infra)
end procedure file_is_open_infra
module procedure get_file_info_infra
  call get_file_info(handle%handle_infra, ndim, nvar, ntime)
end procedure get_file_info_infra
module procedure get_file_fields_infra
  type(fieldtype), allocatable :: fields_infra(:)
  integer :: i
  character(len=64) :: label
  allocate(fields_infra(size(fields)))
  call get_file_fields(handle%handle_infra, fields_infra)

  do i = 1, size(fields)
    call get_field_atts(fields_infra(i), name=label)
    call handle%fields%append(fields_infra(i), trim(label))
    fields(i)%label = trim(label)
  enddo
end procedure get_file_fields_infra
module procedure get_file_times_infra
  call get_file_times(handle%handle_infra, time_values, ntime=ntime)
end procedure get_file_times_infra
module procedure get_field_atts_infra
  type(fieldtype) :: field_infra
  field_infra = handle%fields%get(field%label)
  call get_field_atts(field_infra, name, units, longname, checksum)
end procedure get_field_atts_infra
module procedure read_field_chksum_infra
  type(fieldtype) :: field_infra
  field_infra = handle%fields%get(field%label)
  call read_field_chksum(field_infra, chksum, valid_chksum)
end procedure read_field_chksum_infra
module procedure get_file_fieldtypes
  type(field_node_infra), pointer :: node
  integer :: i
  node => handle%fields%head
  do i = 1, size(fields)
    if (.not. associated(node%next)) &
      call MOM_error(FATAL, 'fields(:) size exceeds number of registered fields.')
    fields(i) = node%field
    node => node%next
  enddo
end procedure get_file_fieldtypes
module procedure open_file_nc
  if (.not. present(MOM_domain) .and. .not. is_root_PE()) return

  call open_netcdf_file(handle%handle_nc, filename, action)
  handle%is_open = .true.

  if (present(MOM_domain)) then
    handle%domain_decomposed = .true.

    ! Input files use unrotated indexing.
    if (associated(MOM_domain%domain_in)) then
      call hor_index_init(MOM_domain%domain_in, handle%HI)
    else
      call hor_index_init(MOM_domain, handle%HI)
    endif
  endif

  call handle%axes%init()
  call handle%fields%init()
end procedure open_file_nc
module procedure close_file_nc
  if (.not. handle%domain_decomposed .and. .not. is_root_PE()) return

  handle%is_open = .false.
  call close_netcdf_file(handle%handle_nc)
end procedure close_file_nc
module procedure flush_file_nc
  if (.not. is_root_PE()) return

  call flush_netcdf_file(handle%handle_nc)
end procedure flush_file_nc
module procedure register_axis_nc
  type(netcdf_axis) :: axis_nc
  if (is_root_PE()) then
    axis_nc = register_netcdf_axis(handle%handle_nc, label, units, longname, &
        data, cartesian, sense)

    call handle%axes%append(axis_nc, label)
  endif
  axis%label = label
end procedure register_axis_nc
module procedure register_field_nc
  type(netcdf_field) :: field_nc
  type(netcdf_axis), allocatable :: axes_nc(:)
  integer :: i
  if (is_root_PE()) then
    allocate(axes_nc(size(axes)))
    do i = 1, size(axes)
      axes_nc(i) = handle%axes%get(axes(i)%label)
    enddo

    field_nc = register_netcdf_field(handle%handle_nc, label, axes_nc, longname, units)

    call handle%fields%append(field_nc, label)
  endif
  field%label = label
  field%conversion = 1.0 ; if (present(conversion)) field%conversion = conversion
end procedure register_field_nc
module procedure write_attribute_nc
  if (.not. is_root_PE()) return

  call write_netcdf_attribute(handle%handle_nc, name, attribute)
end procedure write_attribute_nc
module procedure write_field_4d_nc
  type(netcdf_field) :: field_nc
  real, allocatable :: unscaled_field(:,:,:,:) ! An unscaled version of field for output [a]
  if (.not. is_root_PE()) return

  field_nc = handle%fields%get(field_md%label)
  if (field_md%conversion == 1.0) then
    call write_netcdf_field(handle%handle_nc, field_nc, field, time=tstamp)
  else
    allocate(unscaled_field, source=field)
    unscaled_field(:,:,:,:) = field_md%conversion * field(:,:,:,:)
    call write_netcdf_field(handle%handle_nc, field_nc, unscaled_field, time=tstamp)
    deallocate(unscaled_field)
  endif
end procedure write_field_4d_nc
module procedure write_field_3d_nc
  type(netcdf_field) :: field_nc
  real, allocatable :: unscaled_field(:,:,:) ! An unscaled version of field for output [a]
  if (.not. is_root_PE()) return

  field_nc = handle%fields%get(field_md%label)
  if (field_md%conversion == 1.0) then
    call write_netcdf_field(handle%handle_nc, field_nc, field, time=tstamp)
  else
    allocate(unscaled_field, source=field)
    unscaled_field(:,:,:) = field_md%conversion * field(:,:,:)
    call write_netcdf_field(handle%handle_nc, field_nc, unscaled_field, time=tstamp)
    deallocate(unscaled_field)
  endif
end procedure write_field_3d_nc
module procedure write_field_2d_nc
  type(netcdf_field) :: field_nc
  real, allocatable :: unscaled_field(:,:) ! An unscaled version of field for output [a]
  if (.not. is_root_PE()) return

  field_nc = handle%fields%get(field_md%label)
  if (field_md%conversion == 1.0) then
    call write_netcdf_field(handle%handle_nc, field_nc, field, time=tstamp)
  else
    allocate(unscaled_field, source=field)
    unscaled_field(:,:) = field_md%conversion * field(:,:)
    call write_netcdf_field(handle%handle_nc, field_nc, unscaled_field, time=tstamp)
    deallocate(unscaled_field)
  endif
end procedure write_field_2d_nc
module procedure write_field_1d_nc
  type(netcdf_field) :: field_nc
  real, allocatable :: unscaled_field(:) ! An unscaled version of field for output [a]
  if (.not. is_root_PE()) return

  field_nc = handle%fields%get(field_md%label)
  if (field_md%conversion == 1.0) then
    call write_netcdf_field(handle%handle_nc, field_nc, field, time=tstamp)
  else
    allocate(unscaled_field, source=field)
    unscaled_field(:) = field_md%conversion * field(:)
    call write_netcdf_field(handle%handle_nc, field_nc, unscaled_field, time=tstamp)
    deallocate(unscaled_field)
  endif
end procedure write_field_1d_nc
module procedure write_field_0d_nc
  type(netcdf_field) :: field_nc
  real :: unscaled_field ! An unscaled version of field for output [a]
  if (.not. is_root_PE()) return

  field_nc = handle%fields%get(field_md%label)
  unscaled_field = field_md%conversion * field
  call write_netcdf_field(handle%handle_nc, field_nc, unscaled_field, time=tstamp)
end procedure write_field_0d_nc
module procedure write_field_axis_nc
  type(netcdf_axis) :: axis_nc
  if (.not. is_root_PE()) return

  axis_nc = handle%axes%get(axis%label)
  call write_netcdf_axis(handle%handle_nc, axis_nc)
end procedure write_field_axis_nc
module procedure file_is_open_nc
  file_is_open_nc = handle%is_open
end procedure file_is_open_nc
module procedure get_file_info_nc
  integer :: ndim_nc, nvar_nc
  if (.not. is_root_PE()) return

  call get_netcdf_size(handle%handle_nc, ndims=ndim_nc, nvars=nvar_nc, nsteps=ntime)

  ! MOM I/O follows legacy FMS behavior and excludes axes from field count
  if (present(ndim)) ndim = ndim_nc
  if (present(nvar)) nvar = nvar_nc - ndim_nc
end procedure get_file_info_nc
module procedure update_file_contents_nc
  type(netcdf_axis), allocatable :: axes_nc(:)
  type(netcdf_field), allocatable :: fields_nc(:)
  integer :: i
  if (.not. handle%domain_decomposed .and. .not. is_root_PE()) return

  call get_netcdf_fields(handle%handle_nc, axes_nc, fields_nc)

  do i = 1, size(axes_nc)
    call handle%axes%append(axes_nc(i), axes_nc(i)%label)
  enddo

  do i = 1, size(fields_nc)
    call handle%fields%append(fields_nc(i), fields_nc(i)%label)
  enddo
end procedure update_file_contents_nc
module procedure get_file_fields_nc
  type(field_node_nc), pointer :: node => null()
  integer :: n
  if (.not. is_root_PE()) return

  ! Generate the manifest of axes and fields
  call handle%update()

  n = 0
  node => handle%fields%head
  do while (associated(node%next))
    n = n + 1
    fields(n)%label = trim(node%label)
    node => node%next
  enddo
end procedure get_file_fields_nc
module procedure get_field_atts_nc
  call MOM_error(FATAL, 'get_field_atts over netCDF is not yet implemented.')
end procedure get_field_atts_nc
module procedure read_field_chksum_nc
  call MOM_error(FATAL, 'read_field_chksum over netCDF is not yet implemented.')
  chksum = -1_int64
  valid_chksum = .false.
end procedure read_field_chksum_nc
module procedure get_field_nc
  logical :: data_domain
  logical :: compute_domain
  type(netcdf_field) :: field_nc
  integer :: isc, iec, jsc, jec
  integer :: isd, ied, jsd, jed
  integer :: iscl, iecl, jscl, jecl
  integer :: bounds(2,2)
  real, allocatable :: values_c(:,:)
  isc = handle%HI%isc
  iec = handle%HI%iec
  jsc = handle%HI%jsc
  jec = handle%HI%jec

  isd = handle%HI%isd
  ied = handle%HI%ied
  jsd = handle%HI%jsd
  jed = handle%HI%jed

  data_domain = all(shape(values) == [ied-isd+1, jed-jsd+1])
  compute_domain = all(shape(values) == [iec-isc+1, jec-jsc+1])

  ! NOTE: Data on face and vertex points is not yet supported.  This is a
  ! temporary check to detect such cases, but may be removed in the future.
  if (.not. (compute_domain .or. data_domain)) &
    call MOM_error(FATAL, 'get_field_nc trying to read '//trim(label)//' from '//&
                   trim(get_netcdf_filename(handle%handle_nc))//&
                   ': Only compute and data domains are currently supported.')

  field_nc = handle%fields%get(label)

  if (data_domain) &
    allocate(values_c(1:iec-isc+1,1:jec-jsc+1))

  if (handle%domain_decomposed) then
    bounds(1,:) = [isc, jsc] + [handle%HI%idg_offset, handle%HI%jdg_offset]
    bounds(2,:) = [iec, jec] + [handle%HI%idg_offset, handle%HI%jdg_offset]
    if (data_domain) then
      call read_netcdf_field(handle%handle_nc, field_nc, values_c, bounds=bounds)
    else
      call read_netcdf_field(handle%handle_nc, field_nc, values, bounds=bounds)
    endif
  else
    if (data_domain) then
      call read_netcdf_field(handle%handle_nc, field_nc, values_c)
    else
      call read_netcdf_field(handle%handle_nc, field_nc, values)
    endif
  endif

  if (data_domain) then
    iscl = isc - isd + 1
    iecl = iec - isd + 1
    jscl = jsc - jsd + 1
    jecl = jec - jsd + 1

    values(iscl:iecl,jscl:jecl) = values_c(:,:)
  else
    iscl = 1
    iecl = iec - isc + 1
    jscl = 1
    jecl = jec - jsc + 1
  endif

  ! NOTE: It is more efficient to do the rescale in-place while copying
  ! values_c(:,:) to values(:,:).  But since rescale is only present for
  ! debugging, we can probably disregard this impact on performance.
  if (present(rescale)) then
    if (rescale /= 1.0) then
      values(iscl:iecl,jscl:jecl) = rescale * values(iscl:iecl,jscl:jecl)
    endif
  endif
end procedure get_field_nc
end submodule MOM_io_file_s
