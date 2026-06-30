submodule (MOM_interp_infra) MOM_interp_infra_s
  implicit none
contains
module procedure horizontal_interp_init
  call horiz_interp_init()
end procedure horizontal_interp_init
module procedure time_interp_extern_init
  call time_interp_external_init()
end procedure time_interp_extern_init
module procedure horiz_interp_from_weights_field2d
  call horiz_interp(Interp, data_in, data_out, verbose, &
                    mask_in, mask_out, missing_value, missing_permit, &
                    err_msg, new_missing_handle=.true. )

end procedure horiz_interp_from_weights_field2d
module procedure horiz_interp_from_weights_field3d
  call horiz_interp(Interp, data_in, data_out, verbose, mask_in, mask_out, &
                    missing_value, missing_permit, err_msg)

end procedure horiz_interp_from_weights_field3d
module procedure build_horiz_interp_weights_2d_to_2d
  call horiz_interp_new(Interp, lon_in, lat_in, lon_out, lat_out, &
                        verbose, interp_method, num_nbrs, max_dist, &
                        src_modulo, mask_in, mask_out, &
                        is_latlon_in, is_latlon_out)

end procedure build_horiz_interp_weights_2d_to_2d
module procedure get_extern_field_size
  get_extern_field_size = get_external_field_size(index)

end procedure get_extern_field_size
module procedure get_extern_field_axes
  integer :: ndims
  integer, allocatable :: dims(:)
  character(len=256) :: dim_name
  integer :: dim_len
  integer :: var_dim
  real, allocatable :: axis_points(:)
  integer :: ncid
  integer :: varid
  integer :: rc
  integer :: nc_start(1)
  integer :: nc_count(1)
  integer :: d
  character(len=2) :: d_str
  rc = nf90_open(trim(field%filename), NF90_NOWRITE, ncid)
  if (rc /= NF90_NOERR) &
    call MOM_error(FATAL, "Error opening file " // trim(field%filename) // ".")

  ! Use field%label to get the netCDF varid
  rc = nf90_inq_varid(ncid, trim(field%label), varid)
  if (rc /= NF90_NOERR) &
    call MOM_error(FATAL, "Error finding variable " // trim(field%label) &
        // " in " // trim(field%filename) // ".")

  ! Use the varid to get the number of dims (ndims) and their IDs (dims(:))
  !   Verify that ndims >= 3
  rc = nf90_inquire_variable(ncid, varid, ndims=ndims)
  if (rc /= NF90_NOERR) &
    call MOM_error(FATAL, "Error querying variable " // trim(field%label) &
        // " in " // trim(field%filename) // ".")

  if (ndims < 3) &
    call MOM_error(FATAL, trim(field%label) // " in " // trim(field%filename) &
        // " has too few dimensions to be read as a 3D array.")

  allocate(dims(ndims))

  rc = nf90_inquire_variable(ncid, varid, dimids=dims)
  if (rc /= NF90_NOERR) &
    call MOM_error(FATAL, "Error querying variable " // trim(field%label) &
        // " in " // trim(field%filename) // ".")

  do d=1,ndims
    ! Determine the name of each dimension
    rc = nf90_inquire_dimension(ncid, dims(d), dim_name, len=dim_len)
    if (rc /= NF90_NOERR) then
      write(d_str, '(i0)') d
      call MOM_error(FATAL, "Error querying dimension " // trim(d_str) &
          // " of " // trim(field%label) // " in " // trim(field%filename) &
          // ".")
    endif

    ! Now locate a variable with the same name as the dimension (e.g. "x")
    rc = nf90_inq_varid(ncid, dim_name, var_dim)
    if (rc /= NF90_NOERR) &
      call MOM_error(FATAL, "Error finding dimension variable " &
          // trim(dim_name) // " of " // trim(field%label) // " in " &
          // trim(field%filename))

    allocate(axis_points(dim_len))

    ! Get the dimensional axis values
    nc_start(1) = 1
    nc_count(1) = dim_len
    rc = nf90_get_var(ncid, var_dim, axis_points, nc_start, nc_count)
    if (rc /= NF90_NOERR) &
      call MOM_error(FATAL, "Error reading dimension " // trim(dim_name) &
          // " axis data of " // trim(field%label) // " in " &
          // trim(field%filename))

    ! write via set_axis_info() equivalent for axistype
    call set_axis_data(axes(d), dim_name, axis_points)

    deallocate(axis_points)
  enddo

  deallocate(dims)

  ! Close external file
  rc = nf90_close(ncid)
  if (rc /= NF90_NOERR) &
    call MOM_error(FATAL, "Error closing file "//trim(field%filename)//".")

end procedure get_extern_field_axes
module procedure get_extern_field_missing
  get_extern_field_missing = get_external_field_missing(index)

end procedure get_extern_field_missing
module procedure get_external_field_info
  if (present(size)) then
    size(:) = get_extern_field_size(field%id)
  endif

  if (present(axes)) then
    axes(:) = get_extern_field_axes(field)
  endif

  if (present(missing)) then
    missing = get_extern_field_missing(field%id)
  endif
end procedure get_external_field_info
module procedure time_interp_extern_0d
  call time_interp_external(field%id, time, data_in, verbose=verbose)
end procedure time_interp_extern_0d
module procedure time_interp_extern_2d
  call time_interp_external(field%id, time, data_in, interp=interp, verbose=verbose, &
                            horz_interp=horz_interp, mask_out=mask_out)
end procedure time_interp_extern_2d
module procedure time_interp_extern_3d
  call time_interp_external(field%id, time, data_in, interp=interp, verbose=verbose, &
                            horz_interp=horz_interp, mask_out=mask_out)
end procedure time_interp_extern_3d
module procedure init_extern_field
  type(FmsNetcdfFile_t) :: extern_file
  integer :: num_fields
  character(len=256), allocatable :: extern_fieldnames(:)
  character(len=:), allocatable :: label
  logical :: rc
  integer :: i
  field%filename = file

  ! FMS2's init_external_field is case sensitive, so we must replicate the
  !   case-insensitivity of FMS1.  This requires opening the file twice.

  rc = netcdf_file_open(extern_file, file, 'read')
  if (.not. rc) then
    call MOM_error(FATAL, 'init_extern_file: file ' // trim(file) &
        // ' could not be opened.')
  endif

  ! TODO: broadcast = .false.?
  num_fields = get_num_variables(extern_file)

  allocate(extern_fieldnames(num_fields))
  call get_variable_names(extern_file, extern_fieldnames)

  do i = 1, num_fields
    if (lowercase(extern_fieldnames(i)) == lowercase(fieldname)) then
      field%label = extern_fieldnames(i)
      exit
    endif
  enddo

  call netcdf_file_close(extern_file)

  if (.not. allocated(field%label)) then
    call MOM_error(FATAL, 'init_extern_field: field ' // trim(fieldname) &
        // ' not found in ' // trim(file) // '.')
  endif

  ! Pass to FMS2 implementation of init_external_field

  ! NOTE: external fields are currently assumed to be on-grid, which holds
  ! across the current codebase.  In the future, we may need to either enforce
  ! this or somehow relax this requirement.

  if (present(MOM_Domain)) then
    field%id = init_external_field(file, field%label, domain=MOM_domain%mpp_domain, &
             verbose=verbose, ierr=ierr, ignore_axis_atts=ignore_axis_atts, &
             correct_leap_year_inconsistency=correct_leap_year_inconsistency, &
             ongrid=.true.)
  else
    field%id = init_external_field(file, field%label, domain=domain, &
             verbose=verbose, ierr=ierr, ignore_axis_atts=ignore_axis_atts, &
             correct_leap_year_inconsistency=correct_leap_year_inconsistency, &
             ongrid=.true.)
  endif
end procedure init_extern_field
end submodule MOM_interp_infra_s
