submodule (MOM_IS_diag_mediator) MOM_IS_diag_mediator_s
#define DIAG_ALLOC_CHUNK_SIZE 15
  implicit none
contains
module procedure set_IS_axes_info
  integer :: id_xq, id_yq, id_xh, id_yh, id_null
  integer :: i, j
  character(len=80) :: set_name
  real, allocatable, dimension(:) :: IaxB, iax ! Index-based integer and half-integer i-axis labels [nondim]
  real, allocatable, dimension(:) :: JaxB, jax ! Index-based integer and half-integer j-axis labels [nondim]
  set_name = "ice_shelf" ; if (present(axes_set_name)) set_name = trim(axes_set_name)

  if (diag_cs%index_space_axes) then
    allocate(IaxB(G%IsgB:G%IegB))
    do I=G%IsgB,G%IegB
      Iaxb(I) = real(I)
    enddo
    allocate(iax(G%isg:G%ieg))
    do i=G%isg,G%ieg
      iax(i) = real(i)-0.5
    enddo
    allocate(JaxB(G%JsgB:G%JegB))
    do J=G%JsgB,G%JegB
      JaxB(J) = real(J)
    enddo
    allocate(jax(G%jsg:G%jeg))
    do j=G%jsg,G%jeg
      jax(j) = real(j)-0.5
    enddo
  endif

  ! Horizontal axes for the native grids.
  if (diag_cs%index_space_axes) then
    if (G%symmetric) then
      id_xq = MOM_diag_axis_init('Iq', IaxB(G%IsgB:G%IegB), 'none', 'x', &
          'Boundary (q) point grid-space longitude', G%Domain, position=EAST, set_name=set_name)
      id_yq = MOM_diag_axis_init('Jq', JaxB(G%JsgB:G%JegB), 'none', 'y', &
          'Boundary (q) point grid-space latitude', G%Domain, position=NORTH, set_name=set_name)
    else
      id_xq = MOM_diag_axis_init('Iq', IaxB(G%isg:G%ieg), 'none', 'x', &
          'Boundary (q) point grid-space longitude', G%Domain, position=EAST, set_name=set_name)
      id_yq = MOM_diag_axis_init('Jq', JaxB(G%jsg:G%jeg), 'none', 'y', &
          'Boundary (q) point grid-space latitude', G%Domain, position=NORTH, set_name=set_name)
    endif

    id_xh = MOM_diag_axis_init('ih', iax, 'none', 'x', &
        'Tracer (h) point grid-space longitude', G%Domain, set_name=set_name)
    id_yh = MOM_diag_axis_init('jh', jax, 'none', 'y', &
        'Tracer (h) point grid-space latitude', G%Domain, set_name=set_name)
  else
    if (G%symmetric) then
      id_xq = MOM_diag_axis_init('xB', G%gridLonB(G%isgB:G%iegB), G%x_axis_units, 'x', &
          'Boundary point nominal longitude', G%Domain, position=EAST, set_name=set_name)
      id_yq = MOM_diag_axis_init('yB', G%gridLatB(G%jsgB:G%jegB), G%y_axis_units, 'y', &
          'Boundary point nominal latitude', G%Domain, position=NORTH, set_name=set_name)

    else
      id_xq = MOM_diag_axis_init('xB', G%gridLonB(G%isg:G%ieg), G%x_axis_units, 'x', &
          'Boundary point nominal longitude', G%Domain, position=EAST, set_name=set_name)
      id_yq = MOM_diag_axis_init('yB', G%gridLatB(G%jsg:G%jeg), G%y_axis_units, 'y', &
          'Boundary point nominal latitude', G%Domain, position=NORTH, set_name=set_name)

    endif
    id_xh = MOM_diag_axis_init('xT', G%gridLonT(G%isg:G%ieg), G%x_axis_units, 'x', &
        'Tracer point nominal longitude', G%Domain, set_name=set_name)
    id_yh = MOM_diag_axis_init('yT', G%gridLatT(G%jsg:G%jeg), G%y_axis_units, 'y', &
        'Tracer point nominal latitude', G%Domain, set_name=set_name)
  endif

  ! Axis groupings for 2-D arrays
  call define_axes_group(diag_cs, (/id_xh, id_yh/), diag_cs%axesT1, &
       x_cell_method='mean', y_cell_method='mean', is_h_point=.true.)
  call define_axes_group(diag_cs, (/id_xq, id_yq/), diag_cs%axesB1, &
       x_cell_method='point', y_cell_method='point', is_q_point=.true.)
  call define_axes_group(diag_cs, (/id_xq, id_yh/), diag_cs%axesCu1, &
       x_cell_method='point', y_cell_method='mean', is_u_point=.true.)
  call define_axes_group(diag_cs, (/id_xh, id_yq/), diag_cs%axesCv1, &
       x_cell_method='mean', y_cell_method='point', is_v_point=.true.)

  ! Axis group for special null axis for scalars from diag manager.
  id_null = MOM_diag_axis_init('scalar_axis', (/0./), 'none', 'N', 'none', null_axis=.true.)
  call define_axes_group(diag_cs, (/ id_null /), diag_cs%axesNull)

  if (diag_cs%index_space_axes) then
    deallocate(IaxB, iax, JaxB, jax)
  endif

end procedure set_IS_axes_info
module procedure diag_register_area_ids
  integer :: fms_id, i
  if (present(id_area_t)) then
    fms_id = diag_cs%diags(id_area_t)%fms_diag_id
    diag_cs%axesT1%id_area = fms_id
  endif
  if (present(id_area_q)) then
    fms_id = diag_cs%diags(id_area_q)%fms_diag_id
    diag_cs%axesB1%id_area = fms_id
  endif
end procedure diag_register_area_ids
module procedure define_axes_group
  integer :: n
  n = size(handles)
  if (n<1 .or. n>2) call MOM_error(FATAL, "define_axes_group: wrong size for list of handles!")
  allocate( axes%handles(n) )
  axes%id = ints_to_string(handles, max(n,2)) ! Identifying string
  axes%rank = n
  axes%handles(:) = handles(:)
  axes%diag_cs => diag_cs ! A (circular) link back to the diag_ctrl structure

  if ((axes%rank<2) .and. (present(x_cell_method) .or. present(x_cell_method))) &
    call MOM_error(FATAL, 'define_axes_group: Can not set x_cell_method or y_cell_method for rank<2.')
  axes%x_cell_method = '' ; if (present(x_cell_method)) axes%x_cell_method = trim(x_cell_method)
  axes%y_cell_method = '' ; if (present(y_cell_method)) axes%y_cell_method = trim(y_cell_method)

  if (present(is_h_point)) axes%is_h_point = is_h_point
  if (present(is_q_point)) axes%is_q_point = is_q_point
  if (present(is_u_point)) axes%is_u_point = is_u_point
  if (present(is_v_point)) axes%is_v_point = is_v_point

  ! Setup masks for this axes group
  axes%mask2d => null()
  if (axes%rank==2) then
    if (axes%is_h_point) axes%mask2d => diag_cs%mask2dT
    if (axes%is_h_point) axes%mask2d_comp => diag_cs%mask2dT_comp
    if (axes%is_u_point) axes%mask2d => diag_cs%mask2dCu
    if (axes%is_v_point) axes%mask2d => diag_cs%mask2dCv
    if (axes%is_q_point) axes%mask2d => diag_cs%mask2dBu
  endif

end procedure define_axes_group
module procedure set_IS_diag_mediator_grid
  diag_cs%is = G%isc - (G%isd-1) ; diag_cs%ie = G%iec - (G%isd-1)
  diag_cs%js = G%jsc - (G%jsd-1) ; diag_cs%je = G%jec - (G%jsd-1)
  diag_cs%isd = G%isd ; diag_cs%ied = G%ied
  diag_cs%jsd = G%jsd ; diag_cs%jed = G%jed

end procedure set_IS_diag_mediator_grid
module procedure post_IS_data_0d
  real :: locfield ! The field being offered in arbitrary unscaled units [a]
  logical :: used, is_stat
  type(diag_type), pointer :: diag => null()
  integer :: time_days
  integer :: time_seconds
  character(len=300) :: debug_mesg
  if (id_clock_diag_mediator>0) call cpu_clock_begin(id_clock_diag_mediator)
  is_stat = .false. ; if (present(is_static)) is_stat = is_static

  ! Iterate over list of diag 'variants', e.g. CMOR aliases, call send_data
  ! for each one.
  call assert(diag_field_id < diag_cs%next_free_diag_id, &
              'post_IS_data_0d: Unregistered diagnostic id')
  diag => diag_cs%diags(diag_field_id)

  do while (associated(diag))
    locfield = field
    if (diag%conversion_factor /= 0.) &
      locfield = locfield * diag%conversion_factor

    if (diag_cs%diag_as_chksum) then
      ! Append timestep to mesg
      call get_time(diag_cs%time_end, time_seconds, days=time_days)
      write(debug_mesg, '(a, 1x, i0, 1x, i0)') &
          trim(diag%debug_str), time_days, time_seconds

      call chksum0(locfield, debug_mesg, logunit=diag_cs%chksum_iounit)
    elseif (is_stat) then
      used = send_data_infra(diag%fms_diag_id, locfield)
    elseif (diag_cs%ave_enabled) then
      used = send_data_infra(diag%fms_diag_id, locfield, diag_cs%time_end)
    endif

    diag => diag%next
  enddo

  if (id_clock_diag_mediator>0) call cpu_clock_end(id_clock_diag_mediator)
end procedure post_IS_data_0d
module procedure post_IS_data_2d
  type(diag_type), pointer :: diag => NULL()
  if (id_clock_diag_mediator>0) call cpu_clock_begin(id_clock_diag_mediator)

  ! Iterate over list of diag 'variants' (e.g. CMOR aliases) and post each.
  call assert(diag_field_id < diag_cs%next_free_diag_id, &
              'post_IS_data_2d: Unregistered diagnostic id')
  diag => diag_cs%diags(diag_field_id)
  do while (associated(diag))
    call post_data_2d_low(diag, field, diag_cs, is_static, mask)
    diag => diag%next
  enddo

  if (id_clock_diag_mediator>0) call cpu_clock_end(id_clock_diag_mediator)
end procedure post_IS_data_2d
module procedure post_data_2d_low
  real, dimension(:,:), pointer :: locfield ! The field being offered in arbitrary unscaled units [a]
  real, dimension(:,:), pointer :: locmask  ! A pointer to the data mask to use [nondim]
  logical :: used  ! The return value of send_data is not used for anything.
  logical :: is_stat
  logical :: i_data, j_data ! True if the field is on the data domain in the i or j directions.
  integer :: cszi, cszj, dszi, dszj
  integer :: isv, iev, jsv, jev, i, j
  integer :: time_days, time_seconds
  character(len=300) :: mesg
  character(len=300) :: debug_mesg
  locfield => NULL()
  locmask => NULL()
  is_stat = .false. ; if (present(is_static)) is_stat = is_static

  ! Determine the proper array indices, noting that because of the (:,:)
  ! declaration of field, symmetric arrays are using a SW-grid indexing,
  ! but non-symmetric arrays are using a NE-grid indexing.  Send_data
  ! actually only uses the difference between ie and is to determine
  ! the output data size and assumes that halos are symmetric.
  isv = diag_cs%is ; iev = diag_cs%ie ; jsv = diag_cs%js ; jev = diag_cs%je

  cszi = (diag_cs%ie-diag_cs%is) +1 ; dszi = (diag_cs%ied-diag_cs%isd) +1
  cszj = (diag_cs%je-diag_cs%js) +1 ; dszj = (diag_cs%jed-diag_cs%jsd) +1
  if ( size(field,1) == dszi ) then
    isv = diag_cs%is ; iev = diag_cs%ie ; i_data = .true.   ! Data domain
  elseif ( size(field,1) == dszi + 1 ) then
    isv = diag_cs%is ; iev = diag_cs%ie+1 ; i_data = .true. ! Symmetric data domain
  elseif ( size(field,1) == cszi ) then
    isv = 1 ; iev = cszi ; i_data = .false. ! Computational domain
  elseif ( size(field,1) == cszi + 1 ) then
    isv = 1 ; iev = cszi+1 ; i_data = .false. ! Symmetric computational domain
  else
    write (mesg,*) " peculiar size ",size(field,1)," in i-direction\n"//&
       "does not match one of ", cszi, cszi+1, dszi, dszi+1
    call MOM_error(FATAL,"post_IS_data_2d_low: "//trim(diag%debug_str)//trim(mesg))
  endif

  if ( size(field,2) == dszj ) then
    jsv = diag_cs%js ; jev = diag_cs%je ; j_data = .true.   ! Data domain
  elseif ( size(field,2) == dszj + 1 ) then
    jsv = diag_cs%js ; jev = diag_cs%je+1 ; j_data = .true. ! Symmetric data domain
  elseif ( size(field,2) == cszj ) then
    jsv = 1 ; jev = cszj ; j_data = .false. ! Computational domain
  ! This was: elseif ( size(field,1) == cszj + 1 ) then
  elseif ( size(field,2) == cszj + 1 ) then
    jsv = 1 ; jev = cszj+1 ; j_data = .false. ! Symmetric computational domain
  else
    write (mesg,*) " peculiar size ",size(field,2)," in j-direction\n"//&
       "does not match one of ", cszj, cszj+1, dszj, dszj+1
    call MOM_error(FATAL,"post_IS_data_2d_low: "//trim(diag%debug_str)//trim(mesg))
  endif

  if ((diag%conversion_factor /= 0.) .and. (diag%conversion_factor /= 1.)) then
    allocate( locfield( lbound(field,1):ubound(field,1), lbound(field,2):ubound(field,2) ) )
    do j=jsv,jev ; do i=isv,iev
      if (field(i,j) == diag_cs%missing_value) then
        locfield(i,j) = diag_cs%missing_value
      else
        locfield(i,j) = field(i,j) * diag%conversion_factor
      endif
    enddo ; enddo
    locfield(isv:iev,jsv:jev) = field(isv:iev,jsv:jev) * diag%conversion_factor
  else
    locfield => field
  endif

  ! Handle cases where the data and computational domain are the same size.
  if (diag_cs%ied-diag_cs%isd == diag_cs%ie-diag_cs%is) i_data = j_data
  if (diag_cs%jed-diag_cs%jsd == diag_cs%je-diag_cs%js) j_data = i_data
  if ( i_data .NEQV. j_data ) then
    call MOM_error(FATAL, "post_IS_data_2d: post_IS_data called for "//&
                   trim(diag%debug_str)//" with mixed computational and data domain array sizes.")
  endif

  if (present(mask)) then
    locmask => mask
  elseif (.not.is_stat) then  ! Static fields do not have assigned axes.
    if (i_data .and. associated(diag%axes%mask2d)) then
      locmask => diag%axes%mask2d
    elseif ((.not.i_data) .and. associated(diag%axes%mask2d_comp)) then
      locmask => diag%axes%mask2d_comp
    endif
  endif
  if (associated(locmask)) call assert(size(locfield) == size(locmask), &
        'post_data_2d_low: mask size mismatch: '//trim(diag%debug_str))

  if (diag_cs%diag_as_chksum) then
    ! Append timestep to mesg
    call get_time(diag_cs%time_end, time_seconds, days=time_days)
    write(debug_mesg, '(a, 1x, i0, 1x, i0)') &
        trim(diag%debug_str), time_days, time_seconds

    if (diag%axes%is_h_point) then
      call hchksum(locfield, debug_mesg, diag_cs%G%HI, &
                   logunit=diag_cs%chksum_iounit)
    elseif (diag%axes%is_u_point) then
      call uchksum(locfield, debug_mesg, diag_cs%G%HI, &
                   logunit=diag_cs%chksum_iounit)
    elseif (diag%axes%is_v_point) then
      call vchksum(locfield, debug_mesg, diag_cs%G%HI, &
                   logunit=diag_cs%chksum_iounit)
    elseif (diag%axes%is_q_point) then
      call Bchksum(locfield, debug_mesg, diag_cs%G%HI, &
                   logunit=diag_cs%chksum_iounit)
    else
      call MOM_error(FATAL, "post_data_2d_low: unknown axis type.")
    endif
  else
    if (is_stat) then
      if (associated(locmask)) then
        used = send_data_infra(diag%fms_diag_id, locfield, &
                         is_in=isv, ie_in=iev, js_in=jsv, je_in=jev, rmask=locmask)
      else
        used = send_data_infra(diag%fms_diag_id, locfield, &
                         is_in=isv, ie_in=iev, js_in=jsv, je_in=jev)
      endif
    elseif (diag_cs%ave_enabled) then
      if (associated(locmask)) then
        used = send_data_infra(diag%fms_diag_id, locfield, &
                         is_in=isv, ie_in=iev, js_in=jsv, je_in=jev, &
                         time=diag_cs%time_end, weight=diag_cs%time_int, rmask=locmask)
      else
        used = send_data_infra(diag%fms_diag_id, locfield, &
                         is_in=isv, ie_in=iev, js_in=jsv, je_in=jev, &
                         time=diag_cs%time_end, weight=diag_cs%time_int)
      endif
    endif
  endif

  if ((diag%conversion_factor /= 0.) .and. (diag%conversion_factor /= 1.)) deallocate( locfield )

end procedure post_data_2d_low
module procedure enable_averaging
  diag_cs%time_int = time_int_in
  diag_cs%time_end = time_end_in
  diag_cs%ave_enabled = .true.
end procedure enable_averaging
module procedure enable_averages
  if (present(T_to_s)) then
    diag_cs%time_int = time_int*T_to_s
  elseif (associated(diag_CS%US)) then
    diag_cs%time_int = time_int*diag_CS%US%T_to_s
  else
    diag_cs%time_int = time_int
  endif
  diag_cs%time_end = time_end
  diag_cs%ave_enabled = .true.
end procedure enable_averages
module procedure disable_averaging
  diag_cs%time_int = 0.0
  diag_cs%ave_enabled = .false.
end procedure disable_averaging
module procedure query_averaging_enabled
  if (present(time_int)) time_int = diag_cs%time_int
  if (present(time_end)) time_end = diag_cs%time_end
  query_averaging_enabled = diag_cs%ave_enabled
end procedure query_averaging_enabled
module procedure MOM_IS_diag_mediator_infrastructure_init
  call MOM_diag_manager_init(err_msg=err_msg)
end procedure MOM_IS_diag_mediator_infrastructure_init
module procedure get_diag_time_end
  get_diag_time_end = diag_cs%time_end
end procedure get_diag_time_end
module procedure register_MOM_IS_diag_field
  real :: MOM_missing_value ! A value used to indicate missing values in output files, in arbitrary units [a]
  type(diag_ctrl), pointer :: diag_cs => NULL() ! A structure that is used to regulate diagnostic output
  type(axes_grp), pointer :: axes
  integer :: dm_id
  character(len=256) :: msg
  character(len=256) :: cm_string ! A string describing the cell methods returned from attach_cell_methods.
  character(len=256) :: new_module_name
  character(len=480) :: module_list, var_list
  character(len=24)  :: dimensions
  integer :: num_modnm, num_varnm
  logical :: active
  diag_cs => axes_in%diag_cs

  ! Check if the axes match a standard grid axis.
  ! If not, allocate the new axis and copy the contents.
  if (axes_in%id == diag_cs%axesT1%id) then
    axes => diag_cs%axesT1
  elseif (axes_in%id == diag_cs%axesB1%id) then
    axes => diag_cs%axesB1
  elseif (axes_in%id == diag_cs%axesCu1%id) then
    axes => diag_cs%axesCu1
  elseif (axes_in%id == diag_cs%axesCv1%id) then
    axes => diag_cs%axesCv1
  else
    allocate(axes)
    axes = axes_in
  endif

  MOM_missing_value = axes%diag_cs%missing_value
  if (present(missing_value)) MOM_missing_value = missing_value

  diag_cs => axes%diag_cs
  dm_id = -1

  module_list = "{"//trim(module_name)
  num_modnm = 1

  ! Register the native diagnostic
  active = register_diag_field_expand_cmor(dm_id, module_name, field_name, axes, &
             init_time, long_name=long_name, units=units, missing_value=MOM_missing_value, &
             range=range, mask_variant=mask_variant, standard_name=standard_name, &
             verbose=verbose, do_not_log=do_not_log, err_msg=err_msg, &
             interp_method=interp_method, tile_count=tile_count, &
             cmor_field_name=cmor_field_name, cmor_long_name=cmor_long_name, &
             cmor_units=cmor_units, cmor_standard_name=cmor_standard_name, &
             cell_methods=cell_methods, x_cell_method=x_cell_method, y_cell_method=y_cell_method, &
             conversion=conversion)
  num_varnm = 1 ; var_list = "{"//trim(field_name)
  if (present(cmor_field_name)) then
    num_varnm = num_varnm + 1
    var_list = trim(var_list)//","//trim(cmor_field_name)
  endif
  var_list = trim(var_list)//"}"

  dimensions = ""
  if (axes_in%is_h_point)   dimensions = trim(dimensions)//" xh, yh,"
  if (axes_in%is_q_point)   dimensions = trim(dimensions)//" xq, yq,"
  if (axes_in%is_u_point)   dimensions = trim(dimensions)//" xq, yh,"
  if (axes_in%is_v_point)   dimensions = trim(dimensions)//" xh, yq,"
  if (len_trim(dimensions) > 0) dimensions = trim_trailing_commas(dimensions)

  if (is_root_pe() .and. (diag_CS%available_diag_doc_unit > 0)) then
    msg = ''
    if (present(cmor_field_name)) msg = 'CMOR equivalent is "'//trim(cmor_field_name)//'"'
    call attach_cell_methods(-1, axes, cm_string, cell_methods, x_cell_method, y_cell_method)
    module_list = trim(module_list)//"}"
    if (num_modnm <= 1) module_list = module_name
    if (num_varnm <= 1) var_list = ''

    call log_available_diag(dm_id>0, module_list, field_name, cm_string, msg, diag_CS, &
                            long_name, units, standard_name, variants=var_list, dimensions=dimensions)
  endif

  register_diag_field = dm_id

end procedure register_MOM_IS_diag_field
module procedure register_diag_field_expand_cmor
  real :: MOM_missing_value ! A value used to indicate missing values in output files, in arbitrary units [a]
  type(diag_ctrl), pointer :: diag_cs => null()
  type(diag_type), pointer :: this_diag => null()
  integer :: fms_id
  character(len=256) :: posted_cmor_units, posted_cmor_standard_name, posted_cmor_long_name
  character(len=256) :: cm_string ! A string describing the cell methods returned from attach_cell_methods.
  MOM_missing_value = axes%diag_cs%missing_value
  if (present(missing_value)) MOM_missing_value = missing_value

  register_diag_field_expand_cmor = .false.
  diag_cs => axes%diag_cs

  ! Set up the 'primary' diagnostic, first get an underlying FMS id
  fms_id = register_diag_field_expand_axes(module_name, field_name, axes, init_time, &
             long_name=long_name, units=units, missing_value=MOM_missing_value, &
             range=range, mask_variant=mask_variant, standard_name=standard_name, &
             verbose=verbose, do_not_log=do_not_log, err_msg=err_msg, &
             interp_method=interp_method, tile_count=tile_count)
  if (.not. diag_cs%diag_as_chksum) &
    call attach_cell_methods(fms_id, axes, cm_string, cell_methods, x_cell_method, y_cell_method)

  this_diag => null()
  if (fms_id /= DIAG_FIELD_NOT_FOUND) then
    call add_diag_to_list(diag_cs, dm_id, fms_id, this_diag, axes, module_name, field_name)
    if (present(conversion)) this_diag%conversion_factor = conversion
    register_diag_field_expand_cmor = .true.
  endif

  ! For the CMOR variation of the above diagnostic
  if (present(cmor_field_name) .and. .not. diag_cs%diag_as_chksum) then
    ! Fallback values for strings set to "NULL"
    posted_cmor_units = "not provided"         !
    posted_cmor_standard_name = "not provided" ! Values might be able to be replaced with a CS%missing field?
    posted_cmor_long_name = "not provided"     !

    ! If attributes are present for MOM variable names, use them first for the register_MOM_IS_diag_field
    ! call for CMOR verison of the variable
    if (present(units)) posted_cmor_units = units
    if (present(standard_name)) posted_cmor_standard_name = standard_name
    if (present(long_name)) posted_cmor_long_name = long_name

    ! If specified in the call to register_MOM_IS_diag_field, override attributes with the CMOR versions
    if (present(cmor_units)) posted_cmor_units = cmor_units
    if (present(cmor_standard_name)) posted_cmor_standard_name = cmor_standard_name
    if (present(cmor_long_name)) posted_cmor_long_name = cmor_long_name

    fms_id = register_diag_field_expand_axes(module_name, cmor_field_name, axes, init_time,    &
               long_name=trim(posted_cmor_long_name), units=trim(posted_cmor_units),                  &
               missing_value=MOM_missing_value, range=range, mask_variant=mask_variant,               &
               standard_name=trim(posted_cmor_standard_name), verbose=verbose, do_not_log=do_not_log, &
               err_msg=err_msg, interp_method=interp_method, tile_count=tile_count)
    call attach_cell_methods(fms_id, axes, cm_string, cell_methods, x_cell_method, y_cell_method)

    this_diag => null()
    if (fms_id /= DIAG_FIELD_NOT_FOUND) then
      call add_diag_to_list(diag_cs, dm_id, fms_id, this_diag, axes, module_name, field_name)
      if (present(conversion)) this_diag%conversion_factor = conversion
      register_diag_field_expand_cmor = .true.
    endif
  endif

end procedure register_diag_field_expand_cmor
module procedure register_diag_field_expand_axes
  integer :: fms_id, area_id
  area_id = axes%id_area

  ! Get the FMS diagnostic id
  if (axes%diag_cs%diag_as_chksum) then
    fms_id = axes%diag_cs%num_chksum_diags + 1
    axes%diag_cs%num_chksum_diags = fms_id
  elseif (present(interp_method) .or. axes%is_h_point) then
    ! If interp_method is provided we must use it
    if (area_id>0) then
      fms_id = register_diag_field_infra(module_name, field_name, axes%handles, &
               init_time, long_name=long_name, units=units, missing_value=missing_value, &
               range=range, mask_variant=mask_variant, standard_name=standard_name, &
               verbose=verbose, do_not_log=do_not_log, err_msg=err_msg, &
               interp_method=interp_method, tile_count=tile_count, area=area_id)
    else
      fms_id = register_diag_field_infra(module_name, field_name, axes%handles, &
               init_time, long_name=long_name, units=units, missing_value=missing_value, &
               range=range, mask_variant=mask_variant, standard_name=standard_name, &
               verbose=verbose, do_not_log=do_not_log, err_msg=err_msg, &
               interp_method=interp_method, tile_count=tile_count)
    endif
  else
    ! If interp_method is not provided and the field is not at an h-point then interp_method='none'
    if (area_id>0) then
      fms_id = register_diag_field_infra(module_name, field_name, axes%handles, &
               init_time, long_name=long_name, units=units, missing_value=missing_value, &
               range=range, mask_variant=mask_variant, standard_name=standard_name, &
               verbose=verbose, do_not_log=do_not_log, err_msg=err_msg, &
               interp_method='none', tile_count=tile_count, area=area_id)
    else
      fms_id = register_diag_field_infra(module_name, field_name, axes%handles, &
               init_time, long_name=long_name, units=units, missing_value=missing_value, &
               range=range, mask_variant=mask_variant, standard_name=standard_name, &
               verbose=verbose, do_not_log=do_not_log, err_msg=err_msg, &
               interp_method='none', tile_count=tile_count)
    endif
  endif

  register_diag_field_expand_axes = fms_id

end procedure register_diag_field_expand_axes
module procedure add_diag_to_list
  if (dm_id == -1) dm_id = get_new_diag_id(diag_cs)
  ! Create a new diag_type to store links in
  call alloc_diag_with_id(dm_id, diag_cs, this_diag)
  call assert(associated(this_diag), 'add_diag_to_list: allocation failed for '//trim(field_name))
  ! Record FMS id, masks and conversion factor, in diag_type
  this_diag%fms_diag_id = fms_id
  this_diag%debug_str = trim(module_name)//"-"//trim(field_name)
  this_diag%axes => axes

end procedure add_diag_to_list
module procedure attach_cell_methods
  character(len=9) :: axis_name
  logical :: x_mean, y_mean, x_sum, y_sum
  x_mean = .false.
  y_mean = .false.
  x_sum = .false.
  y_sum = .false.

  ostring = ''
  if (present(cell_methods)) then
    if (present(x_cell_method) .or. present(y_cell_method)) then
      call MOM_error(FATAL, "attach_cell_methods: " // &
           'Individual direction cell method was specified along with a "cell_methods" string.')
    endif
    if (len(trim(cell_methods))>0) then
      call MOM_diag_field_add_attribute(id, 'cell_methods', trim(cell_methods))
      ostring = trim(cell_methods)
    endif
  else
    if (present(x_cell_method)) then
      if (len(trim(x_cell_method))>0) then
        call get_MOM_diag_axis_name(axes%handles(1), axis_name)
        call MOM_diag_field_add_attribute(id, 'cell_methods', trim(axis_name)//':'//trim(x_cell_method))
        ostring = trim(adjustl(ostring))//' '//trim(axis_name)//':'//trim(x_cell_method)
        if (trim(x_cell_method)=='mean') x_mean=.true.
        if (trim(x_cell_method)=='sum') x_sum=.true.
      endif
    else
      if (len(trim(axes%x_cell_method))>0) then
        call get_MOM_diag_axis_name(axes%handles(1), axis_name)
        call MOM_diag_field_add_attribute(id, 'cell_methods', trim(axis_name)//':'//trim(axes%x_cell_method))
        ostring = trim(adjustl(ostring))//' '//trim(axis_name)//':'//trim(axes%x_cell_method)
        if (trim(axes%x_cell_method)=='mean') x_mean=.true.
        if (trim(axes%x_cell_method)=='sum') x_sum=.true.
      endif
    endif
    if (present(y_cell_method)) then
      if (len(trim(y_cell_method))>0) then
        call get_MOM_diag_axis_name(axes%handles(2), axis_name)
        call MOM_diag_field_add_attribute(id, 'cell_methods', trim(axis_name)//':'//trim(y_cell_method))
        ostring = trim(adjustl(ostring))//' '//trim(axis_name)//':'//trim(y_cell_method)
        if (trim(y_cell_method)=='mean') y_mean=.true.
        if (trim(y_cell_method)=='sum') y_sum=.true.
      endif
    else
      if (len(trim(axes%y_cell_method))>0) then
        call get_MOM_diag_axis_name(axes%handles(2), axis_name)
        call MOM_diag_field_add_attribute(id, 'cell_methods', trim(axis_name)//':'//trim(axes%y_cell_method))
        ostring = trim(adjustl(ostring))//' '//trim(axis_name)//':'//trim(axes%y_cell_method)
        if (trim(axes%y_cell_method)=='mean') y_mean=.true.
        if (trim(axes%y_cell_method)=='sum') y_sum=.true.
      endif
    endif
    if (x_mean .and. y_mean) then
      call MOM_diag_field_add_attribute(id, 'cell_methods', 'area:mean')
      ostring = trim(adjustl(ostring))//' area:mean'
    elseif (x_sum .and. y_sum) then
      call MOM_diag_field_add_attribute(id, 'cell_methods', 'area:sum')
      ostring = trim(adjustl(ostring))//' area:sum'
    endif
  endif
  ostring = adjustl(ostring)
end procedure attach_cell_methods
module procedure register_scalar_field_axes
  register_scalar_field = register_scalar_field_CS(module_name, field_name, init_time, axes%diag_cs, &
            long_name, units, missing_value, range, standard_name, &
            do_not_log, err_msg, interp_method, cmor_field_name, &
            cmor_long_name, cmor_units, cmor_standard_name, conversion)

end procedure register_scalar_field_axes
module procedure register_scalar_field_CS
  real :: MOM_missing_value ! A value used to indicate missing values in output files, in arbitrary units [a]
  integer :: dm_id, fms_id
  type(diag_type), pointer :: diag => null(), cmor_diag => null()
  character(len=256) :: posted_cmor_units, posted_cmor_standard_name, posted_cmor_long_name
  character(len=16)  :: dimensions
  MOM_missing_value = diag_cs%missing_value
  if (present(missing_value)) MOM_missing_value = missing_value

  dm_id = -1
  diag => null()
  cmor_diag => null()

  if (diag_cs%diag_as_chksum) then
    fms_id = diag_cs%num_chksum_diags + 1
    diag_cs%num_chksum_diags = fms_id
  else
    fms_id = register_diag_field_infra(module_name, field_name, init_time, &
                long_name=long_name, units=units, missing_value=MOM_missing_value, &
                range=range, standard_name=standard_name, do_not_log=do_not_log, &
                err_msg=err_msg)
  endif

  if (fms_id /= DIAG_FIELD_NOT_FOUND) then
    dm_id = get_new_diag_id(diag_cs)
    call alloc_diag_with_id(dm_id, diag_cs, diag)
    call assert(associated(diag), 'register_scalar_field: diag allocation failed')
    diag%fms_diag_id = fms_id
    diag%debug_str = trim(module_name)//"-"//trim(field_name)
    if (present(conversion)) diag%conversion_factor = conversion
  endif

  if (present(cmor_field_name)) then
    ! Fallback values for strings set to "not provided"
    posted_cmor_units = "not provided"
    posted_cmor_standard_name = "not provided"
    posted_cmor_long_name = "not provided"

    ! If attributes are present for MOM variable names, use them as defaults for the
    ! register_diag_field_infra call for CMOR verison of the variable
    if (present(units)) posted_cmor_units = units
    if (present(standard_name)) posted_cmor_standard_name = standard_name
    if (present(long_name)) posted_cmor_long_name = long_name

    ! If specified in the call to register_MOM_IS_scalar_field, override attributes with the CMOR versions
    if (present(cmor_units)) posted_cmor_units = cmor_units
    if (present(cmor_standard_name)) posted_cmor_standard_name = cmor_standard_name
    if (present(cmor_long_name)) posted_cmor_long_name = cmor_long_name

    fms_id = register_diag_field_infra(module_name, cmor_field_name, init_time, &
           long_name=trim(posted_cmor_long_name), units=trim(posted_cmor_units), &
           missing_value=MOM_missing_value, range=range, &
           standard_name=trim(posted_cmor_standard_name), do_not_log=do_not_log, err_msg=err_msg)
    if (fms_id /= DIAG_FIELD_NOT_FOUND) then
      if (dm_id == -1) then
        dm_id = get_new_diag_id(diag_cs)
      endif
      call alloc_diag_with_id(dm_id, diag_cs, cmor_diag)
      cmor_diag%fms_diag_id = fms_id
      cmor_diag%debug_str = trim(module_name)//"-"//trim(cmor_field_name)
      if (present(conversion)) cmor_diag%conversion_factor = conversion
    endif
  endif

  dimensions = "scalar"

  ! Document diagnostics in list of available diagnostics
  if (is_root_pe() .and. diag_CS%available_diag_doc_unit > 0) then
    if (present(cmor_field_name)) then
      call log_available_diag(associated(diag), module_name, field_name, '', '', diag_CS, &
                              long_name, units, standard_name, &
                              variants="{"//trim(field_name)//","//trim(cmor_field_name)//"}", &
                              dimensions=dimensions)
    else
      call log_available_diag(associated(diag), module_name, field_name, '', '', diag_CS, &
                              long_name, units, standard_name, dimensions=dimensions)
    endif
  endif

  register_scalar_field = dm_id

end procedure register_scalar_field_CS
module procedure register_MOM_IS_static_field
  real :: MOM_missing_value ! A value used to indicate missing values in output files, in arbitrary units [a]
  type(diag_ctrl), pointer :: diag_cs => null() !< A structure that is used to regulate diagnostic output
  type(diag_type), pointer :: diag => null(), cmor_diag => null()
  integer :: dm_id, fms_id
  character(len=256) :: posted_cmor_units, posted_cmor_standard_name, posted_cmor_long_name
  character(len=9) :: axis_name
  character(len=24) :: dimensions
  MOM_missing_value = axes%diag_cs%missing_value
  if (present(missing_value)) MOM_missing_value = missing_value

  diag_cs => axes%diag_cs
  dm_id = -1
  diag => null()
  cmor_diag => null()

  if (diag_cs%diag_as_chksum) then
    fms_id = diag_cs%num_chksum_diags + 1
    diag_cs%num_chksum_diags = fms_id
  else
    fms_id = register_static_field_infra(module_name, field_name, axes%handles, &
           long_name=long_name, units=units, missing_value=MOM_missing_value, &
           range=range, mask_variant=mask_variant, standard_name=standard_name, &
           do_not_log=do_not_log, &
           interp_method=interp_method, tile_count=tile_count, area=area)
  endif

  if (fms_id /= DIAG_FIELD_NOT_FOUND) then
    dm_id = get_new_diag_id(diag_cs)
    call alloc_diag_with_id(dm_id, diag_cs, diag)
    call assert(associated(diag), 'register_static_field: diag allocation failed')
    diag%fms_diag_id = fms_id
    diag%debug_str = trim(module_name)//"-"//trim(field_name)
    if (present(conversion)) diag%conversion_factor = conversion

    if (diag_cs%diag_as_chksum) then
      diag%axes => axes
    else
      if (present(x_cell_method)) then
        call get_MOM_diag_axis_name(axes%handles(1), axis_name)
        call MOM_diag_field_add_attribute(fms_id, 'cell_methods', &
            trim(axis_name)//':'//trim(x_cell_method))
      endif
      if (present(y_cell_method)) then
        call get_MOM_diag_axis_name(axes%handles(2), axis_name)
        call MOM_diag_field_add_attribute(fms_id, 'cell_methods', &
            trim(axis_name)//':'//trim(y_cell_method))
      endif
      if (present(area_cell_method)) then
        call MOM_diag_field_add_attribute(fms_id, 'cell_methods', &
            'area:'//trim(area_cell_method))
      endif
    endif
  endif

  if (present(cmor_field_name) .and. .not. diag_cs%diag_as_chksum) then
    ! Fallback values for strings set to "not provided"
    posted_cmor_units = "not provided"
    posted_cmor_standard_name = "not provided"
    posted_cmor_long_name = "not provided"

    ! If attributes are present for MOM variable names, use them first for the register_static_field
    ! call for CMOR verison of the variable
    if (present(units)) posted_cmor_units = units
    if (present(standard_name)) posted_cmor_standard_name = standard_name
    if (present(long_name)) posted_cmor_long_name = long_name

    ! If specified in the call to register_static_field, override attributes with the CMOR versions
    if (present(cmor_units)) posted_cmor_units = cmor_units
    if (present(cmor_standard_name)) posted_cmor_standard_name = cmor_standard_name
    if (present(cmor_long_name)) posted_cmor_long_name = cmor_long_name

    fms_id = register_static_field_infra(module_name, cmor_field_name, axes%handles, &
                long_name=trim(posted_cmor_long_name), units=trim(posted_cmor_units), &
                missing_value=MOM_missing_value, range=range, mask_variant=mask_variant, &
                standard_name=trim(posted_cmor_standard_name), do_not_log=do_not_log, &
                interp_method=interp_method, tile_count=tile_count, area=area)
    if (fms_id /= DIAG_FIELD_NOT_FOUND) then
      if (dm_id == -1) then
        dm_id = get_new_diag_id(diag_cs)
      endif
      call alloc_diag_with_id(dm_id, diag_cs, cmor_diag)
      cmor_diag%fms_diag_id = fms_id
      cmor_diag%debug_str = trim(module_name)//"-"//trim(cmor_field_name)
      if (present(conversion)) cmor_diag%conversion_factor = conversion
      if (present(x_cell_method)) then
        call get_MOM_diag_axis_name(axes%handles(1), axis_name)
        call MOM_diag_field_add_attribute(fms_id, 'cell_methods', trim(axis_name)//':'//trim(x_cell_method))
      endif
      if (present(y_cell_method)) then
        call get_MOM_diag_axis_name(axes%handles(2), axis_name)
        call MOM_diag_field_add_attribute(fms_id, 'cell_methods', trim(axis_name)//':'//trim(y_cell_method))
      endif
      if (present(area_cell_method)) then
        call MOM_diag_field_add_attribute(fms_id, 'cell_methods', 'area:'//trim(area_cell_method))
      endif
    endif
  endif

  dimensions = ""
  if (axes%is_h_point)   dimensions = trim(dimensions)//" xh, yh,"
  if (axes%is_q_point)   dimensions = trim(dimensions)//" xq, yq,"
  if (axes%is_u_point)   dimensions = trim(dimensions)//" xq, yh,"
  if (axes%is_v_point)   dimensions = trim(dimensions)//" xh, yq,"
  if (len_trim(dimensions) > 0) dimensions = trim_trailing_commas(dimensions)

  ! Document diagnostics in list of available diagnostics
  if (is_root_pe() .and. diag_CS%available_diag_doc_unit > 0) then
    if (present(cmor_field_name)) then
      call log_available_diag(associated(diag), module_name, field_name, '', '', diag_CS, &
                              long_name, units, standard_name, &
                              variants="{"//trim(field_name)//","//trim(cmor_field_name)//"}", &
                              dimensions=dimensions)
    else
      call log_available_diag(associated(diag), module_name, field_name, '', '', diag_CS, &
                              long_name, units, standard_name, dimensions=dimensions)
    endif
  endif

  register_static_field = dm_id

end procedure register_MOM_IS_static_field
module procedure describe_option
  character(len=480) :: mesg
  integer :: len_ind
  len_ind = len_trim(value)

  mesg = "    ! "//trim(opt_name)//": "//trim(value)
  write(diag_CS%available_diag_doc_unit, '(a)') trim(mesg)
end procedure describe_option
module procedure MOM_IS_diag_mediator_init
  integer :: ios, i, new_unit
  logical :: opened, new_file
  character(len=8)   :: this_pe
  character(len=240) :: doc_file, doc_file_dflt, doc_path
  character(len=40)  :: doc_file_param
# include "version_variable.h"
  character(len=40) :: mdl = "MOM_IS_diag_mediator" ! This module's name.
  character(len=32) :: filename_appendix = '' !fms appendix to filename for ensemble runs
  call MOM_diag_manager_init(err_msg=err_msg)

  id_clock_diag_mediator = cpu_clock_id('(Ice shelf diagnostics framework)', grain=CLOCK_MODULE)

  ! Allocate and initialize list of all diagnostics (and variants)
  allocate(diag_cs%diags(DIAG_ALLOC_CHUNK_SIZE))
  diag_cs%next_free_diag_id = 1
  do i=1, DIAG_ALLOC_CHUNK_SIZE
    call initialize_diag_type(diag_cs%diags(i))
  enddo

  diag_cs%show_call_tree = callTree_showQuery()

  ! Read all relevant parameters and write them to the model log.
  call log_version(param_file, mdl, version, "")

  call get_param(param_file, mdl, 'USE_INDEX_DIAGNOSTIC_AXES', diag_cs%index_space_axes, &
                 'If true, use a grid index coordinate convention for diagnostic axes. ',&
                 default=.false.)

  call get_param(param_file, mdl, 'DIAG_MISVAL', diag_cs%missing_value, &
                 'Set the default missing value to use for diagnostics.', &
                 units="various", default=-1.e34)
  call get_param(param_file, mdl, 'DIAG_AS_CHKSUM', diag_cs%diag_as_chksum, &
                 'Instead of writing diagnostics to the diag manager, write '//&
                 'a text file containing the checksum (bitcount) of the array.',  &
                 default=.false.)

  if (diag_cs%diag_as_chksum) &
    diag_cs%num_chksum_diags = 0

  ! Keep pointers to the grid for diagnostic checksums
  diag_cs%G => G
  diag_cs%US => US

  diag_cs%is = G%isc - (G%isd-1) ; diag_cs%ie = G%iec - (G%isd-1)
  diag_cs%js = G%jsc - (G%jsd-1) ; diag_cs%je = G%jec - (G%jsd-1)
  diag_cs%isd = G%isd ; diag_cs%ied = G%ied
  diag_cs%jsd = G%jsd ; diag_cs%jed = G%jed

  ! Initialize available diagnostic log file
  if (is_root_pe() .and. (diag_CS%available_diag_doc_unit < 0)) then
    if (present(component)) then
      doc_file_dflt = trim(component)//".available_diags"
      doc_file_param = trim(uppercase(component))//"_AVAILABLE_DIAGS_FILE"
    else
      write(this_pe,'(i6.6)') PE_here()
      doc_file_dflt = "MOM_IS.available_diags."//this_pe
      doc_file_param = "AVAILABLE_MOM_IS_DIAGS_FILE"
    endif
    call get_param(param_file, mdl, trim(doc_file_param), doc_file, &
                 "A file into which to write a list of all available "//&
                 "ice shelf diagnostics that can be included in a diag_table.", &
                 default=doc_file_dflt, do_not_log=(diag_CS%available_diag_doc_unit/=-1))
    if (len_trim(doc_file) > 0) then
      new_file = .true. ; if (diag_CS%available_diag_doc_unit /= -1) new_file = .false.
    ! Find an unused unit number.
      do new_unit=512,42,-1
        inquire( new_unit, opened=opened)
        if (.not.opened) exit
      enddo
      if (opened) call MOM_error(FATAL, &
          "diag_mediator_init failed to find an unused unit number.")

      doc_path = doc_file
      if (present(doc_file_dir)) then ; if (len_trim(doc_file_dir) > 0) then
        doc_path = trim(slasher(doc_file_dir))//trim(doc_file)
      endif ; endif

      diag_CS%available_diag_doc_unit = new_unit

      if (new_file) then
        open(diag_CS%available_diag_doc_unit, file=trim(doc_path), access='SEQUENTIAL', form='FORMATTED', &
             action='WRITE', status='REPLACE', iostat=ios)
      else ! This file is being reopened, and should be appended.
        open(diag_CS%available_diag_doc_unit, file=trim(doc_path), access='SEQUENTIAL', form='FORMATTED', &
             action='WRITE', status='OLD', position='APPEND', iostat=ios)
      endif
      inquire(diag_CS%available_diag_doc_unit, opened=opened)
      if ((.not.opened) .or. (ios /= 0)) then
        call MOM_error(FATAL, "Failed to open available diags file "//trim(doc_path)//".")
      endif
    endif
  endif

  if (is_root_pe() .and. (diag_CS%chksum_iounit < 0) .and. diag_CS%diag_as_chksum) then
    !write(this_pe,'(i6.6)') PE_here()
    !doc_file_dflt = "chksum_diag."//this_pe
    doc_file_dflt = "chksum_diag"
    call get_param(param_file, mdl, "CHKSUM_DIAG_FILE", doc_file, &
                 "A file into which to write all checksums of the "//&
                 "diagnostics listed in the diag_table.", &
                 default=doc_file_dflt, do_not_log=(diag_CS%chksum_iounit/=-1))

    call get_filename_appendix(filename_appendix)
    if (len_trim(filename_appendix) > 0) then
      doc_file = trim(doc_file) //'.'//trim(filename_appendix)
    endif
#ifdef STATSLABEL
    doc_file = trim(doc_file)//"."//trim(adjustl(STATSLABEL))
#endif

    if (len_trim(doc_file) > 0) then
      new_file = .true. ; if (diag_CS%chksum_iounit /= -1) new_file = .false.
    ! Find an unused unit number.
      do new_unit=512,42,-1
        inquire( new_unit, opened=opened)
        if (.not.opened) exit
      enddo
      if (opened) call MOM_error(FATAL, &
          "diag_mediator_init failed to find an unused unit number.")

      doc_path = doc_file
      if (present(doc_file_dir)) then ; if (len_trim(doc_file_dir) > 0) then
        doc_path = trim(slasher(doc_file_dir))//trim(doc_file)
      endif ; endif

      diag_CS%chksum_iounit = new_unit

      if (new_file) then
        open(diag_CS%chksum_iounit, file=trim(doc_path), access='SEQUENTIAL', form='FORMATTED', &
             action='WRITE', status='REPLACE', iostat=ios)
      else ! This file is being reopened, and should be appended.
        open(diag_CS%chksum_iounit, file=trim(doc_path), access='SEQUENTIAL', form='FORMATTED', &
             action='WRITE', status='OLD', position='APPEND', iostat=ios)
      endif
      inquire(diag_CS%chksum_iounit, opened=opened)
      if ((.not.opened) .or. (ios /= 0)) then
        call MOM_error(FATAL, "Failed to open checksum diags file "//trim(doc_path)//".")
      endif
    endif
  endif

  call diag_masks_set(G, diag_cs%missing_value, diag_cs)

end procedure MOM_IS_diag_mediator_init
module procedure diag_masks_set
  integer :: i, j
  diag_cs%mask2dT  => G%mask2dT
  diag_cs%mask2dBu => G%mask2dBu
  diag_cs%mask2dCu => G%mask2dCu
  diag_cs%mask2dCv => G%mask2dCv

  allocate(diag_cs%mask2dT_comp(G%isc:G%iec,G%jsc:G%jec))
  do j=G%jsc,G%jec ; do i=G%isc,G%iec
    diag_cs%mask2dT_comp(i,j) = diag_cs%mask2dT(i,j)
  enddo ; enddo

  diag_cs%missing_value = missing_value

end procedure diag_masks_set
module procedure MOM_IS_diag_mediator_close_registration
  if (diag_CS%available_diag_doc_unit > -1) then
    close(diag_CS%available_diag_doc_unit) ; diag_CS%available_diag_doc_unit = -2
  endif

end procedure MOM_IS_diag_mediator_close_registration
module procedure MOM_IS_diag_mediator_end
  type(diag_type), pointer :: diag, next_diag
  integer :: i
  if (diag_CS%available_diag_doc_unit > -1) then
    close(diag_CS%available_diag_doc_unit) ; diag_CS%available_diag_doc_unit = -3
  endif
  if (diag_CS%chksum_iounit > -1) then
    close(diag_CS%chksum_iounit) ; diag_CS%chksum_iounit = -3
  endif

  do i=1, diag_cs%next_free_diag_id - 1
    if (associated(diag_cs%diags(i)%next)) then
      next_diag => diag_cs%diags(i)%next
      do while (associated(next_diag))
        diag => next_diag
        next_diag => diag%next
        deallocate(diag)
      enddo
    endif
  enddo

  deallocate(diag_cs%diags)

  ! These points to arrays in the grid type, so they can not be deallocated here.
  if (associated(diag_cs%mask2dT))  diag_cs%mask2dT => NULL()
  if (associated(diag_cs%mask2dBu)) diag_cs%mask2dBu => NULL()
  if (associated(diag_cs%mask2dCu)) diag_cs%mask2dCu => NULL()
  if (associated(diag_cs%mask2dCv)) diag_cs%mask2dCv => NULL()
  if (associated(diag_cs%mask2dT_comp)) deallocate(diag_cs%mask2dT_comp)

end procedure MOM_IS_diag_mediator_end
module procedure get_new_diag_id
  type(diag_type), dimension(:), allocatable :: tmp
  integer :: i
  if (diag_cs%next_free_diag_id > size(diag_cs%diags)) then
    call assert(diag_cs%next_free_diag_id - size(diag_cs%diags) == 1, &
                'get_new_diag_id: inconsistent diag id')

    ! Increase the size of diag_cs%diags and copy data over.
    ! Do not use move_alloc() because it is not supported by Fortran 90
    allocate(tmp(size(diag_cs%diags)))
    tmp(:) = diag_cs%diags(:)
    deallocate(diag_cs%diags)
    allocate(diag_cs%diags(size(tmp) + DIAG_ALLOC_CHUNK_SIZE))
    diag_cs%diags(1:size(tmp)) = tmp(:)
    deallocate(tmp)

    ! Initialize new part of the diag array.
    do i=diag_cs%next_free_diag_id, size(diag_cs%diags)
      call initialize_diag_type(diag_cs%diags(i))
    enddo
  endif

  get_new_diag_id = diag_cs%next_free_diag_id
  diag_cs%next_free_diag_id = diag_cs%next_free_diag_id + 1

end procedure get_new_diag_id
module procedure initialize_diag_type
  diag%in_use = .false.
  diag%fms_diag_id = -1
  diag%axes => null()
  diag%next => null()
  diag%conversion_factor = 0.

end procedure initialize_diag_type
module procedure alloc_diag_with_id
  type(diag_type), pointer :: tmp => NULL()
  if (.not. diag_cs%diags(diag_id)%in_use) then
    diag => diag_cs%diags(diag_id)
  else
    allocate(diag)
    tmp => diag_cs%diags(diag_id)%next
    diag_cs%diags(diag_id)%next => diag
    diag%next => tmp
  endif
  diag%in_use = .true.

end procedure alloc_diag_with_id
module procedure log_available_diag
  character(len=240) :: mesg
  if (used) then
    mesg = '"'//trim(field_name)//'"  [Used]'
  else
    mesg = '"'//trim(field_name)//'"  [Unused]'
  endif
  if (len(trim((comment)))>0) then
    write(diag_CS%available_diag_doc_unit, '(a,1x,"(",a,")")') trim(mesg),trim(comment)
  else
    write(diag_CS%available_diag_doc_unit, '(a)') trim(mesg)
  endif
  call describe_option("modules", module_name, diag_CS)
  if (present(dimensions)) then ; if (len(trim(dimensions)) > 0) then
    call describe_option("dimensions", dimensions, diag_CS)
  endif ; endif
  if (present(long_name)) call describe_option("long_name", long_name, diag_CS)
  if (present(units)) call describe_option("units", units, diag_CS)
  if (present(standard_name)) &
    call describe_option("standard_name", standard_name, diag_CS)
  if (len(trim((cell_methods_string)))>0) &
    call describe_option("cell_methods", trim(cell_methods_string), diag_CS)
  if (present(variants)) then ; if (len(trim(variants)) > 0) then
    call describe_option("variants", variants, diag_CS)
  endif ; endif
end procedure log_available_diag
module procedure log_chksum_diag
  write(docunit, '(a,1x,i9.8)') description, chksum
  flush(docunit)

end procedure log_chksum_diag
module procedure found_in_diagtable
  integer :: handle ! Integer handle returned from diag_manager
  handle = register_static_field_infra('ice_shelf_model', varName, diag%axesT1%handles)

  found_in_diagtable = (handle>0)

end procedure found_in_diagtable
module procedure MOM_IS_diag_send_complete
  call diag_send_complete_infra()
end procedure MOM_IS_diag_send_complete
end submodule MOM_IS_diag_mediator_s
