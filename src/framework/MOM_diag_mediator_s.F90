submodule (MOM_diag_mediator) MOM_diag_mediator_s
#define IMPLIES(A, B) ((.not. (A)) .or. (B))
#define DIAG_ALLOC_CHUNK_SIZE 100
  implicit none
contains
module procedure set_axes_info
  integer :: id_xq, id_yq, id_zl, id_zi, id_xh, id_yh, id_null
  integer :: id_zl_native, id_zi_native
  integer :: i, j, nz
  real :: zlev(GV%ke)     ! Numerical values for layer vertical coordinates, in unscaled units
  real :: zinter(GV%ke+1) ! Numerical values for interface vertical coordinates, in unscaled units
  logical :: set_vert
  real, allocatable, dimension(:) :: IaxB, iax ! Index-based integer and half-integer i-axis labels [nondim]
  real, allocatable, dimension(:) :: JaxB, jax ! Index-based integer and half-integer j-axis labels [nondim]
  set_vert = .true. ; if (present(set_vertical)) set_vert = set_vertical

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

  ! Horizontal axes for the native grids
  if (diag_cs%index_space_axes) then
    if (G%symmetric) then
      id_xq = diag_axis_init('Iq', IaxB(G%IsgB:G%IegB), 'none', 'x', &
          'Boundary (q) point grid-space longitude', G%Domain, position=EAST)
      id_yq = diag_axis_init('Jq', JaxB(G%JsgB:G%JegB), 'none', 'y', &
          'Boundary (q) point grid-space latitude', G%Domain, position=NORTH)
    else
      id_xq = diag_axis_init('Iq', IaxB(G%isg:G%ieg), 'none', 'x', &
          'Boundary (q) point grid-space longitude', G%Domain, position=EAST)
      id_yq = diag_axis_init('Jq', JaxB(G%jsg:G%jeg), 'none', 'y', &
          'Boundary (q) point grid-space latitude', G%Domain, position=NORTH)
    endif
    id_xh = diag_axis_init('ih', iax(G%isg:G%ieg), 'none', 'x', &
        'Tracer (h) point grid-space longitude', G%Domain)
    id_yh = diag_axis_init('jh', jax(G%jsg:G%jeg), 'none', 'y', &
        'Tracer (h) point grid-space latitude', G%Domain)
  else
    if (G%symmetric) then
      id_xq = diag_axis_init('xq', G%gridLonB(G%IsgB:G%IegB), G%x_axis_units, 'x', &
          'q point nominal longitude', G%Domain, position=EAST)
      id_yq = diag_axis_init('yq', G%gridLatB(G%JsgB:G%JegB), G%y_axis_units, 'y', &
          'q point nominal latitude', G%Domain, position=NORTH)
    else
      id_xq = diag_axis_init('xq', G%gridLonB(G%isg:G%ieg), G%x_axis_units, 'x', &
          'q point nominal longitude', G%Domain, position=EAST)
      id_yq = diag_axis_init('yq', G%gridLatB(G%jsg:G%jeg), G%y_axis_units, 'y', &
          'q point nominal latitude', G%Domain, position=NORTH)
    endif
    id_xh = diag_axis_init('xh', G%gridLonT(G%isg:G%ieg), G%x_axis_units, 'x', &
        'h point nominal longitude', G%Domain)
    id_yh = diag_axis_init('yh', G%gridLatT(G%jsg:G%jeg), G%y_axis_units, 'y', &
        'h point nominal latitude', G%Domain)
  endif

  if (set_vert) then
    nz = GV%ke
    zinter(1:nz+1) = GV%sInterface(1:nz+1)
    zlev(1:nz) = GV%sLayer(1:nz)
    id_zl = diag_axis_init('zl', zlev, trim(GV%zAxisUnits), 'z', &
                           'Layer '//trim(GV%zAxisLongName), direction=GV%direction)
    id_zi = diag_axis_init('zi', zinter, trim(GV%zAxisUnits), 'z', &
                           'Interface '//trim(GV%zAxisLongName), direction=GV%direction)
  else
    id_zl = -1 ; id_zi = -1
  endif
  id_zl_native = id_zl ; id_zi_native = id_zi
  ! Vertical axes for the interfaces and layers
  call define_axes_group(diag_cs, (/ id_zi /), diag_cs%axesZi, &
       v_cell_method='point', is_interface=.true.)
  call define_axes_group(diag_cs, (/ id_zL /), diag_cs%axesZL, &
       v_cell_method='mean', is_layer=.true.)

  ! Axis groupings for the model layers
  call define_axes_group(diag_cs, (/ id_xh, id_yh, id_zL /), diag_cs%axesTL, &
       x_cell_method='mean', y_cell_method='mean', v_cell_method='mean', &
       is_h_point=.true., is_layer=.true., xyave_axes=diag_cs%axesZL)
  call define_axes_group(diag_cs, (/ id_xq, id_yq, id_zL /), diag_cs%axesBL, &
       x_cell_method='point', y_cell_method='point', v_cell_method='mean', &
       is_q_point=.true., is_layer=.true.)
  call define_axes_group(diag_cs, (/ id_xq, id_yh, id_zL /), diag_cs%axesCuL, &
       x_cell_method='point', y_cell_method='mean', v_cell_method='mean', &
       is_u_point=.true., is_layer=.true., xyave_axes=diag_cs%axesZL)
  call define_axes_group(diag_cs, (/ id_xh, id_yq, id_zL /), diag_cs%axesCvL, &
       x_cell_method='mean', y_cell_method='point', v_cell_method='mean', &
       is_v_point=.true., is_layer=.true., xyave_axes=diag_cs%axesZL)

  ! Axis groupings for the model interfaces
  call define_axes_group(diag_cs, (/ id_xh, id_yh, id_zi /), diag_cs%axesTi, &
       x_cell_method='mean', y_cell_method='mean', v_cell_method='point', &
       is_h_point=.true., is_interface=.true., xyave_axes=diag_cs%axesZi)
  call define_axes_group(diag_cs, (/ id_xq, id_yq, id_zi /), diag_cs%axesBi, &
       x_cell_method='point', y_cell_method='point', v_cell_method='point', &
       is_q_point=.true., is_interface=.true.)
  call define_axes_group(diag_cs, (/ id_xq, id_yh, id_zi /), diag_cs%axesCui, &
       x_cell_method='point', y_cell_method='mean', v_cell_method='point', &
       is_u_point=.true., is_interface=.true., xyave_axes=diag_cs%axesZi)
  call define_axes_group(diag_cs, (/ id_xh, id_yq, id_zi /), diag_cs%axesCvi, &
       x_cell_method='mean', y_cell_method='point', v_cell_method='point', &
       is_v_point=.true., is_interface=.true., xyave_axes=diag_cs%axesZi)

  ! Axis groupings for 2-D arrays
  call define_axes_group(diag_cs, (/ id_xh, id_yh /), diag_cs%axesT1, &
       x_cell_method='mean', y_cell_method='mean', is_h_point=.true.)
  call define_axes_group(diag_cs, (/ id_xq, id_yq /), diag_cs%axesB1, &
       x_cell_method='point', y_cell_method='point', is_q_point=.true.)
  call define_axes_group(diag_cs, (/ id_xq, id_yh /), diag_cs%axesCu1, &
       x_cell_method='point', y_cell_method='mean', is_u_point=.true.)
  call define_axes_group(diag_cs, (/ id_xh, id_yq /), diag_cs%axesCv1, &
       x_cell_method='mean', y_cell_method='point', is_v_point=.true.)

  ! Define array extents for all piecemeal buffers
  call set_piecemeal_extents(diag_cs)

  ! Axis group for special null axis for scalars from diag manager.
  id_null = diag_axis_init('scalar_axis', (/0./), 'none', 'N', 'none', null_axis=.true.)
  call define_axes_group(diag_cs, (/ id_null /), diag_cs%axesNull)

  ! Set axis groups for non-native, non-downsampled grids
  if (diag_cs%num_diag_coords>0) then
    allocate(diag_cs%remap_axesZL(diag_cs%num_diag_coords))
    allocate(diag_cs%remap_axesTL(diag_cs%num_diag_coords))
    allocate(diag_cs%remap_axesBL(diag_cs%num_diag_coords))
    allocate(diag_cs%remap_axesCuL(diag_cs%num_diag_coords))
    allocate(diag_cs%remap_axesCvL(diag_cs%num_diag_coords))
    allocate(diag_cs%remap_axesZi(diag_cs%num_diag_coords))
    allocate(diag_cs%remap_axesTi(diag_cs%num_diag_coords))
    allocate(diag_cs%remap_axesBi(diag_cs%num_diag_coords))
    allocate(diag_cs%remap_axesCui(diag_cs%num_diag_coords))
    allocate(diag_cs%remap_axesCvi(diag_cs%num_diag_coords))
  endif

  do i=1, diag_cs%num_diag_coords
    ! For each possible diagnostic coordinate
    call diag_remap_configure_axes(diag_cs%diag_remap_cs(i), G, GV, US, param_file)

    ! Allocate these arrays since the size of the diagnostic array is now known
    allocate(diag_cs%diag_remap_cs(i)%h(G%isd:G%ied,G%jsd:G%jed, diag_cs%diag_remap_cs(i)%nz))
    allocate(diag_cs%diag_remap_cs(i)%h_extensive(G%isd:G%ied,G%jsd:G%jed, diag_cs%diag_remap_cs(i)%nz))

    ! This vertical coordinate has been configured so can be used.
    if (diag_remap_axes_configured(diag_cs%diag_remap_cs(i))) then

      ! This fetches the 1D-axis id for layers and interfaces and overwrite
      ! id_zl and id_zi from above. It also returns the number of layers.
      call diag_remap_get_axes_info(diag_cs%diag_remap_cs(i), nz, id_zL, id_zi)

      ! Axes for z layers
      call define_axes_group(diag_cs, (/ id_zL /), diag_cs%remap_axesZL(i), &
           nz=nz, vertical_coordinate_number=i, &
           v_cell_method='mean', &
           is_h_point=.true., is_layer=.true., is_native=.false., needs_remapping=.true.)
      call define_axes_group(diag_cs, (/ id_xh, id_yh, id_zL /), diag_cs%remap_axesTL(i), &
           nz=nz, vertical_coordinate_number=i, &
           x_cell_method='mean', y_cell_method='mean', v_cell_method='mean', &
           is_h_point=.true., is_layer=.true., is_native=.false., needs_remapping=.true., &
           xyave_axes=diag_cs%remap_axesZL(i))

       !! \note Remapping for B points is not yet implemented so needs_remapping is not
       !! provided for remap_axesBL
      call define_axes_group(diag_cs, (/ id_xq, id_yq, id_zL /), diag_cs%remap_axesBL(i), &
           nz=nz, vertical_coordinate_number=i, &
           x_cell_method='point', y_cell_method='point', v_cell_method='mean', &
           is_q_point=.true., is_layer=.true., is_native=.false.)

      call define_axes_group(diag_cs, (/ id_xq, id_yh, id_zL /), diag_cs%remap_axesCuL(i), &
           nz=nz, vertical_coordinate_number=i, &
           x_cell_method='point', y_cell_method='mean', v_cell_method='mean', &
           is_u_point=.true., is_layer=.true., is_native=.false., needs_remapping=.true., &
           xyave_axes=diag_cs%remap_axesZL(i))

      call define_axes_group(diag_cs, (/ id_xh, id_yq, id_zL /), diag_cs%remap_axesCvL(i), &
           nz=nz, vertical_coordinate_number=i, &
           x_cell_method='mean', y_cell_method='point', v_cell_method='mean', &
           is_v_point=.true., is_layer=.true., is_native=.false., needs_remapping=.true., &
           xyave_axes=diag_cs%remap_axesZL(i))

      ! Axes for z interfaces
      call define_axes_group(diag_cs, (/ id_zi /), diag_cs%remap_axesZi(i), &
           nz=nz, vertical_coordinate_number=i, &
           v_cell_method='point', &
           is_h_point=.true., is_interface=.true., is_native=.false., needs_interpolating=.true.)
      call define_axes_group(diag_cs, (/ id_xh, id_yh, id_zi /), diag_cs%remap_axesTi(i), &
           nz=nz, vertical_coordinate_number=i, &
           x_cell_method='mean', y_cell_method='mean', v_cell_method='point', &
           is_h_point=.true., is_interface=.true., is_native=.false., needs_interpolating=.true., &
           xyave_axes=diag_cs%remap_axesZi(i))

      !! \note Remapping for B points is not yet implemented so needs_remapping is not provided for remap_axesBi
      call define_axes_group(diag_cs, (/ id_xq, id_yq, id_zi /), diag_cs%remap_axesBi(i), &
           nz=nz, vertical_coordinate_number=i, &
           x_cell_method='point', y_cell_method='point', v_cell_method='point', &
           is_q_point=.true., is_interface=.true., is_native=.false.)

      call define_axes_group(diag_cs, (/ id_xq, id_yh, id_zi /), diag_cs%remap_axesCui(i), &
           nz=nz, vertical_coordinate_number=i, &
           x_cell_method='point', y_cell_method='mean', v_cell_method='point', &
           is_u_point=.true., is_interface=.true., is_native=.false., &
           needs_interpolating=.true., xyave_axes=diag_cs%remap_axesZi(i))

      call define_axes_group(diag_cs, (/ id_xh, id_yq, id_zi /), diag_cs%remap_axesCvi(i), &
           nz=nz, vertical_coordinate_number=i, &
           x_cell_method='mean', y_cell_method='point', v_cell_method='point', &
           is_v_point=.true., is_interface=.true., is_native=.false., &
           needs_interpolating=.true., xyave_axes=diag_cs%remap_axesZi(i))
    endif
  enddo

  if (diag_cs%index_space_axes) then
    deallocate(IaxB, iax, JaxB, jax)
  endif
  ! Define the downsampled axes
  call set_axes_info_dsamp(G, GV, param_file, diag_cs, id_zl_native, id_zi_native)

  call diag_grid_storage_init(diag_CS%diag_grid_temp, G, GV, diag_CS)

end procedure set_axes_info
module procedure set_axes_info_dsamp
  integer :: id_xq, id_yq, id_zl, id_zi, id_xh, id_yh
  integer :: i, j, nz, dl, dlfac
  real, dimension(:), pointer :: gridLonT_dsamp =>NULL() ! The longitude of downsampled T points for labeling
  real, dimension(:), pointer :: gridLatT_dsamp =>NULL() ! The latitude of downsampled T points for labeling
  real, dimension(:), pointer :: gridLonB_dsamp =>NULL() ! The longitude of downsampled B points for labeling
  real, dimension(:), pointer :: gridLatB_dsamp =>NULL() ! The latitude of downsampled B points for labeling
  do dl=1, diag_cs%num_diag_dsamp_levels
    dlfac = diag_cs%diag_dsamp_levels(dl) ! The actual downsampling factor for this level
    if (G%symmetric) then
      allocate(gridLonB_dsamp(diag_cs%dsamp(dl)%isgB:diag_cs%dsamp(dl)%iegB))
      allocate(gridLatB_dsamp(diag_cs%dsamp(dl)%jsgB:diag_cs%dsamp(dl)%jegB))
      do i=diag_cs%dsamp(dl)%isgB,diag_cs%dsamp(dl)%iegB ; gridLonB_dsamp(i) = G%gridLonB(G%isgB+dlfac*i) ; enddo
      do j=diag_cs%dsamp(dl)%jsgB,diag_cs%dsamp(dl)%jegB ; gridLatB_dsamp(j) = G%gridLatB(G%jsgB+dlfac*j) ; enddo
      id_xq = diag_axis_init('xq', gridLonB_dsamp, G%x_axis_units, 'x', &
            'q point nominal longitude', G%Domain, coarsen=dl)
      id_yq = diag_axis_init('yq', gridLatB_dsamp, G%y_axis_units, 'y', &
            'q point nominal latitude', G%Domain, coarsen=dl)
      deallocate(gridLonB_dsamp, gridLatB_dsamp)
    else
      allocate(gridLonB_dsamp(diag_cs%dsamp(dl)%isg:diag_cs%dsamp(dl)%ieg))
      allocate(gridLatB_dsamp(diag_cs%dsamp(dl)%jsg:diag_cs%dsamp(dl)%jeg))
      do i=diag_cs%dsamp(dl)%isg,diag_cs%dsamp(dl)%ieg ; gridLonB_dsamp(i) = G%gridLonB(G%isg+dlfac*i-2) ; enddo
      do j=diag_cs%dsamp(dl)%jsg,diag_cs%dsamp(dl)%jeg ; gridLatB_dsamp(j) = G%gridLatB(G%jsg+dlfac*j-2) ; enddo
      id_xq = diag_axis_init('xq', gridLonB_dsamp, G%x_axis_units, 'x', &
            'q point nominal longitude', G%Domain, coarsen=dl)
      id_yq = diag_axis_init('yq', gridLatB_dsamp, G%y_axis_units, 'y', &
            'q point nominal latitude', G%Domain, coarsen=dl)
      deallocate(gridLonB_dsamp, gridLatB_dsamp)
    endif

    allocate(gridLonT_dsamp(diag_cs%dsamp(dl)%isg:diag_cs%dsamp(dl)%ieg))
    allocate(gridLatT_dsamp(diag_cs%dsamp(dl)%jsg:diag_cs%dsamp(dl)%jeg))
    do i=diag_cs%dsamp(dl)%isg,diag_cs%dsamp(dl)%ieg ; gridLonT_dsamp(i) = G%gridLonT(G%isg+dlfac*i-2) ; enddo
    do j=diag_cs%dsamp(dl)%jsg,diag_cs%dsamp(dl)%jeg ; gridLatT_dsamp(j) = G%gridLatT(G%jsg+dlfac*j-2) ; enddo
    id_xh = diag_axis_init('xh', gridLonT_dsamp, G%x_axis_units, 'x', &
          'h point nominal longitude', G%Domain, coarsen=dl)
    id_yh = diag_axis_init('yh', gridLatT_dsamp, G%y_axis_units, 'y', &
          'h point nominal latitude', G%Domain, coarsen=dl)

    deallocate(gridLonT_dsamp, gridLatT_dsamp)

    ! Axis groupings for the model layers
    id_zl = id_zl_native ; id_zi = id_zi_native

    call define_axes_group_dsamp(diag_cs, (/ id_xh, id_yh, id_zL /), diag_cs%dsamp(dl)%axesTL, dl, &
            x_cell_method='mean', y_cell_method='mean', v_cell_method='mean', &
            is_h_point=.true., is_layer=.true., xyave_axes=diag_cs%axesZL)
    call define_axes_group_dsamp(diag_cs, (/ id_xq, id_yq, id_zL /), diag_cs%dsamp(dl)%axesBL, dl, &
            x_cell_method='point', y_cell_method='point', v_cell_method='mean', &
            is_q_point=.true., is_layer=.true.)
    call define_axes_group_dsamp(diag_cs, (/ id_xq, id_yh, id_zL /), diag_cs%dsamp(dl)%axesCuL, dl, &
            x_cell_method='point', y_cell_method='mean', v_cell_method='mean', &
            is_u_point=.true., is_layer=.true., xyave_axes=diag_cs%axesZL)
    call define_axes_group_dsamp(diag_cs, (/ id_xh, id_yq, id_zL /), diag_cs%dsamp(dl)%axesCvL, dl, &
            x_cell_method='mean', y_cell_method='point', v_cell_method='mean', &
            is_v_point=.true., is_layer=.true., xyave_axes=diag_cs%axesZL)

    ! Axis groupings for the model interfaces
    call define_axes_group_dsamp(diag_cs, (/ id_xh, id_yh, id_zi /), diag_cs%dsamp(dl)%axesTi, dl, &
            x_cell_method='mean', y_cell_method='mean', v_cell_method='point', &
            is_h_point=.true., is_interface=.true., xyave_axes=diag_cs%axesZi)
    call define_axes_group_dsamp(diag_cs, (/ id_xq, id_yq, id_zi /), diag_cs%dsamp(dl)%axesBi, dl, &
            x_cell_method='point', y_cell_method='point', v_cell_method='point', &
            is_q_point=.true., is_interface=.true.)
    call define_axes_group_dsamp(diag_cs, (/ id_xq, id_yh, id_zi /), diag_cs%dsamp(dl)%axesCui, dl, &
            x_cell_method='point', y_cell_method='mean', v_cell_method='point', &
            is_u_point=.true., is_interface=.true., xyave_axes=diag_cs%axesZi)
    call define_axes_group_dsamp(diag_cs, (/ id_xh, id_yq, id_zi /), diag_cs%dsamp(dl)%axesCvi, dl, &
            x_cell_method='mean', y_cell_method='point', v_cell_method='point', &
            is_v_point=.true., is_interface=.true., xyave_axes=diag_cs%axesZi)

    ! Axis groupings for 2-D arrays
    call define_axes_group_dsamp(diag_cs, (/ id_xh, id_yh /), diag_cs%dsamp(dl)%axesT1, dl, &
            x_cell_method='mean', y_cell_method='mean', is_h_point=.true.)
    call define_axes_group_dsamp(diag_cs, (/ id_xq, id_yq /), diag_cs%dsamp(dl)%axesB1, dl, &
            x_cell_method='point', y_cell_method='point', is_q_point=.true.)
    call define_axes_group_dsamp(diag_cs, (/ id_xq, id_yh /), diag_cs%dsamp(dl)%axesCu1, dl, &
            x_cell_method='point', y_cell_method='mean', is_u_point=.true.)
    call define_axes_group_dsamp(diag_cs, (/ id_xh, id_yq /), diag_cs%dsamp(dl)%axesCv1, dl, &
            x_cell_method='mean', y_cell_method='point', is_v_point=.true.)

    ! Axis groupings with a non-native vertical coordinate
    if (diag_cs%num_diag_coords>0) then
      allocate(diag_cs%dsamp(dl)%remap_axesTL(diag_cs%num_diag_coords))
      allocate(diag_cs%dsamp(dl)%remap_axesBL(diag_cs%num_diag_coords))
      allocate(diag_cs%dsamp(dl)%remap_axesCuL(diag_cs%num_diag_coords))
      allocate(diag_cs%dsamp(dl)%remap_axesCvL(diag_cs%num_diag_coords))
      allocate(diag_cs%dsamp(dl)%remap_axesTi(diag_cs%num_diag_coords))
      allocate(diag_cs%dsamp(dl)%remap_axesBi(diag_cs%num_diag_coords))
      allocate(diag_cs%dsamp(dl)%remap_axesCui(diag_cs%num_diag_coords))
      allocate(diag_cs%dsamp(dl)%remap_axesCvi(diag_cs%num_diag_coords))
    endif

    do i=1, diag_cs%num_diag_coords
      ! For each possible diagnostic coordinate
      ! call diag_remap_configure_axes(diag_cs%diag_remap_cs(i), G, GV, param_file)

      ! This vertical coordinate has been configured so can be used.
      if (diag_remap_axes_configured(diag_cs%diag_remap_cs(i))) then

        ! This fetches the 1D-axis id for layers and interfaces and overwrite
        ! id_zl and id_zi from above. It also returns the number of layers.
        call diag_remap_get_axes_info(diag_cs%diag_remap_cs(i), nz, id_zL, id_zi)

        ! Axes for z layers
        call define_axes_group_dsamp(diag_cs, (/ id_xh, id_yh, id_zL /), diag_cs%dsamp(dl)%remap_axesTL(i), dl, &
                nz=nz, vertical_coordinate_number=i, &
                x_cell_method='mean', y_cell_method='mean', v_cell_method='mean', &
                is_h_point=.true., is_layer=.true., is_native=.false., needs_remapping=.true., &
                xyave_axes=diag_cs%remap_axesZL(i))

        !! \note Remapping for B points is not yet implemented so needs_remapping is not
        !! provided for remap_axesBL
        call define_axes_group_dsamp(diag_cs, (/ id_xq, id_yq, id_zL /), diag_cs%dsamp(dl)%remap_axesBL(i), dl, &
                nz=nz, vertical_coordinate_number=i, &
                x_cell_method='point', y_cell_method='point', v_cell_method='mean', &
                is_q_point=.true., is_layer=.true., is_native=.false.)

        call define_axes_group_dsamp(diag_cs, (/ id_xq, id_yh, id_zL /), diag_cs%dsamp(dl)%remap_axesCuL(i), dl, &
                nz=nz, vertical_coordinate_number=i, &
                x_cell_method='point', y_cell_method='mean', v_cell_method='mean', &
                is_u_point=.true., is_layer=.true., is_native=.false., needs_remapping=.true., &
                xyave_axes=diag_cs%remap_axesZL(i))

        call define_axes_group_dsamp(diag_cs, (/ id_xh, id_yq, id_zL /), diag_cs%dsamp(dl)%remap_axesCvL(i), dl, &
                nz=nz, vertical_coordinate_number=i, &
                x_cell_method='mean', y_cell_method='point', v_cell_method='mean', &
                is_v_point=.true., is_layer=.true., is_native=.false., needs_remapping=.true., &
                xyave_axes=diag_cs%remap_axesZL(i))

        ! Axes for z interfaces
        call define_axes_group_dsamp(diag_cs, (/ id_xh, id_yh, id_zi /), diag_cs%dsamp(dl)%remap_axesTi(i), dl, &
                nz=nz, vertical_coordinate_number=i, &
                x_cell_method='mean', y_cell_method='mean', v_cell_method='point', &
                is_h_point=.true., is_interface=.true., is_native=.false., needs_interpolating=.true., &
                xyave_axes=diag_cs%remap_axesZi(i))

        !! \note Remapping for B points is not yet implemented so needs_remapping is not provided for remap_axesBi
        call define_axes_group_dsamp(diag_cs, (/ id_xq, id_yq, id_zi /), diag_cs%dsamp(dl)%remap_axesBi(i), dl, &
                nz=nz, vertical_coordinate_number=i, &
                x_cell_method='point', y_cell_method='point', v_cell_method='point', &
                is_q_point=.true., is_interface=.true., is_native=.false.)

        call define_axes_group_dsamp(diag_cs, (/ id_xq, id_yh, id_zi /), diag_cs%dsamp(dl)%remap_axesCui(i), dl, &
                nz=nz, vertical_coordinate_number=i, &
                x_cell_method='point', y_cell_method='mean', v_cell_method='point', &
                is_u_point=.true., is_interface=.true., is_native=.false., &
                needs_interpolating=.true., xyave_axes=diag_cs%remap_axesZi(i))

        call define_axes_group_dsamp(diag_cs, (/ id_xh, id_yq, id_zi /), diag_cs%dsamp(dl)%remap_axesCvi(i), dl, &
                nz=nz, vertical_coordinate_number=i, &
                x_cell_method='mean', y_cell_method='point', v_cell_method='point', &
                is_v_point=.true., is_interface=.true., is_native=.false., &
                needs_interpolating=.true., xyave_axes=diag_cs%remap_axesZi(i))
      endif
    enddo
  enddo

end procedure set_axes_info_dsamp
module procedure set_masks_for_axes
  integer :: c, nk, i, j, k
  type(axes_grp), pointer :: axes => NULL(), h_axes => NULL() ! Current axes, for convenience
  do c=1, diag_cs%num_diag_coords
    ! This vertical coordinate has been configured so can be used.
    if (diag_remap_axes_configured(diag_cs%diag_remap_cs(c))) then

      ! Level/layer h-points in diagnostic coordinate
      axes => diag_cs%remap_axesTL(c)
      nk = axes%nz
      allocate( axes%mask3d(G%isd:G%ied,G%jsd:G%jed,nk), source=0. )
      call diag_remap_calc_hmask(diag_cs%diag_remap_cs(c), G, axes%mask3d)

      h_axes => diag_cs%remap_axesTL(c) ! Use the h-point masks to generate the u-, v- and q- masks

      ! Level/layer u-points in diagnostic coordinate
      axes => diag_cs%remap_axesCuL(c)
      call assert(axes%nz == nk, 'set_masks_for_axes: vertical size mismatch at u-layers')
      call assert(.not. associated(axes%mask3d), 'set_masks_for_axes: already associated')
      allocate( axes%mask3d(G%IsdB:G%IedB,G%jsd:G%jed,nk), source=0. )
      do k = 1, nk ; do j=G%jsc,G%jec ; do I=G%isc-1,G%iec
        if (h_axes%mask3d(i,j,k) + h_axes%mask3d(i+1,j,k) > 0.) axes%mask3d(I,j,k) = 1.
      enddo ; enddo ; enddo

      ! Level/layer v-points in diagnostic coordinate
      axes => diag_cs%remap_axesCvL(c)
      call assert(axes%nz == nk, 'set_masks_for_axes: vertical size mismatch at v-layers')
      call assert(.not. associated(axes%mask3d), 'set_masks_for_axes: already associated')
      allocate( axes%mask3d(G%isd:G%ied,G%JsdB:G%JedB,nk), source=0. )
      do k = 1, nk ; do J=G%jsc-1,G%jec ; do i=G%isc,G%iec
        if (h_axes%mask3d(i,j,k) + h_axes%mask3d(i,j+1,k) > 0.) axes%mask3d(i,J,k) = 1.
      enddo ; enddo ; enddo

      ! Level/layer q-points in diagnostic coordinate
      axes => diag_cs%remap_axesBL(c)
      call assert(axes%nz == nk, 'set_masks_for_axes: vertical size mismatch at q-layers')
      call assert(.not. associated(axes%mask3d), 'set_masks_for_axes: already associated')
      allocate( axes%mask3d(G%IsdB:G%IedB,G%JsdB:G%JedB,nk), source=0. )
      do k = 1, nk ; do J=G%jsc-1,G%jec ; do I=G%isc-1,G%iec
        if (h_axes%mask3d(i,j,k) + h_axes%mask3d(i+1,j+1,k) + &
            h_axes%mask3d(i+1,j,k) + h_axes%mask3d(i,j+1,k) > 0.) axes%mask3d(I,J,k) = 1.
      enddo ; enddo ; enddo

      ! Interface h-points in diagnostic coordinate (w-point)
      axes => diag_cs%remap_axesTi(c)
      call assert(axes%nz == nk, 'set_masks_for_axes: vertical size mismatch at h-interfaces')
      call assert(.not. associated(axes%mask3d), 'set_masks_for_axes: already associated')
      allocate( axes%mask3d(G%isd:G%ied,G%jsd:G%jed,nk+1), source=0. )
      do J=G%jsc-1,G%jec+1 ; do i=G%isc-1,G%iec+1
        if (h_axes%mask3d(i,j,1) > 0.) axes%mask3d(i,J,1) = 1.
        do K = 2, nk
          if (h_axes%mask3d(i,j,k-1) + h_axes%mask3d(i,j,k) > 0.) axes%mask3d(i,J,k) = 1.
        enddo
        if (h_axes%mask3d(i,j,nk) > 0.) axes%mask3d(i,J,nk+1) = 1.
      enddo ; enddo

      h_axes => diag_cs%remap_axesTi(c) ! Use the w-point masks to generate the u-, v- and q- masks

      ! Interface u-points in diagnostic coordinate
      axes => diag_cs%remap_axesCui(c)
      call assert(axes%nz == nk, 'set_masks_for_axes: vertical size mismatch at u-interfaces')
      call assert(.not. associated(axes%mask3d), 'set_masks_for_axes: already associated')
      allocate( axes%mask3d(G%IsdB:G%IedB,G%jsd:G%jed,nk+1), source=0. )
      do k = 1, nk+1 ; do j=G%jsc,G%jec ; do I=G%isc-1,G%iec
        if (h_axes%mask3d(i,j,k) + h_axes%mask3d(i+1,j,k) > 0.) axes%mask3d(I,j,k) = 1.
      enddo ; enddo ; enddo

      ! Interface v-points in diagnostic coordinate
      axes => diag_cs%remap_axesCvi(c)
      call assert(axes%nz == nk, 'set_masks_for_axes: vertical size mismatch at v-interfaces')
      call assert(.not. associated(axes%mask3d), 'set_masks_for_axes: already associated')
      allocate( axes%mask3d(G%isd:G%ied,G%JsdB:G%JedB,nk+1), source=0. )
      do k = 1, nk+1 ; do J=G%jsc-1,G%jec ; do i=G%isc,G%iec
        if (h_axes%mask3d(i,j,k) + h_axes%mask3d(i,j+1,k) > 0.) axes%mask3d(i,J,k) = 1.
      enddo ; enddo ; enddo

      ! Interface q-points in diagnostic coordinate
      axes => diag_cs%remap_axesBi(c)
      call assert(axes%nz == nk, 'set_masks_for_axes: vertical size mismatch at q-interfaces')
      call assert(.not. associated(axes%mask3d), 'set_masks_for_axes: already associated')
      allocate( axes%mask3d(G%IsdB:G%IedB,G%JsdB:G%JedB,nk+1), source=0. )
      do k = 1, nk ; do J=G%jsc-1,G%jec ; do I=G%isc-1,G%iec
        if (h_axes%mask3d(i,j,k) + h_axes%mask3d(i+1,j+1,k) + &
            h_axes%mask3d(i+1,j,k) + h_axes%mask3d(i,j+1,k) > 0.) axes%mask3d(I,J,k) = 1.
      enddo ; enddo ; enddo
    endif
  enddo

  ! Allocate and initialize the downsampled masks for the axes
  call set_masks_for_axes_dsamp(G, diag_cs)

end procedure set_masks_for_axes
module procedure set_masks_for_axes_dsamp
  integer :: c, dl, dlfac
  type(axes_grp), pointer :: axes => NULL() ! Current axes, for convenience
  do dl=1, diag_cs%num_diag_dsamp_levels
    dlfac = diag_cs%diag_dsamp_levels(dl) ! The actual downsampling factor for this level
    do c=1, diag_cs%num_diag_coords
      ! Level/layer h-points in diagnostic coordinate
      axes => diag_cs%remap_axesTL(c)
      call downsample_mask(axes%mask3d, diag_cs%dsamp(dl)%remap_axesTL(c)%dsamp(dl)%mask3d, &
              dlfac, xyz_method(axes), G%isc, G%jsc, G%isd, G%jsd, &
              G%HId(dl)%isc, G%HId(dl)%iec, G%HId(dl)%jsc, G%HId(dl)%jec, G%HId(dl)%isd, G%HId(dl)%ied, &
              G%HId(dl)%jsd, G%HId(dl)%jed)
      diag_cs%dsamp(dl)%remap_axesTL(c)%mask3d => axes%mask3d ! Set a pointer to the non-downsampled mask
      ! Level/layer u-points in diagnostic coordinate
      axes => diag_cs%remap_axesCuL(c)
      call downsample_mask(axes%mask3d, diag_cs%dsamp(dl)%remap_axesCuL(c)%dsamp(dl)%mask3d, &
              dlfac, xyz_method(axes), G%IscB, G%jsc, G%IsdB, G%jsd, &
              G%HId(dl)%IscB, G%HId(dl)%IecB, G%HId(dl)%jsc, G%HId(dl)%jec, G%HId(dl)%IsdB, G%HId(dl)%IedB, &
              G%HId(dl)%jsd, G%HId(dl)%jed)
      diag_cs%dsamp(dl)%remap_axesCul(c)%mask3d => axes%mask3d ! Set a pointer to the non-downsampled mask
      ! Level/layer v-points in diagnostic coordinate
      axes => diag_cs%remap_axesCvL(c)
      call downsample_mask(axes%mask3d, diag_cs%dsamp(dl)%remap_axesCvL(c)%dsamp(dl)%mask3d, &
              dlfac, xyz_method(axes), G%isc, G%JscB, G%isd, G%JsdB, &
              G%HId(dl)%isc, G%HId(dl)%iec, G%HId(dl)%JscB, G%HId(dl)%JecB, G%HId(dl)%isd, G%HId(dl)%ied, &
              G%HId(dl)%JsdB, G%HId(dl)%JedB)
      diag_cs%dsamp(dl)%remap_axesCvL(c)%mask3d => axes%mask3d ! Set a pointer to the non-downsampled mask
      ! Level/layer q-points in diagnostic coordinate
      axes => diag_cs%remap_axesBL(c)
      call downsample_mask(axes%mask3d, diag_cs%dsamp(dl)%remap_axesBL(c)%dsamp(dl)%mask3d, &
              dlfac, xyz_method(axes), G%IscB, G%JscB, G%IsdB, G%JsdB, &
              G%HId(dl)%IscB, G%HId(dl)%IecB, G%HId(dl)%JscB, G%HId(dl)%JecB, G%HId(dl)%IsdB, G%HId(dl)%IedB, &
              G%HId(dl)%JsdB, G%HId(dl)%JedB)
      diag_cs%dsamp(dl)%remap_axesBL(c)%mask3d => axes%mask3d ! Set a pointer to the non-downsampled mask
      ! Interface h-points in diagnostic coordinate (w-point)
      axes => diag_cs%remap_axesTi(c)
      call downsample_mask(axes%mask3d, diag_cs%dsamp(dl)%remap_axesTi(c)%dsamp(dl)%mask3d,  &
              dlfac, xyz_method(axes), G%isc, G%jsc, G%isd, G%jsd, &
              G%HId(dl)%isc, G%HId(dl)%iec, G%HId(dl)%jsc, G%HId(dl)%jec, G%HId(dl)%isd, G%HId(dl)%ied, &
              G%HId(dl)%jsd, G%HId(dl)%jed)
      diag_cs%dsamp(dl)%remap_axesTi(c)%mask3d => axes%mask3d ! Set a pointer to the non-downsampled mask
      ! Interface u-points in diagnostic coordinate
      axes => diag_cs%remap_axesCui(c)
      call downsample_mask(axes%mask3d, diag_cs%dsamp(dl)%remap_axesCui(c)%dsamp(dl)%mask3d,  &
              dlfac, xyz_method(axes), G%IscB, G%jsc, G%IsdB, G%jsd, &
              G%HId(dl)%IscB, G%HId(dl)%IecB, G%HId(dl)%jsc, G%HId(dl)%jec, G%HId(dl)%IsdB, G%HId(dl)%IedB, &
              G%HId(dl)%jsd, G%HId(dl)%jed)
      diag_cs%dsamp(dl)%remap_axesCui(c)%mask3d => axes%mask3d ! Set a pointer to the non-downsampled mask
      ! Interface v-points in diagnostic coordinate
      axes => diag_cs%remap_axesCvi(c)
      call downsample_mask(axes%mask3d, diag_cs%dsamp(dl)%remap_axesCvi(c)%dsamp(dl)%mask3d,  &
              dlfac, xyz_method(axes), G%isc, G%JscB, G%isd, G%JsdB, &
              G%HId(dl)%isc, G%HId(dl)%iec, G%HId(dl)%JscB, G%HId(dl)%JecB, G%HId(dl)%isd, G%HId(dl)%ied, &
              G%HId(dl)%JsdB, G%HId(dl)%JedB)
      diag_cs%dsamp(dl)%remap_axesCvi(c)%mask3d => axes%mask3d ! Set a pointer to the non-downsampled mask
      ! Interface q-points in diagnostic coordinate
      axes => diag_cs%remap_axesBi(c)
      call downsample_mask(axes%mask3d, diag_cs%dsamp(dl)%remap_axesBi(c)%dsamp(dl)%mask3d,  &
              dlfac, xyz_method(axes), G%IscB, G%JscB, G%IsdB, G%JsdB, &
              G%HId(dl)%IscB, G%HId(dl)%IecB, G%HId(dl)%JscB, G%HId(dl)%JecB, G%HId(dl)%IsdB, G%HId(dl)%IedB, &
              G%HId(dl)%JsdB, G%HId(dl)%JedB)
      diag_cs%dsamp(dl)%remap_axesBi(c)%mask3d => axes%mask3d ! Set a pointer to the non-downsampled mask
    enddo
  enddo
end procedure set_masks_for_axes_dsamp
module procedure diag_register_area_ids
  integer :: fms_id, i
  if (present(id_area_t)) then
    fms_id = diag_cs%diags(id_area_t)%fms_diag_id
    diag_cs%axesT1%id_area = fms_id
    diag_cs%axesTi%id_area = fms_id
    diag_cs%axesTL%id_area = fms_id
    do i=1, diag_cs%num_diag_coords
      diag_cs%remap_axesTL(i)%id_area = fms_id
      diag_cs%remap_axesTi(i)%id_area = fms_id
    enddo
  endif
  if (present(id_area_q)) then
    fms_id = diag_cs%diags(id_area_q)%fms_diag_id
    diag_cs%axesB1%id_area = fms_id
    diag_cs%axesBi%id_area = fms_id
    diag_cs%axesBL%id_area = fms_id
    do i=1, diag_cs%num_diag_coords
      diag_cs%remap_axesBL(i)%id_area = fms_id
      diag_cs%remap_axesBi(i)%id_area = fms_id
    enddo
  endif
end procedure diag_register_area_ids
module procedure register_cell_measure
  integer :: id
  id = register_diag_field('ocean_model', 'volcello', diag%axesTL, &
                           Time, 'Ocean grid-cell volume', units='m3', conversion=1.0, &
                           standard_name='ocean_volume', v_extensive=.true., &
                           x_cell_method='sum', y_cell_method='sum')
  call diag_associate_volume_cell_measure(diag, id)

end procedure register_cell_measure
module procedure diag_associate_volume_cell_measure
  type(diag_type), pointer :: tmp => NULL()
  if (id_h_volume<=0) return ! Do nothing
  diag_cs%volume_cell_measure_dm_id = id_h_volume ! Record for diag_get_volume_cell_measure_dm_id()

  ! Set the cell measure for this axes group to the FMS id in this coordinate system
  diag_cs%diags(id_h_volume)%axes%id_volume = diag_cs%diags(id_h_volume)%fms_diag_id

  tmp => diag_cs%diags(id_h_volume)%next ! First item in the list, if any
  do while (associated(tmp))
    ! Set the cell measure for this axes group to the FMS id in this coordinate system
    tmp%axes%id_volume = tmp%fms_diag_id
    tmp => tmp%next ! Move to next axes group for this field
  enddo

end procedure diag_associate_volume_cell_measure
module procedure diag_get_volume_cell_measure_dm_id
  diag_get_volume_cell_measure_dm_id = diag_cs%volume_cell_measure_dm_id

end procedure diag_get_volume_cell_measure_dm_id
module procedure define_axes_group
  integer :: n
  n = size(handles)
  if (n<1 .or. n>3) call MOM_error(FATAL, "define_axes_group: wrong size for list of handles!")
  allocate( axes%handles(n) )
  axes%id = ints_to_string(handles, max(n,3)) ! Identifying string
  axes%rank = n
  axes%handles(:) = handles(:)
  axes%diag_cs => diag_cs ! A (circular) link back to the diag_cs structure
  if (present(x_cell_method)) then
    if (axes%rank<2) call MOM_error(FATAL, 'define_axes_group: ' // &
                                           'Can not set x_cell_method for rank<2.')
    axes%x_cell_method = trim(x_cell_method)
  else
    axes%x_cell_method = ''
  endif
  if (present(y_cell_method)) then
    if (axes%rank<2) call MOM_error(FATAL, 'define_axes_group: ' // &
                                           'Can not set y_cell_method for rank<2.')
    axes%y_cell_method = trim(y_cell_method)
  else
    axes%y_cell_method = ''
  endif
  if (present(v_cell_method)) then
    if (axes%rank/=1 .and. axes%rank/=3) call MOM_error(FATAL, 'define_axes_group: ' // &
                                           'Can not set v_cell_method for rank<>1 or 3.')
    axes%v_cell_method = trim(v_cell_method)
  else
    axes%v_cell_method = ''
  endif

  if (present(nz)) axes%nz = nz
  if (present(vertical_coordinate_number)) axes%vertical_coordinate_number = vertical_coordinate_number
  if (present(is_h_point)) axes%is_h_point = is_h_point
  if (present(is_q_point)) axes%is_q_point = is_q_point
  if (present(is_u_point)) axes%is_u_point = is_u_point
  if (present(is_v_point)) axes%is_v_point = is_v_point
  if (present(is_layer)) axes%is_layer = is_layer
  if (present(is_interface)) axes%is_interface = is_interface
  if (present(is_native)) axes%is_native = is_native
  if (present(needs_remapping)) axes%needs_remapping = needs_remapping
  if (present(needs_interpolating)) axes%needs_interpolating = needs_interpolating
  if (present(xyave_axes)) axes%xyave_axes => xyave_axes

  ! Setup masks for this axes group
  axes%mask2d => null()
  if (axes%rank==2) then
    if (axes%is_h_point) axes%mask2d => diag_cs%mask2dT
    if (axes%is_u_point) axes%mask2d => diag_cs%mask2dCu
    if (axes%is_v_point) axes%mask2d => diag_cs%mask2dCv
    if (axes%is_q_point) axes%mask2d => diag_cs%mask2dBu
  endif
  ! A static 3d mask for non-native coordinates can only be setup when a grid is available
  axes%mask3d => null()
  if (axes%rank==3 .and. axes%is_native) then
    ! Native variables can/should use the native masks copied into diag_cs
    if (axes%is_layer) then
      if (axes%is_h_point) axes%mask3d => diag_cs%mask3dTL
      if (axes%is_u_point) axes%mask3d => diag_cs%mask3dCuL
      if (axes%is_v_point) axes%mask3d => diag_cs%mask3dCvL
      if (axes%is_q_point) axes%mask3d => diag_cs%mask3dBL
    elseif (axes%is_interface) then
      if (axes%is_h_point) axes%mask3d => diag_cs%mask3dTi
      if (axes%is_u_point) axes%mask3d => diag_cs%mask3dCui
      if (axes%is_v_point) axes%mask3d => diag_cs%mask3dCvi
      if (axes%is_q_point) axes%mask3d => diag_cs%mask3dBi
    endif
  endif


end procedure define_axes_group
module procedure define_axes_group_dsamp
  integer :: n
  n = size(handles)
  if (n<1 .or. n>3) call MOM_error(FATAL, "define_axes_group: wrong size for list of handles!")
  allocate( axes%handles(n) )
  axes%id = ints_to_string(handles, max(n,3)) ! Identifying string
  axes%rank = n
  axes%handles(:) = handles(:)
  axes%diag_cs => diag_cs ! A (circular) link back to the diag_cs structure
  if (present(x_cell_method)) then
    if (axes%rank<2) call MOM_error(FATAL, 'define_axes_group: ' // &
                                           'Can not set x_cell_method for rank<2.')
    axes%x_cell_method = trim(x_cell_method)
  else
    axes%x_cell_method = ''
  endif
  if (present(y_cell_method)) then
    if (axes%rank<2) call MOM_error(FATAL, 'define_axes_group: ' // &
                                           'Can not set y_cell_method for rank<2.')
    axes%y_cell_method = trim(y_cell_method)
  else
    axes%y_cell_method = ''
  endif
  if (present(v_cell_method)) then
    if (axes%rank/=1 .and. axes%rank/=3) call MOM_error(FATAL, 'define_axes_group: ' // &
                                           'Can not set v_cell_method for rank<>1 or 3.')
    axes%v_cell_method = trim(v_cell_method)
  else
    axes%v_cell_method = ''
  endif
  axes%downsample_level_index = dl
  axes%downsample_level_factor = diag_cs%diag_dsamp_levels(dl)
  if (present(nz)) axes%nz = nz
  if (present(vertical_coordinate_number)) axes%vertical_coordinate_number = vertical_coordinate_number
  if (present(is_h_point)) axes%is_h_point = is_h_point
  if (present(is_q_point)) axes%is_q_point = is_q_point
  if (present(is_u_point)) axes%is_u_point = is_u_point
  if (present(is_v_point)) axes%is_v_point = is_v_point
  if (present(is_layer)) axes%is_layer = is_layer
  if (present(is_interface)) axes%is_interface = is_interface
  if (present(is_native)) axes%is_native = is_native
  if (present(needs_remapping)) axes%needs_remapping = needs_remapping
  if (present(needs_interpolating)) axes%needs_interpolating = needs_interpolating
  if (present(xyave_axes)) axes%xyave_axes => xyave_axes

  ! Setup masks for this axes group

  axes%mask2d => null()
  if (axes%rank==2) then
    if (axes%is_h_point) axes%mask2d => diag_cs%mask2dT
    if (axes%is_u_point) axes%mask2d => diag_cs%mask2dCu
    if (axes%is_v_point) axes%mask2d => diag_cs%mask2dCv
    if (axes%is_q_point) axes%mask2d => diag_cs%mask2dBu
  endif
  ! A static 3d mask for non-native coordinates can only be setup when a grid is available
  axes%mask3d => null()
  if (axes%rank==3 .and. axes%is_native) then
    ! Native variables can/should use the native masks copied into diag_cs
    if (axes%is_layer) then
      if (axes%is_h_point) axes%mask3d => diag_cs%mask3dTL
      if (axes%is_u_point) axes%mask3d => diag_cs%mask3dCuL
      if (axes%is_v_point) axes%mask3d => diag_cs%mask3dCvL
      if (axes%is_q_point) axes%mask3d => diag_cs%mask3dBL
    elseif (axes%is_interface) then
      if (axes%is_h_point) axes%mask3d => diag_cs%mask3dTi
      if (axes%is_u_point) axes%mask3d => diag_cs%mask3dCui
      if (axes%is_v_point) axes%mask3d => diag_cs%mask3dCvi
      if (axes%is_q_point) axes%mask3d => diag_cs%mask3dBi
    endif
  endif

  if (.Not. allocated(axes%dsamp)) allocate(axes%dsamp(diag_cs%num_diag_dsamp_levels))
  axes%dsamp(dl)%mask2d => null()
  if (axes%rank==2) then
    if (axes%is_h_point) axes%dsamp(dl)%mask2d => diag_cs%dsamp(dl)%mask2dT
    if (axes%is_u_point) axes%dsamp(dl)%mask2d => diag_cs%dsamp(dl)%mask2dCu
    if (axes%is_v_point) axes%dsamp(dl)%mask2d => diag_cs%dsamp(dl)%mask2dCv
    if (axes%is_q_point) axes%dsamp(dl)%mask2d => diag_cs%dsamp(dl)%mask2dBu
  endif
  ! A static 3d mask for non-native coordinates can only be setup when a grid is available
  axes%dsamp(dl)%mask3d => null()
  if (axes%rank==3 .and. axes%is_native) then
    ! Native variables can/should use the native masks copied into diag_cs
    if (axes%is_layer) then
      if (axes%is_h_point) axes%dsamp(dl)%mask3d => diag_cs%dsamp(dl)%mask3dTL
      if (axes%is_u_point) axes%dsamp(dl)%mask3d => diag_cs%dsamp(dl)%mask3dCuL
      if (axes%is_v_point) axes%dsamp(dl)%mask3d => diag_cs%dsamp(dl)%mask3dCvL
      if (axes%is_q_point) axes%dsamp(dl)%mask3d => diag_cs%dsamp(dl)%mask3dBL
    elseif (axes%is_interface) then
      if (axes%is_h_point) axes%dsamp(dl)%mask3d => diag_cs%dsamp(dl)%mask3dTi
      if (axes%is_u_point) axes%dsamp(dl)%mask3d => diag_cs%dsamp(dl)%mask3dCui
      if (axes%is_v_point) axes%dsamp(dl)%mask3d => diag_cs%dsamp(dl)%mask3dCvi
      if (axes%is_q_point) axes%dsamp(dl)%mask3d => diag_cs%dsamp(dl)%mask3dBi
    endif
  endif

end procedure define_axes_group_dsamp
module procedure set_diag_mediator_grid
  diag_cs%is = G%isc - (G%isd-1) ; diag_cs%ie = G%iec - (G%isd-1)
  diag_cs%js = G%jsc - (G%jsd-1) ; diag_cs%je = G%jec - (G%jsd-1)
  diag_cs%isd = G%isd ; diag_cs%ied = G%ied
  diag_cs%jsd = G%jsd ; diag_cs%jed = G%jed

end procedure set_diag_mediator_grid
module procedure post_data_0d
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
              'post_data_0d: Unregistered diagnostic id')
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
end procedure post_data_0d
module procedure post_data_1d_k
  logical :: used  ! The return value of send_data is not used for anything.
  real, dimension(:), pointer :: locfield => NULL() ! The field being offered in arbitrary unscaled units [a]
  logical :: is_stat
  integer :: k, ks, ke
  type(diag_type), pointer :: diag => null()
  integer :: time_days
  integer :: time_seconds
  character(len=300) :: debug_mesg
  if (id_clock_diag_mediator>0) call cpu_clock_begin(id_clock_diag_mediator)
  is_stat = .false. ; if (present(is_static)) is_stat = is_static

  ! Iterate over list of diag 'variants', e.g. CMOR aliases.
  call assert(diag_field_id < diag_cs%next_free_diag_id, &
              'post_data_1d_k: Unregistered diagnostic id')
  diag => diag_cs%diags(diag_field_id)
  do while (associated(diag))

    if ((diag%conversion_factor /= 0.) .and. (diag%conversion_factor /= 1.)) then
      ks = lbound(field,1) ; ke = ubound(field,1)
      allocate( locfield( ks:ke ) )

      do k=ks,ke
        locfield(k) = field(k) * diag%conversion_factor
      enddo
    else
      locfield => field
    endif

    if (diag_cs%diag_as_chksum) then
      ! Append timestep to mesg
      call get_time(diag_cs%time_end, time_seconds, days=time_days)
      write(debug_mesg, '(a, 1x, i0, 1x, i0)') &
          trim(diag%debug_str), time_days, time_seconds

      call zchksum(locfield, debug_mesg, logunit=diag_cs%chksum_iounit)
    elseif (is_stat) then
      used = send_data_infra(diag%fms_diag_id, locfield)
    elseif (diag_cs%ave_enabled) then
      used = send_data_infra(diag%fms_diag_id, locfield, time=diag_cs%time_end, weight=diag_cs%time_int)
    endif
    if ((diag%conversion_factor /= 0.) .and. (diag%conversion_factor /= 1.)) deallocate( locfield )

    diag => diag%next
  enddo

  if (id_clock_diag_mediator>0) call cpu_clock_end(id_clock_diag_mediator)
end procedure post_data_1d_k
module procedure post_data_2d
  type(diag_type), pointer :: diag => null()
  if (id_clock_diag_mediator>0) call cpu_clock_begin(id_clock_diag_mediator)

  ! Iterate over list of diag 'variants' (e.g. CMOR aliases) and post each.
  call assert(diag_field_id < diag_cs%next_free_diag_id, &
              'post_data_2d: Unregistered diagnostic id')
  diag => diag_cs%diags(diag_field_id)
  do while (associated(diag))
    call post_data_2d_low(diag, field, diag_cs, is_static, mask)
    diag => diag%next
  enddo

  if (id_clock_diag_mediator>0) call cpu_clock_end(id_clock_diag_mediator)
end procedure post_data_2d
module procedure post_data_2d_low
  real, dimension(:,:), pointer :: locfield ! The field being offered in arbitrary unscaled units [a]
  real, dimension(:,:), pointer :: locmask  ! A pointer to the data mask to use [nondim]
  logical :: used  ! The return value of send_data is not used for anything.
  logical :: is_stat, not_static
  integer :: cszi, cszj, dszi, dszj
  integer :: isv, iev, jsv, jev, i, j, isv_o, jsv_o
  real, dimension(:,:), allocatable, target :: locfield_dsamp ! A downsampled version of locfield [a]
  real, dimension(:,:), allocatable, target :: locmask_dsamp  ! A downsampled version of locmask [nondim]
  real, dimension(:,:), pointer :: ones => NULL() ! An array of ones for testing where masks do not apply [nondim]
  real, dimension(:,:), pointer :: mask_in => NULL() ! A pointer to the input mask [nondim]
  integer :: dl, dlfac
  integer :: time_days
  integer :: time_seconds
  character(len=300) :: mesg
  character(len=300) :: debug_mesg
  locfield => NULL()
  locmask => NULL()
  is_stat = .false. ; if (present(is_static)) is_stat = is_static
  not_static = .not. is_stat

  ! Determine the proper array indices, noting that because of the (:,:)
  ! declaration of field, symmetric arrays are using a SW-grid indexing,
  ! but non-symmetric arrays are using a NE-grid indexing.  Send_data
  ! actually only uses the difference between ie and is to determine
  ! the output data size and assumes that halos are symmetric.
  isv = diag_cs%is ; iev = diag_cs%ie ; jsv = diag_cs%js ; jev = diag_cs%je

  cszi = diag_cs%ie-diag_cs%is +1 ; dszi = diag_cs%ied-diag_cs%isd +1
  cszj = diag_cs%je-diag_cs%js +1 ; dszj = diag_cs%jed-diag_cs%jsd +1
  if ( size(field,1) == dszi ) then
    isv = diag_cs%is ; iev = diag_cs%ie     ! Data domain
  elseif ( size(field,1) == dszi + 1 ) then
    isv = diag_cs%is ; iev = diag_cs%ie+1   ! Symmetric data domain
  elseif ( size(field,1) == cszi) then
    isv = 1 ; iev = cszi                    ! Computational domain
  elseif ( size(field,1) == cszi + 1 ) then
    isv = 1 ; iev = cszi+1                  ! Symmetric computational domain
  else
    write (mesg,*) " peculiar size ",size(field,1)," in i-direction\n"//&
       "does not match one of ", cszi, cszi+1, dszi, dszi+1
    call MOM_error(FATAL,"post_data_2d_low: "//trim(diag%debug_str)//trim(mesg))
  endif

  if ( size(field,2) == dszj ) then
    jsv = diag_cs%js ; jev = diag_cs%je     ! Data domain
  elseif ( size(field,2) == dszj + 1 ) then
    jsv = diag_cs%js ; jev = diag_cs%je+1   ! Symmetric data domain
  elseif ( size(field,2) == cszj ) then
    jsv = 1 ; jev = cszj                    ! Computational domain
  elseif ( size(field,2) == cszj+1 ) then
    jsv = 1 ; jev = cszj+1                  ! Symmetric computational domain
  else
    write (mesg,*) " peculiar size ",size(field,2)," in j-direction\n"//&
       "does not match one of ", cszj, cszj+1, dszj, dszj+1
    call MOM_error(FATAL,"post_data_2d_low: "//trim(diag%debug_str)//trim(mesg))
  endif

  if ((diag%conversion_factor /= 0.) .and. (diag%conversion_factor /= 1.)) then
    allocate( locfield( lbound(field,1):ubound(field,1), lbound(field,2):ubound(field,2) ) )
    do j=jsv,jev ; do i=isv,iev
      locfield(i,j) = field(i,j) * diag%conversion_factor
    enddo ; enddo
  else
    locfield => field
  endif

  if (present(mask)) then
    locmask => mask
  elseif (not_static .and. associated(diag%axes)) then
  ! If we were to decide to allow masking of static diagnostics, we could do so by changing the line above to
  ! elseif (associated(diag%axes) .and. (diag_CS%mask_static_diags .or. not_static)) then
    if (associated(diag%axes%mask2d)) locmask => diag%axes%mask2d
  endif

  dlfac = 1
  if (not_static .and. associated(diag%axes)) &
    dlfac = diag%axes%downsample_level_factor ! Static field downsampling is not supported yet.
  ! Downsample the diag field and mask as appropriate.
  if (dlfac > 1) then
    dl = diag%axes%downsample_level_index
    isv_o = isv ; jsv_o = jsv
    call downsample_diag_field(locfield, locfield_dsamp, dl, diag_cs, diag, isv, iev, jsv, jev, mask)
    if ((diag%conversion_factor /= 0.) .and. (diag%conversion_factor /= 1.)) deallocate( locfield )
    locfield => locfield_dsamp
    if (present(mask)) then
      ! Replicate the downsampling of other fields to find unmasked points.
      allocate(ones, mold=locmask) ; ones(:,:) = 1.0
      mask_in => mask
      call downsample_field_2d(ones, locmask_dsamp, dlfac, diag%xyz_method, mask_in, diag_cs, diag, &
                               isv_o, jsv_o, isv, iev, jsv, jev)
      deallocate(ones)
      where (abs(locmask_dsamp) > 0.0) locmask_dsamp = 1.0
      locmask => locmask_dsamp
    elseif (associated(diag%axes%dsamp(dl)%mask2d)) then
      locmask => diag%axes%dsamp(dl)%mask2d
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
  if ((diag%conversion_factor /= 0.) .and. (diag%conversion_factor /= 1.) .and. dlfac<2) &
    deallocate( locfield )
end procedure post_data_2d_low
module procedure post_data_3d
  type(diag_type), pointer :: diag => null()
  real, dimension(:,:,:), allocatable :: remapped_field !< The vertically remapped diagnostic [A ~> a]
  logical :: staggered_in_x, staggered_in_y, dz_diag_needed, dz_begin_needed
  real, dimension(:,:,:), pointer :: h_diag => NULL() !< A pointer to the thickness to use for vertically
  real, dimension(diag_cs%G%isd:diag_cS%G%ied, diag_cs%G%jsd:diag_cS%G%jed, diag_cs%GV%ke) :: &
    dz_diag     ! Layer vertical extents for remapping [Z ~> m]
  if (id_clock_diag_mediator>0) call cpu_clock_begin(id_clock_diag_mediator)

  ! For intensive variables only, we can choose to use a different diagnostic grid to map to
  if (present(alt_h)) then
    h_diag => alt_h
  else
    h_diag => diag_cs%h
  endif

  ! Iterate over list of diag 'variants', e.g. CMOR aliases, different vertical
  ! grids, and post each.
  call assert(diag_field_id < diag_cs%next_free_diag_id, &
              'post_data_3d: Unregistered diagnostic id')

  if (diag_cs%show_call_tree) &
    call callTree_enter("post_data_3d("//trim(diag_cs%diags(diag_field_id)%debug_str)//")")

  ! Find out whether there are any z-based diagnostics
  diag => diag_cs%diags(diag_field_id)
  dz_diag_needed = .false.
  do while (associated(diag))
    if (diag%axes%needs_remapping .or. diag%axes%needs_interpolating) then
      if (diag_cs%diag_remap_cs(diag%axes%vertical_coordinate_number)%Z_based_coord) &
        dz_diag_needed = .true.
    endif
    diag => diag%next
  enddo

  ! Determine the diagnostic grid spacing in height units, if it is needed.
  if (dz_diag_needed) then
    call thickness_to_dz(h_diag, diag_cs%tv, dz_diag, diag_cs%G, diag_cs%GV, diag_cs%US, halo_size=1)
  endif

  diag => diag_cs%diags(diag_field_id)
  do while (associated(diag))
    call assert(associated(diag%axes), 'post_data_3d: axes is not associated')

    staggered_in_x = diag%axes%is_u_point .or. diag%axes%is_q_point
    staggered_in_y = diag%axes%is_v_point .or. diag%axes%is_q_point

    if (diag%v_extensive .and. .not.diag%axes%is_native) then
      ! The field is vertically integrated and needs to be re-gridded
      if (present(mask)) then
        call MOM_error(FATAL,"post_data_3d: no mask for regridded field.")
      endif

      if (id_clock_diag_remap>0) call cpu_clock_begin(id_clock_diag_remap)
      allocate(remapped_field(size(field,1), size(field,2), diag%axes%nz))
      if (diag_cs%diag_remap_cs(diag%axes%vertical_coordinate_number)%Z_based_coord) then
        call vertically_reintegrate_diag_field(                                    &
                diag_cs%diag_remap_cs(diag%axes%vertical_coordinate_number), diag_cs%G, &
                diag_cs%dz_begin, diag_cs%diag_remap_cs(diag%axes%vertical_coordinate_number)%h_extensive, &
                diag_cs%OBC_u, diag_cs%OBC_v, staggered_in_x, staggered_in_y, diag%axes%mask3d, field, remapped_field)
      else
        call vertically_reintegrate_diag_field(                                    &
                diag_cs%diag_remap_cs(diag%axes%vertical_coordinate_number), diag_cs%G, &
                diag_cs%h_begin, diag_cs%diag_remap_cs(diag%axes%vertical_coordinate_number)%h_extensive, &
                diag_cs%OBC_u, diag_cs%OBC_v, staggered_in_x, staggered_in_y, diag%axes%mask3d, field, remapped_field)
      endif
      if (id_clock_diag_remap>0) call cpu_clock_end(id_clock_diag_remap)
      if (associated(diag%axes%mask3d)) then
        ! Since 3d masks do not vary in the vertical, just use as much as is
        ! needed.
        call post_data_3d_low(diag, remapped_field, diag_cs, is_static, &
                              mask=diag%axes%mask3d)
      else
        call post_data_3d_low(diag, remapped_field, diag_cs, is_static)
      endif
      if (id_clock_diag_remap>0) call cpu_clock_begin(id_clock_diag_remap)
      deallocate(remapped_field)
      if (id_clock_diag_remap>0) call cpu_clock_end(id_clock_diag_remap)
    elseif (diag%axes%needs_remapping) then
      ! Remap this field to another vertical coordinate.
      if (present(mask)) then
        call MOM_error(FATAL,"post_data_3d: no mask for regridded field.")
      endif

      if (id_clock_diag_remap>0) call cpu_clock_begin(id_clock_diag_remap)
      allocate(remapped_field(size(field,1), size(field,2), diag%axes%nz))
      if (diag_cs%diag_remap_cs(diag%axes%vertical_coordinate_number)%Z_based_coord) then
        call diag_remap_do_remap(diag_cs%diag_remap_cs(diag%axes%vertical_coordinate_number), &
                diag_cs%G, diag_cs%GV, diag_cs%US, dz_diag, diag_cs%OBC_u, diag_cs%OBC_v, &
                staggered_in_x, staggered_in_y, diag%axes%mask3d, field, remapped_field)
      else
        call diag_remap_do_remap(diag_cs%diag_remap_cs(diag%axes%vertical_coordinate_number), &
                diag_cs%G, diag_cs%GV, diag_cs%US, h_diag, diag_cs%OBC_u, diag_cs%OBC_v, &
                staggered_in_x, staggered_in_y, diag%axes%mask3d, field, remapped_field)
      endif
      if (id_clock_diag_remap>0) call cpu_clock_end(id_clock_diag_remap)
      if (associated(diag%axes%mask3d)) then
        ! Since 3d masks do not vary in the vertical, just use as much as is
        ! needed.
        call post_data_3d_low(diag, remapped_field, diag_cs, is_static, &
                              mask=diag%axes%mask3d)
      else
        call post_data_3d_low(diag, remapped_field, diag_cs, is_static)
      endif
      if (id_clock_diag_remap>0) call cpu_clock_begin(id_clock_diag_remap)
      deallocate(remapped_field)
      if (id_clock_diag_remap>0) call cpu_clock_end(id_clock_diag_remap)
    elseif (diag%axes%needs_interpolating) then
      ! Interpolate this field to another vertical coordinate.
      if (present(mask)) then
        call MOM_error(FATAL,"post_data_3d: no mask for regridded field.")
      endif

      if (id_clock_diag_remap>0) call cpu_clock_begin(id_clock_diag_remap)
      allocate(remapped_field(size(field,1), size(field,2), diag%axes%nz+1))
      if (diag_cs%diag_remap_cs(diag%axes%vertical_coordinate_number)%Z_based_coord) then
        call vertically_interpolate_diag_field(diag_cs%diag_remap_cs(diag%axes%vertical_coordinate_number), &
                diag_cs%G, dz_diag, diag_cs%OBC_u, diag_cs%OBC_v, staggered_in_x, staggered_in_y, &
                diag%axes%mask3d, field, remapped_field)
      else
        call vertically_interpolate_diag_field(diag_cs%diag_remap_cs(diag%axes%vertical_coordinate_number), &
                diag_cs%G, h_diag, diag_cs%OBC_u, diag_cs%OBC_v, staggered_in_x, staggered_in_y, &
                diag%axes%mask3d, field, remapped_field)
      endif
      if (id_clock_diag_remap>0) call cpu_clock_end(id_clock_diag_remap)
      if (associated(diag%axes%mask3d)) then
        ! Since 3d masks do not vary in the vertical, just use as much as is needed.
        call post_data_3d_low(diag, remapped_field, diag_cs, is_static, &
                              mask=diag%axes%mask3d)
      else
        call post_data_3d_low(diag, remapped_field, diag_cs, is_static)
      endif
      if (id_clock_diag_remap>0) call cpu_clock_begin(id_clock_diag_remap)
      deallocate(remapped_field)
      if (id_clock_diag_remap>0) call cpu_clock_end(id_clock_diag_remap)
    else
      call post_data_3d_low(diag, field, diag_cs, is_static, mask)
    endif
    diag => diag%next
  enddo
  if (id_clock_diag_mediator>0) call cpu_clock_end(id_clock_diag_mediator)

  if (diag_cs%show_call_tree) &
    call callTree_leave("post_data_3d("//trim(diag_cs%diags(diag_field_id)%debug_str)//")")

end procedure post_data_3d
module procedure post_data_3d_low
  real, dimension(:,:,:), pointer :: locfield ! The field being offered in arbitrary unscaled units [a]
  real, dimension(:,:,:), pointer :: locmask  ! A pointer to the data mask to use [nondim]
  character(len=300) :: mesg
  logical :: used  ! The return value of send_data is not used for anything.
  logical :: staggered_in_x, staggered_in_y
  logical :: is_stat, not_static
  integer :: cszi, cszj, dszi, dszj
  integer :: isv, iev, jsv, jev, ks, ke, i, j, k, isv_c, jsv_c, isv_o, jsv_o
  real, dimension(:,:,:), allocatable, target :: locfield_dsamp ! A downsampled version of locfield [a]
  real, dimension(:,:,:), allocatable, target :: locmask_dsamp  ! A downsampled version of locmask [nondim]
  real, dimension(:,:,:), pointer :: ones => NULL() ! An array of ones for testing where masks do not apply [nondim]
  real, dimension(:,:,:), pointer :: mask_in => NULL() ! A pointer to the input mask [nondim]
  integer :: dl, dlfac
  integer :: time_days
  integer :: time_seconds
  character(len=300) :: debug_mesg
  locfield => NULL()
  locmask => NULL()
  is_stat = .false. ; if (present(is_static)) is_stat = is_static
  not_static = .not. is_stat

  ! Determine the proper array indices, noting that because of the (:,:)
  ! declaration of field, symmetric arrays are using a SW-grid indexing,
  ! but non-symmetric arrays are using a NE-grid indexing.  Send_data
  ! actually only uses the difference between ie and is to determine
  ! the output data size and assumes that halos are symmetric.
  !isv = diag_cs%is ; iev = diag_cs%ie ; jsv = diag_cs%js ; jev = diag_cs%je

  cszi = (diag_cs%ie-diag_cs%is) +1 ; dszi = (diag_cs%ied-diag_cs%isd) +1
  cszj = (diag_cs%je-diag_cs%js) +1 ; dszj = (diag_cs%jed-diag_cs%jsd) +1
  if ( size(field,1) == dszi ) then
    isv = diag_cs%is ; iev = diag_cs%ie     ! Data domain
  elseif ( size(field,1) == dszi + 1 ) then
    isv = diag_cs%is ; iev = diag_cs%ie+1   ! Symmetric data domain
  elseif ( size(field,1) == cszi) then
    isv = 1 ; iev = cszi                    ! Computational domain
  elseif ( size(field,1) == cszi + 1 ) then
    isv = 1 ; iev = cszi+1                  ! Symmetric computational domain
  else
    write (mesg,*) " peculiar size ",size(field,1)," in i-direction\n"//&
       "does not match one of ", cszi, cszi+1, dszi, dszi+1
    call MOM_error(FATAL,"post_data_3d_low: "//trim(diag%debug_str)//trim(mesg))
  endif

  if ( size(field,2) == dszj ) then
    jsv = diag_cs%js ; jev = diag_cs%je     ! Data domain
  elseif ( size(field,2) == dszj + 1 ) then
    jsv = diag_cs%js ; jev = diag_cs%je+1   ! Symmetric data domain
  elseif ( size(field,2) == cszj ) then
    jsv = 1 ; jev = cszj                    ! Computational domain
  elseif ( size(field,2) == cszj+1 ) then
    jsv = 1 ; jev = cszj+1                  ! Symmetric computational domain
  else
    write (mesg,*) " peculiar size ",size(field,2)," in j-direction\n"//&
       "does not match one of ", cszj, cszj+1, dszj, dszj+1
    call MOM_error(FATAL,"post_data_3d_low: "//trim(diag%debug_str)//trim(mesg))
  endif

  ks = lbound(field,3) ; ke = ubound(field,3)
  if ((diag%conversion_factor /= 0.) .and. (diag%conversion_factor /= 1.)) then
    allocate( locfield( lbound(field,1):ubound(field,1), lbound(field,2):ubound(field,2), ks:ke ) )
    ! locfield(:,:,:) = 0.0  ! Zeroing out this array would be a good idea, but it appears
                             ! not to be necessary.
    isv_c = isv ; jsv_c = jsv
    if (diag%fms_xyave_diag_id>0) then
      staggered_in_x = diag%axes%is_u_point .or. diag%axes%is_q_point
      staggered_in_y = diag%axes%is_v_point .or. diag%axes%is_q_point
      ! When averaging a staggered field, edge points are always required.
      if (staggered_in_x) isv_c = iev - (diag_cs%ie - diag_cs%is) - 1
      if (staggered_in_y) jsv_c = jev - (diag_cs%je - diag_cs%js) - 1
      if (isv_c < lbound(locfield,1)) call MOM_error(FATAL, &
        "It is an error to average a staggered diagnostic field that does not "//&
        "have i-direction space to represent the symmetric computational domain.")
      if (jsv_c < lbound(locfield,2)) call MOM_error(FATAL, &
        "It is an error to average a staggered diagnostic field that does not "//&
        "have j-direction space to represent the symmetric computational domain.")
    endif

    do k=ks,ke ; do j=jsv,jev ; do i=isv,iev
      locfield(i,j,k) = field(i,j,k) * diag%conversion_factor
    enddo ; enddo ; enddo
  else
    locfield => field
  endif

  if (present(mask)) then
    locmask => mask
  elseif (associated(diag%axes) .and. (not_static)) then
  ! If we were to decide to allow masking of static diagnostics, we could do so by changing the line above to
  ! elseif (associated(diag%axes) .and. (diag_CS%mask_static_diags .or. not_static)) then
    if (associated(diag%axes%mask3d)) locmask => diag%axes%mask3d
  endif

  dlfac = 1
  if (not_static .and. associated(diag%axes)) &
    dlfac = diag%axes%downsample_level_factor ! Static field downsampling is not supported yet.
  ! Downsample the diag field and mask as appropriate.
  if (dlfac > 1) then
    dl = diag%axes%downsample_level_index
    isv_o = isv ; jsv_o = jsv
    !Note that downsample_diag_field_3d takes the downsampling index
    call downsample_diag_field(locfield, locfield_dsamp, dl, diag_cs, diag, isv, iev, jsv, jev, mask)
    if ((diag%conversion_factor /= 0.) .and. (diag%conversion_factor /= 1.)) deallocate( locfield )
    locfield => locfield_dsamp
    if (present(mask)) then
      ! Replicate the downsampling of other fields to find unmasked points.
      allocate(ones, mold=locmask) ; ones(:,:,:) = 1.0
      mask_in => mask
      call downsample_field_3d(ones, locmask_dsamp, dlfac, diag%xyz_method, mask_in, diag_cs, diag, &
                               isv_o, jsv_o, isv, iev, jsv, jev)
      deallocate(ones)
      where (abs(locmask_dsamp) > 0.0) locmask_dsamp = 1.0
      locmask => locmask_dsamp
    elseif (associated(diag%axes%dsamp(dl)%mask3d)) then
      locmask => diag%axes%dsamp(dl)%mask3d
    endif
  endif
  if (associated(locmask)) call assert(size(locfield) == size(locmask), &
        'post_data_3d_low: mask size mismatch: '//trim(diag%debug_str))

  if (diag%fms_diag_id>0) then
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
        call MOM_error(FATAL, "post_data_3d_low: unknown axis type.")
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
  endif

  if (diag%fms_xyave_diag_id>0 .and. dlfac<2) then
    call post_xy_average(diag_cs, diag, locfield)
  endif

  if ((diag%conversion_factor /= 0.) .and. (diag%conversion_factor /= 1.) .and. dl<2) &
    deallocate( locfield )

end procedure post_data_3d_low
module procedure post_data_3d_by_column
  type(diag_type), pointer :: diag => null()
  integer :: buffer_slot
  diag => diag_cs%diags(diag_field_id)
  buffer_slot = diag%axes%piecemeal_3d%check_capacity_by_id(diag_field_id)
  diag%axes%piecemeal_3d%buffer(buffer_slot)%field(i,j,:) = field(:)
end procedure post_data_3d_by_column
module procedure post_data_3d_by_point
  type(diag_type), pointer :: diag => null()
  integer :: buffer_slot
  diag => diag_cs%diags(diag_field_id)
  buffer_slot = diag%axes%piecemeal_3d%check_capacity_by_id(diag_field_id)
  diag%axes%piecemeal_3d%buffer(buffer_slot)%field(i,j,k) = field
end procedure post_data_3d_by_point
module procedure post_data_3d_final
  type(diag_type), pointer :: diag => null()
  integer :: buffer_slot
  diag => diag_cs%diags(diag_field_id)
  buffer_slot = diag%axes%piecemeal_3d%find_buffer_slot(diag_field_id)
  ! Only perform an action if the buffer slot was actually used
  if (buffer_slot>0) then
    call post_data(diag_field_id, diag%axes%piecemeal_3d%buffer(buffer_slot)%field(:,:,:), diag_CS)
    call diag%axes%piecemeal_3d%mark_available(diag_field_id)
  endif
end procedure post_data_3d_final
module procedure post_product_u
  real, dimension(G%IsdB:G%IedB, G%jsd:G%jed, nz) :: u_prod ! The product of u_a and u_b [A B]
  integer :: i, j, k
  if (id <= 0) return

  do k=1,nz ; do j=G%jsc,G%jec ; do I=G%IscB,G%IecB
    u_prod(I,j,k) = u_a(I,j,k) * u_b(I,j,k)
  enddo ; enddo ; enddo
  call post_data(id, u_prod, diag, mask=mask, alt_h=alt_h)

end procedure post_product_u
module procedure post_product_sum_u
  real, dimension(G%IsdB:G%IedB, G%jsd:G%jed) :: u_sum  ! The vertical sum of the product of u_a and u_b [A B]
  integer :: i, j, k
  if (id <= 0) return

  u_sum(:,:) = 0.0
  do k=1,nz ; do j=G%jsc,G%jec ; do I=G%IscB,G%IecB
    u_sum(I,j) = u_sum(I,j) + u_a(I,j,k) * u_b(I,j,k)
  enddo ; enddo ; enddo
  call post_data(id, u_sum, diag)

end procedure post_product_sum_u
module procedure post_product_v
  real, dimension(G%isd:G%ied, G%JsdB:G%JedB, nz) :: v_prod ! The product of v_a and v_b [A B]
  integer :: i, j, k
  if (id <= 0) return

  do k=1,nz ; do J=G%JscB,G%JecB ; do i=G%isc,G%iec
    v_prod(i,J,k) = v_a(i,J,k) * v_b(i,J,k)
  enddo ; enddo ; enddo
  call post_data(id, v_prod, diag, mask=mask, alt_h=alt_h)

end procedure post_product_v
module procedure post_product_sum_v
  real, dimension(G%isd:G%ied, G%JsdB:G%JedB) :: v_sum ! The vertical sum of the product of v_a and v_b [A B]
  integer :: i, j, k
  if (id <= 0) return

  v_sum(:,:) = 0.0
  do k=1,nz ; do J=G%JscB,G%JecB ; do i=G%isc,G%iec
    v_sum(i,J) = v_sum(i,J) + v_a(i,J,k) * v_b(i,J,k)
  enddo ; enddo ; enddo
  call post_data(id, v_sum, diag)

end procedure post_product_sum_v
module procedure post_xy_average
  real, dimension(size(field,3)) :: averaged_field ! The horizontally averaged field [A ~> a]
  logical, dimension(size(field,3)) :: averaged_mask
  logical :: staggered_in_x, staggered_in_y, used
  integer :: nz, remap_nz, coord
  integer :: time_days
  integer :: time_seconds
  character(len=300) :: debug_mesg
  if (.not. diag_cs%ave_enabled) then
    return
  endif

  staggered_in_x = diag%axes%is_u_point .or. diag%axes%is_q_point
  staggered_in_y = diag%axes%is_v_point .or. diag%axes%is_q_point

  if (diag%axes%is_native) then
    call horizontally_average_diag_field(diag_cs%G, diag_cs%GV, diag_cs%h, &
                                         staggered_in_x, staggered_in_y, &
                                         diag%axes%is_layer, diag%v_extensive, &
                                         field, averaged_field, averaged_mask)
  else
    nz = size(field, 3)
    coord = diag%axes%vertical_coordinate_number
    remap_nz = diag_cs%diag_remap_cs(coord)%nz

    call assert(diag_cs%diag_remap_cs(coord)%initialized, &
                'post_xy_average: remap_cs not initialized.')

    call assert(IMPLIES(diag%axes%is_layer, nz == remap_nz), &
              'post_xy_average: layer field dimension mismatch.')
    call assert(IMPLIES(.not. diag%axes%is_layer, nz == remap_nz+1), &
              'post_xy_average: interface field dimension mismatch.')

    call horizontally_average_diag_field(diag_cs%G, diag_cs%GV, &
                                         diag_cs%diag_remap_cs(coord)%h, &
                                         staggered_in_x, staggered_in_y, &
                                         diag%axes%is_layer, diag%v_extensive, &
                                         field, averaged_field, averaged_mask)
  endif

  if (diag_cs%diag_as_chksum) then
    ! Append timestep to mesg
    call get_time(diag_cs%time_end, time_seconds, days=time_days)
    write(debug_mesg, '(a, 1x, i0, 1x, i0)') &
        trim(diag%debug_str)//'_xyave', time_days, time_seconds

    call zchksum(averaged_field, debug_mesg, logunit=diag_CS%chksum_iounit)
  else
    used = send_data_infra(diag%fms_xyave_diag_id, averaged_field, &
                           time=diag_cs%time_end, weight=diag_cs%time_int, mask=averaged_mask)
  endif
end procedure post_xy_average
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
module procedure get_diag_time_end
  get_diag_time_end = diag_cs%time_end
end procedure get_diag_time_end
module procedure register_diag_field
  real :: MOM_missing_value ! A value used to indicate missing values in output files, in arbitrary units [a]
  type(diag_ctrl), pointer :: diag_cs => NULL() ! A structure that is used to regulate diagnostic output
  type(axes_grp), pointer :: remap_axes
  type(axes_grp), pointer :: axes
  type(axes_grp), pointer :: axes_d2
  integer :: dm_id, i, dl
  character(len=256) :: msg, cm_string
  character(len=256) :: new_module_name
  character(len=480) :: module_list, var_list
  character(len=24)  :: dimensions
  integer :: num_modnm, num_varnm
  logical :: active
  character(len=2)   :: dl_str
  diag_cs => axes_in%diag_cs

  ! Check if the axes match a standard grid axis.
  ! If not, allocate the new axis and copy the contents.
  if (axes_in%id == diag_cs%axesTL%id) then
    axes => diag_cs%axesTL
  elseif (axes_in%id == diag_cs%axesBL%id) then
    axes => diag_cs%axesBL
  elseif (axes_in%id == diag_cs%axesCuL%id) then
    axes => diag_cs%axesCuL
  elseif (axes_in%id == diag_cs%axesCvL%id) then
    axes => diag_cs%axesCvL
  elseif (axes_in%id == diag_cs%axesTi%id) then
    axes => diag_cs%axesTi
  elseif (axes_in%id == diag_cs%axesBi%id) then
    axes => diag_cs%axesBi
  elseif (axes_in%id == diag_cs%axesCui%id) then
    axes => diag_cs%axesCui
  elseif (axes_in%id == diag_cs%axesCvi%id) then
    axes => diag_cs%axesCvi
  elseif (axes_in%id == diag_cs%axesT1%id) then
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
             cell_methods=cell_methods, x_cell_method=x_cell_method, &
             y_cell_method=y_cell_method, v_cell_method=v_cell_method, &
             conversion=conversion, v_extensive=v_extensive)
  if (associated(axes%xyave_axes)) then
    num_varnm = 2 ; var_list = "{"//trim(field_name)//","//trim(field_name)//"_xyave"
  else
    num_varnm = 1 ; var_list = "{"//trim(field_name)
  endif
  if (present(cmor_field_name)) then
    if (associated(axes%xyave_axes)) then
      num_varnm = num_varnm + 2
      var_list = trim(var_list)//","//trim(cmor_field_name)//","//trim(cmor_field_name)//"_xyave"
    else
      num_varnm = num_varnm + 1
      var_list = trim(var_list)//","//trim(cmor_field_name)
    endif
  endif
  var_list = trim(var_list)//"}"

  ! For each diagnostic coordinate register the diagnostic again under a different module name
  do i=1,diag_cs%num_diag_coords
    new_module_name = trim(module_name)//'_'//trim(diag_cs%diag_remap_cs(i)%diag_module_suffix)

    ! Register diagnostics remapped to z vertical coordinate
    if (axes_in%rank == 3) then
      remap_axes => null()
      if ((axes_in%id == diag_cs%axesTL%id)) then
        remap_axes => diag_cs%remap_axesTL(i)
      elseif (axes_in%id == diag_cs%axesBL%id) then
        remap_axes => diag_cs%remap_axesBL(i)
      elseif (axes_in%id == diag_cs%axesCuL%id ) then
        remap_axes => diag_cs%remap_axesCuL(i)
      elseif (axes_in%id == diag_cs%axesCvL%id) then
        remap_axes => diag_cs%remap_axesCvL(i)
      elseif (axes_in%id == diag_cs%axesTi%id) then
        remap_axes => diag_cs%remap_axesTi(i)
      elseif (axes_in%id == diag_cs%axesBi%id) then
        remap_axes => diag_cs%remap_axesBi(i)
      elseif (axes_in%id == diag_cs%axesCui%id ) then
        remap_axes => diag_cs%remap_axesCui(i)
      elseif (axes_in%id == diag_cs%axesCvi%id) then
        remap_axes => diag_cs%remap_axesCvi(i)
      endif
      ! When the MOM_diag_to_Z module has been obsoleted we can assume remap_axes will
      ! always exist but in the mean-time we have to do this check:
      ! call assert(associated(remap_axes), 'register_diag_field: remap_axes not set')
      if (associated(remap_axes)) then
        if (remap_axes%needs_remapping .or. remap_axes%needs_interpolating) then
          active = register_diag_field_expand_cmor(dm_id, new_module_name, field_name, remap_axes, &
                     init_time, long_name=long_name, units=units, missing_value=MOM_missing_value, &
                     range=range, mask_variant=mask_variant, standard_name=standard_name, &
                     verbose=verbose, do_not_log=do_not_log, err_msg=err_msg, &
                     interp_method=interp_method, tile_count=tile_count, &
                     cmor_field_name=cmor_field_name, cmor_long_name=cmor_long_name, &
                     cmor_units=cmor_units, cmor_standard_name=cmor_standard_name, &
                     cell_methods=cell_methods, x_cell_method=x_cell_method, &
                     y_cell_method=y_cell_method, v_cell_method=v_cell_method, &
                     conversion=conversion, v_extensive=v_extensive)
          if (active) then
            call diag_remap_set_active(diag_cs%diag_remap_cs(i))
          endif
          module_list = trim(module_list)//","//trim(new_module_name)
          num_modnm = num_modnm + 1
        endif ! remap_axes%needs_remapping
      endif ! associated(remap_axes)
    endif ! axes%rank == 3
  enddo ! i

  ! Register downsampled diagnostics
  do dl=1, diag_cs%num_diag_dsamp_levels
    ! Do not attempt to checksum the downsampled diagnostics
    if (diag_cs%diag_as_chksum) cycle

    write(dl_str, '(i0)') diag_cs%diag_dsamp_levels(dl)
    new_module_name = trim(module_name)//'_d'//trim(dl_str)

    axes_d2 => null()
    if (axes_in%rank == 3 .or. axes_in%rank == 2 ) then
      if (axes_in%id == diag_cs%axesTL%id) then
        axes_d2 => diag_cs%dsamp(dl)%axesTL
      elseif (axes_in%id == diag_cs%axesBL%id) then
        axes_d2 => diag_cs%dsamp(dl)%axesBL
      elseif (axes_in%id == diag_cs%axesCuL%id ) then
        axes_d2 => diag_cs%dsamp(dl)%axesCuL
      elseif (axes_in%id == diag_cs%axesCvL%id) then
        axes_d2 => diag_cs%dsamp(dl)%axesCvL
      elseif (axes_in%id == diag_cs%axesTi%id) then
        axes_d2 => diag_cs%dsamp(dl)%axesTi
      elseif (axes_in%id == diag_cs%axesBi%id) then
        axes_d2 => diag_cs%dsamp(dl)%axesBi
      elseif (axes_in%id == diag_cs%axesCui%id ) then
        axes_d2 => diag_cs%dsamp(dl)%axesCui
      elseif (axes_in%id == diag_cs%axesCvi%id) then
        axes_d2 => diag_cs%dsamp(dl)%axesCvi
      elseif (axes_in%id == diag_cs%axesT1%id) then
        axes_d2 => diag_cs%dsamp(dl)%axesT1
      elseif (axes_in%id == diag_cs%axesB1%id) then
        axes_d2 => diag_cs%dsamp(dl)%axesB1
      elseif (axes_in%id == diag_cs%axesCu1%id ) then
        axes_d2 => diag_cs%dsamp(dl)%axesCu1
      elseif (axes_in%id == diag_cs%axesCv1%id) then
        axes_d2 => diag_cs%dsamp(dl)%axesCv1
      else
        !Niki: Should we worry about these, e.g., diag_to_Z_CS?
        call MOM_error(WARNING,"register_diag_field: Could not find a proper axes for " &
              //trim(new_module_name)//"-"//trim(field_name))
      endif
    endif

    ! Register the native diagnostic
    if (associated(axes_d2)) then
      active = register_diag_field_expand_cmor(dm_id, new_module_name, field_name, axes_d2, &
                init_time, long_name=long_name, units=units, missing_value=MOM_missing_value, &
                range=range, mask_variant=mask_variant, standard_name=standard_name, &
                verbose=verbose, do_not_log=do_not_log, err_msg=err_msg, &
                interp_method=interp_method, tile_count=tile_count, &
                cmor_field_name=cmor_field_name, cmor_long_name=cmor_long_name, &
                cmor_units=cmor_units, cmor_standard_name=cmor_standard_name, &
                cell_methods=cell_methods, x_cell_method=x_cell_method, &
                y_cell_method=y_cell_method, v_cell_method=v_cell_method, &
                conversion=conversion, v_extensive=v_extensive)
      module_list = trim(module_list)//","//trim(new_module_name)
      num_modnm = num_modnm + 1
    endif

    ! For each diagnostic coordinate register the diagnostic again under a different module name
    do i=1,diag_cs%num_diag_coords
      new_module_name = trim(module_name)//'_'//trim(diag_cs%diag_remap_cs(i)%diag_module_suffix)//'_d'//trim(dl_str)

      ! Register diagnostics remapped to z vertical coordinate
      if (axes_in%rank == 3) then
        remap_axes => null()
        if ((axes_in%id == diag_cs%axesTL%id)) then
          remap_axes => diag_cs%dsamp(dl)%remap_axesTL(i)
        elseif (axes_in%id == diag_cs%axesBL%id) then
          remap_axes => diag_cs%dsamp(dl)%remap_axesBL(i)
        elseif (axes_in%id == diag_cs%axesCuL%id ) then
          remap_axes => diag_cs%dsamp(dl)%remap_axesCuL(i)
        elseif (axes_in%id == diag_cs%axesCvL%id) then
          remap_axes => diag_cs%dsamp(dl)%remap_axesCvL(i)
        elseif (axes_in%id == diag_cs%axesTi%id) then
          remap_axes => diag_cs%dsamp(dl)%remap_axesTi(i)
        elseif (axes_in%id == diag_cs%axesBi%id) then
          remap_axes => diag_cs%dsamp(dl)%remap_axesBi(i)
        elseif (axes_in%id == diag_cs%axesCui%id ) then
          remap_axes => diag_cs%dsamp(dl)%remap_axesCui(i)
        elseif (axes_in%id == diag_cs%axesCvi%id) then
          remap_axes => diag_cs%dsamp(dl)%remap_axesCvi(i)
        endif

        ! When the MOM_diag_to_Z module has been obsoleted we can assume remap_axes will
        ! always exist but in the mean-time we have to do this check:
        ! call assert(associated(remap_axes), 'register_diag_field: remap_axes not set')
        if (associated(remap_axes)) then
          if (remap_axes%needs_remapping .or. remap_axes%needs_interpolating) then
            active = register_diag_field_expand_cmor(dm_id, new_module_name, field_name, remap_axes, &
                    init_time, long_name=long_name, units=units, missing_value=MOM_missing_value, &
                    range=range, mask_variant=mask_variant, standard_name=standard_name, &
                    verbose=verbose, do_not_log=do_not_log, err_msg=err_msg, &
                    interp_method=interp_method, tile_count=tile_count, &
                    cmor_field_name=cmor_field_name, cmor_long_name=cmor_long_name, &
                    cmor_units=cmor_units, cmor_standard_name=cmor_standard_name, &
                    cell_methods=cell_methods, x_cell_method=x_cell_method, &
                    y_cell_method=y_cell_method, v_cell_method=v_cell_method, &
                    conversion=conversion, v_extensive=v_extensive)
            if (active) then
              call diag_remap_set_active(diag_cs%diag_remap_cs(i))
            endif
            module_list = trim(module_list)//","//trim(new_module_name)
            num_modnm = num_modnm + 1
          endif ! remap_axes%needs_remapping
        endif ! associated(remap_axes)
      endif ! axes%rank == 3
    enddo ! i
  enddo ! dl

  dimensions = ""
  if (axes_in%is_h_point)   dimensions = trim(dimensions)//" xh, yh,"
  if (axes_in%is_q_point)   dimensions = trim(dimensions)//" xq, yq,"
  if (axes_in%is_u_point)   dimensions = trim(dimensions)//" xq, yh,"
  if (axes_in%is_v_point)   dimensions = trim(dimensions)//" xh, yq,"
  if (axes_in%is_layer)     dimensions = trim(dimensions)//" zl,"
  if (axes_in%is_interface) dimensions = trim(dimensions)//" zi,"
  if (len_trim(dimensions) > 0) dimensions = trim_trailing_commas(dimensions)

  if (is_root_pe() .and. (diag_CS%available_diag_doc_unit > 0)) then
    msg = ''
    if (present(cmor_field_name)) msg = 'CMOR equivalent is "'//trim(cmor_field_name)//'"'
    call attach_cell_methods(-1, axes, cm_string, cell_methods, &
                             x_cell_method, y_cell_method, v_cell_method, &
                             v_extensive=v_extensive)
    module_list = trim(module_list)//"}"
    if (num_modnm <= 1) module_list = module_name
    if (num_varnm <= 1) var_list = ''

    call log_available_diag(dm_id>0, module_list, field_name, cm_string, msg, diag_CS, &
                            long_name, units, standard_name, variants=var_list, dimensions=dimensions)
  endif

  register_diag_field = dm_id

end procedure register_diag_field
module procedure register_diag_field_expand_cmor
  real :: MOM_missing_value ! A value used to indicate missing values in output files, in arbitrary units [a]
  type(diag_ctrl), pointer :: diag_cs => null()
  type(diag_type), pointer :: this_diag => null()
  integer :: fms_id, fms_xyave_id
  character(len=256) :: posted_cmor_units, posted_cmor_standard_name, posted_cmor_long_name, cm_string
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
    call attach_cell_methods(fms_id, axes, cm_string, cell_methods, &
                             x_cell_method, y_cell_method, v_cell_method, &
                             v_extensive=v_extensive)
  ! Associated horizontally area-averaged diagnostic
  fms_xyave_id = DIAG_FIELD_NOT_FOUND
  if (associated(axes%xyave_axes)) then
    fms_xyave_id = register_diag_field_expand_axes(module_name, trim(field_name)//'_xyave', &
             axes%xyave_axes, init_time, &
             long_name=long_name, units=units, missing_value=MOM_missing_value, &
             range=range, mask_variant=mask_variant, standard_name=standard_name, &
             verbose=verbose, do_not_log=do_not_log, err_msg=err_msg, &
             interp_method=interp_method, tile_count=tile_count)
    if (.not. diag_cs%diag_as_chksum) &
      call attach_cell_methods(fms_xyave_id, axes%xyave_axes, cm_string, &
                               cell_methods, v_cell_method, v_extensive=v_extensive)
  endif
  this_diag => null()
  if (fms_id /= DIAG_FIELD_NOT_FOUND .or. fms_xyave_id /= DIAG_FIELD_NOT_FOUND) then
    call add_diag_to_list(diag_cs, dm_id, fms_id, this_diag, axes, module_name, field_name)
    this_diag%fms_xyave_diag_id = fms_xyave_id
    ! Encode and save the cell methods for this diagnostic
    this_diag%xyz_method = xyz_method(axes, x_cell_method, y_cell_method, v_cell_method, v_extensive)
    if (present(v_extensive)) this_diag%v_extensive = v_extensive
    if (present(conversion)) this_diag%conversion_factor = conversion
    register_diag_field_expand_cmor = .true.
  endif

  ! For the CMOR variation of the above diagnostic
  if (present(cmor_field_name) .and. .not. diag_cs%diag_as_chksum) then
    ! Fallback values for strings set to "NULL"
    posted_cmor_units = "not provided"         !
    posted_cmor_standard_name = "not provided" ! Values might be able to be replaced with a CS%missing field?
    posted_cmor_long_name = "not provided"     !

    ! If attributes are present for MOM variable names, use them first for the register_diag_field
    ! call for CMOR version of the variable
    if (present(units)) posted_cmor_units = units
    if (present(standard_name)) posted_cmor_standard_name = standard_name
    if (present(long_name)) posted_cmor_long_name = long_name

    ! If specified in the call to register_diag_field, override attributes with the CMOR versions
    if (present(cmor_units)) posted_cmor_units = cmor_units
    if (present(cmor_standard_name)) posted_cmor_standard_name = cmor_standard_name
    if (present(cmor_long_name)) posted_cmor_long_name = cmor_long_name

    fms_id = register_diag_field_expand_axes(module_name, cmor_field_name, axes, init_time,    &
               long_name=trim(posted_cmor_long_name), units=trim(posted_cmor_units),                  &
               missing_value=MOM_missing_value, range=range, mask_variant=mask_variant,               &
               standard_name=trim(posted_cmor_standard_name), verbose=verbose, do_not_log=do_not_log, &
               err_msg=err_msg, interp_method=interp_method, tile_count=tile_count)
    call attach_cell_methods(fms_id, axes, cm_string, &
                             cell_methods, x_cell_method, y_cell_method, v_cell_method, &
                             v_extensive=v_extensive)
    ! Associated horizontally area-averaged diagnostic
    fms_xyave_id = DIAG_FIELD_NOT_FOUND
    if (associated(axes%xyave_axes)) then
      fms_xyave_id = register_diag_field_expand_axes(module_name, trim(cmor_field_name)//'_xyave', &
               axes%xyave_axes, init_time, &
               long_name=trim(posted_cmor_long_name), units=trim(posted_cmor_units),                  &
               missing_value=MOM_missing_value, range=range, mask_variant=mask_variant,               &
               standard_name=trim(posted_cmor_standard_name), verbose=verbose, do_not_log=do_not_log, &
               err_msg=err_msg, interp_method=interp_method, tile_count=tile_count)
      call attach_cell_methods(fms_xyave_id, axes%xyave_axes, cm_string, &
                               cell_methods, v_cell_method, v_extensive=v_extensive)
    endif
    this_diag => null()
    if (fms_id /= DIAG_FIELD_NOT_FOUND .or. fms_xyave_id /= DIAG_FIELD_NOT_FOUND) then
      call add_diag_to_list(diag_cs, dm_id, fms_id, this_diag, axes, module_name, field_name)
      this_diag%fms_xyave_diag_id = fms_xyave_id
      ! Encode and save the cell methods for this diagnostic
      this_diag%xyz_method = xyz_method(axes, x_cell_method, y_cell_method, v_cell_method, v_extensive)
      if (present(v_extensive)) this_diag%v_extensive = v_extensive
      if (present(conversion)) this_diag%conversion_factor = conversion
      register_diag_field_expand_cmor = .true.
    endif
  endif

end procedure register_diag_field_expand_cmor
module procedure register_diag_field_expand_axes
  integer :: fms_id, area_id, volume_id
  area_id = axes%id_area
  volume_id = axes%id_volume

  ! Get the FMS diagnostic id
  if (axes%diag_cs%diag_as_chksum) then
    fms_id = axes%diag_cs%num_chksum_diags + 1
    axes%diag_cs%num_chksum_diags = fms_id
  elseif (present(interp_method) .or. axes%is_h_point) then
    ! If interp_method is provided we must use it
    if (area_id>0) then
      if (volume_id>0) then
        fms_id = register_diag_field_infra(module_name, field_name, axes%handles, &
                 init_time, long_name=long_name, units=units, missing_value=missing_value, &
                 range=range, mask_variant=mask_variant, standard_name=standard_name, &
                 verbose=verbose, do_not_log=do_not_log, err_msg=err_msg, &
                 interp_method=interp_method, tile_count=tile_count, area=area_id, volume=volume_id)
      else
        fms_id = register_diag_field_infra(module_name, field_name, axes%handles, &
                 init_time, long_name=long_name, units=units, missing_value=missing_value, &
                 range=range, mask_variant=mask_variant, standard_name=standard_name, &
                 verbose=verbose, do_not_log=do_not_log, err_msg=err_msg, &
                 interp_method=interp_method, tile_count=tile_count, area=area_id)
      endif
    else
      if (volume_id>0) then
        fms_id = register_diag_field_infra(module_name, field_name, axes%handles, &
                 init_time, long_name=long_name, units=units, missing_value=missing_value, &
                 range=range, mask_variant=mask_variant, standard_name=standard_name, &
                 verbose=verbose, do_not_log=do_not_log, err_msg=err_msg, &
                 interp_method=interp_method, tile_count=tile_count, volume=volume_id)
      else
        fms_id = register_diag_field_infra(module_name, field_name, axes%handles, &
                 init_time, long_name=long_name, units=units, missing_value=missing_value, &
                 range=range, mask_variant=mask_variant, standard_name=standard_name, &
                 verbose=verbose, do_not_log=do_not_log, err_msg=err_msg, &
                 interp_method=interp_method, tile_count=tile_count)
      endif
    endif
  else
    ! If interp_method is not provided and the field is not at an h-point then interp_method='none'
    if (area_id>0) then
      if (volume_id>0) then
        fms_id = register_diag_field_infra(module_name, field_name, axes%handles, &
                 init_time, long_name=long_name, units=units, missing_value=missing_value, &
                 range=range, mask_variant=mask_variant, standard_name=standard_name, &
                 verbose=verbose, do_not_log=do_not_log, err_msg=err_msg, &
                 interp_method='none', tile_count=tile_count, area=area_id, volume=volume_id)
      else
        fms_id = register_diag_field_infra(module_name, field_name, axes%handles, &
                 init_time, long_name=long_name, units=units, missing_value=missing_value, &
                 range=range, mask_variant=mask_variant, standard_name=standard_name, &
                 verbose=verbose, do_not_log=do_not_log, err_msg=err_msg, &
                 interp_method='none', tile_count=tile_count, area=area_id)
      endif
    else
      if (volume_id>0) then
        fms_id = register_diag_field_infra(module_name, field_name, axes%handles, &
                 init_time, long_name=long_name, units=units, missing_value=missing_value, &
                 range=range, mask_variant=mask_variant, standard_name=standard_name, &
                 verbose=verbose, do_not_log=do_not_log, err_msg=err_msg, &
                 interp_method='none', tile_count=tile_count, volume=volume_id)
      else
        fms_id = register_diag_field_infra(module_name, field_name, axes%handles, &
                 init_time, long_name=long_name, units=units, missing_value=missing_value, &
                 range=range, mask_variant=mask_variant, standard_name=standard_name, &
                 verbose=verbose, do_not_log=do_not_log, err_msg=err_msg, &
                 interp_method='none', tile_count=tile_count)
      endif
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
module procedure xyz_method
  character(len=9) :: mstr
  xyz_method = 111

  mstr = axes%v_cell_method
  if (present(v_extensive)) then
    if (present(v_cell_method)) call MOM_error(FATAL, "xyz_method: " // &
       'Vertical cell method was specified along with the vertically extensive flag.')
    if (v_extensive) then
      mstr='sum'
    else
      mstr='mean'
    endif
  elseif (present(v_cell_method)) then
    mstr = v_cell_method
  endif
  if (trim(mstr)=='sum') then
    xyz_method = xyz_method + 1
  elseif (trim(mstr)=='mean') then
    xyz_method = xyz_method + 2
  endif

  mstr = axes%y_cell_method
  if (present(y_cell_method)) mstr = y_cell_method
  if (trim(mstr)=='sum') then
    xyz_method = xyz_method + 10
  elseif (trim(mstr)=='mean') then
    xyz_method = xyz_method + 20
  endif

  mstr = axes%x_cell_method
  if (present(x_cell_method)) mstr = x_cell_method
  if (trim(mstr)=='sum') then
    xyz_method = xyz_method + 100
  elseif (trim(mstr)=='mean') then
    xyz_method = xyz_method + 200
  endif

end procedure xyz_method
module procedure attach_cell_methods
  character(len=9) :: axis_name
  logical :: x_mean, y_mean, x_sum, y_sum
  x_mean = .false.
  y_mean = .false.
  x_sum = .false.
  y_sum = .false.

  ostring = ''
  if (present(cell_methods)) then
    if (present(x_cell_method) .or. present(y_cell_method) .or. present(v_cell_method) &
        .or. present(v_extensive)) then
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
    if (present(v_cell_method)) then
      if (present(v_extensive)) call MOM_error(FATAL, "attach_cell_methods: " // &
           'Vertical cell method was specified along with the vertically extensive flag.')
      if (len(trim(v_cell_method))>0) then
        if (axes%rank==1) then
          call get_MOM_diag_axis_name(axes%handles(1), axis_name)
        elseif (axes%rank==3) then
          call get_MOM_diag_axis_name(axes%handles(3), axis_name)
        endif
        call MOM_diag_field_add_attribute(id, 'cell_methods', trim(axis_name)//':'//trim(v_cell_method))
        ostring = trim(adjustl(ostring))//' '//trim(axis_name)//':'//trim(v_cell_method)
      endif
    elseif (present(v_extensive)) then
      if (v_extensive) then
        if (axes%rank==1) then
          call get_MOM_diag_axis_name(axes%handles(1), axis_name)
        elseif (axes%rank==3) then
          call get_MOM_diag_axis_name(axes%handles(3), axis_name)
        endif
        call MOM_diag_field_add_attribute(id, 'cell_methods', trim(axis_name)//':sum')
        ostring = trim(adjustl(ostring))//' '//trim(axis_name)//':sum'
      endif
    else
      if (len(trim(axes%v_cell_method))>0) then
        if (axes%rank==1) then
          call get_MOM_diag_axis_name(axes%handles(1), axis_name)
        elseif (axes%rank==3) then
          call get_MOM_diag_axis_name(axes%handles(3), axis_name)
        endif
        call MOM_diag_field_add_attribute(id, 'cell_methods', trim(axis_name)//':'//trim(axes%v_cell_method))
        ostring = trim(adjustl(ostring))//' '//trim(axis_name)//':'//trim(axes%v_cell_method)
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
    ! register_diag_field_infra call for CMOR version of the variable
    if (present(units)) posted_cmor_units = units
    if (present(standard_name)) posted_cmor_standard_name = standard_name
    if (present(long_name)) posted_cmor_long_name = long_name

    ! If specified in the call to register_scalar_field, override attributes with the CMOR versions
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
module procedure register_static_field
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
    ! call for CMOR version of the variable
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
  if (axes%is_layer)     dimensions = trim(dimensions)//" zl,"
  if (axes%is_interface) dimensions = trim(dimensions)//" zi,"
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

end procedure register_static_field
module procedure describe_option
  character(len=480) :: mesg
  integer :: len_ind
  len_ind = len_trim(value)  ! Add error handling for long values?

  mesg = "    ! "//trim(opt_name)//": "//trim(value)
  write(diag_CS%available_diag_doc_unit, '(a)') trim(mesg)
end procedure describe_option
module procedure ocean_register_diag
  character(len=64) :: var_name         ! A variable's name.
  character(len=48) :: units            ! A variable's units.
  character(len=240) :: longname        ! A variable's longname.
  character(len=8) :: hor_grid, z_grid  ! Variable grid info.
  real :: conversion ! A multiplicative factor for unit conversions for output,
  type(axes_grp), pointer :: axes => NULL()
  call query_vardesc(var_desc, units=units, longname=longname, hor_grid=hor_grid, &
                     z_grid=z_grid, conversion=conversion, caller="ocean_register_diag")

  ! Use the hor_grid and z_grid components of vardesc to determine the
  ! desired axes to register the diagnostic field for.
  select case (z_grid)

    case ("L")
      select case (hor_grid)
        case ("q")  ; axes => diag_cs%axesBL
        case ("h")  ; axes => diag_cs%axesTL
        case ("u")  ; axes => diag_cs%axesCuL
        case ("v")  ; axes => diag_cs%axesCvL
        case ("Bu") ; axes => diag_cs%axesBL
        case ("T")  ; axes => diag_cs%axesTL
        case ("Cu") ; axes => diag_cs%axesCuL
        case ("Cv") ; axes => diag_cs%axesCvL
        case ("z")  ; axes => diag_cs%axeszL
        case default ; call MOM_error(FATAL, "ocean_register_diag: " // &
                                      "unknown hor_grid component "//trim(hor_grid))
      end select

    case ("i")
      select case (hor_grid)
        case ("q")  ; axes => diag_cs%axesBi
        case ("h")  ; axes => diag_cs%axesTi
        case ("u")  ; axes => diag_cs%axesCui
        case ("v")  ; axes => diag_cs%axesCvi
        case ("Bu") ; axes => diag_cs%axesBi
        case ("T")  ; axes => diag_cs%axesTi
        case ("Cu") ; axes => diag_cs%axesCui
        case ("Cv") ; axes => diag_cs%axesCvi
        case ("z")  ; axes => diag_cs%axeszi
        case default ; call MOM_error(FATAL, "ocean_register_diag: " // &
                                      "unknown hor_grid component "//trim(hor_grid))
      end select

    case ("1")
      select case (hor_grid)
        case ("q")  ; axes => diag_cs%axesB1
        case ("h")  ; axes => diag_cs%axesT1
        case ("u")  ; axes => diag_cs%axesCu1
        case ("v")  ; axes => diag_cs%axesCv1
        case ("Bu") ; axes => diag_cs%axesB1
        case ("T")  ; axes => diag_cs%axesT1
        case ("Cu") ; axes => diag_cs%axesCu1
        case ("Cv") ; axes => diag_cs%axesCv1
        case default ; call MOM_error(FATAL, "ocean_register_diag: " // &
                                      "unknown hor_grid component "//trim(hor_grid))
      end select

    case default
      call MOM_error(FATAL,&
        "ocean_register_diag: unknown z_grid component "//trim(z_grid))
  end select

  ocean_register_diag = register_diag_field("ocean_model", trim(var_name), axes, day, &
          trim(longname), units=trim(units), conversion=conversion, missing_value=-1.0e+34)

end procedure ocean_register_diag
module procedure diag_mediator_infrastructure_init
  call MOM_diag_manager_init(err_msg=err_msg)
end procedure diag_mediator_infrastructure_init
module procedure diag_mediator_init
  integer :: ios, i, new_unit, dl, dlfac
  logical :: opened, new_file
  integer :: remap_answer_date    ! The vintage of the order of arithmetic and expressions to use
  integer :: default_answer_date  ! The default setting for the various ANSWER_DATE flags.
  logical :: om4_remap_via_sub_cells ! Use the OM4-era ramap_via_sub_cells for diagnostics
  logical :: dz_diag_needed       ! Logical set True if we need to store dz_begin for reintegrating
  character(len=8)   :: this_pe
  character(len=240) :: doc_file, doc_file_dflt, doc_path
  character(len=240), allocatable :: diag_coords(:)
# include "version_variable.h"
  character(len=40) :: mdl = "MOM_diag_mediator" ! This module's name.
  character(len=32) :: filename_appendix = '' ! FMS appendix to filename for ensemble runs
  character(len=16) :: dsamp_domain_name
  id_clock_diag_mediator = cpu_clock_id('(Ocean diagnostics framework)', grain=CLOCK_MODULE)
  id_clock_diag_remap = cpu_clock_id('(Ocean diagnostics remapping)', grain=CLOCK_ROUTINE)
  id_clock_diag_grid_updates = cpu_clock_id('(Ocean diagnostics grid updates)', grain=CLOCK_ROUTINE)

  ! Allocate and initialize list of all diagnostics (and variants)
  allocate(diag_cs%diags(DIAG_ALLOC_CHUNK_SIZE))
  diag_cs%next_free_diag_id = 1
  do i=1, DIAG_ALLOC_CHUNK_SIZE
    call initialize_diag_type(diag_cs%diags(i))
  enddo

  diag_cs%show_call_tree = callTree_showQuery()

  ! Read all relevant parameters and write them to the model log.
  call log_version(param_file, mdl, version, "")

  call get_param(param_file, mdl, 'NUM_DIAG_COORDS', diag_cs%num_diag_coords, &
                 'The number of diagnostic vertical coordinates to use. '//&
                 'For each coordinate, an entry in DIAG_COORDS must be provided.', &
                 default=1)
  call get_param(param_file, mdl, "DEFAULT_ANSWER_DATE", default_answer_date, &
                 "This sets the default value for the various _ANSWER_DATE parameters.", &
                 default=99991231)
  call get_param(param_file, mdl, "REMAPPING_USE_OM4_SUBCELLS", om4_remap_via_sub_cells, &
                 do_not_log=.true., default=.true.)
  call get_param(param_file, mdl, "DIAG_REMAPPING_USE_OM4_SUBCELLS", om4_remap_via_sub_cells, &
                 "If true, use the OM4 remapping-via-subcells algorithm for diagnostics. "//&
                 "See REMAPPING_USE_OM4_SUBCELLS for details. "//&
                 "We recommend setting this option to false.", default=om4_remap_via_sub_cells)
  call get_param(param_file, mdl, "REMAPPING_ANSWER_DATE", remap_answer_date, &
                 "The vintage of the expressions and order of arithmetic to use for remapping.  "//&
                 "Values below 20190101 result in the use of older, less accurate expressions "//&
                 "that were in use at the end of 2018.  Higher values result in the use of more "//&
                 "robust and accurate forms of mathematically equivalent expressions.", &
                 default=default_answer_date, do_not_log=.not.GV%Boussinesq)
  if (.not.GV%Boussinesq) remap_answer_date = max(remap_answer_date, 20230701)
  call get_param(param_file, mdl, 'USE_INDEX_DIAGNOSTIC_AXES', diag_cs%index_space_axes, &
                 'If true, use a grid index coordinate convention for diagnostic axes. ', &
                 default=.false.)
  call get_param(param_file, mdl, 'SYMMETRIC_DOWNSAMPLE_SUMS', diag_cs%symmetric_downsample_sums, &
                 'If true, use rotationally symmetric sums when downsampling diagnostics.', &
                 default=.false.)


  dz_diag_needed = .false.
  if (diag_cs%num_diag_coords>0) then
    allocate(diag_coords(diag_cs%num_diag_coords))
    if (diag_cs%num_diag_coords==1) then ! The default is to provide just one instance of Z*
      call get_param(param_file, mdl, 'DIAG_COORDS', diag_coords, &
                 'A list of string tuples associating diag_table modules to '//&
                 'a coordinate definition used for diagnostics. Each string '//&
                 'is of the form "MODULE_SUFFIX PARAMETER_SUFFIX COORDINATE_NAME".', &
                 default='z Z ZSTAR')
    else ! If using more than 1 diagnostic coordinate, all must be explicitly defined
      call get_param(param_file, mdl, 'DIAG_COORDS', diag_coords, &
                 'A list of string tuples associating diag_table modules to '//&
                 'a coordinate definition used for diagnostics. Each string '//&
                 'is of the form "MODULE_SUFFIX,PARAMETER_SUFFIX,COORDINATE_NAME".', &
                 fail_if_missing=.true.)
    endif
    allocate(diag_cs%diag_remap_cs(diag_cs%num_diag_coords))
    ! Initialize each diagnostic vertical coordinate
    do i=1, diag_cs%num_diag_coords
      call diag_remap_init(diag_cs%diag_remap_cs(i), diag_coords(i), om4_remap_via_sub_cells, remap_answer_date, GV)
      if (diag_cs%diag_remap_cs(i)%Z_based_coord) dz_diag_needed = .true.
    enddo
    deallocate(diag_coords)
  endif

  call get_param(param_file, mdl, 'DIAG_MISVAL', diag_cs%missing_value, &
                 'Set the default missing value to use for diagnostics.', &
                 units="various", default=1.e20)
  call get_param(param_file, mdl, 'DIAG_AS_CHKSUM', diag_cs%diag_as_chksum, &
                 'Instead of writing diagnostics to the diag manager, write '//&
                 'a text file containing the checksum (bitcount) of the array.',  &
                 default=.false.)

  if (diag_cs%diag_as_chksum) &
    diag_cs%num_chksum_diags = 0

  ! Keep pointers to the grid, h, T, S needed for diagnostic remapping
  diag_cs%G => G
  diag_cs%GV => GV
  diag_cs%US => US
  diag_cs%h => null()
  diag_cs%T => null()
  diag_cs%S => null()
  diag_cs%eqn_of_state => null()
  diag_cs%tv => null()

  allocate(diag_cs%h_begin(G%isd:G%ied,G%jsd:G%jed,nz))
  if (dz_diag_needed) allocate(diag_cs%dz_begin(G%isd:G%ied,G%jsd:G%jed,nz))
#if defined(DEBUG) || defined(__DO_SAFETY_CHECKS__)
  allocate(diag_cs%h_old(G%isd:G%ied,G%jsd:G%jed,nz), source=0.0)
#endif
  allocate(diag_cs%OBC_u(G%IsdB:G%IedB,G%jsd:G%jed), source=0)
  allocate(diag_cs%OBC_v(G%isd:G%ied,G%JsdB:G%JedB), source=0)

  diag_cs%is = G%isc - (G%isd-1) ; diag_cs%ie = G%iec - (G%isd-1)
  diag_cs%js = G%jsc - (G%jsd-1) ; diag_cs%je = G%jec - (G%jsd-1)
  diag_cs%isd = G%isd ; diag_cs%ied = G%ied
  diag_cs%jsd = G%jsd ; diag_cs%jed = G%jed

  ! In this code design
  ! diag_cs%num_diag_dsamp_levels is the number of downsampling levels requested in the parameters
  ! diag_cs%diag_dsamp_levels(dl) is the actual downsampling factor for each level,
  ! which is used to as the division factor to define the axes for that level.
  ! Note that the downsampling axes and domains are created at initialization based on what is
  ! requested in the parameter files (default is none) regardless of whether
  ! any downsampled diagnostics are present in the diag_table.
  ! Are downsampled diagnostics requested?
  call get_param(param_file, mdl, 'NUM_DIAG_DOWNSAMP_LEV', diag_cs%num_diag_dsamp_levels, &
                 'The number of diagnostic downsample levels to use. '//&
                 'For each level, an entry in DIAG_DOWNSAMP_LEV must be provided.', &
                 default=0)
  if (diag_cs%num_diag_dsamp_levels > 0) then
    allocate(diag_cs%diag_dsamp_levels(diag_cs%num_diag_dsamp_levels))
    call get_param(param_file, mdl, 'DIAG_DOWNSAMP_LEVS', diag_cs%diag_dsamp_levels, &
                  'A comma separated list of diagnostic downsample levels to be used. ', &
                  fail_if_missing=.true.)

    allocate(diag_cs%dsamp(diag_cs%num_diag_dsamp_levels))
    ! Initialize the global grid extents for all requested levels of diagnostics coarsening.
    allocate(G%HId(diag_cs%num_diag_dsamp_levels))
    !Allocate downsampling domains
    allocate(G%Domain%mpp_domain_d(diag_cs%num_diag_dsamp_levels))
    !Create and populated the downsampling domains and grids
    do dl=1, diag_cs%num_diag_dsamp_levels
      dlfac = diag_cs%diag_dsamp_levels(dl)
      !Create the auxiliary mpp_domain for this level of downsampled diagnostics
      !Downsample diagnostics calculations do not need halos.
      write(dsamp_domain_name, '(a,i0)') trim("MOM_domain_d"),dlfac
      call clone_MOM_domain(G%Domain, G%Domain%mpp_domain_d(dl), coarsen=dlfac, & !halo_size=0, &
                            domain_name=dsamp_domain_name)

      !Set the grid extents for this level of downsampling.
      call get_domain_extent(G%Domain, G%HId(dl)%isc, G%HId(dl)%iec, G%HId(dl)%jsc, G%HId(dl)%jec, &
                             G%HId(dl)%isd, G%HId(dl)%ied, G%HId(dl)%jsd, G%HId(dl)%jed, &
                             G%HId(dl)%isg, G%HId(dl)%ieg, G%HId(dl)%jsg, G%HId(dl)%jeg, &
                             coarsen=dl)

      ! Set array sizes for fields that are discretized at tracer cell boundaries.
      G%HId(dl)%IscB = G%HId(dl)%isc ; G%HId(dl)%JscB = G%HId(dl)%jsc
      G%HId(dl)%IsdB = G%HId(dl)%isd ; G%HId(dl)%JsdB = G%HId(dl)%jsd
      G%HId(dl)%IsgB = G%HId(dl)%isg ; G%HId(dl)%JsgB = G%HId(dl)%jsg
      if (G%symmetric) then
        G%HId(dl)%IscB = G%HId(dl)%isc-1 ; G%HId(dl)%JscB = G%HId(dl)%jsc-1
        G%HId(dl)%IsdB = G%HId(dl)%isd-1 ; G%HId(dl)%JsdB = G%HId(dl)%jsd-1
        G%HId(dl)%IsgB = G%HId(dl)%isg-1 ; G%HId(dl)%JsgB = G%HId(dl)%jsg-1
      endif
      G%HId(dl)%IecB = G%HId(dl)%iec ; G%HId(dl)%JecB = G%HId(dl)%jec
      G%HId(dl)%IedB = G%HId(dl)%ied ; G%HId(dl)%JedB = G%HId(dl)%jed
      G%HId(dl)%IegB = G%HId(dl)%ieg ; G%HId(dl)%JegB = G%HId(dl)%jeg

      !Downsample indices for diagnostics that are on a coarser grid than the model grid.
      diag_cs%dsamp(dl)%isc = G%HId(dl)%isc - (G%HId(dl)%isd-1)
      diag_cs%dsamp(dl)%iec = G%HId(dl)%iec - (G%HId(dl)%isd-1)
      diag_cs%dsamp(dl)%jsc = G%HId(dl)%jsc - (G%HId(dl)%jsd-1)
      diag_cs%dsamp(dl)%jec = G%HId(dl)%jec - (G%HId(dl)%jsd-1)
      diag_cs%dsamp(dl)%isd = G%HId(dl)%isd ; diag_cs%dsamp(dl)%ied = G%HId(dl)%ied
      diag_cs%dsamp(dl)%jsd = G%HId(dl)%jsd ; diag_cs%dsamp(dl)%jed = G%HId(dl)%jed
      diag_cs%dsamp(dl)%isg = G%HId(dl)%isg ; diag_cs%dsamp(dl)%ieg = G%HId(dl)%ieg
      diag_cs%dsamp(dl)%jsg = G%HId(dl)%jsg ; diag_cs%dsamp(dl)%jeg = G%HId(dl)%jeg
      diag_cs%dsamp(dl)%isgB = G%HId(dl)%isgB ; diag_cs%dsamp(dl)%iegB = G%HId(dl)%iegB
      diag_cs%dsamp(dl)%jsgB = G%HId(dl)%jsgB ; diag_cs%dsamp(dl)%jegB = G%HId(dl)%jegB
    enddo
  endif
  ! Initialze available diagnostic log file
  if (is_root_pe() .and. (diag_CS%available_diag_doc_unit < 0)) then
    write(this_pe,'(i6.6)') PE_here()
    doc_file_dflt = "available_diags."//this_pe
    call get_param(param_file, mdl, "AVAILABLE_DIAGS_FILE", doc_file, &
                 "A file into which to write a list of all available "//&
                 "ocean diagnostics that can be included in a diag_table.", &
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
    ! write(this_pe,'(i6.6)') PE_here()
    ! doc_file_dflt = "chksum_diag."//this_pe
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

end procedure diag_mediator_init
module procedure diag_mediator_set_OBC_info
  integer :: i, j
  do j=G%jsd,G%jed ; do i=G%IsdB,G%IedB
    diag_cs%OBC_u(I,j) = 0
    if (OBC_seg_u(I,j) > 0) diag_cs%OBC_u(I,j) = 1
    if (OBC_seg_u(I,j) < 0) diag_cs%OBC_u(I,j) = -1
  enddo ; enddo

  do J=G%JsdB,G%JedB ; do i=G%isd,G%ied
    diag_cs%OBC_v(i,J) = 0.0
    if (OBC_seg_v(i,J) > 0) diag_cs%OBC_v(i,J) = 1
    if (OBC_seg_v(i,J) < 0) diag_cs%OBC_v(i,J) = -1
  enddo ; enddo

end procedure diag_mediator_set_OBC_info
module procedure diag_set_state_ptrs
  diag_cs%h => h
  diag_cs%T => tv%T
  diag_cs%S => tv%S
  diag_cs%eqn_of_state => tv%eqn_of_state
  diag_cs%tv => tv

end procedure diag_set_state_ptrs
module procedure diag_update_remap_grids
  integer :: m
  real, dimension(:,:,:), pointer :: h_diag => NULL() ! The layer thicknesses for diagnostics [H ~> m or kg m-2]
  real, dimension(:,:,:), pointer :: T_diag => NULL() ! The layer temperatures for diagnostics [C ~> degC]
  real, dimension(:,:,:), pointer :: S_diag => NULL() ! The layer salinities for diagnostics [S ~> ppt]
  real, dimension(diag_cs%G%isd:diag_cS%G%ied, diag_cs%G%jsd:diag_cS%G%jed, diag_cs%GV%ke) :: &
    dz_diag     ! Layer vertical extents for remapping [Z ~> m]
  logical :: update_intensive_local, update_extensive_local, dz_diag_needed
  if (diag_cs%show_call_tree) call callTree_enter("diag_update_remap_grids()")

  ! Set values based on optional input arguments
  if (present(alt_h)) then
    h_diag => alt_h
  else
    h_diag => diag_cs%h
  endif

  if (present(alt_T)) then
    T_diag => alt_T
  else
    T_diag => diag_CS%T
  endif

  if (present(alt_S)) then
    S_diag => alt_S
  else
    S_diag => diag_CS%S
  endif

  ! Defaults here are based on wanting to update intensive quantities frequently as soon as the model state changes.
  ! Conversely, for extensive quantities, in an effort to close budgets and to be consistent with the total time
  ! tendency, we construct the diagnostic grid at the beginning of the baroclinic timestep and remap all extensive
  ! quantities to the same grid
  update_intensive_local = .true.
  if (present(update_intensive)) update_intensive_local = update_intensive
  update_extensive_local = .false.
  if (present(update_extensive)) update_extensive_local = update_extensive

  if (id_clock_diag_grid_updates>0) call cpu_clock_begin(id_clock_diag_grid_updates)

  if (diag_cs%diag_grid_overridden) then
    call MOM_error(FATAL, "diag_update_remap_grids was called, but current grids in "// &
                          "diagnostic structure have been overridden")
  endif

  ! Determine the diagnostic grid spacing in height units, if it is needed.
  dz_diag_needed = .false.
  if (update_intensive_local .or. update_extensive_local) then
    do m=1, diag_cs%num_diag_coords
      if (diag_cs%diag_remap_cs(m)%Z_based_coord) dz_diag_needed = .true.
    enddo
  endif
  if (dz_diag_needed) then
    call thickness_to_dz(h_diag, diag_cs%tv, dz_diag, diag_cs%G, diag_cs%GV, diag_cs%US, halo_size=1)
  endif

  if (update_intensive_local) then
    do m=1, diag_cs%num_diag_coords
      if (diag_cs%diag_remap_cs(m)%Z_based_coord) then
        call diag_remap_update(diag_cs%diag_remap_cs(m), diag_cs%G, diag_cs%GV, diag_cs%US, dz_diag, T_diag, S_diag, &
                               diag_cs%eqn_of_state, diag_cs%diag_remap_cs(m)%h)
      else
        call diag_remap_update(diag_cs%diag_remap_cs(m), diag_cs%G, diag_cs%GV, diag_cs%US, h_diag, T_diag, S_diag, &
                               diag_cs%eqn_of_state, diag_cs%diag_remap_cs(m)%h)
      endif
    enddo
  endif
  if (update_extensive_local) then
    diag_cs%h_begin(:,:,:) = diag_cs%h(:,:,:)
    if (dz_diag_needed) diag_cs%dz_begin(:,:,:) = dz_diag(:,:,:)
    do m=1, diag_cs%num_diag_coords
      if (diag_cs%diag_remap_cs(m)%Z_based_coord) then
        call diag_remap_update(diag_cs%diag_remap_cs(m), diag_cs%G, diag_cs%GV, diag_cs%US, dz_diag, T_diag, S_diag, &
                               diag_cs%eqn_of_state, diag_cs%diag_remap_cs(m)%h_extensive)
      else
        call diag_remap_update(diag_cs%diag_remap_cs(m), diag_cs%G, diag_cs%GV, diag_cs%US, h_diag, T_diag, S_diag, &
                               diag_cs%eqn_of_state, diag_cs%diag_remap_cs(m)%h_extensive)
      endif
    enddo
  endif

#if defined(DEBUG) || defined(__DO_SAFETY_CHECKS__)
  ! Keep a copy of H - used to check whether grids are up-to-date
  ! when doing remapping.
  diag_cs%h_old(:,:,:) = diag_cs%h(:,:,:)
#endif

  if (id_clock_diag_grid_updates>0) call cpu_clock_end(id_clock_diag_grid_updates)

  if (diag_cs%show_call_tree) call callTree_leave("diag_update_remap_grids()")

end procedure diag_update_remap_grids
module procedure diag_masks_set
  integer :: k
  diag_cs%mask2dT  => G%mask2dT
  diag_cs%mask2dBu => G%mask2dBu
  diag_cs%mask2dCu => G%mask2dCu
  diag_cs%mask2dCv => G%mask2dCv

  ! 3d native masks are needed by diag_manager but the native variables
  ! can only be masked 2d - for ocean points, all layers exists.
  allocate(diag_cs%mask3dTL(G%isd:G%ied,G%jsd:G%jed,1:nz))
  allocate(diag_cs%mask3dBL(G%IsdB:G%IedB,G%JsdB:G%JedB,1:nz))
  allocate(diag_cs%mask3dCuL(G%IsdB:G%IedB,G%jsd:G%jed,1:nz))
  allocate(diag_cs%mask3dCvL(G%isd:G%ied,G%JsdB:G%JedB,1:nz))
  do k=1,nz
    diag_cs%mask3dTL(:,:,k) = diag_cs%mask2dT(:,:)
    diag_cs%mask3dBL(:,:,k) = diag_cs%mask2dBu(:,:)
    diag_cs%mask3dCuL(:,:,k) = diag_cs%mask2dCu(:,:)
    diag_cs%mask3dCvL(:,:,k) = diag_cs%mask2dCv(:,:)
  enddo
  allocate(diag_cs%mask3dTi(G%isd:G%ied,G%jsd:G%jed,1:nz+1))
  allocate(diag_cs%mask3dBi(G%IsdB:G%IedB,G%JsdB:G%JedB,1:nz+1))
  allocate(diag_cs%mask3dCui(G%IsdB:G%IedB,G%jsd:G%jed,1:nz+1))
  allocate(diag_cs%mask3dCvi(G%isd:G%ied,G%JsdB:G%JedB,1:nz+1))
  do k=1,nz+1
    diag_cs%mask3dTi(:,:,k) = diag_cs%mask2dT(:,:)
    diag_cs%mask3dBi(:,:,k) = diag_cs%mask2dBu(:,:)
    diag_cs%mask3dCui(:,:,k) = diag_cs%mask2dCu(:,:)
    diag_cs%mask3dCvi(:,:,k) = diag_cs%mask2dCv(:,:)
  enddo

  ! Allocate and initialize the downsampled masks
  call downsample_diag_masks_set(G, nz, diag_cs)

end procedure diag_masks_set
module procedure set_piecemeal_extents
  call diag_cs%axesT1%piecemeal_2d%set_extents_from_array(diag_cs%mask2dT, diag_cs%missing_value)
  call diag_cs%axesB1%piecemeal_2d%set_extents_from_array(diag_cs%mask2dBu, diag_cs%missing_value)
  call diag_cs%axesCu1%piecemeal_2d%set_extents_from_array(diag_cs%mask2dCu, diag_cs%missing_value)
  call diag_cs%axesCv1%piecemeal_2d%set_extents_from_array(diag_cs%mask2dCv, diag_cs%missing_value)

  ! Piecemeal buffers for 3d axes
  call diag_cs%axesTL%piecemeal_3d%set_extents_from_array(diag_cs%mask3dTL, diag_cs%missing_value)
  call diag_cs%axesBL%piecemeal_3d%set_extents_from_array(diag_cs%mask3dBL, diag_cs%missing_value)
  call diag_cs%axesCuL%piecemeal_3d%set_extents_from_array(diag_cs%mask3dCuL, diag_cs%missing_value)
  call diag_cs%axesCvL%piecemeal_3d%set_extents_from_array(diag_cs%mask3dCvL, diag_cs%missing_value)
  call diag_cs%axesTi%piecemeal_3d%set_extents_from_array(diag_cs%mask3dTi, diag_cs%missing_value)
  call diag_cs%axesBi%piecemeal_3d%set_extents_from_array(diag_cs%mask3dBi, diag_cs%missing_value)
  call diag_cs%axesCui%piecemeal_3d%set_extents_from_array(diag_cs%mask3dCui, diag_cs%missing_value)
  call diag_cs%axesCvi%piecemeal_3d%set_extents_from_array(diag_cs%mask3dCvi, diag_cs%missing_value)

end procedure set_piecemeal_extents
module procedure diag_mediator_close_registration
  integer :: i
  if (diag_CS%available_diag_doc_unit > -1) then
    close(diag_CS%available_diag_doc_unit) ; diag_CS%available_diag_doc_unit = -2
  endif

  do i=1, diag_cs%num_diag_coords
    call diag_remap_diag_registration_closed(diag_cs%diag_remap_cs(i))
  enddo

end procedure diag_mediator_close_registration
module procedure axes_grp_end
  deallocate(axes%handles)
  if (associated(axes%mask2d)) deallocate(axes%mask2d)
  if (associated(axes%mask3d)) deallocate(axes%mask3d)
end procedure axes_grp_end
module procedure diag_mediator_end
  type(diag_type), pointer :: diag, next_diag
  integer :: i, dl
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

  do i=1, diag_cs%num_diag_coords
    call diag_remap_end(diag_cs%diag_remap_cs(i))
  enddo

  call diag_grid_storage_end(diag_cs%diag_grid_temp)
  if (associated(diag_cs%mask3dTL))  deallocate(diag_cs%mask3dTL)
  if (associated(diag_cs%mask3dBL))  deallocate(diag_cs%mask3dBL)
  if (associated(diag_cs%mask3dCuL)) deallocate(diag_cs%mask3dCuL)
  if (associated(diag_cs%mask3dCvL)) deallocate(diag_cs%mask3dCvL)
  if (associated(diag_cs%mask3dTi))  deallocate(diag_cs%mask3dTi)
  if (associated(diag_cs%mask3dBi))  deallocate(diag_cs%mask3dBi)
  if (associated(diag_cs%mask3dCui)) deallocate(diag_cs%mask3dCui)
  if (associated(diag_cs%mask3dCvi)) deallocate(diag_cs%mask3dCvi)
  do dl=1, diag_cs%num_diag_dsamp_levels
    if (associated(diag_cs%dsamp(dl)%mask2dT))   deallocate(diag_cs%dsamp(dl)%mask2dT)
    if (associated(diag_cs%dsamp(dl)%mask2dBu))  deallocate(diag_cs%dsamp(dl)%mask2dBu)
    if (associated(diag_cs%dsamp(dl)%mask2dCu))  deallocate(diag_cs%dsamp(dl)%mask2dCu)
    if (associated(diag_cs%dsamp(dl)%mask2dCv))  deallocate(diag_cs%dsamp(dl)%mask2dCv)
    if (associated(diag_cs%dsamp(dl)%mask3dTL))  deallocate(diag_cs%dsamp(dl)%mask3dTL)
    if (associated(diag_cs%dsamp(dl)%mask3dBL))  deallocate(diag_cs%dsamp(dl)%mask3dBL)
    if (associated(diag_cs%dsamp(dl)%mask3dCuL)) deallocate(diag_cs%dsamp(dl)%mask3dCuL)
    if (associated(diag_cs%dsamp(dl)%mask3dCvL)) deallocate(diag_cs%dsamp(dl)%mask3dCvL)
    if (associated(diag_cs%dsamp(dl)%mask3dTi))  deallocate(diag_cs%dsamp(dl)%mask3dTi)
    if (associated(diag_cs%dsamp(dl)%mask3dBi))  deallocate(diag_cs%dsamp(dl)%mask3dBi)
    if (associated(diag_cs%dsamp(dl)%mask3dCui)) deallocate(diag_cs%dsamp(dl)%mask3dCui)
    if (associated(diag_cs%dsamp(dl)%mask3dCvi)) deallocate(diag_cs%dsamp(dl)%mask3dCvi)

    do i=1,diag_cs%num_diag_coords
      if (associated(diag_cs%dsamp(dl)%remap_axesTL(i)%dsamp(dl)%mask3d)) &
        deallocate(diag_cs%dsamp(dl)%remap_axesTL(i)%dsamp(dl)%mask3d)
      if (associated(diag_cs%dsamp(dl)%remap_axesCuL(i)%dsamp(dl)%mask3d)) &
        deallocate(diag_cs%dsamp(dl)%remap_axesCuL(i)%dsamp(dl)%mask3d)
      if (associated(diag_cs%dsamp(dl)%remap_axesCvL(i)%dsamp(dl)%mask3d)) &
        deallocate(diag_cs%dsamp(dl)%remap_axesCvL(i)%dsamp(dl)%mask3d)
      if (associated(diag_cs%dsamp(dl)%remap_axesBL(i)%dsamp(dl)%mask3d)) &
        deallocate(diag_cs%dsamp(dl)%remap_axesBL(i)%dsamp(dl)%mask3d)
      if (associated(diag_cs%dsamp(dl)%remap_axesTi(i)%dsamp(dl)%mask3d)) &
        deallocate(diag_cs%dsamp(dl)%remap_axesTi(i)%dsamp(dl)%mask3d)
      if (associated(diag_cs%dsamp(dl)%remap_axesCui(i)%dsamp(dl)%mask3d)) &
        deallocate(diag_cs%dsamp(dl)%remap_axesCui(i)%dsamp(dl)%mask3d)
      if (associated(diag_cs%dsamp(dl)%remap_axesCvi(i)%dsamp(dl)%mask3d)) &
        deallocate(diag_cs%dsamp(dl)%remap_axesCvi(i)%dsamp(dl)%mask3d)
      if (associated(diag_cs%dsamp(dl)%remap_axesBi(i)%dsamp(dl)%mask3d)) &
        deallocate(diag_cs%dsamp(dl)%remap_axesBi(i)%dsamp(dl)%mask3d)
    enddo
  enddo

  ! axes_grp masks may point to diag_cs masks, so do these after mask dealloc
  do i=1, diag_cs%num_diag_coords
    call axes_grp_end(diag_cs%remap_axesZL(i))
    call axes_grp_end(diag_cs%remap_axesZi(i))
    call axes_grp_end(diag_cs%remap_axesTL(i))
    call axes_grp_end(diag_cs%remap_axesTi(i))
    call axes_grp_end(diag_cs%remap_axesBL(i))
    call axes_grp_end(diag_cs%remap_axesBi(i))
    call axes_grp_end(diag_cs%remap_axesCuL(i))
    call axes_grp_end(diag_cs%remap_axesCui(i))
    call axes_grp_end(diag_cs%remap_axesCvL(i))
    call axes_grp_end(diag_cs%remap_axesCvi(i))
  enddo

  if (diag_cs%num_diag_coords > 0) then
    deallocate(diag_cs%remap_axesZL)
    deallocate(diag_cs%remap_axesZi)
    deallocate(diag_cs%remap_axesTL)
    deallocate(diag_cs%remap_axesTi)
    deallocate(diag_cs%remap_axesBL)
    deallocate(diag_cs%remap_axesBi)
    deallocate(diag_cs%remap_axesCuL)
    deallocate(diag_cs%remap_axesCui)
    deallocate(diag_cs%remap_axesCvL)
    deallocate(diag_cs%remap_axesCvi)
  endif

  do dl=1, diag_cs%num_diag_dsamp_levels
    if (allocated(diag_cs%dsamp(dl)%remap_axesTL)) &
      deallocate(diag_cs%dsamp(dl)%remap_axesTL)
    if (allocated(diag_cs%dsamp(dl)%remap_axesTi)) &
      deallocate(diag_cs%dsamp(dl)%remap_axesTi)
    if (allocated(diag_cs%dsamp(dl)%remap_axesBL)) &
      deallocate(diag_cs%dsamp(dl)%remap_axesBL)
    if (allocated(diag_cs%dsamp(dl)%remap_axesBi)) &
      deallocate(diag_cs%dsamp(dl)%remap_axesBi)
    if (allocated(diag_cs%dsamp(dl)%remap_axesCuL)) &
      deallocate(diag_cs%dsamp(dl)%remap_axesCuL)
    if (allocated(diag_cs%dsamp(dl)%remap_axesCui)) &
      deallocate(diag_cs%dsamp(dl)%remap_axesCui)
    if (allocated(diag_cs%dsamp(dl)%remap_axesCvL)) &
      deallocate(diag_cs%dsamp(dl)%remap_axesCvL)
    if (allocated(diag_cs%dsamp(dl)%remap_axesCvi)) &
      deallocate(diag_cs%dsamp(dl)%remap_axesCvi)
  enddo


#if defined(DEBUG) || defined(__DO_SAFETY_CHECKS__)
  deallocate(diag_cs%h_old)
#endif

  if (present(end_diag_manager)) then
    if (end_diag_manager) call MOM_diag_manager_end(time)
  endif

end procedure diag_mediator_end
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
  if (present(dimensions)) then
    if (len(trim(dimensions)) > 0) then
      call describe_option("dimensions", dimensions, diag_CS)
    endif
  endif
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
module procedure diag_grid_storage_init
  integer :: m, nz
  grid_storage%num_diag_coords = diag%num_diag_coords

  ! Don't do anything else if there are no remapped coordinates
  if (grid_storage%num_diag_coords < 1) return

  ! Allocate memory for the native space
  allocate( grid_storage%h_state(G%isd:G%ied, G%jsd:G%jed, GV%ke))
  ! Allocate diagnostic remapping structures
  allocate(grid_storage%diag_grids(diag%num_diag_coords))
  ! Loop through and allocate memory for the grid on each target coordinate
  do m = 1, diag%num_diag_coords
    nz = diag%diag_remap_cs(m)%nz
    allocate(grid_storage%diag_grids(m)%h(G%isd:G%ied,G%jsd:G%jed, nz))
  enddo

end procedure diag_grid_storage_init
module procedure diag_copy_diag_to_storage
  integer :: m
  if (grid_storage%num_diag_coords < 1) return

  grid_storage%h_state(:,:,:) = h_state(:,:,:)
  do m = 1,grid_storage%num_diag_coords
    if (diag%diag_remap_cs(m)%nz > 0) &
      grid_storage%diag_grids(m)%h(:,:,:) = diag%diag_remap_cs(m)%h(:,:,:)
  enddo

end procedure diag_copy_diag_to_storage
module procedure diag_copy_storage_to_diag
  integer :: m
  if (grid_storage%num_diag_coords < 1) return

  diag%diag_grid_overridden = .true.
  do m = 1,grid_storage%num_diag_coords
    if (diag%diag_remap_cs(m)%nz > 0) &
      diag%diag_remap_cs(m)%h(:,:,:) = grid_storage%diag_grids(m)%h(:,:,:)
  enddo

end procedure diag_copy_storage_to_diag
module procedure diag_save_grids
  integer :: m
  if (diag%num_diag_coords < 1) return

  do m = 1,diag%num_diag_coords
    if (diag%diag_remap_cs(m)%nz > 0) &
      diag%diag_grid_temp%diag_grids(m)%h(:,:,:) = diag%diag_remap_cs(m)%h(:,:,:)
  enddo

end procedure diag_save_grids
module procedure diag_restore_grids
  integer :: m
  if (diag%num_diag_coords < 1) return

  diag%diag_grid_overridden = .false.
  do m = 1,diag%num_diag_coords
    if (diag%diag_remap_cs(m)%nz > 0) &
      diag%diag_remap_cs(m)%h(:,:,:) = diag%diag_grid_temp%diag_grids(m)%h(:,:,:)
  enddo

end procedure diag_restore_grids
module procedure diag_grid_storage_end
  integer :: m
  if (grid_storage%num_diag_coords < 1) return

  ! Deallocate memory for the native space
  deallocate(grid_storage%h_state)
  ! Loop through and deallocate memory for the grid on each target coordinate
  do m = 1, grid_storage%num_diag_coords
    deallocate(grid_storage%diag_grids(m)%h)
  enddo
  ! Deallocate diagnostic remapping structures
  deallocate(grid_storage%diag_grids)
end procedure diag_grid_storage_end
module procedure downsample_diag_masks_set
  integer :: k, dl, dlfac
  do dl=1, diag_cs%num_diag_dsamp_levels
    dlfac = diag_cs%diag_dsamp_levels(dl) !Actual downsampling factor for this level
    ! 2d mask
    call downsample_mask(G%mask2dT, diag_cs%dsamp(dl)%mask2dT,  dlfac, MMP, G%isc, G%jsc, G%isd, G%jsd, &
                         G%HId(dl)%isc, G%HId(dl)%iec, G%HId(dl)%jsc, G%HId(dl)%jec, G%HId(dl)%isd, G%HId(dl)%ied, &
                         G%HId(dl)%jsd, G%HId(dl)%jed)
    call downsample_mask(G%mask2dBu, diag_cs%dsamp(dl)%mask2dBu, dlfac, PPP,G%IscB, G%JscB, G%IsdB, G%JsdB, &
                         G%HId(dl)%IscB,G%HId(dl)%IecB, G%HId(dl)%JscB,G%HId(dl)%JecB,G%HId(dl)%IsdB,G%HId(dl)%IedB, &
                         G%HId(dl)%JsdB,G%HId(dl)%JedB)
    call downsample_mask(G%mask2dCu, diag_cs%dsamp(dl)%mask2dCu, dlfac, PMP, G%IscB, G%jsc, G%IsdB, G%jsd, &
                         G%HId(dl)%IscB,G%HId(dl)%IecB, G%HId(dl)%jsc, G%HId(dl)%jec,G%HId(dl)%IsdB,G%HId(dl)%IedB, &
                         G%HId(dl)%jsd, G%HId(dl)%jed)
    call downsample_mask(G%mask2dCv, diag_cs%dsamp(dl)%mask2dCv, dlfac, MPP, G %isc ,G%JscB, G%isd, G%JsdB, &
                         G%HId(dl)%isc ,G%HId(dl)%iec, G%HId(dl)%JscB,G%HId(dl)%JecB,G%HId(dl)%isd ,G%HId(dl)%ied, &
                         G%HId(dl)%JsdB,G%HId(dl)%JedB)
    ! 3d native masks are needed by diag_manager but the native variables
    ! can only be masked 2d - for ocean points, all layers exists.
    allocate(diag_cs%dsamp(dl)%mask3dTL(G%HId(dl)%isd:G%HId(dl)%ied,G%HId(dl)%jsd:G%HId(dl)%jed,1:nz))
    allocate(diag_cs%dsamp(dl)%mask3dBL(G%HId(dl)%IsdB:G%HId(dl)%IedB,G%HId(dl)%JsdB:G%HId(dl)%JedB,1:nz))
    allocate(diag_cs%dsamp(dl)%mask3dCuL(G%HId(dl)%IsdB:G%HId(dl)%IedB,G%HId(dl)%jsd:G%HId(dl)%jed,1:nz))
    allocate(diag_cs%dsamp(dl)%mask3dCvL(G%HId(dl)%isd:G%HId(dl)%ied,G%HId(dl)%JsdB:G%HId(dl)%JedB,1:nz))
    do k=1,nz
      diag_cs%dsamp(dl)%mask3dTL(:,:,k) = diag_cs%dsamp(dl)%mask2dT(:,:)
      diag_cs%dsamp(dl)%mask3dBL(:,:,k) = diag_cs%dsamp(dl)%mask2dBu(:,:)
      diag_cs%dsamp(dl)%mask3dCuL(:,:,k) = diag_cs%dsamp(dl)%mask2dCu(:,:)
      diag_cs%dsamp(dl)%mask3dCvL(:,:,k) = diag_cs%dsamp(dl)%mask2dCv(:,:)
    enddo
    allocate(diag_cs%dsamp(dl)%mask3dTi(G%HId(dl)%isd:G%HId(dl)%ied,G%HId(dl)%jsd:G%HId(dl)%jed,1:nz+1))
    allocate(diag_cs%dsamp(dl)%mask3dBi(G%HId(dl)%IsdB:G%HId(dl)%IedB,G%HId(dl)%JsdB:G%HId(dl)%JedB,1:nz+1))
    allocate(diag_cs%dsamp(dl)%mask3dCui(G%HId(dl)%IsdB:G%HId(dl)%IedB,G%HId(dl)%jsd:G%HId(dl)%jed,1:nz+1))
    allocate(diag_cs%dsamp(dl)%mask3dCvi(G%HId(dl)%isd:G%HId(dl)%ied,G%HId(dl)%JsdB:G%HId(dl)%JedB,1:nz+1))
    do k=1,nz+1
      diag_cs%dsamp(dl)%mask3dTi(:,:,k) = diag_cs%dsamp(dl)%mask2dT(:,:)
      diag_cs%dsamp(dl)%mask3dBi(:,:,k) = diag_cs%dsamp(dl)%mask2dBu(:,:)
      diag_cs%dsamp(dl)%mask3dCui(:,:,k) = diag_cs%dsamp(dl)%mask2dCu(:,:)
      diag_cs%dsamp(dl)%mask3dCvi(:,:,k) = diag_cs%dsamp(dl)%mask2dCv(:,:)
    enddo
  enddo
end procedure downsample_diag_masks_set
module procedure downsample_diag_indices_get
  integer :: dszi, cszi, dszj, cszj, f1, f2, dlfac
  character(len=500) :: mesg
  logical, save :: first_check = .true.
  dlfac = diag_cs%diag_dsamp_levels(dl) !Actual downsampling factor for this level
  if (first_check) then
    if (mod(diag_cs%ie-diag_cs%is+1, dlfac) /= 0 .OR. &
        mod(diag_cs%je-diag_cs%js+1, dlfac) /= 0) then
      write (mesg,*) "Non-commensurate downsampled domain is not supported. "//&
             "Please choose a layout such that NIGLOBAL/Layout_X and NJGLOBAL/Layout_Y are both divisible by dL=", &
             dlfac,&
             " Current domain extents: ", diag_cs%is,diag_cs%ie, diag_cs%js,diag_cs%je
      call MOM_error(FATAL,"downsample_diag_indices_get: "//trim(mesg))
    endif
    first_check = .false.
  endif

  cszi = diag_cs%dsamp(dl)%iec-diag_cs%dsamp(dl)%isc +1 ; dszi = diag_cs%dsamp(dl)%ied-diag_cs%dsamp(dl)%isd +1
  cszj = diag_cs%dsamp(dl)%jec-diag_cs%dsamp(dl)%jsc +1 ; dszj = diag_cs%dsamp(dl)%jed-diag_cs%dsamp(dl)%jsd +1
  !isv = diag_cs%dsamp(dl)%isc ; iev = diag_cs%dsamp(dl)%iec
  !jsv = diag_cs%dsamp(dl)%jsc ; jev = diag_cs%dsamp(dl)%jec

  f1 = fo1/dlfac
  f2 = fo2/dlfac
  ! Correction for the symmetric case
  if (diag_cs%G%symmetric) then
    f1 = f1 + mod(fo1,dlfac)
    f2 = f2 + mod(fo2,dlfac)
  endif

  ! Find the range of indices in the downsampled computational domain.
  if ( f1 == dszi ) then
    isv = diag_cs%dsamp(dl)%isc ; iev = diag_cs%dsamp(dl)%iec   ! Field on Data domain, take compute domain indices
  elseif ( f1 == dszi + 1 ) then
    isv = diag_cs%dsamp(dl)%isc ; iev = diag_cs%dsamp(dl)%iec+1   ! Symmetric data domain
  elseif ( f1 == cszi) then
    isv = 1 ; iev = (diag_cs%dsamp(dl)%iec-diag_cs%dsamp(dl)%isc) +1  ! Computational domain
  elseif ( f1 == cszi + 1 ) then
    isv = 1 ; iev = (diag_cs%dsamp(dl)%iec-diag_cs%dsamp(dl)%isc) +2  ! Symmetric computational domain
  else
    write (mesg,*) " dl =",dl,",dL =",dlfac,",fo1 =",fo1," f1 =",f1," peculiar size for diag field in i-direction\n"//&
          "does not match one of ", cszi, cszi+1, dszi, dszi+1
    call MOM_error(FATAL,"downsample_diag_indices_get: "//trim(mesg))
  endif
  if ( f2 == dszj ) then
    jsv = diag_cs%dsamp(dl)%jsc ; jev = diag_cs%dsamp(dl)%jec     ! Data domain
  elseif ( f2 == dszj + 1 ) then
    jsv = diag_cs%dsamp(dl)%jsc ; jev = diag_cs%dsamp(dl)%jec+1   ! Symmetric data domain
  elseif ( f2 == cszj) then
    jsv = 1 ; jev = (diag_cs%dsamp(dl)%jec-diag_cs%dsamp(dl)%jsc) +1  ! Computational domain
  elseif ( f2 == cszj + 1 ) then
    jsv = 1 ; jev = (diag_cs%dsamp(dl)%jec-diag_cs%dsamp(dl)%jsc) +2  ! Symmetric computational domain
  else
    write (mesg,*) " dl =",dl,",dL =",dlfac,",fo2 =",fo2," f2 =",f2," peculiar size for diag field in j-direction\n"//&
          "does not match one of ", cszj, cszj+1, dszj, dszj+1
    call MOM_error(FATAL,"downsample_diag_indices_get: "//trim(mesg))
  endif
end procedure downsample_diag_indices_get
module procedure downsample_diag_field_3d
  real, dimension(:,:,:), pointer :: locmask ! A pointer to the mask [nondim]
  integer :: f1, f2, isv_o, jsv_o
  locmask => NULL()
  ! Get the correct indices corresponding to input field based on its shape.
  f1 = size(locfield, 1)
  f2 = size(locfield, 2)
  ! Save the extents of the original (fine) domain
  isv_o = isv ; jsv_o = jsv
  ! Get the shape of the downsampled field and overwrite isv, iev, jsv and jev with them
  call downsample_diag_indices_get(f1, f2, dl, diag_cs, isv, iev, jsv, jev)
  ! Set the pointer to the non-downsampled mask, which must be associated and initialized
  if (present(mask)) then
    locmask => mask
  elseif (associated(diag%axes%mask3d)) then
    locmask => diag%axes%mask3d
  else
    call MOM_error(FATAL, "downsample_diag_field_3d: Cannot downsample without a mask!!! ")
  endif
  call downsample_field(locfield, locfield_dsamp, diag_cs%diag_dsamp_levels(dl), diag%xyz_method, &
                        locmask, diag_cs, diag, isv_o, jsv_o, isv, iev, jsv, jev)

end procedure downsample_diag_field_3d
module procedure downsample_diag_field_2d
  real, dimension(:,:), pointer :: locmask ! A pointer to the mask [nondim]
  integer :: f1, f2, isv_o, jsv_o
  locmask => NULL()
  ! Get the correct indices corresponding to input field based on its shape.
  f1 = size(locfield,1)
  f2 = size(locfield,2)
  ! Save the extents of the original (fine) domain
  isv_o = isv ; jsv_o = jsv
  ! Get the shape of the downsampled field and overwrite isv, iev, jsv and jev with them
  call downsample_diag_indices_get(f1, f2, dl, diag_cs, isv, iev, jsv, jev)
  ! Set the non-downsampled mask, it must be associated and initialized
  if (present(mask)) then
    locmask => mask
  elseif (associated(diag%axes%mask2d)) then
    locmask => diag%axes%mask2d
  else
    call MOM_error(FATAL, "downsample_diag_field_2d: Cannot downsample without a mask!!! ")
  endif

  call downsample_field(locfield, locfield_dsamp, diag_cs%diag_dsamp_levels(dl), diag%xyz_method, &
                        locmask, diag_cs, diag, isv_o, jsv_o, isv, iev, jsv, jev)

end procedure downsample_diag_field_2d
module procedure downsample_field_3d
  character(len=240) :: mesg
  integer :: i, j, k, i_dn, j_dn, ks, ke, i0, j0, f1, f2, f_in1, f_in2
  integer :: ii, jj  ! The index locations on the full grid that contribute to the averages.
  integer :: i0_off, j0_off  ! The starting point offsets between full array and reduced array
  real :: wt(dL,dL) ! The nondimensional, area-, volume- or mass-based weight for an input
  real :: wtd_field(dL,dL) ! The weighted field to sum, in [A ~> a], [A L2 ~> a m2],
  real :: wt_1d(dL) ! The nondimensional, area-, volume- or mass-based weight for an input
  real :: wtd_field_1d(dL) ! The weighted field to sum, in [A ~> a], [A L2 ~> a m2],
  real :: ave       ! The running sum of the average, in [A ~> a], [A L2 ~> a m2],
  real :: weight    ! The nondimensional, area-, volume- or mass-based weight for an input
  real :: total_weight ! The sum of weights contributing to a point [nondim], [L2 ~> m2],
  real :: h_face    ! The thickness at a velocity face [H ~> m or kg m-2]
  real :: eps_vol   ! A negligibly small volume or mass [H L2 ~> m3 or kg]
  real :: eps_area  ! A negligibly small area [L2 ~> m2]
  real :: eps_face  ! A negligibly small face area [H L ~> m2 or kg m-1]
  logical :: naive  ! If true, use naive rotatially variant sums to reproduct previous answers.
  ks = 1 ; ke = size(field_in,3)

  ! It would be better to use a max with eps_vol instead of adding it into the denominator.
  eps_face = 1.0e-20 * diag_cs%G%US%m_to_L * diag_cs%GV%m_to_H
  eps_area = 1.0e-20 * diag_cs%G%US%m_to_L**2
  eps_vol = 1.0e-20 * diag_cs%G%US%m_to_L**2 * diag_cs%GV%m_to_H

  naive = .not.diag_CS%symmetric_downsample_sums

  ! Allocate the down sampled field on the down sampled data domain
!  allocate(field_out(diag_cs%dsamp(dl)%isd:diag_cs%dsamp(dl)%ied,diag_cs%dsamp(dl)%jsd:diag_cs%dsamp(dl)%jed,ks:ke))
!  allocate(field_out(1:size(field_in,1)/dl,1:size(field_in,2)/dl,ks:ke))
  f_in1 = size(field_in, 1)
  f_in2 = size(field_in, 2)
  f1 = f_in1 / dL
  f2 = f_in2 / dL
  ! Correction for the symmetric case
  if (diag_cs%G%symmetric) then
    f1 = f1 + mod(f_in1, dL)
    f2 = f2 + mod(f_in2, dL)
  endif
  allocate(field_out(1:f1,1:f2,ks:ke), source=0.0)

  ! These are the starting point offsets between full array and reduced array indices when i or j is 0.
  i0_off = (isv_o-1) - dL*isv_d
  j0_off = (jsv_o-1) - dL*jsv_d

  ! Fill the down sampled field on the down sampled diagnostics (almost always compute) domain
  if (method == MMM) then
    do k=ks,ke ; do j=jsv_d,jev_d ; do i=isv_d,iev_d
      do j_dn=1,dL ; do i_dn=1,dL
        ! ii and jj are the index locations on the full grid that contribute to the averages.
        jj = j_dn + (dL*j + j0_off) ; ii = i_dn + (dL*i + i0_off)
        wt(i_dn,j_dn) = mask(ii,jj,k) * diag_cs%G%areaT(ii,jj) * diag_cs%h(ii,jj,k)
        wtd_field(i_dn,j_dn) = field_in(ii,jj,k) * wt(i_dn,j_dn)
      enddo ; enddo
      field_out(i,j,k) = square_sum(wtd_field(1:dL,1:dL), dL, naive) / &
                        (square_sum(wt(1:dL,1:dL), dL, naive) + eps_vol) ! Eps_vol avoids division by 0.
    enddo ; enddo ; enddo
  elseif (method == SSS) then   ! e.g., volcello
    do k=ks,ke ; do j=jsv_d,jev_d ; do i=isv_d,iev_d
      do j_dn=1,dL ; do i_dn=1,dL
        jj = j_dn + (dL*j + j0_off) ; ii = i_dn + (dL*i + i0_off)
        wtd_field(i_dn,j_dn) = field_in(ii,jj,k) * mask(ii,jj,k)
      enddo ; enddo
      field_out(i,j,k)  = square_sum(wtd_field(1:dL,1:dL), dL, naive) ! This is a masked sum.
    enddo ; enddo ; enddo
  elseif (method == MMP .or. method == MMS) then   ! e.g., T_advection_xy
    do k=ks,ke ; do j=jsv_d,jev_d ; do i=isv_d,iev_d
      do j_dn=1,dL ; do i_dn=1,dL
        ! ii and jj are the index locations on the full grid that contribute to the averages.
        jj = j_dn + (dL*j + j0_off) ; ii = i_dn + (dL*i + i0_off)
        wt(i_dn,j_dn) = mask(ii,jj,k) * diag_cs%G%areaT(ii,jj)
        wtd_field(i_dn,j_dn) = field_in(ii,jj,k) * wt(i_dn,j_dn)
      enddo ; enddo
      field_out(i,j,k) = square_sum(wtd_field(1:dL,1:dL), dL, naive) / &
                        (square_sum(wt(1:dL,1:dL), dL, naive) + eps_area) ! Eps_area avoids division by 0.
    enddo ; enddo ; enddo
  elseif (method == PMM) then
    do k=ks,ke ; do j=jsv_d,jev_d ; do i=isv_d,iev_d
      II = dL*I + I0_off + (dL-1)
      do j_dn=1,dL
        jj = j_dn + (dL*j + j0_off)
        if (diag_cs%OBC_u(II,jj) == 0) then    ! This is not an OBC face.
          h_face = 0.5*(diag_cs%h(ii,jj,k) + diag_cs%h(ii+1,jj,k))
        elseif (diag_cs%OBC_u(II,jj) < 0) then ! This is a western OBC face.
          h_face = diag_cs%h(ii+1,jj,k)
        else  ! (diag_cs%OBC_u(II,jj) > 0)     ! This is an eastern OBC face.
          h_face = diag_cs%h(ii,jj,k)
        endif
        wt_1d(j_dn) = mask(II,jj,k) * diag_cs%G%dyCu(II,jj) * h_face
        wtd_field_1d(j_dn) = field_in(II,jj,k) * wt_1d(j_dn)
      enddo
      field_out(I,j,k)  = sum_1d(wtd_field_1d(1:dL), dL) / &
                         (sum_1d(wt_1d(1:dL), dL) + eps_face)  ! Eps_face avoids division by 0.
    enddo ; enddo ; enddo
  elseif (method == PSS) then    ! e.g. umo
    do k=ks,ke ; do j=jsv_d,jev_d ; do i=isv_d,iev_d
      II = dL*I + I0_off + (dL-1)
      do j_dn=1,dL
        jj = j_dn + (dL*j + j0_off)
        wtd_field_1d(j_dn) = field_in(II,jj,k) * mask(II,jj,k)
      enddo
      field_out(I,j,k) = sum_1d(wtd_field_1d(1:dL), dL)   ! This is a masked sum.
    enddo ; enddo ; enddo
  elseif (method == SPS) then   ! e.g. vmo
    do k=ks,ke ; do J=jsv_d,jev_d ; do i=isv_d,iev_d
      JJ = dL*J + J0_off + (dL-1)
      do i_dn=1,dL
        ii = i_dn + (dL*i + i0_off)
        wtd_field_1d(i_dn) = field_in(ii,JJ,k) * mask(ii,JJ,k)
      enddo
      field_out(i,J,k) = sum_1d(wtd_field_1d(1:dL), dL)   ! This is a masked sum.
    enddo ; enddo ; enddo
  elseif (method == MPM) then
    do k=ks,ke ; do J=jsv_d,jev_d ; do i=isv_d,iev_d
      JJ = dL*J + J0_off + (dL-1)
      do i_dn=1,dL
        ii = i_dn + (dL*i + i0_off)
        if (diag_cs%OBC_v(ii,JJ) == 0) then    ! This is not an OBC face.
          h_face = 0.5*(diag_cs%h(ii,jj,k) + diag_cs%h(ii,jj+1,k))
        elseif (diag_cs%OBC_v(ii,JJ) < 0) then ! This is a southern OBC face.
          h_face = diag_cs%h(ii,jj+1,k)
        else  ! (diag_cs%OBC_v(ii,JJ) > 0)     ! This is a northern OBC face.
          h_face = diag_cs%h(ii,jj,k)
        endif
        wt_1d(i_dn) = mask(ii,JJ,k) * diag_cs%G%dxCv(ii,JJ) * h_face
        wtd_field_1d(i_dn) = field_in(ii,JJ,k) * wt_1d(i_dn)
      enddo
      field_out(i,J,k)  = sum_1d(wtd_field_1d(1:dL), dL) / &
                         (sum_1d(wt_1d(1:dL), dL) + eps_face)  ! Eps_face avoids division by 0.
    enddo ; enddo ; enddo
  else
    write (mesg,*) " unknown sampling method: ",method
    call MOM_error(FATAL, "downsample_field_3d: "//trim(mesg)//" "//trim(diag%debug_str))
  endif

end procedure downsample_field_3d
module procedure downsample_field_2d
  character(len=240) :: mesg
  integer :: i, j, i_dn, j_dn, i0, j0, f1, f2, f_in1, f_in2
  integer :: ii, jj  ! The index locations on the full grid that contribute to the averages.
  integer :: i0_off, j0_off  ! The starting point offsets between full array and reduced array
  real :: wt(dL,dL) ! The nondimensional, area-, volume- or mass-based weight for an input
  real :: wtd_field(dL,dL) ! The weighted field to sum, in [A ~> a], [A L2 ~> a m2],
  real :: wt_1d(dL) ! The nondimensional, area-, volume- or mass-based weight for an input
  real :: wtd_field_1d(dL) ! The weighted field to sum, in [A ~> a], [A L2 ~> a m2],
  real :: ave       ! The running sum of the average, in [A ~> a] or [A L2 ~> a m2]
  real :: weight    ! The nondimensional or area-weighted weight for an input value [nondim] or [L2 ~> m2]
  real :: total_weight ! The sum of weights contributing to a point [nondim] or [L2 ~> m2]
  real :: eps_area  ! A negligibly small area [L2 ~> m2]
  real :: eps_len   ! A negligibly small horizontal length [L ~> m]
  logical :: naive  ! If true, use naive rotatially variant sums to reproduct previous answers.
  eps_len = 1.0e-20 * diag_cs%G%US%m_to_L
  eps_area = 1.0e-20 * diag_cs%G%US%m_to_L**2

  naive = .not.diag_CS%symmetric_downsample_sums

  ! Allocate the down sampled field on the down sampled data domain
!  allocate(field_out(diag_cs%dsamp(dl)%isd:diag_cs%dsamp(dl)%ied,diag_cs%dsamp(dl)%jsd:diag_cs%dsamp(dl)%jed))
!  allocate(field_out(1:size(field_in,1)/dl,1:size(field_in,2)/dl))
  ! Fill the down sampled field on the down sampled diagnostics (almost always compute) domain
  f_in1 = size(field_in,1)
  f_in2 = size(field_in,2)
  f1 = f_in1 / dL
  f2 = f_in2 / dL
  ! Correction for the symmetric case
  if (diag_cs%G%symmetric) then
    f1 = f1 + mod(f_in1,dl)
    f2 = f2 + mod(f_in2,dl)
  endif
  allocate(field_out(1:f1,1:f2))

  ! These are the starting point offsets between full array and reduced array indices when i or j is 0.
  i0_off = (isv_o-1) - dL*isv_d
  j0_off = (jsv_o-1) - dL*jsv_d

  if (method == MMP) then
    do j=jsv_d,jev_d ; do i=isv_d,iev_d
      do j_dn=1,dL ; do i_dn=1,dL
        ! ii and jj are the index locations on the full grid that contribute to the averages.
        jj = j_dn + (dL*j + j0_off) ; ii = i_dn + (dL*i + i0_off)
        wt(i_dn,j_dn) = mask(ii,jj) * diag_cs%G%areaT(ii,jj)
        wtd_field(i_dn,j_dn) = field_in(ii,jj) * wt(i_dn,j_dn)
      enddo ; enddo
      field_out(i,j) = square_sum(wtd_field(1:dL,1:dL), dL, naive) / &
                      (square_sum(wt(1:dL,1:dL), dL, naive) + eps_area) ! Eps_area avoids division by 0.
    enddo ; enddo
  elseif (method == SSP) then    ! e.g., T_dfxy_cont_tendency_2d
    do j=jsv_d,jev_d ; do i=isv_d,iev_d
      do j_dn=1,dL ; do i_dn=1,dL
        jj = j_dn + (dL*j + j0_off) ; ii = i_dn + (dL*i + i0_off)
        wtd_field(i_dn,j_dn) = field_in(ii,jj) * mask(ii,jj)
      enddo ; enddo
      field_out(i,j)  = square_sum(wtd_field(1:dL,1:dL), dL, naive) ! This is a masked sum.
    enddo ; enddo
  elseif (method == PSP) then   ! e.g., umo_2d
    do j=jsv_d,jev_d ; do I=isv_d,iev_d
      II = dL*I + i0_off + (dL-1)
      do j_dn=1,dL
        jj = j_dn + (dL*j + j0_off)
        wtd_field_1d(j_dn) = field_in(II,jj) * mask(II,jj)
      enddo
      field_out(I,j) = sum_1d(wtd_field_1d(1:dL), dL)   ! This is a masked sum.
    enddo ; enddo
  elseif (method == SPP) then   ! e.g., vmo_2d
    do J=jsv_d,jev_d ; do i=isv_d,iev_d
      JJ = dL*J + J0_off + (dL-1)
      do i_dn=1,dL
        ii = i_dn + (dL*i + i0_off)
        wtd_field_1d(i_dn) = field_in(ii,JJ) * mask(ii,JJ)
      enddo
      field_out(i,J) = sum_1d(wtd_field_1d(1:dL), dL)   ! This is a masked sum.
    enddo ; enddo
  elseif (method == PMP) then
    do j=jsv_d,jev_d ; do I=isv_d,iev_d
      II = dL*I + I0_off + (dL-1)
      do j_dn=1,dL
        jj = j_dn + (dL*j + j0_off)
        ! Should this weight include the total thickness interpolated to velocity points?
        wt_1d(j_dn) = mask(II,jj) * diag_cs%G%dyCu(II,jj)
        wtd_field_1d(j_dn) = field_in(II,jj) * wt_1d(j_dn)
      enddo
      field_out(I,j) = sum_1d(wtd_field_1d(1:dL), dL) / &
                      (sum_1d(wt_1d(1:dL), dL) + eps_len)  ! Eps_len avoids division by 0.
    enddo ; enddo
  elseif (method == MPP) then
    do J=jsv_d,jev_d ; do i=isv_d,iev_d
      JJ = dL*J + J0_off + (dL-1)
      do i_dn=1,dL
        ii = i_dn + (dL*i + i0_off)
        ! Should this weight include the total thickness interpolated to velocity points?
        wt_1d(i_dn) = mask(ii,JJ) * diag_cs%G%dxCv(ii,JJ)
        wtd_field_1d(i_dn) = field_in(ii,JJ) * wt_1d(i_dn)
      enddo
      field_out(i,J) = sum_1d(wtd_field_1d(1:dL), dL) / &
                      (sum_1d(wt_1d(1:dL), dL) + Eps_len)  ! Eps_len avoids division by 0.
    enddo ; enddo
  else
    write (mesg,*) " unknown sampling method: ",method
    call MOM_error(FATAL, "downsample_field_2d: "//trim(mesg)//" "//trim(diag%debug_str))
  endif

end procedure downsample_field_2d
module procedure sum_1d
  integer :: i, sz_2
  if (sz == 2) then      ! The order of arithmetic does not matter.
    sum = field(1) + field(2)
  elseif (sz == 3) then  ! Use simpler code that has the same order of sums as the general case.
    sum = field(2) + (field(1) + field(3))
  else
    ! This is a copy of the general code from symmetric_sum_1d in MOM_array_transform
    sz_2 = sz / 2 ! Note that for an odd number sz_2 is rounded down.
    sum = 0.0
    if (2*sz_2 < sz) sum = field(sz_2+1)
    ! Add pairs of values, working from the inside out.
    do i=sz_2,1,-1
      sum = sum + (field(i) + field(sz+1-i))
    enddo
  endif
end procedure sum_1d
module procedure square_sum
  integer :: i, j
  logical :: simple_sum
  simple_sum = .true. ; if (present(naive_sum)) simple_sum = naive_sum

  if (sz == 1) then
    sum = field(1,1)
  elseif (simple_sum) then
    ! This non-rotationally symmetric sum is here to reproduce previous results.
    sum = 0.0
    do j=1,sz ; do i=1,sz ; sum = sum + field(i,j) ; enddo ; enddo
  elseif (sz == 2) then
    ! This copy of code from symmetric_sum may facilitate inlining in a common case.
    sum = (field(1,1) + field(2,2)) + (field(2,1) + field(1,2))
  elseif (sz == 3) then
    ! This copy of code from symmetric_sum may facilitate inlining in a common case.
    sum = (field(2,2) + ((field(1,2) + field(3,2)) + (field(2,1) + field(2,3)))) + &
          ((field(1,1) + field(3,3)) + (field(3,1) + field(1,3)))
  else
    sum = symmetric_sum(field(1:sz,1:sz))
  endif

end procedure square_sum
module procedure downsample_mask_2d
  real    :: tot_non_zero    ! The sum of mask values in the down-scaled cell or face [nondim]
  character(len=8) :: method_str
  integer :: i, j, i_dn, j_dn
  integer :: ii, jj  ! The index locations on the full grid that contribute to the averages.
  integer :: i0_off, j0_off  ! The starting point offsets between full array and reduced array
  allocate(mask_dsamp(isd_d:ied_d, jsd_d:jed_d), source=0.0)

  i0_off = ((isc_o-1) - dL*isc_d)
  j0_off = ((jsc_o-1) - dL*jsc_d)
  if ((method == MMM) .or. (method == MMP) .or. (method == MMS) .or. (method == SSS)) then
    ! This applies at tracer points.
    do j=jsc_d,jec_d ; do i=isc_d,iec_d
      tot_non_zero = 0.0
      do j_dn=1,dL ; do i_dn=1,dL
        ii = i_dn + (dL*i + i0_off)
        jj = j_dn + (dL*j + j0_off)
        tot_non_zero = tot_non_zero + abs(mask_in(ii,jj))
      enddo ; enddo
      if (tot_non_zero > 0.0) mask_dsamp(i,j) = 1.0
    enddo ; enddo
  elseif ((method == PMM) .or. (method == PSP) .or. (method == PMP) .or. (method == PSS)) then
    ! This applies at u-velocity points.
    do j=jsc_d,jec_d ; do I=isc_d,iec_d
      tot_non_zero = 0.0
      II = (dL*I + I0_off) + (dL-1)
      do j_dn=1,dL
        jj = j_dn + (dL*j + j0_off)
        tot_non_zero = tot_non_zero + abs(mask_in(II,jj))
      enddo
      if (tot_non_zero > 0.0) mask_dsamp(I,j) = 1.0
    enddo ; enddo
  elseif ((method == MPM) .or. (method == SPP) .or. (method == MPP) .or. (method == SPS)) then
    ! This applies at v-velocity points.
    do J=jsc_d,jec_d ; do i=isc_d,iec_d
      tot_non_zero = 0.0
      JJ = (dL*J + J0_off) + (dL-1)
      do i_dn=1,dL
        ii = i_dn + (dL*i + i0_off)
        tot_non_zero = tot_non_zero + abs(mask_in(ii,JJ))
      enddo
      if (tot_non_zero > 0.0) mask_dsamp(i,J) = 1.0
    enddo ; enddo
  elseif ((method == PPP) .or. (method == PPM)) then
    ! This applies at corner (vorticity) points.
    do j=jsc_d,jec_d ; do I=isc_d,iec_d
      II = (dL*I + I0_off) + (dL-1)
      JJ = (dL*J + J0_off) + (dL-1)
      if (abs(mask_in(II,JJ)) > 0.0) mask_dsamp(I,J) = 1.0
    enddo ; enddo
  else
    write(method_str, '(I0)') method
    call MOM_error(FATAL, "downsample_mask_2d: unknown sampling method "//trim(method_str))
  endif

end procedure downsample_mask_2d
module procedure downsample_mask_3d
  real    :: tot_non_zero    ! The sum of mask values in the down-scaled cell or face [nondim]
  character(len=8) :: method_str
  integer :: i, j, i_dn, j_dn, k, ks, ke
  integer :: ii, jj  ! The index locations on the full grid that contribute to the averages.
  integer :: i0_off, j0_off  ! The starting point offsets between full array and reduced array
  ks = lbound(mask_in, 3) ; ke = ubound(mask_in, 3)
  allocate(mask_dsamp(isd_d:ied_d, jsd_d:jed_d, ks:ke), source=0.0)

  i0_off = ((isc_o-1) - dL*isc_d)
  j0_off = ((jsc_o-1) - dL*jsc_d)
  if ((method == MMM) .or. (method == MMP) .or. (method == MMS) .or. (method == SSS)) then
    ! This applies at tracer points.
    do k=ks,ke ; do j=jsc_d,jec_d ; do i=isc_d,iec_d
      tot_non_zero = 0.0
      do j_dn=1,dL ; do i_dn=1,dL
        ii = i_dn + (dL*i + i0_off)
        jj = j_dn + (dL*j + j0_off)
        tot_non_zero = tot_non_zero + abs(mask_in(ii,jj,k))
      enddo ; enddo
      if (tot_non_zero > 0.0) mask_dsamp(i,j,k) = 1.0
    enddo ; enddo ; enddo
  elseif ((method == PMM) .or. (method == PSP) .or. (method == PMP) .or. (method == PSS)) then
    ! This applies at u-velocity points.
    do k=ks,ke ; do j=jsc_d,jec_d ; do I=isc_d,iec_d
      tot_non_zero = 0.0
      II = (dL*I + I0_off) + (dL-1)
      do j_dn=1,dL
        jj = j_dn + (dL*j + j0_off)
        tot_non_zero = tot_non_zero + abs(mask_in(II,jj,k))
      enddo
      if (tot_non_zero > 0.0) mask_dsamp(I,j,k) = 1.0
    enddo ; enddo ; enddo
  elseif ((method == MPM) .or. (method == SPP) .or. (method == MPP) .or. (method == SPS)) then
    ! This applies at v-velocity points.
    do k=ks,ke ; do J=jsc_d,jec_d ; do i=isc_d,iec_d
      tot_non_zero = 0.0
      JJ = (dL*J + J0_off) + (dL-1)
      do i_dn=1,dL
        ii = i_dn + (dL*i + i0_off)
        tot_non_zero = tot_non_zero + abs(mask_in(ii,JJ,k))
      enddo
      if (tot_non_zero > 0.0) mask_dsamp(i,J,k) = 1.0
    enddo ; enddo ; enddo
  elseif ((method == PPP) .or. (method == PPM)) then
    ! This applies at corner (vorticity) points.
    do k=ks,ke ; do j=jsc_d,jec_d ; do I=isc_d,iec_d
      II = (dL*I + I0_off) + (dL-1)
      JJ = (dL*J + J0_off) + (dL-1)
      if (abs(mask_in(II,JJ,k)) > 0.0) mask_dsamp(I,J,k) = 1.0
    enddo ; enddo ; enddo
  else
    write(method_str, '(I0)') method
    call MOM_error(FATAL, "downsample_mask_3d: unknown sampling method "//trim(method_str))
  endif

end procedure downsample_mask_3d
module procedure found_in_diagtable
  integer :: handle ! Integer handle returned from diag_manager
  handle = register_static_field_infra('ocean_model', varName, diag%axesT1%handles)

  found_in_diagtable = (handle>0)

end procedure found_in_diagtable
module procedure MOM_diag_send_complete
  call diag_send_complete_infra()
end procedure MOM_diag_send_complete
end submodule MOM_diag_mediator_s
