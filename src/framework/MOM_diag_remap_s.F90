submodule (MOM_diag_remap) MOM_diag_remap_s
#include "MOM_memory.h"
  implicit none
contains
module procedure diag_remap_init
  remap_cs%diag_module_suffix = trim(extractWord(coord_tuple, 1))
  remap_cs%diag_coord_name = trim(extractWord(coord_tuple, 2))
  remap_cs%vertical_coord_name = trim(extractWord(coord_tuple, 3))
  remap_cs%vertical_coord = coordinateMode(remap_cs%vertical_coord_name)
  remap_cs%Z_based_coord = .false.
  if (.not.(GV%Boussinesq .or. GV%semi_Boussinesq) .and. &
      ((remap_cs%vertical_coord == coordinateMode('ZSTAR')) .or. &
       (remap_cs%vertical_coord == coordinateMode('SIGMA')) .or. &
       (remap_cs%vertical_coord == coordinateMode('RHO'))) ) &
    remap_cs%Z_based_coord = .true.

  remap_cs%configured = .false.
  remap_cs%initialized = .false.
  remap_cs%used = .false.
  remap_cs%om4_remap_via_sub_cells = om4_remap_via_sub_cells
  remap_cs%answer_date = answer_date
  remap_cs%nz = 0

end procedure diag_remap_init
module procedure diag_remap_end
  if (allocated(remap_cs%h)) deallocate(remap_cs%h)

  remap_cs%configured = .false.
  remap_cs%initialized = .false.
  remap_cs%used = .false.
  remap_cs%nz = 0

end procedure diag_remap_end
module procedure diag_remap_diag_registration_closed
  if (.not. remap_cs%used) then
    call diag_remap_end(remap_cs)
    call end_regridding(remap_cs%regrid_cs)
  endif

end procedure diag_remap_diag_registration_closed
module procedure diag_remap_set_active
  remap_cs%used = .true.

end procedure diag_remap_set_active
module procedure diag_remap_configure_axes
  character(len=40)  :: mod  = "MOM_diag_remap" ! This module's name.
  character(len=8)   :: units
  character(len=34)  :: longname
  real, allocatable, dimension(:) :: &
    interfaces, & ! Numerical values for interface vertical coordinates, in unscaled units
                  ! that might be [m], [kg m-3] or [nondim], depending on the coordinate.
    layers        ! Numerical values for layer vertical coordinates, in unscaled units
                  ! that might be [m], [kg m-3] or [nondim], depending on the coordinate.

  call initialize_regridding(remap_cs%regrid_cs, G, GV, US, GV%max_depth, param_file, mod, &
           trim(remap_cs%vertical_coord_name), "DIAG_COORD", trim(remap_cs%diag_coord_name))
  call set_regrid_params(remap_cs%regrid_cs, min_thickness=0., integrate_downward_for_e=.false.)

  remap_cs%nz = get_regrid_size(remap_cs%regrid_cs)

  if (remap_cs%vertical_coord == coordinateMode('SIGMA')) then
    units = 'nondim'
    longname = 'Fraction'
  elseif (remap_cs%vertical_coord == coordinateMode('RHO')) then
    units = 'kg m-3'
    longname = 'Target Potential Density'
  else
    units = 'meters'
    longname = 'Depth'
  endif

  ! Make axes objects
  allocate(interfaces(remap_cs%nz+1))
  allocate(layers(remap_cs%nz))

  interfaces(:) = getCoordinateInterfaces(remap_cs%regrid_cs, undo_scaling=.true.)
  layers(:) = 0.5 * ( interfaces(1:remap_cs%nz) + interfaces(2:remap_cs%nz+1) )

  remap_cs%interface_axes_id = MOM_diag_axis_init(lowercase(trim(remap_cs%diag_coord_name))//'_i', &
                                              interfaces, trim(units), 'z', &
                                              trim(longname)//' at interface', direction=-1)
  remap_cs%layer_axes_id = MOM_diag_axis_init(lowercase(trim(remap_cs%diag_coord_name))//'_l', &
                                          layers, trim(units), 'z', &
                                          trim(longname)//' at cell center', direction=-1, &
                                          edges=remap_cs%interface_axes_id)

  ! Axes have now been configured.
  remap_cs%configured = .true.

  deallocate(interfaces)
  deallocate(layers)

end procedure diag_remap_configure_axes
module procedure diag_remap_get_axes_info
  nz = remap_cs%nz
  id_layer = remap_cs%layer_axes_id
  id_interface = remap_cs%interface_axes_id

end procedure diag_remap_get_axes_info
module procedure diag_remap_axes_configured
  diag_remap_axes_configured = remap_cs%configured

end procedure diag_remap_axes_configured
module procedure diag_remap_update
  real, dimension(remap_cs%nz + 1) :: zInterfaces ! Interface positions [H ~> m or kg m-2] or [Z ~> m]
  real :: h_neglect, h_neglect_edge ! Negligible thicknesses [H ~> m or kg m-2] or [Z ~> m]
  real :: bottom_depth(SZI_(G),SZJ_(G)) ! The depth of the bathymetry in [H ~> m or kg m-2] or [Z ~> m]
  real :: h_tot(SZI_(G),SZJ_(G))        ! The total thickness of the water column [H ~> m or kg m-2] or [Z ~> m]
  real :: Z_unit_scale   ! A conversion factor from Z-units the internal work units in this routine,
  integer :: i, j, k, is, ie, js, je, nz
  if (.not. remap_cs%configured) return

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec

  ! Set the bottom depth and negligible thicknesses used in the coordinate remapping in the right units.
  if (remap_cs%Z_based_coord) then
    h_neglect = set_dz_neglect(GV, US, remap_cs%answer_date, h_neglect_edge)
    Z_unit_scale = 1.0
    do j=js-1,je+1 ; do i=is-1,ie+1
      bottom_depth(i,j) = G%bathyT(i,j) + G%Z_ref
    enddo ; enddo
  else
    h_neglect = set_h_neglect(GV, remap_cs%answer_date, h_neglect_edge)
    Z_unit_scale = GV%Z_to_H  ! This branch is not used in fully non-Boussinesq mode.
    do j=js-1,je+1 ; do i=is-1,ie+1
      bottom_depth(i,j) = GV%Z_to_H * (G%bathyT(i,j) + G%Z_ref)
    enddo ; enddo
  endif

  if (.not. remap_cs%initialized) then
    ! Initialize remapping and regridding on the first call
    call initialize_remapping(remap_cs%remap_cs, 'PPM_IH4', boundary_extrapolation=.false., &
                              om4_remap_via_sub_cells=remap_cs%om4_remap_via_sub_cells, &
                              answer_date=remap_cs%answer_date, &
                              h_neglect=h_neglect, h_neglect_edge=h_neglect_edge)
    remap_cs%initialized = .true.
  endif

  ! Calculate the total thickness of the water column, if it is needed,
  if ((remap_cs%vertical_coord == coordinateMode('ZSTAR')) .or. &
      (remap_cs%vertical_coord == coordinateMode('SIGMA'))) then
    if (remap_CS%answer_date >= 20240201) then
      ! Avoid using sum to have a specific order for the vertical sums.
      ! For some compilers, the explicit expression gives the same answers as the sum function.
      h_tot(:,:) = 0.0
      do k=1,GV%ke ; do j=js-1,je+1 ; do i=is-1,ie+1
        h_tot(i,j) = h_tot(i,j) + h(i,j,k)
      enddo ; enddo ; enddo
    else
      do j=js-1,je+1 ; do i=is-1,ie+1
        h_tot(i,j) = sum(h(i,j,:))
      enddo ; enddo
    endif
  endif

  ! Calculate remapping thicknesses for different target grids based on
  ! nominal/target interface locations. This happens for every call on the
  ! assumption that h, T, S has changed.
  h_target(:,:,:) = 0.0

  nz = remap_cs%nz
  if (remap_cs%vertical_coord == coordinateMode('ZSTAR')) then
    do j=js-1,je+1 ; do i=is-1,ie+1 ; if (G%mask2dT(i,j) > 0.0) then
      ! This function call can work with the last 4 arguments all in units of [Z ~> m] or [H ~> kg m-2].
      call build_zstar_column(get_zlike_CS(remap_cs%regrid_cs), &
                              bottom_depth(i,j), h_tot(i,j), zInterfaces, zScale=Z_unit_scale)
      do k=1,nz ; h_target(i,j,k) = zInterfaces(K) - zInterfaces(K+1) ; enddo
    endif ; enddo ; enddo
  elseif (remap_cs%vertical_coord == coordinateMode('SIGMA')) then
    do j=js-1, je+1 ; do i=is-1,ie+1 ; if (G%mask2dT(i,j) > 0.0) then
      ! This function call can work with the last 3 arguments all in units of [Z ~> m] or [H ~> kg m-2].
      call build_sigma_column(get_sigma_CS(remap_cs%regrid_cs), &
                              bottom_depth(i,j), h_tot(i,j), zInterfaces)
      do k=1,nz ; h_target(i,j,k) = zInterfaces(K) - zInterfaces(K+1) ; enddo
    endif ; enddo ; enddo
  elseif (remap_cs%vertical_coord == coordinateMode('RHO')) then
    do j=js-1,je+1 ; do i=is-1,ie+1 ; if (G%mask2dT(i,j) > 0.0) then
      ! This function call can work with 5 arguments in units of [Z ~> m] or [H ~> kg m-2].
      call build_rho_column(get_rho_CS(remap_cs%regrid_cs), GV%ke, &
                            bottom_depth(i,j), h(i,j,:), T(i,j,:), S(i,j,:), &
                            eqn_of_state, zInterfaces, h_neglect=h_neglect, h_neglect_edge=h_neglect_edge)
      do k=1,nz ; h_target(i,j,k) = zInterfaces(K) - zInterfaces(K+1) ; enddo
    endif ; enddo ; enddo
  elseif (remap_cs%vertical_coord == coordinateMode('HYCOM1')) then
    call MOM_error(FATAL,"diag_remap_update: HYCOM1 coordinate not coded for diagnostics yet!")
!    do j=js-1,je+1 ; do i=is-1,ie+1 ; if (G%mask2dT(i,j) > 0.0) then
!      call build_hycom1_column(remap_cs%regrid_cs, nz, &
!                           bottom_depth(i,j), h_tot(i,j), zInterfaces)
!      do k=1,nz ; h_target(i,j,k) = zInterfaces(K) - zInterfaces(K+1) ; enddo
!    endif ; enddo ; enddo
  endif

end procedure diag_remap_update
module procedure diag_remap_do_remap
  integer :: isdf, jsdf !< The starting i- and j-indices in memory for field
  call assert(remap_cs%initialized, 'diag_remap_do_remap: remap_cs not initialized.')
  call assert(size(field, 3) == size(h, 3), &
              'diag_remap_do_remap: Remap field and thickness z-axes do not match.')

  isdf = G%isd ; if (staggered_in_x) Isdf = G%IsdB
  jsdf = G%jsd ; if (staggered_in_y) Jsdf = G%JsdB

  if (associated(mask)) then
    call do_remap(remap_cs, G, GV, US, isdf, jsdf, h, OBC_u, OBC_v, staggered_in_x, staggered_in_y, &
                  field, remapped_field, mask(:,:,1))
  else
    call do_remap(remap_cs, G, GV, US, isdf, jsdf, h, OBC_u, OBC_v, staggered_in_x, staggered_in_y, &
                  field, remapped_field)
  endif

end procedure diag_remap_do_remap
module procedure do_remap
  real, dimension(remap_cs%nz) :: h_dest ! Destination thicknesses [H ~> m or kg m-2] or [Z ~> m]
  real, dimension(size(h,3)) :: h_src    ! A column of source thicknesses [H ~> m or kg m-2] or [Z ~> m]
  integer :: nz_src, nz_dest        ! The number of layers on the native and remapped grids
  integer :: i, j                   ! Grid index
  nz_src = size(field,3)
  nz_dest = remap_cs%nz
  remapped_field(:,:,:) = 0.

  if (staggered_in_x .and. .not. staggered_in_y) then
    ! U-points
    if (present(mask)) then
      do j=G%jsc,G%jec ; do I=G%IscB,G%IecB ; if (mask(I,j) > 0.) then
        if (OBC_u(I,j) == 0) then    ! This is not an OBC face.
          h_src(:) = 0.5*(h(i,j,:) + h(i+1,j,:))
          h_dest(:) = 0.5*(remap_cs%h(i,j,:) + remap_cs%h(i+1,j,:))
        elseif (OBC_u(I,j) < 0) then ! This is a western OBC face.
          h_src(:) = h(i+1,j,:)
          h_dest(:) = remap_cs%h(i+1,j,:)
        else    ! (OBC_u(I,j) > 0)   ! This is a eastern OBC face.
          h_src(:) = h(i,j,:)
          h_dest(:) = remap_cs%h(i,j,:)
        endif
        call remapping_core_h(remap_cs%remap_cs, nz_src, h_src(:), field(I,j,:), &
                              nz_dest, h_dest(:), remapped_field(I,j,:))
      endif ; enddo ; enddo
    else
      do j=G%jsc,G%jec ; do I=G%IscB,G%IecB
        if (OBC_u(I,j) == 0) then    ! This is not an OBC face.
          h_src(:) = 0.5*(h(i,j,:) + h(i+1,j,:))
          h_dest(:) = 0.5*(remap_cs%h(i,j,:) + remap_cs%h(i+1,j,:))
        elseif (OBC_u(I,j) < 0) then ! This is a western OBC face.
          h_src(:) = h(i+1,j,:)
          h_dest(:) = remap_cs%h(i+1,j,:)
        else    ! (OBC_u(I,j) > 0)   ! This is a eastern OBC face.
          h_src(:) = h(i,j,:)
          h_dest(:) = remap_cs%h(i,j,:)
        endif
        call remapping_core_h(remap_cs%remap_cs, nz_src, h_src(:), field(I,j,:), &
                              nz_dest, h_dest(:), remapped_field(I,j,:))
      enddo ; enddo
    endif
  elseif (staggered_in_y .and. .not. staggered_in_x) then
    ! V-points
    if (present(mask)) then
      do J=G%jscB,G%jecB ; do i=G%isc,G%iec ; if (mask(i,j) > 0.) then
        if (OBC_v(i,J) == 0) then    ! This is not an OBC face.
          h_src(:) = 0.5*(h(i,j,:) + h(i,j+1,:))
          h_dest(:) = 0.5*(remap_cs%h(i,j,:) + remap_cs%h(i,j+1,:))
        elseif (OBC_v(i,J) < 0) then ! This is a southern OBC face
          h_src(:) = h(i,j+1,:)
          h_dest(:) = remap_cs%h(i,j+1,:)
        else    ! (OBC_v(i,J) > 0)   ! This is a northern OBC face
          h_src(:) = h(i,j,:)
          h_dest(:) = remap_cs%h(i,j,:)
        endif
        call remapping_core_h(remap_cs%remap_cs, nz_src, h_src(:), field(i,J,:), &
                              nz_dest, h_dest(:), remapped_field(i,J,:))
      endif ; enddo ; enddo
    else
      do J=G%jscB,G%jecB ; do i=G%isc,G%iec
        if (OBC_v(i,J) == 0) then    ! This is not an OBC face.
          h_src(:) = 0.5*(h(i,j,:) + h(i,j+1,:))
          h_dest(:) = 0.5*(remap_cs%h(i,j,:) + remap_cs%h(i,j+1,:))
        elseif (OBC_v(i,J) < 0) then ! This is a southern OBC face
          h_src(:) = h(i,j+1,:)
          h_dest(:) = remap_cs%h(i,j+1,:)
        else    ! (OBC_v(i,J) > 0)   ! This is a northern OBC face
          h_src(:) = h(i,j,:)
          h_dest(:) = remap_cs%h(i,j,:)
        endif
        call remapping_core_h(remap_cs%remap_cs, nz_src, h_src(:), field(i,J,:), &
                              nz_dest, h_dest(:), remapped_field(i,J,:))
      enddo ; enddo
    endif
  elseif ((.not. staggered_in_x) .and. (.not. staggered_in_y)) then
    ! H-points
    if (present(mask)) then
      do j=G%jsc,G%jec ; do i=G%isc,G%iec ; if (mask(i,j) > 0.) then
        call remapping_core_h(remap_cs%remap_cs, nz_src, h(i,j,:), field(i,j,:), &
                              nz_dest, remap_cs%h(i,j,:), remapped_field(i,j,:))
      endif ; enddo ; enddo
    else
      do j=G%jsc,G%jec ; do i=G%isc,G%iec
        call remapping_core_h(remap_cs%remap_cs, nz_src, h(i,j,:), field(i,j,:), &
                              nz_dest, remap_cs%h(i,j,:), remapped_field(i,j,:))
      enddo ; enddo
    endif
  else
    call assert(.false., 'diag_remap_do_remap: Unsupported axis combination')
  endif

end procedure do_remap
module procedure diag_remap_calc_hmask
  real, dimension(remap_cs%nz) :: h_dest ! Destination thicknesses [H ~> m or kg m-2] or [Z ~> m]
  integer :: i, j, k
  logical :: mask_vanished_layers
  real :: h_tot      ! Sum of all thicknesses [H ~> m or kg m-2] or [Z ~> m]
  real :: h_err      ! An estimate of a negligible thickness [H ~> m or kg m-2] or [Z ~> m]
  call assert(remap_cs%initialized, 'diag_remap_calc_hmask: remap_cs not initialized.')

  ! Only z*-like diagnostic coordinates should have a 3d mask
  mask_vanished_layers = (remap_cs%vertical_coord == coordinateMode('ZSTAR'))
  mask(:,:,:) = 0.

  do j=G%jsc-1,G%jec+1 ; do i=G%isc-1,G%iec+1
    if (G%mask2dT(i,j)>0.) then
      if (mask_vanished_layers) then
        h_dest(:) = remap_cs%h(i,j,:)
        h_tot = 0.
        h_err = 0.
        do k=1, remap_cs%nz
          h_tot = h_tot + h_dest(k)
          ! This is an overestimate of how thick a vanished layer might be, that
          ! appears due to round-off.
          h_err = h_err + epsilon(h_tot) * h_tot
          ! Mask out vanished layers
          if (h_dest(k)<=8.*h_err) then
            mask(i,j,k) = 0.
          else
            mask(i,j,k) = 1.
          endif
        enddo
      else ! all layers might contain data
        mask(i,j,:) = 1.
      endif
    endif
  enddo ; enddo

end procedure diag_remap_calc_hmask
module procedure vertically_reintegrate_diag_field
  integer :: isdf, jsdf !< The starting i- and j-indices in memory for field
  call assert(remap_cs%initialized, 'vertically_reintegrate_diag_field: remap_cs not initialized.')
  call assert(size(field, 3) == size(h, 3), &
              'vertically_reintegrate_diag_field: Remap field and thickness z-axes do not match.')

  isdf = G%isd ; if (staggered_in_x) Isdf = G%IsdB
  jsdf = G%jsd ; if (staggered_in_y) Jsdf = G%JsdB

  if (associated(mask)) then
    call vertically_reintegrate_field(remap_cs, G, isdf, jsdf, h, h_target, OBC_u, OBC_v, &
                                      staggered_in_x, staggered_in_y, field, reintegrated_field, mask(:,:,1))
  else
    call vertically_reintegrate_field(remap_cs, G, isdf, jsdf, h, h_target, OBC_u, OBC_v, &
                                      staggered_in_x, staggered_in_y, field, reintegrated_field)
  endif

end procedure vertically_reintegrate_diag_field
module procedure vertically_reintegrate_field
  real, dimension(remap_cs%nz) :: h_dest ! Destination thicknesses [H ~> m or kg m-2] or [Z ~> m]
  real, dimension(size(h,3)) :: h_src    ! A column of source thicknesses [H ~> m or kg m-2] or [Z ~> m]
  integer :: nz_src, nz_dest        ! The number of layers on the native and remapped grids
  integer :: i, j                   ! Grid index
  nz_src = size(field,3)
  nz_dest = remap_cs%nz
  reintegrated_field(:,:,:) = 0.

  if (staggered_in_x .and. .not. staggered_in_y) then
    ! U-points
    if (present(mask)) then
      do j=G%jsc,G%jec ; do I=G%IscB,G%IecB ; if (mask(I,j) > 0.0) then
        if (OBC_u(I,j) == 0) then    ! This is not an OBC face.
          h_src(:) = 0.5*(h(i,j,:) + h(i+1,j,:))
          h_dest(:) = 0.5*(h_target(i,j,:) + h_target(i+1,j,:))
        elseif (OBC_u(I,j) < 0) then ! This is a western OBC face
          h_src(:) = h(i+1,j,:)
          h_dest(:) = h_target(i+1,j,:)
        else    ! (OBC_u(I,j) > 0)   ! This is an eastern OBC face
          h_src(:) = h(i,j,:)
          h_dest(:) = h_target(i,j,:)
        endif
        call reintegrate_column(nz_src, h_src, field(I,j,:), &
                                nz_dest, h_dest, reintegrated_field(I,j,:))
      endif ; enddo ; enddo
    else
      do j=G%jsc,G%jec ; do I=G%IscB,G%IecB
        if (OBC_u(I,j) == 0) then    ! This is not an OBC face.
          h_src(:) = 0.5*(h(i,j,:) + h(i+1,j,:))
          h_dest(:) = 0.5*(h_target(i,j,:) + h_target(i+1,j,:))
        elseif (OBC_u(I,j) < 0) then ! This is a western OBC face
          h_src(:) = h(i+1,j,:)
          h_dest(:) = h_target(i+1,j,:)
        else    ! (OBC_u(I,j) > 0)   ! This is an eastern OBC face
          h_src(:) = h(i,j,:)
          h_dest(:) = h_target(i,j,:)
        endif
        call reintegrate_column(nz_src, h_src, field(I,j,:), &
                                nz_dest, h_dest, reintegrated_field(I,j,:))
      enddo ; enddo
    endif
  elseif (staggered_in_y .and. .not. staggered_in_x) then
    ! V-points
    if (present(mask)) then
      do J=G%jscB,G%jecB ; do i=G%isc,G%iec ; if (mask(i,J) > 0.0) then
        if (OBC_v(i,J) == 0) then    ! This is not an OBC face.
          h_src(:) = 0.5*(h(i,j,:) + h(i,j+1,:))
          h_dest(:) = 0.5*(h_target(i,j,:) + h_target(i,j+1,:))
        elseif (OBC_v(i,J) < 0) then ! This is a southern OBC face
          h_src(:) = h(i,j+1,:)
          h_dest(:) = h_target(i,j+1,:)
        else    ! (OBC_v(i,J) > 0)   ! This is a northern OBC face
          h_src(:) = h(i,j,:)
          h_dest(:) = h_target(i,j,:)
        endif
        call reintegrate_column(nz_src, h_src, field(i,J,:), &
                                nz_dest, h_dest, reintegrated_field(i,J,:))
      endif ; enddo ; enddo
    else
      do J=G%jscB,G%jecB ; do i=G%isc,G%iec
        if (OBC_v(i,J) == 0) then    ! This is not an OBC face.
          h_src(:) = 0.5*(h(i,j,:) + h(i,j+1,:))
          h_dest(:) = 0.5*(h_target(i,j,:) + h_target(i,j+1,:))
        elseif (OBC_v(i,J) < 0) then ! This is a southern OBC face
          h_src(:) = h(i,j+1,:)
          h_dest(:) = h_target(i,j+1,:)
        else    ! (OBC_v(i,J) > 0)   ! This is a northern OBC face
          h_src(:) = h(i,j,:)
          h_dest(:) = h_target(i,j,:)
        endif
        call reintegrate_column(nz_src, h_src, field(i,J,:), &
                                nz_dest, h_dest, reintegrated_field(i,J,:))
      enddo ; enddo
    endif
  elseif ((.not. staggered_in_x) .and. (.not. staggered_in_y)) then
    ! H-points
    if (present(mask)) then
      do j=G%jsc,G%jec ; do i=G%isc,G%iec ; if (mask(i,J) > 0.0) then
        call reintegrate_column(nz_src, h(i,j,:), field(i,j,:), &
                                nz_dest, h_target(i,j,:), reintegrated_field(i,j,:))
      endif ; enddo ; enddo
    else
      do j=G%jsc,G%jec ; do i=G%isc,G%iec
        call reintegrate_column(nz_src, h(i,j,:), field(i,j,:), &
                                nz_dest, h_target(i,j,:), reintegrated_field(i,j,:))
      enddo ; enddo
    endif
  else
    call assert(.false., 'vertically_reintegrate_diag_field: Q point remapping is not coded yet.')
  endif

end procedure vertically_reintegrate_field
module procedure vertically_interpolate_diag_field
  integer :: isdf, jsdf !< The starting i- and j-indices in memory for field
  call assert(remap_cs%initialized, 'vertically_interpolate_diag_field: remap_cs not initialized.')
  call assert(size(field, 3) == size(h, 3)+1, &
              'vertically_interpolate_diag_field: Remap field and thickness z-axes do not match.')

  isdf = G%isd ; if (staggered_in_x) Isdf = G%IsdB
  jsdf = G%jsd ; if (staggered_in_y) Jsdf = G%JsdB

  if (associated(mask)) then
    call vertically_interpolate_field(remap_cs, G, isdf, jsdf, h, OBC_u, OBC_v, staggered_in_x, staggered_in_y, &
                                      field, interpolated_field, mask(:,:,1))
  else
    call vertically_interpolate_field(remap_cs, G, isdf, jsdf, h, OBC_u, OBC_v, staggered_in_x, staggered_in_y, &
                                      field, interpolated_field)
  endif

end procedure vertically_interpolate_diag_field
module procedure vertically_interpolate_field
  real, dimension(remap_cs%nz) :: h_dest ! Destination thicknesses [H ~> m or kg m-2] or [Z ~> m]
  real, dimension(size(h,3)) :: h_src    ! A column of source thicknesses [H ~> m or kg m-2] or [Z ~> m]
  integer :: nz_src, nz_dest        ! The number of layers on the native and remapped grids
  integer :: i, j                   !< Grid index
  interpolated_field(:,:,:) = 0.

  nz_src = size(h,3)
  nz_dest = remap_cs%nz

  if (staggered_in_x .and. .not. staggered_in_y) then
    ! U-points
    if (present(mask)) then
      do j=G%jsc,G%jec ; do I=G%IscB,G%IecB ; if (mask(I,j) > 0.0) then
        if (OBC_u(I,j) == 0) then    ! This is not an OBC face.
          h_src(:) = 0.5*(h(i,j,:) + h(i+1,j,:))
          h_dest(:) = 0.5*(remap_cs%h(i,j,:) + remap_cs%h(i+1,j,:))
        elseif (OBC_u(I,j) < 0) then ! This is a western OBC face.
          h_src(:) = h(i+1,j,:)
          h_dest(:) = remap_cs%h(i+1,j,:)
        else    ! (OBC_u(I,j) > 0)   ! This is a eastern OBC face.
          h_src(:) = h(i,j,:)
          h_dest(:) = remap_cs%h(i,j,:)
        endif
        call interpolate_column(nz_src, h_src, field(I,j,:), &
                                nz_dest, h_dest, interpolated_field(I,j,:), .true.)
      endif ; enddo ; enddo
    else
      do j=G%jsc,G%jec ; do I=G%IscB,G%IecB
        if (OBC_u(I,j) == 0) then    ! This is not an OBC face.
          h_src(:) = 0.5*(h(i,j,:) + h(i+1,j,:))
          h_dest(:) = 0.5*(remap_cs%h(i,j,:) + remap_cs%h(i+1,j,:))
        elseif (OBC_u(I,j) < 0) then ! This is a western OBC face.
          h_src(:) = h(i+1,j,:)
          h_dest(:) = remap_cs%h(i+1,j,:)
        else    ! (OBC_u(I,j) > 0)   ! This is a eastern OBC face.
          h_src(:) = h(i,j,:)
          h_dest(:) = remap_cs%h(i,j,:)
        endif
        call interpolate_column(nz_src, h_src, field(I,j,:), &
                                nz_dest, h_dest, interpolated_field(I,j,:), .true.)
      enddo ; enddo
    endif
  elseif (staggered_in_y .and. .not. staggered_in_x) then
    ! V-points
    if (present(mask)) then
      do J=G%jscB,G%jecB ; do i=G%isc,G%iec ; if (mask(I,j) > 0.0) then
        if (OBC_v(i,J) == 0) then    ! This is not an OBC face.
          h_src(:) = 0.5*(h(i,j,:) + h(i,j+1,:))
          h_dest(:) = 0.5*(remap_cs%h(i,j,:) + remap_cs%h(i,j+1,:))
        elseif (OBC_v(i,J) < 0) then ! This is a southern OBC face
          h_src(:) = h(i,j+1,:)
          h_dest(:) = remap_cs%h(i,j+1,:)
        else    ! (OBC_v(i,J) > 0)   ! This is a northern OBC face
          h_src(:) = h(i,j,:)
          h_dest(:) = remap_cs%h(i,j,:)
        endif
        call interpolate_column(nz_src, h_src, field(i,J,:), &
                                nz_dest, h_dest, interpolated_field(i,J,:), .true.)
      endif ; enddo ; enddo
    else
      do J=G%jscB,G%jecB ; do i=G%isc,G%iec
        if (OBC_v(i,J) == 0) then    ! This is not an OBC face.
          h_src(:) = 0.5*(h(i,j,:) + h(i,j+1,:))
          h_dest(:) = 0.5*(remap_cs%h(i,j,:) + remap_cs%h(i,j+1,:))
        elseif (OBC_v(i,J) < 0) then ! This is a southern OBC face
          h_src(:) = h(i,j+1,:)
          h_dest(:) = remap_cs%h(i,j+1,:)
        else    ! (OBC_v(i,J) > 0)   ! This is a northern OBC face
          h_src(:) = h(i,j,:)
          h_dest(:) = remap_cs%h(i,j,:)
        endif
        call interpolate_column(nz_src, h_src, field(i,J,:), &
                                nz_dest, h_dest, interpolated_field(i,J,:), .true.)
      enddo ; enddo
    endif
  elseif ((.not. staggered_in_x) .and. (.not. staggered_in_y)) then
    ! H-points
    if (present(mask)) then
      do j=G%jsc,G%jec ; do i=G%isc,G%iec ; if (mask(i,j) > 0.0) then
        call interpolate_column(nz_src, h(i,j,:), field(i,j,:), &
                                nz_dest, remap_cs%h(i,j,:), interpolated_field(i,j,:), .true.)
      endif ; enddo ; enddo
    else
      do j=G%jsc,G%jec ; do i=G%isc,G%iec
        call interpolate_column(nz_src, h(i,j,:), field(i,j,:), &
                                nz_dest, remap_cs%h(i,j,:), interpolated_field(i,j,:), .true.)
      enddo ; enddo
    endif
  else
    call assert(.false., 'vertically_interpolate_diag_field: Q point remapping is not coded yet.')
  endif

end procedure vertically_interpolate_field
module procedure horizontally_average_diag_field
  integer :: isdf, jsdf !< The starting i- and j-indices in memory for field
  isdf = G%isd ; if (staggered_in_x) Isdf = G%IsdB
  jsdf = G%jsd ; if (staggered_in_y) Jsdf = G%JsdB

  call horizontally_average_field(G, GV, isdf, jsdf, h, staggered_in_x, staggered_in_y, &
                                  is_layer, is_extensive, field, averaged_field, averaged_mask)

end procedure horizontally_average_diag_field
module procedure horizontally_average_field
  real :: volume(G%isc:G%iec, G%jsc:G%jec, size(field,3)) ! The area [L2 ~> m2], volume [L2 m ~> m3]
  real :: stuff(G%isc:G%iec, G%jsc:G%jec, size(field,3))  ! The area, volume or mass-weighted integral of the
  real, dimension(size(field, 3)) :: vol_sum   ! The global sum of the areas [m2], volumes [m3] or mass [kg]
  real, dimension(size(field, 3)) :: stuff_sum ! The global sum of the weighted field in all cells, in
  type(EFP_type), dimension(2*size(field,3)) :: sums_EFP ! Sums of volume or stuff by layer
  real :: height  ! An average thickness attributed to an velocity point [H ~> m or kg m-2]
  integer :: i, j, k, nz
  nz = size(field, 3)

  ! TODO: These averages could potentially be modified to use the function in
  !       the MOM_spatial_means module.
  ! NOTE: Reproducible sums must be computed in the original MKS units

  if (staggered_in_x .and. .not. staggered_in_y) then
    if (is_layer) then
      ! U-points
      do k=1,nz
        vol_sum(k) = 0.
        stuff_sum(k) = 0.
        if (is_extensive) then
          do j=G%jsc, G%jec ; do I=G%isc, G%iec
            volume(I,j,k) = G%areaCu(I,j) * G%mask2dCu(I,j)
            stuff(I,j,k) = volume(I,j,k) * field(I,j,k)
          enddo ; enddo
        else ! Intensive
          do j=G%jsc, G%jec ; do I=G%isc, G%iec
            height = 0.5 * (h(i,j,k) + h(i+1,j,k))
            volume(I,j,k) = G%areaCu(I,j)  * (GV%H_to_MKS * height) * G%mask2dCu(I,j)
            stuff(I,j,k) = volume(I,j,k) * field(I,j,k)
          enddo ; enddo
        endif
      enddo
    else ! Interface
      do k=1,nz
        do j=G%jsc, G%jec ; do I=G%isc, G%iec
          volume(I,j,k) = G%areaCu(I,j) * G%mask2dCu(I,j)
          stuff(I,j,k) = volume(I,j,k) * field(I,j,k)
        enddo ; enddo
      enddo
    endif
  elseif (staggered_in_y .and. .not. staggered_in_x) then
    if (is_layer) then
      ! V-points
      do k=1,nz
        if (is_extensive) then
          do J=G%jsc, G%jec ; do i=G%isc, G%iec
            volume(i,J,k) = G%areaCv(i,J) * G%mask2dCv(i,J)
            stuff(i,J,k) = volume(i,J,k) * field(i,J,k)
          enddo ; enddo
        else ! Intensive
          do J=G%jsc, G%jec ; do i=G%isc, G%iec
            height = 0.5 * (h(i,j,k) + h(i,j+1,k))
            volume(i,J,k) = G%areaCv(i,J) * (GV%H_to_MKS * height) * G%mask2dCv(i,J)
            stuff(i,J,k) = volume(i,J,k) * field(i,J,k)
          enddo ; enddo
        endif
      enddo
    else ! Interface
      do k=1,nz
        do J=G%jsc, G%jec ; do i=G%isc, G%iec
          volume(i,J,k) = G%areaCv(i,J) * G%mask2dCv(i,J)
          stuff(i,J,k) = volume(i,J,k) * field(i,J,k)
        enddo ; enddo
      enddo
    endif
  elseif ((.not. staggered_in_x) .and. (.not. staggered_in_y)) then
    if (is_layer) then
      ! H-points
      do k=1,nz
        if (is_extensive) then
          do j=G%jsc, G%jec ; do i=G%isc, G%iec
            if (h(i,j,k) > 0.) then
              volume(i,j,k) = G%areaT(i,j) * G%mask2dT(i,j)
              stuff(i,j,k) = volume(i,j,k) * field(i,j,k)
            else
              volume(i,j,k) = 0.
              stuff(i,j,k) = 0.
            endif
          enddo ; enddo
        else ! Intensive
          do j=G%jsc, G%jec ; do i=G%isc, G%iec
            volume(i,j,k) = G%areaT(i,j) * (GV%H_to_MKS * h(i,j,k)) * G%mask2dT(i,j)
            stuff(i,j,k) = volume(i,j,k) * field(i,j,k)
          enddo ; enddo
        endif
      enddo
    else ! Interface
      do k=1,nz
        do j=G%jsc, G%jec ; do i=G%isc, G%iec
          volume(i,j,k) = G%areaT(i,j) * G%mask2dT(i,j)
          stuff(i,j,k) = volume(i,j,k) * field(i,j,k)
        enddo ; enddo
      enddo
    endif
  else
    call assert(.false., 'horizontally_average_diag_field: Q point averaging is not coded yet.')
  endif

  ! Packing the sums into a single array with a single call to sum across PEs saves reduces
  ! the costs of communication.
  do k=1,nz
    sums_EFP(2*k-1) = reproducing_sum_EFP(volume(:,:,k), only_on_PE=.true., unscale=G%US%L_to_m**2)
    sums_EFP(2*k)   = reproducing_sum_EFP(stuff(:,:,k), only_on_PE=.true., unscale=G%US%L_to_m**2)
  enddo
  call EFP_sum_across_PEs(sums_EFP, 2*nz)
  do k=1,nz
    vol_sum(k) = EFP_to_real(sums_EFP(2*k-1))
    stuff_sum(k) = EFP_to_real(sums_EFP(2*k))
  enddo

  averaged_mask(:) = .true.
  do k=1,nz
    if (vol_sum(k) > 0.) then
      averaged_field(k) = stuff_sum(k) / vol_sum(k)
    else
      averaged_field(k) = 0.
      averaged_mask(k) = .false.
    endif
  enddo

end procedure horizontally_average_field
end submodule MOM_diag_remap_s
