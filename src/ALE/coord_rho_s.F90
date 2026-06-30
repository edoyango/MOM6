submodule (coord_rho) coord_rho_s
  implicit none
contains
module procedure init_coord_rho
  if (associated(CS)) call MOM_error(FATAL, "init_coord_rho: CS already associated!")
  allocate(CS)
  allocate(CS%target_density(nk+1))

  CS%nk                = nk
  CS%ref_pressure      = ref_pressure
  CS%target_density(:) = target_density(:)
  CS%interp_CS         = interp_CS

end procedure init_coord_rho
module procedure end_coord_rho
  if (.not. associated(CS)) return
  deallocate(CS%target_density)
  deallocate(CS)
end procedure end_coord_rho
module procedure set_rho_params
  if (.not. associated(CS)) call MOM_error(FATAL, "set_rho_params: CS not associated")

  if (present(min_thickness)) CS%min_thickness = min_thickness
  if (present(integrate_downward_for_e)) CS%integrate_downward_for_e = integrate_downward_for_e
  if (present(interp_CS)) CS%interp_CS = interp_CS
  if (present(ref_pressure)) CS%ref_pressure = ref_pressure
end procedure set_rho_params
module procedure build_rho_column
  integer :: k, count_nonzero_layers
  integer, dimension(nz) :: mapping
  real, dimension(nz) :: pres     ! Pressures used to calculate density [R L2 T-2 ~> Pa]
  real, dimension(nz) :: h_nv     ! Thicknesses of non-vanishing layers [H ~> m or kg m-2]
  real, dimension(nz) :: densities ! Layer density [R ~> kg m-3]
  real, dimension(nz+1) :: xTmp   ! Temporary positions [H ~> m or kg m-2]
  real, dimension(CS%nk) :: h_new ! New thicknesses [H ~> m or kg m-2]
  real, dimension(CS%nk+1) :: x1  ! Interface heights [H ~> m or kg m-2]
  call copy_finite_thicknesses(nz, h, CS%min_thickness, count_nonzero_layers, h_nv, mapping)

  if (count_nonzero_layers > 1) then
    xTmp(1) = 0.0
    do k = 1,count_nonzero_layers
      xTmp(k+1) = xTmp(k) + h_nv(k)
    enddo

    ! Compute densities on source column
    pres(:) = CS%ref_pressure
    call calculate_density(T, S, pres, densities, eqn_of_state)
    do k = 1,count_nonzero_layers
      densities(k) = densities(mapping(k))
    enddo

    ! Based on source column density profile, interpolate to generate a new grid
    call build_and_interpolate_grid(CS%interp_CS, densities, count_nonzero_layers, &
                                    h_nv, xTmp, CS%target_density, CS%nk, h_new, &
                                    x1, h_neglect, h_neglect_edge)

    ! Inflate vanished layers
    call old_inflate_layers_1d(CS%min_thickness, CS%nk, h_new)

    ! Comment: The following adjustment of h_new, and re-calculation of h_new via x1 needs to be removed
    x1(1) = 0.0 ; do k = 1,CS%nk ; x1(k+1) = x1(k) + h_new(k) ; enddo
    do k = 1,CS%nk
      h_new(k) = x1(k+1) - x1(k)
    enddo

  else ! count_nonzero_layers <= 1
    if (nz == CS%nk) then
      h_new(:) = h(:) ! This keeps old behavior
    else
      h_new(:) = 0.
      h_new(1) = h(1)
    endif
  endif

  ! Return interface positions
  if (CS%integrate_downward_for_e) then
    ! Remapping is defined integrating from zero
    z_interface(1) = 0.
    do k = 1,CS%nk
      z_interface(k+1) = z_interface(k) - h_new(k)
    enddo
  else
    ! The rest of the model defines grids integrating up from the bottom
    z_interface(CS%nk+1) = -depth
    do k = CS%nk,1,-1
      z_interface(k) = z_interface(k+1) + h_new(k)
    enddo
  endif

end procedure build_rho_column
module procedure build_rho_column_iteratively
  real, dimension(nz+1) :: x0, x1, xTmp ! Temporary interface heights [Z ~> m]
  real, dimension(nz) :: pres       ! The pressure used in the equation of state [R L2 T-2 ~> Pa].
  real, dimension(nz) :: densities  ! Layer densities [R ~> kg m-3]
  real, dimension(nz) :: T_tmp, S_tmp ! A temporary profile of temperature [C ~> degC] and salinity [S ~> ppt].
  real, dimension(nz) :: h0, h1, hTmp ! Temporary thicknesses [Z ~> m]
  real :: deviation            ! When iterating to determine the final grid, this is the
  real :: deviation_tol        ! Deviation tolerance between succesive grids in
  real :: threshold            ! The minimum thickness for a layer to be considered to exist [Z ~> m]
  integer, dimension(nz) :: mapping ! The indices of the massive layers in the initial column.
  integer :: k, m, count_nonzero_layers
  integer, parameter :: NB_REGRIDDING_ITERATIONS = 1
  threshold = CS%min_thickness
  pres(:) = CS%ref_pressure
  T_tmp(:) = T(:)
  S_tmp(:) = S(:)
  h0(:) = h(:)

  ! Start iterations to build grid
  m = 1
  deviation_tol = 1.0e-15*depth ; if (present(dev_tol)) deviation_tol = dev_tol

  do m=1,NB_REGRIDDING_ITERATIONS

    ! Construct column with vanished layers removed
    call copy_finite_thicknesses(nz, h0, threshold, count_nonzero_layers, hTmp, mapping)
    if ( count_nonzero_layers <= 1 ) then
      h1(:) = h0(:)
      exit  ! stop iterations here
    endif

    xTmp(1) = 0.0
    do k = 1,count_nonzero_layers
      xTmp(k+1) = xTmp(k) + hTmp(k)
    enddo

    ! Compute densities within current water column
    call calculate_density(T_tmp, S_tmp, pres, densities, eqn_of_state)

    do k = 1,count_nonzero_layers
      densities(k) = densities(mapping(k))
    enddo

    ! One regridding iteration
    ! Based on global density profile, interpolate to generate a new grid
    call build_and_interpolate_grid(CS%interp_CS, densities, count_nonzero_layers, &
         hTmp, xTmp, CS%target_density, nz, h1, x1, h_neglect, h_neglect_edge)

    call old_inflate_layers_1d( CS%min_thickness, nz, h1 )
    x1(1) = 0.0 ; do k = 1,nz ; x1(k+1) = x1(k) + h1(k) ; enddo

    ! Remap T and S from previous grid to new grid
    do k = 1,nz
      h1(k) = x1(k+1) - x1(k)
    enddo

    call remapping_core_h(remapCS, nz, h0, S, nz, h1, S_tmp)

    call remapping_core_h(remapCS, nz, h0, T, nz, h1, T_tmp)

    ! Compute the deviation between two successive grids
    deviation = 0.0
    x0(1) = 0.0
    x1(1) = 0.0
    do k = 2,nz
      x0(k) = x0(k-1) + h0(k-1)
      x1(k) = x1(k-1) + h1(k-1)
      deviation = deviation + (x0(k)-x1(k))**2
    enddo
    deviation = sqrt( deviation / (nz-1) )

    if ( deviation <= deviation_tol ) exit

    ! Copy final grid onto start grid for next iteration
    h0(:) = h1(:)
  enddo ! end regridding iterations

  if (CS%integrate_downward_for_e) then
    zInterface(1) = 0.
    do k = 1,nz
      zInterface(k+1) = zInterface(k) - h1(k)
      ! Adjust interface position to accommodate inflating layers
      ! without disturbing the interface above
    enddo
  else
    ! The rest of the model defines grids integrating up from the bottom
    zInterface(nz+1) = -depth
    do k = nz,1,-1
      zInterface(k) = zInterface(k+1) + h1(k)
      ! Adjust interface position to accommodate inflating layers
      ! without disturbing the interface above
    enddo
  endif

end procedure build_rho_column_iteratively
module procedure copy_finite_thicknesses
  integer :: k, k_thickest
  real :: thickness_in_vanished ! Summed thicknesses in discarded layers [H ~> m or kg m-2] or [Z ~> m]
  real :: thickest_h_out        ! Thickness of the thickest layer [H ~> m or kg m-2] or [Z ~> m]
  nout = 0
  thickness_in_vanished = 0.0
  thickest_h_out = h_in(1)
  k_thickest = 1
  do k = 1, nk
    mapping(k) = nout ! Note k>=nout always
    h_out(k) = 0.  ! Make sure h_out is set everywhere
    if (h_in(k) > thresh) then
      ! For non-vanished layers
      nout = nout + 1
      mapping(nout) = k
      h_out(nout) = h_in(k)
      if (h_out(nout) > thickest_h_out) then
        thickest_h_out = h_out(nout)
        k_thickest = nout
      endif
    else
      ! Add up mass in vanished layers
      thickness_in_vanished = thickness_in_vanished + h_in(k)
    endif
  enddo

  ! No finite layers
  if (nout <= 1) return

  ! Adjust for any lost volume in vanished layers
  h_out(k_thickest) = h_out(k_thickest) + thickness_in_vanished

end procedure copy_finite_thicknesses
module procedure old_inflate_layers_1d
  integer   :: k
  integer   :: k_found
  integer   :: count_nonzero_layers
  real      :: delta         ! An increase to a layer to increase it to the minimum thickness in the
  real      :: correction    ! The accumulated correction that will be applied to the thickest layer
  real      :: maxThickness  ! The thickness of the thickest layer in the same units as h, often [H ~> m or kg m-2]
  count_nonzero_layers = 0
  do k = 1,nk
    if ( h(k) > min_thickness ) then
      count_nonzero_layers = count_nonzero_layers + 1
    endif
  enddo

  ! If all layer thicknesses are greater than the threshold, exit routine
  if ( count_nonzero_layers == nk ) return

  ! If all thicknesses are zero, inflate them all and exit
  if ( count_nonzero_layers == 0 ) then
    do k = 1,nk
      h(k) = min_thickness
    enddo
    return
  endif

  ! Inflate zero layers
  correction = 0.0
  do k = 1,nk
    if ( h(k) <= min_thickness ) then
      delta = min_thickness - h(k)
      correction = correction + delta
      h(k) = h(k) + delta
    endif
  enddo

  ! Modify thicknesses of nonzero layers to ensure volume conservation
  maxThickness = h(1)
  k_found = 1
  do k = 1,nk
    if ( h(k) > maxThickness ) then
      maxThickness = h(k)
      k_found = k
    endif
  enddo

  h(k_found) = h(k_found) - correction

end procedure old_inflate_layers_1d
end submodule coord_rho_s
