submodule (coord_hycom) coord_hycom_s
#include <MOM_memory.h>
  implicit none
contains
module procedure init_coord_hycom
  if (associated(CS)) call MOM_error(FATAL, "init_coord_hycom: CS already associated!")
  allocate(CS)
  allocate(CS%coordinateResolution(nk))
  allocate(CS%target_density(nk+1))

  CS%nk                      = nk
  CS%coordinateResolution(:) = coordinateResolution(:)
  CS%target_density(:)       = target_density(:)
  CS%use_3d                  = .false.
  CS%interp_CS               = interp_CS

  if (is_root_pe()) call MOM_error(NOTE, "init_coord_hycom: use_3d = .false.")

end procedure init_coord_hycom
module procedure init_3d_coord_hycom
  integer   :: i,j,k
  if (associated(CS)) call MOM_error(FATAL, "init_3d_coord_hycom: CS already associated!")

  allocate(CS)
  allocate(CS%coordinateResolution_3d(nk,SZI_(G),SZJ_(G)), source=0.0)
  allocate(CS%target_density_3d(nk+1,SZI_(G),SZJ_(G)), source=0.0)

  CS%nk        = nk
  CS%use_3d    = .true.
  CS%interp_CS = interp_CS

  do i=G%isc-1,G%iec+1 ; do j=G%jsc-1,G%jec+1
    if (G%mask2dT(i,j)>0.) then
      do k= 1,nk
        CS%coordinateResolution_3d(k,i,j) = coordinateResolution(i,j,k)
        CS%target_density_3d(k,i,j) = target_density(i,j,k)
      enddo
      CS%target_density_3d(nk+1,i,j) = target_density(i,j,nk+1)
    endif !mask2dT
  enddo ; enddo

  if (is_root_pe()) call MOM_error(NOTE, "init_3d_coord_hycom: use_3d = .true.")

end procedure init_3d_coord_hycom
module procedure end_coord_hycom
  if (.not. associated(CS)) return

  if (allocated(CS%coordinateResolution)) deallocate(CS%coordinateResolution)
  if (allocated(CS%target_density)) deallocate(CS%target_density)
  if (allocated(CS%coordinateResolution_3d)) deallocate(CS%coordinateResolution_3d)
  if (allocated(CS%target_density_3d)) deallocate(CS%target_density_3d)
  if (allocated(CS%max_interface_depths)) deallocate(CS%max_interface_depths)
  if (allocated(CS%max_layer_thickness)) deallocate(CS%max_layer_thickness)
  deallocate(CS)
end procedure end_coord_hycom
module procedure set_hycom_params
  if (.not. associated(CS)) call MOM_error(FATAL, "set_hycom_params: CS not associated")

  if (present(max_interface_depths)) then
    if (size(max_interface_depths) /= CS%nk+1) &
        call MOM_error(FATAL, "set_hycom_params: max_interface_depths inconsistent size")
    allocate(CS%max_interface_depths(CS%nk+1))
    CS%max_interface_depths(:) = max_interface_depths(:)
  endif

  if (present(max_layer_thickness)) then
    if (size(max_layer_thickness) /= CS%nk) &
        call MOM_error(FATAL, "set_hycom_params: max_layer_thickness inconsistent size")
    allocate(CS%max_layer_thickness(CS%nk))
    CS%max_layer_thickness(:) = max_layer_thickness(:)
  endif

  if (present(only_improves)) CS%only_improves = only_improves

  if (present(interp_CS)) CS%interp_CS = interp_CS
end procedure set_hycom_params
module procedure build_hycom1_column
  integer   :: k
  real, dimension(nz)      :: rho_col   ! Layer densities in a column [R ~> kg m-3]
  real, dimension(CS%nk)   :: h_col_new ! New layer thicknesses [H ~> m or kg m-2]
  real, dimension(CS%nk)   :: r_col_new ! New layer densities [R ~> kg m-3]
  real, dimension(CS%nk)   :: T_col_new ! New layer temperatures [C ~> degC]
  real, dimension(CS%nk)   :: S_col_new ! New layer salinities [S ~> ppt]
  real, dimension(CS%nk)   :: p_col_new ! New layer pressure [R L2 T-2 ~> Pa]
  real, dimension(CS%nk+1) :: RiA_ini   ! Initial nk+1 interface density anomaly w.r.t. the
  real, dimension(CS%nk+1) :: RiA_new   ! New interface density anomaly w.r.t. the
  real :: z_1, z_nz  ! mid point of 1st and last layers [H ~> m or kg m-2]
  real :: z_scale    ! A scaling factor from the input thicknesses to the target thicknesses,
  real :: stretching ! z* stretching, converts z* to z [nondim].
  real :: nominal_z ! Nominal depth of interface when using z* [H ~> m or kg m-2]
  logical :: maximum_depths_set ! If true, the maximum depths of interface have been set.
  logical :: maximum_h_set      ! If true, the maximum layer thicknesses have been set.
  maximum_depths_set = allocated(CS%max_interface_depths)
  maximum_h_set = allocated(CS%max_layer_thickness)

  z_scale = 1.0 ; if (present(zScale)) z_scale = zScale

  if (CS%only_improves .and. nz == CS%nk) then
    call build_hycom1_target_anomaly(CS, remapCS, eqn_of_state, CS%nk, ix, jy, depth, &
        h, T, S, p_col, rho_col, RiA_ini, h_neglect, h_neglect_edge)
  else
    ! Work bottom recording potential density
    call calculate_density(T, S, p_col, rho_col, eqn_of_state)
    ! This ensures the potential density profile is monotonic
    ! although not necessarily single valued.
    do k = nz-1, 1, -1
      rho_col(k) = min( rho_col(k), rho_col(k+1) )
    enddo
  endif

  ! Interpolates for the target interface position with the rho_col profile
  ! Based on global density profile, interpolate to generate a new grid
  if (CS%use_3d) then
    call build_and_interpolate_grid(CS%interp_CS, rho_col, nz, h(:), z_col, &
             CS%target_density_3d(:,ix,jy), CS%nk, h_col_new, z_col_new, h_neglect, h_neglect_edge)
  else
    call build_and_interpolate_grid(CS%interp_CS, rho_col, nz, h(:), z_col, &
             CS%target_density, CS%nk, h_col_new, z_col_new, h_neglect, h_neglect_edge)
  endif
  if (CS%only_improves .and. nz == CS%nk) then
    ! Only move an interface if it improves the density fit
    z_1 = 0.5 * ( z_col(1) + z_col(2) )
    z_nz  = 0.5 * ( z_col(nz) + z_col(nz+1) )
    do k = 1,CS%nk
      p_col_new(k) = p_col(1) + ( 0.5 * ( z_col_new(K) + z_col_new(K+1) ) - z_1 ) &
                                / ( z_nz - z_1 ) * ( p_col(nz) - p_col(1) )
    enddo
    ! Remap from original h and T,S to get T,S_col_new
    call remapping_core_h(remapCS, nz, h(:), T, CS%nk, h_col_new, T_col_new)
    call remapping_core_h(remapCS, nz, h(:), S, CS%nk, h_col_new, S_col_new)
    call build_hycom1_target_anomaly(CS, remapCS, eqn_of_state, CS%nk, ix, jy, depth, &
        h_col_new, T_col_new, S_col_new, p_col_new, r_col_new, RiA_new, h_neglect, h_neglect_edge)
    do k= 2,CS%nk
      if     ( abs(RiA_ini(K)) <= abs(RiA_new(K)) .and. z_col(K) > z_col_new(K-1) .and. &
               z_col(K) < z_col_new(K+1)) then
        z_col_new(K) = z_col(K)
      endif
    enddo
  endif !only_improves

  ! Sweep down the interfaces and make sure that the interface is at least
  ! as deep as a nominal target z* grid
  nominal_z = 0.
  stretching = z_col(nz+1) / depth ! Stretches z* to z
  if (CS%use_3d) then
    do k = 2, CS%nk+1
      nominal_z = nominal_z + (z_scale * CS%coordinateResolution_3d(k-1,ix,jy)) * stretching
      z_col_new(k) = max( z_col_new(k), nominal_z )
      z_col_new(k) = min( z_col_new(k), z_col(nz+1) )
    enddo
  else
    do k = 2, CS%nk+1
      nominal_z = nominal_z + (z_scale * CS%coordinateResolution(k-1)) * stretching
      z_col_new(k) = max( z_col_new(k), nominal_z )
      z_col_new(k) = min( z_col_new(k), z_col(nz+1) )
    enddo
  endif

  if (maximum_depths_set .and. maximum_h_set) then ; do k=2,CS%nk
    ! The loop bounds are 2 & nz so the top and bottom interfaces do not move.
    ! Recall that z_col_new is positive downward.
    z_col_new(K) = min(z_col_new(K), CS%max_interface_depths(K), &
                       z_col_new(K-1) + CS%max_layer_thickness(k-1))
  enddo ; elseif (maximum_depths_set) then ; do K=2,CS%nk
    z_col_new(K) = min(z_col_new(K), CS%max_interface_depths(K))
  enddo ; elseif (maximum_h_set) then ; do k=2,CS%nk
    z_col_new(K) = min(z_col_new(K), z_col_new(K-1) + CS%max_layer_thickness(k-1))
  enddo ; endif
end procedure build_hycom1_column
module procedure build_hycom1_target_anomaly
  integer   :: degree,k
  real, dimension(nz)   :: rho_col ! Layer densities in a column [R ~> kg m-3]
  real, dimension(nz,2) :: ppoly_E ! Polynomial edge values [R ~> kg m-3]
  real, dimension(nz,2) :: ppoly_S ! Polynomial edge slopes [R H-1]
  real, dimension(nz,DEGREE_MAX+1) :: ppoly_C ! Polynomial interpolant coeficients on the local 0-1 grid [R ~> kg m-3]
  call calculate_density(T, S, p_col, rho_col, eqn_of_state)
  ! This ensures the potential density profile is monotonic
  ! although not necessarily single valued.
  do k = nz-1, 1, -1
    rho_col(k) = min( rho_col(k), rho_col(k+1) )
  enddo

  call regridding_set_ppolys(CS%interp_CS, rho_col, nz, h, ppoly_E, ppoly_S, ppoly_C, &
                             degree, h_neglect, h_neglect_edge)

  if (CS%use_3d) then
    R(1) = rho_col(1)
    RiAnom(1) = ppoly_E(1,1) - CS%target_density_3d(1,ix,jy)
    do k= 2,nz
      R(k) = rho_col(k)
      if (ppoly_E(k-1,2) > CS%target_density_3d(k,ix,jy)) then
        RiAnom(k) = ppoly_E(k-1,2) - CS%target_density_3d(k,ix,jy)  !interface is heavier than target
      elseif (ppoly_E(k,1) < CS%target_density_3d(k,ix,jy)) then
        RiAnom(k) = ppoly_E(k,1)   - CS%target_density_3d(k,ix,jy)  !interface is lighter than target
      else
        RiAnom(k) = 0.0  !interface spans the target
      endif
    enddo
    RiAnom(nz+1) = ppoly_E(nz,2) - CS%target_density_3d(nz+1,ix,jy)
  else
    R(1) = rho_col(1)
    RiAnom(1) = ppoly_E(1,1) - CS%target_density(1)
    do k= 2,nz
      R(k) = rho_col(k)
      if (ppoly_E(k-1,2) > CS%target_density(k)) then
        RiAnom(k) = ppoly_E(k-1,2) - CS%target_density(k)  !interface is heavier than target
      elseif (ppoly_E(k,1) < CS%target_density(k)) then
        RiAnom(k) = ppoly_E(k,1)   - CS%target_density(k)  !interface is lighter than target
      else
        RiAnom(k) = 0.0  !interface spans the target
      endif
    enddo
    RiAnom(nz+1) = ppoly_E(nz,2) - CS%target_density(nz+1)
  endif !use_3d:else

end procedure build_hycom1_target_anomaly
end submodule coord_hycom_s
