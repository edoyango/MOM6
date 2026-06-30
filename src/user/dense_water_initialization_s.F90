submodule (dense_water_initialization) dense_water_initialization_s
#include <MOM_memory.h>
  implicit none
contains
module procedure dense_water_initialize_topography
  real, dimension(5) :: domain_params ! nondimensional widths of all domain sections [nondim]
  real :: sill_frac     ! Depth of the sill separating downslope from upslope, as a fraction of
  real :: shelf_frac    ! Depth of the shelf region accumulating dense water for overflow,
  real :: x             ! Horizontal position normalized by the domain width [nondim]
  integer :: i, j
  call get_param(param_file, mdl, "DENSE_WATER_DOMAIN_PARAMS", domain_params, &
       "Fractional widths of all the domain sections for the dense water experiment.\n"//&
       "As a 5-element vector:\n"//&
       "  - open ocean, the section at maximum depth\n"//&
       "  - downslope, the downward overflow slope\n"//&
       "  - sill separating downslope from upslope\n"//&
       "  - upslope, the upward slope accumulating dense water\n"//&
       "  - the shelf in the dense formation region.", &
       units="nondim", fail_if_missing=.true.)
  call get_param(param_file, mdl, "DENSE_WATER_SILL_DEPTH", sill_frac, &
       "Depth of the sill separating downslope from upslope, as fraction of basin depth.", &
       units="nondim", default=default_sill)
  call get_param(param_file, mdl, "DENSE_WATER_SHELF_DEPTH", shelf_frac, &
       "Depth of the shelf region accumulating dense water for overflow, as fraction of basin depth.", &
       units="nondim", default=default_shelf)

  do i = 2, 5
    ! turn widths into positions
    domain_params(i) = domain_params(i-1) + domain_params(i)
  enddo

  do j = G%jsc,G%jec
    do i = G%isc,G%iec
      ! compute normalised zonal coordinate
      x = (G%geoLonT(i,j) - G%west_lon) / G%len_lon

      if (x <= domain_params(1)) then
        ! open ocean region
        D(i,j) = max_depth
      elseif (x <= domain_params(2)) then
        ! downslope region, linear
        D(i,j) = max_depth - (1.0 - sill_frac) * max_depth * &
             (x - domain_params(1)) / (domain_params(2) - domain_params(1))
      elseif (x <= domain_params(3)) then
        ! sill region
        D(i,j) = sill_frac * max_depth
      elseif (x <= domain_params(4)) then
        ! upslope region
        D(i,j) = sill_frac * max_depth + (shelf_frac - sill_frac) * max_depth * &
             (x - domain_params(3)) / (domain_params(4) - domain_params(3))
      else
        ! shelf region
        D(i,j) = shelf_frac * max_depth
      endif
    enddo
  enddo

end procedure dense_water_initialize_topography
module procedure dense_water_initialize_TS
  real :: mld             ! The initial mixed layer depth as a fraction of the maximum depth [nondim]
  real :: S_ref, S_range  ! The reference salinity and its range in the initial conditions [S ~> ppt]
  real :: T_ref           ! The reference temperature [C ~> degC]
  real :: zi, zmid        ! Depths from the surface nondimensionalized by the maximum depth [nondim]
  integer :: i, j, k, nz
  nz = GV%ke

  call get_param(param_file, mdl, "DENSE_WATER_MLD", mld, &
       "Depth of unstratified mixed layer as a fraction of the water column.", &
       units="nondim", default=default_mld, do_not_log=just_read)
  call get_param(param_file, mdl, "S_REF", S_ref, 'Reference salinity', &
                 default=35.0, units="ppt", scale=US%ppt_to_S, do_not_log=just_read)
  call get_param(param_file, mdl,"T_REF", T_ref, 'Reference temperature', &
                units='degC', scale=US%degC_to_C, fail_if_missing=.not.just_read, do_not_log=just_read)
  call get_param(param_file, mdl,"S_RANGE", S_range, 'Initial salinity range', &
                units="ppt", default=2.0, scale=US%ppt_to_S, do_not_log=just_read)

  if (just_read) return ! All run-time parameters have been read, so return.

  ! uniform temperature everywhere
  T(:,:,:) = T_ref

  do j = G%jsc,G%jec
    do i = G%isc,G%iec
      zi = 0.
      do k = 1,nz
        ! nondimensional middle of layer
        zmid = zi + 0.5 * h(i,j,k) / G%max_depth

        if (zmid < mld) then
          ! use reference salinity in the mixed layer
          S(i,j,k) = S_ref
        else
          ! linear between bottom of mixed layer and bottom
          S(i,j,k) = S_ref + S_range * (zmid - mld) / (1.0 - mld)
        endif

        zi = zi + h(i,j,k) / G%max_depth
      enddo
    enddo
  enddo
end procedure dense_water_initialize_TS
module procedure dense_water_initialize_sponges
  real :: west_sponge_time_scale, east_sponge_time_scale ! Sponge timescales [T ~> s]
  real :: west_sponge_width ! The fraction of the domain in which the western (outflow) sponge is active [nondim]
  real :: east_sponge_width ! The fraction of the domain in which the eastern (outflow) sponge is active [nondim]
  real, dimension(SZI_(G),SZJ_(G)) :: Idamp ! inverse damping timescale [T-1 ~> s-1]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)) :: dz ! sponge layer thicknesses in height units [Z ~> m]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)) :: T  ! sponge temperature [C ~> degC]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)) :: S  ! sponge salinity [S ~> ppt]
  real, dimension(SZK_(GV)+1) :: e0, eta1D ! interface positions for ALE sponge [Z ~> m]
  real :: x         ! Horizontal position normalized by the domain width [nondim]
  real :: zi, zmid  ! Depths from the surface nondimensionalized by the maximum depth [nondim]
  real :: dist      ! Distance from the edge of a sponge normalized by the width of that sponge [nondim]
  real :: mld       ! The initial mixed layer depth as a fraction of the maximum depth [nondim]
  real :: S_ref, S_range  ! The reference salinity and its range in the initial conditions [S ~> ppt]
  real :: S_dense   ! The salinity of the dense water being formed on the shelf [S ~> ppt]
  real :: T_ref     ! The reference temperature [C ~> degC]
  real :: sill_frac ! Fractional depths of the sill, relative to the maximum depth [nondim]
  integer :: i, j, k, nz
  nz = GV%ke

  call get_param(param_file, mdl, "DENSE_WATER_WEST_SPONGE_TIME_SCALE", west_sponge_time_scale, &
                 "The time scale on the west (outflow) of the domain for restoring. "//&
                 "If zero, the sponge is disabled.", units="s", default=0., scale=US%s_to_T)
  call get_param(param_file, mdl, "DENSE_WATER_WEST_SPONGE_WIDTH", west_sponge_width, &
                 "The fraction of the domain in which the western (outflow) sponge is active.", &
                 units="nondim", default=0.1)
  call get_param(param_file, mdl, "DENSE_WATER_EAST_SPONGE_TIME_SCALE", east_sponge_time_scale, &
                 "The time scale on the east (outflow) of the domain for restoring. "//&
                 "If zero, the sponge is disabled.", units="s", default=0., scale=US%s_to_T)
  call get_param(param_file, mdl, "DENSE_WATER_EAST_SPONGE_WIDTH", east_sponge_width, &
                 "The fraction of the domain in which the eastern (outflow) sponge is active.", &
                 units="nondim", default=0.1)
  call get_param(param_file, mdl, "DENSE_WATER_EAST_SPONGE_SALT", S_dense, &
                 "Salt anomaly of the dense water being formed in the overflow region.", &
                 units="ppt", default=4.0, scale=US%ppt_to_S)

  call get_param(param_file, mdl, "DENSE_WATER_MLD", mld, &
                 units="nondim", default=default_mld, do_not_log=.true.)
  call get_param(param_file, mdl, "DENSE_WATER_SILL_DEPTH", sill_frac, &
                 units="nondim", default=default_sill, do_not_log=.true.)

  call get_param(param_file, mdl, "S_REF", S_ref, &
                 units="ppt", default=35.0, scale=US%ppt_to_S, do_not_log=.true.)
  call get_param(param_file, mdl, "S_RANGE", S_range, &
                 units="ppt", default=2.0, scale=US%ppt_to_S, do_not_log=.true.)
  call get_param(param_file, mdl, "T_REF", T_ref, &
                 units='degC', scale=US%degC_to_C, fail_if_missing=.true., do_not_log=.true.)

  ! no active sponges
  if (west_sponge_time_scale <= 0. .and. east_sponge_time_scale <= 0.) return

  ! everywhere is initially unsponged
  Idamp(:,:) = 0.0

  do j = G%jsc, G%jec
    do i = G%isc,G%iec
      if (G%mask2dT(i,j) > 0.) then
        ! nondimensional x position
        x = (G%geoLonT(i,j) - G%west_lon) / G%len_lon

        if (west_sponge_time_scale > 0. .and. x < west_sponge_width) then
          dist = 1. - x / west_sponge_width
          ! scale restoring by depth into sponge
          Idamp(i,j) = 1. / west_sponge_time_scale * max(0., min(1., dist))
        elseif (east_sponge_time_scale > 0. .and. x > (1. - east_sponge_width)) then
          dist = 1. - (1. - x) / east_sponge_width
          Idamp(i,j) = 1. / east_sponge_time_scale * max(0., min(1., dist))
        endif
      endif
    enddo
  enddo

  if (use_ALE) then
    ! construct a uniform grid for the sponge
    do k = 1,nz
      e0(k) = -G%max_depth * (real(k - 1) / real(nz))
    enddo
    e0(nz+1) = -G%max_depth

    do j = G%jsc,G%jec
      do i = G%isc,G%iec
        eta1D(nz+1) = -depth_tot(i,j)
        do k = nz,1,-1
          eta1D(k) = e0(k)

          if (eta1D(k) < (eta1D(k+1) + GV%Angstrom_Z)) then
            ! is this layer vanished?
            eta1D(k) = eta1D(k+1) + GV%Angstrom_Z
            dz(i,j,k) = GV%Angstrom_Z
          else
            dz(i,j,k) = eta1D(k) - eta1D(k+1)
          endif
        enddo
      enddo
    enddo

    ! construct temperature and salinity for the sponge
    ! start with initial condition
    T(:,:,:) = T_ref
    S(:,:,:) = S_ref

    do j = G%jsc,G%jec
      do i = G%isc,G%iec
        zi = 0.
        x = (G%geoLonT(i,j) - G%west_lon) / G%len_lon
        do k = 1,nz
          ! nondimensional middle of layer
          zmid = zi + 0.5 * dz(i,j,k) / G%max_depth

          if (x > (1. - east_sponge_width)) then
            !if (zmid >= 0.9 * sill_frac) &
              S(i,j,k) = S_ref + S_dense
          else
            ! linear between bottom of mixed layer and bottom
            if (zmid >= mld) &
              S(i,j,k) = S_ref + S_range * (zmid - mld) / (1.0 - mld)
          endif

          zi = zi + dz(i,j,k) / G%max_depth
        enddo
      enddo
    enddo

    ! This call sets up the damping rates and interface heights in the sponges.
    call initialize_ALE_sponge(Idamp, G, GV, param_file, ACSp, dz, nz, data_h_is_Z=.true.)

    if ( associated(tv%T) ) call set_up_ALE_sponge_field(T, G, GV, tv%T, ACSp, 'temp', &
        sp_long_name='temperature', sp_unit='degC s-1')
    if ( associated(tv%S) ) call set_up_ALE_sponge_field(S, G, GV, tv%S, ACSp, 'salt', &
        sp_long_name='salinity', sp_unit='g kg-1 s-1')
  else
    call MOM_error(FATAL, "dense_water_initialize_sponges: trying to use non ALE sponge")
  endif
end procedure dense_water_initialize_sponges
end submodule dense_water_initialization_s
