submodule (baroclinic_zone_initialization) baroclinic_zone_initialization_s
#include <MOM_memory.h>
  implicit none
contains
module procedure bcz_params
  if (.not.just_read) &
    call log_version(param_file, mdl, version, 'Initialization of an analytic baroclinic zone')
  call openParameterBlock(param_file,'BCZIC')
  call get_param(param_file, mdl, "S_REF", S_ref, 'Reference salinity', &
                 units='ppt', default=35., scale=US%ppt_to_S, do_not_log=just_read)
  call get_param(param_file, mdl, "DSDZ", dSdz, 'Salinity stratification', &
                 units='ppt m-1', default=0.0, scale=US%ppt_to_S*US%Z_to_m, do_not_log=just_read)
  call get_param(param_file, mdl, "DELTA_S",delta_S, 'Salinity difference across baroclinic zone', &
                 units='ppt', default=0.0, scale=US%ppt_to_S, do_not_log=just_read)
  call get_param(param_file, mdl, "DSDX", dSdx,'Meridional salinity difference', &
                 units='ppt '//trim(G%x_ax_unit_short)//'-1', default=0.0, scale=US%ppt_to_S, do_not_log=just_read)
  call get_param(param_file, mdl, "T_REF", T_ref, 'Reference temperature', &
                 units='degC', default=10., scale=US%degC_to_C, do_not_log=just_read)
  call get_param(param_file, mdl, "DTDZ", dTdz, 'Temperature stratification', &
                 units='degC m-1', default=0.0, scale=US%degC_to_C*US%Z_to_m, do_not_log=just_read)
  call get_param(param_file, mdl, "DELTA_T", delta_T,'Temperature difference across baroclinic zone', &
                 units='degC', default=0.0, scale=US%degC_to_C, do_not_log=just_read)
  call get_param(param_file, mdl, "DTDX", dTdx,'Meridional temperature difference', &
                 units='degC '//trim(G%x_ax_unit_short)//'-1', default=0.0, scale=US%degC_to_C, do_not_log=just_read)
  call get_param(param_file, mdl, "L_ZONE", L_zone, 'Width of baroclinic zone', &
                 units=G%y_ax_unit_short, default=0.5*G%len_lat, do_not_log=just_read)
  call closeParameterBlock(param_file)

end procedure bcz_params
module procedure baroclinic_zone_init_temperature_salinity
  integer   :: i, j, k, is, ie, js, je, nz
  real      :: T_ref, delta_T ! Parameters describing temperature distribution [C ~> degC]
  real      :: dTdz           ! Vertical temperature gradients [C Z-1 ~> degC m-1]
  real      :: dTdx           ! Zonal temperature gradients [C axis_units-1 ~> degC axis_units-1]
  real      :: S_ref, delta_S ! Parameters describing salinity distribution [S ~> ppt]
  real      :: dSdz           ! Vertical salinity gradients [S Z-1 ~> ppt m-1]
  real      :: dSdx           ! Zonal salinity gradients [S axis_units-1 ~> ppt axis_units-1]
  real      :: L_zone         ! Width of baroclinic zone, often in [km] or [degrees_N], depending
  real      :: zc, zi         ! Depths in depth units [Z ~> m]
  real      :: x              ! X-position relative to the domain center [degrees_E] or [km] or [m]
  real      :: y              ! Y-position relative to the domain center [degrees_N] or [km] or [m]
  real      :: fn             ! A smooth function based on the position in the baroclinic zone [nondim]
  real      :: xs, xd, yd     ! Fractional x- and y-positions relative to the domain extent [nondim]
  real      :: PI             ! 3.1415926... calculated as 4*atan(1) [nondim]
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke

  call bcz_params(G, GV, US, param_file, S_ref, dSdz, delta_S, dSdx, T_ref, dTdz, &
                  delta_T, dTdx, L_zone, just_read)

  if (just_read) return ! All run-time parameters have been read, so return.

  T(:,:,:) = 0.
  S(:,:,:) = 0.
  PI = 4.*atan(1.)

  do j = G%jsc,G%jec ; do i = G%isc,G%iec
    zi = -depth_tot(i,j)
    x = G%geoLonT(i,j) - (G%west_lon + 0.5*G%len_lon) ! Relative to center of domain
    xd = x / G%len_lon ! -1/2 < xd 1/2
    y = G%geoLatT(i,j) - (G%south_lat + 0.5*G%len_lat) ! Relative to center of domain
    yd = y / G%len_lat ! -1/2 < yd 1/2
    if (L_zone/=0.) then
      xs = min(1., max(-1., x/L_zone)) ! -1 < ys < 1
      fn = sin((0.5*PI)*xs)
    else
      xs = sign(1., x) ! +/- 1
      fn = xs
    endif
    do k = nz, 1, -1
      zc = zi + 0.5*h(i,j,k)          ! Position of middle of cell
      zi = zi + h(i,j,k)              ! Top interface position
      T(i,j,k) = T_ref + dTdz * zc  & ! Linear temperature stratification
                 + dTdx * x         & ! Linear gradient
                 + delta_T * fn       ! Smooth fn of width L_zone
      S(i,j,k) = S_ref + dSdz * zc  & ! Linear temperature stratification
                 + dSdx * x         & ! Linear gradient
                 + delta_S * fn       ! Smooth fn of width L_zone
    enddo
  enddo ; enddo

end procedure baroclinic_zone_init_temperature_salinity
end submodule baroclinic_zone_initialization_s
