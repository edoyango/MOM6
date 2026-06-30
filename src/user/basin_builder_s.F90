submodule (basin_builder) basin_builder_s
#include <MOM_memory.h>
  implicit none
contains
module procedure basin_builder_topography
  character(len=17) :: pname1, pname2 ! For construction of parameter names
  character(len=20) :: funcs ! Basin build function
  real, dimension(20) :: pars ! Parameters for each function [various]
  real :: lon ! Longitude [degrees_E]
  real :: lat ! Latitude [degrees_N]
  integer :: i, j, n, n_funcs
  call MOM_mesg("  basin_builder.F90, basin_builder_topography: setting topography", 5)
  call log_version(param_file, mdl, version, "")

  do j=G%jsc,G%jec ; do i=G%isc,G%iec
    D(i,j) = 1.0
  enddo ; enddo

  call get_param(param_file, mdl, "BBUILDER_N", n_funcs, &
                 "Number of pieces of topography to use.", fail_if_missing=.true.)

  do n=1,n_funcs
    write( pname1, "('BBUILDER_',i3.3,'_FUNC')" ) n
    write( pname2, "('BBUILDER_',i3.3,'_PARS')" ) n
    call get_param(param_file, mdl, pname1, funcs, &
                   "The basin builder function to apply with parameters "//&
                   trim(pname2)//". Choices are: NS_COAST, EW_COAST, "//&
                   "CIRC_CONIC_RIDGE, NS_CONIC_RIDGE, CIRC_SCURVE_RIDGE, "//&
                   "NS_SCURVE_RIDGE.", &
                   fail_if_missing=.true.)
    pars(:) = 0.
    if (trim(lowercase(funcs)) == 'ns_coast') then
      call get_param(param_file, mdl, pname2, pars(1:5), &
                     "NS_COAST parameters: longitude, starting latitude, "//&
                     "ending latitude, footprint radius, shelf depth.", &
                     units="degrees_E,degrees_N,degrees_N,degrees,m", &
                     fail_if_missing=.true.)
      pars(5) = pars(5) / max_depth
      do j=G%jsc,G%jec ; do i=G%isc,G%iec
        lon = G%geoLonT(i,j)
        lat = G%geoLatT(i,j)
        D(i,j) = min( D(i,j), NS_coast(lon, lat, pars(1), pars(2), pars(3), pars(4), pars(5)) )
      enddo ; enddo
    elseif (trim(lowercase(funcs)) == 'ns_conic_ridge') then
      call get_param(param_file, mdl, pname2, pars(1:5), &
                     "NS_CONIC_RIDGE parameters: longitude, starting latitude, "//&
                     "ending latitude, footprint radius, ridge height.", &
                     units="degrees_E,degrees_N,degrees_N,degrees,m", &
                     fail_if_missing=.true.)
      pars(5) = pars(5) / max_depth
      do j=G%jsc,G%jec ; do i=G%isc,G%iec
        lon = G%geoLonT(i,j)
        lat = G%geoLatT(i,j)
        D(i,j) = min( D(i,j), NS_conic_ridge(lon, lat, pars(1), pars(2), pars(3), pars(4), pars(5)) )
      enddo ; enddo
    elseif (trim(lowercase(funcs)) == 'ns_scurve_ridge') then
      call get_param(param_file, mdl, pname2, pars(1:5), &
                     "NS_SCURVE_RIDGE parameters: longitude, starting latitude, "//&
                     "ending latitude, footprint radius, ridge height.", &
                     units="degrees_E,degrees_N,degrees_N,degrees,m", &
                     fail_if_missing=.true.)
      pars(5) = pars(5) / max_depth
      do j=G%jsc,G%jec ; do i=G%isc,G%iec
        lon = G%geoLonT(i,j)
        lat = G%geoLatT(i,j)
        D(i,j) = min( D(i,j), NS_scurve_ridge(lon, lat, pars(1), pars(2), pars(3), pars(4), pars(5)) )
      enddo ; enddo
    elseif (trim(lowercase(funcs)) == 'angled_coast') then
      call get_param(param_file, mdl, pname2, pars(1:4), &
                     "ANGLED_COAST parameters: longitude intersection with Equator, "//&
                     "latitude intersection with Prime Meridian, footprint radius, shelf depth.", &
                     units="degrees_E,degrees_N,degrees,m", &
                     fail_if_missing=.true.)
      pars(4) = pars(4) / max_depth
      do j=G%jsc,G%jec ; do i=G%isc,G%iec
        lon = G%geoLonT(i,j)
        lat = G%geoLatT(i,j)
        D(i,j) = min( D(i,j), angled_coast(lon, lat, pars(1), pars(2), pars(3), pars(4)) )
      enddo ; enddo
    elseif (trim(lowercase(funcs)) == 'ew_coast') then
      call get_param(param_file, mdl, pname2, pars(1:5), &
                     "EW_COAST parameters: latitude, starting longitude, "//&
                     "ending longitude, footprint radius, shelf depth.", &
                     units="degrees_N,degrees_E,degrees_E,degrees,m", &
                     fail_if_missing=.true.)
      pars(5) = pars(5) / max_depth
      do j=G%jsc,G%jec ; do i=G%isc,G%iec
        lon = G%geoLonT(i,j)
        lat = G%geoLatT(i,j)
        D(i,j) = min( D(i,j), EW_coast(lon, lat, pars(1), pars(2), pars(3), pars(4), pars(5)) )
      enddo ; enddo
    elseif (trim(lowercase(funcs)) == 'circ_conic_ridge') then
      call get_param(param_file, mdl, pname2, pars(1:5), &
                     "CIRC_CONIC_RIDGE parameters: center longitude, center latitude, "//&
                     "ring radius, footprint radius, ridge height.", &
                     units="degrees_E,degrees_N,degrees,degrees,m", &
                     fail_if_missing=.true.)
      pars(5) = pars(5) / max_depth
      do j=G%jsc,G%jec ; do i=G%isc,G%iec
        lon = G%geoLonT(i,j)
        lat = G%geoLatT(i,j)
        D(i,j) = min( D(i,j), circ_conic_ridge(lon, lat, pars(1), pars(2), pars(3), pars(4), pars(5)) )
      enddo ; enddo
    elseif (trim(lowercase(funcs)) == 'circ_scurve_ridge') then
      call get_param(param_file, mdl, pname2, pars(1:5), &
                     "CIRC_SCURVe_RIDGE parameters: center longitude, center latitude, "//&
                     "ring radius, footprint radius, ridge height.", &
                     units="degrees_E,degrees_N,degrees,degrees,m", &
                     fail_if_missing=.true.)
      pars(5) = pars(5) / max_depth
      do j=G%jsc,G%jec ; do i=G%isc,G%iec
        lon = G%geoLonT(i,j)
        lat = G%geoLatT(i,j)
        D(i,j) = min( D(i,j), circ_scurve_ridge(lon, lat, pars(1), pars(2), pars(3), pars(4), pars(5)) )
      enddo ; enddo
    else
      call MOM_error(FATAL, "basin_builder.F90, basin_builer_topography:\n"//&
                     "Unrecognized function "//trim(funcs))
    endif

  enddo ! n

  do j=G%jsc,G%jec ; do i=G%isc,G%iec
    ! Dimensionalize by scaling 1 to max_depth
    D(i,j) = D(i,j) * max_depth
  enddo ; enddo

end procedure basin_builder_topography
module procedure cone
  cone = max( 0., 1. - abs(x - x0) / L )
  if (present(clip)) cone = min(clip, cone)
end procedure cone
module procedure scurve
  real :: s ! A rescaled position [nondim]
  s = max( 0., min( 1.,( x - x0 ) / L ) )
  scurve = ( 3. - 2.*s ) * ( s * s )
end procedure scurve
module procedure cstprof
  real :: s ! A rescaled position [nondim]
  s = max( 0., min( 1.,( x - x0 ) / L ) )
  cstprof = sh * scurve(s-lf,0.,bf) + (1.-sh) * scurve(s - (1.-sf),0.,sf)
end procedure cstprof
module procedure dist_line_fixed_x
  real :: dx, yr, dy ! Relative positions in arbitrary units [A]
  dx = x - x0
  yr = min( max(y0,y1), max( min(y0,y1), y ) ) ! bound y by y0,y1
  dy = y - yr ! =0 within y0<y<y1, =y0-y for y<y0, =y-y1 for y>y1
  dist_line_fixed_x = sqrt( (dx*dx) + (dy*dy) )
end procedure dist_line_fixed_x
module procedure dist_line_fixed_y
  dist_line_fixed_y = dist_line_fixed_x(y, x, y0, x0, x1)
end procedure dist_line_fixed_y
module procedure angled_coast
  real :: r ! A relative position [degrees]
  real :: I_dr ! The inverse of a distance [degrees-1]
  I_dr = 1/sqrt( lat_mer*lat_mer + lon_eq*lon_eq )
  r = I_dr * ( lat_mer*lon + lon_eq*lat - lon_eq*lat_mer)
  angled_coast = cstprof(r, 0., dr, 0.125, 0.125, 0.5, sh)
end procedure angled_coast
module procedure NS_coast
  real :: r  ! A relative position [degrees]
  r = dist_line_fixed_x( lon, lat, lonC, lat0, lat1 )
  NS_coast = cstprof(r, 0., dlon, 0.125, 0.125, 0.5, sh)
end procedure NS_coast
module procedure EW_coast
  real :: r  ! A relative position [degrees]
  r = dist_line_fixed_y( lon, lat, lon0, lon1, latC )
  EW_coast = cstprof(r, 0., dlat, 0.125, 0.125, 0.5, sh)
end procedure EW_coast
module procedure NS_conic_ridge
  real :: r  ! A relative position [degrees]
  r = dist_line_fixed_x( lon, lat, lonC, lat0, lat1 )
  NS_conic_ridge = 1. - rh * cone(r, 0., dlon)
end procedure NS_conic_ridge
module procedure NS_scurve_ridge
  real :: r  ! A relative position [degrees]
  r = dist_line_fixed_x( lon, lat, lonC, lat0, lat1 )
  NS_scurve_ridge = 1. - rh * (1. - scurve(r, 0., dlon) )
end procedure NS_scurve_ridge
module procedure circ_conic_ridge
  real :: r  ! A relative position [degrees]
  real :: frac_ht ! The fractional height of the topography [nondim]
  r = sqrt( ((lon - lon0)**2) + ((lat - lat0)**2) ) ! Pseudo-distance from a point
  r = abs( r - ring_radius) ! Pseudo-distance from a circle
  frac_ht = cone(r, 0., ring_thickness, ridge_height) ! 0 .. frac_ridge_height
  circ_conic_ridge = 1. - frac_ht ! nondim depths (1-frac_ridge_height) .. 1
end procedure circ_conic_ridge
module procedure circ_scurve_ridge
  real :: r  ! A relative position [degrees]
  real :: s  ! A function of the normalized position [nondim]
  real :: frac_ht ! The fractional height of the topography [nondim]
  r = sqrt( ((lon - lon0)**2) + ((lat - lat0)**2) ) ! Pseudo-distance from a point
  r = abs( r - ring_radius) ! Pseudo-distance from a circle
  s = 1. - scurve(r, 0., ring_thickness) ! 0 .. 1
  frac_ht = s * ridge_height ! 0 .. frac_ridge_height
  circ_scurve_ridge = 1. - frac_ht ! nondim depths (1-frac_ridge_height) .. 1
end procedure circ_scurve_ridge
end submodule basin_builder_s
