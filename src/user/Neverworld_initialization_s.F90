submodule (Neverworld_initialization) Neverworld_initialization_s
#include <MOM_memory.h>
  implicit none
contains
module procedure Neverworld_initialize_topography
  real :: PI                   ! 3.1415926... calculated as 4*atan(1) [nondim]
  real :: x, y ! Lateral positions normalized by the domain size [nondim]
# include "version_variable.h"
  character(len=40)  :: mdl = "Neverworld_initialize_topography" ! This subroutine's name.
  real :: nl_top_amp       ! Amplitude of large-scale topographic features as a fraction of the maximum depth [nondim]
  real :: nl_roughness_amp ! Amplitude of topographic roughness as a fraction of the maximum depth [nondim]
  integer :: i, j, is, ie, js, je, isd, ied, jsd, jed
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed

  call MOM_mesg("  Neverworld_initialization.F90, Neverworld_initialize_topography: setting topography", 5)

  call log_version(param_file, mdl, version, "")
  call get_param(param_file, mdl, "NL_ROUGHNESS_AMP", nl_roughness_amp, &
                 "Amplitude of wavy signal in bathymetry.", units="nondim", default=0.05)
  call get_param(param_file, mdl, "NL_CONTINENT_AMP", nl_top_amp, &
                 "Scale factor for topography - 0.0 for no continents.", units="nondim", default=1.0)

  PI = 4.0*atan(1.0)

!  Calculate the depth of the bottom.
  do j=js,je ; do i=is,ie
    x = (G%geoLonT(i,j)-G%west_lon) / G%len_lon
    y = (G%geoLatT(i,j)-G%south_lat) / G%len_lat
!  This sets topography that has a reentrant channel to the south.
    D(i,j) = 1.0 - 1.1 * spike(y-1,0.12) - 1.1 * spike(y,0.12) - & !< The great northern wall and Antarctica
              nl_top_amp*( &
                (1.2 * spike(x,0.2) + 1.2 * spike(x-1.0,0.2)) * spike(MIN(0.0,y-0.3),0.2) & !< South America
              +  1.2 * spike(x-0.5,0.2) * spike(MIN(0.0,y-0.55),0.2)       & !< Africa
              +  1.2 * (spike(x,0.12)  + spike(x-1,0.12)) * spike(MAX(0.0,y-0.06),0.12)    & !< Antarctic Peninsula
              +  0.1 * (cosbell(x,0.1) + cosbell(x-1,0.1))                 & !< Drake Passage ridge
              +  0.5 * cosbell(x-0.16,0.05) * (cosbell(y-0.18,0.13)**0.4)  & !< Scotia Arc East
              +  0.4 * (cosbell(x-0.09,0.08)**0.4) * cosbell(y-0.26,0.05)  & !< Scotia Arc North
              +  0.4 * (cosbell(x-0.08,0.08)**0.4) * cosbell(y-0.1,0.05))   & !< Scotia Arc South
              -  nl_roughness_amp * cos(14*PI*x) * sin(14*PI*y)            & !< roughness
              -  nl_roughness_amp * cos(20*PI*x) * cos(20*PI*y)              !< roughness
    if (D(i,j) < 0.0) D(i,j) = 0.0
    D(i,j) = D(i,j) * max_depth
  enddo ; enddo

end procedure Neverworld_initialize_topography
module procedure cosbell
  real              :: PI      !< 3.1415926... calculated as 4*atan(1) [nondim]
  PI      = 4.0*atan(1.0)
  cosbell = 0.5 * (1 + cos(PI*MIN(ABS(x/L),1.0)))
end procedure cosbell
module procedure spike
  real              :: PI      !< 3.1415926... calculated as 4*atan(1) [nondim]
  PI    = 4.0*atan(1.0)
  spike = (1 - sin(PI*MIN(ABS(x/L),0.5)))
end procedure spike
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
module procedure NS_coast
  real :: r  ! A relative position [nondim]
  r = dist_line_fixed_x( lon, lat, lon0, lat0, lat1 )
  NS_coast = cstprof(r, 0., dlon, 0.125, 0.125, 0.5, sh)
end procedure NS_coast
module procedure EW_coast
  real :: r  ! A relative position [nondim]
  r = dist_line_fixed_y( lon, lat, lon0, lon1, lat0 )
  EW_coast = cstprof(r, 0., dlat, 0.125, 0.125, 0.5, sh)
end procedure EW_coast
module procedure NS_ridge
  real :: r ! A distance from a point [degrees]
  r = dist_line_fixed_x( lon, lat, lon0, lat0, lat1 )
  NS_ridge = 1. - rh * cone(r, 0., dlon)
end procedure NS_ridge
module procedure circ_ridge
  real :: r ! A relative position [degrees]
  real :: frac_ht ! The fractional height of the topography [nondim]
  r = sqrt( ((lon - lon0)**2) + ((lat - lat0)**2) ) ! Pseudo-distance from a point
  r = abs( r - ring_radius) ! Pseudo-distance from a circle
  frac_ht = cone(r, 0., ring_thickness, ridge_height) ! 0 .. frac_ridge_height
  circ_ridge = 1. - frac_ht ! Fractional depths (1-frac_ridge_height) .. 1
end procedure circ_ridge
module procedure Neverworld_initialize_thickness
  real :: e0(SZK_(GV)+1)    ! The resting interface heights, in depth units [Z ~> m],
  real, dimension(SZK_(GV)) :: h_profile ! Vector of initial thickness profile [Z ~> m].
  real :: e_interface ! Current interface position [Z ~> m].
  real :: x, y    ! horizontal coordinates for computation of the initial perturbation normalized
  real :: r1, r2  ! radial coordinates for computation of initial perturbation, normalized
  real :: pert_amp ! Amplitude of perturbations as a fraction of layer thicknesses [nondim]
  real :: h_noise ! Amplitude of noise to scale h by [nondim]
  real :: noise   ! Fractional noise in the layer thicknesses [nondim]
  type(randomNumberStream) :: rns ! Random numbers for stochastic tidal parameterization
  character(len=40)  :: mdl = "Neverworld_initialize_thickness" ! This subroutine's name.
  integer :: i, j, k, is, ie, js, je, nz
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke

  call MOM_mesg("  Neverworld_initialization.F90, Neverworld_initialize_thickness: setting thickness", 5)
  call get_param(param_file, mdl, "INIT_THICKNESS_PROFILE", h_profile, &
                 "Profile of initial layer thicknesses.", units="m", scale=US%m_to_Z, &
                 fail_if_missing=.true.)
  call get_param(param_file, mdl, "NL_THICKNESS_PERT_AMP", pert_amp, &
                 "Amplitude of finite scale perturbations as fraction of depth.", &
                 units="nondim", default=0.)
  call get_param(param_file, mdl, "NL_THICKNESS_NOISE_AMP", h_noise, &
                 "Amplitude of noise to scale layer by.", units="nondim", default=0.)

  ! e0 is the notional position of interfaces
  e0(1) = 0. ! The surface
  do k=1,nz
    e0(k+1) = e0(k) - h_profile(k)
  enddo

  do j=js,je ; do i=is,ie
    e_interface = -depth_tot(i,j)
    do k=nz,2,-1
      h(i,j,k) = e0(k) - e_interface ! Nominal thickness
      x = (G%geoLonT(i,j)-G%west_lon)/G%len_lon
      y = (G%geoLatT(i,j)-G%south_lat)/G%len_lat
      r1 = sqrt(((x-0.7)**2) + ((y-0.2)**2))
      r2 = sqrt(((x-0.3)**2) + ((y-0.25)**2))
      h(i,j,k) = h(i,j,k) + pert_amp * (e0(k) - e0(nz+1)) * &
                            (spike(r1,0.15)-spike(r2,0.15)) ! Prescribed perturbation
      if (h_noise /= 0.) then
        rns = initializeRandomNumberStream( int( 4096*(x + (y+1.)) ) )
        call getRandomNumbers(rns, noise) ! x will be in (0,1)
        noise = h_noise * 2. * ( noise - 0.5 ) ! range -h_noise to h_noise
        h(i,j,k) = ( 1. + noise ) * h(i,j,k)
      endif
      h(i,j,k) = max( GV%Angstrom_Z, h(i,j,k) ) ! Limit to non-negative
      e_interface = e_interface + h(i,j,k) ! Actual position of upper interface
    enddo
    h(i,j,1) = e0(1) - e_interface ! Nominal thickness
    h(i,j,1) = max( GV%Angstrom_Z, h(i,j,1) ) ! Limit to non-negative
  enddo ; enddo

end procedure Neverworld_initialize_thickness
end submodule Neverworld_initialization_s
