submodule (circle_obcs_initialization) circle_obcs_initialization_s
#include <MOM_memory.h>
  implicit none
contains
module procedure circle_obcs_initialize_thickness
  real :: e0(SZK_(GV)+1)   ! The resting interface heights, in depth units [Z ~> m], usually
  real :: eta1D(SZK_(GV)+1)! Interface height relative to the sea surface
  real :: IC_amp           ! The amplitude of the initial height displacement [Z ~> m].
  real :: diskrad          ! Radius of the elevated disk [km] or [degrees] or [m]
  real :: rad              ! Distance from the center of the elevated disk [km] or [degrees] or [m]
  real :: lonC             ! The x-position of a point [km] or [degrees] or [m]
  real :: latC             ! The y-position of a point [km] or [degrees] or [m]
  real :: xOffset          ! The x-offset of the elevated disc center relative to the domain
# include "version_variable.h"
  character(len=40)  :: mdl = "circle_obcs_initialization"   ! This module's name.
  integer :: i, j, k, is, ie, js, je, nz
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke

  if (.not.just_read) &
    call MOM_mesg("  circle_obcs_initialization.F90, circle_obcs_initialize_thickness: setting thickness", 5)

  if (.not.just_read) call log_version(param_file, mdl, version, "")
  ! Parameters read by cartesian grid initialization
  call get_param(param_file, mdl, "DISK_RADIUS", diskrad, &
                 "The radius of the initially elevated disk in the "//&
                 "circle_obcs test case.", units=G%x_ax_unit_short, &
                 fail_if_missing=.not.just_read, do_not_log=just_read)
  call get_param(param_file, mdl, "DISK_X_OFFSET", xOffset, &
                 "The x-offset of the initially elevated disk in the "//&
                 "circle_obcs test case.", units=G%x_ax_unit_short, &
                 default=0.0, do_not_log=just_read)
  call get_param(param_file, mdl, "DISK_IC_AMPLITUDE", IC_amp, &
                 "Initial amplitude of interface height displacements "//&
                 "in the circle_obcs test case.", &
                 units='m', default=5.0, scale=US%m_to_Z, do_not_log=just_read)

  if (just_read) return ! All run-time parameters have been read, so return.

  do k=1,nz
    e0(K) = -G%max_depth * real(k-1) / real(nz)
  enddo

  ! Uniform thicknesses for base state
  do j=js,je ; do i=is,ie
    eta1D(nz+1) = -depth_tot(i,j)
    do k=nz,1,-1
      eta1D(K) = e0(K)
      if (eta1D(K) < (eta1D(K+1) + GV%Angstrom_Z)) then
        eta1D(K) = eta1D(K+1) + GV%Angstrom_Z
        h(i,j,k) = GV%Angstrom_Z
      else
        h(i,j,k) = eta1D(K) - eta1D(K+1)
      endif
    enddo
  enddo ; enddo

  ! Perturb base state by circular anomaly in center
  k=nz
  latC = G%south_lat + 0.5*G%len_lat
  lonC = G%west_lon + 0.5*G%len_lon + xOffset
  do j=js,je ; do i=is,ie
    rad = sqrt(((G%geoLonT(i,j)-lonC)**2) + ((G%geoLatT(i,j)-latC)**2)) / diskrad
    ! if (rad <= 6.*diskrad) h(i,j,k) = h(i,j,k)+10.0*exp( -0.5*( rad**2 ) )
    rad = min( rad, 1. ) ! Flatten outside radius of diskrad
    rad = rad*(2.*asin(1.)) ! Map 0-1 to 0-pi
    if (nz==1) then
      ! The model is barotropic
      h(i,j,k) = h(i,j,k) + IC_amp * 0.5*(1.+cos(rad)) ! cosine bell
    else
      ! The model is baroclinic
      do k = 1, nz
        h(i,j,k) = h(i,j,k) - 0.5*(1.+cos(rad)) * IC_amp * real( 2*k-nz )
      enddo
    endif
  enddo ; enddo

end procedure circle_obcs_initialize_thickness
end submodule circle_obcs_initialization_s
