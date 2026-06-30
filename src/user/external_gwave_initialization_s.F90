submodule (external_gwave_initialization) external_gwave_initialization_s
#include <MOM_memory.h>
  implicit none
contains
module procedure external_gwave_initialize_thickness
  real :: eta1D(SZK_(GV)+1)  ! Interface height relative to the sea surface
  real :: ssh_anomaly_height ! Vertical height of ssh anomaly [Z ~> m]
  real :: ssh_anomaly_width  ! Lateral width of anomaly, often in [km] or [degrees_E]
  character(len=40)  :: mdl = "external_gwave_initialize_thickness" ! This subroutine's name.
# include "version_variable.h"
  integer :: i, j, k, is, ie, js, je, nz
  real :: PI       ! The ratio of the circumference of a circle to its diameter [nondim]
  real :: Xnondim  ! A normalized x position [nondim]
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke

  if (.not.just_read) &
    call MOM_mesg("  external_gwave_initialization.F90, external_gwave_initialize_thickness: setting thickness", 5)

  if (.not.just_read) call log_version(param_file, mdl, version, "")
  call get_param(param_file, mdl, "SSH_ANOMALY_HEIGHT", ssh_anomaly_height, &
                 "The vertical displacement of the SSH anomaly. ", &
                 units="m", scale=US%m_to_Z, fail_if_missing=.not.just_read, do_not_log=just_read)
  call get_param(param_file, mdl, "SSH_ANOMALY_WIDTH", ssh_anomaly_width, &
                 "The lateral width of the SSH anomaly. ", &
                 units=G%x_ax_unit_short, fail_if_missing=.not.just_read, do_not_log=just_read)

  if (just_read) return ! All run-time parameters have been read, so return.

  PI = 4.0*atan(1.0)
  do j=G%jsc,G%jec ; do i=G%isc,G%iec
    Xnondim = (G%geoLonT(i,j)-G%west_lon-0.5*G%len_lon) / ssh_anomaly_width
    Xnondim = min(1., abs(Xnondim))
    eta1D(1) = ssh_anomaly_height * 0.5 * ( 1. + cos(PI*Xnondim) ) ! Cosine bell
    do k=2,nz
      eta1D(K) = -G%max_depth & ! Stretch interior interfaces with SSH
              + (eta1D(1)+G%max_depth) * ( real(nz+1-k)/real(nz) ) ! Stratification
    enddo
    eta1D(nz+1) = -G%max_depth ! Force bottom interface to bottom
    do k=1,nz
      h(i,j,k) = eta1D(K) - eta1D(K+1)
    enddo
  enddo ; enddo

end procedure external_gwave_initialize_thickness
end submodule external_gwave_initialization_s
