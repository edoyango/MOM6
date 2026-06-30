submodule (lock_exchange_initialization) lock_exchange_initialization_s
#include <MOM_memory.h>
  implicit none
contains
module procedure lock_exchange_initialize_thickness
  real :: eta1D(SZK_(GV)+1)! Interface height relative to the sea surface
  real :: front_displacement ! Vertical displacement across front [Z ~> m]
  real :: thermocline_thickness ! Thickness of stratified region [Z ~> m]
# include "version_variable.h"
  character(len=40)  :: mdl = "lock_exchange_initialize_thickness" ! This subroutine's name.
  integer :: i, j, k, is, ie, js, je, nz
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke

  if (.not.just_read) &
    call MOM_mesg("  lock_exchange_initialization.F90, lock_exchange_initialize_thickness: setting thickness", 5)

  if (.not.just_read) call log_version(param_file, mdl, version, "")
  call get_param(param_file, mdl, "FRONT_DISPLACEMENT", front_displacement, &
                 "The vertical displacement of interfaces across the front. "//&
                 "A value larger in magnitude that MAX_DEPTH is truncated,", &
                 units="m", fail_if_missing=.not.just_read, do_not_log=just_read, scale=US%m_to_Z)
  call get_param(param_file, mdl, "THERMOCLINE_THICKNESS", thermocline_thickness, &
                 "The thickness of the thermocline in the lock exchange "//&
                 "experiment.  A value of zero creates a two layer system "//&
                 "with vanished layers in between the two inflated layers.", &
                 default=0., units="m", do_not_log=just_read, scale=US%m_to_Z)

  if (just_read) return ! All run-time parameters have been read, so return.

  do j=G%jsc,G%jec ; do i=G%isc,G%iec
    do k=2,nz
      eta1D(K) = -0.5 * G%max_depth & ! Middle of column
              - thermocline_thickness * ( (real(k-1))/real(nz) -0.5 ) ! Stratification
      if (G%geoLonT(i,j)-G%west_lon < 0.5 * G%len_lon) then
        eta1D(K) = eta1D(K) + 0.5 * front_displacement
      elseif (G%geoLonT(i,j)-G%west_lon > 0.5 * G%len_lon) then
        eta1D(K) = eta1D(K) - 0.5 * front_displacement
      endif
    enddo
    eta1D(nz+1) = -G%max_depth ! Force bottom interface to bottom
    do k=nz,2,-1 ! Make sure interfaces increase upwards
      eta1D(K) = max( eta1D(K), eta1D(K+1) + GV%Angstrom_Z )
    enddo
    eta1D(1) = 0. ! Force bottom interface to bottom
    do k=2,nz ! Make sure interfaces decrease downwards
      eta1D(K) = min( eta1D(K), eta1D(K-1) - GV%Angstrom_Z )
    enddo
    do k=nz,1,-1
      h(i,j,k) = eta1D(K) - eta1D(K+1)
    enddo
  enddo ; enddo

end procedure lock_exchange_initialize_thickness
end submodule lock_exchange_initialization_s
