submodule (soliton_initialization) soliton_initialization_s
#include <MOM_memory.h>
  implicit none
contains
module procedure soliton_initialize_thickness
  real    :: max_depth  ! Maximum depth of the model bathymetry [Z ~> m]
  real    :: cg_max     ! The external wave speed based on max_depth [L T-1 ~> m s-1]
  real    :: beta       ! The meridional gradient of the Coriolis parameter [T-1 L-1 ~> s-1 m-1]
  real    :: L_eq       ! The equatorial deformation radius used in nondimensionalizing this problem [L ~> m]
  real    :: scale_pos  ! A conversion factor to nondimensionalize the axis units, usually [m-1]
  real    :: x0    ! Initial x-position of the soliton in the same units as geoLonT, often [m].
  real    :: y0    ! Initial y-position of the soliton in the same units as geoLatT, often [m].
  real    :: x, y  ! Nondimensionalized positions [nondim]
  real    :: I_nz  ! The inverse of the number of layers [nondim]
  real    :: val1  ! A nondimensionlized zonal decay scale [nondim]
  real    :: val2  ! An overall surface height anomaly amplitude [L T-1 ~> m s-1]
  real    :: val3  ! A decay factor [nondim]
  real    :: val4  ! The local velocity amplitude [L T-1 ~> m s-1]
# include "version_variable.h"
  integer :: i, j, k, is, ie, js, je, nz
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke

  if (.not.just_read) &
    call MOM_mesg("soliton_initialization.F90, soliton_initialize_thickness: setting thickness")

  if (.not.just_read) call log_version(param_file, mdl, version, "")
  call get_param(param_file, mdl, "MAXIMUM_DEPTH", max_depth, &
                 units="m", default=-1.e9, scale=US%m_to_Z, do_not_log=.true.)
  call get_param(param_file, mdl, "BETA", beta, &
                 "The northward gradient of the Coriolis parameter with the betaplane option.", &
                 units="m-1 s-1", default=0.0, scale=US%T_to_s*US%L_to_m, do_not_log=.true.)

  if (just_read) return ! All run-time parameters have been read, so return.

  if (max_depth <= 0.0) call MOM_error(FATAL, &
      "soliton_initialization, soliton_initialize_thickness: "//&
      "This module requires a positive value of MAXIMUM_DEPTH.")
  if (abs(beta) <= 0.0) call MOM_error(FATAL, &
      "soliton_initialization, soliton_initialize_thickness: "//&
      "This module requires a non-zero value of BETA.")
  if (G%grid_unit_to_L <= 0.) call MOM_error(FATAL, "soliton_initialization.F90: "//&
          "soliton_initialize_thickness() is only set to work with Cartesian axis units.")

  cg_max = sqrt(GV%g_Earth * max_depth)
  L_eq = sqrt(cg_max / abs(beta))
  scale_pos = G%grid_unit_to_L / L_eq
  I_nz = 1.0 / real(nz)

  x0 = 2.0*G%len_lon/3.0
  y0 = 0.0
  val1 = 0.395
  val2 = max_depth * 0.771*(val1*val1)

  do j = G%jsc,G%jec ; do i = G%isc,G%iec
    do k = 1, nz
      x = (G%geoLonT(i,j)-x0) * scale_pos
      y = (G%geoLatT(i,j)-y0) * scale_pos
      val3 = exp(-val1*x)
      val4 = val2 * ( 2.0*val3 / (1.0 + (val3*val3)) )**2
      h(i,j,k) = (0.25*val4*(6.0*y*y + 3.0) * exp(-0.5*y*y) + depth_tot(i,j)) * I_nz
    enddo
  enddo ; enddo

end procedure soliton_initialize_thickness
module procedure soliton_initialize_velocity
  real    :: max_depth  ! Maximum depth of the model bathymetry [Z ~> m]
  real    :: cg_max     ! The external wave speed based on max_depth [L T-1 ~> m s-1]
  real    :: beta       ! The meridional gradient of the Coriolis parameter [T-1 L-1 ~> s-1 m-1]
  real    :: L_eq       ! The equatorial deformation radius used in nondimensionalizing this problem [L ~> m]
  real    :: scale_pos  ! A conversion factor to nondimensionalize the axis units, usually [m-1]
  real    :: x0    ! Initial x-position of the soliton in the same units as geoLonT, often [m].
  real    :: y0    ! Initial y-position of the soliton in the same units as geoLatT, often [m].
  real    :: x, y  ! Nondimensionalized positions [nondim]
  real    :: val1  ! A nondimensionlized zonal decay scale [nondim]
  real    :: val2  ! An overall velocity amplitude [L T-1 ~> m s-1]
  real    :: val3  ! A decay factor [nondim]
  real    :: val4  ! The local velocity amplitude [L T-1 ~> m s-1]
  integer :: i, j, k, is, ie, js, je, nz
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke

  if (.not.just_read) &
    call MOM_mesg("soliton_initialization.F90, soliton_initialize_thickness: setting thickness")

  call get_param(param_file, mdl, "MAXIMUM_DEPTH", max_depth, &
                 units="m", default=-1.e9, scale=US%m_to_Z, do_not_log=.true.)
  call get_param(param_file, mdl, "BETA", beta, &
                 "The northward gradient of the Coriolis parameter with the betaplane option.", &
                 units="m-1 s-1", default=0.0, scale=US%T_to_s*US%L_to_m, do_not_log=.true.)

  if (just_read) return ! All run-time parameters have been read, so return.

  if (max_depth <= 0.0) call MOM_error(FATAL, &
      "soliton_initialization, soliton_initialize_velocity: "//&
      "This module requires a positive value of MAXIMUM_DEPTH.")
  if (abs(beta) <= 0.0) call MOM_error(FATAL, &
      "soliton_initialization, soliton_initialize_velocity: "//&
      "This module requires a non-zero value of BETA.")
  if (G%grid_unit_to_L <= 0.) call MOM_error(FATAL, "soliton_initialization.F90: "//&
          "soliton_initialize_velocity() is only set to work with Cartesian axis units.")

  cg_max = sqrt(GV%g_Earth * max_depth)
  L_eq = sqrt(cg_max / abs(beta))
  scale_pos = G%grid_unit_to_L / L_eq

  x0 = 2.0*G%len_lon/3.0
  y0 = 0.0
  val1 = 0.395
  val2 = cg_max * 0.771*(val1*val1)

  v(:,:,:) = 0.0
  u(:,:,:) = 0.0

  do j = G%jsc,G%jec ; do I = G%isc-1,G%iec+1
    do k = 1, nz
      x = (0.5*(G%geoLonT(i+1,j)+G%geoLonT(i,j))-x0) * scale_pos
      y = (0.5*(G%geoLatT(i+1,j)+G%geoLatT(i,j))-y0) * scale_pos
      val3 = exp(-val1*x)
      val4 = val2*((2.0*val3/(1.0+(val3*val3)))**2)
      u(I,j,k) = 0.25*val4*(6.0*y*y-9.0) * exp(-0.5*y*y)
    enddo
  enddo ; enddo
  do j = G%jsc-1,G%jec+1 ; do I = G%isc,G%iec
    do k = 1, nz
      x = 0.5*(G%geoLonT(i,j+1)+G%geoLonT(i,j))-x0
      y = 0.5*(G%geoLatT(i,j+1)+G%geoLatT(i,j))-y0
      val3 = exp(-val1*x)
      val4 = val2*((2.0*val3/(1.0+(val3*val3)))**2)
      v(i,J,k) = 2.0*val4*y*(-2.0*val1*tanh(val1*x)) * exp(-0.5*y*y)
    enddo
  enddo ; enddo

end procedure soliton_initialize_velocity
end submodule soliton_initialization_s
