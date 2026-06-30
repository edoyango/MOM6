submodule (sloshing_initialization) sloshing_initialization_s
#include <MOM_memory.h>
  implicit none
contains
module procedure sloshing_initialize_topography
  integer   :: i, j
  do i=G%isc,G%iec ; do j=G%jsc,G%jec
    D(i,j) = max_depth
  enddo ; enddo

end procedure sloshing_initialize_topography
module procedure sloshing_initialize_thickness
  real    :: displ(SZK_(GV)+1)  ! The interface displacement [Z ~> m].
  real    :: z_unif(SZK_(GV)+1) ! Fractional uniform interface heights [nondim].
  real    :: z_inter(SZK_(GV)+1) ! Interface heights [Z ~> m]
  real    :: a0                 ! The displacement amplitude [Z ~> m].
  real    :: weight_z           ! A depth-space weighting [nondim].
  real    :: x1, y1, x2, y2     ! Dimensonless parameters specifying the depth profile [nondim]
  real    :: x, t               ! Dimensionless depth coordinates scales [nondim]
  logical :: use_IC_bug         ! If true, set the initial conditions retaining an old bug.
# include "version_variable.h"
  character(len=40)  :: mdl = "sloshing_initialization" !< This module's name.
  integer :: i, j, k, is, ie, js, je, nz
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke

  if (.not.just_read) call log_version(param_file, mdl, version, "")
  call get_param(param_file, mdl, "SLOSHING_IC_AMPLITUDE", a0, &
                 "Initial amplitude of sloshing internal interface height "//&
                 "displacements it the sloshing test case.", &
                 units='m', default=75.0, scale=US%m_to_Z, do_not_log=just_read)
  call get_param(param_file, mdl, "SLOSHING_IC_BUG", use_IC_bug, &
                 "If true, use code with a bug to set the sloshing initial conditions.", &
                 default=.false., do_not_log=just_read)

  if (just_read) return ! All run-time parameters have been read, so return.

  ! Define thicknesses
  do j=G%jsc,G%jec ; do i=G%isc,G%iec

    ! Define uniform interfaces
    do k = 0,nz
      z_unif(k+1) = -real(k)/real(nz)
    enddo

    ! 1. Define stratification
    do k = 1,nz+1

      ! Thin pycnocline in the middle
      !z_inter(k) = (2.0**(n-1)) * (z_unif(k) + 0.5)**n - 0.5

      ! Thin pycnocline in the middle (piecewise linear profile)
      x1 = 0.30 ; y1 = 0.48 ; x2 = 0.70 ; y2 = 0.52

      x = -z_unif(k)

      if ( x <= x1 ) then
        t = y1*x/x1
      elseif ( (x > x1 ) .and. ( x < x2 )) then
        t = y1 + (y2-y1) * (x-x1) / (x2-x1)
      else
        t = y2 + (1.0-y2) * (x-x2) / (1.0-x2)
      endif

      t = - z_unif(k)

      z_inter(k) = -t * G%max_depth

    enddo

    ! 2. Define displacement
    ! a0 is set via get_param; by default a0 is a 75m Displacement amplitude in depth units.
    do k = 1,nz+1

      weight_z = - 4.0 * ( z_unif(k) + 0.5 )**2 + 1.0

      x = G%geoLonT(i,j) / G%len_lon
      if (use_IC_bug) then
        displ(k) = a0 * cos(acos(-1.0)*x) + weight_z * US%m_to_Z ! There is a flag to fix this bug.
      else
        displ(k) = a0 * cos(acos(-1.0)*x) * weight_z
      endif

      if ( k == 1 ) then
        displ(k) = 0.0
      endif

      if ( k == nz+1 ) then
        displ(k) = 0.0
      endif

      z_inter(k) = z_inter(k) + displ(k)

    enddo

    ! 3. The last interface must coincide with the seabed
    z_inter(nz+1) = -depth_tot(i,j)
    ! Modify interface heights to make sure all thicknesses are strictly positive
    do k = nz,1,-1
      if ( z_inter(k) < (z_inter(k+1) + GV%Angstrom_Z) ) then
        z_inter(k) = z_inter(k+1) + GV%Angstrom_Z
      endif
    enddo

    ! 4. Define layers
    do k = 1,nz
      h(i,j,k) = z_inter(k) - z_inter(k+1)
    enddo

  enddo ; enddo

end procedure sloshing_initialize_thickness
module procedure sloshing_initialize_temperature_salinity
  real    :: delta_T            ! Temperature difference between layers [C ~> degC]
  real    :: S_ref, T_ref       ! Reference salinity [S ~> ppt] and temperature [C ~> degC] within surface layer
  real    :: S_range, T_range   ! Range of salinities [S ~> ppt] and temperatures [C ~> degC] over the vertical
  real    :: S_surf             ! Initial surface salinity [S ~> ppt]
  real    :: T_pert             ! A perturbed temperature [C ~> degC]
  integer :: kdelta             ! Half the number of layers with the temperature perturbation
  real    :: deltah             ! Thickness of each layer [Z ~> m]
  real    :: xi0, xi1           ! Fractional vertical positions [nondim]
  character(len=40)  :: mdl = "sloshing_initialization" ! This module's name.
  integer :: i, j, k, is, ie, js, je, nz
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke

  call get_param(param_file, mdl, "S_REF", S_ref, 'Reference value for salinity', &
                 default=35.0, units="ppt", scale=US%ppt_to_S, do_not_log=just_read)
  call get_param(param_file, mdl, "T_REF", T_ref, 'Reference value for temperature', &
                 units='degC', scale=US%degC_to_C, fail_if_missing=.not.just_read, do_not_log=just_read)

  ! The default is to assume an increase by 2 ppt for the salinity and a uniform temperature.
  call get_param(param_file, mdl, "S_RANGE", S_range, 'Initial salinity range.', &
                 units="ppt", default=2.0, scale=US%ppt_to_S, do_not_log=just_read)
  call get_param(param_file, mdl, "T_RANGE", T_range, 'Initial temperature range', &
                 units='degC', default=0.0, scale=US%degC_to_C, do_not_log=just_read)
  call get_param(param_file, mdl, "INITIAL_SSS", S_surf, "Initial surface salinity", &
                 units="ppt", default=34.0, scale=US%ppt_to_S, do_not_log=just_read)
  call get_param(param_file, mdl, "SLOSHING_T_PERT", T_pert, &
                 'A mid-column temperature perturbation in the sloshing test case', &
                 units='degC', default=1.0, scale=US%degC_to_C, do_not_log=just_read)

  if (just_read) return ! All run-time parameters have been read, so return.

  ! Prescribe salinity
  !delta_S = S_range / ( GV%ke - 1.0 )

  !S(:,:,1) = S_ref
  !do k = 2,GV%ke
  !  S(:,:,k) = S(:,:,k-1) + delta_S
  !enddo

  deltah = G%max_depth / nz
  do j=js,je ; do i=is,ie
    xi0 = 0.0
    do k = 1,nz
      xi1 = xi0 + deltah / G%max_depth ! =  xi0 + 1.0 / real(nz)
      S(i,j,k) = S_surf + 0.5 * S_range * (xi0 + xi1)
      xi0 = xi1
    enddo
  enddo ; enddo

  ! Prescribe temperature
  delta_T = T_range / ( GV%ke - 1.0 )

  T(:,:,1) = T_ref
  do k = 2,GV%ke
    T(:,:,k) = T(:,:,k-1) + delta_T
  enddo
  kdelta = 2
  ! Perhaps the following lines should instead assign T() = T_pert + T_ref
  T(:,:,GV%ke/2 - (kdelta-1):GV%ke/2 + kdelta) = T_pert

end procedure sloshing_initialize_temperature_salinity
end submodule sloshing_initialization_s
