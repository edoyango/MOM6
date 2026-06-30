submodule (shelfwave_initialization) shelfwave_initialization_s
#include <MOM_memory.h>
  implicit none
contains
module procedure register_shelfwave_OBC
  real :: PI      ! The ratio of the circumference of a circle to its diameter [nondim]
  character(len=32)  :: casename = "shelfwave"       !< This case's name.
  real :: jj      ! Cross-shore wave mode [nondim]
  real :: f0      ! Coriolis parameter [T-1 ~> s-1]
  real :: Lx      ! Long-shore length scale of bathymetry [km] or [m]
  real :: Ly      ! Cross-shore length scale [km] or [m]
  real :: default_amp  ! The default velocity amplitude [m s-1] or amplitude scaling factor [nondim]
  PI = 4.0*atan(1.0)

  if (associated(CS)) then
    call MOM_error(WARNING, "register_shelfwave_OBC called with an "// &
                            "associated control structure.")
    return
  endif
  allocate(CS)

  ! Register the tracer for horizontal advection & diffusion.
  call register_OBC(casename, param_file, OBC_Reg)
  call get_param(param_file, mdl, "F_0", f0, &
                 default=0.0, units="s-1", scale=US%T_to_s, do_not_log=.true.)
  call get_param(param_file, mdl,"SHELFWAVE_X_WAVELENGTH", Lx, &
                 "Length scale of shelfwave in x-direction.",&
                 units=G%x_ax_unit_short, default=100.)
  call get_param(param_file, mdl, "SHELFWAVE_Y_LENGTH_SCALE", Ly, &
                 "Length scale of exponential dropoff of topography in the y-direction.", &
                 units=G%y_ax_unit_short, default=50.)
  call get_param(param_file, mdl, "SHELFWAVE_Y_MODE", jj, &
                 "Cross-shore wave mode.",               &
                 units="nondim", default=1.)
  call get_param(param_file, mdl, "SHELFWAVE_CORRECT_AMPLITUDE", CS%shelfwave_correct_amplitude, &
                 "If true, SHELFWAVE_AMPLITUDE gives the actual inflow velocity, rather than giving "//&
                 "an overall scaling factor for the flow.", default=.true.)
  default_amp = 1.0 ; if (CS%shelfwave_correct_amplitude) default_amp = 0.1
  call get_param(param_file, mdl, "SHELFWAVE_AMPLITUDE", CS%my_amp, &
                 "Amplitude of the open boundary current inflows in the shelfwave configuration.", &
                 units="m s-1", default=default_amp, scale=US%m_s_to_L_T)

  CS%alpha = 1. / Ly
  CS%ll = 2. * PI / Lx
  CS%kk = jj * PI / G%len_lat
  CS%omega = 2 * CS%alpha * f0 * CS%ll / &
             (CS%kk*CS%kk + CS%alpha*CS%alpha + CS%ll*CS%ll)
  register_shelfwave_OBC = .true.

end procedure register_shelfwave_OBC
module procedure shelfwave_OBC_end
  if (associated(CS)) then
    deallocate(CS)
  endif
end procedure shelfwave_OBC_end
module procedure shelfwave_initialize_topography
  real      :: y    ! Position relative to the southern boundary [km] or [m] or [degrees_N]
  real      :: rLy  ! Exponential decay rate of the topography [km-1] or [m-1] or [degrees_N-1]
  real      :: Ly   ! Exponential decay lengthscale of the topography [km] or [m] or [degrees_N]
  real      :: H0   ! The minimum depth of the ocean [Z ~> m]
  integer   :: i, j
  call get_param(param_file, mdl,"SHELFWAVE_Y_LENGTH_SCALE", Ly, &
                 units=G%y_ax_unit_short, default=50., do_not_log=.true.)
  call get_param(param_file, mdl,"MINIMUM_DEPTH", H0, &
                 units="m", default=10., scale=US%m_to_Z, do_not_log=.true.)

  rLy = 0. ; if (Ly>0.) rLy = 1. / Ly

  do j=G%jsc,G%jec ; do i=G%isc,G%iec
    ! Compute normalized zonal coordinates (x,y=0 at center of domain)
    y = ( G%geoLatT(i,j) - G%south_lat )
    D(i,j) = H0 * exp(2 * rLy * y)
  enddo ; enddo

end procedure shelfwave_initialize_topography
module procedure shelfwave_set_OBC_data
  real :: time_sec ! The time in the run [T ~> s]
  real :: cos_wt, sin_wt ! Cosine and sine associated with the propagating x-direction structure [nondim]
  real :: cos_ky, sin_ky ! Cosine and sine associated with the y-direction structure [nondim]
  real :: x   ! Position relative to the western boundary [km] or [m] or [degrees_E]
  real :: y   ! Position relative to the southern boundary [km] or [m] or [degrees_N]
  real :: I_yscale  ! A factor to give the correct inflow velocity [km-1] or [m-1] or [degrees_N-1] or
  real :: my_amp    ! Amplitude of the open boundary current inflows, including sign changes
  integer :: i, j, is, ie, js, je, n
  integer :: turns    ! Number of index quarter turns
  type(OBC_segment_type), pointer :: segment => NULL()
  if (.not.associated(OBC)) return

  turns = modulo(G%HI%turns, 4)
  my_amp = CS%my_amp ; if ((turns==2) .or. (turns==3)) my_amp = -CS%my_amp

  time_sec = time_to_real(Time, scale=US%s_to_T)
  if (CS%shelfwave_correct_amplitude) then
    ! This makes the units and edge value of normal_vel_bt the same as my_amp.
    I_yscale = 1.0 / CS%kk
  else ! This preserves the previous answers.
    if (G%grid_unit_to_L == 0.0) call MOM_error(FATAL, &
          "shelfwave_set_OBC_data requires the use of Cartesian coordinates.")
    I_yscale = (1.0e3 * US%m_to_L) / G%grid_unit_to_L
  endif
  do n = 1, OBC%number_of_segments
    segment => OBC%segment(n)
    if (.not. segment%on_pe) cycle
    if (rotate_OBC_segment_direction(segment%direction, -turns) /= OBC_DIRECTION_W) cycle

    if (segment%is_E_or_W) then
      ! segment thicknesses are defined at cell face centers.
      is = segment%HI%isdB ; ie = segment%HI%iedB
      js = segment%HI%jsd ; je = segment%HI%jed
    else
      is = segment%HI%isd ; ie = segment%HI%ied
      js = segment%HI%jsdB ; je = segment%HI%jedB
    endif

    do j=js,je ; do I=is,ie
      if (segment%is_E_or_W) then
        x = G%geoLonCu(I,j) - G%west_lon
        y = G%geoLatCu(I,j) - G%south_lat
      else
        x = G%geoLonCv(i,J) - G%west_lon
        y = G%geoLatCv(i,J) - G%south_lat
      endif
      sin_wt = sin(CS%ll*x - CS%omega*time_sec)
      cos_wt = cos(CS%ll*x - CS%omega*time_sec)
      sin_ky = sin(CS%kk * y)
      cos_ky = cos(CS%kk * y)
      segment%normal_vel_bt(I,j) = my_amp * exp(- CS%alpha * y) * cos_wt * &
           ((CS%alpha * sin_ky + CS%kk * cos_ky) * I_yscale)
!     segment%tangential_vel_bt(I,j) = my_amp * (CS%ll * I_yscale) * exp(- CS%alpha * y) * sin_wt * sin_ky
!     segment%vorticity_bt(I,j) = my_amp * exp(- CS%alpha * y) * cos_wt * sin_ky * &
!           ((CS%ll**2 + CS%kk**2 + CS%alpha**2) * (I_yscale / G%grid_unit_to_L))
    enddo ; enddo
  enddo

end procedure shelfwave_set_OBC_data
end submodule shelfwave_initialization_s
