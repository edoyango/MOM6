submodule (tidal_bay_initialization) tidal_bay_initialization_s
#include <MOM_memory.h>
  implicit none
contains
module procedure register_tidal_bay_OBC
  character(len=32)  :: casename = "tidal bay"       !< This case's name.
  character(len=40)  :: mdl = "tidal_bay_initialization" ! This module's name.
  call get_param(param_file, mdl, "TIDAL_BAY_FLOW", CS%tide_flow, &
                 "Maximum total tidal volume flux.", &
                 units="m3 s-1", default=3.0e6, scale=US%m_s_to_L_T*US%m_to_L*US%m_to_Z)
  call get_param(param_file, mdl, "TIDAL_BAY_PERIOD", CS%tide_period, &
                 "Period of the inflow in the tidal bay configuration.", &
                 units="s", default=12.0*3600.0, scale=US%s_to_T)
  call get_param(param_file, mdl, "TIDAL_BAY_SSH_ANOM", CS%tide_ssh_amp, &
                 "Magnitude of the sea surface height anomalies at the inflow with the "//&
                 "tidal bay configuration.", &
                 units="m", default=0.1, scale=US%m_to_Z)

  ! Register the open boundaries.
  call register_OBC(casename, param_file, OBC_Reg)
  register_tidal_bay_OBC = .true.

end procedure register_tidal_bay_OBC
module procedure tidal_bay_set_OBC_data
  real :: time_sec    ! Elapsed model time [T ~> s]
  real :: cff_eta     ! The sea surface height anomalies associated with the inflow [Z ~> m]
  real :: my_flux     ! The volume flux through the face [L2 Z T-1 ~> m3 s-1]
  real :: total_area  ! The total face area of the OBCs [L Z ~> m2]
  real :: normal_vel  ! The normal velocity through the inflow face [L T-1 ~> m s-1]
  real :: PI          ! The ratio of the circumference of a circle to its diameter [nondim]
  real, allocatable :: my_area(:,:) ! The total OBC inflow area [L Z ~> m2]
  integer :: turns    ! Number of index quarter turns
  integer :: i, j, k, is, ie, js, je, isd, ied, jsd, jed, nz, n
  integer :: IsdB, IedB, JsdB, JedB
  type(OBC_segment_type), pointer :: segment => NULL()
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  IsdB = G%IsdB ; IedB = G%IedB ; JsdB = G%JsdB ; JedB = G%JedB

  PI = 4.0*atan(1.0)

  turns = modulo(G%HI%turns, 4)

  if (.not.associated(OBC)) return

  time_sec = time_to_real(Time, scale=US%s_to_T)
  cff_eta = CS%tide_ssh_amp * sin(2.0*PI*time_sec / CS%tide_period)

  segment => OBC%segment(1)

  if (turns == 0) then
    allocate(my_area(1:1,js:je), source=0.0)
    do j=segment%HI%jsc,segment%HI%jec ; do I=segment%HI%IscB,segment%HI%IecB
      if (OBC%segnum_u(I,j) > 0) then ! (segment%direction == OBC_DIRECTION_E)
        do k=1,nz
          my_area(1,j) = my_area(1,j) + h(i,j,k)*(GV%H_to_m*US%m_to_Z)*G%dyCu(I,j)
        enddo
      endif
    enddo ; enddo
  elseif (turns == 1) then
    allocate(my_area(is:ie,1:1), source=0.0)
    do J=segment%HI%JscB,segment%HI%JecB ; do i=segment%HI%isc,segment%HI%iec
      if (OBC%segnum_v(i,J) > 0) then ! (segment%direction == OBC_DIRECTION_N)
        do k=1,nz
          my_area(i,1) = my_area(i,1) + h(i,j,k)*(GV%H_to_m*US%m_to_Z)*G%dxCv(i,J)
        enddo
      endif
    enddo ; enddo
  elseif (turns == 2) then
    allocate(my_area(1:1,js:je), source=0.0)
    do j=segment%HI%jsc,segment%HI%jec ; do I=segment%HI%IscB,segment%HI%IecB
      if (OBC%segnum_u(I,j) < 0) then ! (segment%direction == OBC_DIRECTION_W)
        do k=1,nz
          my_area(1,j) = my_area(1,j) + h(i+1,j,k)*(GV%H_to_m*US%m_to_Z)*G%dyCu(I,j)
        enddo
      endif
    enddo ; enddo
  elseif (turns == 3) then
    allocate(my_area(is:ie,1:1), source=0.0)
    do J=segment%HI%JscB,segment%HI%JecB ; do i=segment%HI%isc,segment%HI%iec
      if (OBC%segnum_v(i,J) < 0) then ! (segment%direction == OBC_DIRECTION_S)
        do k=1,nz
          my_area(i,1) = my_area(i,1) + h(i,j+1,k)*(GV%H_to_m*US%m_to_Z)*G%dxCv(i,J)
        enddo
      endif
    enddo ; enddo
  endif

  total_area = reproducing_sum(my_area, unscale=US%Z_to_m*US%L_to_m)
  my_flux = - CS%tide_flow * SIN(2.0*PI*time_sec / CS%tide_period)
  normal_vel = my_flux / total_area
  if ((turns==2) .or. (turns==3)) normal_vel = -1.0 * normal_vel

  do n = 1, OBC%number_of_segments
    segment => OBC%segment(n)
    if (.not. segment%on_pe) cycle

    segment%normal_vel_bt(:,:) = normal_vel
    segment%SSH(:,:) = cff_eta

  enddo ! end segment loop

end procedure tidal_bay_set_OBC_data
end submodule tidal_bay_initialization_s
