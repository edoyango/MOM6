submodule (dyed_channel_initialization) dyed_channel_initialization_s
#include <MOM_memory.h>
  implicit none
contains
module procedure register_dyed_channel_OBC
  logical :: enable_bugs  ! If true, the defaults for recently added bug-fix flags are set to
  character(len=32)  :: casename = "dyed channel"     ! This case's name.
  character(len=40)  :: mdl = "register_dyed_channel_OBC" ! This subroutine's name.
  if (associated(CS)) then
    call MOM_error(WARNING, "register_dyed_channel_OBC called with an "// &
                            "associated control structure.")
    return
  endif
  allocate(CS)

  call get_param(param_file, mdl, "CHANNEL_MEAN_FLOW", CS%zonal_flow, &
                 "Mean zonal flow imposed at upstream open boundary.", &
                 units="m/s", default=8.57, scale=US%m_s_to_L_T)
  call get_param(param_file, mdl, "CHANNEL_TIDAL_AMP", CS%tidal_amp, &
                 "Sloshing amplitude imposed at upstream open boundary.", &
                 units="m/s", default=0.0, scale=US%m_s_to_L_T)
  call get_param(param_file, mdl, "CHANNEL_FLOW_FREQUENCY", CS%frequency, &
                 "Frequency of oscillating zonal flow.", &
                 units="s-1", default=0.0, scale=US%T_to_s)
  call get_param(param_file, mdl, "ENABLE_BUGS_BY_DEFAULT", enable_bugs, &
                 default=.true., do_not_log=.true.)  ! This is logged from MOM.F90.
  call get_param(param_file, mdl, "CHANNEL_FLOW_OBC_TRANSPORT_BUG", CS%OBC_transport_bug, &
                 "If true and specified open boundary conditions are being used, use a 1 m "//&
                 "(if Boussienesq) or 1 kg m-2 layer thickness instead of the actual thickness.", &
                 default=enable_bugs)

  ! Register the open boundaries.
  call register_OBC(casename, param_file, OBC_Reg)
  register_dyed_channel_OBC = .true.

end procedure register_dyed_channel_OBC
module procedure dyed_channel_OBC_end
  if (associated(CS)) then
    deallocate(CS)
  endif
end procedure dyed_channel_OBC_end
module procedure dyed_channel_set_OBC_tracer_data
  character(len=40)  :: mdl = "dyed_channel_set_OBC_tracer_data" ! This subroutine's name.
  character(len=80)  :: name, longname
  integer :: m, n, ntr_id
  real :: dye  ! Inflow dye concentrations [arbitrary]
  type(tracer_type), pointer      :: tr_ptr => NULL()
  if (.not.associated(OBC)) call MOM_error(FATAL, 'dyed_channel_initialization.F90: '// &
        'dyed_channel_set_OBC_data() was called but OBC type was not initialized!')

  call get_param(param_file, mdl, "NUM_DYE_TRACERS", ntr, &
                 "The number of dye tracers in this run. Each tracer "//&
                 "should have a separate boundary segment.", default=0,   &
                 do_not_log=.true.)

  if (OBC%number_of_segments < ntr) then
    call MOM_error(WARNING, "Error in dyed_obc segment setup")
    return   !!! Need a better error message here
  endif

! ! Set the inflow values of the dyes, one per segment.
! ! We know the order: north, south, east, west
  do m=1,ntr
    write(name,'("dye_",I2.2)') m
    write(longname,'("Concentration of dyed_obc Tracer ",I2.2, " on segment ",I2.2)') m, m
    call tracer_name_lookup(tr_Reg, ntr_id, tr_ptr, name)

    do n=1,OBC%number_of_segments
      if (n == m) then
        dye = 1.0
      else
        dye = 0.0
      endif
      call register_segment_tracer(tr_ptr, ntr_id, param_file, GV, &
                                   OBC%segment(n), OBC_scalar=dye)
    enddo
  enddo

end procedure dyed_channel_set_OBC_tracer_data
module procedure dyed_channel_update_flow
  real :: flow      ! The OBC velocity [L T-1 ~> m s-1]
  real :: PI        ! 3.1415926535... [nondim]
  real :: time_sec  ! The elapsed time since the start of the calendar [T ~> s]
  real :: fixed_thickness ! A fixed layer thickness, hard-coded to 1 mks unit, that is used to
  logical :: cross_channel  ! True if the segment runs across the channel
  integer :: turns    ! Number of index quarter turns
  integer :: i, j, k, l_seg, isd, ied, jsd, jed
  integer :: IsdB, IedB, JsdB, JedB, is, ie, js, je
  type(OBC_segment_type), pointer :: segment => NULL()
  if (.not.associated(OBC)) call MOM_error(FATAL, 'dyed_channel_initialization.F90: '// &
        'dyed_channel_update_flow() was called but OBC type was not initialized!')

  time_sec = time_to_real(Time, scale=US%s_to_T)
  PI = 4.0*atan(1.0)

  turns = modulo(G%HI%turns, 4)

  do l_seg=1, OBC%number_of_segments
    segment => OBC%segment(l_seg)
    if (.not. segment%on_pe) cycle
    if (segment%gradient) cycle
    if (segment%oblique .and. (.not. segment%nudged) .and. (.not. segment%Flather)) cycle

    if (CS%frequency == 0.0) then
      flow = CS%zonal_flow
    else
      flow = CS%zonal_flow + CS%tidal_amp * cos(2 * PI * CS%frequency * time_sec)
    endif
    if ((turns==2) .or. (turns==3)) flow = -1.0 * flow

    isd = segment%HI%isd ; ied = segment%HI%ied
    jsd = segment%HI%jsd ; jed = segment%HI%jed
    IsdB = segment%HI%IsdB ; IedB = segment%HI%IedB
    JsdB = segment%HI%JsdB ; JedB = segment%HI%JedB
    if (segment%is_E_or_W) then
      is = IsdB ; ie = IedB ; js = jsd ; je = jed
    else
      is = isd ; ie = ied ; js = JsdB ; je = JedB
    endif
    cross_channel = ((segment%is_E_or_W .and. ((turns==0) .or. (turns==2))) .or. &
                     (segment%is_N_or_S .and. ((turns==1) .or. (turns==3))))

    if ((segment%specified .or. segment%nudged) .and. cross_channel) then
      do k=1,GV%ke ; do j=js,je ; do I=is,ie
        segment%normal_vel(I,j,k) = flow
      enddo ; enddo ; enddo
    endif

    if (segment%specified .and. cross_channel) then
      if (CS%OBC_transport_bug) then
        fixed_thickness = 1.0 / GV%H_to_mks  ! This replicates the prevoius answers without rescaling.
        if ((segment%direction == OBC_DIRECTION_W) .or. (segment%direction == OBC_DIRECTION_E)) then
          do k=1,GV%ke ; do j=jsd,jed ; do I=IsdB,IedB
            segment%normal_trans(I,j,k) = flow * G%dyCu(I,j) * fixed_thickness
          enddo ; enddo ; enddo
        elseif ((segment%direction == OBC_DIRECTION_S) .or. (segment%direction == OBC_DIRECTION_N)) then
          do k=1,GV%ke ; do J=JsdB,JedB ; do i=isd,ied
            segment%normal_trans(i,J,k) = flow * G%dxCv(i,J) * fixed_thickness
          enddo ; enddo ; enddo
        endif
      else
        if (segment%direction == OBC_DIRECTION_W) then
          do k=1,GV%ke ; do j=jsd,jed ; do I=IsdB,IedB
            segment%normal_trans(I,j,k) = flow * G%dyCu(I,j) * h(i+1,j,k)
          enddo ; enddo ; enddo
        elseif (segment%direction == OBC_DIRECTION_E) then
          do k=1,GV%ke ; do j=jsd,jed ; do I=IsdB,IedB
            segment%normal_trans(I,j,k) = flow * G%dyCu(I,j) * h(i,j,k)
          enddo ; enddo ; enddo
        elseif (segment%direction == OBC_DIRECTION_S) then
          do k=1,GV%ke ; do J=JsdB,JedB ; do i=isd,ied
            segment%normal_trans(i,J,k) = flow * G%dxCv(i,J) * h(i,j+1,k)
          enddo ; enddo ; enddo
        elseif (segment%direction == OBC_DIRECTION_N) then
          do k=1,GV%ke ; do J=JsdB,JedB ; do i=isd,ied
            segment%normal_trans(i,J,k) = flow * G%dxCv(i,J) * h(i,j,k)
          enddo ; enddo ; enddo
        endif
      endif
    endif

    if (cross_channel) then
      do j=js,je ; do I=is,ie
        segment%normal_vel_bt(I,j) = flow
      enddo ; enddo
    else
      do J=js,je ; do i=is,ie
        segment%normal_vel_bt(i,J) = 0.0
      enddo ; enddo
    endif

  enddo

end procedure dyed_channel_update_flow
end submodule dyed_channel_initialization_s
