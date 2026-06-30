submodule (MOM_cap_outputlog) MOM_cap_outputlog_s
  implicit none
contains
module procedure outputlog_init
  rc = ESMF_SUCCESS
end procedure outputlog_init
module procedure outputlog_run
  rc = ESMF_SUCCESS
end procedure outputlog_run
module procedure outputlog_restart
  rc = ESMF_SUCCESS
end procedure outputlog_restart
module procedure outputlog_init
  type(ESMF_Time)         :: mcurrTime
  type(ESMF_TimeInterval) :: alarmoffset
  logical                 :: isPresent, isSet
  integer                 :: n
  integer                 :: year, month, day, hour
  character(len=3)        :: chour
  character(len=256)      :: msgString
  character(len=256)      :: value
  character(len=256)      :: subname='MOM_cap:(outputlog_init)'
  rc = ESMF_SUCCESS
  call ESMF_GridCompGet(gcomp, vm=vm, rc=rc)
  if (ChkErr(rc,__LINE__,u_FILE_u)) return

  call NUOPC_CompAttributeGet(gcomp, name="mom6_restart_dir", value=value, &
       isPresent=isPresent, isSet=isSet, rc=rc)
  if (ChkErr(rc,__LINE__,u_FILE_u)) return
  if (isPresent .and. isSet) then
    restartdir = trim(value)
  else
    restartdir = './'
  endif
  if (restartdir(len_trim(restartdir):len_trim(restartdir)) /= '/') then
    restartdir = trim(restartdir)//'/'
  endif
  write(msgString,'(A)')'MOM_cap:MOM6 restart directory = '//trim(restartdir)
  call ESMF_LogWrite(trim(msgString), ESMF_LOGMSG_INFO)

  call NUOPC_CompAttributeGet(gcomp, name="mom6_output_dir", value=value, &
       isPresent=isPresent, isSet=isSet, rc=rc)
  if (ChkErr(rc,__LINE__,u_FILE_u)) return
  if (isPresent .and. isSet) then
    outputdir = trim(value)
  else
    outputdir = './'
  endif
  if (outputdir(len_trim(outputdir):len_trim(outputdir)) /= '/') then
    outputdir = trim(outputdir)//'/'
  endif
  write(msgString,'(A)')'MOM_cap:MOM6 output directory = '//trim(outputdir)
  call ESMF_LogWrite(trim(msgString), ESMF_LOGMSG_INFO)

  call NUOPC_CompAttributeGet(gcomp, name="mom6_output_fh", value=value, &
       isPresent=isPresent, isSet=isSet, rc=rc)
  if (ChkErr(rc,__LINE__,u_FILE_u)) return
  if (isPresent .and. isSet) then
    if (len_trim(value) == 1) then
      output_fh = '0'//trim(value)
    else
      output_fh = trim(value)
    endif
  else
    output_fh = '06'
  endif
  write(msgString,'(A)')'MOM_cap:MOM6 output frequency = '//trim(output_fh)
  call ESMF_LogWrite(trim(msgString), ESMF_LOGMSG_INFO)

  debug = .false.
  call NUOPC_CompAttributeGet(gcomp, name="debug_outputlog", value=value, &
       isPresent=isPresent, isSet=isSet, rc=rc)
  if (ChkErr(rc,__LINE__,u_FILE_u)) return
  if (isPresent .and. isSet) debug=(trim(value)=="true")
  if (debug) call ESMF_LogWrite('MOM_cap:MOM6 output debug ON', ESMF_LOGMSG_INFO)

  call ESMF_ClockGet(mclock, currTime=mcurrTime, rc=rc)
  if (ChkErr(rc,__LINE__,u_FILE_u)) return
  call ESMF_TimeIntervalSet(tincrement, m=1, rc=rc)
  if (ChkErr(rc,__LINE__,u_FILE_u)) return

  ! get start hour time offset (ie, fhrot)
  call ESMF_TimeGet(mcurrTime, yy=year, mm=month, dd=day, h=hour, rc=rc)
  if (ChkErr(rc,__LINE__,u_FILE_u)) return
  if (mod(hour,6) /= 0) then
    toffset = hour - 6
  else
    toffset = 0
  endif
  if (debug .and. is_root_pe()) then
    print '(A,i8)',trim(subname)//' toffset = ',toffset
  endif
  ! initialize
  lastrestart = mcurrTime

  do n = 1,n_freq
    write(chour,'(I2.2,A)')freq(n),'h'
    olog(n)%alarm_name          = 'output_alarm'//trim(chour)
    olog(n)%opt_n               = freq(n)
    olog(n)%chkfile_nextAdvance = .false.
    olog(n)%use_filesize        = .false.
    olog(n)%filename            = ''
    olog(n)%createsize            = 0
    olog(n)%time_lastrestart    = lastrestart
    olog(n)%fhoffset            = 60*freq(n)*tincrement
    olog(n)%filename_fhoffset   = 90*freq(n)*tincrement

    ! the time offset in hours required to ensure the alarm rings at multiples of 6
    if (freq(n) >= 6) then
      alarmoffset = toffset*60*tincrement
    else
      alarmoffset = 0*tincrement
    endif

    call AlarmInit(mclock,                  &
         alarm     = olog(n)%alarm,         &
         option    = 'nhours',              &
         opt_n     = olog(n)%opt_n,         &
         opt_ymd   = -999,                  &
         RefTime   = mcurrTime+alarmoffset, &
         alarmname = olog(n)%alarm_name, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_AlarmSet(olog(n)%alarm, clock=mclock, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    write(msgString,'(A)')trim(subname)//' Output alarm '//trim(olog(n)%alarm_name)//' Created & Set'
    call ESMF_LogWrite(trim(msgString), ESMF_LOGMSG_INFO)
    if (debug .and. is_root_pe()) then
      call ESMF_TimeIntervalPrint(olog(n)%filename_fhoffset, options="string", rc=rc)
      if (ChkErr(rc,__LINE__,u_FILE_u)) return
    endif
  enddo
end procedure outputlog_init
module procedure outputlog_run
  type(ESMF_Time)    :: nextTime, currTime, startTime, prevRing
  logical            :: lstop
  logical            :: filecomplete
  integer            :: n, nlen(1), fsize(1)
  character(len=3)   :: chour
  character(len=40)  :: importexport
  character(len=16)  :: timestr
  character(len=256) :: fname
  character(len=256) :: subname='MOM_cap:(outputlog_run)'
  rc = ESMF_SUCCESS

  call ESMF_ClockGet(mclock, startTime=startTime, currTime=currTime, rc=rc)
  if (ChkErr(rc,__LINE__,u_FILE_u)) return
  call ESMF_ClockGetNextTime(mclock, nextTime, rc=rc)
  if (ChkErr(rc,__LINE__,u_FILE_u)) return
  importexport = get_importexport(currTime, nextTime, rc=rc)
  if (ChkErr(rc,__LINE__,u_FILE_u)) return

  lstop = .false.
  if (present(atStopTime)) then
    lstop = atStopTime
  endif

  filecomplete = .false.
  fsize(1) = nf90_fill_int
  nlen(1)  = nf90_fill_int

  do n = 1,n_freq
    write(chour,'(I2.2,A)')freq(n),'h'
    if (chour(1:2) == output_fh(1:2)) then
      call ESMF_ClockGetAlarm(mclock, alarmname=trim(olog(n)%alarm_name), alarm=olog(n)%alarm, rc=rc)
      if (ChkErr(rc,__LINE__,u_FILE_u)) return
      ! when the alarm rings, set file check on next advance and construct the filename
      if (ESMF_AlarmIsRinging(olog(n)%alarm, rc=rc)) then
        if (ChkErr(rc,__LINE__,u_FILE_u)) return
        call ESMF_AlarmRingerOff(olog(n)%alarm, rc=rc )
        if (ChkErr(rc,__LINE__,u_FILE_u)) return
        olog(n)%chkfile_nextAdvance = .true.

        timestr = get_timestr(nextTime-olog(n)%filename_fhoffset, rc=rc)
        if (ChkErr(rc,__LINE__,u_FILE_u)) return
        write(olog(n)%filename,'(A)')trim(outputdir)//'ocn_'//trim(timestr)//'.nc'

        fname = trim(olog(n)%filename)
        inquire(file=fname, exist=existflag)
        if (existflag) then
          if (is_root_pe()) then
            nlen(1) = get_unlimited_len(trim(fname))
            inquire(file=fname, size=fsize(1))
          endif
          call ESMF_VMBroadCast(vm, nlen, 1, 0, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          call ESMF_VMBroadCast(vm, fsize, 1, 0, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          olog(n)%createsize = fsize(1)

          if (nlen(1) == 0) then
            olog(n)%use_filesize = .false.
          else
            olog(n)%use_filesize = .true.
          endif
        endif
        if (debug .and. is_root_pe()) then
          print '(A,2(A,L),A,2i16)',trim(subname)//' fname '//trim(olog(n)%filename)//'  '//trim(importexport), &
               ' checkflag ',olog(n)%chkfile_nextAdvance,' use_filesize ',olog(n)%use_filesize,                 &
               '  ',olog(n)%createsize,nlen(1)
        endif
      endif

      if (olog(n)%chkfile_nextAdvance) then
        fname = trim(olog(n)%filename)
        filecomplete = file_is_complete(fname, olog(n)%use_filesize, olog(n)%createsize, rc)
        if (ChkErr(rc,__LINE__,u_FILE_u)) return

        if (filecomplete) then
          olog(n)%chkfile_nextAdvance = .false.
          olog(n)%time_lastrestart = lastrestart
          if (is_root_pe()) then
            call log_restart_fh(currTime-olog(n)%fhoffset, startTime, 'mom6.'//chour, prefixtime=.true., &
                 lastrestart=olog(n)%time_lastrestart, lastoutput=olog(n)%filename, rc=rc)
            if (ChkErr(rc,__LINE__,u_FILE_u)) return
          endif
        endif
      endif
      if (debug .and. is_root_pe()) call debug_info(trim(subname)//'  ',trim(olog(n)%filename), &
           olog(n)%chkfile_nextAdvance, olog(n)%createsize, importexport)

      if (lstop) then
        ! use prevRing in place of currTime to allow for stopping between averaging intervals
        ! prevring == currTime if stopping on intervals
        call ESMF_AlarmGet(olog(n)%alarm, prevRingTime=prevring, rc=rc)
        if (ChkErr(rc,__LINE__,u_FILE_u)) return

        timestr = get_timestr(prevring-30*freq(n)*tincrement, rc=rc)
        if (ChkErr(rc,__LINE__,u_FILE_u)) return
        write(olog(n)%filename,'(A)')trim(outputdir)//'ocn_'//trim(timestr)//'.nc'

        fname = trim(olog(n)%filename)
        filecomplete = file_is_complete(fname, olog(n)%use_filesize, olog(n)%createsize, rc)
        if (ChkErr(rc,__LINE__,u_FILE_u)) return

        if (filecomplete) then
          olog(n)%chkfile_nextAdvance = .false.
          olog(n)%time_lastrestart = lastrestart
          if (is_root_pe()) then
            call log_restart_fh(prevring, startTime, 'mom6.lstop.'//chour, prefixtime=.true., &
                 lastrestart=olog(n)%time_lastrestart, lastoutput=olog(n)%filename, rc=rc)
            if (ChkErr(rc,__LINE__,u_FILE_u)) return
          endif
        endif
        if (debug .and. is_root_pe()) call debug_info(trim(subname)//' lstop ',trim(olog(n)%filename), &
             olog(n)%chkfile_nextAdvance, olog(n)%createsize, importexport)

      endif ! lstop
    endif ! chour = output_fh
  enddo
end procedure outputlog_run
module procedure outputlog_restart
  type(ESMF_Time)      :: startTime, currTime, nextTime
  integer              :: n, nlen(1)
  integer              :: year, month, day, hour, minute, seconds
  character(len=256)   :: fname
  character(len=15)    :: timestr
  character(len=40)    :: importexport
  logical, allocatable :: allDone(:)
  character(len=8)     :: suffix
  character(len=256)   :: subname='MOM_cap:(outputlog_restart)'
  rc = ESMF_SUCCESS

  call ESMF_ClockGet(mclock, startTime=startTime, currTime=currTime, rc=rc)
  if (ChkErr(rc,__LINE__,u_FILE_u)) return
  call ESMF_ClockGetNextTime(mclock, nextTime, rc=rc)
  if (ChkErr(rc,__LINE__,u_FILE_u)) return
  importexport = get_importexport(currTime, nextTime, rc=rc)
  if (ChkErr(rc,__LINE__,u_FILE_u)) return

  call ESMF_TimeGet(nextTime, yy=year, mm=month, dd=day, h=hour, m=minute, s=seconds, rc=rc )
  if (ChkErr(rc,__LINE__,u_FILE_u)) return
  write(timestr,'(I4.4,2(I2.2),A,3(I2.2))') year, month, day,".", hour, minute, seconds

  allocate(allDone(1:num_rest_files))
  allDone = .false.

  do n = 1,num_rest_files
    if (n == 1) then
      suffix = ''
    else if (n-1 < 10) then
      write(suffix,'("_",I1)') n-1
    else
      write(suffix,'("_",I2)') n-1
    endif
    if (len_trim(suffix) == 0) then
      fname = trim(restartdir)//trim(timestr)//'.MOM.res.nc'
    else
      fname = trim(restartdir)//trim(timestr)//'.MOM.res'//trim(suffix)//'.nc'
    endif

    ! check if file is written
    inquire(file=trim(fname), exist=existflag)
    if (existflag) then
      if (is_root_pe())then
        nlen(1) = get_unlimited_len(trim(fname))
      endif
      call ESMF_VMBroadCast(vm, nlen, 1, 0, rc=rc)
      if (ChkErr(rc,__LINE__,u_FILE_u)) return

      if (nlen(1) > 0) allDone(n) = .true.
      if (debug .and. is_root_pe()) then
        if (nlen(1) > 0) then
          print '(A)',trim(subname)//' restart '//trim(fname)//'  '//trim(importexport)//' complete'
        else
          print '(A)',trim(subname)//' restart '//trim(fname)//'  '//trim(importexport)//' still 0'
        endif
      endif
    endif
  enddo ! num_rest_files

  if (all(allDone) .eqv. .true.) then
    lastrestart = nextTime
    if (is_root_pe()) then
      call log_restart_fh(nextTime, startTime, 'mom6.res', prefixtime=.true., rc=rc)
      if (ChkErr(rc,__LINE__,u_FILE_u)) return
    endif
  endif
end procedure outputlog_restart
module procedure file_is_complete
  integer :: nlen(1), fsize(1)
  rc = ESMF_SUCCESS

  filecomplete = .false.
  nlen(1) = nf90_fill_int
  fsize(1) = nf90_fill_int

  inquire(file=fname, exist=existflag)
  if (existflag) then
    if (is_root_pe()) then
      nlen(1) = get_unlimited_len(fname)
      inquire(file=fname, size=fsize(1))
    endif
    call ESMF_VMBroadCast(vm, nlen, 1, 0, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_VMBroadCast(vm, fsize, 1, 0, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
  endif

  if (chk4size) then
    filecomplete = (nlen(1) > 0 .and. fsize(1) > createsize)
  else
    filecomplete = (nlen(1) > 0)
  endif
end procedure file_is_complete
module procedure get_unlimited_len
  integer :: ncid, dimid
  unlen = 0
  call nf90_err(nf90_open(trim(fname), nf90_nowrite, ncid), 'nf90_open: '//trim(fname))
  call nf90_err(nf90_inquire(ncid, unlimiteddimid=dimid), 'inquire unlimiteddimid')
  call nf90_err(nf90_inquire_dimension(ncid, dimid, len=unlen), 'inquire unlimited dimension')
  call nf90_err(nf90_close(ncid), 'close: '//trim(fname))
end procedure get_unlimited_len
module procedure get_timestr
  integer :: year, month, day, hour, minute
  rc = ESMF_SUCCESS

  call ESMF_TimeGet(MyTime, yy=year, mm=month, dd=day, h=hour, m=minute, rc=rc)
  if (ChkErr(rc,__LINE__,u_FILE_u)) return
  write(timestr,'(I4.4,4(A,I2.2))')year,'_',month,'_',day,'_',hour,'_',minute
end procedure get_timestr
module procedure get_importexport
  character(len=19) :: import_timestr, export_timestr
  rc = ESMF_SUCCESS

  call ESMF_TimeGet(currTime, timestring=import_timestr, rc=rc)
  if (ChkErr(rc,__LINE__,u_FILE_u)) return
  call ESMF_TimeGet(nextTime, timestring=export_timestr, rc=rc)
  if (ChkErr(rc,__LINE__,u_FILE_u)) return
  importexport = trim(import_timestr)//'  '//trim(export_timestr)
end procedure get_importexport
module procedure debug_info
  integer :: fsize
  character(len=256) :: msgString
  inquire(file=fname, exist=existflag)
  if (existflag) then
    inquire(file=fname, size=fsize)
    write(msgString,'(A)')tag//'  '//fname//' exists '//timestring
    if (chkflag) then
      print '(A,L,2i16)',trim(msgString)//' not complete, chkflag ',chkflag,filesize,fsize
    else
      print '(A,L,2i16)',trim(msgString)//'     complete, chkflag ',chkflag,filesize,fsize
    endif
  else
    write(msgString,'(A)')tag//'  '//fname//' does not exist '//timestring
    print '(A)',trim(msgString)
  endif
end procedure debug_info
module procedure nf90_err
  if (ierr /= nf90_noerr) then
    write(0, '(A)') 'FATAL ERROR: ' // trim(string)// ' : ' // trim(nf90_strerror(ierr))
    ! This fails on WCOSS2 with Intel 19 compiler. See https://community.intel.com/
    ! Search term "STOP and ERROR STOP with variable stop codes"
    ! When WCOSS2 moves to Intel 2020+, uncomment the next line and remove stop 99
    !stop ierr
    stop 99
  endif
end procedure nf90_err
end submodule MOM_cap_outputlog_s
