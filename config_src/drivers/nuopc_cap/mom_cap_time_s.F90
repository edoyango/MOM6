submodule (MOM_cap_time) MOM_cap_time_s
  implicit none
contains
module procedure AlarmInit
  type(ESMF_Calendar)     :: cal              ! calendar
  integer                 :: lymd             ! local ymd
  integer                 :: ltod             ! local tod
  integer                 :: cyy,cmm,cdd,csec ! time info
  integer                 :: nyy,nmm,ndd,nsec ! time info
  character(len=64)       :: lalarmname       ! local alarm name
  logical                 :: update_nextalarm ! update next alarm
  type(ESMF_Time)         :: CurrTime         ! Current Time
  type(ESMF_Time)         :: NextAlarm        ! Next restart alarm time
  type(ESMF_TimeInterval) :: AlarmInterval    ! Alarm interval
  character(len=*), parameter :: subname = '(AlarmInit): '
  rc = ESMF_SUCCESS

  lalarmname = 'alarm_unknown'
  if (present(alarmname)) lalarmname = trim(alarmname)
  ltod = 0
  if (present(opt_tod)) ltod = opt_tod
  lymd = -1
  if (present(opt_ymd)) lymd = opt_ymd

  ! verify parameters
  if (trim(option) == optNSteps    .or. trim(option) == optNStep   .or. &
      trim(option) == optNSeconds  .or. trim(option) == optNSecond .or. &
      trim(option) == optNMinutes  .or. trim(option) == optNMinute .or. &
      trim(option) == optNHours    .or. trim(option) == optNHour   .or. &
      trim(option) == optNDays     .or. trim(option) == optNDay    .or. &
      trim(option) == optNMonths   .or. trim(option) == optNMonth  .or. &
      trim(option) == optNYears    .or. trim(option) == optNYear   .or. &
      trim(option) == optIfdays0) then
    if (.not. present(opt_n)) then
      call ESMF_LogSetError(ESMF_RC_ARG_BAD, &
           msg=subname//trim(option)//' requires opt_n', &
           line=__LINE__, &
           file=__FILE__, rcToReturn=rc)
      return
    endif
    if (opt_n <= 0) then
      call ESMF_LogSetError(ESMF_RC_ARG_BAD, &
           msg=subname//trim(option)//' invalid opt_n', &
           line=__LINE__, &
           file=__FILE__, rcToReturn=rc)
      return
    endif
  endif

  call ESMF_ClockGet(clock, CurrTime=CurrTime, rc=rc)
  if (ChkErr(rc,__LINE__,u_FILE_u)) return

  call ESMF_TimeGet(CurrTime, yy=cyy, mm=cmm, dd=cdd, s=csec, rc=rc )
  if (ChkErr(rc,__LINE__,u_FILE_u)) return

  call ESMF_TimeGet(CurrTime, yy=nyy, mm=nmm, dd=ndd, s=nsec, rc=rc )
  if (ChkErr(rc,__LINE__,u_FILE_u)) return

  ! initial guess of next alarm, this will be updated below
  if (present(RefTime)) then
    NextAlarm = RefTime
  else
    NextAlarm = CurrTime
  endif

  ! Determine calendar
  call ESMF_ClockGet(clock, calendar=cal, rc=rc)
  if (ChkErr(rc,__LINE__,u_FILE_u)) return

  ! Determine inputs for call to create alarm
  selectcase (trim(option))

  case (optNONE, optNever)
    call ESMF_TimeIntervalSet(AlarmInterval, yy=9999, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_TimeSet( NextAlarm, yy=9999, mm=12, dd=1, s=0, calendar=cal, rc=rc )
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    update_nextalarm  = .false.

  case (optDate)
    if (.not. present(opt_ymd)) then
      call ESMF_LogSetError(ESMF_RC_ARG_BAD, &
           msg=subname//trim(option)//' requires opt_ymd', &
           line=__LINE__, &
           file=__FILE__, rcToReturn=rc)
      return
    endif
    if (lymd < 0 .or. ltod < 0) then
      call ESMF_LogSetError(ESMF_RC_ARG_BAD, &
           msg=subname//trim(option)//'opt_ymd, opt_tod invalid', &
           line=__LINE__, &
           file=__FILE__, rcToReturn=rc)
      return
    endif
    call ESMF_TimeIntervalSet(AlarmInterval, yy=9999, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call TimeInit(NextAlarm, lymd, cal, tod=ltod, desc="optDate", rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    update_nextalarm  = .false.

  case (optIfdays0)
    if (.not. present(opt_ymd)) then
      call ESMF_LogSetError(ESMF_RC_ARG_BAD, &
           msg=subname//trim(option)//' requires opt_ymd', &
           line=__LINE__, &
           file=__FILE__, rcToReturn=rc)
      return
    endif
    call ESMF_TimeIntervalSet(AlarmInterval, mm=1, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_TimeSet( NextAlarm, yy=cyy, mm=cmm, dd=opt_n, s=0, calendar=cal, rc=rc )
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    update_nextalarm  = .true.

  case (optNSteps, optNStep)
    call ESMF_ClockGet(clock, TimeStep=AlarmInterval, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    AlarmInterval = AlarmInterval * opt_n
    update_nextalarm  = .true.

  case (optNSeconds, optNSecond)
    call ESMF_TimeIntervalSet(AlarmInterval, s=1, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    AlarmInterval = AlarmInterval * opt_n
    update_nextalarm  = .true.

  case (optNMinutes, optNMinute)
    call ESMF_TimeIntervalSet(AlarmInterval, s=60, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    AlarmInterval = AlarmInterval * opt_n
    update_nextalarm  = .true.

  case (optNHours, optNHour)
    call ESMF_TimeIntervalSet(AlarmInterval, s=3600, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    AlarmInterval = AlarmInterval * opt_n
    update_nextalarm  = .true.

  case (optNDays, optNDay)
    call ESMF_TimeIntervalSet(AlarmInterval, d=1, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    AlarmInterval = AlarmInterval * opt_n
    update_nextalarm  = .true.

  case (optNMonths, optNMonth)
    call ESMF_TimeIntervalSet(AlarmInterval, mm=1, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    AlarmInterval = AlarmInterval * opt_n
    update_nextalarm  = .true.

  case (optMonthly)
    call ESMF_TimeIntervalSet(AlarmInterval, mm=1, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_TimeSet( NextAlarm, yy=cyy, mm=cmm, dd=1, s=0, calendar=cal, rc=rc )
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    update_nextalarm  = .true.

  case (optNYears, optNYear)
    call ESMF_TimeIntervalSet(AlarmInterval, yy=1, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    AlarmInterval = AlarmInterval * opt_n
    update_nextalarm  = .true.

  case (optYearly)
    call ESMF_TimeIntervalSet(AlarmInterval, yy=1, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_TimeSet( NextAlarm, yy=cyy, mm=1, dd=1, s=0, calendar=cal, rc=rc )
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    update_nextalarm  = .true.

  case default
    call ESMF_LogSetError(ESMF_RC_ARG_BAD, &
         msg=subname//' unknown option: '//trim(option), &
         line=__LINE__, &
         file=__FILE__, rcToReturn=rc)
    return

  end select

  ! --------------------------------------------------------------------------------
  ! --- AlarmInterval and NextAlarm should be set ---
  ! --------------------------------------------------------------------------------

  ! --- advance Next Alarm so it won't ring on first timestep for
  ! --- most options above. go back one alarminterval just to be careful

  if (update_nextalarm) then
    NextAlarm = NextAlarm - AlarmInterval
    do while (NextAlarm <= CurrTime)
      NextAlarm = NextAlarm + AlarmInterval
    enddo
  endif

  alarm = ESMF_AlarmCreate( name=lalarmname, clock=clock, ringTime=NextAlarm, ringInterval=AlarmInterval, rc=rc)
  if (ChkErr(rc,__LINE__,u_FILE_u)) return

end procedure AlarmInit
module procedure TimeInit
  integer                     :: yr, mon, day ! Year, month, day as integers
  integer                     :: ltod         ! local tod
  character(len=256)          :: ldesc        ! local desc
  character(len=*), parameter :: subname = '(TimeInit) '
  ltod = 0
  if (present(tod)) ltod = tod
  ldesc = ''
  if (present(desc)) ldesc = desc

  if ( (ymd < 0) .or. (ltod < 0) .or. (ltod > SecPerDay) )then
    if (present(logunit)) then
      write(logunit,*) subname//': ERROR yymmdd is a negative number or '// &
            'time-of-day out of bounds', ymd, ltod
    endif
    call ESMF_LogSetError(ESMF_RC_ARG_BAD, &
         msg=subname//' yymmdd is negative or time-of-day out of bounds ', &
         line=__LINE__, &
         file=__FILE__, rcToReturn=rc)
    return
  endif

  call date2ymd (ymd,yr,mon,day)

  call ESMF_TimeSet( Time, yy=yr, mm=mon, dd=day, s=ltod, calendar=cal, rc=rc )
  if (ChkErr(rc,__LINE__,u_FILE_u)) return

end procedure TimeInit
module procedure date2ymd
  integer :: tdate   ! temporary date
  character(*),parameter :: subName = "(date2ymd)"
  tdate = abs(date)
  year = int(tdate/10000)
  if (date < 0) then
    year = -year
  endif
  month = int( mod(tdate,10000)/  100)
  day = mod(tdate,  100)

end procedure date2ymd
end submodule MOM_cap_time_s
