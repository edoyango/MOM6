submodule (time_utils_mod) time_utils_mod_s
  implicit none
contains
module procedure fms2esmf_cal_c
  select case( uppercase(trim(calendar)) )
  case( 'GREGORIAN' )
    fms2esmf_cal_c = ESMF_CALKIND_GREGORIAN
  case( 'JULIAN' )
    fms2esmf_cal_c = ESMF_CALKIND_JULIAN
  case( 'NOLEAP' )
    fms2esmf_cal_c = ESMF_CALKIND_NOLEAP
  case( 'THIRTY_DAY' )
    fms2esmf_cal_c = ESMF_CALKIND_360DAY
  case( 'NO_CALENDAR' )
    fms2esmf_cal_c = ESMF_CALKIND_NOCALENDAR
  case default
    call mpp_error(FATAL, &
    'ocean_solo: ocean_solo_nml entry calendar must be one of GREGORIAN|JULIAN|NOLEAP|THIRTY_DAY|NO_CALENDAR.' )
  end select
end procedure fms2esmf_cal_c
module procedure fms2esmf_cal_i
  select case(calendar)
    case(THIRTY_DAY_MONTHS)
      fms2esmf_cal_i = ESMF_CALKIND_360DAY
    case(GREGORIAN)
      fms2esmf_cal_i = ESMF_CALKIND_GREGORIAN
    case(JULIAN)
      fms2esmf_cal_i = ESMF_CALKIND_JULIAN
    case(NOLEAP)
      fms2esmf_cal_i = ESMF_CALKIND_NOLEAP
    case(NO_CALENDAR)
      fms2esmf_cal_i = ESMF_CALKIND_NOCALENDAR
  end select
end procedure fms2esmf_cal_i
module procedure esmf2fms_time_t
  integer                            :: yy, mm, dd, h, m, s
  type(ESMF_CALKIND_FLAG)            :: calkind
  integer                            :: rc
  call ESMF_TimeGet(time, yy=yy, mm=mm, dd=dd, h=h, m=m, s=s, &
      calkindflag=calkind, rc=rc)
  if (ChkErr(rc,__LINE__,u_FILE_u)) return

  esmf2fms_time_t = set_date(yy, mm, dd, h, m, s)

end procedure esmf2fms_time_t
module procedure esmf2fms_timestep
  integer                            :: s
  type(ESMF_CALKIND_FLAG)            :: calkind
  integer                            :: rc
  call ESMF_TimeIntervalGet(timestep, s=s, calkindflag=calkind, rc=rc)
  if (ChkErr(rc,__LINE__,u_FILE_u)) return

  esmf2fms_timestep = set_time(s, 0)

end procedure esmf2fms_timestep
module procedure fms2esmf_time
  integer                            :: yy, mm, d, h, m, s
  type(ESMF_CALKIND_FLAG)            :: l_calkind
  integer                            :: rc
  if (present(calkind)) then
    l_calkind = calkind
  else
    l_calkind = fms2esmf_cal(fms_get_calendar_type())
  endif

  call get_date(time, yy, mm, d, h, m, s)

  call ESMF_TimeSet(fms2esmf_time, yy=yy, mm=mm, d=d, h=h, m=m, s=s, &
      calkindflag=l_calkind, rc=rc)
  if (ChkErr(rc,__LINE__,u_FILE_u)) return

end procedure fms2esmf_time
module procedure string_to_date
  integer                                 :: yr,mon,day,hr,min,sec
  if (present(rc)) rc = ESMF_SUCCESS

  read(string, '(I4.4,I2.2,I2.2,".",I2.2,I2.2,I2.2)') yr, mon, day, hr, min, sec
  string_to_date = set_date(yr, mon, day, hr, min, sec)

end procedure string_to_date
end submodule time_utils_mod_s
