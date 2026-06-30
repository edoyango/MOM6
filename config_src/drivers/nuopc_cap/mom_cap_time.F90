! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> This was originally share code in CIME, but required CIME as a
!! dependency to build the MOM cap.  The options here for setting
!! a restart alarm are useful for all caps, so a second step is to
!! determine if/how these could be offered more generally in a
!! shared library.  For now we really want the MOM cap to only
!! depend on MOM and ESMF/NUOPC.
module MOM_cap_time

! !USES:
use ESMF                  , only : ESMF_Time, ESMF_Clock, ESMF_Calendar, ESMF_Alarm
use ESMF                  , only : ESMF_TimeGet, ESMF_TimeSet
use ESMF                  , only : ESMF_TimeInterval, ESMF_TimeIntervalSet
use ESMF                  , only : ESMF_ClockGet, ESMF_AlarmCreate
use ESMF                  , only : ESMF_SUCCESS, ESMF_LogWrite, ESMF_LOGMSG_INFO
use ESMF                  , only : ESMF_LogSetError, ESMF_LogFoundError, ESMF_LOGERR_PASSTHRU
use ESMF                  , only : ESMF_RC_ARG_BAD
use ESMF                  , only : operator(<), operator(/=), operator(+), operator(-), operator(*) , operator(>=)
use ESMF                  , only : operator(<=), operator(>), operator(==)
use MOM_cap_methods       , only : ChkErr

implicit none ; private

public  :: AlarmInit  ! initialize an alarm

private :: TimeInit
private :: date2ymd

! Clock and alarm options
character(len=*), private, parameter :: &
    optNONE           = "none"      , &
    optNever          = "never"     , &
    optNSteps         = "nsteps"    , &
    optNStep          = "nstep"     , &
    optNSeconds       = "nseconds"  , &
    optNSecond        = "nsecond"   , &
    optNMinutes       = "nminutes"  , &
    optNMinute        = "nminute"   , &
    optNHours         = "nhours"    , &
    optNHour          = "nhour"     , &
    optNDays          = "ndays"     , &
    optNDay           = "nday"      , &
    optNMonths        = "nmonths"   , &
    optNMonth         = "nmonth"    , &
    optNYears         = "nyears"    , &
    optNYear          = "nyear"     , &
    optMonthly        = "monthly"   , &
    optYearly         = "yearly"    , &
    optDate           = "date"      , &
    optIfdays0        = "ifdays0"   , &
    optGLCCouplingPeriod = "glc_coupling_period"

! Module data
integer, parameter          :: SecPerDay = 86400 ! Seconds per day
character(len=*), parameter :: u_FILE_u = &
    __FILE__


  interface
module subroutine AlarmInit( clock, alarm, option, &
                      opt_n, opt_ymd, opt_tod, RefTime, alarmname, rc)
  type(ESMF_Clock)            , intent(inout) :: clock     !< ESMF clock
  type(ESMF_Alarm)            , intent(inout) :: alarm     !< ESMF alarm
  character(len=*)            , intent(in)    :: option    !< alarm option
  integer          , optional , intent(in)    :: opt_n     !< alarm freq
  integer          , optional , intent(in)    :: opt_ymd   !< alarm ymd
  integer          , optional , intent(in)    :: opt_tod   !< alarm tod (sec)
  type(ESMF_Time)  , optional , intent(in)    :: RefTime   !< ref time
  character(len=*) , optional , intent(in)    :: alarmname !< alarm name
  integer                     , intent(inout) :: rc        !< Return code

  ! local variables
  !-------------------------------------------------------------------------------

end subroutine AlarmInit
module subroutine TimeInit( Time, ymd, cal, tod, desc, logunit, rc)
  type(ESMF_Time)     , intent(inout)         :: Time !< ESMF time
  integer             , intent(in)            :: ymd  !< year, month, day YYYYMMDD
  type(ESMF_Calendar) , intent(in)            :: cal  !< ESMF calendar
  integer             , intent(in),  optional :: tod  !< time of day in [sec]
  character(len=*)    , intent(in),  optional :: desc !< description of time to set
  integer             , intent(in),  optional :: logunit!< Unit for stdout output
  integer             , intent(out), optional :: rc   !< Return code

  ! local varaibles
  !-------------------------------------------------------------------------------

end subroutine TimeInit
module subroutine date2ymd (date, year, month, day)
  integer, intent(in)  :: date             !< coded-date (yyyymmdd)
  integer, intent(out) :: year,month,day   !< calendar year,month,day

  ! local variables
  !-------------------------------------------------------------------------------

end subroutine date2ymd
  end interface

end module MOM_cap_time
