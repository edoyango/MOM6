! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> Set of time utilities for converting between FMS and ESMF time type.
module time_utils_mod

! FMS
use fms_mod,            only: uppercase
use mpp_mod,            only: mpp_error, FATAL
use time_manager_mod,   only: time_type, set_time, set_date, get_date
use time_manager_mod,   only: GREGORIAN, JULIAN, NOLEAP, THIRTY_DAY_MONTHS, NO_CALENDAR
use time_manager_mod,   only: fms_get_calendar_type => get_calendar_type
! ESMF
use ESMF,               only: ESMF_CALKIND_FLAG, ESMF_CALKIND_GREGORIAN
use ESMF,               only: ESMF_CALKIND_JULIAN, ESMF_CALKIND_NOLEAP
use ESMF,               only: ESMF_CALKIND_360DAY, ESMF_CALKIND_NOCALENDAR
use ESMF,               only: ESMF_Time, ESMF_TimeGet, ESMF_LogFoundError
use ESMF,               only: ESMF_LOGERR_PASSTHRU,ESMF_TimeInterval
use ESMF,               only: ESMF_TimeIntervalGet, ESMF_TimeSet, ESMF_SUCCESS
use MOM_cap_methods,    only: ChkErr

implicit none ; private

!> Converts calendar from FMS to ESMF format
interface fms2esmf_cal
  module procedure fms2esmf_cal_c
  module procedure fms2esmf_cal_i
end interface fms2esmf_cal

!> Converts time from FMS to ESMF format
interface esmf2fms_time
  module procedure esmf2fms_time_t
  module procedure esmf2fms_timestep
end interface esmf2fms_time

public fms2esmf_cal
public esmf2fms_time
public fms2esmf_time
public string_to_date

character(len=*),parameter :: u_FILE_u = &
     __FILE__


  interface
module function fms2esmf_cal_c(calendar)
  type(ESMF_CALKIND_FLAG)            :: fms2esmf_cal_c !< ESMF calendar type
  character(len=*), intent(in)       :: calendar       !< Type of calendar

end function fms2esmf_cal_c
module function fms2esmf_cal_i(calendar)
  type(ESMF_CALKIND_FLAG)            :: fms2esmf_cal_i !< ESMF calendar structure
  integer, intent(in)                :: calendar       !< Type of calendar

end function fms2esmf_cal_i
module function esmf2fms_time_t(time)
  type(Time_type)                    :: esmf2fms_time_t !< FMS time structure
  type(ESMF_Time), intent(in)        :: time            !< ESMF time structure

  ! Local Variables


end function esmf2fms_time_t
module function esmf2fms_timestep(timestep)
  type(Time_type)                    :: esmf2fms_timestep !< FMS time structure
  type(ESMF_TimeInterval), intent(in):: timestep          !< time-interval following
                                                          !! ESMF format [s]
  ! Local Variables


end function esmf2fms_timestep
module function fms2esmf_time(time, calkind)
  type(ESMF_Time)                    :: fms2esmf_time !< ESMF time structure
  type(time_type), intent(in)        :: time          !< FMS time structure
  type(ESMF_CALKIND_FLAG), intent(in), optional :: calkind !< ESMF calendar structure

  ! Local Variables


end function fms2esmf_time
module function string_to_date(string, rc)
  character(len=15), intent(in)           :: string        !< String representing a date
  integer, intent(out), optional          :: rc            !< ESMF error handler
  type(time_type)                         :: string_to_date!< FMS time structure

  ! Local variables

end function string_to_date
  end interface

end module time_utils_mod
