submodule (user_initialization) user_initialization_s
#include <MOM_memory.h>
  implicit none
contains
module procedure USER_set_coord
  call MOM_error(FATAL, &
    "USER_initialization.F90, USER_set_coord: " // &
    "Unmodified user routine called - you must edit the routine to use it")
  Rlay(:) = 0.0
  g_prime(:) = 0.0

  if (first_call) call write_user_log(param_file)

end procedure USER_set_coord
module procedure USER_initialize_topography
  call MOM_error(FATAL, &
    "USER_initialization.F90, USER_initialize_topography: " // &
    "Unmodified user routine called - you must edit the routine to use it")

  D(:,:) = 0.0

  if (first_call) call write_user_log(param_file)

end procedure USER_initialize_topography
module procedure USER_initialize_thickness
  call MOM_error(FATAL, &
    "USER_initialization.F90, USER_initialize_thickness: " // &
    "Unmodified user routine called - you must edit the routine to use it")

  if (just_read) return ! All run-time parameters have been read, so return.

  h(:,:,1:GV%ke) = 0.0 ! h should be set in [Z ~> m].  It will be converted to thickness units
                       ! [H ~> m or kg m-2] once the temperatures and salinities are known.

  if (first_call) call write_user_log(param_file)

end procedure USER_initialize_thickness
module procedure USER_initialize_velocity
  call MOM_error(FATAL, &
    "USER_initialization.F90, USER_initialize_velocity: " // &
    "Unmodified user routine called - you must edit the routine to use it")

  if (just_read) return ! All run-time parameters have been read, so return.

  u(:,:,1) = 0.0
  v(:,:,1) = 0.0

  if (first_call) call write_user_log(param_file)

end procedure USER_initialize_velocity
module procedure USER_init_temperature_salinity
  call MOM_error(FATAL, &
    "USER_initialization.F90, USER_init_temperature_salinity: " // &
    "Unmodified user routine called - you must edit the routine to use it")

  if (just_read) return ! All run-time parameters have been read, so return.

  T(:,:,1) = 0.0
  S(:,:,1) = 0.0

  if (first_call) call write_user_log(param_file)

end procedure USER_init_temperature_salinity
module procedure USER_initialize_sponges
  call MOM_error(FATAL, &
    "USER_initialization.F90, USER_initialize_sponges: " // &
    "Unmodified user routine called - you must edit the routine to use it")

  if (first_call) call write_user_log(param_file)

end procedure USER_initialize_sponges
module procedure USER_set_OBC_data
  if (first_call) call write_user_log(param_file)

end procedure USER_set_OBC_data
module procedure USER_set_rotation
  call MOM_error(FATAL, &
    "USER_initialization.F90, USER_set_rotation: " // &
    "Unmodified user routine called - you must edit the routine to use it")

  if (first_call) call write_user_log(param_file)

end procedure USER_set_rotation
module procedure write_user_log
# include "version_variable.h"
  character(len=40)  :: mdl = "user_initialization" ! This module's name.
  call log_version(param_file, mdl, version)
  first_call = .false.

end procedure write_user_log
end submodule user_initialization_s
