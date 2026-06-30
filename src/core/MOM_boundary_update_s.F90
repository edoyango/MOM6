submodule (MOM_boundary_update) MOM_boundary_update_s
#include <MOM_memory.h>
  implicit none
contains
module procedure call_OBC_register
  logical :: debug
  character(len=200) :: config
  character(len=40)  :: mdl = "MOM_boundary_update" ! This module's name.
# include "version_variable.h"
  if (associated(CS)) then
    call MOM_error(WARNING, "call_OBC_register called with an associated "// &
                            "control structure.")
    return
  else ; allocate(CS) ; endif

  call log_version(param_file, mdl, version, "")

  call get_param(param_file, mdl, "OBC_VALUE_UPDATE_BUG", CS%value_update_bug, &
                 "If true, recover a bug that OBC segment data does not update if all segments "//&
                 "use 'value' and none uses 'file'.", default=.true.)
  call get_param(param_file, mdl, "USE_FILE_OBC", CS%use_files, &
                 "If true, use external files for the open boundary.", &
                 default=.false.)
  call get_param(param_file, mdl, "USE_TIDAL_BAY_OBC", CS%use_tidal_bay, &
                 "If true, use the tidal_bay open boundary.", &
                 default=.false.)
  call get_param(param_file, mdl, "USE_KELVIN_WAVE_OBC", CS%use_Kelvin, &
                 "If true, use the Kelvin wave open boundary.", &
                 default=.false.)
  call get_param(param_file, mdl, "USE_SHELFWAVE_OBC", CS%use_shelfwave, &
                 "If true, use the shelfwave open boundary.", &
                 default=.false.)
  call get_param(param_file, mdl, "USE_DYED_CHANNEL_OBC", CS%use_dyed_channel, &
                 "If true, use the dyed channel open boundary.", &
                 default=.false.)
  call get_param(param_file, mdl, "OBC_USER_CONFIG", config, &
               "A string that sets how the user code is invoked to set open boundary data: \n"//&
               "   DOME - specified inflow on northern boundary\n"//&
               "   dyed_channel - supercritical with dye on the inflow boundary\n"//&
               "   dyed_obcs - circle_obcs with dyes on the open boundaries\n"//&
               "   Kelvin - barotropic Kelvin wave forcing on the western boundary\n"//&
               "   shelfwave - Flather with shelf wave forcing on western boundary\n"//&
               "   supercritical - now only needed here for the allocations\n"//&
               "   tidal_bay - Flather with tidal forcing on eastern boundary\n"//&
               "   USER - user specified", default="none", do_not_log=.true.)
  call get_param(param_file, mdl, "DEBUG", debug, &
                 "If true, write out verbose debugging data.", &
                 default=.false., debuggingParam=.true.)
  call get_param(param_file, mdl, "DEBUG_OBCS", CS%debug_OBCs, &
                 "If true, write out verbose debugging data about OBCs.", &
                 default=.false., debuggingParam=.true.)
  call get_param(param_file, mdl, "NK_OBC_DEBUG", CS%nk_OBC_debug, &
                 "The number of layers of OBC segment data to write out in full "//&
                 "when DEBUG_OBCS is true.", &
                 default=0, debuggingParam=.true., do_not_log=.not.CS%debug_OBCs)

  if (CS%use_files) CS%use_files = &
    register_file_OBC(param_file, CS%file_OBC_CSp, US, &
               OBC%OBC_Reg)

  if (trim(config) == "DOME") then
    call register_DOME_OBC(param_file, US, OBC, tr_Reg)
!  elseif (trim(config) == "tidal_bay") then
!  elseif (trim(config) == "Kelvin") then
!  elseif (trim(config) == "shelfwave") then
!  elseif (trim(config) == "dyed_channel") then
  endif

  if (CS%use_tidal_bay) CS%use_tidal_bay = &
    register_tidal_bay_OBC(param_file, CS%tidal_bay_OBC, US, &
               OBC%OBC_Reg)
  if (CS%use_Kelvin) CS%use_Kelvin = &
    register_Kelvin_OBC(param_file, CS%Kelvin_OBC_CSp, US, &
               OBC%OBC_Reg)
  if (CS%use_shelfwave) CS%use_shelfwave = &
    register_shelfwave_OBC(param_file, CS%shelfwave_OBC_CSp, G, US, &
               OBC%OBC_Reg)
  if (CS%use_dyed_channel) CS%use_dyed_channel = &
    register_dyed_channel_OBC(param_file, CS%dyed_channel_OBC_CSp, US, &
               OBC%OBC_Reg)

end procedure call_OBC_register
module procedure update_OBC_data
  if (CS%use_tidal_bay) &
      call tidal_bay_set_OBC_data(OBC, CS%tidal_bay_OBC, G, GV, US, h, Time)
  if (CS%use_Kelvin)  &
      call Kelvin_set_OBC_data(OBC, CS%Kelvin_OBC_CSp, G, GV, US, h, Time)
  if (CS%use_shelfwave) &
      call shelfwave_set_OBC_data(OBC, CS%shelfwave_OBC_CSp, G, GV, US, h, Time)
  if (CS%use_dyed_channel) &
      call dyed_channel_update_flow(OBC, CS%dyed_channel_OBC_CSp, G, GV, US, h, Time)

  if (.not. OBC%user_BCs_set_globally) then
    if (OBC%any_needs_IO_for_data) call read_OBC_segment_data(G, GV, US, OBC, tv, h, Time)
    if ((.not.CS%value_update_bug) .or. (OBC%any_needs_IO_for_data .or. OBC%add_tide_constituents)) &
      call update_OBC_segment_data(G, GV, US, OBC, h, Time)
  endif

  if (CS%debug_OBCs) call chksum_OBC_segments(OBC, G, GV, US, CS%nk_OBC_debug)

end procedure update_OBC_data
module procedure OBC_register_end
  if (CS%use_files) call file_OBC_end(CS%file_OBC_CSp)
  if (CS%use_Kelvin) call Kelvin_OBC_end(CS%Kelvin_OBC_CSp)

  if (associated(CS)) deallocate(CS)
end procedure OBC_register_end
end submodule MOM_boundary_update_s
