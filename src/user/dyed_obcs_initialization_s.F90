submodule (dyed_obcs_initialization) dyed_obcs_initialization_s
#include <MOM_memory.h>
  implicit none
contains
module procedure dyed_obcs_set_OBC_data
  character(len=40)  :: mdl = "dyed_obcs_set_OBC_data" ! This subroutine's name.
  character(len=80)  :: name, longname
  integer :: is, ie, js, je, isd, ied, jsd, jed, m, n, nz, ntr_id
  integer :: IsdB, IedB, JsdB, JedB
  integer :: n_dye ! Number of regionsl dye tracers
  real :: dye ! Inflow dye concentration [arbitrary]
  type(tracer_type), pointer      :: tr_ptr => NULL()
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  IsdB = G%IsdB ; IedB = G%IedB ; JsdB = G%JsdB ; JedB = G%JedB

  if (.not.associated(OBC)) return

  call get_param(param_file, mdl, "NUM_DYED_TRACERS", ntr, &
                 "The number of dyed_obc tracers in this run. Each tracer "//&
                 "should have a separate boundary segment.  "//&
                 "If not present, use NUM_DYE_TRACERS.", default=-1, do_not_log=.true.)
  if (ntr == -1) then
    !for backward compatibility
    call get_param(param_file, mdl, "NUM_DYE_TRACERS", ntr, &
                   "The number of dye tracers in this run. Each tracer "//&
                   "should have a separate boundary segment.", default=0, do_not_log=.true.)
    n_dye = 0
  else
    call get_param(param_file, mdl, "NUM_DYE_TRACERS", n_dye, &
                   "The number of dye tracers in this run. Each tracer "//&
                   "should have a separate region.", default=0, do_not_log=.true.)
  endif

  call get_param(param_file, mdl, "DYE_OBC_INFLOW", dye_obc_inflow, &
                 "The OBC inflow value of dye tracers.", units="kg kg-1", &
                 default=1.0)

  if (OBC%number_of_segments < ntr) then
    call MOM_error(WARNING, "Error in dyed_obc segment setup")
    return   !!! Need a better error message here
  endif

! ! Set the inflow values of the dyes, one per segment.
! ! We know the order: north, south, east, west
  do m=1,ntr
    write(name,'("dye_",I2.2)') m+n_dye  !after regional dye tracers
    write(longname,'("Concentration of dyed_obc Tracer ",I2.2, " on segment ",I2.2)') m, m
    call tracer_name_lookup(tr_Reg, ntr_id, tr_ptr, name)

    do n=1,OBC%number_of_segments
      if (n == m) then
        dye = dye_obc_inflow
      else
        dye = 0.0
      endif
      call register_segment_tracer(tr_ptr, ntr_id, param_file, GV, &
                                   OBC%segment(n), OBC_scalar=dye)
    enddo
  enddo

end procedure dyed_obcs_set_OBC_data
end submodule dyed_obcs_initialization_s
