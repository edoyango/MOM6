! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> This module contains a set of subroutines that are required by NUOPC.

module MOM_cap_mod

use MOM_domains,              only: get_domain_extent
use MOM_io,                   only: stdout, io_infra_end
use mpp_domains_mod,          only: mpp_get_compute_domains
use mpp_domains_mod,          only: mpp_get_ntile_count, mpp_get_pelist, mpp_get_global_domain
use mpp_domains_mod,          only: mpp_get_domain_npes

use MOM_time_manager,         only: set_calendar_type, time_type, set_time, set_date
use MOM_time_manager,         only: GREGORIAN, JULIAN, NOLEAP
use MOM_time_manager,         only: operator( <= ), operator( < ), operator( >= )
use MOM_time_manager,         only: operator( + ),  operator( - ), operator( / )
use MOM_time_manager,         only: operator( * ), operator( /= ), operator( > )
use MOM_domains,              only: MOM_infra_init, MOM_infra_end
use MOM_file_parser,          only: get_param, log_version, param_file_type, close_param_file
use MOM_get_input,            only: get_MOM_input, directories
use MOM_domains,              only: pass_var, pe_here
use MOM_error_handler,        only: MOM_error, FATAL, is_root_pe
use MOM_grid,                 only: ocean_grid_type, get_global_grid_size
use MOM_ocean_model_nuopc,    only: ice_ocean_boundary_type
use MOM_ocean_model_nuopc,    only: ocean_model_restart, ocean_public_type, ocean_state_type
use MOM_ocean_model_nuopc,    only: ocean_model_init_sfc, ocean_model_flux_init
use MOM_ocean_model_nuopc,    only: ocean_model_init, update_ocean_model, ocean_model_end
use MOM_ocean_model_nuopc,    only: get_ocean_grid, get_eps_omesh, query_ocean_state
use MOM_cap_time,             only: AlarmInit
use MOM_cap_methods,          only: mom_import, mom_export, mom_set_geomtype, mod2med_areacor
use MOM_cap_methods,          only: med2mod_areacor, state_diagnose
use MOM_cap_methods,          only: ChkErr
use MOM_ensemble_manager,     only: ensemble_manager_init
use MOM_coms,                 only: sum_across_PEs

! stub routines for CESMCOUPLED
use mom_cap_outputlog,       only: outputlog_init, outputlog_run, outputlog_restart
#ifdef CESMCOUPLED
use shr_log_mod,             only: shr_log_setLogUnit
use nuopc_shr_methods,       only: get_component_instance
#endif
use time_utils_mod,          only: esmf2fms_time

use, intrinsic :: iso_fortran_env, only: output_unit

use ESMF,  only: ESMF_ClockAdvance, ESMF_ClockGet, ESMF_ClockPrint, ESMF_VMget
use ESMF,  only: ESMF_ClockGetAlarm, ESMF_ClockGetNextTime, ESMF_ClockAdvance
use ESMF,  only: ESMF_ClockSet, ESMF_Clock, ESMF_GeomType_Flag, ESMF_LOGMSG_INFO
use ESMF,  only: ESMF_Grid, ESMF_GridCreate, ESMF_GridAddCoord
use ESMF,  only: ESMF_GridGetCoord, ESMF_GridAddItem, ESMF_GridGetItem
use ESMF,  only: ESMF_GridComp, ESMF_GridCompSetEntryPoint, ESMF_GridCompGet
use ESMF,  only: ESMF_LogWrite, ESMF_LogSetError
use ESMF,  only: ESMF_LOGERR_PASSTHRU, ESMF_KIND_R8, ESMF_RC_VAL_WRONG
use ESMF,  only: ESMF_GEOMTYPE_MESH, ESMF_GEOMTYPE_GRID, ESMF_SUCCESS
use ESMF,  only: ESMF_METHOD_INITIALIZE, ESMF_MethodRemove, ESMF_State
use ESMF,  only: ESMF_LOGMSG_INFO, ESMF_RC_ARG_BAD, ESMF_VM, ESMF_Time
use ESMF,  only: ESMF_TimeInterval, ESMF_MAXSTR, ESMF_VMGetCurrent
use ESMF,  only: ESMF_VMGet, ESMF_TimeGet, ESMF_TimeIntervalGet, ESMF_MeshGet
use ESMF,  only: ESMF_MethodExecute, ESMF_Mesh, ESMF_DeLayout, ESMF_Distgrid
use ESMF,  only: ESMF_DistGridConnection, ESMF_StateItem_Flag, ESMF_KIND_I4
use ESMF,  only: ESMF_KIND_I8, ESMF_FAILURE, ESMF_DistGridCreate, ESMF_MeshCreate
use ESMF,  only: ESMF_FILEFORMAT_ESMFMESH, ESMF_DELayoutCreate, ESMF_DistGridConnectionSet
use ESMF,  only: ESMF_DistGridGet, ESMF_STAGGERLOC_CORNER, ESMF_GRIDITEM_MASK
use ESMF,  only: ESMF_TYPEKIND_I4, ESMF_TYPEKIND_R8, ESMF_STAGGERLOC_CENTER
use ESMF,  only: ESMF_GRIDITEM_AREA, ESMF_Field, ESMF_ALARM, ESMF_VMLogMemInfo
use ESMF,  only: ESMF_AlarmIsRinging, ESMF_AlarmRingerOff, ESMF_StateRemove
use ESMF,  only: ESMF_FieldCreate, ESMF_LOGMSG_ERROR, ESMF_LOGMSG_WARNING
use ESMF,  only: ESMF_COORDSYS_SPH_DEG, ESMF_GridCreate, ESMF_INDEX_DELOCAL
use ESMF,  only: ESMF_MESHLOC_ELEMENT, ESMF_RC_VAL_OUTOFRANGE, ESMF_StateGet
use ESMF,  only: ESMF_TimePrint, ESMF_AlarmSet, ESMF_FieldGet, ESMF_Array
use ESMF,  only: ESMF_FieldRegridGetArea
use ESMF,  only: ESMF_ArrayCreate
use ESMF,  only: ESMF_RC_FILE_OPEN, ESMF_RC_FILE_READ, ESMF_RC_FILE_WRITE
use ESMF,  only: ESMF_VMBroadcast, ESMF_VMReduce, ESMF_REDUCE_MAX, ESMF_REDUCE_MIN
use ESMF,  only: ESMF_AlarmCreate, ESMF_ClockGetAlarmList, ESMF_AlarmList_Flag
use ESMF,  only: ESMF_AlarmGet, ESMF_AlarmIsCreated, ESMF_ALARMLIST_ALL, ESMF_AlarmIsEnabled
use ESMF,  only: ESMF_STATEITEM_NOTFOUND, ESMF_FieldWrite
use ESMF,  only: ESMF_END_ABORT, ESMF_Finalize
use ESMF,  only: ESMF_REDUCE_MAX, ESMF_REDUCE_MIN, ESMF_VMAllReduce
use ESMF,  only: operator(==), operator(/=), operator(+), operator(-)

! TODO ESMF_GridCompGetInternalState does not have an explicit Fortran interface.
!! Model does not compile with "use ESMF,  only: ESMF_GridCompGetInternalState"
!! Is this okay?

use NUOPC,       only: NUOPC_CompDerive, NUOPC_CompSetEntryPoint, NUOPC_CompSpecialize
use NUOPC,       only: NUOPC_CompFilterPhaseMap, NUOPC_CompAttributeGet, NUOPC_CompAttributeAdd
use NUOPC,       only: NUOPC_Advertise, NUOPC_SetAttribute, NUOPC_IsUpdated, NUOPC_Write
use NUOPC,       only: NUOPC_IsConnected, NUOPC_Realize, NUOPC_CompAttributeSet
use NUOPC_Model, only: NUOPC_ModelGet
use NUOPC_Model, only: model_routine_SS           => SetServices
use NUOPC_Model, only: model_label_Advance        => label_Advance
use NUOPC_Model, only: model_label_DataInitialize => label_DataInitialize
use NUOPC_Model, only: model_label_SetRunClock    => label_SetRunClock
use NUOPC_Model, only: model_label_Finalize       => label_Finalize
use NUOPC_Model, only: SetVM

use mom_inline_mod, only : mom_inline_init, mom_inline_run
#ifndef CESMCOUPLED
use shr_is_restart_fh_mod, only : init_is_restart_fh, is_restart_fh, is_restart_fh_type
#endif
use mom_cap_profiling, only: cap_profiling_init, cap_profiling

implicit none ; private

public SetServices
public SetVM

!> Internal state type with pointers to three types defined by MOM.
type ocean_internalstate_type
  type(ocean_public_type),       pointer :: ocean_public_type_ptr
  type(ocean_state_type),        pointer :: ocean_state_type_ptr
  type(ice_ocean_boundary_type), pointer :: ice_ocean_boundary_type_ptr
end type

!>  Wrapper-derived type required to associate an internal state instance
!! with the ESMF/NUOPC component
type ocean_internalstate_wrapper
  type(ocean_internalstate_type), pointer :: ptr
end type

!> Contains field information
type fld_list_type
  character(len=64) :: stdname
  character(len=64) :: shortname
  character(len=64) :: transferOffer
  integer :: ungridded_lbound = 0
  integer :: ungridded_ubound = 0
end type fld_list_type

integer,parameter    :: fldsMax = 100
integer              :: fldsToOcn_num = 0
type (fld_list_type) :: fldsToOcn(fldsMax)
integer              :: fldsFrOcn_num = 0
type (fld_list_type) :: fldsFrOcn(fldsMax)

integer              :: dbug = 0
integer              :: import_slice = 1
integer              :: export_slice = 1
character(len=256)   :: tmpstr
logical              :: write_diagnostics = .false.
logical              :: overwrite_timeslice = .false.
logical              :: write_runtimelog = .false.
character(len=32)    :: runtype  !< run type
logical              :: profile_memory = .true.
logical              :: grid_attach_area = .false.
logical              :: use_coldstart = .true.
logical              :: use_mommesh = .true.
logical              :: set_missing_stks_to_zero = .false.
logical              :: restart_eor = .false.
logical              :: use_cdeps_inline = .false.
character(len=128)   :: scalar_field_name = ''
integer              :: scalar_field_count = 0
integer              :: scalar_field_idx_grid_nx = 0
integer              :: scalar_field_idx_grid_ny = 0
character(len=*),parameter :: u_FILE_u = &
     __FILE__

#ifdef CESMCOUPLED
logical :: cesm_coupled = .true.
type(ESMF_GeomType_Flag) :: geomtype = ESMF_GEOMTYPE_MESH
#else
logical :: cesm_coupled = .false.
type(ESMF_GeomType_Flag) :: geomtype
type(is_restart_fh_type) :: restartfh_info     ! For flexible restarts in UFS
#endif
character(len=8)  :: restart_mode = 'alarms'
character(len=16) :: inst_suffix = ''
logical           :: pointer_date = .true. ! append date to rpointer
real(8) :: timere
integer :: localPet = -1


  interface
module subroutine SetServices(gcomp, rc)

  type(ESMF_GridComp)  :: gcomp !< an ESMF_GridComp object
  integer, intent(out) :: rc    !< return code

  ! local variables


end subroutine SetServices
module subroutine InitializeP0(gcomp, importState, exportState, clock, rc)
  type(ESMF_GridComp)   :: gcomp                    !< ESMF_GridComp object
  type(ESMF_State)      :: importState, exportState !< ESMF_State object for
                                                    !! import/export fields
  type(ESMF_Clock)      :: clock                    !< ESMF_Clock object
  integer, intent(out)  :: rc                       !< return code

  ! local variables

end subroutine InitializeP0
module subroutine InitializeAdvertise(gcomp, importState, exportState, clock, rc)
  type(ESMF_GridComp)            :: gcomp                    !< ESMF_GridComp object
  type(ESMF_State)               :: importState, exportState !< ESMF_State object for
                                                             !! import/export fields
  type(ESMF_Clock)               :: clock                    !< ESMF_Clock object
  integer, intent(out)           :: rc                       !< return code

  ! local variables
                                                                 ! (same as restartfile if single restart file)
!--------------------------------

end subroutine InitializeAdvertise
module subroutine InitializeRealize(gcomp, importState, exportState, clock, rc)
  type(ESMF_GridComp)  :: gcomp                    !< ESMF_GridComp object
  type(ESMF_State)     :: importState, exportState !< ESMF_State object for
                                                   !! import/export fields
  type(ESMF_Clock)     :: clock                    !< ESMF_Clock object
  integer, intent(out) :: rc                       !< return code

  ! Local Variables
  !--------------------------------

end subroutine InitializeRealize
module subroutine DataInitialize(gcomp, rc)
  type(ESMF_GridComp)  :: gcomp !< ESMF_GridComp object
  integer, intent(out) :: rc    !< return code

  ! local variables
  !--------------------------------

end subroutine DataInitialize
module subroutine ModelAdvance(gcomp, rc)
  type(ESMF_GridComp)                    :: gcomp !< ESMF_GridComp object
  integer, intent(out)                   :: rc    !< return code

  ! local variables

end subroutine ModelAdvance
module subroutine ModelSetRunClock(gcomp, rc)
  type(ESMF_GridComp)  :: gcomp
  integer, intent(out) :: rc

  ! local variables
  !--------------------------------

end subroutine ModelSetRunClock
module subroutine ocean_model_finalize(gcomp, rc)

  type(ESMF_GridComp)  :: gcomp !< ESMF_GridComp object
  integer, intent(out) :: rc    !< return code

  ! local variables

end subroutine ocean_model_finalize
module subroutine State_SetScalar(value, scalar_id, State, mytask, scalar_name, scalar_count,  rc)
  real(ESMF_KIND_R8),intent(in)     :: value
  integer,           intent(in)     :: scalar_id
  type(ESMF_State),  intent(inout)  :: State
  integer,           intent(in)     :: mytask
  character(len=*),  intent(in)     :: scalar_name
  integer,           intent(in)     :: scalar_count
  integer,           intent(inout)  :: rc           !< return code

  ! local variables
  !--------------------------------------------------------

end subroutine State_SetScalar
module subroutine MOM_RealizeFields(state, nfields, field_defs, tag, ice_ocean_boundary, grid, mesh, rc)
  type(ESMF_State)             , intent(inout)           :: state !< ESMF_State object for
                                                                  !! import/export fields.
  integer                      , intent(in)              :: nfields !< Number of fields.
  type(fld_list_type)          , intent(inout)           :: field_defs(:) !< Structure with field's
                                                                          !! information.
  type(ice_ocean_boundary_type), intent(inout), optional :: ice_ocean_boundary  !< May need to nullify atm_co2
  character(len=*)             , intent(in)              :: tag !< Import or export.
  type(ESMF_Grid)              , intent(in)   , optional :: grid!< ESMF grid.
  type(ESMF_Mesh)              , intent(in)   , optional :: mesh!< ESMF mesh.
  integer                      , intent(inout)           :: rc  !< Return code.

  ! local variables
  !--------------------------------------------------------

end subroutine MOM_RealizeFields
module subroutine fld_list_add(num, fldlist, stdname, transferOffer, shortname, ungridded_lbound, ungridded_ubound)
  integer,                    intent(inout) :: num
  type(fld_list_type),        intent(inout) :: fldlist(:)
  character(len=*),           intent(in)    :: stdname
  character(len=*),           intent(in)    :: transferOffer
  character(len=*), optional, intent(in)    :: shortname
  integer, optional,          intent(in)    :: ungridded_lbound
  integer, optional,          intent(in)    :: ungridded_ubound

  ! local variables

  ! fill in the new entry
end subroutine fld_list_add
module subroutine shr_log_setLogUnit(nunit)
  integer, intent(in) :: nunit
  ! do nothing for this stub - its just here to replace
  ! having cppdefs in the main program
end subroutine shr_log_setLogUnit
  end interface

end module MOM_cap_mod
