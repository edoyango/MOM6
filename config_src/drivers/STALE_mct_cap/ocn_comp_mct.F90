! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> This is the main driver for MOM6 in CIME
module ocn_comp_mct

! mct modules
use ESMF,                only: ESMF_clock, ESMF_time, ESMF_timeInterval
use ESMF,                only: ESMF_ClockGet, ESMF_TimeGet, ESMF_TimeIntervalGet
use seq_cdata_mod,       only: seq_cdata, seq_cdata_setptrs
use seq_flds_mod,        only: seq_flds_x2o_fields, seq_flds_o2x_fields
use mct_mod,             only: mct_gsMap, mct_gsmap_init, mct_gsMap_lsize, &
                               mct_gsmap_orderedpoints
use mct_mod,             only: mct_aVect, mct_aVect_init, mct_aVect_zero, &
                               mct_aVect_nRattr
use mct_mod,             only: mct_gGrid, mct_gGrid_init, mct_gGrid_importRAttr, &
                               mct_gGrid_importIAttr
use seq_infodata_mod,    only: seq_infodata_type, seq_infodata_GetData, &
                               seq_infodata_start_type_start, seq_infodata_start_type_cont,   &
                               seq_infodata_start_type_brnch, seq_infodata_PutData
use seq_comm_mct,        only: seq_comm_name, seq_comm_inst, seq_comm_suffix
use seq_timemgr_mod,     only: seq_timemgr_EClockGetData, seq_timemgr_RestartAlarmIsOn
use perf_mod,            only: t_startf, t_stopf
use shr_file_mod,        only: shr_file_getUnit, shr_file_freeUnit, shr_file_setIO, &
                               shr_file_getLogUnit, shr_file_getLogLevel, &
                               shr_file_setLogUnit, shr_file_setLogLevel

! MOM6 modules
use MOM,                  only: extract_surface_state
use MOM_variables,        only: surface
use MOM_domains,          only: MOM_infra_init
use MOM_restart,          only: save_restart
use MOM_ice_shelf,        only: ice_shelf_save_restart, adjust_ice_sheet_frazil
use MOM_domains,          only: num_pes, root_pe, pe_here
use MOM_grid,             only: ocean_grid_type, get_global_grid_size
use MOM_error_handler,    only: MOM_error, FATAL, is_root_pe, WARNING
use MOM_time_manager,     only: time_type, set_date, set_time, set_calendar_type, NOLEAP
use MOM_time_manager,     only: operator(+), operator(-), operator(*), operator(/)
use MOM_time_manager,     only: operator(==), operator(/=), operator(>), get_time
use MOM_file_parser,      only: get_param, log_version, param_file_type, close_param_file
use MOM_get_input,        only: Get_MOM_Input, directories
use MOM_EOS,              only: gsw_sp_from_sr, gsw_pt_from_ct
use MOM_constants,        only: CELSIUS_KELVIN_OFFSET
use MOM_domains,          only: AGRID, BGRID_NE, CGRID_NE, pass_vector
use mpp_domains_mod,      only: mpp_get_compute_domain
use MOM_io,               only: stdout

! Previously inlined - now in separate modules
use MOM_ocean_model_mct,     only: ocean_public_type, ocean_state_type
use MOM_ocean_model_mct,     only: ocean_model_init , update_ocean_model, ocean_model_end
use MOM_ocean_model_mct,     only: convert_state_to_ocean_type
use MOM_surface_forcing_mct, only: surface_forcing_CS, forcing_save_restart, ice_ocean_boundary_type
use ocn_cap_methods,      only: ocn_import, ocn_export

! MCT indices structure and import and export routines that access mom data
use ocn_cpl_indices,   only : cpl_indices_type, cpl_indices_init

! GFDL coupler modules
use MOM_coupler_types,   only : coupler_type_spawn
use MOM_coupler_types,   only : coupler_type_initialized, coupler_type_copy_data

! By default make data private
implicit none ; private

#include <MOM_memory.h>

! Public member functions
public :: ocn_init_mct
public :: ocn_run_mct
public :: ocn_final_mct

! Private member functions
private :: ocn_SetGSMap_mct
private :: ocn_domain_mct
private :: get_runtype
private :: ocean_model_init_sfc

! Flag for debugging
logical, parameter :: debug=.true.

!> Control structure for this module
type MCT_MOM_Data
  type(ocean_state_type),  pointer :: ocn_state => NULL()  !< The private state of ocean
  type(ocean_public_type), pointer :: ocn_public => NULL() !< The public state of ocean
  type(ocean_grid_type),   pointer :: grid => NULL()       !< The grid structure
  type(seq_infodata_type), pointer :: infodata             !< The input info type
  type(cpl_indices_type)           :: ind                  !< Variable IDs
  logical                          :: sw_decomp            !< Controls whether shortwave is decomposed into 4 components
  real                             :: c1, c2, c3, c4       !< Coeffs. used in the shortwave decomposition  i/o
  character(len=384)               :: pointer_filename     !< Name of the ascii file that contains the path
                                                           !! and filename of the latest restart file.
end type MCT_MOM_Data

type(MCT_MOM_Data)            :: glb !< global structure
type(ice_ocean_boundary_type) :: ice_ocean_boundary

!=======================================================================

  interface
module subroutine ocn_init_mct( EClock, cdata_o, x2o_o, o2x_o, NLFilename )
  type(ESMF_Clock),             intent(inout) :: EClock      !< Time and time step ? \todo Why must this
                                                             !! be intent(inout)?
  type(seq_cdata)             , intent(inout) :: cdata_o     !< Input parameters
  type(mct_aVect)             , intent(inout) :: x2o_o       !< Fluxes from coupler to ocean, computed by ocean
  type(mct_aVect)             , intent(inout) :: o2x_o       !< Fluxes from ocean to coupler, computed by ocean
  character(len=*), optional  , intent(in)    :: NLFilename  !< Namelist filename

  !  local variable
                                                  !! (same as restartfile if a single restart file is to be read in)
                                                        !! in the x and y dir., respectively.
  ! runtime params

  ! mct variables (these are local for now)

  ! time management



  ! TODO: Change the following vars with the corresponding MOM6 vars
  !logical :: lsend_precip_fact      !< If T,send precip_fact to cpl for use in fw balance
                                    !! (partially-coupled option)

  ! set the cdata pointers:
end subroutine ocn_init_mct
module subroutine ocn_run_mct( EClock, cdata_o, x2o_o, o2x_o)
  type(ESMF_Clock), intent(inout) :: EClock  !< Time and time step ? \todo Why must this be intent(inout)?
  type(seq_cdata),  intent(inout) :: cdata_o !< Input parameters
  type(mct_aVect),  intent(inout) :: x2o_o   !< Fluxes from coupler to ocean, computed by ocean
  type(mct_aVect),  intent(inout) :: o2x_o   !< Fluxes from ocean to coupler, computed by ocean
  ! Local variables

  ! reset shr logging to ocn log file:
end subroutine ocn_run_mct
module subroutine ocn_final_mct( EClock, cdata_o, x2o_o, o2x_o)
    type(ESMF_Clock)            , intent(inout) :: EClock
    type(seq_cdata)             , intent(inout) :: cdata_o
    type(mct_aVect)             , intent(inout) :: x2o_o  !< Fluxes from coupler to ocean, computed by ocean
    type(mct_aVect)             , intent(inout) :: o2x_o  !< Fluxes from ocean to coupler, computed by ocean

end subroutine ocn_final_mct
module subroutine ocn_SetGSMap_mct(mpicom_ocn, MOM_MCT_ID, gsMap_ocn, gsMap3d_ocn)
  integer,         intent(in)    :: mpicom_ocn  !< MPI communicator
  integer,         intent(in)    :: MOM_MCT_ID  !< MCT component ID
  type(mct_gsMap), intent(inout) :: gsMap_ocn   !< MCT global segment map for 2d data
  type(mct_gsMap), intent(inout) :: gsMap3d_ocn !< MCT global segment map for 3d data

  ! Local variables

end subroutine ocn_SetGSMap_mct
module subroutine ocn_domain_mct( lsize, gsMap_ocn, dom_ocn)
  integer        , intent(in)    :: lsize      !< Size of attr. vector
  type(mct_gsMap), intent(in)    :: gsMap_ocn  !< MCT global segment map for 2d data
  type(mct_ggrid), intent(inout) :: dom_ocn    !< WHAT IS THIS?

  ! Local Variables

end subroutine ocn_domain_mct
character(32) module function get_runtype()

end function get_runtype
module subroutine ocean_model_init_sfc(OS, Ocean_sfc)
  type(ocean_state_type),  pointer       :: OS
  type(ocean_public_type), intent(inout) :: Ocean_sfc


end subroutine ocean_model_init_sfc
module subroutine IOB_allocate(IOB, isc, iec, jsc, jec)
  type(ice_ocean_boundary_type), intent(inout)    :: IOB    !< An ice-ocean boundary type with fluxes to drive
  integer, intent(in) :: isc, iec, jsc, jec                 !< The ocean's local grid size

end subroutine IOB_allocate
  end interface

end module ocn_comp_mct
