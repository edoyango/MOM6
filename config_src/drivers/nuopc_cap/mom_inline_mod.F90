!> This module contains a set of subroutines that enables inline CDEPS capability

module mom_inline_mod

use NUOPC            , only: NUOPC_CompAttributeGet
use ESMF             , only: ESMF_GridComp, ESMF_Mesh
use ESMF             , only: ESMF_Clock, ESMF_Time, ESMF_TimeGet, ESMF_ClockGet
use ESMF             , only: ESMF_KIND_R8, ESMF_SUCCESS, ESMF_LogFoundError
use ESMF             , only: ESMF_LOGERR_PASSTHRU, ESMF_LOGMSG_INFO, ESMF_LOGWRITE
use ESMF             , only: ESMF_END_ABORT, ESMF_Finalize, ESMF_MAXSTR
use dshr_mod         , only: dshr_pio_init
use dshr_strdata_mod , only: shr_strdata_type, shr_strdata_print
use dshr_strdata_mod , only: shr_strdata_init_from_inline
use dshr_strdata_mod , only: shr_strdata_advance
use dshr_methods_mod , only: dshr_fldbun_getfldptr, dshr_fldbun_Field_diagnose
use dshr_stream_mod  , only: shr_stream_init_from_esmfconfig
use MOM_cap_methods  , only: ChkErr

implicit none
private

public mom_inline_init
public mom_inline_run

type(shr_strdata_type), allocatable :: sdat(:)

integer                    :: logunit      ! the logunit on the root task
character(len=ESMF_MAXSTR) :: stream_name  ! generic identifier

character(len=*), parameter :: u_FILE_u =  __FILE__

  interface
module subroutine mom_inline_init(gcomp, model_clock, model_mesh, mytask, rc)
  type(ESMF_GridComp)    , intent(in)  :: gcomp        !< ESMF_GridComp object
  type(ESMF_Clock)       , intent(in)  :: model_clock  !< ESMF_Clock object
  type(ESMF_Mesh)        , intent(in)  :: model_mesh   !< ESMF mesh
  integer                , intent(in)  :: mytask       !< the current task
  integer                , intent(out) :: rc           !< Return code

  ! local variables


  !----------------------------------------------------------------------

end subroutine mom_inline_init
module subroutine mom_inline_run(clock, ocean_public, ocean_grid, ice_ocean_boundary, dbug, rc)
  use MOM_ocean_model_nuopc,     only: ocean_public_type
  use MOM_surface_forcing_nuopc, only: ice_ocean_boundary_type
  use MOM_grid,                  only: ocean_grid_type
  use mpp_domains_mod,           only: mpp_get_compute_domain

  type(ESMF_Clock) ,              intent(in)    :: clock              !< ESMF_Clock object
  type(ocean_public_type)       , intent(in)    :: ocean_public       !< Ocean surface state
  type(ocean_grid_type)         , intent(in)    :: ocean_grid         !< Ocean model grid
  type(ice_ocean_boundary_type) , intent(inout) :: ice_ocean_boundary !< Ocean boundary forcing
  integer ,                       intent(in)    :: dbug               !< Integer debug flag
  integer ,                       intent(out)   :: rc                 !< Return code

  ! local variables
  !-----------------------------------------------------------------------

end subroutine mom_inline_run
  end interface

end module mom_inline_mod
