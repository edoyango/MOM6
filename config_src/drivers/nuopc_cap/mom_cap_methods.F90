! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> Contains import/export methods for CMEPS.
module MOM_cap_methods

use ESMF,                      only: ESMF_Clock, ESMF_ClockGet, ESMF_time, ESMF_TimeGet
use ESMF,                      only: ESMF_TimeInterval, ESMF_TimeIntervalGet
use ESMF,                      only: ESMF_State, ESMF_StateGet
use ESMF,                      only: ESMF_Field, ESMF_FieldGet, ESMF_FieldCreate
use ESMF,                      only: ESMF_GridComp, ESMF_Mesh, ESMF_MeshGet, ESMF_Grid, ESMF_GridCreate
use ESMF,                      only: ESMF_DistGrid, ESMF_DistGridCreate
use ESMF,                      only: ESMF_KIND_R8, ESMF_SUCCESS, ESMF_LogFoundError
use ESMF,                      only: ESMF_LOGERR_PASSTHRU, ESMF_LOGMSG_INFO, ESMF_LOGWRITE
use ESMF,                      only: ESMF_LogSetError, ESMF_RC_MEM_ALLOCATE
use ESMF,                      only: ESMF_StateItem_Flag, ESMF_STATEITEM_NOTFOUND
use ESMF,                      only: ESMF_GEOMTYPE_FLAG, ESMF_GEOMTYPE_GRID, ESMF_GEOMTYPE_MESH
use ESMF,                      only: ESMF_RC_VAL_OUTOFRANGE, ESMF_INDEX_DELOCAL, ESMF_MESHLOC_ELEMENT
use ESMF,                      only: ESMF_TYPEKIND_R8, ESMF_FIELDSTATUS_COMPLETE
use ESMF,                      only: ESMF_FieldStatus_Flag, ESMF_LOGMSG_ERROR, ESMF_FAILURE, ESMF_MAXSTR
use ESMF,                      only: operator(/=), operator(==)
use MOM_ocean_model_nuopc,     only: ocean_public_type, ocean_state_type
use MOM_surface_forcing_nuopc, only: ice_ocean_boundary_type
use MOM_grid,                  only: ocean_grid_type
use MOM_domains,               only: pass_var
use mpp_domains_mod,           only: mpp_get_compute_domain

! By default make data private
implicit none ; private

! Public member functions
public :: mom_set_geomtype
public :: mom_import
public :: mom_export
public :: state_diagnose
public :: ChkErr

interface State_getImport
   module procedure State_getImport_2d
   module procedure State_getImport_3d ! third dimension being an ungridded dimension
end interface

private :: State_setExport

!> Get field pointer
interface State_GetFldPtr
  module procedure State_GetFldPtr_1d
  module procedure State_GetFldPtr_2d
end interface

integer                  :: import_cnt = 0!< used to skip using the import state
                                          !! at the first count for cesm
type(ESMF_GeomType_Flag) :: geomtype      !< SMF type describing type of
                                          !! geometry (mesh or grid)

! area correction factors for fluxes send and received from mediator
! these actors are ONLY valid for meshes that are read in - so do not need them for
! grids that are calculated internally

real(ESMF_KIND_R8), public, allocatable :: mod2med_areacor(:) ! ratios of model areas to input mesh areas
real(ESMF_KIND_R8), public, allocatable :: med2mod_areacor(:) ! ratios of input mesh areas to model areas
character(len=*),parameter :: u_FILE_u =  __FILE__


  interface
module subroutine mom_set_geomtype(geomtype_in)
  type(ESMF_GeomType_Flag), intent(in)    :: geomtype_in !< ESMF type describing type of
                                                         !! geometry (mesh or grid)

end subroutine mom_set_geomtype
module subroutine mom_import(ocean_public, ocean_grid, importState, ice_ocean_boundary, &
                      set_missing_stks_to_zero, rc)
  type(ocean_public_type)       , intent(in)    :: ocean_public             !< Ocean surface state
  type(ocean_grid_type)         , intent(in)    :: ocean_grid               !< Ocean model grid
  logical                       , intent(in)    :: set_missing_stks_to_zero !< If true, set
                                                                            !! missing stokes drift to zero
  type(ESMF_State)              , intent(inout) :: importState              !< incoming data from mediator
  type(ice_ocean_boundary_type) , intent(inout) :: ice_ocean_boundary       !< Ocean boundary forcing
  integer                       , intent(inout) :: rc                       !< Return code

  ! Local Variables

end subroutine mom_import
module subroutine mom_export(ocean_public, ocean_grid, ocean_state, exportState, clock, rc)
  type(ocean_public_type) , intent(in)    :: ocean_public !< Ocean surface state
  type(ocean_grid_type)   , intent(in)    :: ocean_grid   !< Ocean model grid
  type(ocean_state_type)  , pointer       :: ocean_state  !< Ocean state
  type(ESMF_State)        , intent(inout) :: exportState  !< outgoing data
  type(ESMF_Clock)        , intent(in)    :: clock        !< ESMF clock
  integer                 , intent(inout) :: rc           !< Return code

  ! Local variables

end subroutine mom_export
module subroutine State_GetFldPtr_1d(State, fldname, fldptr, rc)
  type(ESMF_State)            , intent(in)  :: State    !< ESMF state
  character(len=*)            , intent(in)  :: fldname  !< Field name
  real(ESMF_KIND_R8), pointer , intent(in)  :: fldptr(:)!< Pointer to the 1D field
  integer, optional           , intent(out) :: rc       !< Return code

  ! local variables

end subroutine State_GetFldPtr_1d
module subroutine State_GetFldPtr_2d(State, fldname, fldptr, rc)
  type(ESMF_State)            , intent(in)  :: State      !< ESMF state
  character(len=*)            , intent(in)  :: fldname    !< Field name
  real(ESMF_KIND_R8), pointer , intent(in)  :: fldptr(:,:)!< Pointer to the 2D field
  integer, optional           , intent(out) :: rc         !< Return code

  ! local variables

end subroutine State_GetFldPtr_2d
module subroutine State_GetImport_2d(state, fldname, isc, iec, jsc, jec, output, do_sum, areacor, esmf_ind, rc)
  type(ESMF_State)    , intent(in)    :: state   !< ESMF state
  character(len=*)    , intent(in)    :: fldname !< Field name
  integer             , intent(in)    :: isc     !< The start i-index of cell centers within
                                                 !! the computational domain
  integer             , intent(in)    :: iec     !< The end i-index of cell centers within the
                                                 !! computational domain
  integer             , intent(in)    :: jsc     !< The start j-index of cell centers within
                                                 !! the computational domain
  integer             , intent(in)    :: jec     !< The end j-index of cell centers within
                                                 !! the computational domain
  real (ESMF_KIND_R8) , intent(inout) :: output(isc:iec,jsc:jec)!< Output 2D array
  logical, optional   , intent(in)    :: do_sum  !< If true, sums the data
  real (ESMF_KIND_R8), optional,  intent(in) :: areacor(:) !< flux area correction factors
                                                           !! applicable to meshes
  integer,             optional, intent(in) :: esmf_ind
  integer             , intent(out)   :: rc      !< Return code

  ! local variables
  ! ----------------------------------------------

end subroutine State_GetImport_2d
module subroutine State_GetImport_3d(state, fldname, isc, iec, jsc, jec, lbd, ubd, output, do_sum, areacor, rc)
  type(ESMF_State)    , intent(in)    :: state   !< ESMF state
  character(len=*)    , intent(in)    :: fldname !< Field name
  integer             , intent(in)    :: isc     !< The start i-index of cell centers within
                                                 !! the computational domain
  integer             , intent(in)    :: iec     !< The end i-index of cell centers within the
                                                 !! computational domain
  integer             , intent(in)    :: jsc     !< The start j-index of cell centers within
                                                 !! the computational domain
  integer             , intent(in)    :: jec     !< The end j-index of cell centers within
                                                 !! the computational domain
  integer             , intent(in)    :: lbd     !< lower bound of ungridded dimension
  integer             , intent(in)    :: ubd     !< upper bound of ungridded dimension
  real (ESMF_KIND_R8) , intent(inout) :: output(isc:iec,jsc:jec,lbd:ubd)!< Output 3D array
  logical, optional   , intent(in)    :: do_sum  !< If true, sums the data
  real (ESMF_KIND_R8), optional,  intent(in) :: areacor(:) !< flux area correction factors
                                                           !! applicable to meshes
  integer             , intent(out)   :: rc      !< Return code

  ! local variables
  ! ----------------------------------------------

end subroutine State_GetImport_3d
module subroutine State_SetExport(state, fldname, isc, iec, jsc, jec, input, ocean_grid, areacor, rc)
  type(ESMF_State)      , intent(inout) :: state   !< ESMF state
  character(len=*)      , intent(in)    :: fldname !< Field name
  integer             , intent(in)      :: isc     !< The start i-index of cell centers within
                                                   !! the computational domain
  integer             , intent(in)      :: iec     !< The end i-index of cell centers within the
                                                   !! computational domain
  integer             , intent(in)      :: jsc     !< The start j-index of cell centers within
                                                   !! the computational domain
  integer             , intent(in)      :: jec     !< The end j-index of cell centers within
                                                   !! the computational domain
  real (ESMF_KIND_R8)   , intent(in)    :: input(isc:iec,jsc:jec)!< Input 2D array
  type(ocean_grid_type) , intent(in)    :: ocean_grid !< Ocean horizontal grid
  real (ESMF_KIND_R8), optional,  intent(in) :: areacor(:) !< flux area correction factors
                                                           !! applicable to meshes
  integer               , intent(out)   :: rc         !< Return code

  ! local variables
  ! ----------------------------------------------

end subroutine State_SetExport
module subroutine state_diagnose(State, string, rc)

  ! ----------------------------------------------
  ! Diagnose status of State
  ! ----------------------------------------------

  type(ESMF_State), intent(in)    :: state  !< An ESMF State
  character(len=*), intent(in)    :: string !< A string indicating whether the State is an
                                            !! import or export State
  integer         , intent(out)   :: rc     !< Return code

  ! local variables
  ! ----------------------------------------------

end subroutine state_diagnose
module subroutine field_getfldptr(field, fldptr1, fldptr2, rank, abort, rc)

  ! input/output variables
  type(ESMF_Field)  , intent(in)                        :: field        !< An ESMF field
  real(ESMF_KIND_R8), pointer , intent(inout), optional :: fldptr1(:)   !< A pointer to a rank 1 ESMF field
  real(ESMF_KIND_R8), pointer , intent(inout), optional :: fldptr2(:,:) !< A pointer to a rank 2 ESMF field
  integer           , intent(out)            , optional :: rank         !< Field rank
  logical           , intent(in)             , optional :: abort        !< Abort code
  integer           , intent(out)            , optional :: rc           !< Return code

  ! local variables
  ! ----------------------------------------------

end subroutine field_getfldptr
logical module function ChkErr(rc, line, file)
  integer, intent(in) :: rc            !< return code to check
  integer, intent(in) :: line          !< Integer source line number
  character(len=*), intent(in) :: file !< User-provided source file name
end function ChkErr
  end interface

end module MOM_cap_methods
