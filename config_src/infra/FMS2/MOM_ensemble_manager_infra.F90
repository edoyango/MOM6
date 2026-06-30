! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> A simple (very thin) wrapper for managing ensemble member layout information
module MOM_ensemble_manager_infra

use ensemble_manager_mod, only : FMS_ensemble_manager_init => ensemble_manager_init
use ensemble_manager_mod, only : FMS_ensemble_pelist_setup => ensemble_pelist_setup
use ensemble_manager_mod, only : FMS_get_ensemble_id => get_ensemble_id
use ensemble_manager_mod, only : FMS_get_ensemble_size => get_ensemble_size
use ensemble_manager_mod, only : FMS_get_ensemble_pelist => get_ensemble_pelist
use ensemble_manager_mod, only : FMS_get_ensemble_filter_pelist => get_ensemble_filter_pelist
use fms2_io_mod, only : fms2_io_set_filename_appendix=>set_filename_appendix

implicit none ; private

public :: ensemble_manager_init, ensemble_pelist_setup
public :: get_ensemble_id, get_ensemble_size
public :: get_ensemble_pelist, get_ensemble_filter_pelist


  interface
module subroutine ensemble_manager_init(ensemble_suffix)
  character(len=*), optional, intent(in) :: ensemble_suffix !> Ensemble suffix provided by the cap. This may be
                                                            !! provided to bypass FMS ensemble manager.

end subroutine ensemble_manager_init
module subroutine ensemble_pelist_setup(concurrent, atmos_npes, ocean_npes, land_npes, ice_npes, &
                                   Atm_pelist, Ocean_pelist, Land_pelist, Ice_pelist)
  logical,               intent(in)    :: concurrent !< A logical flag, if True, then ocean fast
                                                     !! PEs are run concurrently with
                                                     !! slow PEs within the coupler.
  integer,               intent(in)    :: atmos_npes !< The number of atmospheric (fast) PEs
  integer,               intent(in)    :: ocean_npes !< The number of ocean (slow) PEs
  integer,               intent(in)    :: land_npes  !< The number of land PEs (fast)
  integer,               intent(in)    :: ice_npes   !< The number of ice (fast) PEs
  integer, dimension(:), intent(inout) :: Atm_pelist !< A list of Atm PEs
  integer, dimension(:), intent(inout) :: Ocean_pelist !< A list of Ocean PEs
  integer, dimension(:), intent(inout) :: Land_pelist !< A list of Land PEs
  integer, dimension(:), intent(inout) :: Ice_pelist !< A list of Ice PEs


end subroutine ensemble_pelist_setup
module function get_ensemble_id()
  integer :: get_ensemble_id

end function get_ensemble_id
module function get_ensemble_size()
  integer, dimension(6) :: get_ensemble_size

end function get_ensemble_size
module subroutine get_ensemble_pelist(pelist, name)
  integer,                    intent(inout) :: pelist(:,:) !< A processor list for all ensemble members
  character(len=*), optional, intent(in)    :: name !< An optional component name (atmos, ocean, land, ice)

end subroutine get_ensemble_pelist
module subroutine get_ensemble_filter_pelist(pelist, name)
  integer,          intent(inout) :: pelist(:) !< A processor list for the ensemble filter
  character(len=*), intent(in)    :: name      !< The component name (atmos, ocean, land, ice)

end subroutine get_ensemble_filter_pelist
  end interface

end module MOM_ensemble_manager_infra
