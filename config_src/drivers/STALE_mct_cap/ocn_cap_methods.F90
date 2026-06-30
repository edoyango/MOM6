! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

module ocn_cap_methods

  use ESMF,                    only: ESMF_clock, ESMF_time, ESMF_ClockGet, ESMF_TimeGet
  use MOM_ocean_model_mct,     only: ocean_public_type, ocean_state_type
  use MOM_surface_forcing_mct, only: ice_ocean_boundary_type
  use MOM_grid,                only: ocean_grid_type
  use MOM_domains,             only: pass_var
  use MOM_error_handler,       only: is_root_pe
  use mpp_domains_mod,         only: mpp_get_compute_domain
  use ocn_cpl_indices,         only: cpl_indices_type

  implicit none
  private

  public :: ocn_import
  public :: ocn_export

  logical, parameter :: debug=.false.

!=======================================================================

  interface
module subroutine ocn_import(x2o, ind, grid, ice_ocean_boundary, ocean_public, logunit, Eclock, c1, c2, c3, c4)
  real(kind=8)                  , intent(in)    :: x2o(:,:)           !< incoming data
  type(cpl_indices_type)        , intent(in)    :: ind                !< Structure with MCT attribute vects and indices
  type(ocean_grid_type)         , intent(in)    :: grid               !< Ocean model grid
  type(ice_ocean_boundary_type) , intent(inout) :: ice_ocean_boundary !< Ocean boundary forcing
  type(ocean_public_type)       , intent(in)    :: ocean_public       !< Ocean surface state
  integer                       , intent(in)    :: logunit            !< Unit for stdout output
  type(ESMF_Clock)              , intent(in)    :: EClock             !< Time and time step ? \todo Why must this
  real(kind=8), optional        , intent(in)    :: c1, c2, c3, c4     !< Coeffs. used in the shortwave decomposition

  ! Local variables
  !-----------------------------------------------------------------------

end subroutine ocn_import
module subroutine ocn_export(ind, ocn_public, grid, o2x, dt_int, ncouple_per_day)
  type(cpl_indices_type),  intent(inout) :: ind        !< Structure with coupler indices and vectors
  type(ocean_public_type), intent(in)    :: ocn_public !< Ocean surface state
  type(ocean_grid_type),   intent(in)    :: grid       !< Ocean model grid
  real(kind=8),            intent(inout) :: o2x(:,:)   !< MCT outgoing bugger
  real(kind=8), intent(in)               :: dt_int     !< Amount of time over which to advance the
                                                       !! ocean (ocean_coupling_time_step), in sec
  integer, intent(in)                    :: ncouple_per_day !< Number of ocean coupling calls per day

  ! Local variables

  !-----------------------------------------------------------------------

  ! Use Adcroft's rule of reciprocals; it does the right thing here.
end subroutine ocn_export
  end interface

end module ocn_cap_methods
