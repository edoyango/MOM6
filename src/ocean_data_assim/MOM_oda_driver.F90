! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> Interfaces for MOM6 ensembles and data assimilation.
module MOM_oda_driver_mod

! This file is part of MOM6. see LICENSE.md for the license.

! MOM infrastructure
use MOM_coms, only : PE_here, num_PEs
use MOM_coms, only : set_PElist, set_rootPE, Get_PElist, broadcast
use MOM_domains, only : domain2d, global_field, get_domain_extent
use MOM_domains, only : pass_var, redistribute_array, broadcast_domain
use MOM_diag_mediator, only : register_diag_field, diag_axis_init, post_data
use MOM_diag_mediator, only : enable_averaging, disable_averaging
use MOM_diag_mediator, only : diag_update_remap_grids
use MOM_ensemble_manager, only : get_ensemble_id, get_ensemble_size
use MOM_ensemble_manager, only : get_ensemble_pelist, get_ensemble_filter_pelist
use MOM_error_handler, only : stdout, stdlog, MOM_error
use MOM_io, only : SINGLE_FILE
use MOM_interp_infra, only : init_extern_field
use MOM_interp_infra, only : time_interp_extern
use MOM_interpolate, only : external_field
use MOM_interpolate, only : get_external_field_info
use MOM_remapping,    only : remappingSchemesDoc
use MOM_time_manager, only : time_type, real_to_time, get_date
use MOM_time_manager, only : operator(+), operator(>=), operator(/=)
use MOM_time_manager, only : operator(==), operator(<)
use MOM_cpu_clock, only : cpu_clock_begin, cpu_clock_end, cpu_clock_id
use MOM_horizontal_regridding, only : horiz_interp_and_extrap_tracer
! ODA Modules
use ocean_da_types_mod, only : grid_type, ocean_profile_type, ocean_control_struct
use ocean_da_core_mod, only : ocean_da_core_init, get_profiles
!This preprocessing directive enables the SPEAR online ensemble data assimilation
!configuration. Existing community based APIs for data assimilation are currently
!called offline for forecast applications using information read from a MOM6 state file.
!The SPEAR configuration (https://doi.org/10.1029/2020MS002149) calculated increments
!efficiently online. A community-based set of APIs should be implemented in place
!of the CPP directive when this is available.
#ifdef ENABLE_ECDA
use eakf_oda_mod, only : ensemble_filter
#endif
use kdtree, only : kd_root !# A kd-tree object using JEDI APIs
! MOM Modules
use MOM_io, only : slasher, MOM_read_data
use MOM_diag_mediator, only : diag_ctrl, set_axes_info
use MOM_error_handler, only : FATAL, WARNING, MOM_error, MOM_mesg, is_root_pe
use MOM_get_input, only : get_MOM_input, directories
use MOM_grid, only : ocean_grid_type, MOM_grid_init
use MOM_grid_initialize, only : set_grid_metrics
use MOM_hor_index, only : hor_index_type, hor_index_init
use MOM_dyn_horgrid, only : dyn_horgrid_type, create_dyn_horgrid, destroy_dyn_horgrid
use MOM_transcribe_grid, only : copy_dyngrid_to_MOM_grid, copy_MOM_grid_to_dyngrid
use MOM_fixed_initialization, only : MOM_initialize_topography
use MOM_coord_initialization, only : MOM_initialize_coord
use MOM_file_parser, only : read_param, get_param, param_file_type
use MOM_string_functions, only : lowercase
use MOM_ALE, only : ALE_CS, ALE_initThicknessToCoord, ALE_init, ALE_updateVerticalGridType
use MOM_domains, only : MOM_domains_init, MOM_domain_type, clone_MOM_domain
use MOM_remapping, only : remapping_CS, initialize_remapping, remapping_core_h
use MOM_regridding, only : regridding_CS, initialize_regridding
use MOM_regridding, only : regridding_main, set_regrid_params, set_h_neglect
use MOM_unit_scaling, only : unit_scale_type, unit_scaling_init
use MOM_variables, only : thermo_var_ptrs
use MOM_verticalGrid, only : verticalGrid_type, verticalGridInit

implicit none ; private

public :: init_oda, oda_end, set_prior_tracer, get_posterior_tracer
public :: set_analysis_time, oda, apply_oda_tracer_increments

!>@{ CPU time clock ID
integer :: id_clock_oda_init
integer :: id_clock_bias_adjustment
integer :: id_clock_apply_increments
!>@}

#include <MOM_memory.h>

!> A structure with a pointer to a domain2d, to allow for the creation of arrays of pointers.
type :: ptr_mpp_domain
  type(domain2d), pointer :: mpp_domain => NULL() !< pointer to a domain2d
end type ptr_mpp_domain

!> A structure containing integer handles for bias adjustment of tracers
type :: INC_CS
  integer :: fldno = 0 !< The number of tracers
  type(external_field) :: T  !< The handle for the temperature file
  type(external_field) :: S  !< The handle for the salinity file
end type INC_CS

!> Control structure that contains a transpose of the ocean state across ensemble members.
type, public :: ODA_CS ; private
  type(ocean_control_struct), pointer :: Ocean_prior=> NULL() !< ensemble ocean prior states in DA space
  type(ocean_control_struct), pointer :: Ocean_posterior=> NULL() !< ensemble ocean posterior states
                                                                  !! or increments to prior in DA space
  type(ocean_control_struct), pointer :: Ocean_increment=> NULL() !< A separate structure for
                                                                  !! increment diagnostics
  integer :: nk !< number of vertical layers used for DA
  type(ocean_grid_type), pointer :: Grid => NULL() !< MOM6 grid type and decomposition for the DA
  type(ocean_grid_type), pointer :: G => NULL() !< MOM6 grid type and decomposition for the model
  type(MOM_domain_type), pointer, dimension(:) :: domains => NULL() !< Pointer to mpp_domain objects
                                                                       !! for ensemble members
  type(verticalGrid_type), pointer :: GV => NULL() !< vertical grid for DA
  type(unit_scale_type), pointer :: &
    US => NULL()    !< structure containing various unit conversion factors for DA

  type(domain2d), pointer :: mpp_domain => NULL() !< Pointer to a mpp domain object for DA
  type(grid_type), pointer :: oda_grid !< local tracer grid
  real, pointer, dimension(:,:,:) :: h => NULL() !<layer thicknesses [H ~> m or kg m-2] for DA
  real, pointer, dimension(:,:,:) :: T_tend => NULL() !<layer temperature tendency from DA [C T-1 ~> degC s-1]
  real, pointer, dimension(:,:,:) :: S_tend => NULL() !<layer salinity tendency from DA [S T-1 ~> ppt s-1]
  real, pointer, dimension(:,:,:) :: T_bc_tend => NULL() !< The layer temperature tendency due
                                                         !! to bias adjustment [C T-1 ~> degC s-1]
  real, pointer, dimension(:,:,:) :: S_bc_tend => NULL() !< The layer salinity tendency due
                                                         !! to bias adjustment [S T-1 ~> ppt s-1]
  integer :: ni          !< global i-direction grid size
  integer :: nj          !< global j-direction grid size
  logical :: reentrant_x !< grid is reentrant in the x direction
  logical :: reentrant_y !< grid is reentrant in the y direction
  logical :: tripolar_N !< grid is folded at its north edge
  logical :: symmetric !< Values at C-grid locations are symmetric
  logical :: use_basin_mask !< If true, use a basin file to delineate weakly coupled ocean basins
  logical :: do_bias_adjustment !< If true, use spatio-temporally varying climatological tendency
                                !! adjustment for Temperature and Salinity
  real :: bias_adjustment_multiplier !< A scaling for the bias adjustment [nondim]
  integer :: assim_method !< Method: NO_ASSIM,EAKF_ASSIM or OI_ASSIM
  integer :: ensemble_size !< Size of the ensemble
  integer :: ensemble_id = 0 !< id of the current ensemble member
  integer, pointer, dimension(:,:) :: ensemble_pelist !< PE list for ensemble members
  integer, pointer, dimension(:) :: filter_pelist !< PE list for ensemble members
  real :: assim_interval !< analysis interval [T ~> s]
  ! Profiles local to the analysis domain
  type(ocean_profile_type), pointer :: Profiles => NULL() !< pointer to linked list of all available profiles
  type(ocean_profile_type), pointer :: CProfiles => NULL()!< pointer to linked list of current profiles
  type(kd_root), pointer :: kdroot => NULL() !< A structure for storing nearest neighbors
  type(ALE_CS), pointer :: ALE_CS=>NULL() !< ALE control structure for DA
  logical :: use_ALE_algorithm !< true is using ALE remapping
  type(regridding_CS) :: regridCS !< ALE control structure for regridding
  type(remapping_CS) :: remapCS !< ALE control structure for remapping
  type(time_type) :: Time !< Current Analysis time
  type(diag_ctrl), pointer :: diag_cs=> NULL() !<Pointer to diagnostics control structure
  type(INC_CS) :: INC_CS !< A Structure containing integer file handles for bias adjustment
  integer :: id_inc_t !< A diagnostic handle for the temperature climatological adjustment
  integer :: id_inc_s !< A diagnostic handle for the salinity climatological adjustment
  integer :: answer_date    !< The vintage of the order of arithmetic and expressions in the
                            !! remapping invoked by the ODA driver.  Values below 20190101 recover
                            !! the answers from the end of 2018, while higher values use updated
                            !! and more robust forms of the same expressions.
  logical :: reproduce_2018_nmme !< true if reproducing older NMME answers.
end type ODA_CS


!>@{  DA parameters
integer, parameter :: NO_ASSIM = 0, OI_ASSIM=1, EAKF_ASSIM=2
!>@}
character(len=40)  :: mdl = "MOM_oda_driver" !< This module's name.


  interface
module subroutine init_oda(Time, G, GV, US, diag_CS, CS)

  type(time_type), intent(in) :: Time !< The current model time.
  type(ocean_grid_type), pointer :: G !< domain and grid information for ocean model
  type(verticalGrid_type), intent(in) :: GV   !< The ocean's vertical grid structure
  type(unit_scale_type),   intent(in) :: US   !< A dimensional unit scaling type
  type(diag_ctrl), target, intent(inout) :: diag_CS !< A pointer to a diagnostic control structure
  type(ODA_CS), pointer, intent(inout) :: CS  !< The DA control structure

! Local variables


end subroutine init_oda
module subroutine set_prior_tracer(Time, G, GV, h, tv, CS)
  type(time_type), intent(in)    :: Time !< The current model time
  type(ocean_grid_type), pointer :: G !< domain and grid information for ocean model
  type(verticalGrid_type),               intent(in)    :: GV   !< The ocean's vertical grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in) :: h   !< Layer thicknesses [H ~> m or kg m-2]
  type(thermo_var_ptrs),                 intent(in) :: tv   !< A structure pointing to various thermodynamic variables

  type(ODA_CS), pointer :: CS !< ocean DA control structure

  ! return if not time for analysis
end subroutine set_prior_tracer
module subroutine get_posterior_tracer(Time, CS, increment)
  type(time_type), intent(in) :: Time !< the current model time
  type(ODA_CS), pointer :: CS !< ocean DA control structure
  logical, optional, intent(in) :: increment !< True if returning increment only



  ! return if not analysis time (retain pointers for h and tv)
end subroutine get_posterior_tracer
module subroutine oda(Time, CS)
  type(time_type), intent(in) :: Time !< the current model time
  type(oda_CS), pointer :: CS !< A pointer the ocean DA control structure

end subroutine oda
module subroutine get_bias_correction_tracer(Time, US, CS)
  type(time_type), intent(in) :: Time !< the current model time
  type(unit_scale_type), intent(in) :: US !< A dimensional unit scaling type
  type(ODA_CS), pointer :: CS !< ocean DA control structure

  ! Local variables
                                                    ! and input-file vertical levels [nondim]


end subroutine get_bias_correction_tracer
module subroutine oda_end(CS)
  type(ODA_CS), intent(inout) :: CS !< the ocean DA control structure

end subroutine oda_end
module subroutine init_ocean_ensemble(CS, Grid, GV, ens_size)
  type(ocean_control_struct), pointer :: CS !< Pointer to ODA control structure
  type(ocean_grid_type), pointer :: Grid !< Pointer to ocean analysis grid
  type(verticalGrid_type), pointer :: GV !< Pointer to DA vertical grid
  integer, intent(in) :: ens_size !< ensemble size


end subroutine init_ocean_ensemble
module subroutine set_analysis_time(Time, CS)
  type(time_type), intent(in) :: Time !< the current model time
  type(ODA_CS), pointer, intent(inout) :: CS !< the DA control structure


end subroutine set_analysis_time
module subroutine apply_oda_tracer_increments(dt, Time_end, G, GV, tv, h, CS)
  real,                     intent(in)    :: dt !< The tracer timestep [T ~> s]
  type(time_type), intent(in)             :: Time_end !< Time at the end of the interval
  type(ocean_grid_type),    intent(in)    :: G  !< ocean grid structure
  type(verticalGrid_type),  intent(in)    :: GV !< The ocean's vertical grid structure
  type(thermo_var_ptrs),    intent(inout) :: tv !< A structure pointing to various thermodynamic variables
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                            intent(in)    :: h  !< layer thickness [H ~> m or kg m-2]
  type(ODA_CS), pointer                   :: CS !< the data assimilation structure

  !! local variables
                                                    !! tendency [C T-1 ~> degC s-1]
                                                    !! tendency [S T-1 ~> ppt s-1]
                                                           !! DA [C T-1 ~> degC s-1]
                                                          !! [S T-1 ~> ppt s-1]

end subroutine apply_oda_tracer_increments
  module subroutine set_up_global_tgrid(T_grid, CS, G)
    type(grid_type), pointer :: T_grid !< global tracer grid
    type(ODA_CS), pointer, intent(in) :: CS !< A pointer to DA control structure.
    type(ocean_grid_type), pointer :: G !< domain and grid information for ocean model

    ! local variables
                   ! global domain [H ~> m or kg m-2]

    !    get global grid information from ocean_model
  end subroutine set_up_global_tgrid
  end interface

end module MOM_oda_driver_mod
