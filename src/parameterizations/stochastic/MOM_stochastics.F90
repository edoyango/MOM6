! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> Top-level module for the MOM6 ocean model in coupled mode.
module MOM_stochastics

! This is the top level module for the MOM6 ocean model.  It contains routines
! for initialization, update, and writing restart of stochastic physics. This
! particular version wraps all of the calls for MOM6 in the calls that had
! been used for MOM4.
!
use MOM_debugging,           only : hchksum, uvchksum, qchksum
use MOM_diag_mediator,       only : register_diag_field, diag_ctrl, time_type, post_data
use MOM_diag_mediator,       only : register_static_field, enable_averages, disable_averaging
use MOM_grid,                only : ocean_grid_type
use MOM_variables,           only : thermo_var_ptrs
use MOM_domains,             only : pass_var, pass_vector, CORNER, SCALAR_PAIR
use MOM_verticalGrid,        only : verticalGrid_type
use MOM_error_handler,       only : MOM_error, MOM_mesg, FATAL, WARNING, is_root_pe
use MOM_error_handler,       only : callTree_enter, callTree_leave
use MOM_file_parser,         only : get_param, log_version, close_param_file, param_file_type
use mpp_domains_mod,         only : domain2d, mpp_get_layout, mpp_get_global_domain
use mpp_domains_mod,         only : mpp_define_domains, mpp_get_compute_domain, mpp_get_data_domain
use MOM_domains,             only : root_PE, num_PEs
use MOM_coms,                only : Get_PElist
use MOM_EOS,                 only : calculate_density, EOS_domain
use stochastic_physics,      only : init_stochastic_physics_ocn, run_stochastic_physics_ocn

#include <MOM_memory.h>

implicit none ; private

public stochastics_init, update_stochastics, apply_skeb

!> This control structure holds parameters for the MOM_stochastics module
type, public:: stochastic_CS
  logical :: do_sppt         !< If true, stochastically perturb the diabatic
  logical :: do_skeb         !< If true, stochastically perturb the horizontal velocity
  logical :: skeb_use_gm     !< If true, adds GM work to the amplitude of SKEBS
  logical :: skeb_use_frict  !< If true, adds viscous dissipation rate to the amplitude of SKEBS
  logical :: pert_epbl       !< If true, then randomly perturb the KE dissipation and genration terms
  integer :: id_sppt_wts    = -1 !< Diagnostic id for SPPT
  integer :: id_skeb_wts    = -1 !< Diagnostic id for SKEB
  integer :: id_skebu       = -1 !< Diagnostic id for SKEB
  integer :: id_skebv       = -1 !< Diagnostic id for SKEB
  integer :: id_diss        = -1 !< Diagnostic id for SKEB
  integer :: skeb_npass     = -1 !< number of passes of the 9-point smoother for the dissipation estimate
  integer :: id_psi         = -1 !< Diagnostic id for SPPT
  integer :: id_epbl1_wts   = -1 !< Diagnostic id for epbl generation perturbation
  integer :: id_epbl2_wts   = -1 !< Diagnostic id for epbl dissipation perturbation
  integer :: id_skeb_taperu = -1 !< Diagnostic id for u taper of SKEB velocity increment
  integer :: id_skeb_taperv = -1 !< Diagnostic id for v taper of SKEB velocity increment
  real    :: skeb_gm_coef     !< If skeb_use_gm is true, then skeb_gm_coef * GM_work is added to the
                              !! dissipation rate used to set the amplitude of SKEBS [nondim]
  real    :: skeb_frict_coef  !< If skeb_use_frict is true, then skeb_gm_coef * GM_work is added to the
                              !! dissipation rate used to set the amplitude of SKEBS [nondim]
  real, allocatable :: skeb_diss(:,:,:) !< Dissipation rate used to set amplitude of SKEBS [L2 T-3 ~> m2 s-3]
                                        !! Index into this at h points.
  ! stochastic patterns
  real, allocatable :: sppt_wts(:,:)  !< Random pattern for ocean SPPT
                                      !! tendencies with a number between 0 and 2 [nondim]
  real, allocatable :: skeb_wts(:,:)  !< Random pattern for ocean SKEB [nondim]
  real, allocatable :: epbl1_wts(:,:) !< Random pattern for K.E. generation [nondim]
  real, allocatable :: epbl2_wts(:,:) !< Random pattern for K.E. dissipation [nondim]
  type(time_type), pointer :: Time !< Pointer to model time (needed for sponges)
  type(diag_ctrl), pointer :: diag=>NULL() !< A structure that is used to regulate the

  ! Taper array to smoothly zero out the SKEBS velocity increment near land
  real, allocatable :: taperCu(:,:) !< Taper applied to u component of stochastic
                                    !! velocity increment range [0,1], [nondim]
  real, allocatable :: taperCv(:,:) !< Taper applied to v component of stochastic
                                    !! velocity increment range [0,1], [nondim]

end type stochastic_CS


  interface
module subroutine stochastics_init(dt, grid, GV, CS, param_file, diag, Time)
  real, intent(in)                       :: dt      !< time step [T ~> s]
  type(ocean_grid_type),   intent(in)    :: grid    !< horizontal grid information
  type(verticalGrid_type), intent(in)    :: GV      !< vertical grid structure
  type(stochastic_CS), pointer, intent(inout) :: CS !< stochastic control structure
  type(param_file_type),   intent(in)    :: param_file !< A structure to parse for run-time parameters
  type(diag_ctrl), target, intent(inout) :: diag    !< structure to regulate diagnostic output
  type(time_type), target                :: Time    !< model time

  ! Local variables
                               ! increments to 0 at the boundary.

  ! This include declares and sets the variable "version".

end subroutine stochastics_init
module subroutine update_stochastics(CS)
  type(stochastic_CS),      intent(inout) :: CS        !< diabatic control structure
end subroutine update_stochastics
module subroutine apply_skeb(grid,GV,CS,uc,vc,thickness,tv,dt,Time_end)

  type(ocean_grid_type),   intent(in)    :: grid   !< ocean grid structure
  type(verticalGrid_type), intent(in)    :: GV     !< ocean vertical grid
  type(stochastic_CS),     intent(inout) :: CS     !< stochastic control structure

  real, dimension(SZIB_(grid),SZJ_(grid),SZK_(GV)), intent(inout) :: uc        !< zonal velocity [L T-1 ~> m s-1]
  real, dimension(SZI_(grid),SZJB_(grid),SZK_(GV)), intent(inout) :: vc        !< meridional velocity [L T-1 ~> m s-1]
  real, dimension(SZI_(grid),SZJ_(grid),SZK_(GV)),  intent(in)    :: thickness !< thickness [H ~> m or kg m-2]
  type(thermo_var_ptrs),                            intent(in)    :: tv       !< points to thermodynamic fields
  real,                                       intent(in)    :: dt       !< time increment [T ~> s]
  type(time_type),                            intent(in)    :: Time_end !< Time at the end of the interval
! locals

                                                                   !! [L2 T-1 ~> m2 s-1]
                                                                   !! [L2 T-3 ~> m2 s-2]
                                                                   !! [L2 ~> m2]


end subroutine apply_skeb
module subroutine smooth_x9_uv(G, field_u, field_v, zero_land)
  type(ocean_grid_type),             intent(in)    :: G         !< Ocean grid
  real, dimension(SZIB_(G),SZJ_(G)), intent(inout) :: field_u   !< u-point field to be smoothed[arbitrary]
  real, dimension(SZI_(G),SZJB_(G)), intent(inout) :: field_v   !< v-point field to be smoothed [arbitrary]
  logical,                 optional, intent(in)    :: zero_land !< If present and false, return the average
                                                                !! of the surrounding ocean points when
                                                                !! smoothing, otherwise use a value of 0 for
                                                                !! land points and include them in the averages.

  ! Local variables.

end subroutine smooth_x9_uv
  end interface

end module MOM_stochastics
