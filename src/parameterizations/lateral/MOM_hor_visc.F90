! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> Calculates horizontal viscosity and viscous stresses
module MOM_hor_visc

use MOM_checksums,             only : hchksum, Bchksum, uvchksum
use MOM_coms,                  only : min_across_PEs
use MOM_diag_mediator,         only : post_data, register_diag_field, safe_alloc_ptr
use MOM_diag_mediator,         only : post_product_u, post_product_sum_u
use MOM_diag_mediator,         only : post_product_v, post_product_sum_v
use MOM_diag_mediator,         only : diag_ctrl, time_type
use MOM_domains,               only : pass_var, CORNER, pass_vector, AGRID, BGRID_NE
use MOM_domains,               only : To_All, Scalar_Pair
use MOM_error_handler,         only : MOM_error, FATAL, WARNING, is_root_pe
use MOM_file_parser,           only : get_param, log_version, param_file_type
use MOM_grid,                  only : ocean_grid_type
use MOM_interface_heights,     only : thickness_to_dz
use MOM_lateral_mixing_coeffs, only : VarMix_CS, calc_QG_slopes, calc_QG_Leith_viscosity
use MOM_barotropic,            only : barotropic_CS, barotropic_get_tav
use MOM_thickness_diffuse,     only : thickness_diffuse_CS, thickness_diffuse_get_KH
use MOM_io,                    only : MOM_read_data, slasher
use MOM_MEKE_types,            only : MEKE_type
use MOM_open_boundary,         only : ocean_OBC_type, OBC_DIRECTION_E, OBC_DIRECTION_W
use MOM_open_boundary,         only : OBC_DIRECTION_N, OBC_DIRECTION_S
use MOM_open_boundary,         only : OBC_STRAIN_NONE, OBC_STRAIN_ZERO, OBC_STRAIN_FREESLIP
use MOM_open_boundary,         only : OBC_STRAIN_COMPUTED, OBC_STRAIN_SPECIFIED
use MOM_stochastics,           only : stochastic_CS
use MOM_unit_scaling,          only : unit_scale_type
use MOM_verticalGrid,          only : verticalGrid_type
use MOM_variables,             only : accel_diag_ptrs, thermo_var_ptrs
use MOM_Zanna_Bolton,          only : ZB2020_lateral_stress, ZB2020_init, ZB2020_end
use MOM_Zanna_Bolton,          only : ZB2020_CS, ZB2020_copy_gradient_and_thickness

implicit none ; private

#include <MOM_memory.h>

public horizontal_viscosity, hor_visc_init, hor_visc_end, hor_visc_vel_stencil

!> Control structure for horizontal viscosity
type, public :: hor_visc_CS ; private
  logical :: initialized = .false. !< True if this control structure has been initialized.
  logical :: Laplacian       !< Use a Laplacian horizontal viscosity if true.
  logical :: biharmonic      !< Use a biharmonic horizontal viscosity if true.
  logical :: debug           !< If true, write verbose checksums for debugging purposes.
  logical :: no_slip         !< If true, no slip boundary conditions are used.
                             !! Otherwise free slip boundary conditions are assumed.
                             !! The implementation of the free slip boundary
                             !! conditions on a C-grid is much cleaner than the
                             !! no slip boundary conditions. The use of free slip
                             !! b.c.s is strongly encouraged. The no slip b.c.s
                             !! are not implemented with the biharmonic viscosity.
  logical :: bound_Kh        !< If true, the Laplacian coefficient is locally
                             !! limited to guarantee stability.
  logical :: EY24_EBT_BS     !! If true, use an equivalent barotropic backscatter
                             !! with a stabilizing kill switch in MEKE,
                             !< developed by Yankovsky et al. 2024
  logical :: bound_Ah        !< If true, the biharmonic coefficient is locally
                             !! limited to guarantee stability.
  real    :: Re_Ah           !! If nonzero, the biharmonic coefficient is scaled
                             !< so that the biharmonic Reynolds number is equal to this [nondim].
  real    :: bound_coef      !< The nondimensional coefficient of the ratio of
                             !! the viscosity bounds to the theoretical maximum
                             !! for stability without considering other terms [nondim].
                             !! The default is 0.8.
  real    :: KS_coef         !< A nondimensional coefficient on the biharmonic viscosity that sets the
                             !! kill switch for backscatter. Default is 1.0 [nondim].
  real    :: KS_timescale    !< A timescale for computing CFL limit for turning off backscatter [T ~> s].
  logical :: backscatter_underbound !< If true, the bounds on the biharmonic viscosity are allowed
                             !! to increase where the Laplacian viscosity is negative (due to
                             !! backscatter parameterizations) beyond the largest timestep-dependent
                             !! stable values of biharmonic viscosity when no Laplacian viscosity is
                             !! applied.  The default is true for historical reasons, but this option
                             !! probably should not be used as it can lead to numerical instabilities.
  logical :: Smagorinsky_Kh  !< If true, use Smagorinsky nonlinear eddy
                             !! viscosity. KH is the background value.
  logical :: Smagorinsky_Ah  !< If true, use a biharmonic form of Smagorinsky
                             !! nonlinear eddy viscosity. AH is the background.
  logical :: Leith_Kh        !< If true, use 2D Leith nonlinear eddy
                             !! viscosity. KH is the background value.
  logical :: Modified_Leith  !< If true, use extra component of Leith viscosity
                             !! to damp divergent flow. To use, still set Leith_Kh=.TRUE.
  logical :: use_beta_in_Leith !< If true, includes the beta term in the Leith viscosity
  logical :: Leith_Ah        !< If true, use a biharmonic form of 2D Leith
                             !! nonlinear eddy viscosity. AH is the background.
  logical :: use_Leithy      !< If true, use a biharmonic form of 2D Leith
                             !! nonlinear eddy viscosity with harmonic backscatter.
                             !! Ah is the background. Leithy = Leith+E
  real    :: c_K             !< Fraction of energy dissipated by the biharmonic term
                             !! that gets backscattered in the Leith+E scheme. [nondim]
  logical :: smooth_Ah       !< If true (default), then Ah and m_leithy are smoothed.
                             !! This smoothing requires a lot of blocking communication.
  logical :: use_QG_Leith_visc    !< If true, use QG Leith nonlinear eddy viscosity.
                             !! KH is the background value.
  logical :: bound_Coriolis  !< If true & SMAGORINSKY_AH is used, the biharmonic
                             !! viscosity is modified to include a term that
                             !! scales quadratically with the velocity shears.
  logical :: use_Kh_bg_2d    !< Read 2d background viscosity from a file.
  logical :: Kh_bg_2d_bug    !< If true, retain an answer-changing horizontal indexing bug
                             !! in setting the corner-point viscosities when USE_KH_BG_2D=True.
  real    :: Kh_bg_min       !< The minimum value allowed for Laplacian horizontal
                             !! viscosity [L2 T-1 ~> m2 s-1]. The default is 0.0.
  logical :: FrictWork_bug   !< If true, retain an answer-changing bug in calculating FrictWork,
                             !! which cancels the h in thickness flux and the h at velocity point.
  logical :: OBC_strain_bug  !< If true, recover a bug that specified shear strain option at open
                             !! boundaries cannot be applied.
  logical :: use_land_mask   !< Use the land mask for the computation of thicknesses
                             !! at velocity locations. This eliminates the dependence on
                             !! arbitrary values over land or outside of the domain.
                             !! Default is False to maintain answers with legacy experiments
                             !! but should be changed to True for new experiments.
  logical :: anisotropic     !< If true, allow anisotropic component to the viscosity.
  logical :: add_LES_viscosity!< If true, adds the viscosity from Smagorinsky and Leith to
                             !! the background viscosity instead of taking the maximum.
  real    :: Kh_aniso        !< The anisotropic viscosity [L2 T-1 ~> m2 s-1].
  logical :: dynamic_aniso   !< If true, the anisotropic viscosity is recomputed as a function
                             !! of state. This is set depending on ANISOTROPIC_MODE.
  logical :: res_scale_MEKE  !< If true, the viscosity contribution from MEKE is scaled by
                             !! the resolution function.
  logical :: use_GME         !< If true, use GME backscatter scheme.
  integer :: answer_date     !< The vintage of the order of arithmetic and expressions in the
                             !! horizontal viscosity calculations.  Values below 20190101 recover
                             !! the answers from the end of 2018, while higher values use updated
                             !! and more robust forms of the same expressions.
  real    :: GME_h0          !< The strength of GME tapers quadratically to zero when the bathymetric
                             !! total water column thickness is less than GME_H0 [H ~> m or kg m-2]
  real    :: GME_efficiency  !< The nondimensional prefactor multiplying the GME coefficient [nondim]
  real    :: GME_limiter     !< The absolute maximum value the GME coefficient is allowed to take [L2 T-1 ~> m2 s-1].
  real    :: min_grid_Kh     !< Minimum horizontal Laplacian viscosity used to
                             !! limit the grid Reynolds number [L2 T-1 ~> m2 s-1]
  real    :: min_grid_Ah     !< Minimun horizontal biharmonic viscosity used to
                             !! limit grid Reynolds number [L4 T-1 ~> m4 s-1]
  logical :: use_cont_thick  !< If true, thickness at velocity points adopts h[uv] in BT_cont from continuity solver.
  logical :: use_cont_thick_bug  !< If true, retain an answer-changing bug for thickness at velocity points.
  type(ZB2020_CS) :: ZB2020  !< Zanna-Bolton 2020 control structure.
  logical :: use_ZB2020      !< If true, use Zanna-Bolton 2020 parameterization.
  logical :: use_circulation !< If true, use circulation theorem to compute vorticity (for ZB20 or Leith)

  real ALLOCABLE_, dimension(NIMEM_,NJMEM_) :: Kh_bg_xx
                      !< The background Laplacian viscosity at h points [L2 T-1 ~> m2 s-1].
                      !! The actual viscosity may be the larger of this
                      !! viscosity and the Smagorinsky and Leith viscosities.
  real, allocatable :: Kh_bg_2d(:,:)
                      !< The background Laplacian viscosity at h points [L2 T-1 ~> m2 s-1].
                      !! The actual viscosity may be the larger of this
                      !! viscosity and the Smagorinsky and Leith viscosities.
  real ALLOCABLE_, dimension(NIMEM_,NJMEM_) :: Ah_bg_xx
                      !< The background biharmonic viscosity at h points [L4 T-1 ~> m4 s-1].
                      !! The actual viscosity may be the larger of this
                      !! viscosity and the Smagorinsky and Leith viscosities.
  real ALLOCABLE_, dimension(NIMEM_,NJMEM_) :: reduction_xx
                      !< The amount by which stresses through h points are reduced
                      !! due to partial barriers [nondim].
  real, allocatable :: Kh_Max_xx(:,:)     !< The maximum permitted Laplacian viscosity [L2 T-1 ~> m2 s-1].
  real, allocatable :: Ah_Max_xx(:,:)     !< The maximum permitted biharmonic viscosity [L4 T-1 ~> m4 s-1].
  real, allocatable :: Ah_Max_xx_KS(:,:)  !< The maximum permitted biharmonic viscosity for
                                          !! the kill switch [L4 T-1 ~> m4 s-1].
  real, allocatable :: n1n2_h(:,:)        !< Factor n1*n2 in the anisotropic direction tensor at h-points [nondim]
  real, allocatable :: n1n1_m_n2n2_h(:,:) !< Factor n1**2-n2**2 in the anisotropic direction tensor at h-points [nondim]
  real ALLOCABLE_, dimension(NIMEM_,NJMEM_) :: &
    grid_sp_h2,     & !< Harmonic mean of the squares of the grid [L2 ~> m2]
    grid_sp_h3        !< Harmonic mean of the squares of the grid^(3/2) [L3 ~> m3]
  real ALLOCABLE_, dimension(NIMEMB_PTR_,NJMEMB_PTR_) :: Kh_bg_xy
                      !< The background Laplacian viscosity at q points [L2 T-1 ~> m2 s-1].
                      !! The actual viscosity may be the larger of this
                      !! viscosity and the Smagorinsky and Leith viscosities.
  real ALLOCABLE_, dimension(NIMEMB_PTR_,NJMEMB_PTR_) :: Ah_bg_xy
                      !< The background biharmonic viscosity at q points [L4 T-1 ~> m4 s-1].
                      !! The actual viscosity may be the larger of this
                      !! viscosity and the Smagorinsky and Leith viscosities.
  real ALLOCABLE_, dimension(NIMEMB_PTR_,NJMEMB_PTR_) :: reduction_xy
                      !< The amount by which stresses through q points are reduced
                      !! due to partial barriers [nondim].
  real, allocatable :: Kh_Max_xy(:,:)  !< The maximum permitted Laplacian viscosity [L2 T-1 ~> m2 s-1].
  real, allocatable :: Ah_Max_xy(:,:)     !< The maximum permitted biharmonic viscosity [L4 T-1 ~> m4 s-1].
  real, allocatable :: Ah_Max_xy_KS(:,:)  !< The maximum permitted biharmonic viscosity for
                                          !! the  kill switch [L4 T-1 ~> m4 s-1].
  real, allocatable :: n1n2_q(:,:)        !< Factor n1*n2 in the anisotropic direction tensor at q-points [nondim]
  real, allocatable :: n1n1_m_n2n2_q(:,:) !< Factor n1**2-n2**2 in the anisotropic direction tensor at q-points [nondim]

  real ALLOCABLE_, dimension(NIMEM_,NJMEM_) :: &
    dx2h,           & !< Pre-calculated dx^2 at h points [L2 ~> m2]
    dy2h,           & !< Pre-calculated dy^2 at h points [L2 ~> m2]
    dx_dyT,         & !< Pre-calculated dx/dy at h points [nondim]
    dy_dxT            !< Pre-calculated dy/dx at h points [nondim]
  real, allocatable :: m_const_leithy(:,:) !< Pre-calculated .5*sqrt(c_K)*max{dx,dy} [L ~> m]
  real, allocatable :: m_leithy_max(:,:)   !< Pre-calculated 4./max(dx,dy)^2 at h points [L-2 ~> m-2]
  real ALLOCABLE_, dimension(NIMEMB_PTR_,NJMEMB_PTR_) :: &
    dx2q,    & !< Pre-calculated dx^2 at q points [L2 ~> m2]
    dy2q,    & !< Pre-calculated dy^2 at q points [L2 ~> m2]
    dx_dyBu, & !< Pre-calculated dx/dy at q points [nondim]
    dy_dxBu    !< Pre-calculated dy/dx at q points [nondim]
  real ALLOCABLE_, dimension(NIMEMB_PTR_,NJMEM_) :: &
    Idx2dyCu, & !< 1/(dx^2 dy) at u points [L-3 ~> m-3]
    Idxdy2u     !< 1/(dx dy^2) at u points [L-3 ~> m-3]
  real ALLOCABLE_, dimension(NIMEM_,NJMEMB_PTR_) :: &
    Idx2dyCv, & !< 1/(dx^2 dy) at v points [L-3 ~> m-3]
    Idxdy2v     !< 1/(dx dy^2) at v points [L-3 ~> m-3]

  ! The following variables are precalculated time-invariant combinations of
  ! parameters and metric terms.
  real, allocatable :: Laplac2_const_xx(:,:) !< Laplacian metric-dependent constants [L2 ~> m2]
  real, allocatable :: Biharm6_const_xx(:,:) !< Biharmonic metric-dependent constants [L6 ~> m6]
  real, allocatable :: Laplac3_const_xx(:,:) !< Laplacian metric-dependent constants [L3 ~> m3]
  real, allocatable :: Biharm_const_xx(:,:)  !< Biharmonic metric-dependent constants [L4 ~> m4]
  real, allocatable :: Biharm_const2_xx(:,:) !< Biharmonic metric-dependent constants [T L4 ~> s m4]
  real, allocatable :: Re_Ah_const_xx(:,:)   !< Biharmonic metric-dependent constants [L3 ~> m3]

  real, allocatable :: Laplac2_const_xy(:,:) !< Laplacian metric-dependent constants [L2 ~> m2]
  real, allocatable :: Biharm6_const_xy(:,:) !< Biharmonic metric-dependent constants [L6 ~> m6]
  real, allocatable :: Laplac3_const_xy(:,:) !< Laplacian metric-dependent constants [L3 ~> m3]
  real, allocatable :: Biharm_const_xy(:,:)  !< Biharmonic metric-dependent constants [L4 ~> m4]
  real, allocatable :: Biharm_const2_xy(:,:) !< Biharmonic metric-dependent constants [T L4 ~> s m4]
  real, allocatable :: Re_Ah_const_xy(:,:)   !< Biharmonic metric-dependent constants [L3 ~> m3]

  type(diag_ctrl), pointer :: diag => NULL() !< structure to regulate diagnostics

  ! real, allocatable :: hf_diffu(:,:,:)  ! Zonal horizontal viscous acceleleration times
  !                                       ! fractional thickness [L T-2 ~> m s-2].
  ! real, allocatable :: hf_diffv(:,:,:)  ! Meridional horizontal viscous acceleleration times
  !                                       ! fractional thickness [L T-2 ~> m s-2].
  ! 3D diagnostics hf_diffu(diffv) are commented because there is no clarity on proper remapping grid option.
  ! The code is retained for debugging purposes in the future.

  integer :: num_smooth_gme !< number of smoothing passes for the GME fluxes.
  !>@{
  !! Diagnostic id
  integer :: id_grid_Re_Ah = -1, id_grid_Re_Kh   = -1
  integer :: id_diffu     = -1, id_diffv         = -1
  ! integer :: id_hf_diffu  = -1, id_hf_diffv      = -1
  integer :: id_h_diffu  = -1, id_h_diffv      = -1
  integer :: id_hf_diffu_2d = -1, id_hf_diffv_2d = -1
  integer :: id_intz_diffu_2d = -1, id_intz_diffv_2d = -1
  integer :: id_diffu_visc_rem = -1, id_diffv_visc_rem = -1
  integer :: id_Ah_h      = -1, id_Ah_q          = -1
  integer :: id_Kh_h      = -1, id_Kh_q          = -1
  integer :: id_GME_coeff_h = -1, id_GME_coeff_q = -1
  integer :: id_dudx_bt = -1, id_dvdy_bt = -1
  integer :: id_dudy_bt = -1, id_dvdx_bt = -1
  integer :: id_vort_xy_q = -1, id_div_xx_h      = -1
  integer :: id_sh_xy_q = -1,    id_sh_xx_h      = -1
  integer :: id_FrictWork = -1, id_FrictWorkIntz = -1
  integer :: id_FrictWork_bh = -1, id_FrictWorkIntz_bh = -1
  integer :: id_FrictWork_GME = -1
  integer :: id_normstress = -1, id_shearstress = -1
  integer :: id_visc_limit_h = -1, id_visc_limit_q = -1
  integer :: id_visc_limit_h_flag = -1, id_visc_limit_q_flag = -1
  integer :: id_visc_limit_h_frac = -1, id_visc_limit_q_frac = -1
  integer :: id_BS_coeff_h = -1, id_BS_coeff_q = -1
  !>@}

end type hor_visc_CS


  interface
module subroutine horizontal_viscosity(u, v, h, uh, vh, diffu, diffv, MEKE, VarMix, G, GV, US, &
                                CS, tv, dt, OBC, BT, TD, ADp, hu_cont, hv_cont, STOCH)
  type(ocean_grid_type),         intent(in)  :: G      !< The ocean's grid structure.
  type(verticalGrid_type),       intent(in)  :: GV     !< The ocean's vertical grid structure.
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                                 intent(in)  :: u      !< The zonal velocity [L T-1 ~> m s-1].
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                                 intent(in)  :: v      !< The meridional velocity [L T-1 ~> m s-1].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                                 intent(inout) :: h    !< Layer thicknesses [H ~> m or kg m-2].
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                                 intent(in)  :: uh      !< The zonal volume transport [H L2 T-1 ~> m3 s-1].
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                                 intent(in)  :: vh      !< The meridional volume transport [H L2 T-1 ~> m3 s-1].
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                                 intent(out) :: diffu  !< Zonal acceleration due to convergence of
                                                       !! along-coordinate stress tensor [L T-2 ~> m s-2]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                                 intent(out) :: diffv  !< Meridional acceleration due to convergence
                                                       !! of along-coordinate stress tensor [L T-2 ~> m s-2].
  type(MEKE_type),               intent(inout) :: MEKE !< MEKE fields
                                                       !! related to Mesoscale Eddy Kinetic Energy.
  type(VarMix_CS),               intent(inout) :: VarMix !< Variable mixing control structure
  type(unit_scale_type),         intent(in)    :: US   !< A dimensional unit scaling type
  type(hor_visc_CS),             intent(inout) :: CS   !< Horizontal viscosity control structure
  type(thermo_var_ptrs),         intent(in)    :: tv   !< A structure pointing to various
                                                       !! thermodynamic variables
  real,                          intent(in)    :: dt   !< Time increment [T ~> s]
  type(ocean_OBC_type), optional, pointer      :: OBC  !< Pointer to an open boundary condition type
  type(barotropic_CS), optional, intent(in)    :: BT   !< Barotropic control structure
  type(thickness_diffuse_CS), optional, intent(in) :: TD !< Thickness diffusion control structure
  type(accel_diag_ptrs), optional, intent(in)  :: ADp  !< Acceleration diagnostics
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                    optional, intent(inout) :: hu_cont !< Layer thickness at u-points [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                    optional, intent(inout) :: hv_cont !< Layer thickness at v-points [H ~> m or kg m-2].
  type(stochastic_CS), intent(inout), optional :: STOCH !< Stochastic control structure

  ! Local variables
end subroutine horizontal_viscosity
module subroutine hor_visc_init(Time, G, GV, US, param_file, diag, CS, ADp)
  type(time_type),         intent(in)    :: Time !< Current model time.
  type(ocean_grid_type),   intent(inout) :: G    !< The ocean's grid structure.
  type(verticalGrid_type), intent(in)    :: GV   !< The ocean's vertical grid structure
  type(unit_scale_type),   intent(in)    :: US   !< A dimensional unit scaling type
  type(param_file_type),   intent(in)    :: param_file !< A structure to parse for run-time
                                                 !! parameters.
  type(diag_ctrl), target, intent(inout) :: diag !< Structure to regulate diagnostic output.
  type(hor_visc_CS),       intent(inout) :: CS   !< Horizontal viscosity control structure
  type(accel_diag_ptrs), intent(in), optional :: ADp !< Acceleration diagnostics

  ! u0v is the Laplacian sensitivities to the v velocities at u points, with u0u, v0u, and v0v defined analogously.
                           ! grid spacing, to limit Laplacian viscosity.
                           ! vorticity points around a thickness point [T-1 ~> s-1]
                           ! [T2 L-2 ~> s2 m-2]
                           ! grid spacing, to limit biharmonic viscosity
                           ! the quadratically varying biharmonic viscosity
                           ! balances Coriolis acceleration [L T-1 ~> m s-1]
                           ! If false and USE_GME = True, issue a FATAL error.
                           ! recreate the bugs, or if false bugs are only used if actively selected.
  ! This include declares and sets the variable "version".
end subroutine hor_visc_init
module function hor_visc_vel_stencil(CS) result(stencil)
  type(hor_visc_CS), intent(in) :: CS !< Control structure for horizontal viscosity
  integer ::  stencil !< The horizontal viscosity velocity stencil size with the current settings.

end function hor_visc_vel_stencil
module subroutine align_aniso_tensor_to_grid(CS, n1, n2)
  type(hor_visc_CS), intent(inout) :: CS !< Control structure for horizontal viscosity
  real,              intent(in) :: n1 !< i-component of direction vector [nondim]
  real,              intent(in) :: n2 !< j-component of direction vector [nondim]
  ! Local variables
  ! For normalizing n=(n1,n2) in case arguments are not a unit vector
end subroutine align_aniso_tensor_to_grid
module subroutine smooth_GME(CS, G, GME_flux_h, GME_flux_q)
  type(hor_visc_CS),                            intent(in)    :: CS        !< Control structure
  type(ocean_grid_type),                        intent(in)    :: G         !< Ocean grid
  real, dimension(SZI_(G),SZJ_(G)),   optional, intent(inout) :: GME_flux_h!< GME diffusive flux
                                                              !! at h points [L2 T-2 ~> m2 s-2]
  real, dimension(SZIB_(G),SZJB_(G)), optional, intent(inout) :: GME_flux_q!< GME diffusive flux
                                                              !! at q points [L2 T-2 ~> m2 s-2]
  ! local variables

end subroutine smooth_GME
module subroutine smooth_x9_h(G, field_h, zero_land)
  type(ocean_grid_type),            intent(in)    :: G         !< Ocean grid
  real, dimension(SZI_(G),SZJ_(G)), intent(inout) :: field_h   !< h-point field to be smoothed [arbitrary]
  logical,                optional, intent(in)    :: zero_land !< If present and false, return the average
                                                               !! of the surrounding ocean points when
                                                               !! smoothing, otherwise use a value of 0 for
                                                               !! land points and include them in the averages.

  ! Local variables

end subroutine smooth_x9_h
module subroutine smooth_x9_uv(G, field_u, field_v, zero_land)
  type(ocean_grid_type),             intent(in)    :: G         !< Ocean grid
  real, dimension(SZIB_(G),SZJ_(G)), intent(inout) :: field_u   !< u-point field to be smoothed [arbitrary]
  real, dimension(SZI_(G),SZJB_(G)), intent(inout) :: field_v   !< v-point field to be smoothed [arbitrary]
  logical,                 optional, intent(in)    :: zero_land !< If present and false, return the average
                                                                !! of the surrounding ocean points when
                                                                !! smoothing, otherwise use a value of 0 for
                                                                !! land points and include them in the averages.

  ! Local variables.

end subroutine smooth_x9_uv
module subroutine hor_visc_end(CS)
  type(hor_visc_CS), intent(inout) :: CS !< Horizontal viscosity control structure
end subroutine hor_visc_end
  end interface

end module MOM_hor_visc
