! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> Isopycnal height diffusion (or Gent McWilliams diffusion)
module MOM_thickness_diffuse

use MOM_debugging,             only : hchksum, uvchksum
use MOM_diag_mediator,         only : post_data, query_averaging_enabled, diag_ctrl
use MOM_diag_mediator,         only : register_diag_field, safe_alloc_ptr, time_type
use MOM_diag_mediator,         only : diag_update_remap_grids
use MOM_domains,               only : pass_var, CORNER, pass_vector
use MOM_error_handler,         only : MOM_error, FATAL, WARNING, is_root_pe
use MOM_EOS,                   only : calculate_density, calculate_density_derivs, EOS_domain
use MOM_EOS,                   only : calculate_density_second_derivs
use MOM_file_parser,           only : get_param, log_version, param_file_type
use MOM_grid,                  only : ocean_grid_type
use MOM_io,                    only : MOM_read_data, slasher
use MOM_interface_heights,     only : find_eta, thickness_to_dz
use MOM_isopycnal_slopes,      only : vert_fill_TS
use MOM_lateral_mixing_coeffs, only : VarMix_CS
use MOM_MEKE_types,            only : MEKE_type
use MOM_stochastics,           only : stochastic_CS
use MOM_unit_scaling,          only : unit_scale_type
use MOM_variables,             only : thermo_var_ptrs, cont_diag_ptrs
use MOM_verticalGrid,          only : verticalGrid_type
use MOM_meso_sfn_ANN,          only : meso_sfn_ANN_compute, MESO_SFN_ANN_CS
use MOM_meso_sfn_ANN,          only : meso_sfn_ANN_init, meso_sfn_ANN_end

implicit none ; private

#include <MOM_memory.h>

public thickness_diffuse, thickness_diffuse_init, thickness_diffuse_end
public thickness_diffuse_get_KH

! A note on unit descriptions in comments: MOM6 uses units that can be rescaled for dimensional
! consistency testing. These are noted in comments with units like Z, H, L, and T, along with
! their mks counterparts with notation like "a velocity [Z T-1 ~> m s-1]".  If the units
! vary with the Boussinesq approximation, the Boussinesq variant is given first.

!> Control structure for thickness_diffuse
type, public :: thickness_diffuse_CS ; private
  logical :: initialized = .false. !< True if this control structure has been initialized.
  real    :: Khth                !< Background isopycnal depth diffusivity [L2 T-1 ~> m2 s-1]
  real    :: Khth_Slope_Cff      !< Slope dependence coefficient of Khth [nondim]
  real    :: max_Khth_CFL        !< Maximum value of the diffusive CFL for isopycnal height diffusion [nondim]
  real    :: Khth_Min            !< Minimum value of Khth [L2 T-1 ~> m2 s-1]
  real    :: Khth_Max            !< Maximum value of Khth [L2 T-1 ~> m2 s-1], or 0 for no max
  real    :: Kh_eta_bg           !< Background isopycnal height diffusivity [L2 T-1 ~> m2 s-1]
  real    :: Kh_eta_vel          !< Velocity scale that is multiplied by the grid spacing to give
                                 !! the isopycnal height diffusivity [L T-1 ~> m s-1]
  real    :: slope_max           !< Slopes steeper than slope_max are limited in some way [Z L-1 ~> nondim]
  real    :: kappa_smooth        !< Vertical diffusivity used to interpolate more sensible values
                                 !! of T & S into thin layers [H Z T-1 ~> m2 s-1 or kg m-1 s-1]
  logical :: thickness_diffuse   !< If true, interfaces heights are diffused.
  logical :: full_depth_khth_min !< If true, KHTH_MIN is enforced throughout the whole water column.
                                 !! Otherwise, KHTH_MIN is only enforced at the surface. This parameter
                                 !! is only available when KHTH_USE_EBT_STRUCT=True and KHTH_MIN>0.
  logical :: use_FGNV_streamfn   !< If true, use the streamfunction formulation of
                                 !! Ferrari et al., 2010, which effectively emphasizes
                                 !! graver vertical modes by smoothing in the vertical.
  real    :: FGNV_scale          !< A coefficient scaling the vertical smoothing term in the
                                 !! Ferrari et al., 2010, streamfunction formulation [nondim].
  real    :: FGNV_c_min          !< A minimum wave speed used in the Ferrari et al., 2010,
                                 !! streamfunction formulation [L T-1 ~> m s-1].
  real    :: N2_floor            !< A floor for squared buoyancy frequency in the Ferrari et al., 2010,
                                 !! streamfunction formulation divided by aspect ratio rescaling factors
                                 !! [L2 Z-2 T-2 ~> s-2].
  logical :: detangle_interfaces !< If true, add 3-d structured interface height
                                 !! diffusivities to horizontally smooth jagged layers.
  real    :: detangle_time       !< If detangle_interfaces is true, this is the
                                 !! timescale over which maximally jagged grid-scale
                                 !! thickness variations are suppressed [T ~> s].  This must be
                                 !! longer than DT, or 0 (the default) to use DT.
  integer :: nkml                !< number of layers within mixed layer
  logical :: debug               !< write verbose checksums for debugging purposes
  logical :: use_GME_thickness_diffuse !< If true, passes GM coefficients to MOM_hor_visc for use
                                 !! with GME closure.
  logical :: MEKE_GEOMETRIC      !< If true, uses the GM coefficient formulation from the GEOMETRIC
                                 !! framework (Marshall et al., 2012)
  real    :: MEKE_GEOMETRIC_alpha!< The nondimensional coefficient governing the efficiency of
                                 !! the GEOMETRIC isopycnal height diffusion [nondim]
  real    :: MEKE_GEOMETRIC_epsilon !< Minimum Eady growth rate for the GEOMETRIC thickness
                                 !! diffusivity [T-1 ~> s-1].
  integer :: MEKE_GEOM_answer_date  !< The vintage of the expressions in the MEKE_GEOMETRIC
                                 !! calculation.  Values below 20190101 recover the answers from the
                                 !! original implementation, while higher values use expressions that
                                 !! satisfy rotational symmetry.
  logical :: Use_KH_in_MEKE      !< If true, uses the isopycnal height diffusivity calculated here to diffuse MEKE.
  real    :: MEKE_min_depth_diff !< The minimum total depth over which to average the diffusivity
                                 !! used for MEKE [H ~> m or kg m-2].  When the total depth is less
                                 !! than this, the diffusivity is scaled away.
  logical :: GM_src_alt          !< If true, use the GM energy conversion form S^2*N^2*kappa rather
                                 !! than the streamfunction for the GM source term for MEKE.
  integer :: MEKE_src_answer_date  !< The vintage of the expressions in the GM energy conversion
                                 !! calculation when MEKE_GM_SRC_ALT is true.  Values below 20240601
                                 !! recover the answers from the original implementation, while higher
                                 !! values use expressions that satisfy rotational symmetry.
  logical :: MEKE_src_slope_bug  !< If true, use a bug that limits the positive values, but not the
                                 !! negative values, of the slopes used when MEKE_GM_SRC_ALT is true.
                                 !! When this is true, it breaks rotational symmetry.
  logical :: use_GM_work_bug     !< If true, use the incorrect sign for the
                                 !! top-level work tendency on the top layer.
  real :: Stanley_det_coeff      !< The coefficient correlating SGS temperature variance with the mean
                                 !! temperature gradient in the deterministic part of the Stanley parameterization.
                                 !! Negative values disable the scheme. [nondim]
  logical :: read_khth           !< If true, read a file containing the spatially varying horizontal
                                 !! isopycnal height diffusivity
  logical :: use_stanley_gm      !< If true, also use the Stanley parameterization in MOM_thickness_diffuse

  logical :: use_meso_sfn_ANN  !< If true, use the meso-scale streamfunction ANN parameterization
  type(MESO_SFN_ANN_CS) :: meso_sfn_ANN_CS !< Control structure for the meso-scale streamfunction ANN parameterization

  type(diag_ctrl), pointer :: diag => NULL() !< structure used to regulate timing of diagnostics
  real, allocatable :: GMwork(:,:)        !< Work by isopycnal height diffusion [R Z L2 T-3 ~> W m-2]
  real, allocatable :: diagSlopeX(:,:,:)  !< Diagnostic: zonal neutral slope [Z L-1 ~> nondim]
  real, allocatable :: diagSlopeY(:,:,:)  !< Diagnostic: zonal neutral slope [Z L-1 ~> nondim]

  real, allocatable :: Kh_eta_u(:,:)    !< Isopycnal height diffusivities at u points [L2 T-1 ~> m2 s-1]
  real, allocatable :: Kh_eta_v(:,:)    !< Isopycnal height diffusivities in v points [L2 T-1 ~> m2 s-1]

  real, allocatable :: KH_u_GME(:,:,:)  !< Isopycnal height diffusivities in u-columns [L2 T-1 ~> m2 s-1]
  real, allocatable :: KH_v_GME(:,:,:)  !< Isopycnal height diffusivities in v-columns [L2 T-1 ~> m2 s-1]
  real, allocatable :: khth2d(:,:)      !< 2D isopycnal height diffusivity at h-points [L2 T-1 ~> m2 s-1]

  !>@{
  !! Diagnostic identifier
  integer :: id_uhGM    = -1, id_vhGM    = -1, id_GMwork = -1
  integer :: id_KH_u    = -1, id_KH_v    = -1, id_KH_t   = -1
  integer :: id_KH_u1   = -1, id_KH_v1   = -1, id_KH_t1  = -1
  integer :: id_slope_x = -1, id_slope_y = -1
  integer :: id_sfn_unlim_x = -1, id_sfn_unlim_y = -1, id_sfn_x = -1, id_sfn_y = -1
  !>@}
end type thickness_diffuse_CS


  interface
module subroutine thickness_diffuse(h, uhtr, vhtr, tv, dt, G, GV, US, MEKE, VarMix, CDp, CS, STOCH, u, v)
  type(ocean_grid_type),                      intent(in)    :: G      !< Ocean grid structure
  type(verticalGrid_type),                    intent(in)    :: GV     !< Vertical grid structure
  type(unit_scale_type),                      intent(in)    :: US     !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),  intent(inout) :: h      !< Layer thickness [H ~> m or kg m-2]
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), intent(inout) :: uhtr   !< Accumulated zonal mass flux
                                                                      !! [L2 H ~> m3 or kg]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), intent(inout) :: vhtr   !< Accumulated meridional mass flux
                                                                      !! [L2 H ~> m3 or kg]
  type(thermo_var_ptrs),                      intent(in)    :: tv     !< Thermodynamics structure
  real,                                       intent(in)    :: dt     !< Time increment [T ~> s]
  type(MEKE_type),                            intent(inout) :: MEKE   !< MEKE fields
  type(VarMix_CS), target,                    intent(in)    :: VarMix !< Variable mixing coefficients
  type(cont_diag_ptrs),                       intent(inout) :: CDp    !< Diagnostics for the continuity equation
  type(thickness_diffuse_CS),                 intent(inout) :: CS     !< Control structure for thickness_diffuse
  type(stochastic_CS),                        intent(inout) :: STOCH !< Stochastic control structure
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), intent(in) :: u !< Zonal velocity [L T-1 ~> m s-1].
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), intent(in) :: v !< Meridional velocity [L T-1 ~> m s-1].

  ! Local variables
                                         ! sea level [Z ~> m], positive up.


                  ! weighting of the interface slopes to that calculated also
                  ! using density gradients at u points.  The physically correct
                  ! slopes occur at 0, while 1 is used for numerical closures [nondim].
                  ! weighting of the interface slopes to that calculated also
                  ! using density gradients at v points.  The physically correct
                  ! slopes occur at 0, while 1 is used for numerical closures [nondim].

                    ! in roundoff and can be neglected [H ~> m or kg m-2].
                                    ! to layer centers [L2 T-1 ~> m2 s-1]
                                    ! to layer centers [L2 T-1 ~> m2 s-1]

end subroutine thickness_diffuse
module subroutine thickness_diffuse_full(h, e, Kh_u, Kh_v, tv, uhD, vhD, cg1, dt, G, GV, US, MEKE, &
                                  CS, int_slope_u, int_slope_v, slope_x, slope_y, STOCH, VarMix, &
                                  Sfn_unlim_u_3D, Sfn_unlim_v_3D)
  type(ocean_grid_type),                        intent(in)  :: G     !< Ocean grid structure
  type(verticalGrid_type),                      intent(in)  :: GV    !< Vertical grid structure
  type(unit_scale_type),                        intent(in)  :: US    !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),    intent(in)  :: h     !< Layer thickness [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)+1),  intent(in)  :: e     !< Interface positions [Z ~> m]
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)+1), intent(in)  :: Kh_u  !< Isopycnal height diffusivity
                                                                     !! at u points [L2 T-1 ~> m2 s-1]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)+1), intent(in)  :: Kh_v  !< Isopycnal height diffusivity
                                                                     !! at v points [L2 T-1 ~> m2 s-1]
  type(thermo_var_ptrs),                        intent(in)  :: tv    !< Thermodynamics structure
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)),   intent(out) :: uhD   !< Zonal mass fluxes
                                                                     !! [H L2 T-1 ~> m3 s-1 or kg s-1]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)),   intent(out) :: vhD   !< Meridional mass fluxes
                                                                     !! [H L2 T-1 ~> m3 s-1 or kg s-1]
  real, dimension(:,:),                         pointer     :: cg1   !< Wave speed [L T-1 ~> m s-1]
  real,                                         intent(in)  :: dt    !< Time increment [T ~> s]
  type(MEKE_type),                              intent(inout) :: MEKE !< MEKE fields
  type(thickness_diffuse_CS),                   intent(inout) :: CS  !< Control structure for thickness_diffuse
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)+1), intent(in)  :: int_slope_u !< Ratio that determine how much of
                                                                     !! the isopycnal slopes are taken directly from
                                                                     !! the interface slopes without consideration of
                                                                     !! density gradients [nondim].
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)+1), intent(in)  :: int_slope_v !< Ratio that determine how much of
                                                                     !! the isopycnal slopes are taken directly from
                                                                     !! the interface slopes without consideration of
                                                                     !! density gradients [nondim].
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)+1), optional, intent(in)  :: slope_x !< Isopyc. slope at u [Z L-1 ~> nondim]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)+1), optional, intent(in)  :: slope_y !< Isopyc. slope at v [Z L-1 ~> nondim]
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)+1), optional, intent(in) :: Sfn_unlim_u_3D !< ANN streamfunction
                                                                      !! at u [Z L2 T-1 ~> m3 s-1]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)+1), optional, intent(in) :: Sfn_unlim_v_3D !< ANN streamfunction
                                                                      !! at v [Z L2 T-1 ~> m3 s-1]
  type(stochastic_CS),                       optional, intent(inout)  :: STOCH !< Stochastic control structure
  type(VarMix_CS), target,                      optional, intent(in)  :: VarMix !< Variable mixing coefficents

  ! Local variables
end subroutine thickness_diffuse_full
module subroutine streamfn_solver(nk, c2_h, hN2, sfn)
  integer,               intent(in)    :: nk   !< Number of layers
  real, dimension(nk),   intent(in)    :: c2_h !< Wave speed squared over thickness in layers, rescaled to
                                               !! [H L2 Z-2 T-2 ~> m s-2 or kg m-2 s-2]
  real, dimension(nk+1), intent(in)    :: hN2  !< Thickness times N2 at interfaces times rescaling factors
                                               !! [H L2 Z-2 T-2 ~> m s-2 or kg m-2 s-2]
  real, dimension(nk+1), intent(inout) :: sfn  !< Streamfunction [H L2 T-1 ~> m3 s-1 or kg s-1] or arbitrary units
                                               !! On entry, equals diffusivity times slope.
                                               !! On exit, equals the streamfunction.
  ! Local variables

end subroutine streamfn_solver
module subroutine add_interface_Kh(G, GV, US, CS, Kh_u, Kh_v, Kh_u_CFL, Kh_v_CFL, int_slope_u, int_slope_v)
  type(ocean_grid_type),                        intent(in)    :: G    !< Ocean grid structure
  type(verticalGrid_type),                      intent(in)    :: GV   !< Vertical grid structure
  type(unit_scale_type),                        intent(in)    :: US   !< A dimensional unit scaling type
  type(thickness_diffuse_CS),                   intent(in)    :: CS   !< Control structure for thickness_diffuse
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)+1), intent(inout) :: Kh_u !< Isopycnal height diffusivity
                                                                      !! at u points [L2 T-1 ~> m2 s-1]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)+1), intent(inout) :: Kh_v !< Isopycnal height diffusivity
                                                                      !! at v points [L2 T-1 ~> m2 s-1]
  real, dimension(SZIB_(G),SZJ_(G)),            intent(in)    :: Kh_u_CFL !< Maximum stable isopycnal height
                                                                      !! diffusivity at u points [L2 T-1 ~> m2 s-1]
  real, dimension(SZI_(G),SZJB_(G)),            intent(in)    :: Kh_v_CFL !< Maximum stable isopycnal height
                                                                      !! diffusivity at v points [L2 T-1 ~> m2 s-1]
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)+1), intent(inout) :: int_slope_u !< Ratio that determine how much of
                                                                      !! the isopycnal slopes are taken directly from
                                                                      !! the interface slopes without consideration
                                                                      !! of density gradients [nondim].
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)+1), intent(inout) :: int_slope_v !< Ratio that determine how much of
                                                                      !! the isopycnal slopes are taken directly from
                                                                      !! the interface slopes without consideration
                                                                      !! of density gradients [nondim].

  ! Local variables

end subroutine add_interface_Kh
module subroutine add_detangling_Kh(h, e, Kh_u, Kh_v, KH_u_CFL, KH_v_CFL, tv, dt, G, GV, US, CS, &
                             int_slope_u, int_slope_v)
  type(ocean_grid_type),                        intent(in)    :: G    !< Ocean grid structure
  type(verticalGrid_type),                      intent(in)    :: GV   !< Vertical grid structure
  type(unit_scale_type),                        intent(in)    :: US   !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),    intent(in)    :: h    !< Layer thickness [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)+1),  intent(in)    :: e    !< Interface positions [Z ~> m]
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)+1), intent(inout) :: Kh_u !< Isopycnal height diffusivity
                                                                      !! at u points [L2 T-1 ~> m2 s-1]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)+1), intent(inout) :: Kh_v !< Isopycnal height diffusivity
                                                                      !! at v points [L2 T-1 ~> m2 s-1]
  real, dimension(SZIB_(G),SZJ_(G)),            intent(in)    :: Kh_u_CFL !< Maximum stable isopycnal height
                                                                      !! diffusivity at u points [L2 T-1 ~> m2 s-1]
  real, dimension(SZI_(G),SZJB_(G)),            intent(in)    :: Kh_v_CFL !< Maximum stable isopycnal height
                                                                      !! diffusivity at v points [L2 T-1 ~> m2 s-1]
  type(thermo_var_ptrs),                        intent(in)    :: tv   !< Thermodynamics structure
  real,                                         intent(in)    :: dt   !< Time increment [T ~> s]
  type(thickness_diffuse_CS),                   intent(in)    :: CS   !< Control structure for thickness_diffuse
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)+1), intent(inout) :: int_slope_u !< Ratio that determine how much of
                                                                      !! the isopycnal slopes are taken directly from
                                                                      !! the interface slopes without consideration
                                                                      !! of density gradients [nondim].
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)+1), intent(inout) :: int_slope_v !< Ratio that determine how much of
                                                                      !! the isopycnal slopes are taken directly from
                                                                      !! the interface slopes without consideration
                                                                      !! of density gradients [nondim].
  ! Local variables
               ! region where the detangling is applied [H ~> m or kg m-2].
               ! u points [L2 T-1 ~> m2 s-1].
               ! v points [L2 T-1 ~> m2 s-1].
               ! detangling is applied [H ~> m or kg m-2].
                    ! with the thinner modified near the boundaries to mask out
                    ! thickness variations due to topography, etc.
                    ! from 0 (smooth) to 1 (jagged) [nondim].  This is the difference
                    ! between the arithmetic and harmonic mean thicknesses
                    ! normalized by the arithmetic mean thickness.
                    ! layers [nondim].
                    ! in roundoff and can be neglected [H ~> m or kg m-2].

                    ! above and below [L Z-1 ~> nondim].
                    ! magnitude one [nondim]. 0 <= Rsl <1.
                    ! and the ratio of the face length to the adjacent cell
                    ! areas for comparability with the diffusivities [L Z T-1 ~> m2 s-1].
                    ! the damping timescale [T-1 ~> s-1].
  !   Variables used only in testing code.
  ! real, dimension(SZK_(GV)) :: uh_here ! The transport in a layer [Z L2 T-1 ~> m3 s-1]
  ! real, dimension(SZK_(GV)+1) :: Sfn ! The streamfunction at an interface [Z L T-1 ~> m2 s-1]

end subroutine add_detangling_Kh
module subroutine thickness_diffuse_init(Time, G, GV, US, param_file, diag, CDp, CS)
  type(time_type),         intent(in) :: Time    !< Current model time
  type(ocean_grid_type),   intent(in) :: G       !< Ocean grid structure
  type(verticalGrid_type), intent(in) :: GV      !< Vertical grid structure
  type(unit_scale_type),   intent(in) :: US      !< A dimensional unit scaling type
  type(param_file_type),   intent(in) :: param_file !< Parameter file handles
  type(diag_ctrl), target, intent(inout) :: diag !< Diagnostics control structure
  type(cont_diag_ptrs),    intent(inout) :: CDp  !< Continuity equation diagnostics
  type(thickness_diffuse_CS), intent(inout) :: CS !< Control structure for thickness_diffuse

  ! Local variables
  ! This include declares and sets the variable "version".
                       ! streamfunction formulation, expressed as a fraction of planetary
                       ! rotation divided by an aspect ratio rescaling factor [L Z-1 ~> nondim]
                        ! temperature variance [nondim]
                                 ! as the vertical structure of thickness diffusivity.
                                 ! Used to determine if FULL_DEPTH_KHTH_MIN should be
                                 ! available.

end subroutine thickness_diffuse_init
module subroutine thickness_diffuse_get_KH(CS, KH_u_GME, KH_v_GME, G, GV)
  type(thickness_diffuse_CS),          intent(in)  :: CS   !< Control structure for this module
  type(ocean_grid_type),               intent(in)  :: G    !< Grid structure
  type(verticalGrid_type),             intent(in)  :: GV   !< Vertical grid structure
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)+1), intent(inout) :: KH_u_GME !< Isopycnal height
                                                   !! diffusivities at u-faces [L2 T-1 ~> m2 s-1]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)+1), intent(inout) :: KH_v_GME !< Isopycnal height
                                                   !! diffusivities at v-faces [L2 T-1 ~> m2 s-1]
  ! Local variables

end subroutine thickness_diffuse_get_KH
module subroutine thickness_diffuse_end(CS, CDp)
  type(thickness_diffuse_CS), intent(inout) :: CS !< Control structure for thickness_diffuse
  type(cont_diag_ptrs), intent(inout) :: CDp      !< Continuity diagnostic control structure

end subroutine thickness_diffuse_end
  end interface

end module MOM_thickness_diffuse
