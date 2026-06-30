! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> Implements vertical viscosity (vertvisc)
module MOM_vert_friction

use MOM_domains,       only : pass_var, To_All, Omit_corners
use MOM_domains,       only : pass_vector, Scalar_Pair
use MOM_diag_mediator, only : post_data, register_diag_field, safe_alloc_ptr
use MOM_diag_mediator, only : post_product_u, post_product_sum_u
use MOM_diag_mediator, only : post_product_v, post_product_sum_v
use MOM_diag_mediator, only : diag_ctrl, query_averaging_enabled
use MOM_domains,       only : create_group_pass, do_group_pass, group_pass_type
use MOM_domains,       only : To_North, To_East
use MOM_debugging,     only : uvchksum, hchksum
use MOM_error_handler, only : MOM_error, FATAL, WARNING, NOTE
use MOM_file_parser,   only : get_param, log_param, log_version, param_file_type
use MOM_forcing_type,  only : mech_forcing, find_ustar
use MOM_get_input,     only : directories
use MOM_grid,          only : ocean_grid_type
use MOM_io,            only : MOM_read_data, slasher
use MOM_open_boundary, only : ocean_OBC_type, OBC_NONE, OBC_DIRECTION_E
use MOM_open_boundary, only : OBC_DIRECTION_W, OBC_DIRECTION_N, OBC_DIRECTION_S
use MOM_PointAccel,    only : write_u_accel, write_v_accel, PointAccel_init
use MOM_PointAccel,    only : PointAccel_CS
use MOM_time_manager,  only : time_type, time_minus_signed
use MOM_unit_scaling,  only : unit_scale_type
use MOM_variables,     only : thermo_var_ptrs, vertvisc_type
use MOM_variables,     only : cont_diag_ptrs, accel_diag_ptrs
use MOM_variables,     only : ocean_internal_state
use MOM_verticalGrid,  only : verticalGrid_type
use MOM_wave_interface, only : wave_parameters_CS
use MOM_set_visc,      only : set_v_at_u, set_u_at_v
use MOM_lateral_mixing_coeffs, only : VarMix_CS

use CVMix_kpp,         only : cvmix_kpp_composite_Gshape

implicit none ; private

#include <MOM_memory.h>

public vertvisc, vertvisc_remnant, vertvisc_coef
public vertvisc_limit_vel, vertvisc_init, vertvisc_end
public updateCFLtruncationValue
public vertFPmix

! A note on unit descriptions in comments: MOM6 uses units that can be rescaled for dimensional
! consistency testing. These are noted in comments with units like Z, H, L, and T, along with
! their mks counterparts with notation like "a velocity [Z T-1 ~> m s-1]".  If the units
! vary with the Boussinesq approximation, the Boussinesq variant is given first.

!> The control structure with parameters and memory for the MOM_vert_friction module
type, public :: vertvisc_CS ; private
  logical :: initialized = .false. !< True if this control structure has been initialized.
  real    :: Hmix            !< The mixed layer thickness [Z ~> m].
  real    :: Hmix_stress     !< The mixed layer thickness over which the wind
                             !! stress is applied with direct_stress [H ~> m or kg m-2].
  real    :: Kvml_invZ2      !< The extra vertical viscosity scale in [H Z T-1 ~> m2 s-1 or Pa s] in a
                             !! surface mixed layer with a characteristic thickness given by Hmix,
                             !! and scaling proportional to (Hmix/z)^2, where z is the distance
                             !! from the surface; this can get very large with thin layers.
  real    :: Kv              !< The interior vertical viscosity [H Z T-1 ~> m2 s-1 or Pa s].
  real    :: Hbbl            !< The static bottom boundary layer thickness [Z ~> m].
  real    :: Hbbl_gl90       !< The static bottom boundary layer thickness used for GL90 [Z ~> m].
  real    :: Kv_extra_bbl    !< An extra vertical viscosity in the bottom boundary layer of thickness
                             !! Hbbl when there is not a bottom drag law in use [H Z T-1 ~> m2 s-1 or Pa s].
  real    :: vonKar          !< The von Karman constant as used for mixed layer viscosity [nondim]

  logical :: use_GL90_in_SSW !< If true, use the GL90 parameterization in stacked shallow water mode (SSW).
                             !! The calculation of the GL90 viscosity coefficient uses the fact that in SSW
                             !! we simply have 1/N^2 = h/g^prime, where g^prime is the reduced gravity.
                             !! This identity does not generalize to non-SSW setups.
  logical :: use_GL90_N2     !< If true, use GL90 vertical viscosity coefficient that is depth-independent;
                             !! this corresponds to a kappa_GM that scales as N^2 with depth.
  real    :: kappa_gl90      !< The scalar diffusivity used in the GL90 vertical viscosity scheme
                             !! [L2 H Z-1 T-1 ~> m2 s-1 or Pa s]
  logical :: read_kappa_gl90 !< If true, read a file containing the spatially varying kappa_gl90
  real    :: alpha_gl90      !< Coefficient used to compute a depth-independent GL90 vertical
                             !! viscosity via Kv_gl90 = alpha_gl90 * f^2. Note that the implied
                             !! Kv_gl90 corresponds to a kappa_gl90 that scales as N^2 with depth.
                             !! [H Z T ~> m2 s or kg s m-1]
  real    :: vel_underflow   !< Velocity components smaller than vel_underflow
                             !! are set to 0 [L T-1 ~> m s-1].
  real    :: CFL_trunc       !< Velocity components will be truncated when they
                             !! are large enough that the corresponding CFL number
                             !! exceeds this value [nondim].
  real    :: CFL_report      !< The value of the CFL number that will cause the
                             !! accelerations to be reported [nondim].  CFL_report
                             !! will often equal CFL_trunc.
  real    :: truncRampTime   !< The time-scale over which to ramp up the value of
                             !! CFL_trunc from CFL_truncS to CFL_truncE [T ~> s]
  real    :: CFL_truncS      !< The start value of CFL_trunc [nondim]
  real    :: CFL_truncE      !< The end/target value of CFL_trunc [nondim]
  logical :: CFLrampingIsActivated = .false. !< True if the ramping has been initialized
  type(time_type) :: rampStartTime !< The time at which the ramping of CFL_trunc starts

  real ALLOCABLE_, dimension(NIMEMB_PTR_,NJMEM_,NK_INTERFACE_) :: &
    a_u                !< The u-drag coefficient across an interface [H T-1 ~> m s-1 or Pa s m-1]
  real ALLOCABLE_, dimension(NIMEMB_PTR_,NJMEM_,NK_INTERFACE_) :: &
    a_u_gl90           !< The u-drag coefficient associated with GL90 across an interface [H T-1 ~> m s-1 or Pa s m-1]
  real ALLOCABLE_, dimension(NIMEMB_PTR_,NJMEM_,NKMEM_) :: &
    h_u                !< The effective layer thickness at u-points [H ~> m or kg m-2].
  real ALLOCABLE_, dimension(NIMEM_,NJMEMB_PTR_,NK_INTERFACE_) :: &
    a_v                !< The v-drag coefficient across an interface [H T-1 ~> m s-1 or Pa s m-1]
  real ALLOCABLE_, dimension(NIMEM_,NJMEMB_PTR_,NK_INTERFACE_) :: &
    a_v_gl90           !< The v-drag coefficient associated with GL90 across an interface [H T-1 ~> m s-1 or Pa s m-1]
  real ALLOCABLE_, dimension(NIMEM_,NJMEMB_PTR_,NKMEM_) :: &
    h_v                !< The effective layer thickness at v-points [H ~> m or kg m-2].
  real, pointer, dimension(:,:) :: a1_shelf_u => NULL() !< The u-momentum coupling coefficient under
                           !! ice shelves [H T-1 ~> m s-1 or Pa s m-1]. Retained to determine stress under shelves.
  real, pointer, dimension(:,:) :: a1_shelf_v => NULL() !< The v-momentum coupling coefficient under
                           !! ice shelves [H T-1 ~> m s-1 or Pa s m-1]. Retained to determine stress under shelves.

  logical :: split          !< If true, use the split time stepping scheme.
  logical :: bottomdraglaw  !< If true, the  bottom stress is calculated with a
                            !! drag law c_drag*|u|*u. The velocity magnitude
                            !! may be an assumed value or it may be based on the
                            !! actual velocity in the bottommost HBBL, depending
                            !! on whether linear_drag is true.
  logical :: harmonic_visc  !< If true, the harmonic mean thicknesses are used
                            !! to calculate the viscous coupling between layers
                            !! except near the bottom.  Otherwise the arithmetic
                            !! mean thickness is used except near the bottom.
  real    :: harm_BL_val    !< A scale to determine when water is in the boundary
                            !! layers based solely on harmonic mean thicknesses
                            !! for the purpose of determining the extent to which
                            !! the thicknesses used in the viscosities are upwinded [nondim].
  logical :: direct_stress  !< If true, the wind stress is distributed over the topmost Hmix_stress
                            !! of fluid, and an added mixed layer viscosity or a physically based
                            !! boundary layer turbulence parameterization is not needed for stability.
  logical :: dynamic_viscous_ML  !< If true, use the results from a dynamic
                            !! calculation, perhaps based on a bulk Richardson
                            !! number criterion, to determine the mixed layer
                            !! thickness for viscosity.
  logical :: fixed_LOTW_ML  !< If true, use a Law-of-the-wall prescription for the mixed layer
                            !! viscosity within a boundary layer that is the lesser of Hmix and the
                            !! total depth of the ocean in a column.
  logical :: apply_LOTW_floor !< If true, use a Law-of-the-wall prescription to set a lower bound
                            !! on the viscous coupling between layers within the surface boundary
                            !! layer, based the distance of interfaces from the surface.  This only
                            !! acts when there are large changes in the thicknesses of successive
                            !! layers or when the viscosity is set externally and the wind stress
                            !! has subsequently increased.
  integer :: answer_date    !< The vintage of the order of arithmetic and expressions in the viscous
                            !! calculations.  Values below 20190101 recover the answers from the end
                            !! of 2018, while higher values use expressions that do not use an
                            !! arbitrary and hard-coded maximum viscous coupling coefficient between
                            !! layers.  In non-Boussinesq cases, values below 20230601 recover a
                            !! form of the viscosity within  the mixed layer that breaks up the
                            !! magnitude of the wind stress with BULKMIXEDLAYER, DYNAMIC_VISCOUS_ML
                            !! or FIXED_DEPTH_LOTW_ML, but not LOTW_VISCOUS_ML_FLOOR.
  logical :: debug          !< If true, write verbose checksums for debugging purposes.
  integer :: nkml           !< The number of layers in the mixed layer.
  integer, pointer :: ntrunc !< The number of times the velocity has been
                            !! truncated since the last call to write_energy.
  character(len=200) :: u_trunc_file  !< The complete path to a file in which a column of
                            !! u-accelerations are written if velocity truncations occur.
  character(len=200) :: v_trunc_file !< The complete path to a file in which a column of
                            !! v-accelerations are written if velocity truncations occur.
  logical :: StokesMixing   !< If true, do Stokes drift mixing via the Lagrangian current
                            !! (Eulerian plus Stokes drift).  False by default and set
                            !! via STOKES_MIXING_COMBINED.

  type(diag_ctrl), pointer :: diag !< A structure that is used to regulate the
                                   !! timing of diagnostic output.
  real, allocatable, dimension(:,:) :: kappa_gl90_2d !< 2D kappa_gl90 at h-points [L2 H Z-1 T-1 ~> m2 s-1 or Pa s]

  !>@{ Diagnostic identifiers
  integer :: id_du_dt_visc = -1, id_dv_dt_visc = -1, id_du_dt_visc_gl90 = -1, id_dv_dt_visc_gl90 = -1
  integer :: id_GLwork = -1
  integer :: id_au_vv = -1, id_av_vv = -1, id_au_gl90_vv = -1, id_av_gl90_vv = -1
  integer :: id_du_dt_str = -1, id_dv_dt_str = -1
  integer :: id_h_u = -1, id_h_v = -1, id_hML_u = -1 , id_hML_v = -1
  integer :: id_Omega_w2x = -1, id_FPtau2s  = -1 , id_FPtau2w = -1
  integer :: id_uE_h  = -1, id_vE_h  = -1
  integer :: id_uStk  = -1, id_vStk  = -1
  integer :: id_uStk0 = -1, id_vStk0 = -1
  integer :: id_uInc_h= -1, id_vInc_h= -1
  integer :: id_taux_bot = -1, id_tauy_bot = -1
  integer :: id_Kv_slow = -1, id_Kv_u = -1, id_Kv_v = -1
  integer :: id_Kv_gl90_u = -1, id_Kv_gl90_v = -1
  ! integer :: id_hf_du_dt_visc    = -1, id_hf_dv_dt_visc    = -1
  integer :: id_h_du_dt_visc    = -1, id_h_dv_dt_visc    = -1
  integer :: id_hf_du_dt_visc_2d = -1, id_hf_dv_dt_visc_2d = -1
  integer :: id_h_du_dt_str    = -1, id_h_dv_dt_str    = -1
  integer :: id_du_dt_str_visc_rem = -1, id_dv_dt_str_visc_rem = -1
  !>@}

  type(PointAccel_CS), pointer :: PointAccel_CSp => NULL() !< A pointer to the control structure
                              !! for recording accelerations leading to velocity truncations

  type(group_pass_type) :: pass_KE_uv !< A handle used for group halo passes
end type vertvisc_CS


  interface
module subroutine vertFPmix(ui, vi, uold, vold, hbl_h, h, forces, dt, lpost, Cemp_NL, G, GV, US, CS, OBC, Waves)
  type(ocean_grid_type),   intent(in)    :: G      !< Ocean grid structure
  type(verticalGrid_type), intent(in)    :: GV     !< Ocean vertical grid structure
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                           intent(inout) :: ui     !< Zonal velocity after vertvisc [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                           intent(inout) :: vi     !< Meridional velocity after vertvisc [L T-1 ~> m s-1]
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                           intent(inout) :: uold   !< Old Zonal velocity [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                           intent(inout) :: vold   !< Old Meridional velocity [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJ_(G)), intent(inout) :: hbl_h !<  boundary layer depth [H ~> m]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in) :: h       !< Layer thicknesses [H ~> m or kg m-2]
  type(mech_forcing),      intent(in) :: forces  !< A structure with the driving mechanical forces
  real,                    intent(in) :: dt      !< Time increment [T ~> s]
  real,                    intent(in) :: Cemp_NL !< empirical coefficient of non-local momentum mixing [nondim]
  logical,                 intent(in) :: lpost   !< Compute and make available FPMix diagnostics
  type(unit_scale_type),   intent(in) :: US      !< A dimensional unit scaling type
  type(vertvisc_CS),       pointer    :: CS      !< Vertical viscosity control structure
  type(ocean_OBC_type),    pointer    :: OBC     !< Open boundary condition structure
  type(wave_parameters_CS), &
                   optional, pointer  :: Waves   !< Container for wave/Stokes information

  ! local variables

end subroutine vertFPmix
module subroutine find_coupling_coef_gl90(a_cpl_gl90, hvel, i, j, z_i, G, GV, CS, VarMix, work_on_u)
  type(ocean_grid_type), intent(in) :: G        !< Grid structure.
  type(verticalGrid_type), intent(in) :: GV     !< Vertical grid structure.
  real, dimension(SZK_(GV)), intent(in) :: hvel !< Distance between interfaces
                                                !! at velocity points [Z ~> m]
  integer, intent(in) :: i                      !< Column i-index
  integer, intent(in) :: j                      !< Column j-index
  real, dimension(SZK_(GV)+1), intent(in) :: z_i  !< Estimate of interface heights above the
                                                !! bottom, normalized by the GL90 bottom
                                                !! boundary layer thickness [nondim]
  real, dimension(SZK_(GV)+1),intent(out) :: a_cpl_gl90   !< Coupling coefficient associated
                                                !! with GL90 across interfaces; is not
                                                !! included in a_cpl [H T-1 ~> m s-1 or Pa s m-1].
  type(vertvisc_cs), intent(in) :: CS           !< Vertical viscosity control structure
  type(VarMix_CS), intent(in) :: VarMix         !< Variable mixing coefficients
  logical, intent(in) :: work_on_u              !< If true, u-points are being calculated,
                                                !! otherwise they are v-points.

  ! local variables
                        ! and can be neglected [Z ~> m].

end subroutine find_coupling_coef_gl90
module subroutine vertvisc(u, v, h, forces, visc, dt, OBC, ADp, CDp, G, GV, US, CS, &
                    taux_bot, tauy_bot, fpmix, Waves)
  type(ocean_grid_type),   intent(in)    :: G      !< Ocean grid structure
  type(verticalGrid_type), intent(in)    :: GV     !< Ocean vertical grid structure
  type(unit_scale_type),   intent(in)    :: US     !< A dimensional unit scaling type
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                           intent(inout) :: u      !< Zonal velocity [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                           intent(inout) :: v      !< Meridional velocity [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in)    :: h      !< Layer thickness [H ~> m or kg m-2]
  type(mech_forcing),    intent(in)      :: forces !< A structure with the driving mechanical forces
  type(vertvisc_type),   intent(inout)   :: visc   !< Viscosities and bottom drag
  real,                  intent(in)      :: dt     !< Time increment [T ~> s]
  type(ocean_OBC_type),  pointer         :: OBC    !< Open boundary condition structure
  type(accel_diag_ptrs), intent(inout)   :: ADp    !< Accelerations in the momentum
                                                   !! equations for diagnostics
  type(cont_diag_ptrs),  intent(inout)   :: CDp    !< Continuity equation terms
  type(vertvisc_CS),     pointer         :: CS     !< Vertical viscosity control structure
  real, dimension(SZIB_(G),SZJ_(G)), &
                   optional, intent(out) :: taux_bot !< Zonal bottom stress from ocean to
                                                     !! rock [R L Z T-2 ~> Pa]
  real, dimension(SZI_(G),SZJB_(G)), &
                   optional, intent(out) :: tauy_bot !< Meridional bottom stress from ocean to
                                                     !! rock [R L Z T-2 ~> Pa]
  logical,         optional, intent(in)  :: fpmix !< fpmix along Eulerian shear
  type(wave_parameters_CS), &
                   optional, pointer     :: Waves !< Container for wave/Stokes information

  ! Fields from forces used in this subroutine:
  !   taux: Zonal wind stress [R L Z T-2 ~> Pa].
  !   tauy: Meridional wind stress [R L Z T-2 ~> Pa].

  ! Local variables

    ! A variable used by the tridiagonal solver [H-1 ~> m-1 or m2 kg-1].
    ! A variable used by the tridiagonal solver [nondim].
    ! d1=1-c1 is used by the tridiagonal solver [nondim].
    ! Ray is the Rayleigh-drag velocity [H T-1 ~> m s-1 or Pa s m-1]
    ! The first term in the denominator of b1 [H ~> m or kg m-2].

                           ! is applied with direct_stress [H ~> m or kg m-2].
                           ! in roundoff and can be neglected [H ~> m or kg m-2].

                           ! by the density [H L T-1 ~> m2 s-1 or kg m-1 s-1].
                           ! than this are diagnosed as 0 [L T-2 ~> m s-2].
    ! The same as stress, unless the wind stress is applied as a body force
    ! [H L T-1 ~> m2 s-1 or kg m-1 s-1].
                                                 ! [H L2 T-3 ~> m3 s-3 or W m-2]
                                              ! [H L4 T-3 ~> m5 s-3 or kg m2 s-3]
                                              ! [H L4 T-3 ~> m5 s-3 or kg m2 s-3]


end subroutine vertvisc
module subroutine vertvisc_remnant(visc, visc_rem_u, visc_rem_v, dt, G, GV, US, CS)
  type(ocean_grid_type), intent(in)   :: G    !< Ocean grid structure
  type(verticalGrid_type), intent(in) :: GV   !< Ocean vertical grid structure
  type(vertvisc_type),   intent(in)   :: visc !< Viscosities and bottom drag
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                         intent(inout) :: visc_rem_u !< Fraction of a time-step's worth of a
                                              !! barotropic acceleration that a layer experiences after
                                              !! viscosity is applied in the zonal direction [nondim]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                         intent(inout) :: visc_rem_v !< Fraction of a time-step's worth of a
                                              !! barotropic acceleration that a layer experiences after
                                              !! viscosity is applied in the meridional direction [nondim]
  real,                  intent(in)    :: dt  !< Time increment [T ~> s]
  type(unit_scale_type), intent(in)    :: US  !< A dimensional unit scaling type
  type(vertvisc_CS),     pointer       :: CS  !< Vertical viscosity control structure

  ! Local variables

    ! A variable used by the tridiagonal solver [H-1 ~> m-1 or m2 kg-1].
    ! A variable used by the tridiagonal solver [nondim].
    ! d1=1-c1 is used by the tridiagonal solver [nondim].
    ! Ray is the Rayleigh-drag velocity [H T-1 ~> m s-1 or Pa s m-1]
    ! The first term in the denominator of b1 [H ~> m or kg m-2].

end subroutine vertvisc_remnant
module subroutine vertvisc_coef(u, v, h, dz, forces, visc, tv, dt, G, GV, US, CS, OBC, VarMix)
  type(ocean_grid_type),   intent(in)    :: G      !< Ocean grid structure
  type(verticalGrid_type), intent(in)    :: GV     !< Ocean vertical grid structure
  type(unit_scale_type),   intent(in)    :: US     !< A dimensional unit scaling type
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                           intent(in)    :: u      !< Zonal velocity [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                           intent(in)    :: v      !< Meridional velocity [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in)    :: h      !< Layer thickness [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in)    :: dz     !< Vertical distance across layers [Z ~> m]
  type(mech_forcing),      intent(in)    :: forces !< A structure with the driving mechanical forces
  type(vertvisc_type),     intent(in)    :: visc   !< Viscosities and bottom drag
  type(thermo_var_ptrs),   intent(in)    :: tv     !< A structure containing pointers to any available
                                                   !! thermodynamic fields.
  real,                    intent(in)    :: dt     !< Time increment [T ~> s]
  type(vertvisc_CS),       intent(inout) :: CS     !< Vertical viscosity control structure
  type(ocean_OBC_type),    pointer       :: OBC    !< Open boundary condition structure
  type(VarMix_CS),         intent(in) :: VarMix !< Variable mixing coefficients
  ! Field from forces used in this subroutine:
  !   ustar: the friction velocity [Z T-1 ~> m s-1], used here as the mixing
  !     velocity in the mixed layer if NKML > 1 in a bulk mixed layer.

  ! Local variables

end subroutine vertvisc_coef
module subroutine find_coupling_coef(a_cpl, hvel, i, j, h_harm, bbl_thick, kv_bbl, z_i, h_ml, &
                              dt, G, GV, US, CS, visc, Ustar_2d, tv, work_on_u, OBC, shelf)
  type(ocean_grid_type),     intent(in)  :: G  !< Ocean grid structure
  type(verticalGrid_type),   intent(in)  :: GV !< Ocean vertical grid structure
  type(unit_scale_type),     intent(in)  :: US !< A dimensional unit scaling type
  real, dimension(SZK_(GV)+1), &
                             intent(out) :: a_cpl !< Coupling coefficient across interfaces [H T-1 ~> m s-1 or Pa s m-1]
  real, dimension(SZK_(GV)), &
                             intent(in)  :: hvel !< Distance between interfaces at velocity points [Z ~> m]
  integer,                   intent(in)  :: i    !< Column i-index
  integer,                   intent(in)  :: j    !< Column j-index
  real, dimension(SZK_(GV)), &
                             intent(in)  :: h_harm !< Harmonic mean of thicknesses around a velocity
                                                   !! grid point [Z ~> m]
  real, intent(in)  :: bbl_thick !< Bottom boundary layer thickness [Z ~> m]
  real, intent(in)  :: kv_bbl !< Bottom boundary layer viscosity, exclusive of
                                                   !! any depth-dependent contributions from
                                                   !! visc%Kv_shear [H Z T-1 ~> m2 s-1 or Pa s]
  real, dimension(SZK_(GV)+1), &
                             intent(in)  :: z_i  !< Estimate of interface heights above the bottom,
                                                 !! normalized by the bottom boundary layer thickness [nondim]
  real, intent(out) :: h_ml !< Mixed layer depth [Z ~> m]
  real,                      intent(in)  :: dt   !< Time increment [T ~> s]
  type(vertvisc_CS),         intent(in)  :: CS   !< Vertical viscosity control structure
  type(vertvisc_type),       intent(in)  :: visc !< Structure containing viscosities and bottom drag
  real, dimension(SZI_(G),SZJ_(G)), &
                             intent(in)  :: Ustar_2d !< The wind friction velocity, calculated using
                                                 !! the Boussinesq reference density or the
                                                 !! time-evolving surface density in non-Boussinesq
                                                 !! mode [Z T-1 ~> m s-1]
  type(thermo_var_ptrs),     intent(in)  :: tv   !< A structure containing pointers to any available
                                                 !! thermodynamic fields.
  logical,                   intent(in)  :: work_on_u !< If true, u-points are being calculated,
                                                  !! otherwise they are v-points
  type(ocean_OBC_type),      pointer     :: OBC   !< Open boundary condition structure
  logical,         optional, intent(in)  :: shelf !< If present and true, use a surface boundary
                                                  !! condition appropriate for an ice shelf.

  ! Local variables

end subroutine find_coupling_coef
module subroutine vertvisc_limit_vel(u, v, h, ADp, CDp, forces, visc, dt, G, GV, US, CS)
  type(ocean_grid_type),   intent(in)    :: G      !< Ocean grid structure
  type(verticalGrid_type), intent(in)    :: GV     !< Ocean vertical grid structure
  type(unit_scale_type),   intent(in)    :: US     !< A dimensional unit scaling type
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                           intent(inout) :: u      !< Zonal velocity [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                           intent(inout) :: v      !< Meridional velocity [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in)    :: h      !< Layer thickness [H ~> m or kg m-2]
  type(accel_diag_ptrs),   intent(in)    :: ADp    !< Acceleration diagnostic pointers
  type(cont_diag_ptrs),    intent(in)    :: CDp    !< Continuity diagnostic pointers
  type(mech_forcing),      intent(in)    :: forces !< A structure with the driving mechanical forces
  type(vertvisc_type),     intent(in)    :: visc   !< Viscosities and bottom drag
  real,                    intent(in)    :: dt     !< Time increment [T ~> s]
  type(vertvisc_CS),       pointer       :: CS     !< Vertical viscosity control structure

  ! Local variables
end subroutine vertvisc_limit_vel
module subroutine vertvisc_init(MIS, Time, G, GV, US, param_file, diag, ADp, dirs, &
                          ntrunc, CS, fpmix)
  type(ocean_internal_state), &
                   target, intent(in)    :: MIS    !< The "MOM Internal State", a set of pointers
                                                   !! to the fields and accelerations that make
                                                   !! up the ocean's physical state
  type(time_type), target, intent(in)    :: Time   !< Current model time
  type(ocean_grid_type),   intent(in)    :: G      !< Ocean grid structure
  type(verticalGrid_type), intent(in)    :: GV     !< Ocean vertical grid structure
  type(unit_scale_type),   intent(in)    :: US     !< A dimensional unit scaling type
  type(param_file_type),   intent(in)    :: param_file !< File to parse for parameters
  type(diag_ctrl), target, intent(inout) :: diag   !< Diagnostic control structure
  type(accel_diag_ptrs),   intent(inout) :: ADp    !< Acceleration diagnostic pointers
  type(directories),       intent(in)    :: dirs   !< Relevant directory paths
  integer, target,         intent(inout) :: ntrunc !< Number of velocity truncations
  type(vertvisc_CS),       pointer       :: CS     !< Vertical viscosity control structure
  logical, optional,       intent(in)    :: fpmix  !< Nonlocal momentum mixing

  ! Local variables

  ! This include declares and sets the variable "version".

end subroutine vertvisc_init
module subroutine updateCFLtruncationValue(Time, CS, US, activate)
  type(time_type), target, intent(in)    :: Time     !< Current model time
  type(vertvisc_CS),       pointer       :: CS       !< Vertical viscosity control structure
  type(unit_scale_type),   intent(in)    :: US       !< A dimensional unit scaling type
  logical, optional,       intent(in)    :: activate !< Specify whether to record the value of
                                                     !! Time as the beginning of the ramp period

  ! Local variables

end subroutine updateCFLtruncationValue
module subroutine vertvisc_end(CS)
  type(vertvisc_CS), intent(inout) :: CS  !< Vertical viscosity control structure that
                                          !! will be deallocated in this subroutine.

end subroutine vertvisc_end
  end interface

end module MOM_vert_friction
