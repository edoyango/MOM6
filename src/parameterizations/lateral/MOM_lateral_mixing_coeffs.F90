! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> Variable mixing coefficients
module MOM_lateral_mixing_coeffs

use MOM_debugging,         only : hchksum, uvchksum
use MOM_error_handler,     only : MOM_error, FATAL, WARNING, MOM_mesg
use MOM_diag_mediator,     only : register_diag_field, safe_alloc_ptr, post_data
use MOM_diag_mediator,     only : diag_ctrl, time_type, query_averaging_enabled
use MOM_domains,           only : create_group_pass, do_group_pass
use MOM_domains,           only : group_pass_type, pass_var, pass_vector
use MOM_EOS,               only : calculate_density_derivs, EOS_domain
use MOM_file_parser,       only : get_param, log_version, param_file_type
use MOM_interface_heights, only : find_eta, thickness_to_dz
use MOM_isopycnal_slopes,  only : calc_isoneutral_slopes
use MOM_grid,              only : ocean_grid_type
use MOM_unit_scaling,      only : unit_scale_type
use MOM_variables,         only : thermo_var_ptrs
use MOM_verticalGrid,      only : verticalGrid_type
use MOM_wave_speed,        only : wave_speed, wave_speed_CS, wave_speed_init
use MOM_open_boundary,     only : ocean_OBC_type, OBC_NONE
use MOM_open_boundary,     only : OBC_DIRECTION_E, OBC_DIRECTION_W, OBC_DIRECTION_N, OBC_DIRECTION_S
use MOM_MEKE_types,        only : MEKE_type

implicit none ; private

#include <MOM_memory.h>

!> Variable mixing coefficients
type, public :: VarMix_CS
  logical :: initialized = .false. !< True if this control structure has been initialized.
  logical :: use_variable_mixing  !< If true, use the variable mixing.
  logical :: Resoln_scaling_used  !< If true, a resolution function is used somewhere to scale
                                  !! away one of the viscosities or diffusivities when the
                                  !! deformation radius is well resolved.
  logical :: Resoln_scaled_Kh     !< If true, scale away the Laplacian viscosity
                                  !! when the deformation radius is well resolved.
  logical :: Resoln_scaled_KhTh   !< If true, scale away the thickness diffusivity
                                  !! when the deformation radius is well resolved.
  logical :: Depth_scaled_KhTh    !< If true, KHTH is scaled away when the depth is
                                  !! shallower than a reference depth.
  logical :: Resoln_scaled_KhTr   !< If true, scale away the tracer diffusivity
                                  !! when the deformation radius is well resolved.
  logical :: interpolate_Res_fn   !< If true, interpolate the resolution function
                                  !! to the velocity points from the thickness
                                  !! points; otherwise interpolate the wave
                                  !! speed and calculate the resolution function
                                  !! independently at each point.
  logical :: use_stored_slopes    !< If true, stores isopycnal slopes in this structure.
  logical :: Resoln_use_ebt       !< If true, use the equivalent barotropic wave speed instead of the
                                  !! first baroclinic wave speed for calculating the resolution function.
  logical :: khth_use_ebt_struct  !< If true, uses the equivalent barotropic structure
                                  !! as the vertical structure of thickness diffusivity.
  logical :: kdgl90_use_ebt_struct  !< If true, uses the equivalent barotropic structure
                                  !! as the vertical structure of diffusivity in the GL90 scheme.
  logical :: kdgl90_use_sqg_struct  !< If true, uses the surface quasigeostrophic structure
                                  !! as the vertical structure of diffusivity in the GL90 scheme.
  logical :: khth_use_sqg_struct  !< If true, uses the surface quasigeostrophic structure
                                  !! as the vertical structure of thickness diffusivity.
  logical :: khtr_use_ebt_struct  !< If true, uses the equivalent barotropic structure
                                  !! as the vertical structure of tracer diffusivity.
  logical :: khtr_use_sqg_struct  !< If true, uses the surface quasigeostrophic structure
                                  !! as the vertical structure of tracer diffusivity.
  logical :: calculate_cg1        !< If true, calls wave_speed() to calculate the first
                                  !! baroclinic wave speed and populate CS%cg1.
                                  !! This parameter is set depending on other parameters.
  logical :: calculate_Rd_dx      !< If true, calculates Rd/dx and populate CS%Rd_dx_h.
                                  !! This parameter is set depending on other parameters.
  logical :: calculate_res_fns    !< If true, calculate all the resolution factors.
                                  !! This parameter is set depending on other parameters.
  logical :: calculate_depth_fns !< If true, calculate all the depth factors.
                                  !! This parameter is set depending on other parameters.
  logical :: calculate_Eady_growth_rate !< If true, calculate all the Eady growth rates.
                                  !! This parameter is set depending on other parameters.
  logical :: use_stanley_iso      !< If true, use Stanley parameterization in MOM_isopycnal_slopes
  logical :: use_simpler_Eady_growth_rate !< If true, use a simpler method to calculate the
                                  !! Eady growth rate that avoids division by layer thickness.
                                  !! This parameter is set depending on other parameters.
  logical :: full_depth_Eady_growth_rate !< If true, calculate the Eady growth rate based on an
                                  !! average that includes contributions from sea-level changes
                                  !! in its denominator, rather than just the nominal depth of
                                  !! the bathymetry.  This only applies when using the model
                                  !! interface heights as a proxy for isopycnal slopes.
  logical :: OBC_friendly         !< If true, use only interior data for thickness weighting and
                                  !! to calculate stratification and other fields at open boundary
                                  !! condition faces.
  logical :: res_fn_OBC_bug       !< If false, use only interior data for calculating the resolution
                                  !! functions at open boundary condition faces and vertices.
  real :: cropping_distance       !< Distance from surface or bottom to filter out outcropped or
                                  !! incropped interfaces for the Eady growth rate calc [Z ~> m]
  real :: h_min_N2                !< The minimum vertical distance to use in the denominator of the
                                  !! buoyancy frequency used in the slope calculation [H ~> m or kg m-2]

  real, allocatable :: SN_u(:,:)      !< S*N at u-points [T-1 ~> s-1]
  real, allocatable :: SN_v(:,:)      !< S*N at v-points [T-1 ~> s-1]
  real, allocatable :: L2u(:,:)       !< Length scale^2 at u-points [L2 ~> m2]
  real, allocatable :: L2v(:,:)       !< Length scale^2 at v-points [L2 ~> m2]
  real, allocatable :: cg1(:,:)       !< The first baroclinic gravity wave speed [L T-1 ~> m s-1].
  real, allocatable :: Res_fn_h(:,:)  !< Non-dimensional function of the ratio the first baroclinic
                                      !! deformation radius to the grid spacing at h points [nondim].
  real, allocatable :: Res_fn_q(:,:)  !< Non-dimensional function of the ratio the first baroclinic
                                      !! deformation radius to the grid spacing at q points [nondim].
  real, allocatable :: Res_fn_u(:,:)  !< Non-dimensional function of the ratio the first baroclinic
                                      !! deformation radius to the grid spacing at u points [nondim].
  real, allocatable :: Res_fn_v(:,:)  !< Non-dimensional function of the ratio the first baroclinic
                                      !! deformation radius to the grid spacing at v points [nondim].
  real, allocatable :: Depth_fn_u(:,:) !< Non-dimensional function of the ratio of the depth to
                                      !! a reference depth (maximum 1) at u points [nondim]
  real, allocatable :: Depth_fn_v(:,:) !< Non-dimensional function of the ratio of the depth to
                                      !! a reference depth (maximum 1) at v points [nondim]
  real, allocatable :: beta_dx2_h(:,:) !< The magnitude of the gradient of the Coriolis parameter
                                      !! times the grid spacing squared at h points [L T-1 ~> m s-1].
  real, allocatable :: beta_dx2_q(:,:) !< The magnitude of the gradient of the Coriolis parameter
                                      !! times the grid spacing squared at q points [L T-1 ~> m s-1].
  real, allocatable :: beta_dx2_u(:,:) !< The magnitude of the gradient of the Coriolis parameter
                                      !! times the grid spacing squared at u points [L T-1 ~> m s-1].
  real, allocatable :: beta_dx2_v(:,:) !< The magnitude of the gradient of the Coriolis parameter
                                      !! times the grid spacing squared at v points [L T-1 ~> m s-1].
  real, allocatable :: f2_dx2_h(:,:)  !< The Coriolis parameter squared times the grid
                                      !! spacing squared at h [L2 T-2 ~> m2 s-2].
  real, allocatable :: f2_dx2_q(:,:)  !< The Coriolis parameter squared times the grid
                                      !! spacing squared at q [L2 T-2 ~> m2 s-2].
  real, allocatable :: f2_dx2_u(:,:)  !< The Coriolis parameter squared times the grid
                                      !! spacing squared at u [L2 T-2 ~> m2 s-2].
  real, allocatable :: f2_dx2_v(:,:)  !< The Coriolis parameter squared times the grid
                                      !! spacing squared at v [L2 T-2 ~> m2 s-2].
  real, allocatable :: Rd_dx_h(:,:)   !< Deformation radius over grid spacing [nondim]

  real, allocatable :: slope_x(:,:,:)     !< Zonal isopycnal slope [Z L-1 ~> nondim]
  real, allocatable :: slope_y(:,:,:)     !< Meridional isopycnal slope [Z L-1 ~> nondim]
  real, allocatable :: ebt_struct(:,:,:)  !< EBT vertical structure to scale diffusivities with [nondim]
  real, allocatable :: sqg_struct(:,:,:)  !< SQG vertical structure to scale diffusivities with [nondim]
  real, allocatable :: BS_struct(:,:,:) !< Vertical structure function used in backscatter [nondim]
  real, allocatable :: khth_struct(:,:,:) !< Vertical structure function used in thickness diffusivity [nondim]
  real, allocatable :: khtr_struct(:,:,:) !< Vertical structure function used in tracer diffusivity [nondim]
  real, allocatable :: kdgl90_struct(:,:,:) !< Vertical structure function used in GL90 diffusivity [nondim]
  real :: BS_EBT_power                !< Power to raise EBT vertical structure to. Default 0.0.
  real :: sqg_expo     !< Exponent for SQG vertical structure [nondim]. Default 1.0
  logical :: interpolated_sqg_struct  !< If true, interpolate properties to velocity points and then
                                      !! interpolate the buoyancy frequencies and layer thicknesses
                                      !! back to tracer points when calculating the SQG vertical
                                      !! structure.
  logical :: BS_use_sqg_struct   !< If true, use sqg_stuct for backscatter vertical structure.

  real, allocatable :: Laplac3_const_u(:,:) !< Laplacian metric-dependent constants at u-points [L3 ~> m3]
  real, allocatable :: Laplac3_const_v(:,:) !< Laplacian metric-dependent constants at u-points [L3 ~> m3]
  real, allocatable :: KH_u_QG(:,:,:) !< QG Leith GM coefficient at u-points [L2 T-1 ~> m2 s-1]
  real, allocatable :: KH_v_QG(:,:,:) !< QG Leith GM coefficient at v-points [L2 T-1 ~> m2 s-1]

  ! Parameters
  logical :: use_Visbeck  !< Use Visbeck formulation for thickness diffusivity
  integer :: VarMix_Ktop  !< Top layer to start downward integrals
  real :: Visbeck_L_scale !< Fixed length scale in Visbeck formula [L ~> m], or if negative a scaling
                          !! factor [nondim] relating this length scale squared to the cell area
  real :: Eady_GR_D_scale !< Depth over which to average SN [Z ~> m]
  real :: Res_coef_khth   !< A coefficient [nondim] that determines the function
                          !! of resolution, used for thickness and tracer mixing, as:
                          !!  F = 1 / (1 + (Res_coef_khth*Ld/dx)^Res_fn_power)
  real :: Res_coef_visc   !< A coefficient [nondim] that determines the function
                          !! of resolution, used for lateral viscosity, as:
                          !!  F = 1 / (1 + (Res_coef_visc*Ld/dx)^Res_fn_power)
  real :: depth_scaled_khth_h0 !< The depth above which KHTH is linearly scaled away [Z ~> m]
  real :: depth_scaled_khth_exp !< The exponent used in the depth dependent scaling function for KHTH [nondim]
  real :: kappa_smooth    !< A diffusivity for smoothing T/S in vanished layers [H Z T-1 ~> m2 s-1 or kg m-1 s-1]
  integer :: Res_fn_power_khth !< The power of dx/Ld in the KhTh resolution function.  Any
                               !! positive integer power may be used, but even powers
                               !! and especially 2 are coded to be more efficient.
  integer :: Res_fn_power_visc !< The power of dx/Ld in the Kh resolution function.  Any
                               !! positive integer power may be used, but even powers
                               !! and especially 2 are coded to be more efficient.
  real :: Visbeck_S_max   !< Upper bound on slope used in Eady growth rate [Z L-1 ~> nondim].

  ! Leith parameters
  logical :: use_QG_Leith_GM      !< If true, uses the QG Leith viscosity as the GM coefficient
  logical :: use_beta_in_QG_Leith !< If true, includes the beta term in the QG Leith GM coefficient

  ! Diagnostics
  !>@{
  !! Diagnostic identifier
  integer :: id_SN_u=-1, id_SN_v=-1, id_L2u=-1, id_L2v=-1, id_Res_fn = -1
  integer :: id_N2_u=-1, id_N2_v=-1, id_S2_u=-1, id_S2_v=-1
  integer :: id_dzu=-1, id_dzv=-1, id_dzSxN=-1, id_dzSyN=-1
  integer :: id_Rd_dx=-1, id_KH_u_QG = -1, id_KH_v_QG = -1
  integer :: id_sqg_struct=-1, id_BS_struct=-1, id_khth_struct=-1, id_khtr_struct=-1
  integer :: id_kdgl90_struct=-1
  type(diag_ctrl), pointer :: diag !< A structure that is used to regulate the
                                   !! timing of diagnostic output.
  !>@}

  type(wave_speed_CS) :: wave_speed !< Wave speed control structure
  type(group_pass_type) :: pass_cg1 !< For group halo pass
  logical :: debug      !< If true, write out checksums of data for debugging
end type VarMix_CS

public VarMix_init, VarMix_end, calc_slope_functions, calc_resoln_function
public calc_QG_slopes, calc_QG_Leith_viscosity, calc_depth_function, calc_sqg_struct


  interface
module subroutine calc_depth_function(G, CS)
  type(ocean_grid_type),  intent(in)    :: G  !< Ocean grid structure
  type(VarMix_CS),        intent(inout) :: CS !< Variable mixing control structure

  ! Local variables
end subroutine calc_depth_function
module subroutine calc_resoln_function(h, tv, G, GV, US, CS, MEKE, OBC, dt)
  type(ocean_grid_type),                     intent(inout) :: G  !< Ocean grid structure
  type(verticalGrid_type),                   intent(in)    :: GV !< Vertical grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in)    :: h  !< Layer thickness [H ~> m or kg m-2]
  type(thermo_var_ptrs),                     intent(in)    :: tv !< Thermodynamic variables
  type(unit_scale_type),                     intent(in)    :: US !< A dimensional unit scaling type
  type(VarMix_CS),                           intent(inout) :: CS !< Variable mixing control structure
  type(MEKE_type),                           intent(in)    :: MEKE !< MEKE struct
  type(ocean_OBC_type),                      pointer       :: OBC !< Open boundaries control structure
  real,                                      intent(in)    :: dt !< Time increment [T ~> s]

  ! Local variables
  ! Depending on the power-function being used, dimensional rescaling may be limited, so some
  ! of the following variables have units that depend on that power.
end subroutine calc_resoln_function
module subroutine calc_sqg_struct(h, tv, G, GV, US, CS, dt, MEKE, OBC)
  type(ocean_grid_type),                     intent(inout) :: G  !< Ocean grid structure
  type(verticalGrid_type),                   intent(in)    :: GV !< Vertical grid structure
  type(unit_scale_type),                     intent(in)    :: US !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in)    :: h  !< Layer thickness [H ~> m or kg m-2]
  type(thermo_var_ptrs),                     intent(in)    :: tv !<Thermodynamic variables
   real,                                     intent(in)    :: dt !< Time increment [T ~> s]
  type(VarMix_CS),                           intent(inout) :: CS !< Variable mixing control struct
  type(MEKE_type),                           intent(in)    :: MEKE !< MEKE struct
  type(ocean_OBC_type),                      pointer       :: OBC !< Open boundaries control structure

  ! Local variables

                                    ! length scale [T L-1 ~> s m-1]
                                    ! state interpolated to an interface [R C-1 ~> kg m-3 ppt-1]
                                    ! state interpolated [R C-1 ~> kg m-3 degC-1]
                        ! [L2 Z-1 T-2 R-1 ~> m4 s-2 kg-1]

end subroutine calc_sqg_struct
module subroutine calc_slope_functions(h, tv, dt, G, GV, US, CS, OBC)
  type(ocean_grid_type),                     intent(inout) :: G  !< Ocean grid structure
  type(verticalGrid_type),                   intent(in)    :: GV !< Vertical grid structure
  type(unit_scale_type),                     intent(in)    :: US !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(inout) :: h  !< Layer thickness [H ~> m or kg m-2]
  type(thermo_var_ptrs),                     intent(in)    :: tv !< Thermodynamic variables
  real,                                      intent(in)    :: dt !< Time increment [T ~> s]
  type(VarMix_CS),                           intent(inout) :: CS !< Variable mixing control structure
  type(ocean_OBC_type),                      pointer       :: OBC !< Open boundaries control structure

  ! Local variables

end subroutine calc_slope_functions
module subroutine calc_Visbeck_coeffs_old(h, slope_x, slope_y, N2_u, N2_v, G, GV, US, CS, OBC)
  type(ocean_grid_type),                        intent(inout) :: G  !< Ocean grid structure
  type(verticalGrid_type),                      intent(in)    :: GV !< Vertical grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),    intent(in)    :: h  !< Layer thickness [H ~> m or kg m-2]
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)+1), intent(in)    :: slope_x !< Zonal isoneutral slope [Z L-1 ~> nondim]
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)+1), intent(in)    :: N2_u    !< Buoyancy (Brunt-Vaisala) frequency
                                                                         !! at u-points [L2 Z-2 T-2 ~> s-2]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)+1), intent(in)    :: slope_y !< Meridional isoneutral slope
                                                                         !! [Z L-1 ~> nondim]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)+1), intent(in)    :: N2_v    !< Buoyancy (Brunt-Vaisala) frequency
                                                                         !! at v-points [L2 Z-2 T-2 ~> s-2]
  type(unit_scale_type),                        intent(in)    :: US !< A dimensional unit scaling type
  type(VarMix_CS),                              intent(inout) :: CS !< Variable mixing control structure
  type(ocean_OBC_type),                         pointer       :: OBC  !< Open boundaries control structure.

  ! Local variables

  ! Note that at some points in the code S2_u and S2_v hold the running depth
  ! integrals of the squared slope [H ~> m or kg m-2] before the average is taken.
                                 ! slope [H Z2 L-2 ~> m or kg m-2] and then the average of the
                                 ! squared slope [Z2 L-2 ~> nondim] at u points.
                                 ! slope [H Z2 L-2 ~> m or kg m-2] and then the average of the
                                 ! squared slope [Z2 L-2 ~> nondim] at v points.
                                 ! eastern OBCs, -1 for western OBCs and 0 at points with no OBCs.
                                 ! northern OBCs, -1 for southern OBCs and 0 at points with no OBCs.
                                 ! interface or the inward equivalent with OBCs [H4 ~> m4 or kg4 m-8]
                                 ! interface or the inward equivalent with OBCs [H4 ~> m4 or kg4 m-8]

end subroutine calc_Visbeck_coeffs_old
module subroutine calc_Eady_growth_rate_2D(CS, G, GV, US, h, e, dzu, dzv, dzSxN, dzSyN, SN_u, SN_v)
  type(VarMix_CS),                              intent(inout) :: CS !< Variable mixing coefficients
  type(ocean_grid_type),                        intent(in) :: G   !< Ocean grid structure
  type(verticalGrid_type),                      intent(in) :: GV  !< Vertical grid structure
  type(unit_scale_type),                        intent(in) :: US  !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),    intent(in) :: h   !< Interface height [Z ~> m]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)+1),  intent(in) :: e   !< Interface height [Z ~> m]
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)+1), intent(in) :: dzu !< dz at u-points [Z ~> m]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)+1), intent(in) :: dzv !< dz at v-points [Z ~> m]
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)+1), intent(in) :: dzSxN !< dz Sx N at u-points [Z T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)+1), intent(in) :: dzSyN !< dz Sy N at v-points [Z T-1 ~> m s-1]
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)+1), intent(inout) :: SN_u !< SN at u-points [T-1 ~> s-1]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)+1), intent(inout) :: SN_v !< SN at v-points [T-1 ~> s-1]
  ! Local variables

end subroutine calc_Eady_growth_rate_2D
module subroutine calc_slope_functions_using_just_e(h, G, GV, US, CS, e)
  type(ocean_grid_type),                       intent(inout) :: G  !< Ocean grid structure
  type(verticalGrid_type),                     intent(in)    :: GV !< Vertical grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),   intent(inout) :: h  !< Layer thickness [H ~> m or kg m-2]
  type(unit_scale_type),                       intent(in)    :: US !< A dimensional unit scaling type
  type(VarMix_CS),                             intent(inout) :: CS !< Variable mixing control structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)+1), intent(in)    :: e  !< Interface position [Z ~> m]
  ! type(thermo_var_ptrs),                     intent(in)    :: tv !< Thermodynamic variables
  ! Local variables
  ! real :: dz(SZI_(G),SZJ_(G),SZK_(GV)) ! The vertical distance across each layer [Z ~> m]
                        ! in roundoff and can be neglected [H ~> m or kg m-2].
                        ! the buoyancy frequency squared at u-points [Z T-2 ~> m s-2]
                        ! the buoyancy frequency squared at v-points [Z T-2 ~> m s-2]
                        ! bathymetric depth for certain calculations.

end subroutine calc_slope_functions_using_just_e
module subroutine calc_QG_slopes(h, tv, dt, G, GV, US, slope_x, slope_y, CS, OBC)
  type(ocean_grid_type),                        intent(in)    :: G  !< Ocean grid structure
  type(verticalGrid_type),                      intent(in)    :: GV !< Vertical grid structure
  type(unit_scale_type),                        intent(in)    :: US !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),    intent(in)    :: h  !< Layer thickness [H ~> m or kg m-2]
  type(thermo_var_ptrs),                        intent(in)    :: tv !< Thermodynamic variables
  real,                                         intent(in)    :: dt !< Time increment [T ~> s]
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)+1), intent(inout) :: slope_x !< Isopycnal slope in i-dir [Z L-1 ~> nondim]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)+1), intent(inout) :: slope_y !< Isopycnal slope in j-dir [Z L-1 ~> nondim]
  type(VarMix_CS),                              intent(in)    :: CS !< Variable mixing control structure
  type(ocean_OBC_type),                         pointer       :: OBC !< Open boundaries control structure
  ! Local variables

end subroutine calc_QG_slopes
module subroutine calc_QG_Leith_viscosity(CS, G, GV, US, h, dz, k, div_xx_dx, div_xx_dy, slope_x, slope_y, &
                                   vort_xy_dx, vort_xy_dy)
  type(VarMix_CS),                           intent(inout) :: CS !< Variable mixing coefficients
  type(ocean_grid_type),                     intent(in)    :: G  !< Ocean grid structure
  type(verticalGrid_type),                   intent(in)    :: GV !< The ocean's vertical grid structure.
  type(unit_scale_type),                     intent(in)    :: US !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in)    :: h  !< Layer thickness [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in)    :: dz !< Layer vertical extents [Z ~> m]
  integer,                                   intent(in)    :: k  !< Layer for which to calculate vorticity magnitude
  real, dimension(SZIB_(G),SZJ_(G)),         intent(in)    :: div_xx_dx  !< x-derivative of horizontal divergence
                                                                 !! (d/dx(du/dx + dv/dy)) [L-1 T-1 ~> m-1 s-1]
  real, dimension(SZI_(G),SZJB_(G)),         intent(in)    :: div_xx_dy  !< y-derivative of horizontal divergence
                                                                 !! (d/dy(du/dx + dv/dy)) [L-1 T-1 ~> m-1 s-1]
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)+1), intent(inout) :: slope_x !< Isopycnal slope in i-dir [Z L-1 ~> nondim]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)+1), intent(inout) :: slope_y !< Isopycnal slope in j-dir [Z L-1 ~> nondim]
  real, dimension(SZI_(G),SZJB_(G)),         intent(inout) :: vort_xy_dx !< x-derivative of vertical vorticity
                                                                 !! (d/dx(dv/dx - du/dy)) [L-1 T-1 ~> m-1 s-1]
  real, dimension(SZIB_(G),SZJ_(G)),         intent(inout) :: vort_xy_dy !< y-derivative of vertical vorticity
                                                                 !! (d/dy(dv/dx - du/dy)) [L-1 T-1 ~> m-1 s-1]
  ! Local variables

                  ! mass-weighted average specific volumes around an interface [H Z-1 ~> nondim or kg m-3]

end subroutine calc_QG_Leith_viscosity
module subroutine VarMix_init(Time, G, GV, US, param_file, diag, CS)
  type(time_type),            intent(in) :: Time !< Current model time
  type(ocean_grid_type),      intent(in) :: G    !< Ocean grid structure
  type(verticalGrid_type),    intent(in) :: GV   !< The ocean's vertical grid structure
  type(unit_scale_type),      intent(in) :: US   !< A dimensional unit scaling type
  type(param_file_type),      intent(in) :: param_file !< Parameter file handles
  type(diag_ctrl), target, intent(inout) :: diag !< Diagnostics control structure
  type(VarMix_CS),         intent(inout) :: CS   !< Variable mixing coefficients

  ! Local variables
                         ! for the epipycnal tracer diffusivity [nondim]
                         ! for the interface depth diffusivity [nondim]
                   ! of the equatorial deformation radius us used [nondim]
                           ! calculating the first-mode wave speed [H ~> m or kg m-2]
                               ! mixing and interface height mixing [nondim]
             ! default value is roughly (pi / (the age of the universe)).
                                  ! for remapping.  Values below 20190101 recover the remapping
                                  ! answers from 2018, while higher values use more robust
                                  ! forms of the same remapping expressions.
                                  ! buoyancy gradients in boundary layer parameterizations [L ~> m]
                                  ! scaled by the resolution function.
                              ! mode wave speed as the starting point for iterations.
                           ! temperature variance [nondim]
                           ! recreate the bugs, or if false bugs are only used if actively selected.
                           ! lateral mixing coefficient calculations and to calculate stratification
                           ! and other fields at open boundary condition faces.
  ! This include declares and sets the variable "version".
end subroutine VarMix_init
module subroutine VarMix_end(CS)
  type(VarMix_CS), intent(inout) :: CS

end subroutine VarMix_end
  end interface

end module MOM_lateral_mixing_coeffs
