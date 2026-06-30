! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> Shear-dependent mixing following Jackson et al. 2008.
module MOM_kappa_shear

use MOM_cpu_clock,         only : cpu_clock_id, cpu_clock_begin, cpu_clock_end
use MOM_cpu_clock,         only : CLOCK_MODULE_DRIVER, CLOCK_MODULE, CLOCK_ROUTINE
use MOM_diag_mediator,     only : post_data, register_diag_field, safe_alloc_ptr
use MOM_diag_mediator,     only : diag_ctrl, time_type
use MOM_debugging,         only : hchksum, Bchksum
use MOM_error_handler,     only : MOM_error, is_root_pe, FATAL, WARNING, NOTE
use MOM_file_parser,       only : get_param, log_version, param_file_type
use MOM_grid,              only : ocean_grid_type
use MOM_interface_heights, only : thickness_to_dz
use MOM_unit_scaling,      only : unit_scale_type
use MOM_variables,         only : thermo_var_ptrs
use MOM_verticalGrid,      only : verticalGrid_type
use MOM_EOS,               only : calculate_density_derivs
use MOM_EOS,               only : calculate_density, calculate_specific_vol_derivs

implicit none ; private

#include <MOM_memory.h>

public Calculate_kappa_shear, Calc_kappa_shear_vertex, kappa_shear_init
public kappa_shear_is_used, kappa_shear_at_vertex

! A note on unit descriptions in comments: MOM6 uses units that can be rescaled for dimensional
! consistency testing. These are noted in comments with units like Z, H, L, and T, along with
! their mks counterparts with notation like "a velocity [Z T-1 ~> m s-1]".  If the units
! vary with the Boussinesq approximation, the Boussinesq variant is given first.

!> This control structure holds the parameters that regulate shear mixing
type, public :: Kappa_shear_CS ; private
  real    :: RiNo_crit       !< The critical shear Richardson number for
                             !! shear-entrainment [nondim]. The theoretical value is 0.25.
                             !! The values found by Jackson et al. are 0.25-0.35.
  real    :: Shearmix_rate   !< A nondimensional rate scale for shear-driven
                             !! entrainment [nondim].  The value given by Jackson et al.
                             !! is 0.085-0.089.
  real    :: FRi_curvature   !<   A constant giving the curvature of the function
                             !! of the Richardson number that relates shear to
                             !! sources in the kappa equation [nondim].
                             !! The values found by Jackson et al. are -0.97 - -0.89.
  real    :: C_N             !<   The coefficient for the decay of TKE due to
                             !! stratification (i.e. proportional to N*tke) [nondim].
                             !! The values found by Jackson et al. are 0.24-0.28.
  real    :: C_S             !<   The coefficient for the decay of TKE due to
                             !! shear (i.e. proportional to |S|*tke) [nondim].
                             !! The values found by Jackson et al. are 0.14-0.12.
  real    :: lambda          !<   The coefficient for the buoyancy length scale
                             !! in the kappa equation [nondim].
                             !! The values found by Jackson et al. are 0.82-0.81.
  real    :: lambda2_N_S     !<   The square of the ratio of the coefficients of
                             !! the buoyancy and shear scales in the diffusivity
                             !! equation, 0 to eliminate the shear scale [nondim].
  real    :: lz_rescale      !<   A coefficient to rescale the distance to the nearest
                             !! solid boundary. This adjustment is to account for
                             !! regions where 3 dimensional turbulence prevents the
                             !! growth of shear instabilities [nondim].
  real    :: TKE_bg          !<   The background level of TKE [Z2 T-2 ~> m2 s-2].
  real    :: kappa_0         !<   The background diapycnal diffusivity [H Z T-1 ~> m2 s-1 or Pa s]
  real    :: kappa_seed      !<   A moderately large seed value of diapycnal diffusivity that
                             !! is used as a starting turbulent diffusivity in the iterations
                             !! to finding an energetically constrained solution for the
                             !! shear-driven diffusivity [H Z T-1 ~> m2 s-1 or Pa s]
  real    :: kappa_trunc     !< Diffusivities smaller than this are rounded to 0 [H Z T-1 ~> m2 s-1 or Pa s]
  real    :: kappa_tol_err   !<   The fractional error in kappa that is tolerated [nondim].
  real    :: Prandtl_turb    !< Prandtl number used to convert Kd_shear into viscosity [nondim].
  integer :: nkml            !<   The number of layers in the mixed layer, as
                             !! treated in this routine.  If the pieces of the
                             !! mixed layer are not to be treated collectively,
                             !! nkml is set to 1.
  integer :: max_RiNo_it     !< The maximum number of iterations that may be used
                             !! to estimate the instantaneous shear-driven mixing.
  integer :: max_KS_it       !< The maximum number of iterations that may be used
                             !! to estimate the time-averaged diffusivity.
  logical :: dKdQ_iteration_bug !< If true. use an older, dimensionally inconsistent estimate of
                             !! the derivative of diffusivity with energy in the Newton's method
                             !! iteration.  The bug causes under-corrections when dz > 1m.
  logical :: KS_at_vertex    !< If true, do the calculations of the shear-driven mixing
                             !! at the cell vertices (i.e., the vorticity points).
  logical :: eliminate_massless !< If true, massless layers are merged with neighboring
                             !! massive layers in this calculation.
                             !  I can think of no good reason why this should be false. - RWH
  real    :: vel_underflow   !< Velocity components smaller than vel_underflow
                             !! are set to 0 [L T-1 ~> m s-1].
  real    :: kappa_src_max_chg !< The maximum permitted increase in the kappa source within an
                             !! iteration relative to the local source [nondim].  This must be
                             !! greater than 1.  The lower limit for the permitted fractional
                             !! decrease is (1 - 0.5/kappa_src_max_chg).  These limits could
                             !! perhaps be made dynamic with an improved iterative solver.
  real    :: VS_GeoMean_Kdmin !< A minimum diffusivity for computing the horizontal averages
                             !! when using the geometric mean with VERTEX_SHEAR=True.  The model
                             !! is sensitive to this value, which is a drawback of using the
                             !! geometric average as currently implemented.
  logical :: psurf_bug       !< If true, do a simple average of the cell surface pressures to get a
                             !! surface pressure at the corner if VERTEX_SHEAR=True.  Otherwise mask
                             !! out any land points in the average.
  logical :: all_layer_TKE_bug !< If true, report back the latest estimate of TKE instead of the
                             !! time average TKE when there is mass in all layers.  Otherwise always
                             !! report the time-averaged TKE, as is currently done when there
                             !! are some massless layers.
  logical :: VS_viscosity_bug !< If true, use a bug in the calculation of the viscosity that sets
                             !! it to zero for all vertices that are on a coastline.
  logical :: vertex_shear_OBC_bug !< If false, use extra masking when interpolating thicknesses to velocity
                             !! points for setting up the shear velocities at vertices to avoid using
                             !! external thicknesses at open boundaries.  When OBCs are not in use,
                             !! this parameter does not change answers, but true is more efficient.
  logical :: VS_GeometricMean !< If true use geometric averaging for Kd from vertices to tracer points
  logical :: VS_ThicknessMean !< If true use thickness weighting when averaging Kd from vertices to
                             !! tracer points
  logical :: restrictive_tolerance_check !< If false, uses the less restrictive tolerance check to
                             !! determine if a timestep is acceptable for the KS_it outer iteration
                             !! loop, as the code was originally written.  True uses the more
                             !! restrictive check.
!  logical :: layer_stagger = .false. ! If true, do the calculations centered at
                             !  layers, rather than the interfaces.
  logical :: debug = .false. !< If true, write verbose debugging messages.
  type(diag_ctrl), pointer :: diag => NULL() !< A structure that is used to
                             !! regulate the timing of diagnostic output.
  !>@{ Diagnostic IDs
  integer :: id_Kd_shear = -1, id_TKE = -1, id_Kd_vertex = -1, &
             id_S2_init = -1, id_N2_init = -1, id_S2_mean = -1, id_N2_mean = -1
  !>@}
end type Kappa_shear_CS

! integer :: id_clock_project, id_clock_KQ, id_clock_avg, id_clock_setup


  interface
module subroutine Calculate_kappa_shear(u_in, v_in, h, tv, p_surf, kappa_io, tke_io, &
                                 kv_io, dt, G, GV, US, CS)
  type(ocean_grid_type),   intent(in)    :: G      !< The ocean's grid structure.
  type(verticalGrid_type), intent(in)    :: GV     !< The ocean's vertical grid structure.
  type(unit_scale_type),   intent(in)    :: US     !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),   &
                           intent(in)    :: u_in   !< Initial zonal velocity [L T-1 ~> m s-1].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),   &
                           intent(in)    :: v_in   !< Initial meridional velocity [L T-1 ~> m s-1].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),   &
                           intent(in)    :: h      !< Layer thicknesses [H ~> m or kg m-2].
  type(thermo_var_ptrs),   intent(in)    :: tv     !< A structure containing pointers to any
                                                   !! available thermodynamic fields. Absent fields
                                                   !! have NULL ptrs.
  real, dimension(:,:),    pointer       :: p_surf !< The pressure at the ocean surface [R L2 T-2 ~> Pa] (or NULL).
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)+1), &
                           intent(inout) :: kappa_io !< The diapycnal diffusivity at each interface
                                                   !! [H Z T-1 ~> m2 s-1 or kg m-1 s-1].  Initially this
                                                   !! is the value from the previous timestep, which may
                                                   !! accelerate the iteration toward convergence.
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)+1), &
                           intent(out) :: tke_io   !< The turbulent kinetic energy per unit mass at
                                                   !! each interface (not layer!) [Z2 T-2 ~> m2 s-2].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)+1), &
                           intent(inout) :: kv_io  !< The vertical viscosity at each interface
                                                   !! (not layer!) [H Z T-1 ~> m2 s-1 or Pa s]. This discards any
                                                   !! previous value (i.e. it is intent out) and
                                                   !! simply sets Kv = Prandtl * Kd_shear
  real,                    intent(in)    :: dt     !< Time increment [T ~> s].
  type(Kappa_shear_CS),    pointer       :: CS     !< The control structure returned by a previous
                                                   !! call to kappa_shear_init.

  ! Local variables
end subroutine Calculate_kappa_shear
module subroutine Calc_kappa_shear_vertex(u_in, v_in, h, T_in, S_in, tv, p_surf, kappa_io, tke_io, &
                                   kv_io, dt, G, GV, US, CS)
  type(ocean_grid_type),   intent(in)    :: G      !< The ocean's grid structure.
  type(verticalGrid_type), intent(in)    :: GV     !< The ocean's vertical grid structure.
  type(unit_scale_type),    intent(in)   :: US     !< A dimensional unit scaling type
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)),   &
                           intent(in)    :: u_in   !< Initial zonal velocity [L T-1 ~> m s-1].
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)),   &
                           intent(in)    :: v_in   !< Initial meridional velocity [L T-1 ~> m s-1].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),   &
                           intent(in)    :: h      !< Layer thicknesses [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),   &
                           intent(in)    :: T_in   !< Layer potential temperatures [C ~> degC]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),   &
                           intent(in)    :: S_in   !< Layer salinities [S ~> ppt]
  type(thermo_var_ptrs),   intent(in)    :: tv     !< A structure containing pointers to any
                                                   !! available thermodynamic fields. Absent fields
                                                   !! have NULL ptrs.
  real, dimension(:,:),    pointer       :: p_surf !< The pressure at the ocean surface [R L2 T-2 ~> Pa]
                                                   !! (or NULL).
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)+1), &
                           intent(out)   :: kappa_io !< The diapycnal diffusivity at each interface
                                                   !! (not layer!) [H Z T-1 ~> m2 s-1 or kg m-1 s-1].
  real, dimension(SZIB_(G),SZJB_(G),SZK_(GV)+1), &
                           intent(out)   :: tke_io !< The turbulent kinetic energy per unit mass at
                                                   !! each interface (not layer!) [Z2 T-2 ~> m2 s-2].
  real, dimension(SZIB_(G),SZJB_(G),SZK_(GV)+1), &
                           intent(inout) :: kv_io  !< The vertical viscosity at each interface
                                                   !! [H Z T-1 ~> m2 s-1 or Pa s].
                                                   !! The previous value is used to initialize kappa
                                                   !! in the vertex columns as Kappa = Kv/Prandtl
                                                   !! to accelerate the iteration toward convergence.
  real,                    intent(in)    :: dt     !< Time increment [T ~> s].
  type(Kappa_shear_CS),    pointer       :: CS     !< The control structure returned by a previous
                                                   !! call to kappa_shear_init.

  ! Local variables



                        ! allocated and are being used as state variables.

                        ! interfaces and the interfaces with massless layers
                        ! merged into nearby massive layers.
                        ! interpolating back to the original index space [nondim].
                        ! diffusivity is taken using thickness weighted powers [H Z s m-2 T-1 ~> 1]
                        ! or [H Z m s kg-1 T-1 ~> 1]
                        ! thickness-weighted averages [H ~> m or kg m-2]

  ! Diagnostics that should be deleted?
end subroutine Calc_kappa_shear_vertex
module subroutine kappa_shear_column(kappa, tke, dt, nzc, f2, surface_pres, hlay, dz_lay, &
                              u0xdz, v0xdz, T0xdz, S0xdz, kappa_avg, tke_avg, N2_init, S2_init, &
                              N2_mean, S2_mean, tv, CS, GV, US )
  type(verticalGrid_type), intent(in)    :: GV !< The ocean's vertical grid structure.
  real, dimension(SZK_(GV)+1), &
                     intent(inout) :: kappa !< The time-weighted average of kappa [H Z T-1 ~> m2 s-1 or Pa s]
  real, dimension(SZK_(GV)+1), &
                     intent(out)   :: tke  !< The Turbulent Kinetic Energy per unit mass at
                                           !! an interface [Z2 T-2 ~> m2 s-2].
  integer,           intent(in)    :: nzc  !< The number of active layers in the column.
  real,              intent(in)    :: f2   !< The square of the Coriolis parameter [T-2 ~> s-2].
  real,              intent(in)    :: surface_pres  !< The surface pressure [R L2 T-2 ~> Pa].
  real, dimension(SZK_(GV)), &
                     intent(in)    :: hlay  !< The layer thickness [H ~> m or kg m-2]
  real, dimension(SZK_(GV)), &
                     intent(in)    :: dz_lay !< The geometric layer thickness in height units [Z ~> m]
  real, dimension(SZK_(GV)), &
                     intent(in)    :: u0xdz !< The initial zonal velocity times hlay [H L T-1 ~> m2 s-1 or kg m-1 s-1]
  real, dimension(SZK_(GV)), &
                     intent(in)    :: v0xdz !< The initial meridional velocity times the
                                            !! layer thickness [H L T-1 ~> m2 s-1 or kg m-1 s-1]
  real, dimension(SZK_(GV)), &
                     intent(in)    :: T0xdz !< The initial temperature times hlay [C H ~> degC m or degC kg m-2]
  real, dimension(SZK_(GV)), &
                     intent(in)    :: S0xdz !< The initial salinity times hlay [S H ~> ppt m or ppt kg m-2]
  real, dimension(SZK_(GV)+1), &
                     intent(out)   :: kappa_avg !< The time-weighted average of kappa [H Z T-1 ~> m2 s-1 or Pa s]
  real, dimension(SZK_(GV)+1), &
                     intent(out)   :: tke_avg  !< The time-weighted average of TKE [Z2 T-2 ~> m2 s-2].
  real, dimension(SZK_(GV)+1), &
                     intent(out)   :: N2_mean  !< The time-weighted average of N2 [Z2 T-2 ~> m2 s-2].
  real, dimension(SZK_(GV)+1), &
                     intent(out)   :: S2_mean  !< The time-weighted average of S2 [Z2 T-2 ~> m2 s-2].
  real, dimension(SZK_(GV)+1), &
                     intent(out)   :: N2_init  !< The initial value of N2 [Z2 T-2 ~> m2 s-2].
  real, dimension(SZK_(GV)+1), &
                     intent(out)   :: S2_init  !< The initial value of S2 [Z2 T-2 ~> m2 s-2].
  real,                    intent(in)    :: dt !< Time increment [T ~> s].
  type(thermo_var_ptrs),   intent(in)    :: tv !< A structure containing pointers to any
                                               !! available thermodynamic fields. Absent fields
                                               !! have NULL ptrs.
  type(Kappa_shear_CS),    pointer       :: CS !< The control structure returned by a previous
                                               !! call to kappa_shear_init.
  type(unit_scale_type),   intent(in)    :: US !< A dimensional unit scaling type

  ! Local variables

end subroutine kappa_shear_column
module subroutine calculate_projected_state(kappa, u0, v0, T0, S0, dt, nz, dz, I_dz_int, dbuoy_dT, dbuoy_dS, &
                                     vel_under, u, v, T, Sal, N2, S2, GV, US, ks_int, ke_int)
  integer,               intent(in)    :: nz  !< The number of layers (after eliminating massless
                                              !! layers?).
  real, dimension(nz+1), intent(in)    :: kappa !< The diapycnal diffusivity at interfaces,
                                              !! [H Z T-1 ~> m2 s-1 or Pa s].
  real, dimension(nz),   intent(in)    :: u0  !< The initial zonal velocity [L T-1 ~> m s-1].
  real, dimension(nz),   intent(in)    :: v0  !< The initial meridional velocity [L T-1 ~> m s-1].
  real, dimension(nz),   intent(in)    :: T0  !< The initial temperature [C ~> degC].
  real, dimension(nz),   intent(in)    :: S0  !< The initial salinity [S ~> ppt].
  real,                  intent(in)    :: dt  !< The time step [T ~> s].
  real, dimension(nz),   intent(in)    :: dz  !< The layer thicknesses [H ~> m or kg m-2]
  real, dimension(nz+1), intent(in)    :: I_dz_int !< The inverse of the distance between successive
                                              !! layer centers [Z-1 ~> m-1].
  real, dimension(nz+1), intent(in)    :: dbuoy_dT !< The partial derivative of buoyancy with
                                              !! temperature [Z T-2 C-1 ~> m s-2 degC-1].
  real, dimension(nz+1), intent(in)    :: dbuoy_dS !< The partial derivative of buoyancy with
                                              !! salinity [Z T-2 S-1 ~> m s-2 ppt-1].
  real,                  intent(in)    :: vel_under !< Any velocities that are smaller in magnitude
                                              !! than this value are set to 0 [L T-1 ~> m s-1].
  real, dimension(nz),   intent(inout) :: u   !< The zonal velocity after dt [L T-1 ~> m s-1].
  real, dimension(nz),   intent(inout) :: v   !< The meridional velocity after dt [L T-1 ~> m s-1].
  real, dimension(nz),   intent(inout) :: T   !< The temperature after dt [C ~> degC].
  real, dimension(nz),   intent(inout) :: Sal !< The salinity after dt [S ~> ppt].
  real, dimension(nz+1), intent(inout) :: N2  !< The buoyancy frequency squared at interfaces [T-2 ~> s-2].
  real, dimension(nz+1), intent(inout) :: S2  !< The squared shear at interfaces [T-2 ~> s-2].
  type(verticalGrid_type), intent(in)  :: GV  !< The ocean's vertical grid structure.
  type(unit_scale_type), intent(in)    :: US  !< A dimensional unit scaling type
  integer, optional,     intent(in)    :: ks_int !< The topmost k-index with a non-zero diffusivity.
  integer, optional,     intent(in)    :: ke_int !< The bottommost k-index with a non-zero
                                              !! diffusivity.

  ! Local variables

end subroutine calculate_projected_state
module subroutine find_kappa_tke(N2, S2, kappa_in, Idz, h_Int, dz_Int, dz_h_Int, I_L2_bdry, f2, &
                          nz, CS, GV, US, K_Q, tke, kappa, kappa_src, local_src)
  integer,               intent(in)    :: nz  !< The number of layers to work on.
  real, dimension(nz+1), intent(in)    :: N2  !< The buoyancy frequency squared at interfaces [T-2 ~> s-2].
  real, dimension(nz+1), intent(in)    :: S2  !< The squared shear at interfaces [T-2 ~> s-2].
  real, dimension(nz+1), intent(in)    :: kappa_in  !< The initial guess at the diffusivity
                                              !! [H Z T-1 ~> m2 s-1 or Pa s]
  real, dimension(nz+1), intent(in)    :: h_Int !< The thicknesses associated with interfaces
                                              !! [H ~> m or kg m-2]
  real, dimension(nz+1), intent(in)    :: dz_Int !< The vertical distances around interfaces [Z ~> m]
  real, dimension(nz+1), intent(in)    :: dz_h_Int !< The ratio of the vertical distances to the
                                              !! thickness around an interface [Z H-1 ~> nondim or m3 kg-1].
                                              !! In non-Boussinesq mode this is the specific volume.
  real, dimension(nz+1), intent(in)    :: I_L2_bdry !< The inverse of the squared distance to
                                              !! boundaries [H-1 Z-1 ~> m-2 or m kg-1].
  real, dimension(nz),   intent(in)    :: Idz !< The inverse grid spacing of layers [Z-1 ~> m-1].
  real,                  intent(in)    :: f2  !< The squared Coriolis parameter [T-2 ~> s-2].
  type(Kappa_shear_CS),  pointer       :: CS  !< A pointer to this module's control structure.
  type(verticalGrid_type), intent(in)  :: GV  !< The ocean's vertical grid structure.
  type(unit_scale_type), intent(in)    :: US  !< A dimensional unit scaling type
  real, dimension(nz+1), intent(inout) :: K_Q !< The shear-driven diapycnal diffusivity divided by
                                              !! the turbulent kinetic energy per unit mass at
                                              !! interfaces [H T Z-1 ~> s or kg s m-3].
  real, dimension(nz+1), intent(out)   :: tke !< The turbulent kinetic energy per unit mass at
                                              !! interfaces [Z2 T-2 ~> m2 s-2].
  real, dimension(nz+1), intent(out)   :: kappa !< The diapycnal diffusivity at interfaces
                                              !! [H Z T-1 ~> m2 s-1 or Pa s]
  real, dimension(nz+1), optional, &
                         intent(out)   :: kappa_src !< The source term for kappa [T-1 ~> s-1]
  real, dimension(nz+1), optional, &
                         intent(out)   :: local_src !< The sum of all local sources for kappa
                                              !! [T-1 ~> s-1]
  ! This subroutine calculates new, consistent estimates of TKE and kappa.

  ! Local variables
end subroutine find_kappa_tke
module function kappa_shear_init(Time, G, GV, US, param_file, diag, CS)
  type(time_type),         intent(in)    :: Time !< The current model time.
  type(ocean_grid_type),   intent(in)    :: G    !< The ocean's grid structure.
  type(verticalGrid_type), intent(in)    :: GV   !< The ocean's vertical grid structure.
  type(unit_scale_type),   intent(in)    :: US   !< A dimensional unit scaling type
  type(param_file_type),   intent(in)    :: param_file !< A structure to parse for run-time
                                                 !! parameters.
  type(diag_ctrl), target, intent(inout) :: diag !< A structure that is used to regulate diagnostic
                                                 !! output.
  type(Kappa_shear_CS),    pointer       :: CS   !< A pointer that is set to point to the control
                                                 !! structure for this module
  logical :: kappa_shear_init !< True if module is to be used, False otherwise

  ! Local variables
                    ! for setting the default of KD_SMOOTH [Z2 T-1 ~> m2 s-1]
                          ! recreate the bugs, or if false bugs are only used if actively selected.
  ! This include declares and sets the variable "version".

end function kappa_shear_init
logical module function kappa_shear_is_used(param_file)
  type(param_file_type), intent(in) :: param_file !< A structure to parse for run-time parameters

  ! Local variables
  ! This function reads the parameter "USE_JACKSON_PARAM" and returns its value.

end function kappa_shear_is_used
logical module function kappa_shear_at_vertex(param_file)
  type(param_file_type), intent(in) :: param_file !< A structure to parse for run-time parameters

  ! Local variables
  ! This function returns true only if the parameters "USE_JACKSON_PARAM" and "VERTEX_SHEAR" are both true.

end function kappa_shear_at_vertex
  end interface

end module MOM_kappa_shear
