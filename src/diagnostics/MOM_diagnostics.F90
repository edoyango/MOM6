! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> Calculates any requested diagnostic quantities
!! that are not calculated in the various subroutines.
!! Diagnostic quantities are requested by allocating them memory.
module MOM_diagnostics

use MOM_coms,              only : reproducing_sum
use MOM_coupler_types,     only : coupler_type_send_data
use MOM_density_integrals, only : int_density_dz
use MOM_diag_mediator,     only : post_data, get_diag_time_end
use MOM_diag_mediator,     only : post_product_u, post_product_sum_u
use MOM_diag_mediator,     only : post_product_v, post_product_sum_v
use MOM_diag_mediator,     only : register_diag_field, register_scalar_field
use MOM_diag_mediator,     only : register_static_field, diag_register_area_ids
use MOM_diag_mediator,     only : diag_ctrl, time_type, safe_alloc_ptr
use MOM_diag_mediator,     only : diag_get_volume_cell_measure_dm_id
use MOM_diag_mediator,     only : diag_grid_storage
use MOM_diag_mediator,     only : diag_save_grids, diag_restore_grids, diag_copy_storage_to_diag
use MOM_domains,           only : create_group_pass, do_group_pass, group_pass_type
use MOM_domains,           only : To_North, To_East
use MOM_EOS,               only : calculate_density, calculate_density_derivs, EOS_domain
use MOM_EOS,               only : cons_temp_to_pot_temp, pot_temp_to_cons_temp
use MOM_EOS,               only : prac_saln_to_abs_saln, abs_saln_to_prac_saln
use MOM_error_handler,     only : MOM_error, FATAL, WARNING
use MOM_file_parser,       only : get_param, log_version, param_file_type
use MOM_grid,              only : ocean_grid_type
use MOM_interface_heights, only : find_eta, find_dz_for_eta, find_col_mass
use MOM_spatial_means,     only : global_area_mean, global_layer_mean
use MOM_spatial_means,     only : global_volume_mean, global_area_integral
use MOM_tracer_registry,   only : tracer_registry_type, post_tracer_transport_diagnostics
use MOM_unit_scaling,      only : unit_scale_type
use MOM_variables,         only : thermo_var_ptrs, ocean_internal_state, p3d
use MOM_variables,         only : accel_diag_ptrs, cont_diag_ptrs, surface
use MOM_verticalGrid,      only : verticalGrid_type, get_thickness_units, get_flux_units
use MOM_wave_speed,        only : wave_speed, wave_speed_CS, wave_speed_init
use Recon1d_EPPM_CWK,      only : EPPM_CWK

implicit none ; private

#include <MOM_memory.h>

public calculate_diagnostic_fields, register_time_deriv, write_static_fields
public register_surface_diags, post_surface_dyn_diags, post_surface_thermo_diags
public register_transport_diags, post_transport_diagnostics
public MOM_diagnostics_init, MOM_diagnostics_end

! A note on unit descriptions in comments: MOM6 uses units that can be rescaled for dimensional
! consistency testing. These are noted in comments with units like Z, H, L, and T, along with
! their mks counterparts with notation like "a velocity [Z T-1 ~> m s-1]".  If the units
! vary with the Boussinesq approximation, the Boussinesq variant is given first.

!> The control structure for the MOM_diagnostics module
type, public :: diagnostics_CS ; private
  logical :: initialized = .false.     !< True if this control structure has been initialized.
  real :: mono_N2_column_fraction = 0. !< The lower fraction of water column over which N2 is limited as
                                       !! monotonic for the purposes of calculating the equivalent
                                       !! barotropic wave speed [nondim].
  real :: mono_N2_depth = -1.          !< The depth below which N2 is limited as monotonic for the purposes of
                                       !! calculating the equivalent barotropic wave speed [H ~> m or kg m-2].
  logical :: accurate_thick_cello      !< If true, use the same careful integrals to find the diagnosed
                                       !! non-Boussinesq layer thicknesses as are used to find the free
                                       !! surface height, instead of using an approximate thickness
                                       !! based on division by the mid-layer density.

  type(diag_ctrl), pointer :: diag => NULL() !< A structure that is used to
                                       !! regulate the timing of diagnostic output.

  ! following arrays store diagnostics calculated here and unavailable outside.

  ! following fields have nz layers.
  real, allocatable :: du_dt(:,:,:) !< net i-acceleration [L T-2 ~> m s-2]
  real, allocatable :: dv_dt(:,:,:) !< net j-acceleration [L T-2 ~> m s-2]
  real, allocatable :: dh_dt(:,:,:) !< thickness rate of change [H T-1 ~> m s-1 or kg m-2 s-1]

  logical :: KE_term_on !< If true, at least one diagnostic term in the KE budget is in use.

  !>@{ Diagnostic IDs
  integer :: id_u   = -1,   id_v   = -1, id_h = -1
  integer :: id_usq = -1,   id_vsq = -1, id_uv = -1
  integer :: id_e              = -1, id_e_D            = -1
  integer :: id_du_dt          = -1, id_dv_dt          = -1
  ! integer :: id_hf_du_dt       = -1, id_hf_dv_dt       = -1
  integer :: id_h_du_dt       = -1, id_h_dv_dt       = -1
  integer :: id_hf_du_dt_2d    = -1, id_hf_dv_dt_2d    = -1
  integer :: id_col_ht         = -1, id_dh_dt          = -1
  integer :: id_KE             = -1, id_dKEdt          = -1
  integer :: id_PE_to_KE       = -1, id_KE_BT          = -1
  integer :: id_KE_SAL         = -1, id_KE_TIDES       = -1
  integer :: id_KE_BT_PF       = -1, id_KE_BT_CF       = -1
  integer :: id_KE_BT_WD       = -1
  integer :: id_PE_to_KE_btbc  = -1, id_KE_Coradv_btbc = -1
  integer :: id_KE_Coradv      = -1, id_KE_adv         = -1
  integer :: id_KE_visc        = -1, id_KE_stress      = -1
  integer :: id_KE_visc_gl90   = -1
  integer :: id_KE_horvisc     = -1, id_KE_dia         = -1
  integer :: id_uh_Rlay        = -1, id_vh_Rlay        = -1
  integer :: id_uhGM_Rlay      = -1, id_vhGM_Rlay      = -1
  integer :: id_h_Rlay         = -1, id_Rd1            = -1
  integer :: id_Rml            = -1, id_Rcv            = -1
  integer :: id_cg1            = -1, id_cfl_cg1        = -1
  integer :: id_cfl_cg1_x      = -1, id_cfl_cg1_y      = -1
  integer :: id_cg_ebt         = -1, id_Rd_ebt         = -1
  integer :: id_p_ebt          = -1
  integer :: id_temp_int       = -1, id_salt_int       = -1
  integer :: id_absscint       = -1, id_pfscint        = -1
  integer :: id_scint          = -1
  integer :: id_chcint         = -1, id_phcint         = -1
  integer :: id_mass_wt        = -1, id_col_mass       = -1
  integer :: id_masscello      = -1, id_masso          = -1
  integer :: id_volcello       = -1
  integer :: id_Tpot           = -1, id_Sprac          = -1
  integer :: id_tob            = -1, id_sob            = -1
  integer :: id_thetaoga       = -1, id_soga           = -1
  integer :: id_bigthetaoga    = -1, id_abssoga        = -1
  integer :: id_sosga          = -1, id_tosga          = -1
  integer :: id_abssosga       = -1, id_bigtosga       = -1
  integer :: id_temp_layer_ave = -1, id_salt_layer_ave = -1
  integer :: id_bigtemp_layer_ave = -1, id_abssalt_layer_ave = -1
  integer :: id_pbo            = -1
  integer :: id_thkcello       = -1, id_rhoinsitu      = -1
  integer :: id_rhopot0        = -1, id_rhopot2        = -1
  integer :: id_drho_dT        = -1, id_drho_dS        = -1
  integer :: id_h_pre_sync     = -1
  integer :: id_tosq           = -1, id_sosq           = -1
  integer :: id_t20d           = -1, id_t17d           = -1

  !>@}
  type(wave_speed_CS) :: wave_speed  !< Wave speed control struct

  type(p3d) :: var_ptr(MAX_FIELDS_)  !< pointers to variables used in the calculation
                                     !! of time derivatives
  type(p3d) :: deriv(MAX_FIELDS_)    !< Time derivatives of various fields
  type(p3d) :: prev_val(MAX_FIELDS_) !< Previous values of variables used in the calculation
                                     !! of time derivatives
  !< previous values of variables used in calculation of time derivatives
  integer   :: nlay(MAX_FIELDS_) !< The number of layers in each diagnostics
  integer   :: num_time_deriv = 0 !< The number of time derivative diagnostics

  type(group_pass_type) :: pass_KE_uv !< A handle used for group halo passes

end type diagnostics_CS


!> A structure with diagnostic IDs of the surface and integrated variables
type, public :: surface_diag_IDs ; private
  !>@{ Diagnostic IDs for 2-d surface and bottom flux and state fields
  !Diagnostic IDs for 2-d surface and bottom fields
  integer :: id_zos  = -1, id_zossq  = -1
  integer :: id_volo = -1, id_speed  = -1
  integer :: id_ssh  = -1, id_ssh_ga = -1
  integer :: id_sst  = -1, id_sst_sq = -1, id_sstcon = -1
  integer :: id_sss  = -1, id_sss_sq = -1, id_sssabs = -1
  integer :: id_ssu  = -1, id_ssv    = -1
  integer :: id_ssu_east = -1, id_ssv_north = -1

  ! Diagnostic IDs for  heat and salt flux fields
  integer :: id_fraz         = -1
  integer :: id_salt_deficit = -1
  integer :: id_Heat_PmE     = -1
  integer :: id_intern_heat  = -1
  !>@}
end type surface_diag_IDs


!> A structure with diagnostic IDs of mass transport related diagnostics
type, public :: transport_diag_IDs ; private
  !>@{  Diagnostics for tracer horizontal transport
  integer :: id_uhtr = -1, id_umo = -1, id_umo_2d = -1
  integer :: id_vhtr = -1, id_vmo = -1, id_vmo_2d = -1
  integer :: id_dynamics_h = -1, id_dynamics_h_tendency = -1
  !>@}
end type transport_diag_IDs



  interface
module subroutine calculate_diagnostic_fields(u, v, h, uh, vh, tv, ADp, CDp, p_surf, &
                                       dt, diag_pre_sync, G, GV, US, CS)
  type(ocean_grid_type),   intent(inout) :: G    !< The ocean's grid structure.
  type(verticalGrid_type), intent(in)    :: GV   !< The ocean's vertical grid structure.
  type(unit_scale_type),   intent(in)    :: US   !< A dimensional unit scaling type
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                           intent(in)    :: u    !< The zonal velocity [L T-1 ~> m s-1].
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                           intent(in)    :: v    !< The meridional velocity [L T-1 ~> m s-1].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in)    :: h    !< Layer thicknesses [H ~> m or kg m-2].
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                           intent(in)    :: uh   !< Transport through zonal faces = u*h*dy,
                                                 !! [H L2 T-1 ~> m3 s-1 or kg s-1].
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                           intent(in)    :: vh   !< Transport through meridional faces = v*h*dx,
                                                 !! [H L2 T-1 ~> m3 s-1 or kg s-1].
  type(thermo_var_ptrs),   intent(in)    :: tv   !< A structure pointing to various
                                                 !! thermodynamic variables.
  type(accel_diag_ptrs),   intent(in)    :: ADp  !< structure with pointers to
                                                 !! accelerations in momentum equation.
  type(cont_diag_ptrs),    intent(in)    :: CDp  !< structure with pointers to
                                                 !! terms in continuity equation.
  real, dimension(:,:),    pointer       :: p_surf !< A pointer to the surface pressure [R L2 T-2 ~> Pa].
                                                 !! If p_surf is not associated, it is the same
                                                 !! as setting the surface pressure to 0.
  real,                    intent(in)    :: dt   !< The time difference since the last
                                                 !! call to this subroutine [T ~> s].
  type(diag_grid_storage), intent(in)    :: diag_pre_sync !< Target grids from previous timestep
  type(diagnostics_CS),    intent(inout) :: CS   !< Control structure returned by a
                                                 !! previous call to diagnostics_init.

  ! Local variables


                                           ! geopotential or the seafloor [Z ~> m].
                                            ! including [nondim] and [H ~> m or kg m-2].
                                           ! overall grid spacing or just one direction [nondim]

  ! tmp array for surface properties
                   ! a list [nondim], scaled so that wt + wt_p = 1.



end subroutine calculate_diagnostic_fields
module subroutine find_weights(Rlist, R_in, k, nz, wt, wt_p)
  real, dimension(:), &
            intent(in)    :: Rlist !< The list of target densities [R ~> kg m-3]
  real,     intent(in)    :: R_in !< The density being inserted into Rlist [R ~> kg m-3]
  integer,  intent(inout) :: k    !< The value of k such that Rlist(k) <= R_in < Rlist(k+1)
                                  !! The input value is a first guess
  integer,  intent(in)    :: nz   !< The number of layers in Rlist
  real,     intent(out)   :: wt   !< The weight of layer k for interpolation [nondim]
  real,     intent(out)   :: wt_p !< The weight of layer k+1 for interpolation [nondim]

  ! This subroutine finds location of R_in in an increasing ordered
  ! list, Rlist, returning as k the element such that
  !  Rlist(k) <= R_in < Rlist(k+1), and where wt and wt_p are the linear
  ! weights that should be assigned to elements k and k+1.


  ! First, bracket the desired point.
end subroutine find_weights
module subroutine calculate_vertical_integrals(h, tv, p_surf, G, GV, US, CS)
  type(ocean_grid_type),   intent(inout) :: G    !< The ocean's grid structure.
  type(verticalGrid_type), intent(in)    :: GV   !< The ocean's vertical grid structure.
  type(unit_scale_type),   intent(in)    :: US   !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in)    :: h    !< Layer thicknesses [H ~> m or kg m-2].
  type(thermo_var_ptrs),   intent(in)    :: tv   !< A structure pointing to various
                                                 !! thermodynamic variables.
  real, dimension(:,:),    pointer       :: p_surf !< A pointer to the surface pressure [R L2 T-2 ~> Pa].
                                                 !! If p_surf is not associated, it is the same
                                                 !! as setting the surface pressure to 0.
  type(diagnostics_CS),    intent(inout) :: CS   !< Control structure returned by a
                                                 !! previous call to diagnostics_init.
  ! Local variables
end subroutine calculate_vertical_integrals
module subroutine calculate_energy_diagnostics(u, v, h, uh, vh, ADp, CDp, G, GV, US, CS)
  type(ocean_grid_type),   intent(inout) :: G    !< The ocean's grid structure.
  type(verticalGrid_type), intent(in)    :: GV   !< The ocean's vertical grid structure.
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                           intent(in)    :: u    !< The zonal velocity [L T-1 ~> m s-1].
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                           intent(in)    :: v    !< The meridional velocity [L T-1 ~> m s-1].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in)    :: h    !< Layer thicknesses [H ~> m or kg m-2].
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                           intent(in)    :: uh   !< Transport through zonal faces=u*h*dy,
                                                 !! [H L2 T-1 ~> m3 s-1 or kg s-1].
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                           intent(in)    :: vh   !< Transport through merid faces=v*h*dx,
                                                 !! [H L2 T-1 ~> m3 s-1 or kg s-1].
  type(accel_diag_ptrs),   intent(in)    :: ADp  !< Structure pointing to accelerations in momentum equation.
  type(cont_diag_ptrs),    intent(in)    :: CDp  !< Structure pointing to terms in continuity equations.
  type(unit_scale_type),   intent(in)    :: US   !< A dimensional unit scaling type
  type(diagnostics_CS),    intent(inout) :: CS   !< Control structure returned by a previous call to
                                                 !! diagnostics_init.

  ! Local variables
                                 ! [H L2 T-3 ~> m3 s-3 or W m-2]
                                 ! [H L4 T-3 ~> m5 s-3 or W]
                                 ! [H L4 T-3 ~> m5 s-3 or W]
                                 ! [H L2 T-3 ~> m3 s-3 or W m-2]

end subroutine calculate_energy_diagnostics
module subroutine register_time_deriv(lb, f_ptr, deriv_ptr, CS)
  integer, intent(in), dimension(3) :: lb     !< Lower index bound of f_ptr
  real, dimension(lb(1):,lb(2):,:), target :: f_ptr
                                              !< Time derivative operand, in arbitrary units [A ~> a]
  real, dimension(lb(1):,lb(2):,:), target :: deriv_ptr
                                              !< Time derivative of f_ptr, in units derived from
                                              !! the arbitrary units of f_ptr [A T-1 ~> a s-1]
  type(diagnostics_CS), intent(inout) :: CS   !< Control structure returned by previous call to
                                              !! diagnostics_init.

  ! This subroutine registers fields to calculate a diagnostic time derivative.
  ! NOTE: Lower bound is required for grid indexing in calculate_derivs().
  !       We assume that the vertical axis is 1-indexed.


end subroutine register_time_deriv
module subroutine calculate_derivs(dt, G, CS)
  real,                  intent(in)    :: dt   !< The time interval over which differences occur [T ~> s].
  type(ocean_grid_type), intent(inout) :: G    !< The ocean's grid structure.
  type(diagnostics_CS),  intent(inout) :: CS   !< Control structure returned by previous call to
                                               !! diagnostics_init.

! This subroutine calculates all registered time derivatives.

end subroutine calculate_derivs
module subroutine post_surface_dyn_diags(IDs, G, diag, sfc_state, ssh)
  type(surface_diag_IDs),   intent(in) :: IDs !< A structure with the diagnostic IDs.
  type(ocean_grid_type),    intent(in) :: G   !< ocean grid structure
  type(diag_ctrl),          intent(in) :: diag !< regulates diagnostic output
  type(surface),            intent(in) :: sfc_state !< structure describing the ocean surface state
  real, dimension(SZI_(G),SZJ_(G)), &
                            intent(in) :: ssh !< Time mean surface height without corrections
                                              !! for ice displacement [Z ~> m]

  ! Local variables

end subroutine post_surface_dyn_diags
module subroutine post_surface_thermo_diags(IDs, G, GV, US, diag, dt_int, sfc_state, tv, &
                                    ssh, ssh_ibc)
  type(surface_diag_IDs),   intent(in) :: IDs !< A structure with the diagnostic IDs.
  type(ocean_grid_type),    intent(in) :: G   !< ocean grid structure
  type(verticalGrid_type),  intent(in) :: GV  !< ocean vertical grid structure
  type(unit_scale_type),    intent(in) :: US  !< A dimensional unit scaling type
  type(diag_ctrl),          intent(in) :: diag  !< regulates diagnostic output
  real,                     intent(in) :: dt_int !< total time step associated with these diagnostics [T ~> s].
  type(surface),            intent(in) :: sfc_state !< structure describing the ocean surface state
  type(thermo_var_ptrs),    intent(in) :: tv  !< A structure pointing to various thermodynamic variables
  real, dimension(SZI_(G),SZJ_(G)), intent(in) :: ssh !< Time mean surface height without corrections
                                              !! for ice displacement [Z ~> m]
  real, dimension(SZI_(G),SZJ_(G)), intent(in) :: ssh_ibc !< Time mean surface height with corrections
                                              !! for ice displacement and the inverse barometer [Z ~> m]


end subroutine post_surface_thermo_diags
module subroutine post_transport_diagnostics(G, GV, US, uhtr, vhtr, h, IDs, diag_pre_dyn, &
                                      diag, dt_trans, Reg)
  type(ocean_grid_type),    intent(inout) :: G   !< ocean grid structure
  type(verticalGrid_type),  intent(in)    :: GV  !< ocean vertical grid structure
  type(unit_scale_type),    intent(in)    :: US  !< A dimensional unit scaling type
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), intent(in) :: uhtr !< Accumulated zonal thickness fluxes
                                                 !! used to advect tracers [H L2 ~> m3 or kg]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), intent(in) :: vhtr !< Accumulated meridional thickness fluxes
                                                 !! used to advect tracers [H L2 ~> m3 or kg]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                            intent(in)    :: h   !< The updated layer thicknesses [H ~> m or kg m-2]
  type(transport_diag_IDs), intent(in)    :: IDs !< A structure with the diagnostic IDs.
  type(diag_grid_storage),  intent(inout) :: diag_pre_dyn !< Stored grids from before dynamics
  type(diag_ctrl),          intent(inout) :: diag !< regulates diagnostic output
  real,                     intent(in)    :: dt_trans !< total time step associated with the transports [T ~> s].
  type(tracer_registry_type), pointer     :: Reg !< Pointer to the tracer registry

  ! Local variables
                          ! [H T-1 ~> m s-1 or kg m-2 s-1].
                          ! [R Z H-1 T-1 ~> kg m-3 s-1 or s-1].
end subroutine post_transport_diagnostics
module subroutine MOM_diagnostics_init(MIS, ADp, CDp, Time, G, GV, US, param_file, diag, CS, tv)
  type(ocean_internal_state), intent(in)    :: MIS  !< For "MOM Internal State" a set of pointers to
                                                    !! the fields and accelerations that make up the
                                                    !! ocean's internal physical state.
  type(accel_diag_ptrs),      intent(inout) :: ADp  !< Structure with pointers to momentum equation
                                                    !! terms.
  type(cont_diag_ptrs),       intent(inout) :: CDp  !< Structure with pointers to continuity
                                                    !! equation terms.
  type(time_type),            intent(in)    :: Time !< Current model time.
  type(ocean_grid_type),      intent(in)    :: G    !< The ocean's grid structure.
  type(verticalGrid_type),    intent(in)    :: GV   !< The ocean's vertical grid structure.
  type(unit_scale_type),      intent(in)    :: US   !< A dimensional unit scaling type
  type(param_file_type),      intent(in)    :: param_file !< A structure to parse for run-time
                                                    !! parameters.
  type(diag_ctrl), target,    intent(inout) :: diag !< Structure to regulate diagnostic output.
  type(diagnostics_CS),       intent(inout) :: CS   !< Diagnostic control struct
  type(thermo_var_ptrs),      intent(in)    :: tv   !< A structure pointing to various
                                                    !! thermodynamic variables.

  ! Local variables
                              ! MKS units (m or kg m-2) for thicknesses depending on whether the
                              ! Boussinesq approximation is being made [m H-1 ~> 1] or [kg m-2 H-1 ~> 1]
                              ! mode wave speed as the starting point for iterations.
  ! This include declares and sets the variable "version".
                                  ! for remapping.  Values below 20190101 recover the remapping
                                  ! answers from 2018, while higher values use more robust
                                  ! forms of the same remapping expressions.

end subroutine MOM_diagnostics_init
module subroutine register_surface_diags(Time, G, US, IDs, diag, tv)
  type(time_type),         intent(in)    :: Time  !< current model time
  type(ocean_grid_type),   intent(in)    :: G     !< ocean grid structure
  type(unit_scale_type),   intent(in)    :: US   !< A dimensional unit scaling type
  type(surface_diag_IDs),  intent(inout) :: IDs   !< A structure with the diagnostic IDs.
  type(diag_ctrl),         intent(inout) :: diag  !< regulates diagnostic output
  type(thermo_var_ptrs),   intent(in)    :: tv    !< A structure pointing to various thermodynamic variables

  ! Vertically integrated, budget, and surface state diagnostics
end subroutine register_surface_diags
module subroutine register_transport_diags(Time, G, GV, US, IDs, diag)
  type(time_type),          intent(in)    :: Time  !< current model time
  type(ocean_grid_type),    intent(in)    :: G     !< ocean grid structure
  type(verticalGrid_type),  intent(in)    :: GV    !< ocean vertical grid structure
  type(unit_scale_type),    intent(in)    :: US    !< A dimensional unit scaling type
  type(transport_diag_IDs), intent(inout) :: IDs   !< A structure with the diagnostic IDs.
  type(diag_ctrl),          intent(inout) :: diag  !< regulates diagnostic output


end subroutine register_transport_diags
module subroutine write_static_fields(G, GV, US, tv, diag)
  type(ocean_grid_type),   intent(in)    :: G    !< ocean grid structure
  type(verticalGrid_type), intent(in)    :: GV   !< ocean vertical grid structure
  type(unit_scale_type),   intent(in)    :: US   !< A dimensional unit scaling type
  type(thermo_var_ptrs),   intent(in)    :: tv   !< A structure pointing to various thermodynamic variables
  type(diag_ctrl), target, intent(inout) :: diag !< regulates diagnostic output

  ! Local variables

end subroutine write_static_fields
module subroutine set_dependent_diagnostics(MIS, ADp, CDp, G, GV, CS)
  type(ocean_internal_state), intent(in)    :: MIS !< For "MOM Internal State" a set of pointers to
                                                   !! the fields and accelerations making up ocean
                                                   !! internal physical state.
  type(accel_diag_ptrs),      intent(inout) :: ADp !< Structure pointing to accelerations in
                                                   !! momentum equation.
  type(cont_diag_ptrs),       intent(inout) :: CDp !< Structure pointing to terms in continuity
                                                   !! equation.
  type(ocean_grid_type),      intent(in)    :: G   !< The ocean's grid structure.
  type(verticalGrid_type),    intent(in)    :: GV   !< ocean vertical grid structure
  type(diagnostics_CS),       intent(inout) :: CS  !< Pointer to the control structure for this
                                                   !! module.

end subroutine set_dependent_diagnostics
module subroutine MOM_diagnostics_end(CS, ADp, CDp)
  type(diagnostics_CS),  intent(inout) :: CS  !< Control structure returned by a
                                              !! previous call to diagnostics_init.
  type(accel_diag_ptrs), intent(inout) :: ADp !< structure with pointers to
                                              !! accelerations in momentum equation.
  type(cont_diag_ptrs),  intent(inout) :: CDp !< Structure pointing to terms in continuity
                                              !! equation.

end subroutine MOM_diagnostics_end
  end interface

end module MOM_diagnostics
