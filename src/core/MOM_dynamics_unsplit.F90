! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> Time steps the ocean dynamics with an unsplit quasi 3rd order scheme
module MOM_dynamics_unsplit

!********+*********+*********+*********+*********+*********+*********+**
!*                                                                     *
!*  By Robert Hallberg, 1993-2012                                      *
!*                                                                     *
!*    This file contains code that does the time-stepping of the       *
!*  adiabatic dynamic core, in this case with an unsplit third-order   *
!*  Runge-Kutta time stepping scheme for the momentum and a forward-   *
!*  backward coupling between the momentum and continuity equations.   *
!*  This was the orignal unsplit time stepping scheme used in early    *
!*  versions of HIM and its precursor.  While it is very simple and    *
!*  accurate, it is much less efficient that the split time stepping   *
!*  scheme for realistic oceanographic applications.  It has been      *
!*  retained for all of these years primarily to verify that the split *
!*  scheme is giving the right answers, and to debug the failings of   *
!*  the split scheme when it is not.  The split time stepping scheme   *
!*  is now sufficiently robust that it should be first choice for      *
!*  almost any conceivable application, except perhaps from cases      *
!*  with just a few layers for which the exact timing of the high-     *
!*  frequency barotropic gravity waves is of paramount importance.     *
!*  This scheme is slightly more efficient than the other unsplit      *
!*  scheme that can be found in MOM_dynamics_unsplit_RK2.F90.          *
!*                                                                     *
!*    The subroutine step_MOM_dyn_unsplit actually does the time       *
!*  stepping, while register_restarts_dyn_unsplit  sets the fields     *
!*  that are found in a full restart file with this scheme, and        *
!*  initialize_dyn_unsplit  initializes the cpu clocks that are        *                                      *
!*  used in this module.  For largely historical reasons, this module  *
!*  does not have its own control structure, but shares the same       *
!*  control structure with MOM.F90 and the other MOM_dynamics_...      *
!*  modules.                                                           *
!*                                                                     *
!*  Macros written all in capital letters are defined in MOM_memory.h. *
!*                                                                     *
!*     A small fragment of the grid is shown below:                    *
!*                                                                     *
!*    j+1  x ^ x ^ x   At x:  q, CoriolisBu                            *
!*    j+1  > o > o >   At ^:  v, PFv, CAv, vh, diffv, tauy, vbt, vhtr  *
!*    j    x ^ x ^ x   At >:  u, PFu, CAu, uh, diffu, taux, ubt, uhtr  *
!*    j    > o > o >   At o:  h, bathyT, eta, T, S, tr                 *
!*    j-1  x ^ x ^ x                                                   *
!*        i-1  i  i+1                                                  *
!*           i  i+1                                                    *
!*                                                                     *
!*  The boundaries always run through q grid points (x).               *
!*                                                                     *
!********+*********+*********+*********+*********+*********+*********+**

use MOM_variables, only : vertvisc_type, thermo_var_ptrs, porous_barrier_type
use MOM_variables, only : accel_diag_ptrs, ocean_internal_state, cont_diag_ptrs
use MOM_forcing_type, only : mech_forcing
use MOM_checksum_packages, only : MOM_thermo_chksum, MOM_state_chksum, MOM_accel_chksum
use MOM_cpu_clock, only : cpu_clock_id, cpu_clock_begin, cpu_clock_end
use MOM_cpu_clock, only : CLOCK_COMPONENT, CLOCK_SUBCOMPONENT
use MOM_cpu_clock, only : CLOCK_MODULE_DRIVER, CLOCK_MODULE, CLOCK_ROUTINE
use MOM_diag_mediator, only : diag_mediator_init, enable_averages
use MOM_diag_mediator, only : disable_averaging, post_data, safe_alloc_ptr
use MOM_diag_mediator, only : register_diag_field, register_static_field
use MOM_diag_mediator, only : set_diag_mediator_grid, diag_ctrl, diag_update_remap_grids
use MOM_domains, only : pass_var, pass_var_start, pass_var_complete
use MOM_domains, only : pass_vector, pass_vector_start, pass_vector_complete
use MOM_domains, only : To_South, To_West, To_All, CGRID_NE, SCALAR_PAIR
use MOM_error_handler, only : MOM_error, MOM_mesg, FATAL, WARNING, is_root_pe
use MOM_file_parser, only : get_param, log_param, log_version, param_file_type
use MOM_get_input, only : directories
use MOM_time_manager, only : time_type, real_to_time, operator(+)
use MOM_time_manager, only : operator(-), operator(>), operator(*), operator(/)

use MOM_ALE, only : ALE_CS
use MOM_barotropic, only : barotropic_CS
use MOM_boundary_update, only : update_OBC_data, update_OBC_CS
use MOM_continuity, only : continuity, continuity_init, continuity_CS, continuity_stencil
use MOM_CoriolisAdv, only : CorAdCalc, CoriolisAdv_init, CoriolisAdv_CS, CoriolisAdv_stencil
use MOM_debugging, only : check_redundant
use MOM_grid, only : ocean_grid_type
use MOM_hor_index, only : hor_index_type
use MOM_hor_visc, only : horizontal_viscosity, hor_visc_init, hor_visc_CS
use MOM_interface_heights, only : find_eta, thickness_to_dz
use MOM_lateral_mixing_coeffs, only : VarMix_CS
use MOM_MEKE_types, only : MEKE_type
use MOM_open_boundary, only : ocean_OBC_type
use MOM_open_boundary, only : radiation_open_bdry_conds
use MOM_open_boundary, only : open_boundary_zero_normal_flow
use MOM_PressureForce, only : PressureForce, PressureForce_init, PressureForce_CS
use MOM_set_visc, only : set_viscous_ML, set_visc_CS
use MOM_stochastics,   only : stochastic_CS
use MOM_tidal_forcing, only : tidal_forcing_init, tidal_forcing_end, tidal_forcing_CS
use MOM_self_attr_load, only : SAL_init, SAL_end, SAL_CS
use MOM_unit_scaling,  only : unit_scale_type
use MOM_vert_friction, only : vertvisc, vertvisc_coef, vertvisc_init, vertvisc_CS
use MOM_verticalGrid, only : verticalGrid_type, get_thickness_units
use MOM_verticalGrid, only : get_flux_units, get_tr_flux_units
use MOM_wave_interface, only: wave_parameters_CS

implicit none ; private

#include <MOM_memory.h>

!> MOM_dynamics_unsplit module control structure
type, public :: MOM_dyn_unsplit_CS ; private
  real ALLOCABLE_, dimension(NIMEMB_PTR_,NJMEM_,NKMEM_) :: &
    CAu, &    !< CAu = f*v - u.grad(u) [L T-2 ~> m s-2].
    PFu, &    !< PFu = -dM/dx [L T-2 ~> m s-2].
    diffu     !< Zonal acceleration due to convergence of the along-isopycnal stress tensor [L T-2 ~> m s-2].

  real ALLOCABLE_, dimension(NIMEM_,NJMEMB_PTR_,NKMEM_) :: &
    CAv, &    !< CAv = -f*u - u.grad(v) [L T-2 ~> m s-2].
    PFv, &    !< PFv = -dM/dy [L T-2 ~> m s-2].
    diffv     !< Meridional acceleration due to convergence of the along-isopycnal stress tensor [L T-2 ~> m s-2].

  real, pointer, dimension(:,:) :: taux_bot => NULL() !< frictional x-bottom stress from the ocean
                                                      !! to the seafloor [R L Z T-2 ~> Pa]
  real, pointer, dimension(:,:) :: tauy_bot => NULL() !< frictional y-bottom stress from the ocean
                                                      !! to the seafloor [R L Z T-2 ~> Pa]

  logical :: dt_visc_bug    !< If false, use the correct timestep in viscous terms applied in the
                            !! first predictor step and in the calculation of the turbulent mixed
                            !! layer properties for viscosity.  If this is true, an older incorrect
                            !! setting is used.
  logical :: debug          !< If true, write verbose checksums for debugging purposes.
  logical :: calculate_SAL  !< If true, calculate self-attraction and loading.
  logical :: use_tides      !< If true, tidal forcing is enabled.

  logical :: module_is_initialized = .false. !< Record whether this module has been initialized.

  !>@{ Diagnostic IDs
  integer :: id_uh = -1, id_vh = -1
  integer :: id_ueffA = -1, id_veffA = -1
  integer :: id_PFu = -1, id_PFv = -1, id_CAu = -1, id_CAv = -1
  !>@}

  type(diag_ctrl), pointer :: diag => NULL() !< A structure that is used to
                                   !! regulate the timing of diagnostic output.
  type(accel_diag_ptrs), pointer :: ADp => NULL() !< A structure pointing to the
                                   !! accelerations in the momentum equations,
                                   !! which can later be used to calculate
                                   !! derived diagnostics like energy budgets.
  type(cont_diag_ptrs), pointer :: CDp => NULL() !< A structure with pointers to
                                   !! various terms in the continuity equations,
                                   !! which can later be used to calculate
                                   !! derived diagnostics like energy budgets.

  ! The remainder of the structure points to child subroutines' control structures.
  !> A pointer to the horizontal viscosity control structure
  type(hor_visc_CS) :: hor_visc
  !> A pointer to the continuity control structure
  type(continuity_CS) :: continuity_CSp
  !> A pointer to the CoriolisAdv control structure
  type(CoriolisAdv_CS) :: CoriolisAdv
  !> A pointer to the PressureForce control structure
  type(PressureForce_CS) :: PressureForce_CSp
  !> A pointer to the vertvisc control structure
  type(vertvisc_CS), pointer :: vertvisc_CSp => NULL()
  !> A pointer to the set_visc control structure
  type(set_visc_CS), pointer :: set_visc_CSp => NULL()
  !> A pointer to the SAL control structure
  type(SAL_CS) :: SAL_CSp
  !> A pointer to the tidal forcing control structure
  type(tidal_forcing_CS) :: tides_CSp
  !> A pointer to the ALE control structure.
  type(ALE_CS), pointer :: ALE_CSp => NULL()

  type(ocean_OBC_type), pointer :: OBC => NULL() !< A pointer to an open boundary
     ! condition type that specifies whether, where, and  what open boundary
     ! conditions are used.  If no open BCs are used, this pointer stays
     ! nullified.  Flather OBCs use open boundary_CS as well.
  !> A pointer to the update_OBC control structure
  type(update_OBC_CS),    pointer :: update_OBC_CSp => NULL()

end type MOM_dyn_unsplit_CS

public step_MOM_dyn_unsplit, register_restarts_dyn_unsplit
public initialize_dyn_unsplit, end_dyn_unsplit

!>@{ CPU time clock IDs
integer :: id_clock_Cor, id_clock_pres, id_clock_vertvisc
integer :: id_clock_continuity, id_clock_horvisc, id_clock_mom_update
integer :: id_clock_pass, id_clock_pass_init
!>@}


  interface
module subroutine step_MOM_dyn_unsplit(u, v, h, tv, visc, Time_local, dt, forces, &
                  p_surf_begin, p_surf_end, uh, vh, uhtr, vhtr, eta_av, G, GV, US, CS, &
                  VarMix, MEKE, pbv, STOCH, Waves)
  type(ocean_grid_type),   intent(inout) :: G      !< The ocean's grid structure.
  type(verticalGrid_type), intent(in)    :: GV     !< The ocean's vertical grid structure.
  type(unit_scale_type),   intent(in)    :: US     !< A dimensional unit scaling type
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), intent(inout) :: u !< The zonal velocity [L T-1 ~> m s-1].
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), intent(inout) :: v !< The meridional velocity [L T-1 ~> m s-1].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(inout) :: h !< Layer thicknesses [H ~> m or kg m-2].
  type(thermo_var_ptrs),   intent(in)    :: tv     !< A structure pointing to various
                                                   !! thermodynamic variables.
  type(vertvisc_type),     intent(inout) :: visc   !< A structure containing vertical
                                 !! viscosities, bottom drag viscosities, and related fields.
  type(time_type),         intent(in)    :: Time_local   !< The model time at the end
                                                         !! of the time step.
  real,                    intent(in)    :: dt     !< The dynamics time step [T ~> s].
  type(mech_forcing),      intent(in)    :: forces !< A structure with the driving mechanical forces
  real, dimension(:,:),    pointer       :: p_surf_begin !< A pointer (perhaps NULL) to the surface
                                                   !! pressure at the start of this dynamic step [R L2 T-2 ~> Pa].
  real, dimension(:,:),    pointer       :: p_surf_end   !< A pointer (perhaps NULL) to the surface
                                                   !! pressure at the end of this dynamic step [R L2 T-2 ~> Pa].
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), intent(inout) :: uh !< The zonal volume or mass transport
                                                   !! [H L2 T-1 ~> m3 s-1 or kg s-1].
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), intent(inout) :: vh !< The meridional volume or mass
                                                   !! transport [H L2 T-1 ~> m3 s-1 or kg s-1].
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), intent(inout) :: uhtr !< The accumulated zonal volume or mass
                                                   !! transport since the last tracer advection [H L2 ~> m3 or kg].
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), intent(inout) :: vhtr !< The accumulated meridional volume or mass
                                                   !! transport since the last tracer advection [H L2 ~> m3 or kg].
  real, dimension(SZI_(G),SZJ_(G)), intent(out) :: eta_av !< The time-mean free surface height or
                                                   !! column mass [H ~> m or kg m-2].
  type(MOM_dyn_unsplit_CS), pointer      :: CS     !< The control structure set up by
                                                   !! initialize_dyn_unsplit.
  type(VarMix_CS),         intent(inout) :: VarMix !< Variable mixing control structure
  type(MEKE_type),         intent(inout) :: MEKE   !< MEKE fields
  type(porous_barrier_type), intent(in) :: pbv     !< porous barrier fractional cell metrics
  type(stochastic_CS),   intent(inout) :: STOCH    !< Stochastic control structure
  type(wave_parameters_CS), optional, pointer :: Waves !< A pointer to a structure containing
                                 !! fields related to the surface wave conditions

  ! Local variables
end subroutine step_MOM_dyn_unsplit
module subroutine register_restarts_dyn_unsplit(HI, GV, param_file, CS)
  type(hor_index_type),      intent(in) :: HI         !< A horizontal index type structure.
  type(verticalGrid_type),   intent(in) :: GV         !< The ocean's vertical grid structure.
  type(param_file_type),     intent(in) :: param_file !< A structure to parse for
                                                      !! run-time parameters.
  type(MOM_dyn_unsplit_CS),  pointer    :: CS         !< The control structure set up by
                                                      !! initialize_dyn_unsplit.

end subroutine register_restarts_dyn_unsplit
module subroutine initialize_dyn_unsplit(u, v, h, tv, Time, G, GV, US, param_file, diag, CS, &
                                  Accel_diag, Cont_diag, MIS, &
                                  OBC, update_OBC_CSp, ALE_CSp, set_visc, &
                                  visc, dirs, ntrunc, cont_stencil, dyn_h_stencil)
  type(ocean_grid_type),          intent(inout) :: G          !< The ocean's grid structure.
  type(verticalGrid_type),        intent(in)    :: GV         !< The ocean's vertical grid structure.
  type(unit_scale_type),          intent(in)    :: US         !< A dimensional unit scaling type
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                                  intent(inout) :: u          !< The zonal velocity [L T-1 ~> m s-1].
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                                  intent(inout) :: v          !< The meridional velocity [L T-1 ~> m s-1].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                                  intent(inout) :: h          !< Layer thicknesses [H ~> m or kg m-2]
  type(thermo_var_ptrs),          intent(in)    :: tv         !< Thermodynamic type
  type(time_type),        target, intent(in)    :: Time       !< The current model time.
  type(param_file_type),          intent(in)    :: param_file !< A structure to parse
                                                              !! for run-time parameters.
  type(diag_ctrl),        target, intent(inout) :: diag       !< A structure that is used to
                                                              !! regulate diagnostic output.
  type(MOM_dyn_unsplit_CS),       pointer       :: CS         !< The control structure set up
                                                              !! by initialize_dyn_unsplit.
  type(accel_diag_ptrs),  target, intent(inout) :: Accel_diag !< A set of pointers to the various
                                     !! accelerations in the momentum equations, which can be used
                                     !! for later derived diagnostics, like energy budgets.
  type(cont_diag_ptrs),   target, intent(inout) :: Cont_diag  !< A structure with pointers to
                                                              !! various terms in the continuity
                                                              !! equations.
  type(ocean_internal_state),     intent(inout) :: MIS        !< The "MOM6 Internal State"
                                                   !! structure, used to pass around pointers
                                                   !! to various arrays for diagnostic purposes.
  type(ocean_OBC_type),           pointer       :: OBC        !< If open boundary conditions are
                                                       !! used, this points to the ocean_OBC_type
                                                       !! that was set up in MOM_initialization.
  type(update_OBC_CS),            pointer       :: update_OBC_CSp !< If open boundary condition
                                                            !! updates are used, this points to
                                                            !! the appropriate control structure.
  type(ALE_CS),                   pointer       :: ALE_CSp    !< This points to the ALE control
                                                              !! structure.
  type(set_visc_CS),      target, intent(in)    :: set_visc   !< set_visc control structure
  type(vertvisc_type),            intent(inout) :: visc       !< A structure containing vertical
                                                              !! viscosities, bottom drag
                                                              !! viscosities, and related fields.
  type(directories),              intent(in)    :: dirs       !< A structure containing several
                                                              !! relevant directory paths.
  integer, target,                intent(inout) :: ntrunc     !< A target for the variable that
                                                        !! records the number of times the velocity
                                                        !! is truncated (this should be 0).
  integer,                        intent(out)   :: cont_stencil !< The stencil for thickness
                                                              !! from the continuity solver.
  integer,                        intent(out)   :: dyn_h_stencil !< The stencil for thickness
                                                              !! for the dynamics based on the
                                                              !! continuity solver and Coriolis scheme.

  !   This subroutine initializes all of the variables that are used by this
  ! dynamic core, including diagnostics and the cpu clocks.

  ! Local variables
  ! This include declares and sets the variable "version".
end subroutine initialize_dyn_unsplit
module subroutine end_dyn_unsplit(CS)
  type(MOM_dyn_unsplit_CS), pointer :: CS !< unsplit dynamics control structure that
                                          !! will be deallocated in this subroutine.

end subroutine end_dyn_unsplit
  end interface

end module MOM_dynamics_unsplit
