! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> Implements a boundary impulse response tracer to calculate Green's functions
module boundary_impulse_tracer

use MOM_coms,            only : EFP_type
use MOM_coupler_types,   only : set_coupler_type_data, atmos_ocn_coupler_flux
use MOM_diag_mediator,   only : diag_ctrl
use MOM_error_handler,   only : MOM_error, FATAL, WARNING
use MOM_file_parser,     only : get_param, log_param, log_version, param_file_type
use MOM_forcing_type,    only : forcing
use MOM_grid,            only : ocean_grid_type
use MOM_hor_index,       only : hor_index_type
use MOM_io,              only : vardesc, var_desc, query_vardesc
use MOM_open_boundary,   only : ocean_OBC_type
use MOM_restart,         only : register_restart_field, query_initialized, set_initialized, MOM_restart_CS
use MOM_spatial_means,   only : global_mass_int_EFP
use MOM_sponge,          only : set_up_sponge_field, sponge_CS
use MOM_time_manager,    only : time_type
use MOM_tracer_registry, only : register_tracer, tracer_registry_type
use MOM_tracer_diabatic, only : tracer_vertdiff, applyTracerBoundaryFluxesInOut
use MOM_tracer_Z_init,   only : tracer_Z_init
use MOM_unit_scaling,    only : unit_scale_type
use MOM_variables,       only : surface, thermo_var_ptrs
use MOM_verticalGrid,    only : verticalGrid_type

implicit none ; private

#include <MOM_memory.h>

public register_boundary_impulse_tracer, initialize_boundary_impulse_tracer
public boundary_impulse_tracer_column_physics, boundary_impulse_tracer_surface_state
public boundary_impulse_stock, boundary_impulse_tracer_end

!> NTR_MAX is the maximum number of tracers in this module.
integer, parameter :: NTR_MAX = 1

!> The control structure for the boundary impulse tracer package
type, public :: boundary_impulse_tracer_CS ; private
  integer :: ntr=NTR_MAX    !< The number of tracers that are actually used.
  logical :: coupled_tracers = .false. !< These tracers are not offered to the  coupler.
  type(time_type), pointer :: Time => NULL() !< A pointer to the ocean model's clock.
  type(tracer_registry_type), pointer :: tr_Reg => NULL() !< A pointer to the tracer registry
  real, pointer :: tr(:,:,:,:) => NULL() !< The array of tracers used in this subroutine, in [CU ~> conc] (g m-3)?
  logical :: tracers_may_reinit  !< If true, boundary_impulse can be initialized if not found in restart file
  integer, dimension(NTR_MAX) :: ind_tr  !< Indices returned by atmos_ocn_coupler_flux if it is used and the
                                         !! surface tracer concentrations are to be provided to the coupler.

  integer :: nkml !< Number of layers in mixed layer
  real, dimension(NTR_MAX)  :: land_val = -1.0 !< A value to use to fill in tracers over land [CU ~> conc]
  real :: remaining_source_time !< How much longer (same units as the timestep) to
                                !! inject the tracer at the surface [T ~> s]

  type(diag_ctrl), pointer :: diag => NULL() !< A structure that is used to
                                   !! regulate the timing of diagnostic output.
  type(MOM_restart_CS), pointer :: restart_CSp => NULL() !< A pointer to the retart control structure

  type(vardesc) :: tr_desc(NTR_MAX) !< Descriptions and metadata for the tracers
end type boundary_impulse_tracer_CS


  interface
module function register_boundary_impulse_tracer(HI, GV, US, param_file, CS, tr_Reg, restart_CS)
  type(hor_index_type),             intent(in   ) :: HI   !< A horizontal index type structure
  type(verticalGrid_type),          intent(in   ) :: GV   !< The ocean's vertical grid structure
  type(unit_scale_type),            intent(in   ) :: US   !< A dimensional unit scaling type
  type(param_file_type),            intent(in   ) :: param_file !< A structure to parse for run-time parameters
  type(boundary_impulse_tracer_CS), pointer       :: CS   !< The control structure returned by a previous
                                                          !! call to register_boundary_impulse_tracer.
  type(tracer_registry_type),       pointer       :: tr_Reg !< A pointer that is set to point to the control
                                                          !! structure for the tracer advection and
                                                          !! diffusion module
  type(MOM_restart_CS), target, intent(inout) :: restart_CS !< MOM restart control struct

  ! Local variables
                            ! kg(tracer) kg(water)-1 m3 s-1 or kg(tracer) s-1.
  ! This include declares and sets the variable "version".
  logical :: register_boundary_impulse_tracer
end function register_boundary_impulse_tracer
module subroutine initialize_boundary_impulse_tracer(restart, day, G, GV, US, h, diag, OBC, CS, &
                                  sponge_CSp, tv)
  logical,                            intent(in) :: restart !< .true. if the fields have already
                                                         !! been read from a restart file.
  type(time_type),            target, intent(in) :: day  !< Time of the start of the run.
  type(ocean_grid_type),              intent(in) :: G    !< The ocean's grid structure
  type(verticalGrid_type),            intent(in) :: GV   !< The ocean's vertical grid structure
  type(unit_scale_type),              intent(in) :: US   !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                                      intent(in) :: h    !< Layer thicknesses [H ~> m or kg m-2]
  type(diag_ctrl),            target, intent(in) :: diag !< A structure that is used to regulate
                                                         !! diagnostic output.
  type(ocean_OBC_type),               pointer    :: OBC  !< This open boundary condition type specifies
                                                         !! whether, where, and what open boundary
                                                         !! conditions are used.
  type(boundary_impulse_tracer_CS),   pointer    :: CS   !< The control structure returned by a previous
                                                         !! call to register_boundary_impulse_tracer.
  type(sponge_CS),                    pointer    :: sponge_CSp !< Pointer to the control structure for the sponges.
  type(thermo_var_ptrs),              intent(in) :: tv   !< A structure pointing to various
                                                         !! thermodynamic variables
  ! Local variables

end subroutine initialize_boundary_impulse_tracer
module subroutine boundary_impulse_tracer_column_physics(h_old, h_new, ea, eb, fluxes, dt, G, GV, US, CS, &
                     tv, debug, evap_CFL_limit, minimum_forcing_depth)
  type(ocean_grid_type),   intent(in) :: G    !< The ocean's grid structure
  type(verticalGrid_type), intent(in) :: GV   !< The ocean's vertical grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in) :: h_old !< Layer thickness before entrainment [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in) :: h_new !< Layer thickness after entrainment [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in) :: ea   !< an array to which the amount of fluid entrained
                                              !! from the layer above during this call will be
                                              !! added [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in) :: eb   !< an array to which the amount of fluid entrained
                                              !! from the layer below during this call will be
                                              !! added [H ~> m or kg m-2].
  type(forcing),           intent(in) :: fluxes !< A structure containing pointers to thermodynamic
                                              !! and tracer forcing fields.  Unused fields have NULL ptrs.
  real,                    intent(in) :: dt   !< The amount of time covered by this call [T ~> s]
  type(unit_scale_type),   intent(in) :: US   !< A dimensional unit scaling type
  type(boundary_impulse_tracer_CS),  pointer :: CS !< The control structure returned by a previous
                                              !! call to register_boundary_impulse_tracer.
  type(thermo_var_ptrs),   intent(in) :: tv   !< A structure pointing to various
                                              !! thermodynamic variables
  logical,                 intent(in) :: debug !< If true calculate checksums
  real,          optional, intent(in) :: evap_CFL_limit !< Limit on the fraction of the water that can
                                              !! be fluxed out of the top layer in a timestep [nondim]
  real,          optional, intent(in) :: minimum_forcing_depth !< The smallest depth over which
                                              !! fluxes can be applied [H ~> m or kg m-2]

!   This subroutine applies diapycnal diffusion and any other column
! tracer physics or chemistry to the tracers from this file.
! This is a simple example of a set of advected passive tracers.

! The arguments to this subroutine are redundant in that
!     h_new(k) = h_old(k) + ea(k) - eb(k-1) + eb(k) - ea(k+1)

  ! Local variables

end subroutine boundary_impulse_tracer_column_physics
module function boundary_impulse_stock(h, stocks, G, GV, CS, names, units, stock_index)
  type(ocean_grid_type),                    intent(in   ) :: G    !< The ocean's grid structure
  type(verticalGrid_type),                  intent(in   ) :: GV   !< The ocean's vertical grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in   ) :: h    !< Layer thicknesses [H ~> m or kg m-2]
  type(EFP_type), dimension(:),             intent(  out) :: stocks !< The mass-weighted integrated amount of each
                                                                  !! tracer, in kg times concentration units [kg conc]
  type(boundary_impulse_tracer_CS),         pointer       :: CS   !< The control structure returned by a previous
                                                                  !! call to register_boundary_impulse_tracer.
  character(len=*), dimension(:),           intent(  out) :: names  !< The names of the stocks calculated.
  character(len=*), dimension(:),           intent(  out) :: units  !< The units of the stocks calculated.
  integer, optional,                        intent(in   ) :: stock_index !< The coded index of a specific stock
                                                                  !! being sought.
  integer :: boundary_impulse_stock  !< Return value: the number of stocks calculated here.

  ! Local variables

end function boundary_impulse_stock
module subroutine boundary_impulse_tracer_surface_state(sfc_state, h, G, GV, CS)
  type(ocean_grid_type),   intent(in)    :: G  !< The ocean's grid structure.
  type(verticalGrid_type), intent(in)    :: GV !< The ocean's vertical grid structure
  type(surface),           intent(inout) :: sfc_state !< A structure containing fields that
                                               !! describe the surface state of the ocean.
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in)    :: h  !< Layer thickness [H ~> m or kg m-2].
  type(boundary_impulse_tracer_CS), pointer :: CS !< The control structure returned by a previous
                                               !! call to register_boundary_impulse_tracer.

  ! This particular tracer package does not report anything back to the coupler.
  ! The code that is here is just a rough guide for packages that would.

end subroutine boundary_impulse_tracer_surface_state
module subroutine boundary_impulse_tracer_end(CS)
  type(boundary_impulse_tracer_CS), pointer :: CS   !< The control structure returned by a previous
                                                    !! call to register_boundary_impulse_tracer.
end subroutine boundary_impulse_tracer_end
  end interface

end module boundary_impulse_tracer
