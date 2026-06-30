! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> Simulates CFCs using the OCMIP2 protocols
module MOM_OCMIP2_CFC

use MOM_coms,            only : EFP_type
use MOM_coupler_types,   only : extract_coupler_type_data, set_coupler_type_data
use MOM_coupler_types,   only : atmos_ocn_coupler_flux
use MOM_diag_mediator,   only : diag_ctrl
use MOM_error_handler,   only : MOM_error, FATAL, WARNING
use MOM_file_parser,     only : get_param, log_param, log_version, param_file_type
use MOM_forcing_type,    only : forcing
use MOM_hor_index,       only : hor_index_type
use MOM_grid,            only : ocean_grid_type
use MOM_io,              only : file_exists, MOM_read_data, slasher
use MOM_io,              only : vardesc, var_desc, query_vardesc
use MOM_open_boundary,   only : ocean_OBC_type
use MOM_restart,         only : query_initialized, set_initialized, MOM_restart_CS
use MOM_spatial_means,   only : global_mass_int_EFP
use MOM_sponge,          only : set_up_sponge_field, sponge_CS
use MOM_time_manager,    only : time_type
use MOM_tracer_registry, only : register_tracer, tracer_registry_type
use MOM_tracer_diabatic, only : tracer_vertdiff, applyTracerBoundaryFluxesInOut
use MOM_tracer_Z_init,   only : tracer_Z_init
use MOM_unit_scaling,    only : unit_scale_type
use MOM_variables,       only : surface
use MOM_verticalGrid,    only : verticalGrid_type

implicit none ; private

#include <MOM_memory.h>

public register_OCMIP2_CFC, initialize_OCMIP2_CFC, flux_init_OCMIP2_CFC
public OCMIP2_CFC_column_physics, OCMIP2_CFC_surface_state
public OCMIP2_CFC_stock, OCMIP2_CFC_end

!> The control structure for the  OCMPI2_CFC tracer package
type, public :: OCMIP2_CFC_CS ; private
  character(len=200) :: IC_file !< The file in which the CFC initial values can
                                !! be found, or an empty string for internal initilaization.
  logical :: Z_IC_file !< If true, the IC_file is in Z-space.  The default is false..
  type(time_type), pointer :: Time => NULL() !< A pointer to the ocean model's clock.
  type(tracer_registry_type), pointer :: tr_Reg => NULL() !< A pointer to the MOM6 tracer registry
  real, pointer, dimension(:,:,:) :: &
    CFC11 => NULL(), &     !< The CFC11 concentration [mol m-3].
    CFC12 => NULL()        !< The CFC12 concentration [mol m-3].
  ! In the following variables a suffix of _11 refers to CFC11 and _12 to CFC12.
  !>@{ Coefficients used in the CFC11 and CFC12 solubility calculation
  real :: a1_11, a1_12   ! Coefficients for calculating CFC11 and CFC12 Schmidt numbers [nondim]
  real :: a2_11, a2_12   ! Coefficients for calculating CFC11 and CFC12 Schmidt numbers [degC-1]
  real :: a3_11, a3_12   ! Coefficients for calculating CFC11 and CFC12 Schmidt numbers [degC-2]
  real :: a4_11, a4_12   ! Coefficients for calculating CFC11 and CFC12 Schmidt numbers [degC-3]

  real :: d1_11, d1_12   ! Coefficients for calculating CFC11 and CFC12 solubilities [nondim]
  real :: d2_11, d2_12   ! Coefficients for calculating CFC11 and CFC12 solubilities [hectoKelvin-1]
  real :: d3_11, d3_12   ! Coefficients for calculating CFC11 and CFC12 solubilities [log(hectoKelvin)-1]
  real :: d4_11, d4_12   ! Coefficients for calculating CFC11 and CFC12 solubilities [hectoKelvin-2]

  real :: e1_11, e1_12   ! Coefficients for calculating CFC11 and CFC12 solubilities [PSU-1]
  real :: e2_11, e2_12   ! Coefficients for calculating CFC11 and CFC12 solubilities [PSU-1 hectoKelvin-1]
  real :: e3_11, e3_12   ! Coefficients for calculating CFC11 and CFC12 solubilities [PSU-2 hectoKelvin-2]
  !>@}
  real :: CFC11_IC_val = 0.0    !< The initial value assigned to CFC11 [mol m-3].
  real :: CFC12_IC_val = 0.0    !< The initial value assigned to CFC12 [mol m-3].
  real :: CFC11_land_val = -1.0 !< The value of CFC11 used where land is masked out [mol m-3].
  real :: CFC12_land_val = -1.0 !< The value of CFC12 used where land is masked out [mol m-3].
  logical :: tracers_may_reinit !< If true, tracers may be reset via the initialization code
                                !! if they are not found in the restart files.
  character(len=16) :: CFC11_name !< CFC11 variable name
  character(len=16) :: CFC12_name !< CFC12 variable name

  integer :: ind_cfc_11_flux  !< Index returned by atmos_ocn_coupler_flux that is used to
                              !! pack and unpack surface boundary condition arrays.
  integer :: ind_cfc_12_flux  !< Index returned by atmos_ocn_coupler_flux that is used to
                              !! pack and unpack surface boundary condition arrays.

  type(diag_ctrl), pointer :: diag => NULL() !< A structure that is used to regulate
                                             !! the timing of diagnostic output.
  type(MOM_restart_CS), pointer :: restart_CSp => NULL()  !< Model restart control structure

  ! The following vardesc types contain a package of metadata about each tracer.
  type(vardesc) :: CFC11_desc !< A set of metadata for the CFC11 tracer
  type(vardesc) :: CFC12_desc !< A set of metadata for the CFC12 tracer
end type OCMIP2_CFC_CS


  interface
module function register_OCMIP2_CFC(HI, GV, param_file, CS, tr_Reg, restart_CS)
  type(hor_index_type),    intent(in) :: HI         !< A horizontal index type structure.
  type(verticalGrid_type), intent(in) :: GV         !< The ocean's vertical grid structure.
  type(param_file_type),   intent(in) :: param_file !< A structure to parse for run-time parameters.
  type(OCMIP2_CFC_CS),     pointer    :: CS         !< A pointer that is set to point to the control
                                                    !! structure for this module.
  type(tracer_registry_type), &
                           pointer    :: tr_Reg     !< A pointer to the tracer registry.
  type(MOM_restart_CS), target, intent(inout) :: restart_CS !< MOM restart control struct

  ! Local variables
  ! This include declares and sets the variable "version".
  logical :: register_OCMIP2_CFC

end function register_OCMIP2_CFC
module subroutine flux_init_OCMIP2_CFC(CS, verbosity)
  type(OCMIP2_CFC_CS), optional, pointer :: CS !< An optional pointer to the control structure
                                               !! for this module; if not present, the flux indicies
                                               !! are not stored.
  integer,             optional, intent(in) :: verbosity !< A 0-9 integer indicating a level of verbosity.

  ! These can be overridden later in via the field manager?

  ! These calls obtain the indices for the CFC11 and CFC12 flux coupling.  They
  ! can safely be called multiple times.
end subroutine flux_init_OCMIP2_CFC
module subroutine initialize_OCMIP2_CFC(restart, day, G, GV, US, h, diag, OBC, CS, &
                                 sponge_CSp)
  logical,                        intent(in) :: restart    !< .true. if the fields have already been
                                                           !! read from a restart file.
  type(time_type), target,        intent(in) :: day        !< Time of the start of the run.
  type(ocean_grid_type),          intent(in) :: G          !< The ocean's grid structure.
  type(verticalGrid_type),        intent(in) :: GV         !< The ocean's vertical grid structure.
  type(unit_scale_type),          intent(in) :: US         !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                                  intent(in) :: h          !< Layer thicknesses [H ~> m or kg m-2].
  type(diag_ctrl), target,        intent(in) :: diag       !< A structure that is used to regulate
                                                           !! diagnostic output.
  type(ocean_OBC_type),           pointer    :: OBC        !< This open boundary condition type
                                                           !! specifies whether, where, and what
                                                           !! open boundary conditions are used.
  type(OCMIP2_CFC_CS),            pointer    :: CS         !< The control structure returned by a
                                                           !! previous call to register_OCMIP2_CFC.
  type(sponge_CS),                pointer    :: sponge_CSp !< A pointer to the control structure for
                                                           !! the sponges, if they are in use.
                                                           !! Otherwise this may be unassociated.

end subroutine initialize_OCMIP2_CFC
module subroutine init_tracer_CFC(h, tr, name, land_val, IC_val, G, GV, US, CS)
  type(ocean_grid_type),                     intent(in)  :: G    !< The ocean's grid structure
  type(verticalGrid_type),                   intent(in)  :: GV   !< The ocean's vertical grid structure.
  type(unit_scale_type),                     intent(in)  :: US   !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in)  :: h    !< Layer thicknesses [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(out) :: tr   !< The CFC tracer concentration array [mol m-3]
  character(len=*),                          intent(in)  :: name !< The tracer name
  real,                                      intent(in)  :: land_val !< A value the tracer takes over land [mol m-3]
  real,                                      intent(in)  :: IC_val !< The initial condition value for
                                                                 !! the CRC tracer [mol m-3]
  type(OCMIP2_CFC_CS),                       pointer     :: CS   !< The control structure returned by a
                                                                 !! previous call to register_OCMIP2_CFC.

  ! This subroutine initializes a tracer array.

end subroutine init_tracer_CFC
module subroutine OCMIP2_CFC_column_physics(h_old, h_new, ea, eb, fluxes, dt, G, GV, US, CS, &
              evap_CFL_limit, minimum_forcing_depth)
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
  type(OCMIP2_CFC_CS),     pointer    :: CS   !< The control structure returned by a
                                              !! previous call to register_OCMIP2_CFC.
  real,          optional, intent(in) :: evap_CFL_limit !< Limit on the fraction of the water that can
                                              !! be fluxed out of the top layer in a timestep [nondim]
  real,          optional, intent(in) :: minimum_forcing_depth !< The smallest depth over which
                                              !! fluxes can be applied [H ~> m or kg m-2]
!   This subroutine applies diapycnal diffusion and any other column
! tracer physics or chemistry to the tracers from this file.
! CFCs are relatively simple, as they are passive tracers. with only a surface
! flux as a source.

! The arguments to this subroutine are redundant in that
!     h_new(k) = h_old(k) + ea(k) - eb(k-1) + eb(k) - ea(k+1)

  ! Local variables

end subroutine OCMIP2_CFC_column_physics
module function OCMIP2_CFC_stock(h, stocks, G, GV, CS, names, units, stock_index)
  type(ocean_grid_type),           intent(in)    :: G      !< The ocean's grid structure.
  type(verticalGrid_type),         intent(in)    :: GV     !< The ocean's vertical grid structure.
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                                   intent(in)    :: h      !< Layer thicknesses [H ~> m or kg m-2].
  type(EFP_type), dimension(:),    intent(out)   :: stocks !< The mass-weighted integrated amount of each
                                                           !! tracer, in kg times concentration units [kg conc]
  type(OCMIP2_CFC_CS),             pointer       :: CS     !< The control structure returned by a
                                                           !! previous call to register_OCMIP2_CFC.
  character(len=*), dimension(:),  intent(out)   :: names  !< The names of the stocks calculated.
  character(len=*), dimension(:),  intent(out)   :: units  !< The units of the stocks calculated.
  integer, optional,               intent(in)    :: stock_index !< The coded index of a specific
                                                                !! stock being sought.
  integer                                        :: OCMIP2_CFC_stock !< The number of stocks calculated here.


end function OCMIP2_CFC_stock
module subroutine OCMIP2_CFC_surface_state(sfc_state, h, G, GV, US, CS)
  type(ocean_grid_type),   intent(in)    :: G  !< The ocean's grid structure.
  type(verticalGrid_type), intent(in)    :: GV !< The ocean's vertical grid structure
  type(surface),           intent(inout) :: sfc_state !< A structure containing fields that
                                               !! describe the surface state of the ocean.
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in)    :: h  !< Layer thickness [H ~> m or kg m-2].
  type(unit_scale_type),   intent(in)    :: US !< A dimensional unit scaling type
  type(OCMIP2_CFC_CS),     pointer       :: CS !< The control structure returned by a previous
                                               !! call to register_OCMIP2_CFC.

  ! Local variables

end subroutine OCMIP2_CFC_surface_state
module subroutine OCMIP2_CFC_end(CS)
  type(OCMIP2_CFC_CS), pointer :: CS   !< The control structure returned by a
                                       !! previous call to register_OCMIP2_CFC.
!   This subroutine deallocates the memory owned by this module.
! Argument: CS - The control structure returned by a previous call to
!                register_OCMIP2_CFC.

end subroutine OCMIP2_CFC_end
  end interface

end module MOM_OCMIP2_CFC
