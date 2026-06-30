! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

  !> Simulates CFCs using atmospheric pressure, wind speed and sea ice cover
!! provided via cap (only NUOPC cap is implemented so far).
module MOM_CFC_cap

use MOM_coms,            only : EFP_type
use MOM_debugging,       only : hchksum
use MOM_diag_mediator,   only : diag_ctrl, register_diag_field, post_data
use MOM_error_handler,   only : MOM_error, FATAL, WARNING
use MOM_file_parser,     only : get_param, log_param, log_version, param_file_type
use MOM_forcing_type,    only : forcing
use MOM_hor_index,       only : hor_index_type
use MOM_grid,            only : ocean_grid_type
use MOM_CVMix_KPP,       only : KPP_NonLocalTransport, KPP_CS
use MOM_io,              only : file_exists, MOM_read_data, slasher
use MOM_io,              only : vardesc, var_desc, query_vardesc, stdout
use MOM_tracer_registry, only : tracer_type
use MOM_open_boundary,   only : ocean_OBC_type
use MOM_restart,         only : query_initialized, set_initialized, MOM_restart_CS
use MOM_spatial_means,   only : global_mass_int_EFP
use MOM_time_manager,    only : time_type, increment_date
use MOM_interpolate,     only : external_field, init_external_field, time_interp_external
use MOM_tracer_registry, only : register_tracer
use MOM_tracer_types,    only : tracer_registry_type
use MOM_tracer_diabatic, only : tracer_vertdiff, applyTracerBoundaryFluxesInOut
use MOM_tracer_Z_init,   only : tracer_Z_init
use MOM_unit_scaling,    only : unit_scale_type
use MOM_variables,       only : surface, thermo_var_ptrs
use MOM_verticalGrid,    only : verticalGrid_type

implicit none ; private

#include <MOM_memory.h>

public register_CFC_cap, initialize_CFC_cap, CFC_cap_unit_tests
public CFC_cap_column_physics, CFC_cap_set_forcing
public CFC_cap_stock, CFC_cap_end

integer, parameter :: NTR = 2 !< the number of tracers in this module.

!> Contains the concentration array, surface flux, a pointer to Tr in Tr_reg,
!! and some metadata for a single CFC tracer
type, private :: CFC_tracer_data
  type(vardesc) :: desc                     !< A set of metadata for the tracer
  real :: IC_val = 0.0                      !< The initial value assigned to the tracer [mol kg-1].
  real :: land_val = -1.0                   !< The value of the tracer used where land is
                                            !! masked out [mol kg-1].
  character(len=32) :: name                 !< Tracer variable name
  integer :: id_cmor = -1                   !< Diagnostic id
  integer :: id_sfc_flux = -1               !< Surface flux id
  real, pointer, dimension(:,:,:) :: conc   !< The tracer concentration [mol kg-1].
  real, pointer, dimension(:,:) :: sfc_flux !< Surface flux [CU R Z T-1 ~> mol m-2 s-1]
  type(tracer_type), pointer :: tr_ptr      !< pointer to tracer inside Tr_reg
end type CFC_tracer_data

!> The control structure for the CFC_cap tracer package
type, public :: CFC_cap_CS ; private
  logical :: debug              !< If true, write verbose checksums for debugging purposes.
  character(len=200) :: IC_file !< The file in which the CFC initial values can
                                !! be found, or an empty string for internal initilaization.
  logical :: Z_IC_file !< If true, the IC_file is in Z-space.  The default is false.
  type(time_type), pointer :: Time => NULL() !< A pointer to the ocean model's clock.
  type(tracer_registry_type), pointer :: tr_Reg => NULL() !< A pointer to the MOM6 tracer registry
  logical :: tracers_may_reinit !< If true, tracers may be reset via the initialization code
                                !! if they are not found in the restart files.
  type(diag_ctrl), pointer :: diag => NULL() !< A structure that is used to regulate
                                             !! the timing of diagnostic output.
  type(MOM_restart_CS), pointer :: restart_CSp => NULL() !< Model restart control structure

  type(CFC_tracer_data), dimension(NTR) :: CFC_data      !< per-tracer parameters / metadata
  integer :: CFC_BC_year_offset = 0 !< offset to add to model time to get time value used in CFC_BC_file
  type(external_field) :: cfc11_atm_nh_handle !< Handle for time-interpolated CFC11 atm NH
  type(external_field) :: cfc11_atm_sh_handle !< Handle for time-interpolated CFC11 atm SH
  type(external_field) :: cfc12_atm_nh_handle !< Handle for time-interpolated CFC12 atm NH
  type(external_field) :: cfc12_atm_sh_handle !< Handle for time-interpolated CFC12 atm SH
end type CFC_cap_CS


  interface
module function register_CFC_cap(HI, GV, param_file, CS, tr_Reg, restart_CS)
  type(hor_index_type),    intent(in) :: HI         !< A horizontal index type structure.
  type(verticalGrid_type), intent(in) :: GV         !< The ocean's vertical grid structure.
  type(param_file_type),   intent(in) :: param_file !< A structure to parse for run-time parameters.
  type(CFC_cap_CS),        pointer    :: CS         !< A pointer that is set to point to the control
                                                    !! structure for this module.
  type(tracer_registry_type), &
                           pointer    :: tr_Reg     !< A pointer to the tracer registry.
  type(MOM_restart_CS), target, intent(inout) :: restart_CS !< MOM restart control struct

  ! Local variables
  ! This include declares and sets the variable "version".
  logical :: register_CFC_cap

end function register_CFC_cap
module subroutine initialize_CFC_cap(restart, day, G, GV, US, h, diag, OBC, CS)
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
  type(CFC_cap_CS),               pointer    :: CS         !< The control structure returned by a
                                                           !! previous call to register_CFC_cap.

  ! local variables

end subroutine initialize_CFC_cap
module subroutine init_tracer_CFC(h, tr, name, land_val, IC_val, G, GV, US, CS)
  type(ocean_grid_type),                     intent(in)  :: G        !< The ocean's grid structure
  type(verticalGrid_type),                   intent(in)  :: GV       !< The ocean's vertical grid structure.
  type(unit_scale_type),                     intent(in)  :: US       !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in)  :: h        !< Layer thicknesses [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(out) :: tr       !< The tracer concentration array [mol kg-1]
  character(len=*),                          intent(in)  :: name     !< The tracer name
  real,                                      intent(in)  :: land_val !< A value the tracer takes over land [mol kg-1]
  real,                                      intent(in)  :: IC_val   !< The initial condition value for the
                                                                     !! tracer [mol kg-1]
  type(CFC_cap_CS),                          pointer     :: CS       !< The control structure returned by a
                                                                     !! previous call to register_CFC_cap.

  ! local variables
end subroutine init_tracer_CFC
module subroutine CFC_cap_column_physics(h_old, h_new, ea, eb, fluxes, dt, G, GV, US, CS, KPP_CSp, &
                                  nonLocalTrans, evap_CFL_limit, minimum_forcing_depth)
  type(ocean_grid_type),   intent(in) :: G     !< The ocean's grid structure
  type(verticalGrid_type), intent(in) :: GV    !< The ocean's vertical grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in) :: h_old !< Layer thickness before entrainment [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in) :: h_new !< Layer thickness after entrainment [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in) :: ea    !< an array to which the amount of fluid entrained
                                               !! from the layer above during this call will be
                                               !! added [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in) :: eb    !< an array to which the amount of fluid entrained
                                               !! from the layer below during this call will be
                                               !! added [H ~> m or kg m-2].
  type(forcing),           intent(in) :: fluxes!< A structure containing pointers to thermodynamic
                                               !! and tracer forcing fields.  Unused fields have NULL ptrs.
  real,                    intent(in) :: dt    !< The amount of time covered by this call [T ~> s]
  type(unit_scale_type),   intent(in) :: US    !< A dimensional unit scaling type
  type(CFC_cap_CS),        pointer    :: CS    !< The control structure returned by a
                                               !! previous call to register_CFC_cap.
  type(KPP_CS),  optional, pointer    :: KPP_CSp  !< KPP control structure
  real,          optional, intent(in) :: nonLocalTrans(:,:,:) !< Non-local transport [nondim]
  real,          optional, intent(in) :: evap_CFL_limit !< Limit on the fraction of the water that can
                                               !! be fluxed out of the top layer in a timestep [nondim]
  real,          optional, intent(in) :: minimum_forcing_depth !< The smallest depth over which
                                               !! fluxes can be applied [H ~> m or kg m-2]

  ! The arguments to this subroutine are redundant in that
  !     h_new(k) = h_old(k) + ea(k) - eb(k-1) + eb(k) - ea(k+1)

  ! Local variables

end subroutine CFC_cap_column_physics
module function CFC_cap_stock(h, stocks, G, GV, CS, names, units, stock_index)
  type(ocean_grid_type),           intent(in)    :: G      !< The ocean's grid structure.
  type(verticalGrid_type),         intent(in)    :: GV     !< The ocean's vertical grid structure.
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                                   intent(in)    :: h      !< Layer thicknesses [H ~> m or kg m-2].
  type(EFP_type), dimension(:),    intent(out)   :: stocks !< The mass-weighted integrated amount of each
                                                           !! tracer, in kg times concentration units [kg conc]
  type(CFC_cap_CS),                pointer       :: CS     !< The control structure returned by a
                                                           !! previous call to register_CFC_cap.
  character(len=*), dimension(:),  intent(out)   :: names  !< The names of the stocks calculated.
  character(len=*), dimension(:),  intent(out)   :: units  !< The units of the stocks calculated.
  integer, optional,               intent(in)    :: stock_index !< The coded index of a specific
                                                                !! stock being sought.
  integer                                        :: CFC_cap_stock !< The number of stocks calculated here.

  ! Local variables
end function CFC_cap_stock
module subroutine CFC_cap_set_forcing(sfc_state, fluxes, day_start, day_interval, G, US, Rho0, CS)
  type(surface),         intent(in   ) :: sfc_state !< A structure containing fields
                                       !! that describe the surface state of the ocean.
  type(forcing),         intent(inout) :: fluxes !< A structure containing pointers
                                       !! to thermodynamic and tracer forcing fields. Unused fields
                                       !! have NULL ptrs.
  type(time_type),       intent(in)    :: day_start !< Start time of the fluxes.
  type(time_type),       intent(in)    :: day_interval !< Length of time over which these
                                       !! fluxes will be applied.
  type(ocean_grid_type), intent(in)    :: G  !< The ocean's grid structure.
  type(unit_scale_type), intent(in)    :: US !< A dimensional unit scaling type
  real,                  intent(in)    :: Rho0 !< The mean ocean density [R ~> kg m-3]
  type(CFC_cap_CS),      pointer       :: CS !< The control structure returned by a
                                       !! previous call to register_CFC_cap.

  ! Local variables
end subroutine CFC_cap_set_forcing
module subroutine get_solubility(alpha_11, alpha_12, ta, sal , mask)
  real, intent(inout) :: alpha_11 !< The solubility of CFC 11 [mol kg-1 atm-1]
  real, intent(inout) :: alpha_12 !< The solubility of CFC 12 [mol kg-1 atm-1]
  real, intent(in   ) :: ta       !< Absolute sea surface temperature [hectoKelvin]
  real, intent(in   ) :: sal      !< Surface salinity [PSU].
  real, intent(in   ) :: mask     !< ocean mask [nondim]

  ! Local variables

  ! Coefficients for calculating CFC11 solubilities
  ! from Table 5 in Warner and Weiss (1985) DSR, vol 32.



  ! Coefficients for calculating CFC12 solubilities
  ! from Table 5 in Warner and Weiss (1985) DSR, vol 32.




  ! Eq. 9 from Warner and Weiss (1985) DSR, vol 32.
end subroutine get_solubility
module subroutine comp_CFC_schmidt(sst_in, cfc11_sc, cfc12_sc)
  real, intent(in)    :: sst_in   !< The sea surface temperature [degC].
  real, intent(inout) :: cfc11_sc !< Schmidt number of CFC11 [nondim].
  real, intent(inout) :: cfc12_sc !< Schmidt number of CFC12 [nondim].

  !local variables


  ! clip SST to avoid bad values
end subroutine comp_CFC_schmidt
module subroutine CFC_cap_end(CS)
  type(CFC_cap_CS), pointer :: CS !< The control structure returned by a
                                  !! previous call to register_CFC_cap.

  ! local variables

end subroutine CFC_cap_end
logical module function CFC_cap_unit_tests(verbose)
  logical, intent(in) :: verbose !< If true, output additional
                                 !! information for debugging unit tests

  ! Local variables

end function CFC_cap_unit_tests
logical module function compare_values(verbose, test_name, calc, ans, limit)
  logical,             intent(in) :: verbose   !< If true, write results to stdout
  character(len=80),   intent(in) :: test_name !< Brief description of the unit test
  real,                intent(in) :: calc      !< computed value in arbitrary units [A]
  real,                intent(in) :: ans       !< correct value [A]
  real,                intent(in) :: limit     !< value above which test fails [A]

  ! Local variables

end function compare_values
  end interface

end module MOM_CFC_cap
