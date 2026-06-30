! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> Provides the Montgomery potential form of pressure gradient
module MOM_PressureForce_Mont

use MOM_density_integrals, only : int_specific_vol_dp
use MOM_diag_mediator, only : post_data, register_diag_field
use MOM_diag_mediator, only : safe_alloc_ptr, diag_ctrl, time_type
use MOM_error_handler, only : MOM_error, MOM_mesg, FATAL, WARNING, is_root_pe
use MOM_file_parser, only : get_param, log_param, log_version, param_file_type
use MOM_grid, only : ocean_grid_type
use MOM_self_attr_load, only : calc_SAL, SAL_CS
use MOM_tidal_forcing, only : calc_tidal_forcing, tidal_forcing_CS
use MOM_unit_scaling, only : unit_scale_type
use MOM_variables, only : thermo_var_ptrs
use MOM_verticalGrid, only : verticalGrid_type
use MOM_EOS, only : calculate_density, calculate_density_derivs
use MOM_EOS, only : query_compressible

implicit none ; private

#include <MOM_memory.h>

public PressureForce_Mont_Bouss, PressureForce_Mont_nonBouss, Set_pbce_Bouss
public Set_pbce_nonBouss, PressureForce_Mont_init

! A note on unit descriptions in comments: MOM6 uses units that can be rescaled for dimensional
! consistency testing. These are noted in comments with units like Z, H, L, and T, along with
! their mks counterparts with notation like "a velocity [Z T-1 ~> m s-1]".  If the units
! vary with the Boussinesq approximation, the Boussinesq variant is given first.

!> Control structure for the Montgomery potential form of pressure gradient
type, public :: PressureForce_Mont_CS ; private
  logical :: initialized = .false. !< True if this control structure has been initialized.
  logical :: calculate_SAL  !< If true, calculate self-attraction and loading.
  logical :: tides          !< If true, apply tidal momentum forcing.
  real    :: Rho0           !< The density used in the Boussinesq
                            !! approximation [R ~> kg m-3].
  real    :: GFS_scale      !< Ratio between gravity applied to top interface and the
                            !! gravitational acceleration of the planet [nondim].
                            !! Usually this ratio is 1.
  type(time_type), pointer :: Time => NULL() !< A pointer to the ocean model's clock.
  type(diag_ctrl), pointer :: diag => NULL() !< A structure that is used to regulate
                            !! the timing of diagnostic output.
  real, allocatable :: PFu_bc(:,:,:) !< Zonal accelerations due to pressure gradients
                            !! deriving from density gradients within layers [L T-2 ~> m s-2].
  real, allocatable :: PFv_bc(:,:,:) !< Meridional accelerations due to pressure gradients
                            !! deriving from density gradients within layers [L T-2 ~> m s-2].
  !>@{ Diagnostic IDs
  integer :: id_PFu_bc = -1, id_PFv_bc = -1, id_e_sal = -1
  integer :: id_e_tide = -1,  id_e_tide_eq = -1, id_e_tide_sal = -1
  !>@}
  type(SAL_CS), pointer :: SAL_CSp => NULL() !< SAL control structure
  type(tidal_forcing_CS), pointer :: tides_CSp => NULL() !< The tidal forcing control structure
end type PressureForce_Mont_CS


  interface
module subroutine PressureForce_Mont_nonBouss(h, tv, PFu, PFv, G, GV, US, CS, p_atm, pbce, eta)
  type(ocean_grid_type),                      intent(in)  :: G   !< Ocean grid structure.
  type(verticalGrid_type),                    intent(in)  :: GV  !< Vertical grid structure.
  type(unit_scale_type),                      intent(in)  :: US  !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),  intent(in)  :: h   !< Layer thickness, [H ~> kg m-2].
  type(thermo_var_ptrs),                      intent(in)  :: tv  !< Thermodynamic variables.
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), intent(out) :: PFu !< Zonal acceleration due to pressure gradients
                                                                 !! (equal to -dM/dx) [L T-2 ~> m s-2].
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), intent(out) :: PFv !< Meridional acceleration due to pressure gradients
                                                                 !! (equal to -dM/dy) [L T-2 ~> m s-2].
  type(PressureForce_Mont_CS),                intent(inout) :: CS  !< Control structure for Montgomery potential PGF
  real, dimension(:,:),                       pointer     :: p_atm !< The pressure at the ice-ocean or
                                                                 !! atmosphere-ocean [R L2 T-2 ~> Pa].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                                    optional, intent(out) :: pbce !< The baroclinic pressure anomaly in
                                                                 !! each layer due to free surface height anomalies,
                                                                 !! [L2 T-2 H-1 ~> m s-2 or m4 kg-1 s-2].
  real, dimension(SZI_(G),SZJ_(G)), optional, intent(out) :: eta !< The total column mass used to calculate
                                                                 !! PFu and PFv [H ~> kg m-2].
  ! Local variables
                ! p may be adjusted (with a nonlinear equation of state) so that
                ! its derivative compensates for the adiabatic compressibility
                ! in seawater, but p will still be close to the pressure.
end subroutine PressureForce_Mont_nonBouss
module subroutine PressureForce_Mont_Bouss(h, tv, PFu, PFv, G, GV, US, CS, p_atm, pbce, eta)
  type(ocean_grid_type),                      intent(in)  :: G   !< Ocean grid structure.
  type(verticalGrid_type),                    intent(in)  :: GV  !< Vertical grid structure.
  type(unit_scale_type),                      intent(in)  :: US  !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),  intent(in)  :: h   !< Layer thickness [H ~> m].
  type(thermo_var_ptrs),                      intent(in)  :: tv  !< Thermodynamic variables.
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), intent(out) :: PFu !< Zonal acceleration due to pressure gradients
                                                                 !! (equal to -dM/dx) [L T-2 ~> m s-2].
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), intent(out) :: PFv !< Meridional acceleration due to pressure gradients
                                                                 !! (equal to -dM/dy) [L T-2 ~> m s-2].
  type(PressureForce_Mont_CS),                intent(inout) :: CS  !< Control structure for Montgomery potential PGF
  real, dimension(:,:),                       pointer     :: p_atm !< The pressure at the ice-ocean or
                                                                !! atmosphere-ocean [R L2 T-2 ~> Pa].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), optional, intent(out) :: pbce !< The baroclinic pressure anomaly in
                                                                !! each layer due to free surface height anomalies
                                                                !! [L2 T-2 H-1 ~> m s-2].
  real, dimension(SZI_(G),SZJ_(G)),          optional, intent(out) :: eta !< Free surface height [H ~> m].
  ! Local variables
                ! corrected e times (G_Earth/Rho0) [L2 Z-1 T-2 ~> m s-2].
                ! e may be adjusted (with a nonlinear equation of state) so that
                ! its derivative compensates for the adiabatic compressibility
                ! in seawater, but e will still be close to the interface depth.
end subroutine PressureForce_Mont_Bouss
module subroutine Set_pbce_Bouss(e, tv, G, GV, US, Rho0, GFS_scale, pbce, rho_star)
  type(ocean_grid_type),                intent(in)  :: G    !< Ocean grid structure
  type(verticalGrid_type),              intent(in)  :: GV   !< Vertical grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)+1), intent(in) :: e !< Interface height [Z ~> m].
  type(thermo_var_ptrs),                intent(in)  :: tv   !< Thermodynamic variables
  type(unit_scale_type),                intent(in)  :: US   !< A dimensional unit scaling type
  real,                                 intent(in)  :: Rho0 !< The "Boussinesq" ocean density [R ~> kg m-3].
  real,                                 intent(in)  :: GFS_scale !< Ratio between gravity applied to top
                                                            !! interface and the gravitational acceleration of
                                                            !! the planet [nondim]. Usually this ratio is 1.
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                                        intent(out) :: pbce !< The baroclinic pressure anomaly in each layer due
                                                            !! to free surface height anomalies
                                                            !! [L2 T-2 H-1 ~> m s-2].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                              optional, intent(in)  :: rho_star !< The layer densities (maybe compressibility
                                                            !! compensated), times g/rho_0 [L2 Z-1 T-2 ~> m s-2].

  ! Local variables
                             ! an equation of state.
                             ! in roundoff and can be neglected [Z ~> m].

end subroutine Set_pbce_Bouss
module subroutine Set_pbce_nonBouss(p, tv, G, GV, US, GFS_scale, pbce, alpha_star)
  type(ocean_grid_type),                intent(in)  :: G  !< Ocean grid structure
  type(verticalGrid_type),              intent(in)  :: GV !< Vertical grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)+1), intent(in) :: p !< Interface pressures [R L2 T-2 ~> Pa].
  type(thermo_var_ptrs),                intent(in)  :: tv !< Thermodynamic variables
  type(unit_scale_type),                intent(in)  :: US !< A dimensional unit scaling type
  real,                                 intent(in)  :: GFS_scale !< Ratio between gravity applied to top
                                                          !! interface and the gravitational acceleration of
                                                          !! the planet [nondim]. Usually this ratio is 1.
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(out) :: pbce !< The baroclinic pressure anomaly in each
                                                          !! layer due to free surface height anomalies
                                                          !! [L2 H-1 T-2 ~> m4 kg-1 s-2].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), optional, intent(in) :: alpha_star !< The layer specific volumes
                                                          !! (maybe compressibility compensated) [R-1 ~> m3 kg-1].
  ! Local variables
end subroutine Set_pbce_nonBouss
module subroutine PressureForce_Mont_init(Time, G, GV, US, param_file, diag, CS, SAL_CSp, tides_CSp)
  type(time_type), target, intent(in)    :: Time !< Current model time
  type(ocean_grid_type),   intent(in)    :: G  !< ocean grid structure
  type(verticalGrid_type), intent(in)    :: GV !< Vertical grid structure
  type(unit_scale_type),   intent(in)    :: US !< A dimensional unit scaling type
  type(param_file_type),   intent(in)    :: param_file !< Parameter file handles
  type(diag_ctrl), target, intent(inout) :: diag !< Diagnostics control structure
  type(PressureForce_Mont_CS), intent(inout) :: CS !< Montgomery PGF control structure
  type(SAL_CS), intent(in), target, optional :: SAL_CSp !< SAL control structure
  type(tidal_forcing_CS), intent(in), target, optional :: tides_CSp !< Tides control structure

  ! Local variables
  ! This include declares and sets the variable "version".

end subroutine PressureForce_Mont_init
  end interface

end module MOM_PressureForce_Mont
