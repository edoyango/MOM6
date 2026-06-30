! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> Finite volume pressure gradient (integrated by quadrature or analytically)
module MOM_PressureForce_FV

use MOM_debugging, only : hchksum, uvchksum
use MOM_diag_mediator, only : post_data, register_diag_field
use MOM_diag_mediator, only : safe_alloc_ptr, diag_ctrl, time_type
use MOM_error_handler, only : MOM_error, FATAL, WARNING, is_root_pe
use MOM_file_parser, only : get_param, log_param, log_version, param_file_type
use MOM_grid, only : ocean_grid_type
use MOM_PressureForce_Mont, only : set_pbce_Bouss, set_pbce_nonBouss
use MOM_self_attr_load, only : calc_SAL, SAL_CS
use MOM_tidal_forcing, only : calc_tidal_forcing, tidal_forcing_CS
use MOM_tidal_forcing, only : calc_tidal_forcing_legacy
use MOM_unit_scaling, only : unit_scale_type
use MOM_variables, only : thermo_var_ptrs, accel_diag_ptrs
use MOM_verticalGrid, only : verticalGrid_type
use MOM_EOS, only : calculate_density, calculate_spec_vol, EOS_domain
use MOM_density_integrals, only : int_density_dz, int_specific_vol_dp
use MOM_density_integrals, only : int_density_dz_generic_plm, int_density_dz_generic_ppm
use MOM_density_integrals, only : int_spec_vol_dp_generic_plm
use MOM_density_integrals, only : int_density_dz_generic_pcm, int_spec_vol_dp_generic_pcm
use MOM_density_integrals, only : diagnose_mass_weight_Z, diagnose_mass_weight_p
use MOM_ALE, only : TS_PLM_edge_values, TS_PPM_edge_values, TS_PLM_WLS_edge_values, ALE_CS

implicit none ; private

#include <MOM_memory.h>

public PressureForce_FV_init
public PressureForce_FV_Bouss, PressureForce_FV_nonBouss

! A note on unit descriptions in comments: MOM6 uses units that can be rescaled for dimensional
! consistency testing. These are noted in comments with units like Z, H, L, and T, along with
! their mks counterparts with notation like "a velocity [Z T-1 ~> m s-1]".  If the units
! vary with the Boussinesq approximation, the Boussinesq variant is given first.

!> Finite volume pressure gradient control structure
type, public :: PressureForce_FV_CS ; private
  logical :: initialized = .false. !< True if this control structure has been initialized.
  logical :: calculate_SAL = .false. !< If true, calculate self-attraction and loading.
  logical :: sal_use_bpa = .false. !< If true, use bottom pressure anomaly instead of SSH
                                   !! to calculate SAL.
  logical :: tides = .false.       !< If true, apply tidal momentum forcing.
  real    :: rho_ref        !< The reference density that is subtracted off when calculating pressure
                            !! gradient forces [R ~> kg m-3].
  logical :: rho_ref_bug    !< If true, recover a bug that mixes GV%Rho0 and CS%rho_ref in Boussinesq mode.
  real    :: GFS_scale      !< A scaling of the surface pressure gradients to
                            !! allow the use of a reduced gravity model [nondim].
  type(time_type), pointer :: Time !< A pointer to the ocean model's clock.
  type(diag_ctrl), pointer :: diag !< A structure that is used to regulate the
                            !! timing of diagnostic output.
  integer :: MassWghtInterp !< A flag indicating whether and how to use mass weighting in T/S interpolation
  logical :: correction_intxpa !< If true, apply a correction to the value of intxpa at a selected
                            !! interface under ice, using matching at the end values along with a
                            !! 5-point quadrature integral of the hydrostatic pressure or height
                            !! changes along that interface.  The selected interface is either at the
                            !! ocean's surface or in the interior, depending on reset_intxpa_integral.
  logical :: reset_intxpa_integral !< If true and the surface displacement between adjacent cells
                            !! exceeds the vertical grid spacing, reset intxpa at the interface below
                            !! a trusted interior cell.  (This often applies in ice shelf cavities.)
  logical :: MassWghtInterpVanOnly !< If true, don't do mass weighting of T/S interpolation unless vanished
  logical :: reset_intxpa_flattest !< If true, use flattest interface rather than top for reset integral
                                   !! in cases where no best nonvanished interface
  real    :: h_nonvanished  !< A minimal layer thickness that indicates that a layer is thick enough
                            !! to usefully reestimate the pressure integral across the interface
                            !! below it [H ~> m or kg m-2]
  logical :: use_inaccurate_pgf_rho_anom !< If true, uses the older and less accurate
                            !! method to calculate density anomalies, as used prior to
                            !! March 2018.
  logical :: boundary_extrap !< Indicate whether high-order boundary
                            !! extrapolation should be used within boundary cells

  logical :: reconstruct    !< If true, polynomial profiles of T & S will be
                            !! reconstructed and used in the integrals for the
                            !! finite volume pressure gradient calculation.
                            !! The default depends on whether regridding is being used.

  integer :: Recon_Scheme   !< Order of the polynomial of the reconstruction of T & S
                            !! for the finite volume pressure gradient calculation.
                            !! By the default (1) is for a piecewise linear method

  logical :: debug          !< If true, write verbose checksums for debugging purposes.
  logical :: use_SSH_in_Z0p !< If true, adjust the height at which the pressure used in the
                            !! equation of state is 0 to account for the displacement of the sea
                            !! surface including adjustments for atmospheric or sea-ice pressure.
  logical :: use_stanley_pgf  !< If true, turn on Stanley parameterization in the PGF
  logical :: bq_sal_tides = .false. !< If true, use an alternative method for SAL and tides
                                    !! in Boussinesq mode
  integer :: tides_answer_date = 99991231 !< Recover old answers with tides
  integer :: id_e_tide = -1 !< Diagnostic identifier
  integer :: id_e_tidal_eq = -1 !< Diagnostic identifier
  integer :: id_e_tidal_sal = -1 !< Diagnostic identifier
  integer :: id_e_sal = -1 !< Diagnostic identifier
  integer :: id_rho_pgf = -1 !< Diagnostic identifier
  integer :: id_rho_stanley_pgf = -1 !< Diagnostic identifier
  integer :: id_p_stanley = -1 !< Diagnostic identifier
  integer :: id_MassWt_u = -1 !< Diagnostic identifier
  integer :: id_MassWt_v = -1 !< Diagnostic identifier
  integer :: id_sal_u = -1 !< Diagnostic identifier
  integer :: id_sal_v = -1 !< Diagnostic identifier
  integer :: id_tides_u = -1 !< Diagnostic identifier
  integer :: id_tides_v = -1 !< Diagnostic identifier
  type(SAL_CS), pointer :: SAL_CSp => NULL() !< SAL control structure
  type(tidal_forcing_CS), pointer :: tides_CSp => NULL() !< Tides control structure
end type PressureForce_FV_CS


  interface
module subroutine PressureForce_FV_nonBouss(h, tv, PFu, PFv, G, GV, US, CS, ALE_CSp, ADp, p_atm, pbce, eta)
  type(ocean_grid_type),                      intent(in)  :: G   !< Ocean grid structure
  type(verticalGrid_type),                    intent(in)  :: GV  !< Vertical grid structure
  type(unit_scale_type),                      intent(in)  :: US  !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),  intent(in)  :: h   !< Layer thickness [H ~> kg m-2]
  type(thermo_var_ptrs),                      intent(in)  :: tv  !< Thermodynamic variables
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), intent(out) :: PFu !< Zonal acceleration [L T-2 ~> m s-2]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), intent(out) :: PFv !< Meridional acceleration [L T-2 ~> m s-2]
  type(PressureForce_FV_CS),                  intent(in)  :: CS  !< Finite volume PGF control structure
  type(ALE_CS),                               pointer     :: ALE_CSp !< ALE control structure
  type(accel_diag_ptrs),                      pointer     :: ADp !< Acceleration diagnostic pointers
  real, dimension(:,:),                       pointer     :: p_atm !< The pressure at the ice-ocean
                                                           !! or atmosphere-ocean interface [R L2 T-2 ~> Pa].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), optional, intent(out) :: pbce !< The baroclinic pressure
                                                           !! anomaly in each layer due to eta anomalies
                                                           !! [L2 T-2 H-1 ~> m4 s-2 kg-1].
  real, dimension(SZI_(G),SZJ_(G)),          optional, intent(out) :: eta !< The total column mass used to
                                                           !! calculate PFu and PFv [H ~> kg m-2].
  ! Local variables
end subroutine PressureForce_FV_nonBouss
module subroutine PressureForce_FV_Bouss(h, tv, PFu, PFv, G, GV, US, CS, ALE_CSp, ADp, p_atm, pbce, eta)
  type(ocean_grid_type),                      intent(in)  :: G   !< Ocean grid structure
  type(verticalGrid_type),                    intent(in)  :: GV  !< Vertical grid structure
  type(unit_scale_type),                      intent(in)  :: US  !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),  intent(in)  :: h   !< Layer thickness [H ~> m]
  type(thermo_var_ptrs),                      intent(in)  :: tv  !< Thermodynamic variables
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), intent(out) :: PFu !< Zonal acceleration [L T-2 ~> m s-2]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), intent(out) :: PFv !< Meridional acceleration [L T-2 ~> m s-2]
  type(PressureForce_FV_CS),                  intent(in)  :: CS  !< Finite volume PGF control structure
  type(ALE_CS),                               pointer     :: ALE_CSp !< ALE control structure
  type(accel_diag_ptrs),                      pointer     :: ADp !< Acceleration diagnostic pointers
  real, dimension(:,:),                       pointer     :: p_atm !< The pressure at the ice-ocean
                                                         !! or atmosphere-ocean interface [R L2 T-2 ~> Pa].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), optional, intent(out) :: pbce !< The baroclinic pressure
                                                         !! anomaly in each layer due to eta anomalies
                                                         !! [L2 T-2 H-1 ~> m s-2].
  real, dimension(SZI_(G),SZJ_(G)),          optional, intent(out) :: eta !< The sea-surface height used to
                                                         !! calculate PFu and PFv [H ~> m], with any
                                                         !! tidal contributions.
  ! Local variables
end subroutine PressureForce_FV_Bouss
module subroutine PressureForce_FV_init(Time, G, GV, US, param_file, diag, CS, ADp, SAL_CSp, tides_CSp)
  type(time_type), target,    intent(in)    :: Time !< Current model time
  type(ocean_grid_type),      intent(in)    :: G  !< Ocean grid structure
  type(verticalGrid_type),    intent(in)    :: GV !< Vertical grid structure
  type(unit_scale_type),      intent(in)    :: US !< A dimensional unit scaling type
  type(param_file_type),      intent(in)    :: param_file !< Parameter file handles
  type(diag_ctrl), target,    intent(inout) :: diag !< Diagnostics control structure
  type(PressureForce_FV_CS),  intent(inout) :: CS !< Finite volume PGF control structure
  type(accel_diag_ptrs),      pointer       :: ADp !< Acceleration diagnostic pointers
  type(SAL_CS),           intent(in), target, optional :: SAL_CSp !< SAL control structure
  type(tidal_forcing_CS), intent(in), target, optional :: tides_CSp !< Tides control structure

  ! Local variables
                           ! temperature variance [nondim]
                          ! recreate the bugs, or if false bugs are only used if actively selected.
  ! This include declares and sets the variable "version".

end subroutine PressureForce_FV_init
  end interface

end module MOM_PressureForce_FV
