! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> Provides functions for some diabatic processes such as fraxil, brine rejection,
!! tendency due to surface flux divergence.
module MOM_diagnose_mld

use MOM_diag_mediator, only : post_data
use MOM_diag_mediator, only : diag_ctrl
use MOM_EOS,           only : calculate_density, calculate_TFreeze, EOS_domain
use MOM_EOS,           only : calculate_specific_vol_derivs, calculate_density_derivs
use MOM_error_handler, only : MOM_error, FATAL, WARNING
use MOM_grid,          only : ocean_grid_type
use MOM_interface_heights, only : thickness_to_dz
use MOM_unit_scaling,  only : unit_scale_type
use MOM_variables,     only : thermo_var_ptrs
use MOM_verticalGrid,  only : verticalGrid_type

implicit none ; private

#include <MOM_memory.h>

public diagnoseMLDbyEnergy, diagnoseMLDbyDensityDifference

! A note on unit descriptions in comments: MOM6 uses units that can be rescaled for dimensional
! consistency testing. These are noted in comments with units like Z, H, L, and T, along with
! their mks counterparts with notation like "a velocity [Z T-1 ~> m s-1]".  If the units
! vary with the Boussinesq approximation, the Boussinesq variant is given first.


  interface
module subroutine diagnoseMLDbyDensityDifference(id_MLD, h, tv, densityDiff, G, GV, US, diagPtr, &
                                          ref_h_mld, id_ref_z, id_ref_rho, id_N2subML, id_MLDsq, &
                                          dz_subML, MLD_out)
  type(ocean_grid_type),   intent(in) :: G           !< Grid type
  type(verticalGrid_type), intent(in) :: GV          !< ocean vertical grid structure
  type(unit_scale_type),   intent(in) :: US          !< A dimensional unit scaling type
  integer,                 intent(in) :: id_MLD      !< Handle (ID) of MLD diagnostic
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in) :: h           !< Layer thickness [H ~> m or kg m-2]
  type(thermo_var_ptrs),   intent(in) :: tv          !< Structure containing pointers to any
                                                     !! available thermodynamic fields.
  real,                    intent(in) :: densityDiff !< Density difference to determine MLD [R ~> kg m-3]
  type(diag_ctrl),         pointer    :: diagPtr     !< Diagnostics structure
  real,                    intent(in) :: ref_h_mld   !< Depth of the calculated "surface" densisty [Z ~> m]
  integer,                 intent(in) :: id_ref_z    !< Handle (ID) of reference depth diagnostic
  integer,                 intent(in) :: id_ref_rho  !< Handle (ID) of reference density diagnostic
  integer,       optional, intent(in) :: id_N2subML  !< Optional handle (ID) of subML stratification
  integer,       optional, intent(in) :: id_MLDsq    !< Optional handle (ID) of squared MLD
  real,          optional, intent(in) :: dz_subML    !< The distance over which to calculate N2subML
                                                     !! or 50 m if missing [Z ~> m]
  real, dimension(SZI_(G),SZJ_(G)), &
              optional, intent(out)   :: MLD_out     !< Send MLD to other routines [Z ~> m]

  ! Local variables
                                               ! have been stored already.
                           ! reference density [H T-2 R-1 ~> m4 s-2 kg-1 or m s-2].
                                                  ! the MLD. It can be saved as a diagnostic [R ~> kg m-3].

end subroutine diagnoseMLDbyDensityDifference
module subroutine diagnoseMLDbyEnergy(id_MLD, h, tv, G, GV, US, Mixing_Energy, k_bounds, diagPtr, OM4_iteration, MLD_out)
  ! Author: Brandon Reichl
  ! Date: October 2, 2020
  ! //
  ! *Note that gravity is assumed constant everywhere and divided out of all calculations.
  !
  ! This code has been written to step through the columns layer by layer, summing the PE
  ! change inferred by mixing the layer with all layers above.  When the change exceeds a
  ! threshold (determined by input array Mixing_Energy), the code needs to solve for how far
  ! into this layer the threshold PE change occurs (assuming constant density layers).
  ! This is expressed here via solving the function F(X) = 0 where:
  ! F(X) = 0.5 * ( Ca*X^3/(D1+X) + Cb*X^2/(D1+X) + Cc*X/(D1+X) + Dc/(D1+X)
  !                + Ca2*X^2 + Cb2*X + Cc2)
  ! where all coefficients are determined by the previous mixed layer depth, the
  ! density of the previous mixed layer, the present layer thickness, and the present
  ! layer density.  This equation is worked out by computing the total PE assuming constant
  ! density in the mixed layer as well as in the remaining part of the present layer that is
  ! not mixed.
  ! To solve for X in this equation a Newton's method iteration is employed, which
  ! converges extremely quickly (usually 1 guess) since this equation turns out to be rather
  ! linear for PE change with increasing X.
  ! Input parameters:
  integer, dimension(3),   intent(in) :: id_MLD      !< Energy output diagnostic IDs
  type(ocean_grid_type),   intent(in) :: G           !< Grid type
  type(verticalGrid_type), intent(in) :: GV          !< ocean vertical grid structure
  type(unit_scale_type),   intent(in) :: US          !< A dimensional unit scaling type
  real, dimension(3),      intent(in) :: Mixing_Energy !< Energy values for up to 3 MLDs [R Z3 T-2 ~> J m-2]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in) :: h           !< Layer thickness [H ~> m or kg m-2]
  type(thermo_var_ptrs),   intent(in) :: tv          !< Structure containing pointers to any
                                                     !! available thermodynamic fields.
  type(diag_ctrl),         pointer    :: diagPtr     !< Diagnostics structure
  integer, dimension(2), intent(in)   :: k_bounds    !< vertical interface bounds to apply calculations
  logical, optional, intent(in)       :: OM4_iteration !< Uses a legacy version of the MLD iteration
                                                     !! it is kept to reproduce OM4 output
  real, dimension(SZI_(G),SZJ_(G)), &
              optional, intent(out)   :: MLD_out     !< Send MLD to other routines [Z ~> m]

  ! Local variables
                                          ! depth calculation [R L2 T-2 ~> Pa]

                                  ! for the energy used to mix to the diagnosed depth [nondim]
                                  ! of H_ML, divided by the gravitational acceleration [R Z2 ~> kg m-1]
                                  ! of H_ML, divided by the gravitational acceleration [R Z2 ~> kg m-1]
                                  ! of H_ML_TST, divided by the gravitational acceleration [R Z2 ~> kg m-1]

  ! These are all temporary variables used to shorten the expressions in the iterations.


end subroutine diagnoseMLDbyEnergy
  end interface

end module MOM_diagnose_mld
