! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> Provides integrals of density
module MOM_density_integrals

use MOM_EOS,              only : EOS_type
use MOM_EOS,              only : EOS_quadrature, EOS_domain
use MOM_EOS,              only : analytic_int_density_dz
use MOM_EOS,              only : analytic_int_specific_vol_dp
use MOM_EOS,              only : calculate_density
use MOM_EOS,              only : calculate_spec_vol
use MOM_EOS,              only : calculate_specific_vol_derivs
use MOM_EOS,              only : average_specific_vol
use MOM_error_handler,    only : MOM_error, FATAL, WARNING, MOM_mesg
use MOM_hor_index,        only : hor_index_type
use MOM_string_functions, only : uppercase
use MOM_variables,        only : thermo_var_ptrs
use MOM_unit_scaling,     only : unit_scale_type
use MOM_verticalGrid,     only : verticalGrid_type

implicit none ; private

#include <MOM_memory.h>

public int_density_dz
public int_density_dz_generic_pcm
public int_density_dz_generic_plm
public int_density_dz_generic_ppm
public int_specific_vol_dp
public int_spec_vol_dp_generic_pcm
public int_spec_vol_dp_generic_plm
public avg_specific_vol
public find_depth_of_pressure_in_cell
public diagnose_mass_weight_Z, diagnose_mass_weight_p


  interface
module subroutine int_density_dz(T, S, z_t, z_b, rho_ref, rho_0, G_e, HI, EOS, US, dpa, &
                          intz_dpa, intx_dpa, inty_dpa, bathyT, SSH, dz_neglect, MassWghtInterp, Z_0p, &
                          MassWghtInterpVanOnly, h_nv)
  type(hor_index_type), intent(in)  :: HI  !< Ocean horizontal index structures for the arrays
  real, dimension(SZI_(HI),SZJ_(HI)), &
                        intent(in)  :: T   !< Potential temperature referenced to the surface [C ~> degC]
  real, dimension(SZI_(HI),SZJ_(HI)), &
                        intent(in)  :: S   !< Salinity [S ~> ppt]
  real, dimension(SZI_(HI),SZJ_(HI)), &
                        intent(in)  :: z_t !< Height at the top of the layer in depth units [Z ~> m]
  real, dimension(SZI_(HI),SZJ_(HI)), &
                        intent(in)  :: z_b !< Height at the bottom of the layer [Z ~> m]
  real,                 intent(in)  :: rho_ref !< A mean density [R ~> kg m-3], that is
                                           !! subtracted out to reduce the magnitude of each of the
                                           !! integrals.
  real,                 intent(in)  :: rho_0 !< A density [R ~> kg m-3], that is used
                                           !! to calculate the pressure (as p~=-z*rho_0*G_e)
                                           !! used in the equation of state.
  real,                 intent(in)  :: G_e !< The Earth's gravitational acceleration
                                           !! [L2 Z-1 T-2 ~> m s-2]
  type(EOS_type),       intent(in)  :: EOS !< Equation of state structure
  type(unit_scale_type), intent(in) :: US  !< A dimensional unit scaling type
  real, dimension(SZI_(HI),SZJ_(HI)), &
                      intent(inout) :: dpa !< The change in the pressure anomaly
                                           !! across the layer [R L2 T-2 ~> Pa]
  real, dimension(SZI_(HI),SZJ_(HI)), &
            optional, intent(inout) :: intz_dpa !< The integral through the thickness of the
                                           !! layer of the pressure anomaly relative to the
                                           !! anomaly at the top of the layer [R L2 Z T-2 ~> Pa m]
  real, dimension(SZIB_(HI),SZJ_(HI)), &
            optional, intent(inout) :: intx_dpa !< The integral in x of the difference between
                                          !! the pressure anomaly at the top and bottom of the
                                          !! layer divided by the x grid spacing [R L2 T-2 ~> Pa]
  real, dimension(SZI_(HI),SZJB_(HI)), &
            optional, intent(inout) :: inty_dpa !< The integral in y of the difference between
                                          !! the pressure anomaly at the top and bottom of the
                                          !! layer divided by the y grid spacing [R L2 T-2 ~> Pa]
  real, dimension(SZI_(HI),SZJ_(HI)), &
              optional, intent(in)  :: bathyT !< The depth of the bathymetry [Z ~> m]
  real, dimension(SZI_(HI),SZJ_(HI)), &
              optional, intent(in)  :: SSH !< The sea surface height [Z ~> m]
  real,       optional, intent(in)  :: dz_neglect !< A minuscule thickness change [Z ~> m]
  integer,    optional, intent(in)  :: MassWghtInterp !< A flag indicating whether and how to use
                                           !! mass weighting to interpolate T/S in integrals
  logical,    optional, intent(in)  :: MassWghtInterpVanOnly !< If true, does not do mass weighting
                                           !! of T/S unless one side smaller than h_nv (i.e. vanished)
  real,       optional, intent(in)  :: h_nv !< Nonvanished height [Z ~> m]

  real, dimension(HI%isd:HI%ied,HI%jsd:HI%jed), &
              optional, intent(in)  :: Z_0p !< The height at which the pressure is 0 [Z ~> m]

end subroutine int_density_dz
module subroutine int_density_dz_generic_pcm(T, S, z_t, z_b, rho_ref, rho_0, G_e, HI, &
                                      EOS, US, dpa, intz_dpa, intx_dpa, inty_dpa, bathyT, SSH, &
                                      dz_neglect, MassWghtInterp, use_inaccurate_form, Z_0p, &
                                      MassWghtInterpVanOnly, h_nv)
  type(hor_index_type), intent(in)  :: HI  !< Horizontal index type for input variables.
  real, dimension(SZI_(HI),SZJ_(HI)), &
                        intent(in)  :: T  !< Potential temperature of the layer [C ~> degC]
  real, dimension(SZI_(HI),SZJ_(HI)), &
                        intent(in)  :: S  !< Salinity of the layer [S ~> ppt]
  real, dimension(SZI_(HI),SZJ_(HI)), &
                        intent(in)  :: z_t !< Height at the top of the layer in depth units [Z ~> m]
  real, dimension(SZI_(HI),SZJ_(HI)), &
                        intent(in)  :: z_b !< Height at the bottom of the layer [Z ~> m]
  real,                 intent(in)  :: rho_ref !< A mean density [R ~> kg m-3], that is
                                          !! subtracted out to reduce the magnitude
                                          !! of each of the integrals.
  real,                 intent(in)  :: rho_0 !< A density [R ~> kg m-3], that is used
                                          !! to calculate the pressure (as p~=-z*rho_0*G_e)
                                          !! used in the equation of state.
  real,                 intent(in)  :: G_e !< The Earth's gravitational acceleration
                                          !! [L2 Z-1 T-2 ~> m s-2]
  type(EOS_type),       intent(in)  :: EOS !< Equation of state structure
  type(unit_scale_type), intent(in) :: US !< A dimensional unit scaling type
  real, dimension(SZI_(HI),SZJ_(HI)), &
                      intent(inout) :: dpa !< The change in the pressure anomaly
                                          !! across the layer [R L2 T-2 ~> Pa]
  real, dimension(SZI_(HI),SZJ_(HI)), &
            optional, intent(inout) :: intz_dpa !< The integral through the thickness of the
                                          !! layer of the pressure anomaly relative to the
                                          !! anomaly at the top of the layer [R L2 Z T-2 ~> Pa m]
  real, dimension(SZIB_(HI),SZJ_(HI)), &
            optional, intent(inout) :: intx_dpa !< The integral in x of the difference between
                                          !! the pressure anomaly at the top and bottom of the
                                          !! layer divided by the x grid spacing [R L2 T-2 ~> Pa]
  real, dimension(SZI_(HI),SZJB_(HI)), &
            optional, intent(inout) :: inty_dpa !< The integral in y of the difference between
                                          !! the pressure anomaly at the top and bottom of the
                                          !! layer divided by the y grid spacing [R L2 T-2 ~> Pa]
  real, dimension(SZI_(HI),SZJ_(HI)), &
              optional, intent(in)  :: bathyT !< The depth of the bathymetry [Z ~> m]
  real, dimension(SZI_(HI),SZJ_(HI)), &
              optional, intent(in)  :: SSH !< The sea surface height [Z ~> m]
  real,       optional, intent(in)  :: dz_neglect !< A minuscule thickness change [Z ~> m]
  integer,    optional, intent(in)  :: MassWghtInterp !< A flag indicating whether and how to use
                                           !! mass weighting to interpolate T/S in integrals
  logical,    optional, intent(in)  :: MassWghtInterpVanOnly !< If true, does not do mass weighting
                                           !! of T/S unless one side smaller than h_nv (i.e. vanished)
  real,       optional, intent(in)  :: h_nv !< Nonvanished height [Z ~> m]
  logical,    optional, intent(in)  :: use_inaccurate_form !< If true, uses an inaccurate form of
                                          !! density anomalies, as was used prior to March 2018.
  real, dimension(HI%isd:HI%ied,HI%jsd:HI%jed), &
              optional, intent(in)  :: Z_0p  !< The height at which the pressure is 0 [Z ~> m]

  ! Local variables
                     ! with height at the 5 sub-column locations [R L2 T-2 ~> Pa]
                                    ! if at least one side vanished (0 or 1) [nondim]
                         ! of density anomalies.

  ! These array bounds work for the indexing convention of the input arrays, but
  ! on the computational domain defined for the output arrays.
end subroutine int_density_dz_generic_pcm
module subroutine int_density_dz_generic_plm(k, tv, T_t, T_b, S_t, S_b, e, rho_ref, &
                                      rho_0, G_e, dz_subroundoff, bathyT, HI, GV, EOS, US, use_stanley_eos, dpa, &
                                      intz_dpa, intx_dpa, inty_dpa, MassWghtInterp, &
                                      use_inaccurate_form, Z_0p, MassWghtInterpVanOnly, h_nv)
  integer,              intent(in)  :: k   !< Layer index to calculate integrals for
  type(hor_index_type), intent(in)  :: HI  !< Ocean horizontal index structures for the input arrays
  type(verticalGrid_type), intent(in) :: GV !< Vertical grid structure
  type(thermo_var_ptrs), intent(in) :: tv  !< Thermodynamic variables
  real, dimension(SZI_(HI),SZJ_(HI),SZK_(GV)), &
                        intent(in)  :: T_t !< Potential temperature at the cell top [C ~> degC]
  real, dimension(SZI_(HI),SZJ_(HI),SZK_(GV)), &
                        intent(in)  :: T_b !< Potential temperature at the cell bottom [C ~> degC]
  real, dimension(SZI_(HI),SZJ_(HI),SZK_(GV)), &
                        intent(in)  :: S_t !< Salinity at the cell top [S ~> ppt]
  real, dimension(SZI_(HI),SZJ_(HI),SZK_(GV)), &
                        intent(in)  :: S_b !< Salinity at the cell bottom [S ~> ppt]
  real, dimension(SZI_(HI),SZJ_(HI),SZK_(GV)+1), &
                        intent(in)  :: e   !< Height of interfaces [Z ~> m]
  real,                 intent(in)  :: rho_ref !< A mean density [R ~> kg m-3], that is subtracted
                                           !! out to reduce the magnitude of each of the integrals.
  real,                 intent(in)  :: rho_0 !< A density [R ~> kg m-3], that is used to calculate
                                           !! the pressure (as p~=-z*rho_0*G_e) used in the equation of state.
  real,                 intent(in)  :: G_e !< The Earth's gravitational acceleration [L2 Z-1 T-2 ~> m s-2]
  real,                 intent(in)  :: dz_subroundoff !< A minuscule thickness change [Z ~> m]
  real, dimension(SZI_(HI),SZJ_(HI)), &
                        intent(in)  :: bathyT !< The depth of the bathymetry [Z ~> m]
  type(EOS_type),       intent(in)  :: EOS !< Equation of state structure
  type(unit_scale_type), intent(in) :: US !< A dimensional unit scaling type
  logical,              intent(in) :: use_stanley_eos !< If true, turn on Stanley SGS T variance parameterization
  real, dimension(SZI_(HI),SZJ_(HI)), &
                        intent(inout) :: dpa !< The change in the pressure anomaly across the layer [R L2 T-2 ~> Pa]
  real, dimension(SZI_(HI),SZJ_(HI)), &
              optional, intent(inout) :: intz_dpa !< The integral through the thickness of the layer of
                                           !! the pressure anomaly relative to the anomaly at the
                                           !! top of the layer [R L2 Z T-2 ~> Pa m]
  real, dimension(SZIB_(HI),SZJ_(HI)), &
              optional, intent(inout) :: intx_dpa !< The integral in x of the difference between the
                                           !! pressure anomaly at the top and bottom of the layer
                                           !! divided by the x grid spacing [R L2 T-2 ~> Pa]
  real, dimension(SZI_(HI),SZJB_(HI)), &
              optional, intent(inout) :: inty_dpa !< The integral in y of the difference between the
                                           !! pressure anomaly at the top and bottom of the layer
                                           !! divided by the y grid spacing [R L2 T-2 ~> Pa]
  integer,    optional, intent(in)  :: MassWghtInterp !< A flag indicating whether and how to use
                                           !! mass weighting to interpolate T/S in integrals
  logical,    optional, intent(in)  :: use_inaccurate_form !< If true, uses an inaccurate form of
                                           !! density anomalies, as was used prior to March 2018.
  logical,    optional, intent(in)  :: MassWghtInterpVanOnly !< If true, does not do mass weighting
                                           !! of T/S unless one side smaller than h_nv (i.e. vanished)
  real,       optional, intent(in)  :: h_nv !< Nonvanished height [Z ~> m]
  real, dimension(HI%isd:HI%ied,HI%jsd:HI%jed), &
              optional, intent(in)  :: Z_0p !< The height at which the pressure is 0 [Z ~> m]

! This subroutine calculates (by numerical quadrature) integrals of
! pressure anomalies across layers, which are required for calculating the
! finite-volume form pressure accelerations in a Boussinesq model.  The one
! potentially dodgy assumption here is that rho_0 is used both in the denominator
! of the accelerations, and in the pressure used to calculated density (the
! latter being -z*rho_0*G_e).  These two uses could be separated if need be.
!
! It is assumed that the salinity and temperature profiles are linear in the
! vertical. The top and bottom values within each layer are provided and
! a linear interpolation is used to compute intermediate values.

  ! Local variables
                                             ! locations [C2 ~> degC2]
                                             ! locations [C S ~> degC ppt]
                                             ! locations [R ~> kg m-3]
                                             ! (used for inaccurate form) [R ~> kg m-3]
                                                ! locations [C2 ~> degC2]
                                                ! locations [C S ~> degC ppt]
                                                ! locations [S2 ~> ppt2]
                     ! with height at the 5 sub-column locations [R L2 T-2 ~> Pa]
                                    ! if at least one side vanished (0 or 1) [nondim]
                         ! of density anomalies.

end subroutine int_density_dz_generic_plm
module subroutine int_density_dz_generic_ppm(k, tv, T_t, T_b, S_t, S_b, e, &
                                      rho_ref, rho_0, G_e, dz_subroundoff, bathyT, HI, GV, EOS, US, use_stanley_eos, &
                                      dpa, intz_dpa, intx_dpa, inty_dpa, MassWghtInterp, Z_0p, &
                                      MassWghtInterpVanOnly, h_nv)
  integer,              intent(in)  :: k   !< Layer index to calculate integrals for
  type(hor_index_type), intent(in)  :: HI  !< Ocean horizontal index structures for the input arrays
  type(verticalGrid_type), intent(in) :: GV !< Vertical grid structure
  type(thermo_var_ptrs), intent(in) :: tv  !< Thermodynamic variables
  real, dimension(SZI_(HI),SZJ_(HI),SZK_(GV)), &
                        intent(in)  :: T_t !< Potential temperature at the cell top [C ~> degC]
  real, dimension(SZI_(HI),SZJ_(HI),SZK_(GV)), &
                        intent(in)  :: T_b !< Potential temperature at the cell bottom [C ~> degC]
  real, dimension(SZI_(HI),SZJ_(HI),SZK_(GV)), &
                        intent(in)  :: S_t !< Salinity at the cell top [S ~> ppt]
  real, dimension(SZI_(HI),SZJ_(HI),SZK_(GV)), &
                        intent(in)  :: S_b !< Salinity at the cell bottom [S ~> ppt]
  real, dimension(SZI_(HI),SZJ_(HI),SZK_(GV)+1), &
                        intent(in)  :: e   !< Height of interfaces [Z ~> m]
  real,                 intent(in)  :: rho_ref !< A mean density [R ~> kg m-3], that is
                                           !! subtracted out to reduce the magnitude of each of the integrals.
  real,                 intent(in)  :: rho_0 !< A density [R ~> kg m-3], that is used to calculate
                                           !! the pressure (as p~=-z*rho_0*G_e) used in the equation of state.
  real,                 intent(in)  :: G_e !< The Earth's gravitational acceleration [L2 Z-1 T-2 ~> m s-2]
  real,                 intent(in)  :: dz_subroundoff !< A minuscule thickness change [Z ~> m]
  real, dimension(SZI_(HI),SZJ_(HI)), &
                        intent(in)  :: bathyT !< The depth of the bathymetry [Z ~> m]
  type(EOS_type),       intent(in)  :: EOS !< Equation of state structure
  type(unit_scale_type), intent(in) :: US  !< A dimensional unit scaling type
  logical,              intent(in)  :: use_stanley_eos !< If true, turn on Stanley SGS T variance parameterization
  real, dimension(SZI_(HI),SZJ_(HI)), &
                        intent(inout) :: dpa !< The change in the pressure anomaly across the layer [R L2 T-2 ~> Pa]
  real, dimension(SZI_(HI),SZJ_(HI)), &
              optional, intent(inout) :: intz_dpa !< The integral through the thickness of the layer of
                                           !! the pressure anomaly relative to the anomaly at the
                                           !! top of the layer [R L2 Z T-2 ~> Pa m]
  real, dimension(SZIB_(HI),SZJ_(HI)), &
              optional, intent(inout) :: intx_dpa !< The integral in x of the difference between the
                                           !! pressure anomaly at the top and bottom of the layer
                                           !! divided by the x grid spacing [R L2 T-2 ~> Pa]
  real, dimension(SZI_(HI),SZJB_(HI)), &
              optional, intent(inout) :: inty_dpa !< The integral in y of the difference between the
                                           !! pressure anomaly at the top and bottom of the layer
                                           !! divided by the y grid spacing [R L2 T-2 ~> Pa]
  integer,    optional, intent(in)  :: MassWghtInterp !< A flag indicating whether and how to use
                                           !! mass weighting to interpolate T/S in integrals
  logical,    optional, intent(in)  :: MassWghtInterpVanOnly !< If true, does not do mass weighting
                                           !! of T/S unless one side smaller than h_nv (i.e. vanished)
  real,       optional, intent(in)  :: h_nv !< Nonvanished height [Z ~> m]

  real, dimension(HI%isd:HI%ied,HI%jsd:HI%jed), &
              optional, intent(in)  :: Z_0p !< The height at which the pressure is 0 [Z ~> m]

! This subroutine calculates (by numerical quadrature) integrals of
! pressure anomalies across layers, which are required for calculating the
! finite-volume form pressure accelerations in a Boussinesq model.  The one
! potentially dodgy assumption here is that rho_0 is used both in the denominator
! of the accelerations, and in the pressure used to calculated density (the
! latter being -z*rho_0*G_e).  These two uses could be separated if need be.
!
! It is assumed that the salinity and temperature profiles are parabolic in the
! vertical. The top and bottom values within each layer are provided and
! a parabolic interpolation is used to compute intermediate values.

  ! Local variables
                                             ! locations [C2 ~> degC2]
                                             ! locations [C S ~> degC ppt]
                                             ! locations [R ~> kg m-3]
                                                ! locations [C2 ~> degC2]
                                                ! locations [C S ~> degC ppt]
                                                ! locations [S2 ~> ppt2]
                  ! with height at the 5 sub-column locations [R L2 T-2 ~> Pa]
                                    ! if at least one side vanished (0 or 1) [nondim]

end subroutine int_density_dz_generic_ppm
module subroutine int_specific_vol_dp(T, S, p_t, p_b, alpha_ref, HI, EOS, US, &
                               dza, intp_dza, intx_dza, inty_dza, halo_size, &
                               bathyP, P_surf, dP_tiny, MassWghtInterp, &
                               MassWghtInterpVanOnly, p_nv)
  type(hor_index_type), intent(in)  :: HI  !< The horizontal index structure
  real, dimension(SZI_(HI),SZJ_(HI)), &
                        intent(in)  :: T   !< Potential temperature referenced to the surface [C ~> degC]
  real, dimension(SZI_(HI),SZJ_(HI)), &
                        intent(in)  :: S   !< Salinity [S ~> ppt]
  real, dimension(SZI_(HI),SZJ_(HI)), &
                        intent(in)  :: p_t !< Pressure at the top of the layer [R L2 T-2 ~> Pa]
  real, dimension(SZI_(HI),SZJ_(HI)), &
                        intent(in)  :: p_b !< Pressure at the bottom of the layer [R L2 T-2 ~> Pa]
  real,                 intent(in)  :: alpha_ref !< A mean specific volume that is subtracted out
                            !! to reduce the magnitude of each of the integrals [R-1 ~> m3 kg-1]
                            !! The calculation is mathematically identical with different values of
                            !! alpha_ref, but this reduces the effects of roundoff.
  type(EOS_type),       intent(in)  :: EOS !< Equation of state structure
  type(unit_scale_type), intent(in) :: US  !< A dimensional unit scaling type
  real, dimension(SZI_(HI),SZJ_(HI)), &
                        intent(inout) :: dza !< The change in the geopotential anomaly across
                            !! the layer [L2 T-2 ~> m2 s-2]
  real, dimension(SZI_(HI),SZJ_(HI)), &
              optional, intent(inout) :: intp_dza !< The integral in pressure through the layer of the
                            !! geopotential anomaly relative to the anomaly at the bottom of the
                            !! layer [R L4 T-4 ~> Pa m2 s-2]
  real, dimension(SZIB_(HI),SZJ_(HI)), &
              optional, intent(inout) :: intx_dza !< The integral in x of the difference between the
                            !! geopotential anomaly at the top and bottom of the layer divided by
                            !! the x grid spacing [L2 T-2 ~> m2 s-2]
  real, dimension(SZI_(HI),SZJB_(HI)), &
              optional, intent(inout) :: inty_dza !< The integral in y of the difference between the
                            !! geopotential anomaly at the top and bottom of the layer divided by
                            !! the y grid spacing [L2 T-2 ~> m2 s-2]
  integer,    optional, intent(in)  :: halo_size !< The width of halo points on which to calculate dza.
  real, dimension(SZI_(HI),SZJ_(HI)), &
              optional, intent(in)  :: bathyP  !< The pressure at the bathymetry [R L2 T-2 ~> Pa]
  real, dimension(SZI_(HI),SZJ_(HI)), &
              optional, intent(in)  :: P_surf !< The pressure at the ocean surface [R L2 T-2 ~> Pa]
  real,       optional, intent(in)  :: dP_tiny !< A minuscule pressure change with
                            !! the same units as p_t [R L2 T-2 ~> Pa]
  integer,    optional, intent(in)  :: MassWghtInterp !< A flag indicating whether and how to use
                                           !! mass weighting to interpolate T/S in integrals
  logical,    optional, intent(in)  :: MassWghtInterpVanOnly !< If true, does not do mass weighting
                                           !! of T/S unless one side smaller than h_nv (i.e. vanished)
  real,       optional, intent(in)  :: p_nv !< Nonvanished pressure [R L2 T-2 ~> Pa]


end subroutine int_specific_vol_dp
module subroutine int_spec_vol_dp_generic_pcm(T, S, p_t, p_b, alpha_ref, HI, EOS, US, dza, &
                                       intp_dza, intx_dza, inty_dza, halo_size, &
                                       bathyP, P_surf, dP_neglect, MassWghtInterp, &
                                       MassWghtInterpVanOnly, p_nv)
  type(hor_index_type), intent(in)  :: HI !< A horizontal index type structure.
  real, dimension(SZI_(HI),SZJ_(HI)), &
                        intent(in)  :: T  !< Potential temperature of the layer [C ~> degC]
  real, dimension(SZI_(HI),SZJ_(HI)), &
                        intent(in)  :: S  !< Salinity of the layer [S ~> ppt]
  real, dimension(SZI_(HI),SZJ_(HI)), &
                        intent(in)  :: p_t !< Pressure atop the layer [R L2 T-2 ~> Pa]
  real, dimension(SZI_(HI),SZJ_(HI)), &
                        intent(in)  :: p_b !< Pressure below the layer [R L2 T-2 ~> Pa]
  real,                 intent(in)  :: alpha_ref !< A mean specific volume that is subtracted out
                            !! to reduce the magnitude of each of the integrals [R-1 ~> m3 kg-1]
                            !! The calculation is mathematically identical with different values of
                            !! alpha_ref, but alpha_ref alters the effects of roundoff, and
                            !! answers do change.
  type(EOS_type),       intent(in)  :: EOS !< Equation of state structure
  type(unit_scale_type), intent(in) :: US !< A dimensional unit scaling type
  real, dimension(SZI_(HI),SZJ_(HI)), &
                        intent(inout) :: dza !< The change in the geopotential anomaly
                            !! across the layer [L2 T-2 ~> m2 s-2]
  real, dimension(SZI_(HI),SZJ_(HI)), &
              optional, intent(inout) :: intp_dza !< The integral in pressure through the layer of
                            !! the geopotential anomaly relative to the anomaly at the bottom of the
                            !! layer [R L4 T-4 ~> Pa m2 s-2]
  real, dimension(SZIB_(HI),SZJ_(HI)), &
              optional, intent(inout) :: intx_dza  !< The integral in x of the difference between
                            !! the geopotential anomaly at the top and bottom of the layer divided
                            !! by the x grid spacing [L2 T-2 ~> m2 s-2]
  real, dimension(SZI_(HI),SZJB_(HI)), &
              optional, intent(inout) :: inty_dza  !< The integral in y of the difference between
                            !! the geopotential anomaly at the top and bottom of the layer divided
                            !! by the y grid spacing [L2 T-2 ~> m2 s-2]
  integer,    optional, intent(in)  :: halo_size !< The width of halo points on which to calculate dza.
  real, dimension(SZI_(HI),SZJ_(HI)), &
              optional, intent(in)  :: bathyP !< The pressure at the bathymetry [R L2 T-2 ~> Pa]
  real, dimension(SZI_(HI),SZJ_(HI)), &
              optional, intent(in)  :: P_surf !< The pressure at the ocean surface [R L2 T-2 ~> Pa]
  real,       optional, intent(in)  :: dP_neglect !< A minuscule pressure change with
                                           !! the same units as p_t [R L2 T-2 ~> Pa]
  integer,    optional, intent(in)  :: MassWghtInterp !< A flag indicating whether and how to use
                                           !! mass weighting to interpolate T/S in integrals
  logical,    optional, intent(in)  :: MassWghtInterpVanOnly !< If true, does not do mass weighting
                                           !! of T/S unless one side smaller than h_nv (i.e. vanished)
  real,       optional, intent(in)  :: p_nv !< Nonvanished pressure [R L2 T-2 ~> Pa]

!   This subroutine calculates analytical and nearly-analytical integrals in
! pressure across layers of geopotential anomalies, which are required for
! calculating the finite-volume form pressure accelerations in a non-Boussinesq
! model.  There are essentially no free assumptions, apart from the use of
! Boole's rule to do the horizontal integrals, and from a truncation in the
! series for log(1-eps/1+eps) that assumes that |eps| < 0.34.

  ! Local variables
                                           ! locations [R-1 ~> m3 kg-1]
                     ! 5 sub-column locations [L2 T-2 ~> m2 s-2]
                                    ! if at least one side vanished (0 or 1) [nondim]

end subroutine int_spec_vol_dp_generic_pcm
module subroutine int_spec_vol_dp_generic_plm(T_t, T_b, S_t, S_b, p_t, p_b, alpha_ref, &
                             dP_neglect, bathyP, HI, EOS, US, dza, &
                             intp_dza, intx_dza, inty_dza, P_surf, MassWghtInterp, &
                             MassWghtInterpVanOnly, p_nv)
  type(hor_index_type), intent(in)  :: HI !< A horizontal index type structure.
  real, dimension(SZI_(HI),SZJ_(HI)), &
                        intent(in)  :: T_t  !< Potential temperature at the top of the layer [C ~> degC]
  real, dimension(SZI_(HI),SZJ_(HI)), &
                        intent(in)  :: T_b  !< Potential temperature at the bottom of the layer [C ~> degC]
  real, dimension(SZI_(HI),SZJ_(HI)), &
                        intent(in)  :: S_t  !< Salinity at the top the layer [S ~> ppt]
  real, dimension(SZI_(HI),SZJ_(HI)), &
                        intent(in)  :: S_b  !< Salinity at the bottom the layer [S ~> ppt]
  real, dimension(SZI_(HI),SZJ_(HI)), &
                        intent(in)  :: p_t !< Pressure atop the layer [R L2 T-2 ~> Pa]
  real, dimension(SZI_(HI),SZJ_(HI)), &
                        intent(in)  :: p_b !< Pressure below the layer [R L2 T-2 ~> Pa]
  real,                 intent(in)  :: alpha_ref !< A mean specific volume that is subtracted out
                            !! to reduce the magnitude of each of the integrals [R-1 ~> m3 kg-1]
                            !! The calculation is mathematically identical with different values of
                            !! alpha_ref, but alpha_ref alters the effects of roundoff, and
                            !! answers do change.
  real,                 intent(in)  :: dP_neglect !<!< A miniscule pressure change with
                                             !! the same units as p_t [R L2 T-2 ~> Pa]
  real, dimension(SZI_(HI),SZJ_(HI)), &
                        intent(in)  :: bathyP !< The pressure at the bathymetry [R L2 T-2 ~> Pa]
  type(EOS_type),       intent(in)  :: EOS !< Equation of state structure
  type(unit_scale_type), intent(in) :: US !< A dimensional unit scaling type
  real, dimension(SZI_(HI),SZJ_(HI)), &
                        intent(inout) :: dza !< The change in the geopotential anomaly
                            !! across the layer [L2 T-2 ~> m2 s-2]
  real, dimension(SZI_(HI),SZJ_(HI)), &
              optional, intent(inout) :: intp_dza !< The integral in pressure through the layer of
                            !! the geopotential anomaly relative to the anomaly at the bottom of the
                            !! layer [R L4 T-4 ~> Pa m2 s-2]
  real, dimension(SZIB_(HI),SZJ_(HI)), &
              optional, intent(inout) :: intx_dza  !< The integral in x of the difference between
                            !! the geopotential anomaly at the top and bottom of the layer divided
                            !! by the x grid spacing [L2 T-2 ~> m2 s-2]
  real, dimension(SZI_(HI),SZJB_(HI)), &
              optional, intent(inout) :: inty_dza  !< The integral in y of the difference between
                            !! the geopotential anomaly at the top and bottom of the layer divided
                            !! by the y grid spacing [L2 T-2 ~> m2 s-2]
  real, dimension(SZI_(HI),SZJ_(HI)), &
              optional, intent(in)  :: P_surf !< The pressure at the ocean surface [R L2 T-2 ~> Pa]
  integer,    optional, intent(in)  :: MassWghtInterp !< A flag indicating whether and how to use
                            !! mass weighting to interpolate T/S in integrals
  logical,    optional, intent(in)  :: MassWghtInterpVanOnly !< If true, does not do mass weighting
                                           !! of T/S unless one side smaller than h_nv (i.e. vanished)
  real,       optional, intent(in)  :: p_nv !< Nonvanished pressure [R L2 T-2 ~> Pa]

!   This subroutine calculates analytical and nearly-analytical integrals in
! pressure across layers of geopotential anomalies, which are required for
! calculating the finite-volume form pressure accelerations in a non-Boussinesq
! model.  There are essentially no free assumptions, apart from the use of
! Boole's rule to do the horizontal integrals, and from a truncation in the
! series for log(1-eps/1+eps) that assumes that |eps| < 0.34.

                                             ! locations [R-1 ~> m3 kg-1]

                     ! 5 sub-column locations [L2 T-2 ~> m2 s-2]
                                    ! if at least one side vanished (0 or 1) [nondim]

end subroutine int_spec_vol_dp_generic_plm
module subroutine diagnose_mass_weight_Z(z_t, z_b, bathyT, SSH, dz_neglect, MassWghtInterp, HI, &
                                  MassWt_u, MassWt_v, MassWghtInterpVanOnly, h_nv)
  type(hor_index_type), intent(in)  :: HI !< A horizontal index type structure.
  real, dimension(SZI_(HI),SZJ_(HI)), &
                        intent(in)  :: z_t !< Height at the top of the layer in depth units [Z ~> m]
  real, dimension(SZI_(HI),SZJ_(HI)), &
                        intent(in)  :: z_b !< Height at the bottom of the layer [Z ~> m]
  real, dimension(SZI_(HI),SZJ_(HI)), &
                        intent(in)  :: bathyT !< The depth of the bathymetry [Z ~> m]
  real, dimension(SZI_(HI),SZJ_(HI)), &
                        intent(in)  :: SSH !< The sea surface height [Z ~> m]
  real,                 intent(in)  :: dz_neglect !< A minuscule thickness change [Z ~> m]
  integer,              intent(in)  :: MassWghtInterp !< A flag indicating whether and how to use
                                           !! mass weighting to interpolate T/S in integrals
  real, dimension(SZIB_(HI),SZJ_(HI)), &
                        intent(inout) :: MassWt_u  !< The fractional mass weighting at u-points [nondim]
  real, dimension(SZI_(HI),SZJB_(HI)), &
                        intent(inout) :: MassWt_v  !< The fractional mass weighting at v-points [nondim]
  logical,    optional, intent(in)  :: MassWghtInterpVanOnly !< If true, does not do mass weighting
                                           !! of T/S unless one side smaller than h_nv (i.e. vanished)
  real,       optional, intent(in)  :: h_nv !< Nonvanished height [Z ~> m]

  ! Local variables

end subroutine diagnose_mass_weight_Z
module subroutine diagnose_mass_weight_p(p_t, p_b, bathyP, P_surf, dP_neglect, MassWghtInterp, HI, &
                                  MassWt_u, MassWt_v, MassWghtInterpVanOnly, p_nv)
  type(hor_index_type), intent(in)  :: HI !< A horizontal index type structure.
  real, dimension(SZI_(HI),SZJ_(HI)), &
                        intent(in)  :: p_t !< Pressure atop the layer [R L2 T-2 ~> Pa]
  real, dimension(SZI_(HI),SZJ_(HI)), &
                        intent(in)  :: p_b !< Pressure below the layer [R L2 T-2 ~> Pa]
  real,                 intent(in)  :: dP_neglect !<!< A miniscule pressure change with
                                           !! the same units as p_t [R L2 T-2 ~> Pa]
  real, dimension(SZI_(HI),SZJ_(HI)), &
                        intent(in)  :: bathyP !< The pressure at the bathymetry [R L2 T-2 ~> Pa]
  real, dimension(SZI_(HI),SZJ_(HI)), &
                        intent(in)  :: P_surf !< The pressure at the ocean surface [R L2 T-2 ~> Pa]
  integer,              intent(in)  :: MassWghtInterp !< A flag indicating whether and how to use
                                           !! mass weighting to interpolate T/S in integrals
  real, dimension(SZIB_(HI),SZJ_(HI)), &
                        intent(inout) :: MassWt_u  !< The fractional mass weighting at u-points [nondim]
  real, dimension(SZI_(HI),SZJB_(HI)), &
                        intent(inout) :: MassWt_v  !< The fractional mass weighting at v-points [nondim]
  logical,    optional, intent(in)  :: MassWghtInterpVanOnly !< If true, does not do mass weighting
                                           !! of T/S unless one side smaller than h_nv (i.e. vanished)
  real,       optional, intent(in)  :: p_nv !< Nonvanished pressure [R L2 T-2 ~> Pa]

  ! Local variables
                                    ! if at least one side vanished (0 or 1) [nondim]


end subroutine diagnose_mass_weight_p
module subroutine find_depth_of_pressure_in_cell(T_t, T_b, S_t, S_b, z_t, z_b, P_t, P_tgt, &
                       rho_ref, G_e, EOS, US, P_b, z_out, z_tol, frac_dp_bugfix)
  real,                  intent(in)  :: T_t !< Potential temperature at the cell top [C ~> degC]
  real,                  intent(in)  :: T_b !< Potential temperature at the cell bottom [C ~> degC]
  real,                  intent(in)  :: S_t !< Salinity at the cell top [S ~> ppt]
  real,                  intent(in)  :: S_b !< Salinity at the cell bottom [S ~> ppt]
  real,                  intent(in)  :: z_t !< Absolute height of top of cell [Z ~> m]   (Boussinesq ????)
  real,                  intent(in)  :: z_b !< Absolute height of bottom of cell [Z ~> m]
  real,                  intent(in)  :: P_t !< Anomalous pressure of top of cell, relative
                                            !! to g*rho_ref*z_t [R L2 T-2 ~> Pa]
  real,                  intent(in)  :: P_tgt !< Target pressure at height z_out, relative
                                            !! to g*rho_ref*z_out [R L2 T-2 ~> Pa]
  real,                  intent(in)  :: rho_ref !< Reference density with which calculation
                                            !! are anomalous to [R ~> kg m-3]
  real,                  intent(in)  :: G_e !< Gravitational acceleration [L2 Z-1 T-2 ~> m s-2]
  type(EOS_type),        intent(in)  :: EOS !< Equation of state structure
  type(unit_scale_type), intent(in)  :: US  !< A dimensional unit scaling type
  real,                  intent(out) :: P_b !< Pressure at the bottom of the cell [R L2 T-2 ~> Pa]
  real,                  intent(out) :: z_out !< Absolute depth at which anomalous pressure = p_tgt [Z ~> m]
  real,                  intent(in)  :: z_tol !< The tolerance in finding z_out [Z ~> m]
  logical,               intent(in)  :: frac_dp_bugfix !< If true, use bugfix in frac_dp_at_pos

  ! Local variables

end subroutine find_depth_of_pressure_in_cell
module subroutine avg_specific_vol(T, S, p_t, dp, HI, EOS, SpV_avg, halo_size)
  type(hor_index_type), intent(in)  :: HI  !< The horizontal index structure
  real, dimension(SZI_(HI),SZJ_(HI)), &
                        intent(in)  :: T   !< Potential temperature of the layer [C ~> degC]
  real, dimension(SZI_(HI),SZJ_(HI)), &
                        intent(in)  :: S   !< Salinity of the layer [S ~> ppt]
  real, dimension(SZI_(HI),SZJ_(HI)), &
                        intent(in)  :: p_t !< Pressure at the top of the layer [R L2 T-2 ~> Pa]
  real, dimension(SZI_(HI),SZJ_(HI)), &
                        intent(in)  :: dp  !< Pressure change in the layer [R L2 T-2 ~> Pa]
  type(EOS_type),       intent(in)  :: EOS !< Equation of state structure
  real, dimension(SZI_(HI),SZJ_(HI)), &
                        intent(inout) :: SpV_avg !< The vertical average specific volume
                                           !! in the layer [R-1 ~> m3 kg-1]
  integer,    optional, intent(in)  :: halo_size !< The number of halo points in which to work.

  ! Local variables

end subroutine avg_specific_vol
real module function frac_dp_at_pos(T_t, T_b, S_t, S_b, z_t, z_b, rho_ref, G_e, pos, EOS, frac_dp_bugfix)
  real,           intent(in)  :: T_t !< Potential temperature at the cell top [C ~> degC]
  real,           intent(in)  :: T_b !< Potential temperature at the cell bottom [C ~> degC]
  real,           intent(in)  :: S_t !< Salinity at the cell top [S ~> ppt]
  real,           intent(in)  :: S_b !< Salinity at the cell bottom [S ~> ppt]
  real,           intent(in)  :: z_t !< The geometric height at the top of the layer [Z ~> m]
  real,           intent(in)  :: z_b !< The geometric height at the bottom of the layer [Z ~> m]
  real,           intent(in)  :: rho_ref !< A mean density [R ~> kg m-3], that is subtracted out to
                                     !! reduce the magnitude of each of the integrals.
  real,           intent(in)  :: G_e !< The Earth's gravitational acceleration [L2 Z-1 T-2 ~> m s-2]
  real,           intent(in)  :: pos !< The fractional vertical position, 0 to 1 [nondim]
  type(EOS_type), intent(in)  :: EOS !< Equation of state structure
  logical,        intent(in)  :: frac_dp_bugfix !< If true, use bugfix in frac_dp_at_pos

  ! Local variables

end function frac_dp_at_pos
  end interface

end module MOM_density_integrals
