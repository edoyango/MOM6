! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> Functions for calculating interface heights, including free surface height.
module MOM_interface_heights

use MOM_density_integrals, only : int_specific_vol_dp, avg_specific_vol, int_density_dz
use MOM_debugging,     only : hchksum
use MOM_error_handler, only : MOM_error, FATAL
use MOM_EOS,           only : calculate_density, average_specific_vol, EOS_type, EOS_domain
use MOM_file_parser,   only : log_version
use MOM_grid,          only : ocean_grid_type
use MOM_unit_scaling,  only : unit_scale_type
use MOM_variables,     only : thermo_var_ptrs
use MOM_verticalGrid,  only : verticalGrid_type

implicit none ; private

#include <MOM_memory.h>

public find_eta, find_dz_for_eta, dz_to_thickness, thickness_to_dz, dz_to_thickness_simple
public calc_derived_thermo
public convert_MLD_to_ML_thickness
public find_rho_bottom, find_col_avg_SpV, find_col_mass

!> Calculates the heights of the free surface or all interfaces from layer thicknesses.
interface find_eta
  module procedure find_eta_2d, find_eta_3d
end interface find_eta

!> Calculates layer thickness in thickness units from geometric distance between the
!! interfaces around that layer in height units.
interface dz_to_thickness
  module procedure dz_to_thickness_tv, dz_to_thickness_EoS
end interface dz_to_thickness

!> Converts layer thickness in thickness units into the vertical distance between the
!! interfaces around a layer in height units.
interface thickness_to_dz
  module procedure thickness_to_dz_3d, thickness_to_dz_jslice
end interface thickness_to_dz


  interface
module subroutine find_dz_for_eta(h, tv, G, GV, US, dz_lay, halo_size)
  type(ocean_grid_type),                      intent(in)  :: G   !< The ocean's grid structure.
  type(verticalGrid_type),                    intent(in)  :: GV  !< The ocean's vertical grid structure.
  type(unit_scale_type),                      intent(in)  :: US  !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),  intent(in)  :: h   !< Layer thicknesses [H ~> m or kg m-2]
  type(thermo_var_ptrs),                      intent(in)  :: tv  !< A structure pointing to various
                                                                 !! thermodynamic variables.
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),  intent(out) :: dz_lay !< Height change across layers [Z ~> m]
  integer,                          optional, intent(in)  :: halo_size !< width of halo points on
                                                                 !! which to calculate eta.

  ! Local variables
                                  ! the units of thickness to layer mass [Z H-1 ~> nondim or m3 kg-1]
                                  ! rescaling factor derived from eta_to_m [T2 Z L-2 ~> s2 m-1]

end subroutine find_dz_for_eta
module subroutine find_eta_3d(h, tv, G, GV, US, eta, eta_bt, halo_size, dZref)
  type(ocean_grid_type),                      intent(in)  :: G   !< The ocean's grid structure.
  type(verticalGrid_type),                    intent(in)  :: GV  !< The ocean's vertical grid structure.
  type(unit_scale_type),                      intent(in)  :: US  !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),  intent(in)  :: h   !< Layer thicknesses [H ~> m or kg m-2]
  type(thermo_var_ptrs),                      intent(in)  :: tv  !< A structure pointing to various
                                                                 !! thermodynamic variables.
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)+1), intent(out) :: eta !< layer interface heights [Z ~> m]
  real, dimension(SZI_(G),SZJ_(G)), optional, intent(in)  :: eta_bt !< optional barotropic variable
                    !! that gives the "correct" free surface height (Boussinesq) or total water
                    !! column mass per unit area (non-Boussinesq).  This is used to dilate the layer
                    !! thicknesses when calculating interface heights [H ~> m or kg m-2].
                    !! In Boussinesq mode, eta_bt and G%bathyT use the same reference height.
  integer,                          optional, intent(in)  :: halo_size !< width of halo points on
                                                                 !! which to calculate eta.
  real,                             optional, intent(in)  :: dZref !< The difference in the
                    !! reference height between G%bathyT and eta [Z ~> m]. The default is 0.

  ! Local variables
                    ! dZ_ref is 0 unless the optional argument dZref is present.

end subroutine find_eta_3d
module subroutine find_eta_2d(h, tv, G, GV, US, eta, eta_bt, halo_size, dZref)
  type(ocean_grid_type),                      intent(in)  :: G   !< The ocean's grid structure.
  type(verticalGrid_type),                    intent(in)  :: GV  !< The ocean's vertical grid structure.
  type(unit_scale_type),                      intent(in)  :: US  !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),  intent(in)  :: h   !< Layer thicknesses [H ~> m or kg m-2]
  type(thermo_var_ptrs),                      intent(in)  :: tv  !< A structure pointing to various
                                                                 !! thermodynamic variables.
  real, dimension(SZI_(G),SZJ_(G)),           intent(out) :: eta !< free surface height relative to
                                                                 !! mean sea level (z=0) often [Z ~> m].
  real, dimension(SZI_(G),SZJ_(G)), optional, intent(in)  :: eta_bt !< optional barotropic
                    !! variable that gives the "correct" free surface height (Boussinesq) or total
                    !! water column mass per unit area (non-Boussinesq) [H ~> m or kg m-2].
                    !! In Boussinesq mode, eta_bt and G%bathyT use the same reference height.
  integer,                          optional, intent(in)  :: halo_size !< width of halo points on
                                                                 !! which to calculate eta.
  real,                             optional, intent(in)  :: dZref !< The difference in the
                    !! reference height between G%bathyT and eta [Z ~> m]. The default is 0.

  ! Local variables
                    ! dZ_ref is 0 unless the optional argument dZref is present.

end subroutine find_eta_2d
module subroutine calc_derived_thermo(tv, h, G, GV, US, halo, debug)
  type(ocean_grid_type),   intent(in)    :: G  !< The ocean's grid structure
  type(verticalGrid_type), intent(in)    :: GV !< The ocean's vertical grid structure
  type(unit_scale_type),   intent(in)    :: US !< A dimensional unit scaling type
  type(thermo_var_ptrs),   intent(inout) :: tv !< A structure pointing to various
                                               !! thermodynamic variables, some of
                                               !! which will be set here.
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in)    :: h  !< Layer thicknesses [H ~> m or kg m-2].
  integer,       optional, intent(in)    :: halo !< Width of halo within which to
                                               !! calculate thicknesses
  logical,       optional, intent(in)    :: debug !< If present and true, write debugging checksums
  ! Local variables
                                           ! state is used [R-1 ~> m3 kg-1]

end subroutine calc_derived_thermo
module subroutine find_col_avg_SpV(h, SpV_avg, tv, G, GV, US, halo_size)
  type(ocean_grid_type),    intent(in)    :: G    !< The ocean's grid structure
  type(verticalGrid_type),  intent(in)    :: GV   !< The ocean's vertical grid structure
  type(unit_scale_type),    intent(in)    :: US   !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                            intent(in)    :: h    !< Layer thicknesses [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJ_(G)), &
                            intent(inout) :: SpV_avg !< Column average specific volume [R-1 ~> m3 kg-1]
                                                  ! SpV_avg is intent inout to retain excess halo values.
  type(thermo_var_ptrs),    intent(in)    :: tv   !< Structure containing pointers to any available
                                                  !! thermodynamic fields.
  integer,        optional, intent(in)    :: halo_size !< width of halo points on which to work

  ! Local variables
                                ! the layer thicknesses [H R-1 ~> m4 kg-1 or m]

end subroutine find_col_avg_SpV
module subroutine find_col_mass(h, tv, G, GV, US, mass, p_bot, p_surf)
  type(ocean_grid_type),                      intent(in)  :: G    !< The ocean's grid structure.
  type(verticalGrid_type),                    intent(in)  :: GV   !< The ocean's vertical grid structure.
  type(unit_scale_type),                      intent(in)  :: US   !< A dimensional unit scaling type
  type(thermo_var_ptrs),                      intent(in)  :: tv   !< A structure pointing to various
                                                                  !! thermodynamic variables.
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),  intent(in)  :: h    !< Layer thicknesses [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJ_(G)),           intent(out) :: mass !< Integrated mass of the water column
                                                                  !! [R Z ~> kg m-2]
  real, dimension(SZI_(G),SZJ_(G)), optional, intent(out) :: p_bot  !< Bottom pressure = g * mass + psurf
                                                                    !! [R L2 T-2 ~> Pa]
  real, dimension(:,:),             optional, pointer     :: p_surf !< A pointer to surface pressure
                                                                    !! [R L2 T-2 ~> Pa]

  ! Local variables

end subroutine find_col_mass
module subroutine find_rho_bottom(G, GV, US, tv, h, dz, pres_int, dz_avg, j, Rho_bot, h_bot, k_bot)
  type(ocean_grid_type),    intent(in)  :: G    !< The ocean's grid structure
  type(verticalGrid_type),  intent(in)  :: GV   !< The ocean's vertical grid structure
  type(unit_scale_type),    intent(in)  :: US   !< A dimensional unit scaling type
  type(thermo_var_ptrs),    intent(in)  :: tv   !< Structure containing pointers to any available
                                                !! thermodynamic fields.
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                            intent(in)  :: h    !< Layer thicknesses [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZK_(GV)), &
                            intent(in)  :: dz   !< Height change across layers [Z ~> m]
  real, dimension(SZI_(G),SZK_(GV)+1), &
                            intent(in)  :: pres_int !< Pressure at each interface [R L2 T-2 ~> Pa]
  real, dimension(SZI_(G)), intent(in)  :: dz_avg !< The vertical distance over which to average [Z ~> m]
  integer,                  intent(in)  :: j    !< j-index of row to work on
  real, dimension(SZI_(G)), intent(out) :: Rho_bot  !< Near-bottom density [R ~> kg m-3].
  real, dimension(SZI_(G)), intent(out) :: h_bot !< Bottom boundary layer thickness [H ~> m or kg m-2]
  integer, dimension(SZI_(G)), intent(out) :: k_bot !< Bottom boundary layer top layer index

  ! Local variables
                              ! boundary layer [H R-1 ~> m4 kg-1 or m]
                              ! for [Z ~> m]
                              ! boundary layer [H ~> m or kg m-2]
                              ! boundary layer [C ~> degC]
                              ! boundary layer [S ~> ppt]
                              ! of the boundary layer [R L2 T-2 ~> Pa]
                              ! top of the boundary layer [R-1 ~> m3 kg-1]

end subroutine find_rho_bottom
module subroutine dz_to_thickness_tv(dz, tv, h, G, GV, US, halo_size)
  type(ocean_grid_type),   intent(in)    :: G  !< The ocean's grid structure
  type(verticalGrid_type), intent(in)    :: GV !< The ocean's vertical grid structure
  type(unit_scale_type),   intent(in)    :: US !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in)    :: dz !< Geometric layer thicknesses in height units [Z ~> m]
  type(thermo_var_ptrs),   intent(in)    :: tv !< A structure pointing to various
                                               !! thermodynamic variables
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(inout) :: h  !< Output thicknesses in thickness units [H ~> m or kg m-2].
                                               !! This is essentially intent out, but declared as intent
                                               !! inout to preserve any initialized values in halo points.
  integer,         optional, intent(in)  :: halo_size !< Width of halo within which to
                                               !! calculate thicknesses
  ! Local variables

end subroutine dz_to_thickness_tv
module subroutine dz_to_thickness_EOS(dz, Temp, Saln, EoS, h, G, GV, US, halo_size, p_surf)
  type(ocean_grid_type),   intent(in)    :: G  !< The ocean's grid structure
  type(verticalGrid_type), intent(in)    :: GV !< The ocean's vertical grid structure
  type(unit_scale_type),   intent(in)    :: US !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in)    :: dz !< Geometric layer thicknesses in height units [Z ~> m]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in)    :: Temp !< Input layer temperatures [C ~> degC]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in)    :: Saln !< Input layer salinities [S ~> ppt]
  type(EOS_type),          intent(in)    :: EoS  !< Equation of state structure
    real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(inout) :: h  !< Output thicknesses in thickness units [H ~> m or kg m-2].
                                               !! This is essentially intent out, but declared as intent
                                               !! inout to preserve any initialized values in halo points.
  integer,         optional, intent(in)  :: halo_size !< Width of halo within which to
                                               !! calculate thicknesses
  real, dimension(SZI_(G),SZJ_(G)), optional, intent(in)  :: p_surf !< Surface pressures [R L2 T-2 ~> Pa]
  ! Local variables
                                  ! iteration [R L2 T-2 ~> Pa]
                                  ! acceleration [H T2 R-1 L-2 ~> s2 m2 kg-1 or s2 m-1]

end subroutine dz_to_thickness_EOS
module subroutine dz_to_thickness_simple(dz, h, G, GV, US, halo_size, layer_mode)
  type(ocean_grid_type),   intent(in)    :: G  !< The ocean's grid structure
  type(verticalGrid_type), intent(in)    :: GV !< The ocean's vertical grid structure
  type(unit_scale_type),   intent(in)    :: US !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in)    :: dz !< Geometric layer thicknesses in height units [Z ~> m]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(inout) :: h  !< Output thicknesses in thickness units [H ~> m or kg m-2].
                                               !! This is essentially intent out, but declared as intent
                                               !! inout to preserve any initialized values in halo points.
  integer,         optional, intent(in)  :: halo_size !< Width of halo within which to
                                               !! calculate thicknesses
  logical,         optional, intent(in)  :: layer_mode !< If present and true, do the conversion that
                                               !! is appropriate in pure isopycnal layer mode with
                                               !! no state variables or equation of state.  Otherwise
                                               !! use a simple constant rescaling factor and avoid the
                                               !! use of GV%Rlay.
  ! Local variables
                      ! in pure isopycnal layered mode with no state variables or equation of state.

end subroutine dz_to_thickness_simple
module subroutine thickness_to_dz_3d(h, tv, dz, G, GV, US, halo_size)
  type(ocean_grid_type),   intent(in)    :: G  !< The ocean's grid structure
  type(verticalGrid_type), intent(in)    :: GV !< The ocean's vertical grid structure
  type(unit_scale_type),   intent(in)    :: US !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in)    :: h  !< Input thicknesses in thickness units [H ~> m or kg m-2].
  type(thermo_var_ptrs),   intent(in)    :: tv !< A structure pointing to various
                                               !! thermodynamic variables
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(inout) :: dz !< Geometric layer thicknesses in height units [Z ~> m]
                                               !! This is essentially intent out, but declared as intent
                                               !! inout to preserve any initialized values in halo points.
  integer,       optional, intent(in)    :: halo_size !< Width of halo within which to
                                               !! calculate thicknesses
  ! Local variables

end subroutine thickness_to_dz_3d
module subroutine thickness_to_dz_jslice(h, tv, dz, j, G, GV, halo_size)
  type(ocean_grid_type),   intent(in)    :: G  !< The ocean's grid structure
  type(verticalGrid_type), intent(in)    :: GV !< The ocean's vertical grid structure
   real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in)    :: h  !< Input thicknesses in thickness units [H ~> m or kg m-2].
  type(thermo_var_ptrs),   intent(in)    :: tv !< A structure pointing to various
                                               !! thermodynamic variables
  real, dimension(SZI_(G),SZK_(GV)), &
                           intent(inout) :: dz !< Geometric layer thicknesses in height units [Z ~> m]
                                               !! This is essentially intent out, but declared as intent
                                               !! inout to preserve any initialized values in halo points.
  integer,                 intent(in)    :: j  !< The second (j-) index of the input thicknesses to work with
  integer,       optional, intent(in)    :: halo_size !< Width of halo within which to
                                               !! calculate thicknesses
  ! Local variables

end subroutine thickness_to_dz_jslice
module subroutine convert_MLD_to_ML_thickness(MLD_in, h, h_MLD, tv, G, GV, halo)
  type(ocean_grid_type),   intent(in)    :: G  !< The ocean's grid structure
  type(verticalGrid_type), intent(in)    :: GV !< The ocean's vertical grid structure
  real, dimension(SZI_(G),SZJ_(G)), &
                           intent(in)    :: MLD_in !< Input mixed layer depth [Z ~> m].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in)    :: h  !< Layer thicknesses [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJ_(G)), &
                           intent(out)   :: h_MLD !< Thickness of water in the mixed layer [H ~> m or kg m-2]
  type(thermo_var_ptrs),   intent(in)    :: tv !< Structure containing pointers to any available
                                               !! thermodynamic fields.
  integer,       optional, intent(in)    :: halo !< Halo width over which to calculate frazil

  ! Local variables

end subroutine convert_MLD_to_ML_thickness
  end interface

end module MOM_interface_heights
