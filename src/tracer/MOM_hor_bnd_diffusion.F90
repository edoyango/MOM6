! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> Calculates and applies diffusive fluxes as a parameterization of horizontal mixing (non-neutral) by
!! mesoscale eddies near the top and bottom (to be implemented) boundary layers of the ocean.

module MOM_hor_bnd_diffusion

use MOM_cpu_clock,             only : cpu_clock_id, cpu_clock_begin, cpu_clock_end
use MOM_cpu_clock,             only : CLOCK_MODULE
use MOM_checksums,             only : hchksum
use MOM_domains,               only : pass_var
use MOM_diag_mediator,         only : diag_ctrl, time_type
use MOM_diag_mediator,         only : post_data, register_diag_field
use MOM_error_handler,         only : MOM_error, MOM_mesg, FATAL, is_root_pe
use MOM_file_parser,           only : get_param, log_version, param_file_type
use MOM_grid,                  only : ocean_grid_type
use MOM_remapping,             only : remapping_CS, initialize_remapping, reintegrate_column
use MOM_remapping,             only : extract_member_remapping_CS, remapping_core_h
use MOM_remapping,             only : remappingSchemesDoc, remappingDefaultScheme
use MOM_spatial_means,         only : global_mass_integral
use MOM_tracer_registry,       only : tracer_registry_type, tracer_type
use MOM_unit_scaling,          only : unit_scale_type
use MOM_variables,             only : vertvisc_type
use MOM_verticalGrid,          only : verticalGrid_type
use MOM_CVMix_KPP,             only : KPP_get_BLD, KPP_CS
use MOM_energetic_PBL,         only : energetic_PBL_get_MLD, energetic_PBL_CS
use MOM_diabatic_driver,       only : diabatic_CS, extract_diabatic_member
use MOM_io,                    only : stdout, stderr

implicit none ; private

public near_boundary_unit_tests, hor_bnd_diffusion, hor_bnd_diffusion_init
public boundary_k_range, hor_bnd_diffusion_end

! Private parameters to avoid doing string comparisons for bottom or top boundary layer
integer, public, parameter :: SURFACE = -1 !< Set a value that corresponds to the surface boundary
integer, public, parameter :: BOTTOM  = 1  !< Set a value that corresponds to the bottom boundary
#include <MOM_memory.h>

!> Sets parameters for horizontal boundary mixing module.
type, public :: hbd_CS ; private
  logical :: debug           !< If true, write verbose checksums for debugging.
  integer :: deg             !< Degree of polynomial reconstruction.
  integer :: hbd_nk          !< Maximum number of levels in the HBD grid [nondim]
  integer :: surface_boundary_scheme !< Which boundary layer scheme to use
                             !! 1. ePBL; 2. KPP
  logical :: limiter         !< Controls whether a flux limiter is applied in the
                             !! native grid (default is true).
  logical :: limiter_remap   !< Controls whether a flux limiter is applied in the
                             !! remapped grid (default is false).
  logical :: linear          !< If True, apply a linear transition at the base/top of the boundary.
                             !! The flux will be fully applied at k=k_min and zero at k=k_max.
  real    :: H_subroundoff   !< A thickness that is so small that it can be added to a thickness of
                             !! Angstrom or larger without changing it at the bit level [H ~> m or kg m-2].
                             !! If Angstrom is 0 or exceedingly small, this is negligible compared to 1e-17 m.
  ! HBD dynamic grids
  real,    allocatable, dimension(:,:,:) :: hbd_grd_u   !< HBD thicknesses at t-points adjacent to
                                                          !! u-points                     [H ~> m or kg m-2]
  real,    allocatable, dimension(:,:,:) :: hbd_grd_v   !< HBD thicknesses at t-points adjacent to
                                                          !! v-points (left and right)    [H ~> m or kg m-2]
  integer, allocatable, dimension(:,:)   :: hbd_u_kmax  !< Maximum vertical index in hbd_grd_u      [nondim]
  integer, allocatable, dimension(:,:)   :: hbd_v_kmax  !< Maximum vertical index in hbd_grd_v      [nondim]
  type(remapping_CS)              :: remap_CS          !< Control structure to hold remapping configuration.
  type(KPP_CS),           pointer :: KPP_CSp => NULL() !< KPP control structure needed to get BLD.
  type(energetic_PBL_CS), pointer :: energetic_PBL_CSp => NULL()  !< ePBL control structure needed to get BLD.
  type(diag_ctrl), pointer :: diag => NULL() !< A structure that is used to
                                             !! regulate the timing of diagnostic output.
end type hbd_CS

! This include declares and sets the variable "version".
#include "version_variable.h"
character(len=40) :: mdl = "MOM_hor_bnd_diffusion" !< Name of this module
integer :: id_clock_hbd                            !< CPU clock for hbd


  interface
logical module function hor_bnd_diffusion_init(Time, G, GV, US, param_file, diag, diabatic_CSp, CS)
  type(time_type), target,          intent(in)    :: Time          !< Time structure
  type(ocean_grid_type),            intent(in)    :: G             !< Grid structure
  type(verticalGrid_type),          intent(in)    :: GV            !< ocean vertical grid structure
  type(unit_scale_type),            intent(in)    :: US            !< A dimensional unit scaling type
  type(param_file_type),            intent(in)    :: param_file    !< Parameter file structure
  type(diag_ctrl), target,          intent(inout) :: diag          !< Diagnostics control structure
  type(diabatic_CS),                pointer       :: diabatic_CSp  !< KPP control structure needed to get BLD
  type(hbd_CS),                     pointer       :: CS            !< Horizontal boundary mixing control structure

  ! local variables

end function hor_bnd_diffusion_init
module subroutine hor_bnd_diffusion(G, GV, US, h, Coef_x, Coef_y, dt, Reg, visc, CS)
  type(ocean_grid_type),                        intent(inout) :: G      !< Grid type
  type(verticalGrid_type),                      intent(in)    :: GV     !< ocean vertical grid structure
  type(unit_scale_type),                        intent(in)    :: US     !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),    intent(in)    :: h      !< Layer thickness [H ~> m or kg m-2]
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)+1), intent(in)    :: Coef_x !< dt * Kh * dy / dx at u-points [L2 ~> m2]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)+1), intent(in)    :: Coef_y !< dt * Kh * dx / dy at v-points [L2 ~> m2]
  real,                                         intent(in)    :: dt     !< Tracer time step * I_numitts
                                                                        !! (I_numitts in tracer_hordiff) [T ~> s]
  type(tracer_registry_type),                   pointer       :: Reg    !< Tracer registry
  type(vertvisc_type),                          intent(in)    :: visc   !< Structure with vertical viscosities,
                                                                        !! boundary layer properties and related fields
  type(hbd_CS),                                 pointer       :: CS     !< Control structure for this module

  ! Local variables
                                                            !! [conc H L2 ~> conc m3 or conc kg]
                                                            !! [conc H L2 ~> conc m3 or conc kg]
                                                            !! [conc H L2 ~> conc m3 or conc kg]
                                                            !! [H conc T-1 ~> m conc s-1 or kg m-2 conc s-1],
                                                            !! then converted to [conc T-1 ~> conc s-1].
                                                            ! For temperature these units are
                                                            ! [C H T-1 ~> degC m s-1 or degC kg m-2 s-1] and
                                                            ! then [C T-1 ~> degC s-1].
                                                            !! [H conc T-1 ~> m conc s-1 or kg m-2 conc s-1].
                                                            !! For temperature these units are
                                                            !! [C H T-1 ~> degC m s-1 or degC kg m-2 s-1].
                                                            !! only used to compute tendencies [conc].

end subroutine hor_bnd_diffusion
module subroutine hbd_grid(boundary, G, GV, hbl, h, CS)
  integer,                 intent(in   ) :: boundary !< Which boundary layer SURFACE or BOTTOM       [nondim]
  type(ocean_grid_type),   intent(inout) :: G    !< Grid type
  type(verticalGrid_type), intent(in)    :: GV   !< ocean vertical grid structure
  real, dimension(SZI_(G),SZJ_(G)), &
                           intent(in)    :: hbl  !< Boundary layer depth                   [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in)    :: h    !< Layer thickness in the native grid     [H ~> m or kg m-2]
  type(hbd_CS),            pointer       :: CS   !< Horizontal diffusion control structure

  ! Local variables

  ! reset arrays
end subroutine hbd_grid
real module function harmonic_mean(h1,h2)
  real :: h1 !< Scalar quantity [arbitrary]
  real :: h2 !< Scalar quantity [arbitrary]
end function harmonic_mean
integer module function  find_minimum(x, s, e)
  integer, intent(in) :: s              !< start index
  integer, intent(in) :: e              !< end index
  real, dimension(e), intent(in) :: x   !< 1D array to be checked [arbitrary]

  ! local variables

end function find_minimum
module subroutine swap(a, b)
  real, intent(inout) :: a  !< First value to be swapped [arbitrary]
  real, intent(inout) :: b  !< Second value to be swapped [arbitrary]

  ! local variables

end subroutine swap
module subroutine sort(x, n)
  integer,             intent(in   ) :: n        !< Number of points in the array
  real, dimension(n),  intent(inout) :: x        !< 1D array to be sorted [arbitrary]

  ! local variables

end subroutine sort
module subroutine unique(val, n, val_unique, val_max)
  integer,                         intent(in   ) :: n          !< Number of points in the array.
  real, dimension(n),              intent(in   ) :: val        !< 1D array to be checked [arbitrary]
  real, dimension(:), allocatable, intent(inout) :: val_unique !< Returned 1D array with unique values [arbitrary]
  real,                  optional, intent(in   ) :: val_max    !< sets the maximum value in val_unique to
                                                               !! this value [arbitrary]
  ! local variables

end subroutine unique
module subroutine merge_interfaces(nk, h_L, h_R, hbl_L, hbl_R, H_subroundoff, h)
  integer,                         intent(in   ) :: nk     !< Number of layers                        [nondim]
  real, dimension(nk),             intent(in   ) :: h_L    !< Layer thicknesses in the left column    [H ~> m or kg m-2]
  real, dimension(nk),             intent(in   ) :: h_R    !< Layer thicknesses in the right column   [H ~> m or kg m-2]
  real,                            intent(in   ) :: hbl_L  !< Thickness of the boundary layer in the left column
                                                           !!                                         [H ~> m or kg m-2]
  real,                            intent(in   ) :: hbl_R  !< Thickness of the boundary layer in the right column
                                                           !!                                         [H ~> m or kg m-2]
  real,                            intent(in   ) :: H_subroundoff !< GV%H_subroundoff                 [H ~> m or kg m-2]
  real, dimension(:), allocatable, intent(inout) :: h     !< Combined thicknesses                     [H ~> m or kg m-2]

  ! Local variables
                                                 !! plus hbl_L and hbl_R [H ~> m or kg m-2]
                                                 !! hbl_L and hbl_R [H ~> m or kg m-2]

end subroutine merge_interfaces
module subroutine flux_limiter(F_layer, area_L, area_R, phi_L, phi_R, h_L, h_R)
  real, intent(inout) :: F_layer !< Tracer flux to be checked [H L2 conc ~> m3 conc]
  real, intent(in) :: area_L     !< Area of left cell [L2 ~> m2]
  real, intent(in) :: area_R     !< Area of right cell [L2 ~> m2]
  real, intent(in) :: h_L        !< Thickness of left cell [H ~> m or kg m-2]
  real, intent(in) :: h_R        !< Thickness of right cell [H ~> m or kg m-2]
  real, intent(in) :: phi_L      !< Tracer concentration in the left cell [conc]
  real, intent(in) :: phi_R      !< Tracer concentration in the right cell [conc]

  ! local variables
  ! limit the flux to 0.2 of the tracer *gradient*
  ! Why 0.2?
  !  t=0         t=inf
  !   0           .2
  ! 0 1 0       .2.2.2
  !   0           .2
  !
end subroutine flux_limiter
module subroutine boundary_k_range(boundary, nk, h, hbl, k_top, zeta_top, k_bot, zeta_bot)
  integer,             intent(in   ) :: boundary !< SURFACE or BOTTOM                       [nondim]
  integer,             intent(in   ) :: nk       !< Number of layers                        [nondim]
  real, dimension(nk), intent(in   ) :: h        !< Layer thicknesses of the column         [H ~> m or kg m-2]
  real,                intent(in   ) :: hbl      !< Thickness of the boundary layer         [H ~> m or kg m-2]
                                                 !! If surface, with respect to zbl_ref = 0.
                                                 !! If bottom, with respect to zbl_ref = SUM(h)
  integer,             intent(  out) :: k_top    !< Index of the first layer within the boundary
  real,                intent(  out) :: zeta_top !< Distance from the top of a layer to the intersection of the
                                                 !! top extent of the boundary layer (0 at top, 1 at bottom)  [nondim]
  integer,             intent(  out) :: k_bot    !< Index of the last layer within the boundary
  real,                intent(  out) :: zeta_bot !< Distance of the lower layer to the boundary layer depth
                                                 !! (0 at top, 1 at bottom)  [nondim]
  ! Local variables

  ! Surface boundary layer
end subroutine boundary_k_range
module subroutine fluxes_layer_method(boundary, ke, hbl_L, hbl_R, h_L, h_R, phi_L, phi_R, &
                               khtr_u, F_layer, area_L, area_R, nk, dz_top, CS)

  integer,              intent(in   ) :: boundary !< Which boundary layer SURFACE or BOTTOM           [nondim]
  integer,              intent(in   ) :: ke       !< Number of layers in the native grid              [nondim]
  real,                 intent(in   ) :: hbl_L    !< Thickness of the boundary boundary
                                                  !! layer (left)                           [H ~> m or kg m-2]
  real,                 intent(in   ) :: hbl_R    !< Thickness of the boundary boundary
                                                  !! layer (right)                          [H ~> m or kg m-2]
  real, dimension(ke),  intent(in   ) :: h_L      !< Thicknesses in the native grid (left)  [H ~> m or kg m-2]
  real, dimension(ke),  intent(in   ) :: h_R      !< Thicknesses in the native grid (right) [H ~> m or kg m-2]
  real, dimension(ke),  intent(in   ) :: phi_L    !< Tracer values in the native grid (left)            [conc]
  real, dimension(ke),  intent(in   ) :: phi_R    !< Tracer values in the native grid (right)           [conc]
  real, dimension(ke+1),intent(in   ) :: khtr_u   !< Horizontal diffusivities times the time step
                                                  !! at a velocity point and vertical interfaces    [L2 ~> m2]
  real, dimension(ke),  intent(  out) :: F_layer  !< Layerwise diffusive flux at U- or V-point
                                                  !! in the native grid                 [H L2 conc ~> m3 conc]
  real,                 intent(in   ) :: area_L   !< Area of the horizontal grid (left)             [L2 ~> m2]
  real,                 intent(in   ) :: area_R   !< Area of the horizontal grid (right)            [L2 ~> m2]
  integer,              intent(in   ) :: nk       !< Number of layers in the HBD grid                 [nondim]
  real, dimension(nk),  intent(in   ) :: dz_top   !< The HBD z grid                         [H ~> m or kg m-2]
  type(hbd_CS),         pointer       :: CS       !< Horizontal diffusion control structure

  ! Local variables
                                     !! The harmonic mean is used to avoid zero values      [H ~> m or kg m-2]
                                     !! layer depth in the native grid                                [nondim]
                                     !! layer depth in the native grid                                [nondim]

end subroutine fluxes_layer_method
logical module function near_boundary_unit_tests( verbose )
  logical,               intent(in) :: verbose !< If true, output additional information for debugging unit tests

  ! Local variables

end function near_boundary_unit_tests
logical module function test_layer_fluxes(verbose, nk, test_name, F_calc, F_ans)
  logical,                    intent(in) :: verbose   !< If true, write results to stdout
  character(len=80),          intent(in) :: test_name !< Brief description of the unit test
  integer,                    intent(in) :: nk        !< Number of layers
  real, dimension(nk),        intent(in) :: F_calc    !< Fluxes or other quantity from the algorithm [arbitrary]
  real, dimension(nk),        intent(in) :: F_ans     !< Expected value calculated by hand [arbitrary]
  ! Local variables

end function test_layer_fluxes
logical module function test_boundary_k_range(k_top, zeta_top, k_bot, zeta_bot, k_top_ans, zeta_top_ans,&
                                       k_bot_ans, zeta_bot_ans, test_name, verbose)
  integer :: k_top               !< Index of cell containing top of boundary
  real    :: zeta_top            !< Fractional position in the cell of the top boundary [nondim]
  integer :: k_bot               !< Index of cell containing bottom of boundary
  real    :: zeta_bot            !< Fractional position in the cell of the bottom boundary [nondim]
  integer :: k_top_ans           !< Expected index of cell containing top of boundary
  real    :: zeta_top_ans        !< Expected fractional position of the top boundary [nondim]
  integer :: k_bot_ans           !< Expected index of cell containing bottom of boundary
  real    :: zeta_bot_ans        !< Expected fractional position of the bottom boundary [nondim]
  character(len=80) :: test_name !< Name of the unit test
  logical :: verbose             !< If true always print output

end function test_boundary_k_range
module subroutine hbd_grid_test(boundary, hbl_L, hbl_R, h_L, h_R, CS)
  integer,                 intent(in) :: boundary !< Which boundary layer SURFACE or BOTTOM    [nondim]
  real,                    intent(in) :: hbl_L    !< Boundary layer depth, left                [H ~> m or kg m-2]
  real,                    intent(in) :: hbl_R    !< Boundary layer depth, right               [H ~> m or kg m-2]
  real, dimension(2),      intent(in) :: h_L      !< Layer thickness in the native grid, left  [H ~> m or kg m-2]
  real, dimension(2),      intent(in) :: h_R      !< Layer thickness in the native grid, right [H ~> m or kg m-2]
  type(hbd_CS),            pointer    :: CS       !< Horizontal diffusion control structure

  ! Local variables

  ! reset arrays
end subroutine hbd_grid_test
module subroutine hor_bnd_diffusion_end(CS)
  type(hbd_CS), pointer :: CS  !< Horizontal boundary diffusion control structure

end subroutine hor_bnd_diffusion_end
  end interface

end module MOM_hor_bnd_diffusion
