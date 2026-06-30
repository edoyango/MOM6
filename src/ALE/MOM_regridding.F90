! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> Generates vertical grids as part of the ALE algorithm
module MOM_regridding

use MOM_error_handler, only : MOM_error, FATAL, WARNING, NOTE, assert
use MOM_file_parser,   only : param_file_type, get_param, log_param
use MOM_io,            only : file_exists, field_exists, field_size, MOM_read_data
use MOM_io,            only : read_variable
use MOM_io,            only : vardesc, var_desc, SINGLE_FILE
use MOM_io,            only : MOM_netCDF_file, MOM_field
use MOM_io,            only : create_MOM_file, MOM_write_field
use MOM_io,            only : verify_variable_units, slasher
use MOM_unit_scaling,  only : unit_scale_type
use MOM_variables,     only : ocean_grid_type, thermo_var_ptrs
use MOM_verticalGrid,  only : verticalGrid_type
use MOM_EOS,           only : EOS_type, calculate_density
use MOM_domains,       only : max_across_PEs, pass_var
use MOM_string_functions, only : uppercase, extractWord, extract_integer, extract_real

use MOM_remapping, only : remapping_CS
use regrid_consts, only : state_dependent, coordinateUnits
use regrid_consts, only : coordinateMode, DEFAULT_COORDINATE_MODE
use regrid_consts, only : REGRIDDING_LAYER, REGRIDDING_ZSTAR
use regrid_consts, only : REGRIDDING_RHO, REGRIDDING_SIGMA
use regrid_consts, only : REGRIDDING_ARBITRARY, REGRIDDING_SIGMA_SHELF_ZSTAR
use regrid_consts, only : REGRIDDING_HYCOM1, REGRIDDING_HYBGEN, REGRIDDING_ADAPTIVE
use regrid_interp, only : interp_CS_type
use regrid_interp, only : set_interp_scheme, set_interp_extrap, set_interp_answer_date

use coord_zlike,  only : zlike_CS
use coord_zlike,  only : init_coord_zlike, set_zlike_params, build_zstar_column, end_coord_zlike
use coord_sigma,  only : sigma_CS
use coord_sigma,  only : init_coord_sigma, set_sigma_params, build_sigma_column, end_coord_sigma
use coord_rho,    only : init_coord_rho, rho_CS, set_rho_params, build_rho_column, end_coord_rho
use coord_rho,    only : old_inflate_layers_1d
use coord_hycom,  only : hycom_CS
use coord_hycom,  only : init_coord_hycom, set_hycom_params, build_hycom1_column, end_coord_hycom
use coord_hycom,  only : init_3d_coord_hycom
use coord_adapt,  only : adapt_CS
use coord_adapt,  only : init_coord_adapt, set_adapt_params, build_adapt_column, end_coord_adapt
use MOM_hybgen_regrid, only : hybgen_regrid, hybgen_regrid_CS, init_hybgen_regrid, end_hybgen_regrid
use MOM_hybgen_regrid, only : write_Hybgen_coord_file

implicit none ; private

#include <MOM_memory.h>

! A note on unit descriptions in comments: MOM6 uses units that can be rescaled for dimensional
! consistency testing. These are noted in comments with units like Z, H, L, and T, along with
! their mks counterparts with notation like "a velocity [Z T-1 ~> m s-1]".  If the units
! vary with the Boussinesq approximation, the Boussinesq variant is given first.

!> Regridding control structure
type, public :: regridding_CS ; private

  !> This array is set by function setCoordinateResolution()
  !! It contains the "resolution" or delta coordinate of the target
  !! coordinate.  It has the units of the target coordinate, e.g.
  !! [Z ~> m] for z*, [nondim] for sigma, etc.
  real, dimension(:), allocatable :: coordinateResolution

  !> This is a scaling factor that restores coordinateResolution to values in
  !! the natural units for output, perhaps [nondim]
  real :: coord_scale = 1.0

  !> This array is set by function set_target_densities()
  !! This array is the nominal coordinate of interfaces and is the
  !! running sum of coordinateResolution, in [R ~> kg m-3]. i.e.
  !!  target_density(k+1) = coordinateResolution(k) + coordinateResolution(k)
  !! It is only used in "rho" or "Hycom" mode.
  real, dimension(:), allocatable :: target_density

  !> A flag to indicate that the target_density arrays has been filled with data.
  logical :: target_density_set = .false.

  !> Nominal HYCOM1 3D near-surface resolution [Z ~> m]
  real, allocatable, dimension(:,:,:) :: coordinateResolution_3d

  !> Nominal HYCOM1 3D density of interfaces [R ~> kg m-3]
  real, allocatable, dimension(:,:,:) :: target_density_3d

  !> This array is set by function set_regrid_max_depths()
  !! It specifies the maximum depth that every interface is allowed to take [H ~> m or kg m-2].
  real, dimension(:), allocatable :: max_interface_depths

  !> This array is set by function set_regrid_max_thickness()
  !! It specifies the maximum depth that every interface is allowed to take [H ~> m or kg m-2].
  real, dimension(:), allocatable :: max_layer_thickness

  integer :: nk !< Number of layers/levels in generated grid

  !> Indicates which grid to use in the vertical (z*, sigma, target interface
  !! densities)
  integer :: regridding_scheme

  !> Interpolation control structure
  type(interp_CS_type) :: interp_CS

  !> Minimum thickness allowed when building the new grid through regridding [H ~> m or kg m-2].
  real :: min_thickness

  !> If true, call adjust_interface_motion() after initial grid generation
  logical :: use_adjust_interface_motion

  !> Reference pressure for potential density calculations [R L2 T-2 ~> Pa]
  real :: ref_pressure = 2.e7

  !> If true, always pass through the depth-based time filtering that uses CS%old_grid_weight
  !! If false, allows bypassing of the call if CS%old_grid_weight==0
  logical :: use_depth_based_time_filter

  !> Weight given to old coordinate when blending between new and old grids [nondim]
  !! Used only below depth_of_time_filter_shallow, with a cubic variation
  !! from zero to full effect between depth_of_time_filter_shallow and
  !! depth_of_time_filter_deep.
  real :: old_grid_weight = 0.

  !> Depth above which no time-filtering of grid is applied [H ~> m or kg m-2]
  real :: depth_of_time_filter_shallow = 0.

  !> Depth below which time-filtering of grid is applied at full effect [H ~> m or kg m-2]
  real :: depth_of_time_filter_deep = 0.

  !> Fraction (between 0 and 1) of compressibility to add to potential density
  !! profiles when interpolating for target grid positions [nondim]
  real :: compressibility_fraction = 0.

  !> If true, each interface is given a maximum depth based on a rescaling of
  !! the indexing of coordinateResolution.
  logical :: set_maximum_depths = .false.

  !> If true, integrate for interface positions from the top downward.
  !! If false, integrate from the bottom upward, as does the rest of the model.
  logical :: integrate_downward_for_e = .true.

  !> The vintage of the order of arithmetic and expressions to use for remapping.
  !! Values below 20190101 recover the remapping answers from 2018.
  !! Higher values use more robust forms of the same remapping expressions.
  integer :: remap_answer_date = 99991231

  logical :: use_hybgen_unmix = .false.  !< If true, use the hybgen unmixing code before remapping

  type(zlike_CS),  pointer :: zlike_CS  => null() !< Control structure for z-like coordinate generator
  type(sigma_CS),  pointer :: sigma_CS  => null() !< Control structure for sigma coordinate generator
  type(rho_CS),    pointer :: rho_CS    => null() !< Control structure for rho coordinate generator
  type(hycom_CS),  pointer :: hycom_CS  => null() !< Control structure for hybrid coordinate generator
  type(adapt_CS),  pointer :: adapt_CS  => null() !< Control structure for adaptive coordinate generator
  type(hybgen_regrid_CS), pointer :: hybgen_CS => NULL() !< Control structure for hybgen regridding

end type

! The following routines are visible to the outside world
public initialize_regridding, end_regridding, regridding_main
public regridding_preadjust_reqs, convective_adjustment
public inflate_vanished_layers_old, check_grid_column
public set_regrid_params, get_regrid_size, write_regrid_file
public uniformResolution, setCoordinateResolution
public set_target_densities_from_GV, set_target_densities
public set_regrid_max_depths, set_regrid_max_thickness
public getCoordinateResolution, getCoordinateInterfaces
public getCoordinateUnits, getCoordinateShortName, getStaticThickness
public DEFAULT_COORDINATE_MODE
public set_h_neglect, set_dz_neglect
public get_zlike_CS, get_sigma_CS, get_rho_CS

!> Documentation for coordinate options
character(len=*), parameter, public :: regriddingCoordinateModeDoc = &
                 " LAYER - Isopycnal or stacked shallow water layers\n"//&
                 " ZSTAR, Z* - stretched geopotential z*\n"//&
                 " SIGMA_SHELF_ZSTAR - stretched geopotential z* ignoring shelf\n"//&
                 " SIGMA - terrain following coordinates\n"//&
                 " RHO   - continuous isopycnal\n"//&
                 " HYCOM1 - HyCOM-like hybrid coordinate\n"//&
                 " HYBGEN - Hybrid coordinate from the Hycom hybgen code\n"//&
                 " ADAPTIVE - optimize for smooth neutral density surfaces"

!> Documentation for regridding interpolation schemes
character(len=*), parameter, public :: regriddingInterpSchemeDoc = &
                 " P1M_H2     (2nd-order accurate)\n"//&
                 " P1M_H4     (2nd-order accurate)\n"//&
                 " P1M_IH4    (2nd-order accurate)\n"//&
                 " PLM        (2nd-order accurate)\n"//&
                 " PPM_CW     (3rd-order accurate)\n"//&
                 " PPM_H4     (3rd-order accurate)\n"//&
                 " PPM_IH4    (3rd-order accurate)\n"//&
                 " P3M_IH4IH3 (4th-order accurate)\n"//&
                 " P3M_IH6IH5 (4th-order accurate)\n"//&
                 " PQM_IH4IH3 (4th-order accurate)\n"//&
                 " PQM_IH6IH5 (5th-order accurate)"

!> Default interpolation scheme
character(len=*), parameter, public :: regriddingDefaultInterpScheme = "P1M_H2"
!> Default mode for boundary extrapolation
logical, parameter, public :: regriddingDefaultBoundaryExtrapolation = .false.
!> Default minimum thickness for some coordinate generation modes [m]
real, parameter, public :: regriddingDefaultMinThickness = 1.e-3

!> Maximum length of parameters
integer, parameter :: MAX_PARAM_LENGTH = 120

#undef __DO_SAFETY_CHECKS__


  interface
module subroutine initialize_regridding(CS, G, GV, US, max_depth, param_file, mdl, &
                                 coord_mode, param_prefix, param_suffix)
  type(regridding_CS),        intent(inout) :: CS  !< Regridding control structure
  type(ocean_grid_type),      intent(in)    :: G   !< Ocean grid structure
  type(verticalGrid_type),    intent(in)    :: GV  !< Ocean vertical grid structure
  type(unit_scale_type),      intent(in)    :: US  !< A dimensional unit scaling type
  real,                       intent(in)    :: max_depth  !< The maximum depth of the ocean [Z ~> m].
  type(param_file_type),      intent(in)    :: param_file !< Parameter file
  character(len=*),           intent(in)    :: mdl        !< Name of calling module.
  character(len=*),           intent(in)    :: coord_mode !< Coordinate mode
  character(len=*),           intent(in)    :: param_prefix !< String to prefix to parameter names.
                                                            !! If empty, causes main model parameters to be used.
  character(len=*),           intent(in)    :: param_suffix !< String to append to parameter names.

  ! Local variables
                        ! maximum_depth is large [m] (not in Z).
                                            ! or [Z ~> m] or [H ~> m or kg m-2] or [R ~> kg m-3] or other units.
                                            ! or [Z ~> m] or [H ~> m or kg m-2] or [R ~> kg m-3] or other units.
                                            ! or [Z ~> m] or [H ~> m or kg m-2] or [R ~> kg m-3] or other units.
                                            ! units depending on the coordinate
                                            ! [H ~> m or kg m-2] or other units
  ! Thicknesses [m] that give level centers approximately corresponding to table 2 of WOA09
  ! These are approximate because the WOA09 depths are not smoothly spaced. Levels
  ! 1, 4, 5, 9, 12, 24, and 36 are 2.5, 2.5, 1.25 12.5, 37.5 and 62.5 m deeper than WOA09
  ! but all others are identical.
  ! These are the actual spacings [m] between WOA09 depths which, if used for layer thickness, places
  ! the interfaces at the WOA09 depths.
  ! These are the spacings [m] between WOA23 depths from table 3 of
  ! https://www.ncei.noaa.gov/data/oceans/woa/WOA13/DOC/woa13documentation.pdf

end subroutine initialize_regridding
module subroutine end_regridding(CS)
  type(regridding_CS), intent(inout) :: CS !< Regridding control structure

end subroutine end_regridding
module subroutine regridding_main( remapCS, CS, G, GV, US, h, tv, h_new, dzInterface, &
                            frac_shelf_h, PCM_cell)
!------------------------------------------------------------------------------
! This routine takes care of (1) building a new grid and (2) remapping between
! the old grid and the new grid. The creation of the new grid can be based
! on z coordinates, target interface densities, sigma coordinates or any
! arbitrary coordinate system.
!   The MOM6 interface positions are always calculated from the bottom up by
! accumulating the layer thicknesses starting at z=-G%bathyT.  z increases
! upwards (decreasing k-index).
!   The new grid is defined by the change in position of those interfaces in z
!       dzInterface = zNew - zOld.
!   Thus, if the regridding inflates the top layer, hNew(1) > hOld(1), then the
! second interface moves downward, zNew(2) < zOld(2), and dzInterface(2) < 0.
!       hNew(k) = hOld(k) - dzInterface(k+1) + dzInterface(k)
! IMPORTANT NOTE:
!   This is the converse of the sign convention used in the remapping code!
!------------------------------------------------------------------------------

  ! Arguments
  type(remapping_CS),                         intent(in)    :: remapCS !< Remapping parameters and options
  type(regridding_CS),                        intent(in)    :: CS     !< Regridding control structure
  type(ocean_grid_type),                      intent(in)    :: G      !< Ocean grid structure
  type(verticalGrid_type),                    intent(in)    :: GV     !< Ocean vertical grid structure
  type(unit_scale_type),                      intent(in)    :: US     !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),  intent(in)    :: h      !< Current 3D grid obtained after
                                                                      !! the last time step [H ~> m or kg m-2]
  type(thermo_var_ptrs),                      intent(in)    :: tv     !< Thermodynamical variables (T, S, ...)
  real, dimension(SZI_(G),SZJ_(G),CS%nk),     intent(inout) :: h_new  !< New 3D grid consistent with target
                                                                      !! coordinate [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJ_(G),CS%nk+1),   intent(inout) :: dzInterface !< The change in position of each
                                                                      !! interface [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJ_(G)), optional, intent(in   ) :: frac_shelf_h !< Fractional ice shelf coverage [nondim]
  logical, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                                    optional, intent(out  ) :: PCM_cell !< Use PCM remapping in cells where true

  ! Local variables
                  ! parameters to depth units [H Z-1 ~> nondim or kg m-3]

end subroutine regridding_main
module subroutine regridding_preadjust_reqs(CS, do_conv_adj, do_hybgen_unmix, hybgen_CS)

  ! Arguments
  type(regridding_CS), intent(in)  :: CS          !< Regridding control structure
  logical,             intent(out) :: do_conv_adj !< Convective adjustment should be done
  logical,             intent(out) :: do_hybgen_unmix !< Hybgen unmixing should be done
  type(hybgen_regrid_CS), pointer, &
             optional, intent(out) :: hybgen_CS   !< Control structure for hybgen regridding for sharing parameters.


end subroutine regridding_preadjust_reqs
module subroutine calc_h_new_by_dz(CS, G, GV, h, dzInterface, h_new)
  type(regridding_CS),                       intent(in)    :: CS !< Regridding control structure
  type(ocean_grid_type),                     intent(in)    :: G  !< Grid structure
  type(verticalGrid_type),                   intent(in)    :: GV !< Ocean vertical grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in)    :: h  !< Old layer thicknesses [H ~> m or kg m-2]
                                                                 !! or other units
  real, dimension(SZI_(G),SZJ_(G),CS%nk+1),  intent(in)    :: dzInterface !< Change in interface positions
                                                                 !! in the same units as h [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJ_(G),CS%nk),    intent(inout) :: h_new !< New layer thicknesses in the same
                                                                 !! units as h [H ~> m or kg m-2]
  ! Local variables

end subroutine calc_h_new_by_dz
module subroutine check_grid_column( nk, h, dzInterface, msg )
  integer,               intent(in) :: nk !< Number of cells
  real, dimension(nk),   intent(in) :: h  !< Cell thicknesses [Z ~> m] or arbitrary units
  real, dimension(nk+1), intent(in) :: dzInterface !< Change in interface positions (same units as h), often [Z ~> m]
  character(len=*),      intent(in) :: msg !< Message to append to errors
  ! Local variables

end subroutine check_grid_column
module subroutine filtered_grid_motion( CS, nk, z_old, z_new, dz_g )
  type(regridding_CS),      intent(in)    :: CS !< Regridding control structure
  integer,                  intent(in)    :: nk !< Number of cells in source grid
  real, dimension(nk+1),    intent(in)    :: z_old !< Old grid position [H ~> m or kg m-2]
  real, dimension(CS%nk+1), intent(in)    :: z_new !< New grid position before filtering [H ~> m or kg m-2]
  real, dimension(CS%nk+1), intent(inout) :: dz_g  !< Change in interface positions including
                                                   !! the effects of filtering [H ~> m or kg m-2]
  ! Local variables
                  ! filtered grid movement [H ~> m or kg m-2]
                  ! that may be adjusted for numerical accuracy in a solver [H ~> m or kg m-2]
! For debugging:
!  real, dimension(nk+1) :: ddz_g_s, ddz_g_d

end subroutine filtered_grid_motion
module subroutine build_zstar_grid( CS, G, GV, h, nom_depth_H, dzInterface, frac_shelf_h, zScale)

  ! Arguments
  type(regridding_CS),                       intent(in)    :: CS !< Regridding control structure
  type(ocean_grid_type),                     intent(in)    :: G  !< Ocean grid structure
  type(verticalGrid_type),                   intent(in)    :: GV !< ocean vertical grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in)    :: h  !< Layer thicknesses [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJ_(G)),          intent(in)    :: nom_depth_H !< The bathymetric depth of this column
                                                                 !! relative to mean sea level or another locally
                                                                 !! valid reference height, converted to thickness
                                                                 !! units [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJ_(G),CS%nk+1),  intent(inout) :: dzInterface !< The change in interface depth
                                                                 !! [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZJ_(G)), optional,intent(in)    :: frac_shelf_h !< Fractional
                                                                 !! ice shelf coverage [nondim].
  real,                            optional, intent(in)    :: zScale !< Scaling factor from the target coordinate
                                                                 !! resolution in Z to desired units for zInterface,
                                                                 !! usually Z_to_H in which case it is in
                                                                 !! units of [H Z-1 ~> nondim or kg m-3]
  ! Local variables
#ifdef __DO_SAFETY_CHECKS__
#endif

end subroutine build_zstar_grid
module subroutine build_sigma_grid( CS, G, GV, h, nom_depth_H, dzInterface )
!------------------------------------------------------------------------------
! This routine builds a grid based on terrain-following coordinates.
! The module parameter coordinateResolution(:) determines the resolution in
! sigma coordinate, dSigma(:). sigma-coordinates are defined by
!   sigma = (eta-z)/(H+eta)  s.t. sigma=0 at z=eta and sigma=1 at z=-H .
!------------------------------------------------------------------------------

  ! Arguments
  type(regridding_CS),                       intent(in)    :: CS !< Regridding control structure
  type(ocean_grid_type),                     intent(in)    :: G  !< Ocean grid structure
  type(verticalGrid_type),                   intent(in)    :: GV !< ocean vertical grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in)    :: h  !< Layer thicknesses [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJ_(G)),          intent(in)    :: nom_depth_H !< The bathymetric depth of this column
                                                                 !! relative to mean sea level or another locally
                                                                 !! valid reference height, converted to thickness
                                                                 !! units [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJ_(G),CS%nk+1),  intent(inout) :: dzInterface !< The change in interface depth
                                                                 !! [H ~> m or kg m-2]

  ! Local variables
#ifdef __DO_SAFETY_CHECKS__
#endif

end subroutine build_sigma_grid
module subroutine build_rho_grid( G, GV, US, h, nom_depth_H, tv, dzInterface, remapCS, CS, frac_shelf_h )
!------------------------------------------------------------------------------
! This routine builds a new grid based on a given set of target interface
! densities (these target densities are computed by taking the mean value
! of given layer densities). The algorithm operates as follows within each
! column:
! 1. Given T & S within each layer, the layer densities are computed.
! 2. Based on these layer densities, a global density profile is reconstructed
!    (this profile is monotonically increasing and may be discontinuous)
! 3. The new grid interfaces are determined based on the target interface
!    densities.
! 4. T & S are remapped onto the new grid.
! 5. Return to step 1 until convergence or until the maximum number of
!    iterations is reached, whichever comes first.
!------------------------------------------------------------------------------

  ! Arguments
  type(regridding_CS),                        intent(in)    :: CS !< Regridding control structure
  type(ocean_grid_type),                      intent(in)    :: G  !< Ocean grid structure
  type(verticalGrid_type),                    intent(in)    :: GV !< Ocean vertical grid structure
  type(unit_scale_type),                      intent(in)    :: US !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),  intent(in)    :: h  !< Layer thicknesses [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJ_(G)),           intent(in)    :: nom_depth_H !< The bathymetric depth of this column
                                                                  !! relative to mean sea level or another locally
                                                                  !! valid reference height, converted to thickness
                                                                  !! units [H ~> m or kg m-2]
  type(thermo_var_ptrs),                      intent(in)    :: tv !< Thermodynamics structure
  real, dimension(SZI_(G),SZJ_(G),CS%nk+1),   intent(inout) :: dzInterface !< The change in interface depth
                                                                  !! [H ~> m or kg m-2]
  type(remapping_CS),                         intent(in)    :: remapCS !< The remapping control structure
  real, dimension(SZI_(G),SZJ_(G)), optional, intent(in)    :: frac_shelf_h  !< Fractional ice
                                                                  !! shelf coverage [nondim]
  ! Local variables
#ifdef __DO_SAFETY_CHECKS__
#endif

end subroutine build_rho_grid
module subroutine build_grid_HyCOM1( G, GV, US, h, nom_depth_H, tv, h_new, dzInterface, remapCS, CS, &
                              frac_shelf_h, zScale )
  type(ocean_grid_type),                     intent(in)    :: G  !< Grid structure
  type(verticalGrid_type),                   intent(in)    :: GV !< Ocean vertical grid structure
  type(unit_scale_type),                     intent(in)    :: US !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in)    :: h  !< Existing model thickness [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJ_(G)),          intent(in)    :: nom_depth_H !< The bathymetric depth of this column
                                                                 !! relative to mean sea level or another locally
                                                                 !! valid reference height, converted to thickness
                                                                 !! units [H ~> m or kg m-2]
  type(thermo_var_ptrs),                     intent(in)    :: tv !< Thermodynamics structure
  type(remapping_CS),                        intent(in)    :: remapCS !< The remapping control structure
  type(regridding_CS),                       intent(in)    :: CS !< Regridding control structure
  real, dimension(SZI_(G),SZJ_(G),CS%nk),    intent(inout) :: h_new !< New layer thicknesses [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJ_(G),CS%nk+1),  intent(inout) :: dzInterface !< Changes in interface position
                                                                 !! in thickness units [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJ_(G)), optional, intent(in)   :: frac_shelf_h !< Fractional ice shelf
                                                                 !! coverage [nondim]
  real,                            optional, intent(in)    :: zScale !< Scaling factor from the target coordinate
                                                                 !! resolution in Z to desired units for zInterface,
                                                                 !! usually Z_to_H in which case it is in
                                                                 !! units of [H Z-1 ~> nondim or kg m-3]

  ! Local variables
                          ! in thickness units [H ~> m or kg m-2]

end subroutine build_grid_HyCOM1
module subroutine build_grid_adaptive(G, GV, US, h, nom_depth_H, tv, dzInterface, remapCS, CS)
  type(ocean_grid_type),                       intent(in)    :: G    !< The ocean's grid structure
  type(verticalGrid_type),                     intent(in)    :: GV   !< The ocean's vertical grid structure
  type(unit_scale_type),                       intent(in)    :: US   !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),   intent(in)    :: h    !< Layer thicknesses [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJ_(G)),            intent(in)    :: nom_depth_H !< The bathymetric depth of this column
                                                                     !! relative to mean sea level or another locally
                                                                     !! valid reference height, converted to thickness
                                                                     !! units [H ~> m or kg m-2]
  type(thermo_var_ptrs),                       intent(in)    :: tv   !< A structure pointing to various
                                                                     !! thermodynamic variables
  type(regridding_CS),                         intent(in)    :: CS   !< Regridding control structure
  real, dimension(SZI_(G),SZJ_(G),CS%nk+1),    intent(inout) :: dzInterface !< The change in interface depth
                                                                     !! [H ~> m or kg m-2]
  type(remapping_CS),                          intent(in)    :: remapCS !< The remapping control structure

  ! local variables
  ! current interface positions and after tendency term is applied
  ! positive downward

end subroutine build_grid_adaptive
module subroutine adjust_interface_motion( CS, nk, h_old, dz_int )
  type(regridding_CS),      intent(in)    :: CS !< Regridding control structure
  integer,                  intent(in)    :: nk !< Number of layers in h_old
  real, dimension(nk),      intent(in)    :: h_old  !< Layer thicknesses on the old grid [H ~> m or kg m-2]
  real, dimension(CS%nk+1), intent(inout) :: dz_int !< Interface movements, adjusted to keep the thicknesses
                                                    !! thicker than their minimum value [H ~> m or kg m-2]
  ! Local variables
                  ! that can not be explained by roundoff errors [H ~> m or kg m-2]

end subroutine adjust_interface_motion
module subroutine inflate_vanished_layers_old( CS, G, GV, h )
!------------------------------------------------------------------------------
! This routine is called when initializing the regridding options. The
! objective is to make sure all layers are at least as thick as the minimum
! thickness allowed for regridding purposes (this parameter is set in the
! MOM_input file or defaulted to 1.0e-3). When layers are too thin, they
! are inflated up to the minimum thickness.
!------------------------------------------------------------------------------

  ! Arguments
  type(regridding_CS),                       intent(in)    :: CS   !< Regridding control structure
  type(ocean_grid_type),                     intent(in)    :: G    !< The ocean's grid structure
  type(verticalGrid_type),                   intent(in)    :: GV   !< The ocean's vertical grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(inout) :: h    !< Layer thicknesses [H ~> m or kg m-2]

  ! Local variables

end subroutine inflate_vanished_layers_old
module subroutine convective_adjustment(G, GV, h, tv)
  type(ocean_grid_type),   intent(in)    :: G    !< The ocean's grid structure
  type(verticalGrid_type), intent(in)    :: GV   !< The ocean's vertical grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(inout) :: h    !< Layer thicknesses [H ~> m or kg m-2]
  type(thermo_var_ptrs),   intent(inout) :: tv   !< A structure pointing to various thermodynamic variables
!------------------------------------------------------------------------------
! Check each water column to see whether it is stratified. If not, sort the
! layers by successive swappings of water masses (bubble sort algorithm)
!------------------------------------------------------------------------------

  ! Local variables

  !### Doing convective adjustment based on potential densities with zero pressure seems
  !    questionable, although it does avoid ambiguous sorting. -RWH
end subroutine convective_adjustment
module function uniformResolution(nk,coordMode,maxDepth,rhoLight,rhoHeavy)
!------------------------------------------------------------------------------
! Calculate a vector of uniform resolution in the units of the coordinate
!------------------------------------------------------------------------------
  ! Arguments
  integer,          intent(in) :: nk !< Number of cells in source grid
  character(len=*), intent(in) :: coordMode !< A string indicating the coordinate mode.
                                            !! See the documentation for regrid_consts
                                            !! for the recognized values.
  real,             intent(in) :: maxDepth  !< The range of the grid values in some modes, in coordinate
                                            !! dependent units that might be [m] or [kg m-3] or [nondim]
                                            !! or something else.
  real,             intent(in) :: rhoLight  !< The minimum value of the grid in RHO mode [kg m-3]
  real,             intent(in) :: rhoHeavy  !< The maximum value of the grid in RHO mode [kg m-3]

  real                         :: uniformResolution(nk) !< The returned uniform resolution grid, in
                                            !! coordinate dependent units that might be [m] or
                                            !! [kg m-3] or [nondim] or something else.

  ! Local variables

end function uniformResolution
module subroutine initCoord(CS, G, GV, US, coord_mode, param_file)
  type(regridding_CS),     intent(inout) :: CS !< Regridding control structure
  type(ocean_grid_type),   intent(in)    :: G  !< Ocean grid structure
  type(verticalGrid_type), intent(in)    :: GV !< Ocean vertical grid structure
  type(unit_scale_type),   intent(in)    :: US !< A dimensional unit scaling type
  character(len=*),        intent(in)    :: coord_mode !< A string indicating the coordinate mode.
                                               !! See the documentation for regrid_consts
                                               !! for the recognized values.
  type(param_file_type),   intent(in)    :: param_file !< Parameter file

end subroutine initCoord
module subroutine setCoordinateResolution( dz, CS, scale )
  real, dimension(:),  intent(in)    :: dz !< A vector of vertical grid spacings, in arbitrary coordinate
                                           !! dependent units, such as [m] for a z-coordinate or [kg m-3]
                                           !! for a density coordinate.
  type(regridding_CS), intent(inout) :: CS !< Regridding control structure
  real,      optional, intent(in)    :: scale !< A scaling factor converting dz to the internal represetation
                                           !! of coordRes, in various units that depend on the coordinate,
                                           !! such as [Z m-1 ~> 1] for a z-coordinate or [R m3 kg-1 ~> 1] for
                                           !! a density coordinate.

end subroutine setCoordinateResolution
module subroutine setCoordinateResolution_3d( dz_3d, CS, scale )
  real, dimension(:,:,:),  intent(in)    :: dz_3d !< A vector of vertical grid spacings, in arbitrary coordinate
                                           !! dependent units, such as [m] for a z-coordinate or [kg m-3]
                                           !! for a density coordinate.
  type(regridding_CS), intent(inout) :: CS !< Regridding control structure
  real,      optional, intent(in)    :: scale !< A scaling factor converting dz to coordRes [Z m-1 ~> 1]

end subroutine setCoordinateResolution_3d
module subroutine set_target_densities_from_GV( GV, US, CS )
  type(verticalGrid_type), intent(in)    :: GV !< Ocean vertical grid structure
  type(unit_scale_type),   intent(in)    :: US !< A dimensional unit scaling type
  type(regridding_CS),     intent(inout) :: CS !< Regridding control structure
  ! Local variables

end subroutine set_target_densities_from_GV
module subroutine set_target_densities_3d( CS, G, scale, rho_int_3d )
  type(regridding_CS),  intent(inout) :: CS    !< Regridding control structure
  type(ocean_grid_type),intent(in)    :: G     !< Ocean grid structure
  real,                 intent(in)    :: scale !< A scaling factor converting densities [R m3 kg-1 ~> 1]
  real, dimension(SZI_(G),SZJ_(G),CS%nk+1), intent(in) :: rho_int_3d !< Interface densities [kg m-3]

end subroutine set_target_densities_3d
module subroutine set_target_densities( CS, rho_int )
  type(regridding_CS),      intent(inout) :: CS !< Regridding control structure
  real, dimension(CS%nk+1), intent(in)    :: rho_int !< Interface densities [R ~> kg m-3]

end subroutine set_target_densities
module subroutine set_regrid_max_depths( CS, max_depths, units_to_H )
  type(regridding_CS),      intent(inout) :: CS !< Regridding control structure
  real, dimension(CS%nk+1), intent(in)    :: max_depths !< Maximum interface depths, in arbitrary units, often [m]
  real, optional,           intent(in)    :: units_to_H !< A conversion factor for max_depths into H units,
                                                        !! often in [H m-1 ~> 1 or kg m-3]
  ! Local variables
                   ! if units_to_H is present, or [nondim] if it is absent.

end subroutine set_regrid_max_depths
module subroutine set_regrid_max_thickness( CS, max_h, units_to_H )
  type(regridding_CS),      intent(inout) :: CS !< Regridding control structure
  real, dimension(CS%nk+1), intent(in)    :: max_h !< Maximum layer thicknesses, in arbitrary units, often [m]
  real, optional,           intent(in)    :: units_to_H !< A conversion factor for max_h into H units,
                                                        !! often [H m-1 ~> 1 or kg m-3]
  ! Local variables
                   ! if units_to_H is present, or [nondim] if it is absent.

end subroutine set_regrid_max_thickness
module subroutine write_regrid_file( CS, GV, filepath )
  type(regridding_CS),     intent(in) :: CS        !< Regridding control structure
  type(verticalGrid_type), intent(in) :: GV        !< ocean vertical grid structure
  character(len=*),        intent(in) :: filepath  !< The full path to the file to write

                                                 ! in axes in files, in coordinate-dependent units that can
                                                 ! be obtained from getCoordinateUnits [various]

end subroutine write_regrid_file
module function set_h_neglect(GV, remap_answer_date, h_neglect_edge) result(h_neglect)
  type(verticalGrid_type), intent(in)  :: GV   !< Ocean vertical grid structure
  integer,                 intent(in)  :: remap_answer_date !< The vintage of the expressions to use
                                               !! for remapping.  Values below 20190101 recover the
                                               !! remapping answers from 2018. Higher values use more
                                               !! robust forms of the same remapping algorithms.
  real,                    intent(out) :: h_neglect_edge !< A negligibly small thickness used in
                                               !! remapping edge value calculations [H ~> m or kg m-2]
  real                                 :: h_neglect !< A negligibly small thickness used in
                                               !! remapping cell reconstructions [H ~> m or kg m-2]

end function set_h_neglect
module function set_dz_neglect(GV, US, remap_answer_date, dz_neglect_edge) result(dz_neglect)
  type(verticalGrid_type), intent(in)  :: GV   !< Ocean vertical grid structure
  type(unit_scale_type),   intent(in)  :: US   !< A dimensional unit scaling type
  integer,                 intent(in)  :: remap_answer_date !< The vintage of the expressions to use
                                               !! for remapping.  Values below 20190101 recover the
                                               !! remapping answers from 2018. Higher values use more
                                               !! robust forms of the same remapping algorithms.
  real,                    intent(out) :: dz_neglect_edge !< A negligibly small vertical layer extent
                                               !! used in remapping edge value calculations [Z ~> m]
  real                                 :: dz_neglect !< A negligibly small vertical layer extent
                                               !! used in remapping cell reconstructions [Z ~> m]

end function set_dz_neglect
module function getCoordinateResolution( CS, undo_scaling )
  type(regridding_CS), intent(in) :: CS !< Regridding control structure
  logical,   optional, intent(in) :: undo_scaling !< If present and true, undo any internal
                                        !! rescaling of the resolution data.
  real, dimension(CS%nk)          :: getCoordinateResolution !< The resolution or delta of the target coordinate,
                                                             !! in units that depend on the coordinate [various]

end function getCoordinateResolution
module function getCoordinateInterfaces( CS, undo_scaling )
  type(regridding_CS), intent(in) :: CS                      !< Regridding control structure
  logical,   optional, intent(in) :: undo_scaling            !< If present and true, undo any internal
                                                             !! rescaling of the resolution data.
  real, dimension(CS%nk+1)        :: getCoordinateInterfaces !< Interface positions in target coordinate,
                                                             !! in units that depend on the coordinate [various]

end function getCoordinateInterfaces
module function getCoordinateUnits( CS )
  type(regridding_CS), intent(in) :: CS !< Regridding control structure
  character(len=20)               :: getCoordinateUnits

end function getCoordinateUnits
module function getCoordinateShortName( CS )
  type(regridding_CS), intent(in) :: CS !< Regridding control structure
  character(len=20)               :: getCoordinateShortName

end function getCoordinateShortName
module subroutine set_regrid_params( CS, boundary_extrapolation, min_thickness, old_grid_weight, &
             use_depth_based_time_filter, depth_of_time_filter_shallow, depth_of_time_filter_deep, &
             interp_scheme, use_adjust_interface_motion, compress_fraction, ref_pressure, &
             integrate_downward_for_e, remap_answers_2018, remap_answer_date, regrid_answer_date, &
             adaptTimeRatio, adaptZoom, adaptZoomCoeff, adaptBuoyCoeff, &
             adaptAlpha, adaptDoMin, adaptDrho0)
  type(regridding_CS), intent(inout) :: CS !< Regridding control structure
  logical, optional, intent(in) :: boundary_extrapolation !< Extrapolate in boundary cells
  real,    optional, intent(in) :: min_thickness    !< Minimum thickness allowed when building the
                                                    !! new grid [H ~> m or kg m-2]
  real,    optional, intent(in) :: old_grid_weight  !< Weight given to old coordinate when time-filtering grid [nondim]
  logical, optional, intent(in) :: use_depth_based_time_filter !< Allow depth-based time filtering
  real,    optional, intent(in) :: depth_of_time_filter_shallow !< Depth to start cubic [H ~> m or kg m-2]
  real,    optional, intent(in) :: depth_of_time_filter_deep !< Depth to end cubic [H ~> m or kg m-2]
  character(len=*), optional, intent(in) :: interp_scheme !< Interpolation method for state-dependent coordinates
  logical, optional, intent(in) :: use_adjust_interface_motion !< Call adjust_interface_motion()
  real,    optional, intent(in) :: compress_fraction !< Fraction of compressibility to add to potential density [nondim]
  real,    optional, intent(in) :: ref_pressure     !< The reference pressure for density-dependent
                                                    !! coordinates [R L2 T-2 ~> Pa]
  logical, optional, intent(in) :: integrate_downward_for_e !< If true, integrate for interface positions downward
                                                    !! from the top.
  logical, optional, intent(in) :: remap_answers_2018 !< If true, use the order of arithmetic and expressions
                                                    !! that recover the remapping answers from 2018.  Otherwise
                                                    !! use more robust but mathematically equivalent expressions.
  integer, optional, intent(in) :: remap_answer_date !< The vintage of the expressions to use for remapping
  integer, optional, intent(in) :: regrid_answer_date !< The vintage of the expressions to use for regridding
  real,    optional, intent(in) :: adaptTimeRatio   !< Ratio of the ALE timestep to the grid timescale [nondim].
  real,    optional, intent(in) :: adaptZoom        !< Depth of near-surface zooming region [H ~> m or kg m-2].
  real,    optional, intent(in) :: adaptZoomCoeff   !< Coefficient of near-surface zooming diffusivity [nondim].
  real,    optional, intent(in) :: adaptBuoyCoeff   !< Coefficient of buoyancy diffusivity [nondim].
  real,    optional, intent(in) :: adaptAlpha       !< Scaling factor on optimization tendency [nondim].
  logical, optional, intent(in) :: adaptDoMin       !< If true, make a HyCOM-like mixed layer by
                                                    !! preventing interfaces from being shallower than
                                                    !! the depths specified by the regridding coordinate.
  real,    optional, intent(in) :: adaptDrho0       !< Reference density difference for stratification-dependent
                                                    !! diffusion. [R ~> kg m-3]

end subroutine set_regrid_params
integer module function get_regrid_size(CS)
  type(regridding_CS), intent(inout) :: CS !< Regridding control structure

end function get_regrid_size
module function get_zlike_CS(CS)
  type(regridding_CS), intent(in) :: CS !< Regridding control structure
  type(zlike_CS) :: get_zlike_CS

end function get_zlike_CS
module function get_sigma_CS(CS)
  type(regridding_CS), intent(in) :: CS !< Regridding control structure
  type(sigma_CS) :: get_sigma_CS

end function get_sigma_CS
module function get_rho_CS(CS)
  type(regridding_CS), intent(in) :: CS !< Regridding control structure
  type(rho_CS) :: get_rho_CS

end function get_rho_CS
module function getStaticThickness( CS, SSH, depth )
  type(regridding_CS), intent(in) :: CS !< Regridding control structure
  real,                intent(in) :: SSH   !< The sea surface height, in the same units as depth, often [Z ~> m]
  real,                intent(in) :: depth !< The maximum depth of the grid, often [Z ~> m]
  real, dimension(CS%nk)          :: getStaticThickness !< The returned thicknesses in the units of
                                           !! depth, often [Z ~> m]
  ! Local

end function getStaticThickness
module subroutine dz_function1( string, dz )
  character(len=*),   intent(in)    :: string !< String with list of parameters in form
                                              !! dz_min, H_total, power, precision
  real, dimension(:), intent(inout) :: dz     !< Profile of nominal thicknesses [m] or other units
  ! Local variables

end subroutine dz_function1
module function create_coord_param(param_prefix, param_name, param_suffix) result(coord_param)
  character(len=*) :: param_name   !< The base name of the parameter (e.g. the one used for the main coordinate)
  character(len=*) :: param_prefix !< String to prefix to parameter names.
  character(len=*) :: param_suffix !< String to append to parameter names.
  character(len=MAX_PARAM_LENGTH) :: coord_param  !< Parameter name prepended by param_prefix
                                                  !! and appended with param_suffix

end function create_coord_param
integer module function rho_function1( string, rho_target )
  character(len=*),   intent(in)    :: string !< String with list of parameters in form
                                              !! dz_min, H_total, power, precision
  real, dimension(:), allocatable, intent(inout) :: rho_target !< Profile of interface densities [kg m-3]
  ! Local variables
                      ! in subsequent layers [kg m-3]

end function rho_function1
  end interface

end module MOM_regridding
