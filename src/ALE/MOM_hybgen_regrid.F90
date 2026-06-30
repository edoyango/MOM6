! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> This module contains the hybgen regridding routines from HYCOM, with minor
!! modifications to follow the MOM6 coding conventions
module MOM_hybgen_regrid

use MOM_EOS,              only : EOS_type, calculate_density
use MOM_error_handler,    only : MOM_mesg, MOM_error, FATAL, WARNING, assert
use MOM_file_parser,      only : get_param, param_file_type, log_param
use MOM_io,               only : create_MOM_file, file_exists
use MOM_io,               only : MOM_infra_file, MOM_field
use MOM_io,               only : MOM_read_data, MOM_write_field, vardesc, var_desc, SINGLE_FILE
use MOM_string_functions, only : slasher
use MOM_unit_scaling,     only : unit_scale_type
use MOM_variables,        only : ocean_grid_type, thermo_var_ptrs
use MOM_verticalGrid,     only : verticalGrid_type

implicit none ; private

#include <MOM_memory.h>

!> Control structure containing required parameters for the hybgen coordinate generator
type, public :: hybgen_regrid_CS ; private

  real :: min_thickness !< Minimum thickness allowed for layers [H ~> m or kg m-2]

  integer :: nk         !< Number of layers on the target grid

  !> Reference pressure for density calculations [R L2 T-2 ~> Pa]
  real :: ref_pressure

  !> Hybgen uses PCM if layer is within hybiso of target density [R ~> kg m-3]
  real :: hybiso
  !> Number of sigma levels used by HYBGEN
  integer :: nsigma

  real :: dp00i    !< Deep isopycnal spacing minimum thickness [H ~> m or kg m-2]
  real :: qhybrlx  !< Fractional relaxation within a regridding step [nondim]

  real, allocatable, dimension(:) ::  &
    dp0k, & !< minimum deep    z-layer separation [H ~> m or kg m-2]
    ds0k    !< minimum shallow z-layer separation [H ~> m or kg m-2]

  real :: coord_scale = 1.0     !< A scaling factor to restores the depth coordinates to
                                !! values in m [m H-1 ~> 1 or m3 kg-1]
  real :: Rho_coord_scale = 1.0 !< A scaling factor to restores the denesity coordinates to
                                !! values in kg m-3 [kg m-3 R-1 ~> 1]

  real :: dpns  !< depth to start terrain following [H ~> m or kg m-2]
  real :: dsns  !< depth to stop terrain following [H ~> m or kg m-2]
  real :: min_dilate !< The minimum amount of dilation that is permitted when converting target
                     !! coordinates from z to z* [nondim].  This limit applies when wetting occurs.
  real :: max_dilate !< The maximum amount of dilation that is permitted when converting target
                     !! coordinates from z to z* [nondim].  This limit applies when drying occurs.

  real :: thkbot !< Thickness of a bottom boundary layer, within which hybgen does
                 !! something different. [H ~> m or kg m-2]

  !> Shallowest depth for isopycnal layers [H ~> m or kg m-2]
  real :: topiso_const
  ! real, dimension(:,:), allocatable :: topiso

  !> Nominal density of interfaces [R ~> kg m-3]
  real, allocatable, dimension(:) :: target_density

  real :: dp_far_from_sfc  !< A distance that determines when an interface is suffiently far from
                     !! the surface that certain adjustments can be made in the Hybgen regridding
                     !! code [H ~> m or kg m-2].  In Hycom, this is set to tenm (nominally 10 m).
  real :: dp_far_from_bot  !< A distance that determines when an interface is suffiently far from
                     !! the bottom that certain adjustments can be made in the Hybgen regridding
                     !! code [H ~> m or kg m-2].  In Hycom, this is set to onem (nominally 1 m).
  real :: h_thin     !< A layer thickness below which a layer is considered to be too thin for
                     !! certain adjustments to be made in the Hybgen regridding code [H ~> m or kg m-2].
                     !! In Hycom, this is set to onemm (nominally 0.001 m).

  real :: rho_eps    !< A small nonzero density that is used to prevent division by zero
                     !! in several expressions in the Hybgen regridding code [R ~> kg m-3].

  real :: onem       !< Nominally one m in thickness units [H ~> m or kg m-2], used only in
                     !! certain debugging tests.

end type hybgen_regrid_CS


public hybgen_regrid, init_hybgen_regrid, end_hybgen_regrid
public hybgen_column_init, get_hybgen_regrid_params, write_Hybgen_coord_file


  interface
module subroutine init_hybgen_regrid(CS, GV, US, param_file)
  type(hybgen_regrid_CS),  pointer    :: CS  !< Unassociated pointer to hold the control structure
  type(verticalGrid_type), intent(in) :: GV  !< Ocean vertical grid structure
  type(unit_scale_type),   intent(in) :: US  !< A dimensional unit scaling type
  type(param_file_type),   intent(in) :: param_file !< Parameter file


end subroutine init_hybgen_regrid
module subroutine write_Hybgen_coord_file(GV, CS, filepath)
  type(verticalGrid_type), intent(in)  :: GV        !< The ocean's vertical grid structure
  type(hybgen_regrid_CS),  intent(in)  :: CS        !< Control structure for this module
  character(len=*),        intent(in)  :: filepath  !< The full path to the file to write
  ! Local variables

end subroutine write_Hybgen_coord_file
module subroutine end_hybgen_regrid(CS)
  type(hybgen_regrid_CS), pointer :: CS !< Coordinate control structure

  ! nothing to do
end subroutine end_hybgen_regrid
module subroutine get_hybgen_regrid_params(CS, nk, ref_pressure, hybiso, nsigma, dp00i, qhybrlx, &
                                    dp0k, ds0k, dpns, dsns, min_dilate, max_dilate, &
                                    thkbot, topiso_const, target_density)
  type(hybgen_regrid_CS),  pointer    :: CS !< Coordinate regridding control structure
  integer, optional, intent(out) :: nk  !< Number of layers on the target grid
  real,    optional, intent(out) :: ref_pressure !< Reference pressure for density calculations [R L2 T-2 ~> Pa]
  real,    optional, intent(out) :: hybiso  !< Hybgen uses PCM if layer is within hybiso of target density [R ~> kg m-3]
  integer, optional, intent(out) :: nsigma  !< Number of sigma levels used by HYBGEN
  real,    optional, intent(out) :: dp00i   !< Deep isopycnal spacing minimum thickness [H ~> m or kg m-2]
  real,    optional, intent(out) :: qhybrlx !< Fractional relaxation amount per timestep, 0 < qyhbrlx <= 1 [nondim]
  real,    optional, intent(out) :: dp0k(:) !< minimum deep    z-layer separation [H ~> m or kg m-2]
  real,    optional, intent(out) :: ds0k(:) !< minimum shallow z-layer separation [H ~> m or kg m-2]
  real,    optional, intent(out) :: dpns    !< depth to start terrain following [H ~> m or kg m-2]
  real,    optional, intent(out) :: dsns    !< depth to stop terrain following [H ~> m or kg m-2]
  real,    optional, intent(out) :: min_dilate !< The minimum amount of dilation that is permitted when
                                            !! converting target coordinates from z to z* [nondim].
                                            !! This limit applies when wetting occurs.
  real,    optional, intent(out) :: max_dilate !< The maximum amount of dilation that is permitted when
                                            !! converting target coordinates from z to z* [nondim].
                                            !! This limit applies when drying occurs.
  real,    optional, intent(out) :: thkbot  !< Thickness of a bottom boundary layer, within which
                                            !! hybgen does something different. [H ~> m or kg m-2]
  real,    optional, intent(out) :: topiso_const !< Shallowest depth for isopycnal layers [H ~> m or kg m-2]
  ! real, dimension(:,:), allocatable :: topiso
  real,    optional, intent(out) :: target_density(:) !< Nominal density of interfaces [R ~> kg m-3]

end subroutine get_hybgen_regrid_params
module subroutine hybgen_regrid(G, GV, US, dp, nom_depth_H, tv, CS, dzInterface, PCM_cell)
  type(ocean_grid_type),   intent(in)    :: G   !< Ocean grid structure
  type(verticalGrid_type), intent(in)    :: GV  !< Ocean vertical grid structure
  type(unit_scale_type),   intent(in)    :: US  !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in)    :: dp  !< Source grid layer thicknesses [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJ_(G)), &
                           intent(in)    :: nom_depth_H !< The bathymetric depth of this column
                                                !! relative to mean sea level or another locally
                                                !! valid reference height, converted to thickness
                                                !! units [H ~> m or kg m-2]
  type(thermo_var_ptrs),   intent(in)    :: tv  !< Thermodynamics structure
  type(hybgen_regrid_CS),  intent(in)    :: CS  !< hybgen control structure
  real, dimension(SZI_(G),SZJ_(G),CS%nk+1), &
                           intent(inout) :: dzInterface !< The change in height of each interface,
                                                !! using a sign convention opposite to the change
                                                !! in pressure [H ~> m or kg m-2]
  logical, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                 optional, intent(inout) :: PCM_cell !< If true, PCM remapping should be used in a cell.
                                                !! This is effectively intent out, but values in wide
                                                !! halo regions and land points are reused.

  ! --- -------------------------------------
  ! --- hybrid grid generator from HYCOM
  ! --- -------------------------------------

  ! These notes on the parameters for the hybrid grid generator are inhereted from the
  ! Hycom source code for these algorithms.
  !
  ! From blkdat.input (units may have changed from m to pressure):
  !
  ! --- 'nsigma' = number of sigma  levels
  ! --- 'dp0k  ' = layer k deep    z-level spacing minimum thickness (m)
  ! ---              k=1,nk
  ! --- 'ds0k  ' = layer k shallow z-level spacing minimum thickness (m)
  ! ---              k=1,nsigma
  ! --- 'dp00i'  = deep isopycnal spacing minimum thickness (m)
  ! --- 'isotop' = shallowest depth for isopycnal layers     (m)
  !                now in topiso(:,:)
  ! --- 'sigma ' = isopycnal layer target densities (sigma units)
  ! ---            now in theta(:,:,1:nk)
  !
  ! --- the above specifies a vertical coord. that is isopycnal or:
  ! ---  near surface z in    deep water, based on dp0k
  ! ---  near surface z in shallow water, based on ds0k and nsigma
  ! ---   terrain-following between them, based on ds0k and nsigma
  !
  ! --- terrain following starts at depth dpns=sum(dp0k(k),k=1,nsigma) and
  ! --- ends at depth dsns=sum(ds0k(k),k=1,nsigma), and the depth of the
  ! --- k-th layer interface varies linearly with total depth between
  ! --- these two reference depths, i.e. a z-sigma-z fixed coordinate.
  !
  ! --- near the surface (i.e. shallower than isotop), layers are always
  ! --- fixed depth (z or sigma).
  ! --  layer 1 is always fixed, so isotop=0.0 is not realizable.
  ! --- near surface layers can also be forced to be fixed depth
  ! --- by setting target densities (sigma(k)) very small.
  !
  ! --- away from the surface, the minimum layer thickness is dp00i.
  !
  ! --- for fixed depth targets to be:
  ! ---  z-only set nsigma=0,
  ! ---  sigma-z (shallow-deep) use a very small ds0k(:),
  ! ---  sigma-only set nsigma=nk, dp0k large, and ds0k small.

  ! These arrays work with the input column

  ! These arrays are on the target grid.


                            ! each target layer, in the unusual case where the input grid is
                            ! larger than the new grid.  This situation only occurs during certain
                            ! types of initialization or when generating output diagnostics.

end subroutine hybgen_regrid
module subroutine hybgen_column_init(nk, nsigma, dp0k, ds0k, dp00i, topiso_i_j, &
                          qhybrlx, dpns, dsns, h_tot, dilate, h_col, &
                          fixlay, qhrlx, dp0ij, dp0cum)
  integer, intent(in)    :: nk           !< The number of layers in the new grid
  integer, intent(in)    :: nsigma       !< The number of sigma  levels
  real,    intent(in)    :: dp0k(nk)     !< Layer deep z-level spacing minimum thicknesses [H ~> m or kg m-2]
  real,    intent(in)    :: ds0k(nsigma) !< Layer shallow z-level spacing minimum thicknesses [H ~> m or kg m-2]
  real,    intent(in)    :: dp00i        !< Deep isopycnal spacing minimum thickness [H ~> m or kg m-2]
  real,    intent(in)    :: topiso_i_j   !< Shallowest depth for isopycnal layers [H ~> m or kg m-2]
  real,    intent(in)    :: qhybrlx      !< Fractional relaxation amount per timestep, 0 < qyhbrlx <= 1 [nondim]
  real,    intent(in)    :: h_tot        !< The sum of the initial layer thicknesses [H ~> m or kg m-2]
  real,    intent(in)    :: dilate       !< A factor by which to dilate the target positions
                                         !! from z to z* [nondim]
  real,    intent(in)    :: h_col(nk)    !< Initial layer thicknesses [H ~> m or kg m-2]
  real,    intent(in)    :: dpns         !< Vertical sum of dp0k [H ~> m or kg m-2]
  real,    intent(in)    :: dsns         !< Vertical sum of ds0k [H ~> m or kg m-2]
  integer, intent(out)   :: fixlay       !< Deepest fixed coordinate layer
  real,    intent(out)   :: qhrlx(nk+1)  !< Fractional relaxation within a timestep (between 0 and 1) [nondim]
  real,    intent(out)   :: dp0ij(nk)    !< minimum layer thickness [H ~> m or kg m-2]
  real,    intent(out)   :: dp0cum(nk+1) !< minimum interface depth [H ~> m or kg m-2]

  ! --- --------------------------------------------------------------
  ! --- hybrid grid generator, single column - initialization.
  ! --- --------------------------------------------------------------

  ! Local variables

  ! --- dpns = sum(dp0k(k),k=1,nsigma)
  ! --- dsns = sum(ds0k(k),k=1,nsigma)
  ! --- terrain following starts (on the deep side) at depth dpns and ends (on the
  ! --- shallow side) at depth dsns and the depth of the k-th layer interface varies
  ! --- linearly with total depth between these two reference depths.
end subroutine hybgen_column_init
real module function cushn(delp, dp0)
  real, intent(in) :: delp  ! A thickness change [H ~> m or kg m-2]
  real, intent(in) :: dp0   ! A non-negative reference thickness [H ~> m or kg m-2]

  ! These are the nondimensional parameters that define the cushion function.
! real, parameter :: qqmn=-2.0, qqmx=4.0  ! traditional range for cushn [nondim]
! real, parameter :: qqmn=-4.0, qqmx=6.0  ! somewhat wider range for cushn [nondim]
  ! These are derivative nondimensional parameters.
  ! real, parameter :: cusha = qqmn**2 * (qqmx-1.0) / (qqmx-qqmn)**2
  ! real, parameter :: I_qqmn = 1.0 / qqmn

  ! --- if delp >= qqmx*dp0 >>  dp0, cushn returns delp.
  ! --- if delp <= qqmn*dp0 << -dp0, cushn returns dp0.

  ! This is the original version from Hycom.
  ! qq = max(qqmn, min(qqmx, delp/dp0))
  ! cushn = dp0 * (1.0 + cusha * (1.0-I_qqmn*qq)**2) * max(1.0, delp/(dp0*qqmx))

  ! This is mathematically equivalent, has one fewer divide, and works as intended even if dp0 = 0.
end function cushn
module subroutine hybgen_column_regrid(CS, nk, thkbot, Rcv_tgt, &
                                fixlay, qhrlx, dp0ij, dp0cum, Rcv, h_in, dp_int)
  type(hybgen_regrid_CS), intent(in)    :: CS  !< hybgen regridding control structure
  integer, intent(in)    :: nk            !< number of layers
  real,    intent(in)    :: thkbot        !< thickness of bottom boundary layer [H ~> m or kg m-2]
  real,    intent(in)    :: Rcv_tgt(nk)   !< Target potential density [R ~> kg m-3]
  integer, intent(in)    :: fixlay        !< deepest fixed coordinate layer
  real,    intent(in)    :: qhrlx( nk+1)  !< relaxation coefficient per timestep [nondim]
  real,    intent(in)    :: dp0ij( nk)    !< minimum layer thickness [H ~> m or kg m-2]
  real,    intent(in)    :: dp0cum(nk+1)  !< minimum interface depth [H ~> m or kg m-2]
  real,    intent(in)    :: Rcv(nk)       !< Coordinate potential density [R ~> kg m-3]
  real,    intent(in)    :: h_in(nk)      !< Layer thicknesses [H ~> m or kg m-2]
  real,    intent(out)   :: dp_int(nk+1)  !< The change in interface positions [H ~> m or kg m-2]

  ! --- ------------------------------------------------------
  ! --- hybrid grid generator, single column - regrid.
  ! --- ------------------------------------------------------

  ! Local variables
                 ! between layers k and k-1 [H ~> m or kg m-2]

  ! This line needs to be consistent with the parameters set in cushn().
! real, parameter :: qqmn=-2.0, qqmx=4.0  ! traditional range for cushn [nondim]
! real, parameter :: qqmn=-4.0, qqmx=6.0  ! somewhat wider range for cushn [nondim]

end subroutine hybgen_column_regrid
  end interface

end module MOM_hybgen_regrid
