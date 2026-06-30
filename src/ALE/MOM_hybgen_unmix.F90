! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> This module contains the hybgen unmixing routines from HYCOM, with
!! modifications to follow the MOM6 coding conventions and several bugs fixed
module MOM_hybgen_unmix

use MOM_EOS,             only : EOS_type, calculate_density, calculate_density_derivs
use MOM_error_handler,   only : MOM_mesg, MOM_error, FATAL, WARNING
use MOM_file_parser,     only : get_param, param_file_type, log_param
use MOM_hybgen_regrid,   only : hybgen_column_init
use MOM_hybgen_regrid,   only : hybgen_regrid_CS, get_hybgen_regrid_params
use MOM_interface_heights, only : calc_derived_thermo
use MOM_tracer_registry, only : tracer_registry_type, tracer_type, MOM_tracer_chkinv
use MOM_unit_scaling,    only : unit_scale_type
use MOM_variables,       only : ocean_grid_type, thermo_var_ptrs
use MOM_verticalGrid,    only : verticalGrid_type

implicit none ; private

#include <MOM_memory.h>

!> Control structure containing required parameters for the hybgen coordinate generator
type, public :: hybgen_unmix_CS ; private

  integer :: nsigma  !< Number of sigma levels used by HYBGEN
  real :: hybiso     !< Hybgen uses PCM if layer is within hybiso of target density [R ~> kg m-3]

  real :: dp00i   !< Deep isopycnal spacing minimum thickness [H ~> m or kg m-2]
  real :: qhybrlx !< Hybgen relaxation amount per thermodynamic time steps [nondim]

  real, allocatable, dimension(:) ::  &
    dp0k, &     !< minimum deep    z-layer separation [H ~> m or kg m-2]
    ds0k        !< minimum shallow z-layer separation [H ~> m or kg m-2]

  real :: dpns  !< depth to start terrain following [H ~> m or kg m-2]
  real :: dsns  !< depth to stop terrain following [H ~> m or kg m-2]
  real :: min_dilate !< The minimum amount of dilation that is permitted when converting target
                     !! coordinates from z to z* [nondim].  This limit applies when wetting occurs.
  real :: max_dilate !< The maximum amount of dilation that is permitted when converting target
                     !! coordinates from z to z* [nondim].  This limit applies when drying occurs.

  real :: topiso_const !< Shallowest depth for isopycnal layers [H ~> m or kg m-2]
  ! real, dimension(:,:), allocatable :: topiso

  real :: ref_pressure !< Reference pressure for density calculations [R L2 T-2 ~> Pa]
  real, allocatable, dimension(:) :: target_density !< Nominal density of interfaces [R ~> kg m-3]

end type hybgen_unmix_CS

public hybgen_unmix, init_hybgen_unmix, end_hybgen_unmix
public set_hybgen_unmix_params


  interface
module subroutine init_hybgen_unmix(CS, GV, US, param_file, hybgen_regridCS)
  type(hybgen_unmix_CS),   pointer    :: CS  !< Unassociated pointer to hold the control structure
  type(verticalGrid_type), intent(in) :: GV  !< Ocean vertical grid structure
  type(unit_scale_type),   intent(in) :: US  !< A dimensional unit scaling type
  type(param_file_type),   intent(in) :: param_file !< Parameter file
  type(hybgen_regrid_CS),  pointer    :: hybgen_regridCS !< Control structure for hybgen
                                             !! regridding for sharing parameters.

end subroutine init_hybgen_unmix
module subroutine end_hybgen_unmix(CS)
  type(hybgen_unmix_CS), pointer :: CS !< Coordinate control structure

  ! nothing to do
end subroutine end_hybgen_unmix
module subroutine set_hybgen_unmix_params(CS, min_thickness)
  type(hybgen_unmix_CS),  pointer    :: CS !< Coordinate unmixing control structure
  real,    optional, intent(in) :: min_thickness !< Minimum allowed thickness [H ~> m or kg m-2]

end subroutine set_hybgen_unmix_params
module subroutine hybgen_unmix(G, GV, US, CS, tv, Reg, ntr, h)
  type(ocean_grid_type),   intent(in)    :: G   !< Ocean grid structure
  type(verticalGrid_type), intent(in)    :: GV  !< Ocean vertical grid structure
  type(unit_scale_type),   intent(in)    :: US  !< A dimensional unit scaling type
  type(hybgen_unmix_CS),   intent(in)    :: CS  !< hybgen control structure
  type(thermo_var_ptrs),   intent(inout) :: tv  !< Thermodynamics structure
  type(tracer_registry_type), pointer    :: Reg !< Tracer registry structure
  integer,                 intent(in)    :: ntr !< The number of tracers in the registry, or
                                                !! 0 if the registry is not in use.
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(inout) :: h   !< Layer thicknesses [H ~> m or kg m-2]

! --- --------------------------------------------
! --- hybrid grid generator, single j-row (part A).
! --- --------------------------------------------


                            ! vanished layers [H ~> m or kg m-2]



end subroutine hybgen_unmix
module subroutine hybgen_column_unmix(CS, nk, Rcv_tgt, temp, saln, Rcv, eqn_of_state, &
                               ntr, tracer, trcflg, fixlay, qhrlx, h_col, &
                               terrain_following, h_thin)
  type(hybgen_unmix_CS), intent(in) :: CS  !< hybgen unmixing control structure
  integer,        intent(in)    :: nk           !< The number of layers
  integer,        intent(in)    :: fixlay       !< deepest fixed coordinate layer
  real,           intent(in)    :: qhrlx(nk+1)  !< Relaxation fraction per timestep [nondim], < 1.
  real,           intent(in)    :: Rcv_tgt(nk)  !< Target potential density [R ~> kg m-3]
  real,           intent(inout) :: temp(nk)     !< A column of potential temperature [C ~> degC]
  real,           intent(inout) :: saln(nk)     !< A column of salinity [S ~> ppt]
  real,           intent(inout) :: Rcv(nk)      !< Coordinate potential density [R ~> kg m-3]
  type(EOS_type), intent(in)    :: eqn_of_state !< Equation of state structure
  integer,        intent(in)    :: ntr          !< The number of registered passive tracers
  real,           intent(inout) :: tracer(nk, max(ntr,1)) !< Columns of the passive tracers [Conc]
  integer,        intent(in)    :: trcflg(max(ntr,1)) !< Hycom tracer type flag for each tracer
  real,           intent(inout) :: h_col(nk+1)  !< Layer thicknesses [H ~> m or kg m-2]
  logical,        intent(in)    :: terrain_following !< True if this column is terrain following
  real,           intent(in)    :: h_thin       !< A negligibly small thickness to identify
                                                !! essentially vanished layers [H ~> m or kg m-2]

!
! --- ------------------------------------------------------------------
! --- hybrid grid generator, single column - ummix lowest massive layer.
! --- ------------------------------------------------------------------
!
  ! Local variables
                      ! with temperature [R C-1 ~> kg m-3 degC-1]
                      ! with salinity [R S-1 ~> kg m-3 ppt-1]
                      ! layers by which the source layer's property changes by the loss of water
                      ! that matches the destination layers properties via unmixing [nondim].
                      ! used for updating the concentration of passive tracers [nondim]

  ! --- identify the deepest layer kp with significant thickness (> h_thin)
end subroutine hybgen_column_unmix
  end interface

end module MOM_hybgen_unmix
