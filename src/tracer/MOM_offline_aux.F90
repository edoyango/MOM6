! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> Contains routines related to offline transport of tracers. These routines are likely to be called from
!> the MOM_offline_main module
module MOM_offline_aux

use MOM_debugging,        only : check_column_integrals
use MOM_domains,          only : pass_var, pass_vector, To_All
use MOM_diag_mediator,    only : post_data
use MOM_error_handler,    only : callTree_enter, callTree_leave, MOM_error, FATAL, WARNING, is_root_pe
use MOM_file_parser,      only : get_param, log_version, param_file_type
use MOM_forcing_type,     only : forcing
use MOM_grid,             only : ocean_grid_type
use MOM_io,               only : MOM_read_data, MOM_read_vector, CENTER
use MOM_opacity,          only : optics_type
use MOM_time_manager,     only : time_type, operator(-)
use MOM_unit_scaling,     only : unit_scale_type
use MOM_variables,        only : vertvisc_type
use MOM_verticalGrid,     only : verticalGrid_type
use astronomy_mod,        only : orbital_time, diurnal_solar, daily_mean_solar

implicit none ; private

public update_offline_from_files
public update_offline_from_arrays
public update_h_horizontal_flux
public update_h_vertical_flux
public limit_mass_flux_3d
public distribute_residual_uh_barotropic
public distribute_residual_vh_barotropic
public distribute_residual_uh_upwards
public distribute_residual_vh_upwards
public next_modulo_time
public offline_add_diurnal_sw

#include "MOM_memory.h"


  interface
module subroutine update_h_horizontal_flux(G, GV, uhtr, vhtr, h_pre, h_new)
  type(ocean_grid_type),   intent(in)    :: G     !< ocean grid structure
  type(verticalGrid_type), intent(in)    :: GV    !< ocean vertical grid structure
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                           intent(in)    :: uhtr  !< Accumulated mass flux through zonal face [H L2 ~> m3 or kg]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                           intent(in)    :: vhtr  !< Accumulated mass flux through meridional face [H L2 ~> m3 or kg]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in)    :: h_pre !< Previous layer thicknesses [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(inout) :: h_new !< Updated layer thicknesses [H ~> m or kg m-2]

  ! Local variables
  ! Set index-related variables for fields on T-grid
end subroutine update_h_horizontal_flux
module subroutine update_h_vertical_flux(G, GV, ea, eb, h_pre, h_new)
  type(ocean_grid_type),   intent(in)    :: G     !< ocean grid structure
  type(verticalGrid_type), intent(in)    :: GV    !< ocean vertical grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in)    :: ea    !< Mass of fluid entrained from the layer
                                                  !! above within this timestep [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in)    :: eb    !< Mass of fluid entrained from the layer
                                                  !! below within this timestep [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in)    :: h_pre !< Layer thicknesses at the end of the previous
                                                  !! step [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(inout) :: h_new !< Updated layer thicknesses [H ~> m or kg m-2]

  ! Local variables
  ! Set index-related variables for fields on T-grid
end subroutine update_h_vertical_flux
module subroutine limit_mass_flux_3d(G, GV, uh, vh, ea, eb, h_pre)
  type(ocean_grid_type),   intent(in)    :: G     !< ocean grid structure
  type(verticalGrid_type), intent(in)    :: GV    !< ocean vertical grid structure
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                           intent(inout) :: uh    !< Mass flux through zonal face [H L2 ~> m3 or kg]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                           intent(inout) :: vh    !< Mass flux through meridional face [H L2 ~> m3 or kg]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(inout) :: ea    !< Mass of fluid entrained from the layer
                                                  !! above within this timestep [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(inout) :: eb    !< Mass of fluid entrained from the layer
                                                  !! below within this timestep [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in)    :: h_pre !< Layer thicknesses at the end of the previous
                                                  !! step [H ~> m or kg m-2]

  ! Local variables
                                                           ! top [H ~> m or kg m-2]
                                                           ! bottom [H ~> m or kg m-2]

end subroutine limit_mass_flux_3d
module subroutine distribute_residual_uh_barotropic(G, GV, hvol, uh)
  type(ocean_grid_type),   intent(in   ) :: G    !< ocean grid structure
  type(verticalGrid_type), intent(in   ) :: GV   !< ocean vertical grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in   ) :: hvol !< Mass of water in the cells at the end
                                                 !! of the previous timestep [H L2 ~> m3 or kg]
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                           intent(inout) :: uh   !< Zonal mass transport within a timestep [H L2 ~> m3 or kg]

  ! Local variables


  ! Set index-related variables for fields on T-grid
end subroutine distribute_residual_uh_barotropic
module subroutine distribute_residual_vh_barotropic(G, GV, hvol, vh)
  type(ocean_grid_type),   intent(in   ) :: G    !< ocean grid structure
  type(verticalGrid_type), intent(in   ) :: GV   !< ocean vertical grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in   ) :: hvol !< Mass of water in the cells at the end
                                                 !! of the previous timestep [H L2 ~> m3 or kg]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                           intent(inout) :: vh   !< Meridional mass transport within a timestep [H L2 ~> m3 or kg]

  ! Local variables


  ! Set index-related variables for fields on T-grid
end subroutine distribute_residual_vh_barotropic
module subroutine distribute_residual_uh_upwards(G, GV, hvol, uh)
  type(ocean_grid_type),   intent(in   ) :: G     !< ocean grid structure
  type(verticalGrid_type), intent(in   ) :: GV    !< ocean vertical grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in   ) :: hvol  !< Mass of water in the cells at the end
                                                  !! of the previous timestep [H L2 ~> m3 or kg]
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                           intent(inout) :: uh    !< Zonal mass transport within a timestep [H L2 ~> m3 or kg]

  ! Local variables


  ! Set index-related variables for fields on T-grid
end subroutine distribute_residual_uh_upwards
module subroutine distribute_residual_vh_upwards(G, GV, hvol, vh)
  type(ocean_grid_type),   intent(in   ) :: G     !< ocean grid structure
  type(verticalGrid_type), intent(in   ) :: GV    !< ocean vertical grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in   ) :: hvol  !< Mass of water in the cells at the end
                                                  !! of the previous timestep [H L2 ~> m3 or kg]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                           intent(inout) :: vh    !< Meridional mass transport within a timestep [H L2 ~> m3 or kg]

  ! Local variables


  ! Set index-related variables for fields on T-grid
end subroutine distribute_residual_vh_upwards
module subroutine offline_add_diurnal_SW(fluxes, G, Time_start, Time_end)
  type(forcing),         intent(inout) :: fluxes !< The type with atmospheric fluxes to be adjusted.
  type(ocean_grid_type), intent(in)    :: G      !< The ocean lateral grid type.
  type(time_type),       intent(in)    :: Time_start !< The start time for this step.
  type(time_type),       intent(in)    :: Time_end   !< The ending time for this step.

                         ! the orbital ellipse averaged over a day [nondim]
                         ! the orbital ellipse averaged over a timestep [nondim]


end subroutine offline_add_diurnal_SW
module subroutine update_offline_from_files(G, GV, US, nk_input, mean_file, sum_file, snap_file, &
                surf_file, h_end, uhtr, vhtr, temp_mean, salt_mean, mld, Kd, fluxes, &
                ridx_sum, ridx_snap, read_mld, mld_var_name, read_sw, read_ts_uvh, do_ale_in)

  type(ocean_grid_type),   intent(inout) :: G         !< Horizontal grid type
  type(verticalGrid_type), intent(in   ) :: GV        !< Vertical grid type
  type(unit_scale_type),   intent(in   ) :: US        !< A dimensional unit scaling type
  integer,                 intent(in   ) :: nk_input  !< Number of levels in input file
  character(len=*),        intent(in   ) :: mean_file !< Name of file with averages fields
  character(len=*),        intent(in   ) :: sum_file  !< Name of file with summed fields
  character(len=*),        intent(in   ) :: snap_file !< Name of file with snapshot fields
  character(len=*),        intent(in   ) :: surf_file !< Name of file with surface fields
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(inout) :: h_end     !< End of timestep layer thickness [H ~> m or kg m-2]
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                           intent(inout) :: uhtr      !< Zonal mass fluxes [H L2 ~> m3 or kg]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                           intent(inout) :: vhtr      !< Meridional mass fluxes [H L2 ~> m3 or kg]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(inout) :: temp_mean !< Averaged temperature [C ~> degC]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(inout) :: salt_mean !< Averaged salinity [S ~> ppt]
  real, dimension(SZI_(G),SZJ_(G)),          &
                           intent(inout) :: mld       !< Averaged mixed layer depth [Z ~> m]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)+1), &
                           intent(inout) :: Kd        !< Diapycnal diffusivities at interfaces
                                                      !! [H Z T-1 ~> m2 s-1 or kg m-1 s-1]
  type(forcing),           intent(inout) :: fluxes    !< Fields with surface fluxes
  integer,                 intent(in   ) :: ridx_sum  !< Read index for sum, mean, and surf files
  integer,                 intent(in   ) :: ridx_snap !< Read index for snapshot file
  logical,                 intent(in   ) :: read_mld  !< True if reading in MLD
  character(len=*),        intent(in   ) :: mld_var_name !< Name of the mixed layer depth variable
                                                      !! to read from a file.
  logical,                 intent(in   ) :: read_sw   !< True if reading in radiative fluxes
  logical,                 intent(in   ) :: read_ts_uvh !< True if reading in uh, vh, and h
  logical,       optional, intent(in   ) :: do_ale_in !< True if using ALE algorithms

                           ! file to H [H m-1 ~> 1] or [H m2 kg-1 ~> 1]

end subroutine update_offline_from_files
module subroutine update_offline_from_arrays(G, GV, nk_input, ridx_sum, mean_file, sum_file, snap_file, uhtr, vhtr, &
                                      hend, uhtr_all, vhtr_all, hend_all, temp, salt, temp_all, salt_all )
  type(ocean_grid_type),                     intent(inout) :: G         !< Horizontal grid type
  type(verticalGrid_type),                   intent(in   ) :: GV        !< Vertical grid type
  integer,                                   intent(in   ) :: nk_input  !< Number of levels in input file
  integer,                                   intent(in   ) :: ridx_sum  !< Index to read from
  character(len=200),                        intent(in   ) :: mean_file !< Name of file with averages fields
  character(len=200),                        intent(in   ) :: sum_file  !< Name of file with summed fields
  character(len=200),                        intent(in   ) :: snap_file !< Name of file with snapshot fields
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), intent(inout) :: uhtr     !< Zonal mass fluxes [H L2 ~> m3 or kg]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), intent(inout) :: vhtr     !< Meridional mass fluxes [H L2 ~> m3 or kg]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(inout) :: hend      !< End of timestep layer thickness
                                                                        !! [H ~> m or kg m-2]
  real, dimension(:,:,:,:), allocatable,     intent(inout) :: uhtr_all  !< Zonal mass fluxes [H L2 ~> m3 or kg]
  real, dimension(:,:,:,:), allocatable,     intent(inout) :: vhtr_all  !< Meridional mass fluxes [H L2 ~> m3 or kg]
  real, dimension(:,:,:,:), allocatable,     intent(inout) :: hend_all  !< End of timestep layer thickness
                                                                        !! [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(inout) :: temp      !< Temperature array [C ~> degC]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(inout) :: salt      !< Salinity array [S ~> ppt]
  real, dimension(:,:,:,:), allocatable,     intent(inout) :: temp_all  !< Temperature array [C ~> degC]
  real, dimension(:,:,:,:), allocatable,     intent(inout) :: salt_all  !< Salinity array [S ~> ppt]

end subroutine update_offline_from_arrays
module function next_modulo_time(inidx, numtime)
  ! Returns the next time interval to be read
  integer                 :: numtime              ! Number of time levels in input fields
  integer                 :: inidx                ! The current time index

                                                  ! to the current timestep

  integer                 :: next_modulo_time

end function next_modulo_time
  end interface

end module MOM_offline_aux
