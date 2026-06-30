! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> Initialize ice shelf variables
module MOM_ice_shelf_initialize

use MOM_grid, only : ocean_grid_type
use MOM_array_transform,      only : rotate_array
use MOM_hor_index,  only : hor_index_type
use MOM_file_parser, only : get_param, read_param, log_param, param_file_type
use MOM_io, only: MOM_read_data, file_exists, field_exists, slasher, CORNER
use MOM_error_handler, only : MOM_error, MOM_mesg, FATAL, WARNING, is_root_pe
use MOM_unit_scaling, only : unit_scale_type
use user_shelf_init, only: USER_init_ice_thickness

implicit none ; private

#include <MOM_memory.h>

public initialize_ice_thickness
public initialize_ice_shelf_boundary_channel
public initialize_ice_flow_from_file
public initialize_ice_shelf_boundary_from_file
public initialize_ice_C_basal_friction
public initialize_ice_AGlen
public initialize_ice_SMB
! A note on unit descriptions in comments: MOM6 uses units that can be rescaled for dimensional
! consistency testing. These are noted in comments with units like Z, H, L, and T, along with
! their mks counterparts with notation like "a velocity [Z T-1 ~> m s-1]".  If the units
! vary with the Boussinesq approximation, the Boussinesq variant is given first.


  interface
module subroutine initialize_ice_thickness(h_shelf, area_shelf_h, hmask, melt_mask, G, G_in, US, PF, rotate_index, turns)
  type(ocean_grid_type), intent(in)    :: G    !< The ocean's grid structure
  type(ocean_grid_type), intent(in)    :: G_in    !< The ocean's unrotated grid structure
  real, dimension(SZDI_(G),SZDJ_(G)), &
                         intent(inout) :: h_shelf !< The ice shelf thickness [Z ~> m].
  real, dimension(SZDI_(G),SZDJ_(G)), &
                         intent(inout) :: area_shelf_h !< The area per cell covered by the ice shelf [L2 ~> m2].
  real, dimension(SZDI_(G),SZDJ_(G)), &
                         intent(inout) :: hmask !< A mask indicating which tracer points are
                                             !! partly or fully covered by an ice-shelf [nondim]
  real, dimension(SZDI_(G),SZDJ_(G)), &
                         intent(inout) :: melt_mask !< A mask indicating where to allow ice-shelf melting [nondim]
  type(unit_scale_type), intent(in)    :: US !< A structure containing unit conversion factors
  type(param_file_type), intent(in)    :: PF !< A structure to parse for run-time parameters
  logical, intent(in), optional        :: rotate_index !< If true, this is a rotation test
  integer, intent(in), optional        :: turns !< Number of turns for rotation test


end subroutine initialize_ice_thickness
module subroutine initialize_ice_thickness_from_file(h_shelf, area_shelf_h, hmask, melt_mask, G, US, PF)
  type(ocean_grid_type), intent(in)    :: G    !< The ocean's grid structure
  real, dimension(SZDI_(G),SZDJ_(G)), &
                         intent(inout) :: h_shelf !< The ice shelf thickness [Z ~> m].
  real, dimension(SZDI_(G),SZDJ_(G)), &
                         intent(inout) :: area_shelf_h !< The area per cell covered by the ice shelf [L2 ~> m2].
  real, dimension(SZDI_(G),SZDJ_(G)), &
                         intent(inout) :: hmask !< A mask indicating which tracer points are
                                             !! partly or fully covered by an ice-shelf [nondim]
  real, dimension(SZDI_(G),SZDJ_(G)), &
                         intent(inout) :: melt_mask !< A mask indicating where to allow ice-shelf melting [nondim]
  type(unit_scale_type), intent(in)    :: US !< A structure containing unit conversion factors
  type(param_file_type), intent(in)    :: PF !< A structure to parse for run-time parameters

  !  This subroutine reads ice thickness and area from a file and puts it into
  !  h_shelf [Z ~> m] and area_shelf_h [L2 ~> m2] (and dimensionless) and updates hmask

end subroutine initialize_ice_thickness_from_file
module subroutine initialize_ice_thickness_channel(h_shelf, area_shelf_h, hmask, G, US, PF)
  type(ocean_grid_type), intent(in)    :: G    !< The ocean's grid structure
  real, dimension(SZDI_(G),SZDJ_(G)), &
                         intent(inout) :: h_shelf !< The ice shelf thickness [Z ~> m].
  real, dimension(SZDI_(G),SZDJ_(G)), &
                         intent(inout) :: area_shelf_h !< The area per cell covered by the ice shelf [L2 ~> m2].
  real, dimension(SZDI_(G),SZDJ_(G)), &
                         intent(inout) :: hmask !< A mask indicating which tracer points are
                                             !! partly or fully covered by an ice-shelf
  type(unit_scale_type), intent(in)    :: US !< A structure containing unit conversion factors
  type(param_file_type), intent(in)    :: PF !< A structure to parse for run-time parameters


end subroutine initialize_ice_thickness_channel
module subroutine initialize_ice_shelf_boundary_channel(u_face_mask_bdry, v_face_mask_bdry, &
                u_flux_bdry_val, v_flux_bdry_val, u_bdry_val, v_bdry_val, u_shelf, v_shelf, h_bdry_val, &
                hmask,  h_shelf, G, US, PF )

  type(ocean_grid_type), intent(in)    :: G    !< The ocean's grid structure
  real, dimension(SZIB_(G),SZJB_(G)), &
                         intent(inout) :: u_face_mask_bdry !< A boundary-type mask at C-grid u faces

  real, dimension(SZIB_(G),SZJ_(G)), &
                         intent(inout) :: u_flux_bdry_val  !< The boundary thickness flux through
                                                     !! C-grid u faces [L Z T-1 ~> m2 s-1].
  real, dimension(SZIB_(G),SZJB_(G)), &
                         intent(inout) :: v_face_mask_bdry !< A boundary-type mask at C-grid v faces

  real, dimension(SZI_(G),SZJB_(G)), &
                         intent(inout) :: v_flux_bdry_val  !< The boundary thickness flux through
                                                     !! C-grid v faces [L Z T-1 ~> m2 s-1].
  real, dimension(SZIB_(G),SZJB_(G)), &
                         intent(inout) :: u_bdry_val !< The zonal ice shelf velocity at open
                                                      !! boundary vertices [L T-1 ~> m s-1].
  real, dimension(SZIB_(G),SZJB_(G)), &
                         intent(inout) :: v_bdry_val !< The meridional ice shelf velocity at open
  real, dimension(SZIB_(G),SZJB_(G)), &
                         intent(inout) :: u_shelf !< The zonal ice shelf velocity  [L T-1 ~> m s-1].
  real, dimension(SZIB_(G),SZJB_(G)), &
                         intent(inout) :: v_shelf !< The meridional ice shelf velocity  [L T-1 ~> m s-1].
  real, dimension(SZDI_(G),SZDJ_(G)), &
                         intent(inout) :: h_bdry_val !< The ice shelf thickness at open boundaries [Z ~> m]
  real, dimension(SZDI_(G),SZDJ_(G)), &
                         intent(inout) :: hmask !< A mask indicating which tracer points are
                                             !! partly or fully covered by an ice-shelf
  real, dimension(SZDI_(G),SZDJ_(G)), &
                         intent(inout) :: h_shelf !< Ice-shelf thickness [Z ~> m]
  type(unit_scale_type), intent(in)    :: US !< A structure containing unit conversion factors
  type(param_file_type), intent(in)    :: PF !< A structure to parse for run-time parameters


end subroutine initialize_ice_shelf_boundary_channel
module subroutine initialize_ice_flow_from_file(bed_elev,u_shelf, v_shelf,float_cond,&
                                         G, US, PF)
  type(ocean_grid_type), intent(in)    :: G    !< The ocean's grid structure
  real, dimension(SZDI_(G),SZDJ_(G)), &
                         intent(inout) :: bed_elev !< The bed elevation   [Z ~> m].
  real, dimension(SZIB_(G),SZJB_(G)), &
                         intent(inout) :: u_shelf !< The zonal ice shelf velocity  [L T-1 ~> m s-1].
  real, dimension(SZIB_(G),SZJB_(G)), &
                         intent(inout) :: v_shelf !< The meridional ice shelf velocity  [L T-1 ~> m s-1].
  real, dimension(SZDI_(G),SZDJ_(G)), &
                         intent(inout)    :: float_cond !< An array indicating where the ice
                                                !! shelf is floating: 0 if floating, 1 if not. [nondim]
  type(unit_scale_type), intent(in)    :: US !< A structure containing unit conversion factors
  type(param_file_type), intent(in)    :: PF !< A structure to parse for run-time parameters

  !  This subroutine reads ice thickness and area from a file and puts it into
  !  h_shelf [Z ~> m] and area_shelf_h [L2 ~> m2] (and dimensionless) and updates hmask

end subroutine initialize_ice_flow_from_file
module subroutine initialize_ice_shelf_boundary_from_file(u_face_mask_bdry, v_face_mask_bdry, &
                u_bdry_val, v_bdry_val, umask, vmask, h_bdry_val, &
                hmask,  h_shelf, G, US, PF )

  type(ocean_grid_type), intent(in)    :: G    !< The ocean's grid structure
  real, dimension(SZIB_(G),SZJB_(G)), &
                         intent(inout) :: u_face_mask_bdry !< A boundary-type mask at B-grid u faces [nondim]
  real, dimension(SZIB_(G),SZJB_(G)), &
                         intent(inout) :: v_face_mask_bdry !< A boundary-type mask at B-grid v faces [nondim]
  real, dimension(SZIB_(G),SZJB_(G)), &
                         intent(inout) :: u_bdry_val !< The zonal ice shelf velocity at open
                                                      !! boundary vertices [L T-1 ~> m s-1].
  real, dimension(SZIB_(G),SZJB_(G)), &
                         intent(inout) :: v_bdry_val !< The meridional ice shelf velocity at open
                                                      !! boundary vertices [L T-1 ~> m s-1].
  real, dimension(SZDIB_(G),SZDJB_(G)), &
                         intent(inout) :: umask !< A mask for ice shelf velocity [nondim]
  real, dimension(SZDIB_(G),SZDJB_(G)), &
                         intent(inout) :: vmask !< A mask for ice shelf velocity [nondim]
  real, dimension(SZDI_(G),SZDJ_(G)), &
                         intent(inout) :: h_bdry_val !< The ice shelf thickness at open boundaries [Z ~> m]
  real, dimension(SZDI_(G),SZDJ_(G)), &
                         intent(inout) :: hmask !< A mask indicating which tracer points are
                                             !! partly or fully covered by an ice-shelf [nondim]
  real, dimension(SZDI_(G),SZDJ_(G)), &
                         intent(in) :: h_shelf !< Ice-shelf thickness [Z ~> m]
  type(unit_scale_type), intent(in)    :: US !< A structure containing unit conversion factors
  type(param_file_type), intent(in)    :: PF !< A structure to parse for run-time parameters



end subroutine initialize_ice_shelf_boundary_from_file
module subroutine initialize_ice_C_basal_friction(C_basal_friction, G, US, PF)
  type(ocean_grid_type), intent(in)    :: G    !< The ocean's grid structure
  real, dimension(SZDI_(G),SZDJ_(G)), &
                         intent(inout) :: C_basal_friction !< Ice-stream basal friction
                                             !! in units of [R L Z T-2 (s m-1)^n_basal_fric ~> Pa (s m-1)^n_basal_fric]
  type(unit_scale_type), intent(in)    :: US !< A structure containing unit conversion factors
  type(param_file_type), intent(in)    :: PF !< A structure to parse for run-time parameters

!  integer :: i, j
                      ! [R L Z T-2 (s m-1)^n_basal_fric ~> Pa (s m-1)^n_basal_fric]

end subroutine initialize_ice_C_basal_friction
module subroutine initialize_ice_AGlen(AGlen, ice_viscosity_compute, G, US, PF)
  type(ocean_grid_type), intent(in)    :: G    !< The ocean's grid structure
  real, dimension(SZDI_(G),SZDJ_(G)), &
                         intent(inout) :: AGlen !< The ice-stiffness parameter A_Glen, often in [Pa-3 s-1]
  character(len=40) :: ice_viscosity_compute !< Specifies whether the ice viscosity is computed internally
                                             !! according to Glen's flow law; is constant (for debugging purposes)
                                             !! or using observed strain rates and read from a file
  type(unit_scale_type), intent(in)    :: US !< A structure containing unit conversion factors
  type(param_file_type), intent(in)    :: PF !< A structure to parse for run-time parameters


end subroutine initialize_ice_AGlen
module subroutine initialize_ice_SMB(SMB, G, US, PF)
  type(ocean_grid_type), intent(in)    :: G    !< The ocean's grid structure
  real, dimension(SZDI_(G),SZDJ_(G)), &
                         intent(inout) :: SMB !< Ice surface mass balance parameter, often in [R Z T-1 ~> kg m-2 s-1]
  type(unit_scale_type), intent(in)    :: US !< A structure containing unit conversion factors
  type(param_file_type), intent(in)    :: PF !< A structure to parse for run-time parameters


end subroutine initialize_ice_SMB
  end interface

end module MOM_ice_shelf_initialize
