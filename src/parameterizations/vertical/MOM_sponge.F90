! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> Implements sponge regions in isopycnal mode
module MOM_sponge

use MOM_coms,          only : sum_across_PEs
use MOM_diag_mediator, only : post_data, query_averaging_enabled, register_diag_field
use MOM_diag_mediator, only : diag_ctrl
use MOM_error_handler, only : MOM_error, FATAL, NOTE, WARNING, is_root_pe
use MOM_file_parser,   only : get_param, log_param, log_version, param_file_type
use MOM_grid,          only : ocean_grid_type
use MOM_spatial_means, only : global_i_mean
use MOM_time_manager,  only : time_type
use MOM_unit_scaling,  only : unit_scale_type
use MOM_variables,     only : thermo_var_ptrs
use MOM_verticalGrid,  only : verticalGrid_type

! Planned extension:  Support for time varying sponge targets.

implicit none ; private

#include <MOM_memory.h>

public set_up_sponge_field, set_up_sponge_ML_density
public initialize_sponge, apply_sponge, sponge_end, init_sponge_diags

! A note on unit descriptions in comments: MOM6 uses units that can be rescaled for dimensional
! consistency testing. These are noted in comments with units like Z, H, L, and T, along with
! their mks counterparts with notation like "a velocity [Z T-1 ~> m s-1]".  If the units
! vary with the Boussinesq approximation, the Boussinesq variant is given first.

!> A structure for creating arrays of pointers to 3D arrays
type, public :: p3d
  real, dimension(:,:,:), pointer :: p => NULL() !< A pointer to a 3D array [various]
end type p3d
!> A structure for creating arrays of pointers to 2D arrays
type, public :: p2d
  real, dimension(:,:), pointer :: p => NULL() !< A pointer to a 2D array [various]
end type p2d

!> This control structure holds memory and parameters for the MOM_sponge module
type, public :: sponge_CS ; private
  logical :: bulkmixedlayer  !< If true, a refined bulk mixed layer is used with
                       !! nkml sublayers and nkbl buffer layer.
  integer :: nz        !< The total number of layers.
  integer :: num_col   !< The number of sponge points within the computational domain.
  integer :: fldno = 0 !< The number of fields which have already been
                       !! registered by calls to set_up_sponge_field
  integer, pointer :: col_i(:) => NULL() !< Array of the i-indicies of each of the columns being damped.
  integer, pointer :: col_j(:) => NULL() !< Array of the j-indicies of each of the columns being damped.
  real, pointer :: Iresttime_col(:) => NULL() !< The inverse restoring time of each column [T-1 ~> s-1].
  real, pointer :: Rcv_ml_ref(:) => NULL() !< The value toward which the mixed layer
                             !! coordinate-density is being damped [R ~> kg m-3].
  real, pointer :: Ref_eta(:,:) => NULL() !< The value toward which the interface
                             !! heights are being damped [Z ~> m].
  type(p3d) :: var(MAX_FIELDS_) !< Pointers to the fields that are being damped.
  type(p2d) :: Ref_val(MAX_FIELDS_) !< The values to which the fields are damped.

  logical :: do_i_mean_sponge !< If true, apply sponges to the i-mean fields.
  real, pointer :: Iresttime_im(:) => NULL() !< The inverse restoring time of
                             !! each row for i-mean sponges [T-1 ~> s-1].
  real, pointer :: Rcv_ml_ref_im(:) => NULL() !! The value toward which the i-mean
                             !< mixed layer coordinate-density is being damped [R ~> kg m-3].
  real, pointer :: Ref_eta_im(:,:) => NULL() !< The value toward which the i-mean
                             !! interface heights are being damped [Z ~> m].
  type(p2d) :: Ref_val_im(MAX_FIELDS_) !< The values toward which the i-means of
                             !! fields are damped.

  type(diag_ctrl), pointer :: diag => NULL() !< A structure that is used to
                             !! regulate the timing of diagnostic output.
  integer :: id_w_sponge = -1 !< A diagnostic ID
end type sponge_CS


  interface
module subroutine initialize_sponge(Iresttime, int_height, G, param_file, CS, GV, &
                             Iresttime_i_mean, int_height_i_mean)
  type(ocean_grid_type),   intent(in) :: G          !< The ocean's grid structure
  type(verticalGrid_type), intent(in) :: GV         !< The ocean's vertical grid structure
  real, dimension(SZI_(G),SZJ_(G)), &
                           intent(in) :: Iresttime  !< The inverse of the restoring time [T-1 ~> s-1].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)+1), &
                           intent(in) :: int_height !< The interface heights to damp back toward [Z ~> m].
  type(param_file_type),   intent(in) :: param_file !< A structure to parse for run-time parameters
  type(sponge_CS),         pointer    :: CS         !< A pointer that is set to point to the control
                                                    !! structure for this module
  real, dimension(SZJ_(G)), &
                 optional, intent(in) :: Iresttime_i_mean !< The inverse of the restoring time for
                                                          !! the zonal mean properties [T-1 ~> s-1].
  real, dimension(SZJ_(G),SZK_(GV)+1), &
                 optional, intent(in) :: int_height_i_mean !< The interface heights toward which to
                                                           !! damp the zonal mean heights [Z ~> m].


  ! This include declares and sets the variable "version".

end subroutine initialize_sponge
module subroutine init_sponge_diags(Time, G, GV, US, diag, CS)
  type(time_type),       target, intent(in)    :: Time !< The current model time
  type(ocean_grid_type),         intent(in)    :: G    !< The ocean's grid structure
  type(verticalGrid_type),       intent(in)    :: GV   !< The ocean's vertical grid structure
  type(unit_scale_type),         intent(in)    :: US   !< A dimensional unit scaling type
  type(diag_ctrl),       target, intent(inout) :: diag !< A structure that is used to regulate diagnostic output
  type(sponge_CS),               pointer       :: CS   !< A pointer to the control structure for this module that
                                                       !! is set by a previous call to initialize_sponge.

end subroutine init_sponge_diags
module subroutine set_up_sponge_field(sp_val, f_ptr, G, GV, nlay, CS, sp_val_i_mean)
  type(ocean_grid_type),   intent(in) :: G      !< The ocean's grid structure
  type(verticalGrid_type), intent(in) :: GV     !< The ocean's vertical grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in) :: sp_val !< The reference profiles of the quantity being registered [various]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                   target, intent(in) :: f_ptr  !< a pointer to the field which will be damped [various]
  integer,                 intent(in) :: nlay   !< the number of layers in this quantity
  type(sponge_CS),         pointer    :: CS     !< A pointer to the control structure for this module that
                                                !! is set by a previous call to initialize_sponge.
  real, dimension(SZJ_(G),SZK_(GV)),&
               optional, intent(in) :: sp_val_i_mean !< The i-mean reference value for
                                              !! this field with i-mean sponges [various]


end subroutine set_up_sponge_field
module subroutine set_up_sponge_ML_density(sp_val, G, CS, sp_val_i_mean)
  type(ocean_grid_type), intent(in) :: G    !< The ocean's grid structure
  real, dimension(SZI_(G),SZJ_(G)), &
                         intent(in) :: sp_val !< The reference values of the mixed layer density [R ~> kg m-3]
  type(sponge_CS),       pointer    :: CS   !< A pointer to the control structure for this module that is
                                            !! set by a previous call to initialize_sponge.
                                            !  The contents of this structure are intent(inout) here.
  real, dimension(SZJ_(G)), &
               optional, intent(in) :: sp_val_i_mean !< the reference values of the zonal mean mixed
                                            !! layer density [R ~> kg m-3], for use if Iresttime_i_mean > 0.


end subroutine set_up_sponge_ML_density
module subroutine apply_sponge(h, tv, dt, G, GV, US, ea, eb, CS, Rcv_ml)
  type(ocean_grid_type),   intent(inout) :: G   !< The ocean's grid structure
  type(verticalGrid_type), intent(in)    :: GV  !< The ocean's vertical grid structure
  type(unit_scale_type),   intent(in)    :: US  !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(inout) :: h   !< Layer thicknesses [H ~> m or kg m-2]
  type(thermo_var_ptrs),   intent(in)    :: tv  !< A structure pointing to various
                                                !! thermodynamic variables
  real,                    intent(in)    :: dt  !< The amount of time covered by this call [T ~> s].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(inout) :: ea  !< An array to which the amount of fluid entrained
                                                !! from the layer above during this call will be
                                                !! added [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(inout) :: eb  !< An array to which the amount of fluid entrained
                                                !! from the layer below during this call will be
                                                !! added [H ~> m or kg m-2].
  type(sponge_CS),         pointer       :: CS  !< A pointer to the control structure for this module
                                                !! that is set by a previous call to initialize_sponge.
  real, dimension(SZI_(G),SZJ_(G)), &
                 optional, intent(inout) :: Rcv_ml !<  The coordinate density of the mixed layer [R ~> kg m-3].

  ! Local variables
end subroutine apply_sponge
module subroutine sponge_end(CS)
  type(sponge_CS),         pointer    :: CS   !< A pointer to the control structure for this module
                                              !! that is set by a previous call to initialize_sponge.

end subroutine sponge_end
  end interface

end module MOM_sponge
