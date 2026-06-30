! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> Interface height filtering module
module MOM_interface_filter

use MOM_debugging,             only : hchksum, uvchksum
use MOM_diag_mediator,         only : post_data, query_averaging_enabled, diag_ctrl
use MOM_diag_mediator,         only : register_diag_field, safe_alloc_ptr, time_type
use MOM_domains,               only : pass_var, CORNER, pass_vector
use MOM_error_handler,         only : MOM_error, FATAL, WARNING, is_root_pe
use MOM_file_parser,           only : get_param, log_version, param_file_type
use MOM_grid,                  only : ocean_grid_type
use MOM_interface_heights,     only : find_eta
use MOM_unit_scaling,          only : unit_scale_type
use MOM_variables,             only : thermo_var_ptrs, cont_diag_ptrs
use MOM_verticalGrid,          only : verticalGrid_type

implicit none ; private

#include <MOM_memory.h>

public interface_filter, interface_filter_init, interface_filter_end

! A note on unit descriptions in comments: MOM6 uses units that can be rescaled for dimensional
! consistency testing. These are noted in comments with units like Z, H, L, and T, along with
! their mks counterparts with notation like "a velocity [Z T-1 ~> m s-1]".  If the units
! vary with the Boussinesq approximation, the Boussinesq variant is given first.

!> Control structure for interface height filtering
type, public :: interface_filter_CS ; private
  logical :: initialized = .false. !< True if this control structure has been initialized.
  real    :: max_smoothing_CFL   !< Maximum value of the smoothing CFL for interface height filtering [nondim]
  real    :: filter_rate         !< The rate at which grid-scale anomalies are damped away [T-1 ~> s-1]
  integer :: filter_order        !< The even power of the interface height smoothing.
                                 !! At present valid values are 0, 2, or 4.
  logical :: interface_filter    !< If true, interfaces heights are diffused.
  logical :: isotropic_filter    !< If true, use the same filtering lengthscales in both directions,
                                 !! otherwise use filtering lengthscales in each direction that scale
                                 !! with the grid spacing in that direction.
  logical :: debug               !< write verbose checksums for debugging purposes

  type(diag_ctrl), pointer :: diag => NULL() !< structure used to regulate timing of diagnostics

  !>@{
  !! Diagnostic identifier
  integer :: id_uh_sm  = -1, id_vh_sm  = -1
  integer :: id_L2_u  = -1, id_L2_v  = -1
  integer :: id_sfn_x = -1, id_sfn_y = -1
  !>@}
end type interface_filter_CS


  interface
module subroutine interface_filter(h, uhtr, vhtr, tv, dt, G, GV, US, CDp, CS)
  type(ocean_grid_type),                      intent(in)    :: G      !< Ocean grid structure
  type(verticalGrid_type),                    intent(in)    :: GV     !< Vertical grid structure
  type(unit_scale_type),                      intent(in)    :: US     !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),  intent(inout) :: h      !< Layer thickness [H ~> m or kg m-2]
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), intent(inout) :: uhtr   !< Accumulated zonal mass flux
                                                                      !! [L2 H ~> m3 or kg]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), intent(inout) :: vhtr   !< Accumulated meridional mass flux
                                                                      !! [L2 H ~> m3 or kg]
  type(thermo_var_ptrs),                      intent(in)    :: tv     !< Thermodynamics structure
  real,                                       intent(in)    :: dt     !< Time increment [T ~> s]
  type(cont_diag_ptrs),                       intent(inout) :: CDp    !< Diagnostics for the continuity equation
  type(interface_filter_CS),                  intent(inout) :: CS     !< Control structure for interface height
                                                                      !! filtering
  ! Local variables
                                         ! sea level [Z ~> m], positive up.
                                         ! of Laplacian smoothing [Z ~> m], positive downward to avoid
                                         ! having to change other signs in the call to interface_filter.


                                                        ! [H L2 T-1 ~> m3 s-1 or kg s-1]
                                                        ! [H L2 T-1 ~> m3 s-1 or kg s-1]
                    ! in roundoff and can be neglected [H ~> m or kg m-2].

end subroutine interface_filter
module subroutine filter_interface(h, e, Lsm2_u, Lsm2_v, uhD, vhD, tv, G, GV, US, halo_size)
  type(ocean_grid_type),                       intent(in)  :: G     !< Ocean grid structure
  type(verticalGrid_type),                     intent(in)  :: GV    !< Vertical grid structure
  type(unit_scale_type),                       intent(in)  :: US    !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),   intent(in)  :: h     !< Layer thickness [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)+1), intent(in)  :: e     !< Interface positions [Z ~> m]
  real, dimension(SZIB_(G),SZJ_(G)),           intent(in)  :: Lsm2_u !< Interface smoothing lengths squared
                                                                    !! at u points [L2 ~> m2]
  real, dimension(SZI_(G),SZJB_(G)),           intent(in)  :: Lsm2_v !< Interface smoothing lengths squared
                                                                    !! at v points [L2 ~> m2]
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)),  intent(out) :: uhD   !< Zonal mass fluxes
                                                                    !! [H L2 ~> m3 or kg]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)),  intent(out) :: vhD   !< Meridional mass fluxes
                                                                    !! [H L2 ~> m3 or kg]
  type(thermo_var_ptrs),                       intent(in)  :: tv     !< Thermodynamics structure
  integer,                           optional, intent(in)  :: halo_size !< The size of the halo to work on,
                                                                    !! 0 by default.

  ! Local variables
                        ! between -1 and 1 after undoing dimensional scaling, [Z L-1 ~> nondim]
                        ! streamfunction [H L2 ~> m3 or kg].
                        ! in roundoff and can be neglected [H ~> m or kg m-2].

end subroutine filter_interface
module subroutine interface_filter_init(Time, G, GV, US, param_file, diag, CDp, CS)
  type(time_type),         intent(in) :: Time    !< Current model time
  type(ocean_grid_type),   intent(in) :: G       !< Ocean grid structure
  type(verticalGrid_type), intent(in) :: GV      !< Vertical grid structure
  type(unit_scale_type),   intent(in) :: US      !< A dimensional unit scaling type
  type(param_file_type),   intent(in) :: param_file !< Parameter file handles
  type(diag_ctrl), target, intent(inout) :: diag !< Diagnostics control structure
  type(cont_diag_ptrs),    intent(inout) :: CDp  !< Continuity equation diagnostics
  type(interface_filter_CS), intent(inout) :: CS !< Control structure for interface height filtering

  ! Local variables
  ! This include declares and sets the variable "version".

end subroutine interface_filter_init
module subroutine interface_filter_end(CS, CDp)
  type(interface_filter_CS), intent(inout) :: CS !< Control structure for interface height filtering
  type(cont_diag_ptrs), intent(inout) :: CDp      !< Continuity diagnostic control structure

  ! NOTE: [uv]h_smooth are not yet used in diagnostics, but they are here for now for completeness.
end subroutine interface_filter_end
  end interface

end module MOM_interface_filter
