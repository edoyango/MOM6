! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

module MOM_self_attr_load

use MOM_cpu_clock,           only : cpu_clock_id, cpu_clock_begin, cpu_clock_end, CLOCK_MODULE
use MOM_domains,             only : pass_var
use MOM_error_handler,       only : MOM_error, FATAL, WARNING
use MOM_file_parser,         only : get_param, log_param, log_version, param_file_type
use MOM_grid,                only : ocean_grid_type
use MOM_interface_heights,   only : find_col_mass
use MOM_io,                  only : MOM_infra_file, MOM_field, vardesc, slasher
use MOM_io,                  only : create_MOM_file, MOM_read_data, MOM_write_field, var_desc
use MOM_load_love_numbers,   only : Love_Data
use MOM_restart,             only : is_new_run, MOM_restart_CS
use MOM_spherical_harmonics, only : spherical_harmonics_init, spherical_harmonics_end
use MOM_spherical_harmonics, only : spherical_harmonics_forward, spherical_harmonics_inverse
use MOM_spherical_harmonics, only : sht_CS, order2index, calc_lmax
use MOM_string_functions,    only : lowercase
use MOM_unit_scaling,        only : unit_scale_type
use MOM_variables,           only : thermo_var_ptrs
use MOM_verticalGrid,        only : verticalGrid_type

implicit none ; private

public calc_SAL, scalar_SAL_sensitivity, SAL_init, SAL_end

#include <MOM_memory.h>

!> The control structure for the MOM_self_attr_load module
type, public :: SAL_CS ; private
  logical :: use_sal_scalar = .false.
    !< If true, use the scalar approximation to calculate SAL.
  logical :: use_sal_sht = .false.
    !< If true, use online spherical harmonics to calculate SAL
  logical :: use_tidal_sal_prev = .false.
    !< If true, read the tidal SAL from the previous iteration of the tides to
    !! facilitate convergence.
  logical :: use_bpa = .false.
    !< If true, use bottom pressure anomaly instead of SSH to calculate SAL.
  real :: eta_prop
    !< The partial derivative of eta_sal with the local value of eta [nondim].
  real :: linear_scaling
    !< Dimensional coefficients for scalar SAL [nondim] or [Z T2 L-2 R-1 ~> m Pa-1]
  type(sht_CS), allocatable :: sht
    !< Spherical harmonic transforms (SHT) control structure
  integer :: sal_sht_Nd
    !< Maximum degree for spherical harmonic transforms [nondim]
  real, allocatable :: pbot_ref(:,:)
    !< Reference bottom pressure [R L2 T-2 ~> Pa]
  real, allocatable :: Love_scaling(:)
    !< Dimensional coefficients for harmonic SAL, which are functions of Love numbers
    !! [nondim] or [Z T2 L-2 R-1 ~> m Pa-1], depending on the value of use_ppa.
  real, allocatable :: Snm_Re(:), &    !< Real SHT coefficient for SHT SAL [Z ~> m]
                       Snm_Im(:)       !< Imaginary SHT coefficient for SHT SAL [Z ~> m]
end type SAL_CS

integer :: id_clock_SAL   !< CPU clock for self-attraction and loading


  interface
module subroutine calc_SAL(eta, eta_sal, G, CS, tmp_scale)
  type(ocean_grid_type), intent(in) :: G  !< The ocean's grid structure.
  real, dimension(SZI_(G),SZJ_(G)), intent(in)  :: eta     !< The sea surface height anomaly from
              !! a time-mean geoid or total bottom pressure [Z ~> m] or [R L2 T-2 ~> Pa].
  real, dimension(SZI_(G),SZJ_(G)), intent(out) :: eta_sal !< The geopotential height anomaly from
              !! self-attraction and loading [Z ~> m].
  type(SAL_CS), intent(inout) :: CS !< The control structure returned by a previous call to SAL_init.
  real, optional, intent(in)  :: tmp_scale !< A rescaling factor to temporarily convert eta
              !! to MKS units in reproducing sumes [m Z-1 ~> 1]

  ! Local variables

end subroutine calc_SAL
module subroutine scalar_SAL_sensitivity(CS, deta_sal_deta)
  type(SAL_CS), intent(in)  :: CS !< The control structure returned by a previous call to SAL_init.
  real,         intent(out) :: deta_sal_deta !< The partial derivative of eta_sal with
                                             !! the local value of eta [nondim].
end subroutine scalar_SAL_sensitivity
module subroutine calc_love_scaling(rhoW, rhoE, grav, CS)
  real, intent(in) :: rhoW !< The average density of sea water [R ~> kg m-3]
  real, intent(in) :: rhoE !< The average density of Earth [R ~> kg m-3]
  real, intent(in) :: grav !< The gravitational acceleration [L2 Z-1 T-2 ~> m s-2]
  type(SAL_CS), intent(inout) :: CS !< The control structure returned by a previous call to SAL_init.

  ! Local variables
      ! and coef_rhoE = 1.0 / (rhoE * grav) with USE_BPA=True. [nondim] or [Z T2 L-2 R-1 ~> m Pa-1]

end subroutine calc_love_scaling
module subroutine SAL_init(h, tv, G, GV, US, param_file, CS, restart_CS)
  type(ocean_grid_type),          intent(in)    :: G  !< The ocean's grid structure.
  type(verticalGrid_type),        intent(in)    :: GV !< Vertical grid structure
  type(unit_scale_type),          intent(in)    :: US !< A dimensional unit scaling type
  type(param_file_type),          intent(in)    :: param_file !< A structure to parse for run-time parameters.
  type(SAL_CS),                   intent(inout) :: CS !< Self-attraction and loading control structure
  type(thermo_var_ptrs),          intent(in)    :: tv !< A structure pointing to various thermodynamic variables.
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                                  intent(in)    :: h  !< Layer thicknesses [H ~> m or kg m-2]
  type(MOM_restart_CS), optional, intent(in)    :: restart_CS !< MOM restart control structure
  ! Local variables
                                        ! [R Z ~> kg m-2]

end subroutine SAL_init
module subroutine SAL_end(CS)
  type(SAL_CS), intent(inout) :: CS !< The control structure returned by a previous call
                                    !! to SAL_init; it is deallocated here.

end subroutine SAL_end
  end interface

end module MOM_self_attr_load
