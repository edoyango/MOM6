! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> Implements the thermodynamic aspects of ocean / ice-shelf interactions,
!!  along with a crude placeholder for a later implementation of full
!!  ice shelf dynamics, all using the MOM framework and coding style.
module MOM_ice_shelf_state

use MOM_cpu_clock, only : cpu_clock_id, cpu_clock_begin, cpu_clock_end
use MOM_cpu_clock, only : CLOCK_COMPONENT, CLOCK_ROUTINE
use MOM_dyn_horgrid, only : dyn_horgrid_type, create_dyn_horgrid, destroy_dyn_horgrid
use MOM_error_handler, only : MOM_error, MOM_mesg, FATAL, WARNING, is_root_pe
use MOM_file_parser, only : read_param, get_param, log_param, log_version, param_file_type
use MOM_grid, only : MOM_grid_init, ocean_grid_type
use MOM_get_input, only : directories, Get_MOM_input
use MOM_coms, only : reproducing_sum
use MOM_checksums, only : hchksum, qchksum, chksum, uchksum, vchksum, uvchksum

implicit none ; private

public ice_shelf_state_end, ice_shelf_state_init

!> Structure that describes the ice shelf state
type, public :: ice_shelf_state
  real, pointer, dimension(:,:) :: &
    mass_shelf => NULL(), &    !< The mass per unit area of the ice shelf or sheet [R Z ~> kg m-2].
    area_shelf_h => NULL(), &  !< The area per cell covered by the ice shelf [L2 ~> m2].
    melt_mask => NULL(), &     !< Mask is > 0 where melting is allowed [nondim]
    h_shelf => NULL(), &       !< the thickness of the shelf [Z ~> m], redundant with mass but may
                               !! make the code more readable
    dhdt_shelf => NULL(), &       !< the change in thickness of the shelf over time [Z T-1 ~> m s-1]
    hmask => NULL(),&          !< Mask used to indicate ice-covered or partiall-covered cells
                               !! 1: fully covered, solve for velocity here (for now all
                               !!   ice-covered cells are treated the same, this may change)
                               !! 2: partially covered, do not solve for velocity
                               !! 0: no ice in cell.
                               !! 3: bdry condition on thickness set
                               !! -2 : default (out of computational boundary)
                               !! NOTE: hmask will change over time and NEEDS TO BE MAINTAINED
                               !!   otherwise the wrong nodes will be included in velocity calcs.

    tflux_ocn => NULL(), &     !< The downward sensible ocean heat flux at the
                               !! ocean-ice interface [Q R Z T-1 ~> W m-2].
    salt_flux => NULL(), &     !< The downward salt flux at the ocean-ice
                               !! interface [kgSalt kgWater-1 R Z T-1 ~> kgSalt m-2 s-1].
    water_flux => NULL(), &    !< The net downward liquid water flux at the
                               !! ocean-ice interface [R Z T-1 ~> kg m-2 s-1].
    tflux_shelf => NULL(), &   !< The downward diffusive heat flux in the ice
                               !! shelf at the ice-ocean interface [Q R Z T-1 ~> W m-2].
    tfreeze => NULL(), &       !< The freezing point potential temperature
                               !! at the ice-ocean interface [C ~> degC].
    frazil => NULL(), &        !< Accumulated heating [J m-2] from frazil formation in the ocean
                               !! under ice-shelf cells
    !only active when calve_ice_shelf_bergs=true:
    calving => NULL(), &       !< The mass flux per unit area of the ice shelf to convert to
                               !! bergs [R Z T-1 ~> kg m-2 s-1].
    calving_hflx => NULL()     !< Calving heat flux [Q R Z T-1 ~> W m-2].
end type ice_shelf_state


  interface
module subroutine ice_shelf_state_init(ISS, G)
  type(ice_shelf_state), pointer    :: ISS !< A pointer to the ice shelf state structure
  type(ocean_grid_type), intent(in) :: G   !< The grid structure used by the ice shelf.

end subroutine ice_shelf_state_init
module subroutine ice_shelf_state_end(ISS)
  type(ice_shelf_state), pointer :: ISS !< A pointer to the ice shelf state structure

end subroutine ice_shelf_state_end
  end interface

end module MOM_ice_shelf_state
