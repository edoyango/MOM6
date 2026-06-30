! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> Provides checksumming functions for debugging
!!
!! This module contains subroutines that perform various error checking and
!! debugging functions for MOM6.  This routine is similar to it counterpart in
!! the SIS2 code, except for the use of the ocean_grid_type and by keeping them
!! separate we retain the ability to set up MOM6 and SIS2 debugging separately.
module MOM_debugging

use MOM_checksums,     only : hchksum, Bchksum, qchksum, uvchksum, hchksum_pair
use MOM_checksums,     only : is_NaN, chksum, MOM_checksums_init
use MOM_coms,          only : PE_here, root_PE, num_PEs
use MOM_coms,          only : min_across_PEs, max_across_PEs, reproducing_sum
use MOM_domains,       only : pass_vector, pass_var, pe_here
use MOM_domains,       only : BGRID_NE, AGRID, To_All, Scalar_Pair
use MOM_error_handler, only : MOM_error, FATAL, WARNING, is_root_pe
use MOM_file_parser,   only : log_version, param_file_type, get_param
use MOM_grid,          only : ocean_grid_type
use MOM_hor_index,     only : hor_index_type
use MOM_io,            only : stdout
use MOM_unit_scaling,  only : unit_scale_type

implicit none ; private

public :: check_redundant_C, check_redundant_B, check_redundant_T, check_redundant
public :: vec_chksum, vec_chksum_C, vec_chksum_B, vec_chksum_A
public :: MOM_debugging_init, totalStuff, totalTandS
public :: check_column_integral, check_column_integrals
public :: query_debugging_checks

! These interfaces come from MOM_checksums.
public :: hchksum, Bchksum, qchksum, is_NaN, chksum, uvchksum, hchksum_pair

!> Check for consistency between the duplicated points of a C-grid vector
interface check_redundant
  module procedure check_redundant_vC3d, check_redundant_vC2d
end interface check_redundant
!> Check for consistency between the duplicated points of a C-grid vector
interface check_redundant_C
  module procedure check_redundant_vC3d, check_redundant_vC2d
end interface check_redundant_C
!> Check for consistency between the duplicated points of a B-grid vector or scalar
interface check_redundant_B
  module procedure check_redundant_vB3d, check_redundant_vB2d
  module procedure check_redundant_sB3d, check_redundant_sB2d
end interface check_redundant_B
!> Check for consistency between the duplicated points of an A-grid vector or scalar
interface check_redundant_T
  module procedure check_redundant_sT3d, check_redundant_sT2d
  module procedure check_redundant_vT3d, check_redundant_vT2d
end interface check_redundant_T

!> Do checksums on the components of a C-grid vector
interface vec_chksum
  module procedure chksum_vec_C3d, chksum_vec_C2d
end interface vec_chksum
!> Do checksums on the components of a C-grid vector
interface vec_chksum_C
  module procedure chksum_vec_C3d, chksum_vec_C2d
end interface vec_chksum_C
!> Do checksums on the components of a B-grid vector
interface vec_chksum_B
  module procedure chksum_vec_B3d, chksum_vec_B2d
end interface vec_chksum_B
!> Do checksums on the components of an A-grid vector
interface vec_chksum_A
  module procedure chksum_vec_A3d, chksum_vec_A2d
end interface vec_chksum_A

! Note: these parameters are module data but ONLY used when debugging and
!       so can violate the thread-safe requirement of no module/global data.
integer :: max_redundant_prints = 100 !< Maximum number of times to write redundant messages
integer :: redundant_prints(3) = 0 !< Counters for controlling redundant printing
logical :: debug = .false. !< Write out verbose debugging data
logical :: debug_chksums = .true. !< Perform checksums on arrays
logical :: debug_redundant = .true. !< Check redundant values on PE boundaries


  interface
module subroutine MOM_debugging_init(param_file)
  type(param_file_type),   intent(in)    :: param_file !< A structure to parse for run-time parameters
  ! This include declares and sets the variable "version".

end subroutine MOM_debugging_init
module subroutine query_debugging_checks(do_debug, do_chksums, do_redundant)
  logical, optional, intent(out) :: do_debug     !< True if verbose debugging is to be output
  logical, optional, intent(out) :: do_chksums   !< True if checksums are to be output
  logical, optional, intent(out) :: do_redundant !< True if redundant points are to be checked

end subroutine query_debugging_checks
module subroutine check_redundant_vC3d(mesg, u_comp, v_comp, G, is, ie, js, je, &
                                direction, unscale)
  character(len=*),                    intent(in)    :: mesg   !< An identifying message
  type(ocean_grid_type),               intent(inout) :: G      !< The ocean's grid structure
  real, dimension(G%IsdB:,G%jsd:,:),   intent(in)    :: u_comp !< The u-component of the vector to be
                                                               !! checked for consistency in arbitrary,
                                                               !! possibly rescaled units [A ~> a]
  real, dimension(G%isd:,G%JsdB:,:),   intent(in)    :: v_comp !< The u-component of the vector to be
                                                               !! checked for consistency in arbitrary,
                                                               !! possibly rescaled units [A ~> a]
  integer,                   optional, intent(in)    :: is     !< The starting i-index to check
  integer,                   optional, intent(in)    :: ie     !< The ending i-index to check
  integer,                   optional, intent(in)    :: js     !< The starting j-index to check
  integer,                   optional, intent(in)    :: je     !< The ending j-index to check
  integer,                   optional, intent(in)    :: direction !< the direction flag to be
                                                               !! passed to pass_vector
  real,                      optional, intent(in)    :: unscale !< A factor that undoes the scaling for the
                                                               !! arrays to give consistent output [a A-1 ~> 1]

  ! Local variables

end subroutine check_redundant_vC3d
module subroutine check_redundant_vC2d(mesg, u_comp, v_comp, G, is, ie, js, je, &
                                direction, unscale)
  character(len=*),                intent(in)    :: mesg   !< An identifying message
  type(ocean_grid_type),           intent(inout) :: G      !< The ocean's grid structure
  real, dimension(G%IsdB:,G%jsd:), intent(in)    :: u_comp !< The u-component of the vector to be
                                                           !! checked for consistency in arbitrary,
                                                           !! possibly rescaled units [A ~> a]
  real, dimension(G%isd:,G%JsdB:), intent(in)    :: v_comp !< The u-component of the vector to be
                                                           !! checked for consistency in arbitrary,
                                                           !! possibly rescaled units [A ~> a]
  integer,               optional, intent(in)    :: is     !< The starting i-index to check
  integer,               optional, intent(in)    :: ie     !< The ending i-index to check
  integer,               optional, intent(in)    :: js     !< The starting j-index to check
  integer,               optional, intent(in)    :: je     !< The ending j-index to check
  integer,               optional, intent(in)    :: direction !< the direction flag to be
                                                           !! passed to pass_vector
  real,                  optional, intent(in)    :: unscale !< A factor that undoes the scaling for the
                                                           !! arrays to give consistent output [a A-1 ~> 1]
  ! Local variables
  ! In the following comments, [A] is used to indicate the arbitrary, possibly rescaled units
  ! of the input vector while [a] indicates the unscaled (e.g., mks) units to used for output.

end subroutine check_redundant_vC2d
module subroutine check_redundant_sB3d(mesg, array, G, is, ie, js, je, unscale)
  character(len=*),                     intent(in)    :: mesg  !< An identifying message
  type(ocean_grid_type),                intent(inout) :: G     !< The ocean's grid structure
  real, dimension(G%IsdB:,G%JsdB:,:),   intent(in)    :: array !< The array to be checked for consistency in
                                                               !! arbitrary, possibly rescaled units [A ~> a]
  integer,                    optional, intent(in)    :: is    !< The starting i-index to check
  integer,                    optional, intent(in)    :: ie    !< The ending i-index to check
  integer,                    optional, intent(in)    :: js    !< The starting j-index to check
  integer,                    optional, intent(in)    :: je    !< The ending j-index to check
  real,                       optional, intent(in)    :: unscale !< A factor that undoes the scaling for the
                                                               !! arrays to give consistent output [a A-1 ~> 1]

  ! Local variables

end subroutine check_redundant_sB3d
module subroutine check_redundant_sB2d(mesg, array, G, is, ie, js, je, unscale)
  character(len=*),                 intent(in)    :: mesg  !< An identifying message
  type(ocean_grid_type),            intent(inout) :: G     !< The ocean's grid structure
  real, dimension(G%IsdB:,G%JsdB:), intent(in)    :: array !< The array to be checked for consistency in
                                                           !! arbitrary, possibly rescaled units [A ~> a]
  integer,                optional, intent(in)    :: is    !< The starting i-index to check
  integer,                optional, intent(in)    :: ie    !< The ending i-index to check
  integer,                optional, intent(in)    :: js    !< The starting j-index to check
  integer,                optional, intent(in)    :: je    !< The ending j-index to check
  real,                   optional, intent(in)    :: unscale !< A factor that undoes the scaling for the
                                                           !! arrays to give consistent output [a A-1 ~> 1]
  ! Local variables
  ! In the following comments, [A] is used to indicate the arbitrary, possibly rescaled units
  ! of the input array while [a] indicates the unscaled (e.g., mks) units to used for output.

end subroutine check_redundant_sB2d
module subroutine check_redundant_vB3d(mesg, u_comp, v_comp, G, is, ie, js, je, &
                                direction, unscale)
  character(len=*),                    intent(in)    :: mesg   !< An identifying message
  type(ocean_grid_type),               intent(inout) :: G      !< The ocean's grid structure
  real, dimension(G%IsdB:,G%JsdB:,:),  intent(in)    :: u_comp !< The u-component of the vector to be
                                                               !! checked for consistency in arbitrary,
                                                               !! possibly rescaled units [A ~> a]
  real, dimension(G%IsdB:,G%JsdB:,:),  intent(in)    :: v_comp !< The v-component of the vector to be
                                                               !! checked for consistency in arbitrary,
                                                               !! possibly rescaled units [A ~> a]
  integer,                   optional, intent(in)    :: is     !< The starting i-index to check
  integer,                   optional, intent(in)    :: ie     !< The ending i-index to check
  integer,                   optional, intent(in)    :: js     !< The starting j-index to check
  integer,                   optional, intent(in)    :: je     !< The ending j-index to check
  integer,                   optional, intent(in)    :: direction !< the direction flag to be
                                                               !! passed to pass_vector
  real,                      optional, intent(in)    :: unscale !< A factor that undoes the scaling for the
                                                               !! arrays to give consistent output [a A-1 ~> 1]
  ! Local variables

end subroutine check_redundant_vB3d
module subroutine check_redundant_vB2d(mesg, u_comp, v_comp, G, is, ie, js, je, &
                                direction, unscale)
  character(len=*),                 intent(in)    :: mesg   !< An identifying message
  type(ocean_grid_type),            intent(inout) :: G      !< The ocean's grid structure
  real, dimension(G%IsdB:,G%JsdB:), intent(in)    :: u_comp !< The u-component of the vector to be
                                                            !! checked for consistency in arbitrary,
                                                            !! possibly rescaled units [A ~> a]
  real, dimension(G%IsdB:,G%JsdB:), intent(in)    :: v_comp !< The v-component of the vector to be
                                                            !! checked for consistency in arbitrary,
                                                            !! possibly rescaled units [A ~> a]
  integer,                optional, intent(in)    :: is     !< The starting i-index to check
  integer,                optional, intent(in)    :: ie     !< The ending i-index to check
  integer,                optional, intent(in)    :: js     !< The starting j-index to check
  integer,                optional, intent(in)    :: je     !< The ending j-index to check
  integer,                optional, intent(in)    :: direction !< the direction flag to be
                                                            !! passed to pass_vector
  real,                   optional, intent(in)    :: unscale !< A factor that undoes the scaling for the
                                                            !! arrays to give consistent output [a A-1 ~> 1]
  ! Local variables
  ! In the following comments, [A] is used to indicate the arbitrary, possibly rescaled units
  ! of the input vector while [a] indicates the unscaled (e.g., mks) units to used for output.

end subroutine check_redundant_vB2d
module subroutine check_redundant_sT3d(mesg, array, G, is, ie, js, je, unscale)
  character(len=*),                     intent(in)    :: mesg  !< An identifying message
  type(ocean_grid_type),                intent(inout) :: G     !< The ocean's grid structure
  real, dimension(G%isd:,G%jsd:,:),     intent(in)    :: array !< The array to be checked for consistency in
                                                               !! arbitrary, possibly rescaled units [A ~> a]
  integer,                    optional, intent(in)    :: is    !< The starting i-index to check
  integer,                    optional, intent(in)    :: ie    !< The ending i-index to check
  integer,                    optional, intent(in)    :: js    !< The starting j-index to check
  integer,                    optional, intent(in)    :: je    !< The ending j-index to check
  real,                       optional, intent(in)    :: unscale !< A factor that undoes the scaling for the
                                                               !! arrays to give consistent output [a A-1 ~> 1]
  ! Local variables

end subroutine check_redundant_sT3d
module subroutine check_redundant_sT2d(mesg, array, G, is, ie, js, je, unscale)
  character(len=*),                 intent(in)    :: mesg  !< An identifying message
  type(ocean_grid_type),            intent(inout) :: G     !< The ocean's grid structure
  real, dimension(G%isd:,G%jsd:),   intent(in)    :: array !< The array to be checked for consistency in
                                                           !! arbitrary, possibly rescaled units [A ~> a]
  integer,                optional, intent(in)    :: is    !< The starting i-index to check
  integer,                optional, intent(in)    :: ie    !< The ending i-index to check
  integer,                optional, intent(in)    :: js    !< The starting j-index to check
  integer,                optional, intent(in)    :: je    !< The ending j-index to check
  real,                   optional, intent(in)    :: unscale !< A factor that undoes the scaling for the
                                                           !! arrays to give consistent output [a A-1 ~> 1]
  ! Local variables
  ! In the following comments, [A] is used to indicate the arbitrary, possibly rescaled units
  ! of the input array while [a] indicates the unscaled (e.g., mks) units to used for output.

end subroutine check_redundant_sT2d
module subroutine check_redundant_vT3d(mesg, u_comp, v_comp, G, is, ie, js, je, &
                               direction, unscale)
  character(len=*),                    intent(in)    :: mesg   !< An identifying message
  type(ocean_grid_type),               intent(inout) :: G      !< The ocean's grid structure
  real, dimension(G%isd:,G%jsd:,:),    intent(in)    :: u_comp !< The u-component of the vector to be
                                                               !! checked for consistency in arbitrary,
                                                               !! possibly rescaled units [A ~> a]
  real, dimension(G%isd:,G%jsd:,:),    intent(in)    :: v_comp !< The v-component of the vector to be
                                                               !! checked for consistency in arbitrary,
                                                               !! possibly rescaled units [A ~> a]
  integer,                   optional, intent(in)    :: is     !< The starting i-index to check
  integer,                   optional, intent(in)    :: ie     !< The ending i-index to check
  integer,                   optional, intent(in)    :: js     !< The starting j-index to check
  integer,                   optional, intent(in)    :: je     !< The ending j-index to check
  integer,                   optional, intent(in)    :: direction !< the direction flag to be
                                                           !! passed to pass_vector
  real,                      optional, intent(in)    :: unscale !< A factor that undoes the scaling for the
                                                           !! arrays to give consistent output [a A-1 ~> 1]
  ! Local variables

end subroutine check_redundant_vT3d
module subroutine check_redundant_vT2d(mesg, u_comp, v_comp, G, is, ie, js, je, &
                               direction, unscale)
  character(len=*),                intent(in)    :: mesg   !< An identifying message
  type(ocean_grid_type),           intent(inout) :: G      !< The ocean's grid structure
  real, dimension(G%isd:,G%jsd:),  intent(in)    :: u_comp !< The u-component of the vector to be
                                                           !! checked for consistency in arbitrary,
                                                           !! possibly rescaled units [A ~> a]
  real, dimension(G%isd:,G%jsd:),  intent(in)    :: v_comp !< The v-component of the vector to be
                                                           !! checked for consistency in arbitrary,
                                                           !! possibly rescaled units [A ~> a]
  integer,               optional, intent(in)    :: is     !< The starting i-index to check
  integer,               optional, intent(in)    :: ie     !< The ending i-index to check
  integer,               optional, intent(in)    :: js     !< The starting j-index to check
  integer,               optional, intent(in)    :: je     !< The ending j-index to check
  integer,               optional, intent(in)    :: direction !< the direction flag to be
                                                           !! passed to pass_vector
  real,                  optional, intent(in)    :: unscale !< A factor that undoes the scaling for the
                                                           !! arrays to give consistent output [a A-1 ~> 1]
  ! Local variables
  ! In the following comments, [A] is used to indicate the arbitrary, possibly rescaled units
  ! of the input vector while [a] indicates the unscaled (e.g., mks) units to used for output.

end subroutine check_redundant_vT2d
module subroutine chksum_vec_C3d(mesg, u_comp, v_comp, G, halos, scalars, unscale)
  character(len=*),                  intent(in)    :: mesg   !< An identifying message
  type(ocean_grid_type),             intent(inout) :: G      !< The ocean's grid structure
  real, dimension(G%IsdB:,G%jsd:,:), intent(in)    :: u_comp !< The u-component of the vector to be
                                                             !! checked for consistency in arbitrary,
                                                             !! possibly rescaled units [A ~> a]
  real, dimension(G%isd:,G%JsdB:,:), intent(in)    :: v_comp !< The v-component of the vector to be
                                                             !! checked for consistency in arbitrary,
                                                             !! possibly rescaled units [A ~> a]
  integer,                 optional, intent(in)    :: halos  !< The width of halos to check (default 0)
  logical,                 optional, intent(in)    :: scalars !< If true this is a pair of
                                                             !! scalars that are being checked.
  real,                    optional, intent(in)    :: unscale !< A factor that undoes the scaling for the
                                                             !! arrays to give consistent output [a A-1 ~> 1]
  ! Local variables
end subroutine chksum_vec_C3d
module subroutine chksum_vec_C2d(mesg, u_comp, v_comp, G, halos, scalars, unscale)
  character(len=*),                intent(in)    :: mesg   !< An identifying message
  type(ocean_grid_type),           intent(inout) :: G      !< The ocean's grid structure
  real, dimension(G%IsdB:,G%jsd:), intent(in)    :: u_comp !< The u-component of the vector to be
                                                           !! checked for consistency in arbitrary,
                                                           !! possibly rescaled units [A ~> a]
  real, dimension(G%isd:,G%JsdB:), intent(in)    :: v_comp !< The v-component of the vector to be
                                                           !! checked for consistency in arbitrary,
                                                           !! possibly rescaled units [A ~> a]
  integer,               optional, intent(in)    :: halos  !< The width of halos to check (default 0)
  logical,               optional, intent(in)    :: scalars !< If true this is a pair of
                                                           !! scalars that are being checked.
  real,                  optional, intent(in)    :: unscale !< A factor that undoes the scaling for the
                                                           !! arrays to give consistent output [a A-1 ~> 1]
  ! Local variables
end subroutine chksum_vec_C2d
module subroutine chksum_vec_B3d(mesg, u_comp, v_comp, G, halos, scalars, unscale)
  character(len=*),                   intent(in)    :: mesg   !< An identifying message
  type(ocean_grid_type),              intent(inout) :: G      !< The ocean's grid structure
  real, dimension(G%IsdB:,G%JsdB:,:), intent(in)    :: u_comp !< The u-component of the vector to be
                                                              !! checked for consistency in arbitrary,
                                                              !! possibly rescaled units [A ~> a]
  real, dimension(G%IsdB:,G%JsdB:,:), intent(in)    :: v_comp !< The v-component of the vector to be
                                                              !! checked for consistency in arbitrary,
                                                              !! possibly rescaled units [A ~> a]
  integer,                  optional, intent(in)    :: halos  !< The width of halos to check (default 0)
  logical,                  optional, intent(in)    :: scalars !< If true this is a pair of
                                                              !! scalars that are being checked.
  real,                     optional, intent(in)    :: unscale !< A factor that undoes the scaling for the
                                                              !! arrays to give consistent output [a A-1 ~> 1]
  ! Local variables
end subroutine chksum_vec_B3d
module subroutine chksum_vec_B2d(mesg, u_comp, v_comp, G, halos, scalars, symmetric, unscale)
  character(len=*),                 intent(in)    :: mesg   !< An identifying message
  type(ocean_grid_type),            intent(inout) :: G      !< The ocean's grid structure
  real, dimension(G%IsdB:,G%JsdB:), intent(in)    :: u_comp !< The u-component of the vector to be
                                                            !! checked for consistency in arbitrary,
                                                            !! possibly rescaled units [A ~> a]
  real, dimension(G%IsdB:,G%JsdB:), intent(in)    :: v_comp !< The v-component of the vector to be
                                                            !! checked for consistency in arbitrary,
                                                            !! possibly rescaled units [A ~> a]
  integer,                optional, intent(in)    :: halos  !< The width of halos to check (default 0)
  logical,                optional, intent(in)    :: scalars !< If true this is a pair of
                                                            !! scalars that are being checked.
  logical,                optional, intent(in)    :: symmetric !< If true, do the checksums on the
                                                            !! full symmetric computational domain.
  real,                   optional, intent(in)    :: unscale !< A factor that undoes the scaling for the
                                                            !! arrays to give consistent output [a A-1 ~> 1]
  ! Local variables
end subroutine chksum_vec_B2d
module subroutine chksum_vec_A3d(mesg, u_comp, v_comp, G, halos, scalars, unscale)
  character(len=*),                 intent(in)    :: mesg   !< An identifying message
  type(ocean_grid_type),            intent(inout) :: G      !< The ocean's grid structure
  real, dimension(G%isd:,G%jsd:,:), intent(in)    :: u_comp !< The u-component of the vector to be
                                                            !! checked for consistency in arbitrary,
                                                            !! possibly rescaled units [A ~> a]
  real, dimension(G%isd:,G%jsd:,:), intent(in)    :: v_comp !< The v-component of the vector to be
                                                            !! checked for consistency in arbitrary,
                                                            !! possibly rescaled units [A ~> a]
  integer,                optional, intent(in)    :: halos  !< The width of halos to check (default 0)
  logical,                optional, intent(in)    :: scalars !< If true this is a pair of
                                                            !! scalars that are being checked.
  real,                   optional, intent(in)    :: unscale !< A factor that undoes the scaling for the
                                                            !! arrays to give consistent output [a A-1 ~> 1]
  ! Local variables
end subroutine chksum_vec_A3d
module subroutine chksum_vec_A2d(mesg, u_comp, v_comp, G, halos, scalars, unscale)
  character(len=*),               intent(in)    :: mesg   !< An identifying message
  type(ocean_grid_type),          intent(inout) :: G      !< The ocean's grid structure
  real, dimension(G%isd:,G%jsd:), intent(in)    :: u_comp !< The u-component of the vector to be
                                                          !! checked for consistency in arbitrary,
                                                          !! possibly rescaled units [A ~> a]
  real, dimension(G%isd:,G%jsd:), intent(in)    :: v_comp !< The v-component of the vector to be
                                                          !! checked for consistency in arbitrary,
                                                          !! possibly rescaled units [A ~> a]
  integer,              optional, intent(in)    :: halos  !< The width of halos to check (default 0)
  logical,              optional, intent(in)    :: scalars !< If true this is a pair of
                                                          !! scalars that are being checked.
  real,                 optional, intent(in)    :: unscale !< A factor that undoes the scaling for the
                                                          !! arrays to give consistent output [a A-1 ~> 1]
  ! Local variables
end subroutine chksum_vec_A2d
module function totalStuff(HI, hThick, areaT, stuff, unscale)
  type(hor_index_type),               intent(in) :: HI     !< A horizontal index type
  real, dimension(HI%isd:,HI%jsd:,:), intent(in) :: hThick !< The array of thicknesses to use as weights
                                                           !! [H ~> m or kg m-2] or [m] or [kg m-2]
  real, dimension(HI%isd:,HI%jsd:),   intent(in) :: areaT  !< The array of cell areas [L2 ~> m2] or [m2]
  real, dimension(HI%isd:,HI%jsd:,:), intent(in) :: stuff  !< The array of stuff to be summed in arbitrary
                                                           !! units [A ~> a] or [a]
  real,                     optional, intent(in) :: unscale !< A factor that is used to undo scaling of the array
                                                           !! and the cell mass or volume before it is summed in
                                                           !! [a m3 A-1 H-1 L-2 ~> 1] or [a kg A-1 H-1 L-2 ~> 1]
  real                                           :: totalStuff !< the globally integrated amount of stuff
                                                           !! [A H L2 ~> a m3 or a kg] or [a m3]
  ! Local variables
                                                      ! cell [A H L2 ~> a m3 or a kg] or [a m3]

end function totalStuff
module subroutine totalTandS(HI, hThick, areaT, temperature, salinity, mesg, US, H_to_mks)
  type(hor_index_type),               intent(in) :: HI     !< A horizontal index type
  real, dimension(HI%isd:,HI%jsd:,:), intent(in) :: hThick !< The array of thicknesses to use as weights
                                                           !! [H ~> m or kg m-2] or [m] or [kg m-2]
  real, dimension(HI%isd:,HI%jsd:),   intent(in) :: areaT  !< The array of cell areas [L2 ~> m2] or [m2]
  real, dimension(HI%isd:,HI%jsd:,:), intent(in) :: temperature !< The temperature field to sum [C ~> degC] or [degC]
  real, dimension(HI%isd:,HI%jsd:,:), intent(in) :: salinity    !< The salinity field to sum [S ~> ppt] or [ppt]
  character(len=*),                   intent(in) :: mesg        !< An identifying message
  type(unit_scale_type),    optional, intent(in) :: US       !< A dimensional unit scaling type
  real,                     optional, intent(in) :: H_to_MKS !< A constant that translates thickness units to its
                                                             !! MKS units (m or kg m-2) based on whether the model is
                                                             !! Boussinesq [m H-1 ~> 1] or not [kg m-2 H-1 ~> 1]
  ! NOTE: This subroutine uses "save" data which is not thread safe and is purely for
  ! extreme debugging without a proper debugger.
                              ! call [H L2 ~> m3 or kg] or [m3] or [kg]
                              ! call [C H L2 ~> degC m3 or degC kg] or [degC m3] or [degC kg]
                              ! call [S H L2 ~> ppt m3 or ppt kg] or [ppt m3] or [ppt kg]
  ! Local variables
                       ! call [C H L2 ~> degC m3 or degC kg] or [degC m3] or [degC kg]
                       ! call [S H L2 ~> ppt m3 or ppt kg] or [ppt m3] or [ppt kg]
                       ! whether the model is Boussinesq [m H-1 ~> 1] or non-Boussinesq [kg m-2 H-1 ~> 1]
                       ! [degC kg C-1 H-1 L-2 ~> 1]
                       ! [ppt kg S-1 H-1 L-2 ~> 1]

end subroutine totalTandS
logical module function check_column_integral(nk, field, known_answer)
  integer,             intent(in) :: nk           !< Number of levels in column
  real, dimension(nk), intent(in) :: field        !< Field to be summed [arbitrary]
  real, optional,      intent(in) :: known_answer !< If present is the expected sum [arbitrary],
                                                  !! If missing, assumed zero
  ! Local variables

end function check_column_integral
logical module function check_column_integrals(nk_1, field_1, nk_2, field_2, missing_value)
  integer,               intent(in) :: nk_1           !< Number of levels in field 1
  integer,               intent(in) :: nk_2           !< Number of levels in field 2
  real, dimension(nk_1), intent(in) :: field_1        !< First field to be summed [arbitrary]
  real, dimension(nk_2), intent(in) :: field_2        !< Second field to be summed [arbitrary]
  real, optional,        intent(in) :: missing_value  !< If column contains missing values,
                                                      !! mask them from the sum [arbitrary]
  ! Local variables
                            ! from the sums [arbitrary]

  ! Assign missing value
end function check_column_integrals
  end interface

end module MOM_debugging
