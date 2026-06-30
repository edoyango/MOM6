! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> Interface to background mixing schemes, including the Bryan and Lewis (1979)
!! which is applied via CVMix.

module MOM_bkgnd_mixing

use MOM_debugging,       only : hchksum
use MOM_diag_mediator,   only : diag_ctrl, time_type, register_diag_field
use MOM_diag_mediator,   only : post_data
use MOM_error_handler,   only : MOM_error, FATAL, WARNING, NOTE
use MOM_file_parser,     only : get_param, log_param, log_version, param_file_type
use MOM_file_parser,     only : openParameterBlock, closeParameterBlock
use MOM_forcing_type,    only : forcing
use MOM_grid,            only : ocean_grid_type
use MOM_interface_heights, only : thickness_to_dz
use MOM_unit_scaling,    only : unit_scale_type
use MOM_verticalGrid,    only : verticalGrid_type
use MOM_variables,       only : thermo_var_ptrs,  vertvisc_type
use MOM_intrinsic_functions, only : invcosh
use CVMix_background,    only : CVMix_init_bkgnd, CVMix_coeffs_bkgnd

implicit none ; private

#include <MOM_memory.h>

public bkgnd_mixing_init
public bkgnd_mixing_end
public calculate_bkgnd_mixing

! A note on unit descriptions in comments: MOM6 uses units that can be rescaled for dimensional
! consistency testing. These are noted in comments with units like Z, H, L, and T, along with
! their mks counterparts with notation like "a velocity [Z T-1 ~> m s-1]".  If the units
! vary with the Boussinesq approximation, the Boussinesq variant is given first.

!> Control structure including parameters for this module.
type, public :: bkgnd_mixing_cs ; private

  ! Parameters
  real    :: Bryan_Lewis_c1         !< The vertical diffusivity values for  Bryan-Lewis profile
                                    !! at |z|=D [Z2 T-1 ~> m2 s-1]
  real    :: Bryan_Lewis_c2         !< The amplitude of variation in diffusivity for the
                                    !! Bryan-Lewis diffusivity profile [Z2 T-1 ~> m2 s-1]
  real    :: Bryan_Lewis_c3         !< The inverse length scale for transition region in the
                                    !! Bryan-Lewis diffusivity profile [Z-1 ~> m-1]
  real    :: Bryan_Lewis_c4         !< The depth where diffusivity is Bryan_Lewis_bl1 in the
                                    !! Bryan-Lewis profile [Z ~> m]
  real    :: bckgrnd_vdc1           !< Background diffusivity (Ledwell) when
                                    !! horiz_varying_background=.true. [H Z T-1 ~> m2 s-1 or kg m-1 s-1]
  real    :: bckgrnd_vdc_eq         !< Equatorial diffusivity (Gregg) when
                                    !! horiz_varying_background=.true. [H Z T-1 ~> m2 s-1 or kg m-1 s-1]
  real    :: bckgrnd_vdc_psim       !< Max. PSI induced diffusivity (MacKinnon) when
                                    !! horiz_varying_background=.true. [H Z T-1 ~> m2 s-1 or kg m-1 s-1]
  real    :: bckgrnd_vdc_Banda      !< Banda Sea diffusivity (Gordon) when
                                    !! horiz_varying_background=.true. [H Z T-1 ~> m2 s-1 or kg m-1 s-1]
  real    :: Kd_min                 !< minimum diapycnal diffusivity [H Z T-1 ~> m2 s-1 or kg m-1 s-1]
  real    :: Kd                     !< interior diapycnal diffusivity [H Z T-1 ~> m2 s-1 or kg m-1 s-1]
  real    :: omega                  !< The Earth's rotation rate [T-1 ~> s-1].
  real    :: N0_2Omega              !< ratio of the typical Buoyancy frequency to
                                    !! twice the Earth's rotation period, used with the
                                    !! Henyey scaling from the mixing [nondim]
  real    :: Henyey_max_lat         !< A latitude poleward of which the Henyey profile
                                    !! is returned to the minimum diffusivity [degrees_N]
  real    :: prandtl_bkgnd          !< Turbulent Prandtl number used to convert
                                    !! vertical background diffusivity into viscosity [nondim]
  real    :: Kd_tanh_lat_scale      !< A nondimensional scaling for the range of
                                    !! diffusivities with Kd_tanh_lat_fn [nondim]. Valid values
                                    !! are in the range of -2 to 2; 0.4 reproduces CM2M.
  real    :: Kd_tot_ml              !< The mixed layer diapycnal diffusivity [H Z T-1 ~> m2 s-1 or kg m-1 s-1]
                                    !! when no other physically based mixed layer turbulence
                                    !! parameterization is being used.
  real    :: Hmix                   !< mixed layer thickness [H ~> m or kg m-2] when no physically based
                                    !! ocean surface boundary layer parameterization is used.
  logical :: Kd_tanh_lat_fn         !< If true, use the tanh dependence of Kd_sfc on
                                    !! latitude, like GFDL CM2.1/CM2M.  There is no
                                    !! physical justification for this form, and it can
                                    !! not be used with Henyey_IGW_background.
  logical :: Bryan_Lewis_diffusivity!< If true, background vertical diffusivity
                                    !! uses Bryan-Lewis (1979) like tanh profile.
  logical :: horiz_varying_background !< If true, apply vertically uniform, latitude-dependent
                                    !! background diffusivity, as described in Danabasoglu et al., 2012
  logical :: Henyey_IGW_background  !< If true, use a simplified variant of the
             !! Henyey et al, JGR (1986) latitudinal scaling for the background diapycnal diffusivity,
             !! which gives a marked decrease in the diffusivity near the equator.  The simplification
             !! here is to assume that the in-situ stratification is the same as the reference stratificaiton.
  logical :: physical_OBL_scheme !< If true, a physically-based scheme is used to determine mixing in the
                   !! ocean's surface boundary layer, such as ePBL, KPP, or a refined bulk mixed layer scheme.
  logical :: debug !< If true, turn on debugging in this module
  ! Diagnostic handles and pointers
  type(diag_ctrl), pointer :: diag => NULL() !< A structure that regulates diagnostic output

  character(len=40)  :: bkgnd_scheme_str = "none" !< Background scheme identifier

end type bkgnd_mixing_cs

character(len=40)  :: mdl = "MOM_bkgnd_mixing" !< This module's name.


  interface
module subroutine bkgnd_mixing_init(Time, G, GV, US, param_file, diag, CS, physical_OBL_scheme)

  type(time_type),         intent(in)    :: Time       !< The current time.
  type(ocean_grid_type),   intent(in)    :: G          !< Grid structure.
  type(verticalGrid_type), intent(in)    :: GV         !< Vertical grid structure.
  type(unit_scale_type),   intent(in)    :: US         !< A dimensional unit scaling type
  type(param_file_type),   intent(in)    :: param_file !< Run-time parameter file handle
  type(diag_ctrl), target, intent(inout) :: diag       !< Diagnostics control structure.
  type(bkgnd_mixing_cs),   pointer       :: CS         !< This module's control structure.
  logical,                 intent(in)    :: physical_OBL_scheme !< If true, a physically based
                                                       !! parameterization (like KPP or ePBL or a bulk mixed
                                                       !! layer) is used outside of set_diffusivity to
                                                       !! specify the mixing that occurs in the ocean's
                                                       !! surface boundary layer.

  ! Local variables
                                ! number unless it is provided as a parameter
                                ! in setting the default for other diffusivities.

  ! This include declares and sets the variable "version".

end subroutine bkgnd_mixing_init
module subroutine calculate_bkgnd_mixing(h, tv, N2_lay, Kd_lay, Kd_int, Kv_bkgnd, j, G, GV, US, CS)

  type(ocean_grid_type),                     intent(in)    :: G   !< Grid structure.
  type(verticalGrid_type),                   intent(in)    :: GV  !< Vertical grid structure.
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in)    :: h   !< Layer thickness [H ~> m or kg m-2].
  type(thermo_var_ptrs),                     intent(in)    :: tv  !< Thermodynamics structure.
  real, dimension(SZI_(G),SZK_(GV)),         intent(in)    :: N2_lay !< squared buoyancy frequency associated
                                                                  !! with layers [T-2 ~> s-2]
  real, dimension(SZI_(G),SZK_(GV)),         intent(out)   :: Kd_lay !< The background diapycnal diffusivity of each
                                                                  !! layer [H Z T-1 ~> m2 s-1 or kg m-1 s-1]
  real, dimension(SZI_(G),SZK_(GV)+1),       intent(out)   :: Kd_int !< The background diapycnal diffusivity of each
                                                                  !! interface [H Z T-1 ~> m2 s-1 or kg m-1 s-1]
  real, dimension(SZI_(G),SZK_(GV)+1),       intent(out)   :: Kv_bkgnd !< The background vertical viscosity at
                                                                  !! each interface [H Z T-1 ~> m2 s-1 or Pa s]
  integer,                                   intent(in)    :: j   !< Meridional grid index
  type(unit_scale_type),                     intent(in)    :: US  !< A dimensional unit scaling type
  type(bkgnd_mixing_cs),                     pointer       :: CS  !< The control structure returned by
                                                                  !! a previous call to bkgnd_mixing_init.

  ! local variables

end subroutine calculate_bkgnd_mixing
logical module function CVMix_bkgnd_is_used(param_file)
  type(param_file_type), intent(in) :: param_file !< A structure to parse for run-time parameters
end function CVMix_bkgnd_is_used
module subroutine check_bkgnd_scheme(CS, str)
  type(bkgnd_mixing_cs), pointer :: CS  !< Control structure
  character(len=*), intent(in)   :: str !< Background scheme identifier deducted from MOM_input
                                        !! parameters

end subroutine check_bkgnd_scheme
module subroutine bkgnd_mixing_end(CS)
  type(bkgnd_mixing_cs), pointer :: CS !< Control structure for this module that
                                       !! will be deallocated in this subroutine

end subroutine bkgnd_mixing_end
  end interface

end module MOM_bkgnd_mixing
