! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> \brief Parameterization of mixed layer restratification by unresolved mixed-layer eddies.
module MOM_mixed_layer_restrat

use MOM_debugging,     only : hchksum
use MOM_diag_mediator, only : post_data, query_averaging_enabled, diag_ctrl
use MOM_diag_mediator, only : register_diag_field, safe_alloc_ptr, time_type
use MOM_diag_mediator, only : diag_update_remap_grids
use MOM_domains,       only : pass_var, To_West, To_South, Omit_Corners
use MOM_error_handler, only : MOM_error, FATAL, WARNING
use MOM_file_parser,   only : get_param, log_param, log_version, param_file_type
use MOM_file_parser,   only : openParameterBlock, closeParameterBlock
use MOM_forcing_type,  only : mech_forcing, find_ustar
use MOM_grid,          only : ocean_grid_type
use MOM_hor_index,     only : hor_index_type
use MOM_interpolate,   only : init_external_field, time_interp_external, time_interp_external_init
use MOM_interpolate,   only : external_field
use MOM_intrinsic_functions, only : cuberoot
use MOM_io,            only : slasher, MOM_read_data
use MOM_lateral_mixing_coeffs, only : VarMix_CS
use MOM_restart,       only : register_restart_field, query_initialized, MOM_restart_CS
use MOM_unit_scaling,  only : unit_scale_type
use MOM_variables,     only : thermo_var_ptrs
use MOM_verticalGrid,  only : verticalGrid_type, get_thickness_units
use MOM_EOS,           only : calculate_density, calculate_spec_vol, EOS_domain

implicit none ; private

#include <MOM_memory.h>

public mixedlayer_restrat
public mixedlayer_restrat_init
public mixedlayer_restrat_register_restarts
public mixedlayer_restrat_unit_tests

! A note on unit descriptions in comments: MOM6 uses units that can be rescaled for dimensional
! consistency testing. These are noted in comments with units like Z, H, L, and T, along with
! their mks counterparts with notation like "a velocity [Z T-1 ~> m s-1]".  If the units
! vary with the Boussinesq approximation, the Boussinesq variant is given first.

!> Control structure for mom_mixed_layer_restrat
type, public :: mixedlayer_restrat_CS ; private
  logical :: initialized = .false. !< True if this control structure has been initialized.
  real    :: ml_restrat_coef       !< A non-dimensional factor by which the instability is enhanced
                                   !! over what would be predicted based on the resolved gradients
                                   !! [nondim].  This increases with grid spacing^2, up to something
                                   !! of order 500.
  real    :: ml_restrat_coef2      !< As for ml_restrat_coef but using the slow filtered MLD [nondim].
  real    :: front_length          !< If non-zero, is the frontal-length scale [L ~> m] used to calculate the
                                   !! upscaling of buoyancy gradients that is otherwise represented
                                   !! by the parameter FOX_KEMPER_ML_RESTRAT_COEF. If MLE_FRONT_LENGTH is
                                   !! non-zero, it is recommended to set FOX_KEMPER_ML_RESTRAT_COEF=1.0.
  logical :: fl_from_file          !< If true, read the MLE front-length scale from a netCDF file.
  logical :: MLE_use_PBL_MLD       !< If true, use the MLD provided by the PBL parameterization.
                                   !! if false, MLE will calculate a MLD based on a density difference
                                   !! based on the parameter MLE_DENSITY_DIFF.
  logical :: Bodner_detect_MLD        !< If true, detect the MLD based on given density difference criterion
                                   !! (MLE_DENSITY_DIFF) in the Bodner et al. parameterization.
  real    :: vonKar                !< The von Karman constant as used for mixed layer viscosity [nondim]
  real    :: MLE_MLD_decay_time    !< Time-scale to use in a running-mean when MLD is retreating [T ~> s].
  real    :: MLE_MLD_decay_time2   !< Time-scale to use in a running-mean when filtered MLD is retreating [T ~> s].
  real    :: MLE_density_diff      !< Density difference used in detecting mixed-layer depth [R ~> kg m-3].
  real    :: MLE_tail_dh           !< Fraction by which to extend the mixed-layer restratification
                                   !! depth used for a smoother stream function at the base of
                                   !! the mixed-layer [nondim].
  real    :: MLE_MLD_stretch       !< A scaling coefficient for stretching/shrinking the MLD used in
                                   !! the MLE scheme [nondim]. This simply multiplies MLD wherever used.

  ! The following parameters are used in the Bodner et al., 2023, parameterization
  logical :: use_Bodner = .false.  !< If true, use the Bodner et al., 2023, parameterization.
  real    :: Cr                    !< Efficiency coefficient from Bodner et al., 2023 [nondim]
  real    :: mstar                 !< The m* value used to estimate the turbulent vertical momentum flux [nondim]
  real    :: nstar                 !< The n* value used to estimate the turbulent vertical momentum flux [nondim]
  real    :: min_wstar2            !< The minimum lower bound to apply to the vertical momentum flux,
                                   !! w'u', in the Bodner et al., restratification parameterization
                                   !! [Z2 T-2 ~> m2 s-2].  This avoids a division-by-zero in the limit when u*
                                   !! and the buoyancy flux are zero.
  real    :: BLD_growing_Tfilt     !< The time-scale for a running-mean filter applied to the boundary layer
                                   !! depth (BLD) when the BLD is deeper than the running mean [T ~> s].
                                   !! A value of 0 instantaneously sets the running mean to the current value of BLD.
  real    :: BLD_decaying_Tfilt    !< The time-scale for a running-mean filter applied to the boundary layer
                                   !! depth (BLD) when the BLD is shallower than the running mean [T ~> s].
                                   !! A value of 0 instantaneously sets the running mean to the current value of BLD.
  real    :: MLD_decaying_Tfilt    !< The time-scale for a running-mean filter applied to the time-filtered
                                   !! MLD, when the latter is shallower than the running mean [T ~> s].
                                   !! A value of 0 instantaneously sets the running mean to the current value of MLD.
  real    :: MLD_growing_Tfilt     !< The time-scale for a running-mean filter applied to the time-filtered
                                   !! MLD, when the latter is deeper than the running mean [T ~> s].
                                   !! A value of 0 instantaneously sets the running mean to the current value of MLD.
  integer :: answer_date           !< The vintage of the order of arithmetic and expressions in the
                                   !! mixed layer restrat calculations.  Values below 20240201 recover
                                   !! the answers from the end of 2023, while higher values use the new
                                   !! cuberoot function in the Bodner code to avoid needing to undo
                                   !! dimensional rescaling.

  logical :: debug = .false.       !< If true, calculate checksums of fields for debugging.

  type(diag_ctrl), pointer :: diag !< A structure that is used to regulate the
                                   !! timing of diagnostic output.
  type(external_field) :: sbc_fl   !< A handle used in time interpolation of
                                   !! front-length scales read from a file.
  type(time_type), pointer :: Time => NULL() !< A pointer to the ocean model's clock.
  logical :: use_Stanley_ML        !< If true, use the Stanley parameterization of SGS T variance
  real    :: ustar_min             !< A minimum value of ustar in thickness units to avoid numerical
                                   !! problems [H T-1 ~> m s-1 or kg m-2 s-1]
  real    :: Kv_restrat            !< A viscosity that sets a floor on the momentum mixing rate
                                   !! during restratification, rescaled into thickness-based
                                   !! units [H2 T-1 ~> m2 s-1 or kg2 m-4 s-1]
  logical :: MLD_grid              !< If true, read a spacially varying field for MLD_decaying_Tfilt
  logical :: Cr_grid               !< If true, read a spacially varying field for Cr

  real, dimension(:,:), allocatable :: &
         MLD_filtered, &           !< Time-filtered MLD [H ~> m or kg m-2]
         MLD_filtered_slow, &      !< Slower time-filtered MLD [H ~> m or kg m-2]
         wpup_filtered, &          !< Time-filtered vertical momentum flux [H L T-2 ~> m2 s-2 or kg m-1 s-2]
         MLD_Tfilt_space, &        !< Spatially varying time scale for MLD filter [T ~> s]
         Cr_space                  !< Spatially varying Cr coefficient [nondim]

  !>@{
  !! Diagnostic identifier
  integer :: id_urestrat_time = -1
  integer :: id_vrestrat_time = -1
  integer :: id_uhml = -1
  integer :: id_vhml = -1
  integer :: id_MLD  = -1
  integer :: id_BLD  = -1
  integer :: id_Rml  = -1
  integer :: id_uDml = -1
  integer :: id_vDml = -1
  integer :: id_uml  = -1
  integer :: id_vml  = -1
  integer :: id_wpup = -1
  integer :: id_ustar = -1
  integer :: id_bflux = -1
  integer :: id_lfbod = -1
  integer :: id_mle_fl = -1
  !>@}

end type mixedlayer_restrat_CS

character(len=40)  :: mdl = "MOM_mixed_layer_restrat" !< This module's name.


  interface
module subroutine mixedlayer_restrat(h, uhtr, vhtr, tv, forces, dt, MLD, h_MLD, bflux, VarMix, G, GV, US, CS)
  type(ocean_grid_type),                      intent(inout) :: G      !< Ocean grid structure
  type(verticalGrid_type),                    intent(in)    :: GV     !< Ocean vertical grid structure
  type(unit_scale_type),                      intent(in)    :: US     !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),  intent(inout) :: h      !< Layer thickness [H ~> m or kg m-2]
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), intent(inout) :: uhtr   !< Accumulated zonal mass flux
                                                                      !! [H L2 ~> m3 or kg]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), intent(inout) :: vhtr   !< Accumulated meridional mass flux
                                                                      !! [H L2 ~> m3 or kg]
  type(thermo_var_ptrs),                      intent(in)    :: tv     !< Thermodynamic variables structure
  type(mech_forcing),                         intent(in)    :: forces !< A structure with the driving mechanical forces
  real,                                       intent(in)    :: dt     !< Time increment [T ~> s]
  real, dimension(:,:),                       pointer       :: MLD    !< Mixed layer depth provided by the
                                                                      !! planetary boundary layer scheme [Z ~> m]
  real, dimension(:,:),                       pointer       :: h_MLD  !< Mixed layer thickness provided
                                                                      !! by the planetary boundary layer
                                                                      !! scheme [H ~> m or kg m-2]
  real, dimension(:,:),                       pointer       :: bflux  !< Surface buoyancy flux provided by the
                                                                      !! PBL scheme [Z2 T-3 ~> m2 s-3]
  type(VarMix_CS),                            intent(in)    :: VarMix !< Variable mixing control structure
  type(mixedlayer_restrat_CS),                intent(inout) :: CS     !< Module control structure


end subroutine mixedlayer_restrat
module subroutine mixedlayer_restrat_OM4(h, uhtr, vhtr, tv, forces, dt, h_MLD, VarMix, G, GV, US, CS)
  ! Arguments
  type(ocean_grid_type),                      intent(inout) :: G      !< Ocean grid structure
  type(verticalGrid_type),                    intent(in)    :: GV     !< Ocean vertical grid structure
  type(unit_scale_type),                      intent(in)    :: US     !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),  intent(inout) :: h      !< Layer thickness [H ~> m or kg m-2]
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), intent(inout) :: uhtr   !< Accumulated zonal mass flux
                                                                      !!   [H L2 ~> m3 or kg]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), intent(inout) :: vhtr   !< Accumulated meridional mass flux
                                                                      !!   [H L2 ~> m3 or kg]
  type(thermo_var_ptrs),                      intent(in)    :: tv     !< Thermodynamic variables structure
  type(mech_forcing),                         intent(in)    :: forces !< A structure with the driving mechanical forces
  real,                                       intent(in)    :: dt     !< Time increment [T ~> s]
  real, dimension(:,:),                       pointer       :: h_MLD  !< Thickness of water within the
                                                                      !! mixed layer depth provided by
                                                                      !!  the PBL scheme [H ~> m or kg m-2]
  type(VarMix_CS),                            intent(in)    :: VarMix !< Variable mixing control structure
  type(mixedlayer_restrat_CS),                intent(inout) :: CS     !< Module control structure

  ! Local variables
                          ! sublayer of the mixed layer, divided by dt [H L2 T-1 ~> m3 s-1 or kg s-1].
                          ! in non-Boussinesq mode [H T-1 ~> m s-1 or kg m-2 s-1]
end subroutine mixedlayer_restrat_OM4
real module function mu(sigma, dh)
  real, intent(in) :: sigma !< Fractional position within mixed layer [nondim]
                            !! z=0 is surface, z=-1 is the bottom of the mixed layer
  real, intent(in) :: dh    !< Non-dimensional distance over which to extend stream
                            !! function to smooth transport at base [nondim]
  ! Local variables
                        !! to the extended mixed-layer bottom [nondim]
                        !! layer to smooth out the parameterized transport [nondim]

  ! Lower order shape (not used), see eq 10 from FK08b.
  ! Apparently used in CM2G, see eq 14 of FK11.
  !mu = max(0., (1. - (2.*sigma + 1.)**2))

  ! Second order, in Rossby number, shape. See eq 21 from FK08a, eq 9 from FK08b, eq 5 FK11
end function mu
module subroutine mixedlayer_restrat_Bodner(CS, G, GV, US, h, uhtr, vhtr, tv, forces, dt, BLD, h_MLD, bflux)
  ! Arguments
  type(mixedlayer_restrat_CS),                intent(inout) :: CS     !< Module control structure
  type(ocean_grid_type),                      intent(inout) :: G      !< Ocean grid structure
  type(verticalGrid_type),                    intent(in)    :: GV     !< Ocean vertical grid structure
  type(unit_scale_type),                      intent(in)    :: US     !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),  intent(inout) :: h      !< Layer thickness [H ~> m or kg m-2]
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), intent(inout) :: uhtr   !< Accumulated zonal mass flux
                                                                      !!   [H L2 ~> m3 or kg]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), intent(inout) :: vhtr   !< Accumulated meridional mass flux
                                                                      !!   [H L2 ~> m3 or kg]
  type(thermo_var_ptrs),                      intent(in)    :: tv     !< Thermodynamic variables structure
  type(mech_forcing),                         intent(in)    :: forces !< A structure with the driving mechanical forces
  real,                                       intent(in)    :: dt     !< Time increment [T ~> s]
  real, dimension(:,:),                       pointer       :: BLD    !< Active boundary layer depth provided by the
                                                                      !! PBL scheme [Z ~> m] (not H)
  real, dimension(:,:),                       pointer       :: h_MLD  !< Thickness of water within the
                                                                      !! active boundary layer depth provided by
                                                                      !! the PBL scheme [H ~> m or kg m-2]
  real, dimension(:,:),                       pointer       :: bflux  !< Surface buoyancy flux provided by the
                                                                      !! PBL scheme [Z2 T-3 ~> m2 s-3]
  ! Local variables
                          ! each layer, divided by dt [H L2 T-1 ~> m3 s-1 or kg s-1]
end subroutine mixedlayer_restrat_Bodner
real elemental module function rmean2ts(signal, filtered, tau_growing, tau_decaying, dt)
  ! Arguments
  real, intent(in) :: signal       ! Unfiltered signal in arbitrary units [A]
  real, intent(in) :: filtered     ! Current value of running mean in the same arbitrary units [A]
  real, intent(in) :: tau_growing  ! Time scale for growing signal [T ~> s]
  real, intent(in) :: tau_decaying ! Time scale for decaying signal [T ~> s]
  real, intent(in) :: dt           ! Time step [T ~> s]
  ! Local variables

end function rmean2ts
module subroutine mixedlayer_restrat_BML(h, uhtr, vhtr, tv, forces, dt, G, GV, US, CS)
  type(ocean_grid_type),                      intent(in)    :: G      !< Ocean grid structure
  type(verticalGrid_type),                    intent(in)    :: GV     !< Ocean vertical grid structure
  type(unit_scale_type),                      intent(in)    :: US     !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),  intent(inout) :: h      !< Layer thickness [H ~> m or kg m-2]
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), intent(inout) :: uhtr   !< Accumulated zonal mass flux
                                                                      !!   [H L2 ~> m3 or kg]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), intent(inout) :: vhtr   !< Accumulated meridional mass flux
                                                                      !!   [H L2 ~> m3 or kg]
  type(thermo_var_ptrs),                      intent(in)    :: tv     !< Thermodynamic variables structure
  type(mech_forcing),                         intent(in)    :: forces !< A structure with the driving mechanical forces
  real,                                       intent(in)    :: dt     !< Time increment [T ~> s]
  type(mixedlayer_restrat_CS),                intent(inout) :: CS     !< Module control structure

  ! Local variables
                          ! sublayer of the mixed layer, divided by dt [H L2 T-1 ~> m3 s-1 or kg s-1].
                          ! in non-Boussinesq mode [H T-1 ~> m s-1 or kg m-2 s-1]
end subroutine mixedlayer_restrat_BML
module subroutine detect_mld(h, tv, MLD_fast, G, GV, CS)
  type(mixedlayer_restrat_CS),                intent(inout) :: CS     !< Module control structure
  type(ocean_grid_type),                      intent(inout) :: G      !< Ocean grid structure
  type(verticalGrid_type),                    intent(in)    :: GV     !< Ocean vertical grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),  intent(inout) :: h      !< Layer thickness [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJ_(G)),           intent(out)   :: MLD_fast !< detected mixed layer depth [H ~> m or kg m-2]
  type(thermo_var_ptrs),                      intent(in)    :: tv     !< Thermodynamic variables structure

  ! Local variables
                                       ! densities [R L2 T-2 ~> Pa].

end subroutine detect_mld
real module function growth_time(u_star, hBL, absf, h_neg, vonKar, Kv_rest, restrat_coef)
  real, intent(in) :: u_star   !< Surface friction velocity in thickness-based units [H T-1 ~> m s-1 or kg m-2 s-1]
  real, intent(in) :: hBL      !< Boundary layer thickness including at least a negligible
                               !! value to keep it positive definite [H ~> m or kg m-2]
  real, intent(in) :: absf     !< Absolute value of the Coriolis parameter [T-1 ~> s-1]
  real, intent(in) :: h_neg    !< A tiny thickness that is usually lost in roundoff so can be
                               !! neglected [H ~> m or kg m-2]
  real, intent(in) :: Kv_rest  !< The background laminar vertical viscosity used for restratification,
                               !! rescaled into thickness-based units [H2 T-1 ~> m2 s-1 or kg2 m-4 s-1]
  real, intent(in) :: vonKar   !< The von Karman constant, used to scale the turbulent limits
                               !! on the restratification timescales [nondim]
  real, intent(in) :: restrat_coef !< An overall scaling factor for the restratification timescale [nondim]

  ! Local variables

  ! peak ML visc: u_star * von_Karman * (h_ml*u_star)/(absf*h_ml + 4.0*u_star) + Kv_water
  ! momentum mixing rate: pi^2*visc/h_ml^2
end function growth_time
logical module function mixedlayer_restrat_init(Time, G, GV, US, param_file, diag, CS, restart_CS)
  type(time_type), target,     intent(in)    :: Time       !< Current model time
  type(ocean_grid_type),       intent(inout) :: G          !< Ocean grid structure
  type(verticalGrid_type),     intent(in)    :: GV         !< Ocean vertical grid structure
  type(unit_scale_type),       intent(in)    :: US         !< A dimensional unit scaling type
  type(param_file_type),       intent(in)    :: param_file !< Parameter file to parse
  type(diag_ctrl), target,     intent(inout) :: diag       !< Regulate diagnostics
  type(mixedlayer_restrat_CS), intent(inout) :: CS         !< Module control structure
  type(MOM_restart_CS),        intent(in)    :: restart_CS !< MOM restart control structure

  ! Local variables
                           ! temperature variance [nondim]
  ! This include declares and sets the variable "version".
                                    ! when reading from file.


  ! Read all relevant parameters and write them to the model log.
end function mixedlayer_restrat_init
module subroutine mixedlayer_restrat_register_restarts(HI, GV, US, param_file, CS, restart_CS)
  ! Arguments
  type(hor_index_type),        intent(in)    :: HI         !< Horizontal index structure
  type(verticalGrid_type),     intent(in)    :: GV         !< Ocean vertical grid structure
  type(unit_scale_type),       intent(in)    :: US         !< A dimensional unit scaling type
  type(param_file_type),       intent(in)    :: param_file !< Parameter file to parse
  type(mixedlayer_restrat_CS), intent(inout) :: CS         !< Module control structure
  type(MOM_restart_CS),        intent(inout) :: restart_CS !< MOM restart control structure

  ! Local variables

  ! Check to see if this module will be used
end subroutine mixedlayer_restrat_register_restarts
logical module function mixedlayer_restrat_unit_tests(verbose)
  logical, intent(in) :: verbose !< If true, write results to stdout

  ! Local variables

end function mixedlayer_restrat_unit_tests
logical module function test_answer(verbose, u, u_true, label, tol)
  logical,            intent(in) :: verbose !< If true, write results to stdout
  real,               intent(in) :: u      !< Values to test in arbitrary units [A]
  real,               intent(in) :: u_true !< Values to test against (correct answer) [A]
  character(len=*),   intent(in) :: label  !< Message
  real, optional,     intent(in) :: tol    !< The tolerance for differences between u and u_true [A]
  ! Local variables

end function test_answer
  end interface

end module MOM_mixed_layer_restrat
