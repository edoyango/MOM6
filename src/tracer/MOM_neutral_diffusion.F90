! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> A column-wise toolbox for implementing neutral diffusion
module MOM_neutral_diffusion

use MOM_cpu_clock,             only : cpu_clock_id, cpu_clock_begin, cpu_clock_end
use MOM_cpu_clock,             only : CLOCK_MODULE, CLOCK_ROUTINE
use MOM_domains,               only : pass_var
use MOM_diag_mediator,         only : diag_ctrl, time_type
use MOM_diag_mediator,         only : post_data, register_diag_field
use MOM_EOS,                   only : EOS_type, EOS_manual_init, EOS_domain
use MOM_EOS,                   only : calculate_density, calculate_density_derivs
use MOM_EOS,                   only : EOS_LINEAR
use MOM_error_handler,         only : MOM_error, FATAL, WARNING, MOM_mesg, is_root_pe
use MOM_file_parser,           only : get_param, log_version, param_file_type
use MOM_file_parser,           only : openParameterBlock, closeParameterBlock
use MOM_grid,                  only : ocean_grid_type
use MOM_remapping,             only : remapping_CS, initialize_remapping
use MOM_remapping,             only : extract_member_remapping_CS, build_reconstructions_1d
use MOM_remapping,             only : average_value_ppoly, remappingSchemesDoc, remappingDefaultScheme
use MOM_tracer_registry,       only : tracer_registry_type, tracer_type
use MOM_unit_scaling,          only : unit_scale_type
use MOM_variables,             only : vertvisc_type
use MOM_verticalGrid,          only : verticalGrid_type
use polynomial_functions,      only : evaluation_polynomial, first_derivative_polynomial
use PPM_functions,             only : PPM_reconstruction, PPM_boundary_extrapolation
use regrid_edge_values,        only : edge_values_implicit_h4
use MOM_CVMix_KPP,             only : KPP_get_BLD, KPP_CS
use MOM_energetic_PBL,         only : energetic_PBL_get_MLD, energetic_PBL_CS
use MOM_diabatic_driver,       only : diabatic_CS, extract_diabatic_member
use MOM_io,                    only : stdout, stderr
use MOM_hor_bnd_diffusion,     only : boundary_k_range, SURFACE, BOTTOM

implicit none ; private

#include <MOM_memory.h>

public neutral_diffusion, neutral_diffusion_init, neutral_diffusion_end
public neutral_diffusion_calc_coeffs
public neutral_diffusion_unit_tests

!> The control structure for the MOM_neutral_diffusion module
type, public :: neutral_diffusion_CS ; private
  integer :: nkp1     !< Number of interfaces for a column = nk + 1
  integer :: nsurf    !< Number of neutral surfaces
  integer :: deg = 2  !< Degree of polynomial used for reconstructions
  logical :: continuous_reconstruction = .true. !< True if using continuous PPM reconstruction at interfaces
  logical :: debug = .false. !< If true, write verbose debugging messages
  logical :: hard_fail_heff !< Bring down the model if a problem with heff is detected
  integer :: max_iter !< Maximum number of iterations if refine_position is defined
  real :: drho_tol    !< Convergence criterion representing density difference from true neutrality [R ~> kg m-3]
  real :: x_tol       !< Convergence criterion for how small an update of the position can be [nondim]
  real :: ref_pres    !< Reference pressure, negative if using locally referenced neutral
                      !! density [R L2 T-2 ~> Pa]
  logical :: interior_only !< If true, only applies neutral diffusion in the ocean interior.
                      !! That is, the algorithm will exclude the surface and bottom boundary layers.
  logical :: tapering = .false. !< If true, neutral diffusion linearly decays towards zero within a
                      !! transition zone defined using boundary layer depths. Only available when
                      !! interior_only=true.
  logical :: KhTh_use_vert_struct !< If true, uses vertical structure
                                 !! for tracer diffusivity.
  logical :: use_unmasked_transport_bug !< If true, use an older form for the accumulation of
                      !! neutral-diffusion transports that were unmasked, as used prior to Jan 2018.
  real,    allocatable, dimension(:,:)  :: hbl    !< Boundary layer depth [H ~> m or kg m-2]
  ! Coefficients used to apply tapering from neutral to horizontal direction
  real,    allocatable, dimension(:) :: coeff_l   !< Non-dimensional coefficient in the left column,
                                                  !! at cell interfaces [nondim]
  real,    allocatable, dimension(:) :: coeff_r   !< Non-dimensional coefficient in the right column,
                                                  !! at cell interfaces [nondim]
  ! Array used when KhTh_use_vert_struct is true
  real,    allocatable, dimension(:,:,:) :: Coef_h !< Coef_x and Coef_y averaged at t-points [L2 ~> m2]
  ! Positions of neutral surfaces in both the u, v directions
  real,    allocatable, dimension(:,:,:) :: uPoL  !< Non-dimensional position with left layer uKoL-1, u-point [nondim]
  real,    allocatable, dimension(:,:,:) :: uPoR  !< Non-dimensional position with right layer uKoR-1, u-point [nondim]
  integer, allocatable, dimension(:,:,:) :: uKoL  !< Index of left interface corresponding to neutral surface,
                                                  !! at a u-point
  integer, allocatable, dimension(:,:,:) :: uKoR  !< Index of right interface corresponding to neutral surface,
                                                  !! at a u-point
  real,    allocatable, dimension(:,:,:) :: uHeff !< Effective thickness at u-point [H ~> m or kg m-2]
  real,    allocatable, dimension(:,:,:) :: vPoL  !< Non-dimensional position with left layer uKoL-1, v-point [nondim]
  real,    allocatable, dimension(:,:,:) :: vPoR  !< Non-dimensional position with right layer uKoR-1, v-point [nondim]
  integer, allocatable, dimension(:,:,:) :: vKoL  !< Index of left interface corresponding to neutral surface,
                                                  !! at a v-point
  integer, allocatable, dimension(:,:,:) :: vKoR  !< Index of right interface corresponding to neutral surface,
                                                  !! at a v-point
  real,    allocatable, dimension(:,:,:) :: vHeff !< Effective thickness at v-point [H ~> m or kg m-2]
  ! Coefficients of polynomial reconstructions for temperature and salinity
  real,    allocatable, dimension(:,:,:,:) :: ppoly_coeffs_T !< Polynomial coefficients of the
                                                  !! sub-gridscale temperatures [C ~> degC]
  real,    allocatable, dimension(:,:,:,:) :: ppoly_coeffs_S !< Polynomial coefficients of the
                                                  !! sub-gridscale salinity [S ~> ppt]
  ! Variables needed for continuous reconstructions
  real,    allocatable, dimension(:,:,:) :: dRdT !< dRho/dT [R C-1 ~> kg m-3 degC-1] at interfaces
  real,    allocatable, dimension(:,:,:) :: dRdS !< dRho/dS [R S-1 ~> kg m-3 ppt-1] at interfaces
  real,    allocatable, dimension(:,:,:) :: Tint !< Interface T [C ~> degC]
  real,    allocatable, dimension(:,:,:) :: Sint !< Interface S [S ~> ppt]
  real,    allocatable, dimension(:,:,:) :: Pint !< Interface pressure [R L2 T-2 ~> Pa]
  ! Variables needed for discontinuous reconstructions
  real,    allocatable, dimension(:,:,:,:) :: T_i    !< Top edge reconstruction of temperature [C ~> degC]
  real,    allocatable, dimension(:,:,:,:) :: S_i    !< Top edge reconstruction of salinity [S ~> ppt]
  real,    allocatable, dimension(:,:,:,:) :: P_i    !< Interface pressures [R L2 T-2 ~> Pa]
  real,    allocatable, dimension(:,:,:,:) :: dRdT_i !< dRho/dT [R C-1 ~> kg m-3 degC-1] at top edge
  real,    allocatable, dimension(:,:,:,:) :: dRdS_i !< dRho/dS [R S-1 ~> kg m-3 ppt-1] at top edge
  integer, allocatable, dimension(:,:)     :: ns     !< Number of interfaces in a column
  logical, allocatable, dimension(:,:,:) :: stable_cell !< True if the cell is stably stratified wrt to the next cell
  real :: R_to_kg_m3 = 1.0                   !< A rescaling factor translating density to kg m-3 for
                                             !! use in diagnostic messages [kg m-3 R-1 ~> 1].
  type(diag_ctrl), pointer :: diag => NULL() !< A structure that is used to
                                             !! regulate the timing of diagnostic output.
  integer :: neutral_pos_method              !< Method to find the position of a neutral surface within the layer
  character(len=40)  :: delta_rho_form       !< Determine which (if any) approximation is made to the
                                             !! equation describing the difference in density

  integer :: id_uhEff_2d = -1 !< Diagnostic IDs
  integer :: id_vhEff_2d = -1 !< Diagnostic IDs

  type(EOS_type), pointer :: EOS => NULL()  !< Equation of state parameters
  type(remapping_CS) :: remap_CS   !< Remapping control structure used to create sublayers
  integer :: remap_answer_date     !< The vintage of the order of arithmetic and expressions to use
                                   !! for remapping.  Values below 20190101 recover the remapping
                                   !! answers from 2018, while higher values use more robust
                                   !! forms of the same remapping expressions.
  integer :: ndiff_answer_date     !< The vintage of the order of arithmetic to use for the neutral
                                   !! diffusion.  Values of 20240330 or below recover the answers
                                   !! from the original form of this code, while higher values use
                                   !! mathematically equivalent expressions that recover rotational symmetry.
  type(KPP_CS),           pointer :: KPP_CSp => NULL()          !< KPP control structure needed to get BLD
  type(energetic_PBL_CS), pointer :: energetic_PBL_CSp => NULL()!< ePBL control structure needed to get MLD
end type neutral_diffusion_CS

! This include declares and sets the variable "version".
#include "version_variable.h"
character(len=40)  :: mdl = "MOM_neutral_diffusion" !< module name


  interface
logical module function neutral_diffusion_init(Time, G, GV, US, param_file, diag, EOS, diabatic_CSp, CS)
  type(time_type), target,    intent(in)    :: Time       !< Time structure
  type(ocean_grid_type),      intent(in)    :: G          !< Grid structure
  type(verticalGrid_type),    intent(in)    :: GV         !< The ocean's vertical grid structure
  type(unit_scale_type),      intent(in)    :: US         !< A dimensional unit scaling type
  type(diag_ctrl), target,    intent(inout) :: diag       !< Diagnostics control structure
  type(param_file_type),      intent(in)    :: param_file !< Parameter file structure
  type(EOS_type),  target,    intent(in)    :: EOS        !< Equation of state
  type(diabatic_CS),          pointer       :: diabatic_CSp!< diabatic control structure needed to get BLD
  type(neutral_diffusion_CS), pointer       :: CS         !< Neutral diffusion control structure

  ! Local variables
                                  !! extrapolation should be used within boundary cells.

end function neutral_diffusion_init
module subroutine neutral_diffusion_calc_coeffs(G, GV, US, h, T, S, visc, CS, p_surf)
  type(ocean_grid_type),                     intent(in) :: G   !< Ocean grid structure
  type(verticalGrid_type),                   intent(in) :: GV  !< ocean vertical grid structure
  type(unit_scale_type),                     intent(in) :: US  !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in) :: h   !< Layer thickness [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in) :: T   !< Potential temperature [C ~> degC]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in) :: S   !< Salinity [S ~> ppt]
  type(vertvisc_type),                       intent(in) :: visc !< Structure with vertical viscosities,
                                                               !! boundary layer properties and related fields
  type(neutral_diffusion_CS),                pointer    :: CS  !< Neutral diffusion control structure
  real, dimension(SZI_(G),SZJ_(G)), optional, intent(in) :: p_surf !< Surface pressure to include in pressures used
                                                              !! for equation of state calculations [R L2 T-2 ~> Pa]

  ! Local variables
  ! Variables used for reconstructions
                              ! gradient, which for temperature would be [C H-1 ~> degC m-1 or degC m2 kg-1].
                                                   ! top extent of the boundary layer (0 at top, 1 at bottom) [nondim]
                                       ! (H) units [H T2 R-1 Z-2 ~> m Pa-1 or s2 m-1]

end subroutine neutral_diffusion_calc_coeffs
module subroutine neutral_diffusion(G, GV, h, Coef_x, Coef_y, dt, Reg, US, CS)
  type(ocean_grid_type),                        intent(in)    :: G      !< Ocean grid structure
  type(verticalGrid_type),                      intent(in)    :: GV     !< ocean vertical grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),    intent(in)    :: h      !< Layer thickness [H ~> m or kg m-2]
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)+1), intent(in)    :: Coef_x !< dt * Kh * dy / dx at u-points [L2 ~> m2]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)+1), intent(in)    :: Coef_y !< dt * Kh * dx / dy at v-points [L2 ~> m2]
  real,                                         intent(in)    :: dt     !< Tracer time step * I_numitts [T ~> s]
                                                                        !! (I_numitts is in tracer_hordiff)
  type(tracer_registry_type),                   pointer       :: Reg    !< Tracer registry
  type(unit_scale_type),                        intent(in)    :: US     !< A dimensional unit scaling type
  type(neutral_diffusion_CS),                   pointer       :: CS     !< Neutral diffusion control structure

  ! Local variables
                        ! thickness times a concentration ([C H ~> degC m or degC kg m-2] for temperature) or a
                        ! volume or mass times a concentration ([C H L2 ~> degC m3 or degC kg] for temperature),
                        ! depending on the setting of CS%KhTh_use_vert_struct.
                        ! thickness times a concentration ([C H ~> degC m or degC kg m-2] for temperature) or a
                        ! volume or mass times a concentration ([C H L2 ~> degC m3 or degC kg] for temperature),
                        ! depending on the setting of CS%KhTh_use_vert_struct.
                                                              ! [H conc T-1 ~> m conc s-1 or kg m-2 conc s-1]
                                                              ! For temperature these units are
                                                              ! [C H T-1 ~> degC m s-1 or degC kg m-2 s-1].
                                                              ! [H conc T-1 ~> m conc s-1 or kg m-2 conc s-1].
                                                              ! For temperature these units are
                                                              ! [C H T-1 ~> degC m s-1 or degC kg m-2 s-1].
                                                              ! diagnostic.  For temperature this has units of
                                                              ! [C H L2 ~> degC m3 or degC kg].
                                                              ! diagnostic.  For temperature this has units of
                                                              ! [C H L2 ~> degC m3 or degC kg].
                                                              ! [H L2 conc ~> m3 conc or kg conc].  For temperature
                                                              ! these units are [C H L2 ~> degC m3 or degC kg].
                                                              ! into a cell via its logically northern face, in
                                                              ! [H L2 conc ~> m3 conc or kg conc].
                                                              ! into a cell via its logically southern face, in
                                                              ! [H L2 conc ~> m3 conc or kg conc].
                                                              ! into a cell via its logically eastern face, in
                                                              ! [H L2 conc ~> m3 conc or kg conc].
                                                              ! into a cell via its logically western face, in
                                                              ! [H L2 conc ~> m3 conc or kg conc].


end subroutine neutral_diffusion
module subroutine compute_tapering_coeffs(ne, bld_l, bld_r, coeff_l, coeff_r, h_l, h_r)
  integer,               intent(in)    :: ne       !< Number of interfaces
  real,                  intent(in)    :: bld_l    !< Boundary layer depth, left column  [H ~> m or kg m-2]
  real,                  intent(in)    :: bld_r    !< Boundary layer depth, right column [H ~> m or kg m-2]
  real, dimension(ne-1), intent(in)    :: h_l      !< Layer thickness, left column       [H ~> m or kg m-2]
  real, dimension(ne-1), intent(in)    :: h_r      !< Layer thickness, right column      [H ~> m or kg m-2]
  real, dimension(ne),   intent(inout) :: coeff_l  !< Tapering coefficient, left column            [nondim]
  real, dimension(ne),   intent(inout) :: coeff_r  !< Tapering coefficient, right column           [nondim]

  ! Local variables

  ! Initialize coefficients
end subroutine compute_tapering_coeffs
module subroutine interface_scalar(nk, h, S, Si, i_method, h_neglect)
  integer,               intent(in)    :: nk       !< Number of levels
  real, dimension(nk),   intent(in)    :: h        !< Layer thickness [H ~> m or kg m-2]
  real, dimension(nk),   intent(in)    :: S        !< Layer scalar (or concentrations) in arbitrary
                                                   !! concentration units (e.g. [C ~> degC] for temperature)
  real, dimension(nk+1), intent(inout) :: Si       !< Interface scalar (or concentrations) in arbitrary
                                                   !! concentration units (e.g. [C ~> degC] for temperature)
  integer,               intent(in)    :: i_method !< =1 use average of PLM edges
                                                   !! =2 use continuous PPM edge interpolation
  real,                  intent(in)    :: h_neglect !< A negligibly small thickness [H ~> m or kg m-2]
  ! Local variables
                              ! concentration units (e.g. [C ~> degC] for temperature)
                 ! concentration units (e.g. [C ~> degC] for temperature)

end subroutine interface_scalar
real module function ppm_edge(hkm1, hk, hkp1, hkp2,  Ak, Akp1, Pk, Pkp1, h_neglect)
  real, intent(in) :: hkm1 !< Width of cell k-1 in [H ~> m or kg m-2] or other units
  real, intent(in) :: hk   !< Width of cell k in [H ~> m or kg m-2] or other units
  real, intent(in) :: hkp1 !< Width of cell k+1 in [H ~> m or kg m-2] or other units
  real, intent(in) :: hkp2 !< Width of cell k+2 in [H ~> m or kg m-2] or other units
  real, intent(in) :: Ak   !< Average scalar value of cell k in arbitrary concentration
                           !! units (e.g. [C ~> degC] for temperature)
  real, intent(in) :: Akp1 !< Average scalar value of cell k+1 in arbitrary concentration
                           !! units (e.g. [C ~> degC] for temperature)
  real, intent(in) :: Pk   !< PLM slope for cell k in arbitrary concentration
                           !! units (e.g. [C ~> degC] for temperature)
  real, intent(in) :: Pkp1 !< PLM slope for cell k+1 in arbitrary concentration
                           !! units (e.g. [C ~> degC] for temperature)
  real, intent(in) :: h_neglect !< A negligibly small thickness [H ~> m or kg m-2]

  ! Local variables

end function ppm_edge
real module function ppm_ave(xL, xR, aL, aR, aMean)
  real, intent(in) :: xL    !< Fraction position of left bound (0,1) [nondim]
  real, intent(in) :: xR    !< Fraction position of right bound (0,1) [nondim]
  real, intent(in) :: aL    !< Left edge scalar value, at x=0, in arbitrary concentration
                            !! units (e.g. usually [C ~> degC] for temperature)
  real, intent(in) :: aR    !< Right edge scalar value, at x=1 in arbitrary concentration
                            !! units (e.g. usually [C ~> degC] for temperature)
  real, intent(in) :: aMean !< Average scalar value of cell in arbitrary concentration
                            !! units (e.g. usually [C ~> degC] for temperature)

  ! Local variables
                   ! concentration units as aMean (e.g. usually [C ~> degC] for temperature)

end function ppm_ave
real module function signum(a,x)
  real, intent(in) :: a !< The magnitude argument in arbitrary units [arbitrary]
  real, intent(in) :: x !< The sign (or zero) argument [arbitrary]

end function signum
module subroutine PLM_diff(nk, h, S, c_method, b_method, diff)
  integer,             intent(in)    :: nk       !< Number of levels
  real, dimension(nk), intent(in)    :: h        !< Layer thickness [H ~> m or kg m-2] or other units
  real, dimension(nk), intent(in)    :: S        !< Layer salinity (conc, e.g. ppt) or other tracer
                                                 !! concentration in arbitrary units [A ~> a]
  integer,             intent(in)    :: c_method !< Method to use for the centered difference
  integer,             intent(in)    :: b_method !< =1, use PCM in first/last cell, =2 uses linear extrapolation
  real, dimension(nk), intent(inout) :: diff     !< Scalar difference across layer (conc, e.g. ppt)
                                                 !! in the same arbitrary units as S [A ~> a],
                                                 !! determined by the following values for c_method:
                                                 !!   1. Second order finite difference (not recommended)
                                                 !!   2. Second order finite volume (used in original PPM)
                                                 !!   3. Finite-volume weighted least squares linear fit
                                                 !! \todo  The use of c_method to choose a scheme is inefficient
                                                 !! and should eventually be moved up the call tree.

  ! Local variables

end subroutine PLM_diff
real module function fv_diff(hkm1, hk, hkp1, Skm1, Sk, Skp1)
  real, intent(in) :: hkm1 !< Left cell width [H ~> m or kg m-2] or other arbitrary units
  real, intent(in) :: hk   !< Center cell width [H ~> m or kg m-2] or other arbitrary units
  real, intent(in) :: hkp1 !< Right cell width [H ~> m or kg m-2] or other arbitrary units
  real, intent(in) :: Skm1 !< Left cell average value in arbitrary concentration
                           !! units (e.g. [C ~> degC] for temperature)
  real, intent(in) :: Sk   !< Center cell average value in arbitrary concentration
                           !! units (e.g. [C ~> degC] for temperature)
  real, intent(in) :: Skp1 !< Right cell average value in arbitrary concentration
                           !! units (e.g. [C ~> degC] for temperature)

  ! Local variables
                        ! their reciprocals [H-1 ~> m-1 or m2 kg-1]

end function fv_diff
real module function fvlsq_slope(hkm1, hk, hkp1, Skm1, Sk, Skp1)
  real, intent(in) :: hkm1 !< Left cell width [H ~> m or kg m-2] or other arbitrary units
  real, intent(in) :: hk   !< Center cell width [H ~> m or kg m-2] or other arbitrary units
  real, intent(in) :: hkp1 !< Right cell width [H ~> m or kg m-2] or other arbitrary units
  real, intent(in) :: Skm1 !< Left cell average value in arbitrary concentration
                           !! units (e.g. [C ~> degC] for temperature)
  real, intent(in) :: Sk   !< Center cell average value often in arbitrary concentration
                           !! units (e.g. [C ~> degC] for temperature)
  real, intent(in) :: Skp1 !< Right cell average value often in arbitrary concentration
                           !! units (e.g. [C ~> degC] for temperature)

  ! Local variables
                      ! depend on those of Sk (e.g. [C H2 ~> degC m2 or degC kg2 m-4] for temperature)
                      ! those of Sk (e.g. [C H ~> degC m or degC kg m-2] for temperature)

end function fvlsq_slope
module subroutine find_neutral_surface_positions_continuous(nk, Pl, Tl, Sl, dRdTl, dRdSl, Pr, Tr, Sr, &
                                                     dRdTr, dRdSr, PoL, PoR, KoL, KoR, hEff, bl_kl, bl_kr, bl_zl, bl_zr)
  integer,                    intent(in)    :: nk    !< Number of levels
  real, dimension(nk+1),      intent(in)    :: Pl    !< Left-column interface pressure [R L2 T-2 ~> Pa] or other units
  real, dimension(nk+1),      intent(in)    :: Tl    !< Left-column interface potential temperature [C ~> degC]
  real, dimension(nk+1),      intent(in)    :: Sl    !< Left-column interface salinity [S ~> ppt]
  real, dimension(nk+1),      intent(in)    :: dRdTl !< Left-column dRho/dT [R C-1 ~> kg m-3 degC-1]
  real, dimension(nk+1),      intent(in)    :: dRdSl !< Left-column dRho/dS [R S-1 ~> kg m-3 ppt-1]
  real, dimension(nk+1),      intent(in)    :: Pr    !< Right-column interface pressure [R L2 T-2 ~> Pa] or other units
  real, dimension(nk+1),      intent(in)    :: Tr    !< Right-column interface potential temperature [C ~> degC]
  real, dimension(nk+1),      intent(in)    :: Sr    !< Right-column interface salinity [S ~> ppt]
  real, dimension(nk+1),      intent(in)    :: dRdTr !< Left-column dRho/dT [R C-1 ~> kg m-3 degC-1]
  real, dimension(nk+1),      intent(in)    :: dRdSr !< Left-column dRho/dS [R S-1 ~> kg m-3 ppt-1]
  real, dimension(2*nk+2),    intent(inout) :: PoL   !< Fractional position of neutral surface within
                                                     !! layer KoL of left column [nondim]
  real, dimension(2*nk+2),    intent(inout) :: PoR   !< Fractional position of neutral surface within
                                                     !! layer KoR of right column [nondim]
  integer, dimension(2*nk+2), intent(inout) :: KoL   !< Index of first left interface above neutral surface
  integer, dimension(2*nk+2), intent(inout) :: KoR   !< Index of first right interface above neutral surface
  real, dimension(2*nk+1),    intent(inout) :: hEff  !< Effective thickness between two neutral surfaces
                                                     !! [R L2 T-2 ~> Pa] or other units following Pl and Pr.
  integer, optional,          intent(in)    :: bl_kl !< Layer index of the boundary layer (left)
  integer, optional,          intent(in)    :: bl_kr !< Layer index of the boundary layer (right)
  real, optional,             intent(in)    :: bl_zl !< Fractional position of the boundary layer (left) [nondim]
  real, optional,             intent(in)    :: bl_zr !< Fractional position of the boundary layer (right) [nondim]

  ! Local variables

end subroutine find_neutral_surface_positions_continuous
real module function interpolate_for_nondim_position(dRhoNeg, Pneg, dRhoPos, Ppos)
  real, intent(in) :: dRhoNeg !< Negative density difference [R ~> kg m-3]
  real, intent(in) :: Pneg    !< Position of negative density difference [R L2 T-2 ~> Pa] or [nondim]
  real, intent(in) :: dRhoPos !< Positive density difference [R ~> kg m-3]
  real, intent(in) :: Ppos    !< Position of positive density difference [R L2 T-2 ~> Pa] or [nondim]


end function interpolate_for_nondim_position
module subroutine find_neutral_surface_positions_discontinuous(CS, nk, &
                   Pres_l, hcol_l, Tl, Sl, ppoly_T_l, ppoly_S_l, stable_l, &
                   Pres_r, hcol_r, Tr, Sr, ppoly_T_r, ppoly_S_r, stable_r, &
                   PoL, PoR, KoL, KoR, hEff, zeta_bot_L, zeta_bot_R, k_bot_L, k_bot_R, hard_fail_heff)

  type(neutral_diffusion_CS),     intent(inout) :: CS        !< Neutral diffusion control structure
  integer,                        intent(in)    :: nk        !< Number of levels
  real, dimension(nk,2),          intent(in)    :: Pres_l    !< Left-column interface pressure [R L2 T-2 ~> Pa]
  real, dimension(nk),            intent(in)    :: hcol_l    !< Left-column layer thicknesses [H ~> m or kg m-2]
                                                             !! or other units
  real, dimension(nk,2),          intent(in)    :: Tl        !< Left-column top interface potential
                                                             !! temperature [C ~> degC]
  real, dimension(nk,2),          intent(in)    :: Sl        !< Left-column top interface salinity [S ~> ppt]
  real, dimension(:,:),           intent(in)    :: ppoly_T_l !< Left-column coefficients of T reconstruction [C ~> degC]
  real, dimension(:,:),           intent(in)    :: ppoly_S_l !< Left-column coefficients of S reconstruction [S ~> ppt]
  logical, dimension(nk),         intent(in)    :: stable_l  !< True where the left-column is stable
  real, dimension(nk,2),          intent(in)    :: Pres_r    !< Right-column interface pressure [R L2 T-2 ~> Pa]
  real, dimension(nk),            intent(in)    :: hcol_r    !< Left-column layer thicknesses [H ~> m or kg m-2]
                                                             !! or other units
  real, dimension(nk,2),          intent(in)    :: Tr        !< Right-column top interface potential
                                                             !! temperature [C ~> degC]
  real, dimension(nk,2),          intent(in)    :: Sr        !< Right-column top interface salinity [S ~> ppt]
  real, dimension(:,:),           intent(in)    :: ppoly_T_r !< Right-column coefficients of T
                                                             !! reconstruction [C ~> degC]
  real, dimension(:,:),           intent(in)    :: ppoly_S_r !< Right-column coefficients of S reconstruction [S ~> ppt]
  logical, dimension(nk),         intent(in)    :: stable_r  !< True where the right-column is stable
  real, dimension(4*nk),          intent(inout) :: PoL       !< Fractional position of neutral surface within
                                                             !! layer KoL of left column [nondim]
  real, dimension(4*nk),          intent(inout) :: PoR       !< Fractional position of neutral surface within
                                                             !! layer KoR of right column [nondim]
  integer, dimension(4*nk),       intent(inout) :: KoL       !< Index of first left interface above neutral surface
  integer, dimension(4*nk),       intent(inout) :: KoR       !< Index of first right interface above neutral surface
  real, dimension(4*nk-1),        intent(inout) :: hEff      !< Effective thickness between two neutral surfaces
                                                             !! [H ~> m or kg m-2] or other units taken from hcol_l
  real, optional,                 intent(in)    :: zeta_bot_L!< Non-dimensional distance to where the boundary layer
                                                             !! intersects the cell (left) [nondim]
  real, optional,                 intent(in)    :: zeta_bot_R!< Non-dimensional distance to where the boundary layer
                                                             !! intersects the cell (right) [nondim]

  integer, optional,              intent(in)    :: k_bot_L   !< k-index for the boundary layer (left) [nondim]
  integer, optional,              intent(in)    :: k_bot_R   !< k-index for the boundary layer (right) [nondim]
  logical, optional,              intent(in)    :: hard_fail_heff !< If true (default) bring down the model if the
                                                             !! neutral surfaces ever cross
  ! Local variables
                                    ! is true, but it can take its value from hard_fail_heff.
  ! Initialize variables for the search
end subroutine find_neutral_surface_positions_discontinuous
module subroutine mark_unstable_cells(CS, nk, T, S, P, stable_cell)
  type(neutral_diffusion_CS), intent(inout) :: CS      !< Neutral diffusion control structure
  integer,                intent(in)    :: nk          !< Number of levels in a column
  real, dimension(nk,2),  intent(in)    :: T           !< Temperature at interfaces [C ~> degC]
  real, dimension(nk,2),  intent(in)    :: S           !< Salinity at interfaces [S ~> ppt]
  real, dimension(nk,2),  intent(in)    :: P           !< Pressure at interfaces [R L2 T-2 ~> Pa]
  logical, dimension(nk), intent(  out) :: stable_cell !< True if this cell is unstably stratified


end subroutine mark_unstable_cells
real module function search_other_column(CS, ksurf, pos_last, T_from, S_from, P_from, T_top, S_top, P_top, &
                                  T_bot, S_bot, P_bot, T_poly, S_poly ) result(pos)
  type(neutral_diffusion_CS), intent(in   ) :: CS       !< Neutral diffusion control structure
  integer,                    intent(in   ) :: ksurf    !< Current index of neutral surface
  real,                       intent(in   ) :: pos_last !< Last position within the current layer, used as the lower
                                                        !! bound in the root finding algorithm [nondim]
  real,                       intent(in   ) :: T_from   !< Temperature at the searched from interface [C ~> degC]
  real,                       intent(in   ) :: S_from   !< Salinity    at the searched from interface [S ~> ppt]
  real,                       intent(in   ) :: P_from   !< Pressure at the searched from interface [R L2 T-2 ~> Pa]
  real,                       intent(in   ) :: T_top    !< Temperature at the searched to top interface [C ~> degC]
  real,                       intent(in   ) :: S_top    !< Salinity    at the searched to top interface [S ~> ppt]
  real,                       intent(in   ) :: P_top    !< Pressure at the searched to top interface [R L2 T-2 ~> Pa]
                                                        !! interface [R L2 T-2 ~> Pa]
  real,                       intent(in   ) :: T_bot    !< Temperature at the searched to bottom interface [C ~> degC]
  real,                       intent(in   ) :: S_bot    !< Salinity    at the searched to bottom interface [S ~> ppt]
  real,                       intent(in   ) :: P_bot    !< Pressure at the searched to bottom
                                                        !! interface [R L2 T-2 ~> Pa]
  real, dimension(:),         intent(in   ) :: T_poly   !< Temperature polynomial reconstruction
                                                        !! coefficients [C ~> degC]
  real, dimension(:),         intent(in   ) :: S_poly   !< Salinity    polynomial reconstruction
                                                        !! coefficients [S ~> ppt]
  ! Local variables

  ! Calculate the difference in density at the tops or the bottom
end function search_other_column
module subroutine increment_interface(nk, kl, ki, reached_bottom, searching_this_column, searching_other_column)
  integer, intent(in   )                :: nk                     !< Number of vertical levels
  integer, intent(inout)                :: kl                     !< Current layer (potentially updated)
  integer, intent(inout)                :: ki                     !< Current interface
  logical, intent(inout)                :: reached_bottom         !< Updated when kl == nk and ki == 2
  logical, intent(inout)                :: searching_this_column  !< Updated when kl == nk and ki == 2
  logical, intent(inout)                :: searching_other_column !< Updated when kl == nk and ki == 2

end subroutine increment_interface
module function find_neutral_pos_linear( CS, z0, T_ref, S_ref, dRdT_ref, dRdS_ref, &
                                  dRdT_top, dRdS_top, dRdT_bot, dRdS_bot, ppoly_T, ppoly_S ) result( z )
  type(neutral_diffusion_CS),intent(in) :: CS        !< Control structure with parameters for this module
  real,                      intent(in) :: z0        !< Lower bound of position, also serves as the
                                                     !! initial guess [nondim]
  real,                      intent(in) :: T_ref     !< Temperature at the searched from interface [C ~> degC]
  real,                      intent(in) :: S_ref     !< Salinity at the searched from interface [S ~> ppt]
  real,                      intent(in) :: dRdT_ref  !< dRho/dT at the searched from interface
                                                     !! [R C-1 ~> kg m-3 degC-1]
  real,                      intent(in) :: dRdS_ref  !< dRho/dS at the searched from interface
                                                     !! [R S-1 ~> kg m-3 ppt-1]
  real,                      intent(in) :: dRdT_top  !< dRho/dT at top of layer being searched
                                                     !! [R C-1 ~> kg m-3 degC-1]
  real,                      intent(in) :: dRdS_top  !< dRho/dS at top of layer being searched
                                                     !! [R S-1 ~> kg m-3 ppt-1]
  real,                      intent(in) :: dRdT_bot  !< dRho/dT at bottom of layer being searched
                                                     !! [R C-1 ~> kg m-3 degC-1]
  real,                      intent(in) :: dRdS_bot  !< dRho/dS at bottom of layer being searched
                                                     !! [R S-1 ~> kg m-3 ppt-1]
  real, dimension(:),        intent(in) :: ppoly_T   !< Coefficients of the polynomial reconstruction of T within
                                                     !! the layer to be searched [C ~> degC].
  real, dimension(:),        intent(in) :: ppoly_S   !< Coefficients of the polynomial reconstruction of S within
                                                     !! the layer to be searched [S ~> ppt].
  real                                  :: z         !< Position where drho = 0 [nondim]
  ! Local variables
                     ! layer [R C-1 ~> kg m-3 degC-1]
                     ! layer [R S-1 ~> kg m-3 ppt-1]

end function find_neutral_pos_linear
module function find_neutral_pos_full( CS, z0, T_ref, S_ref, P_ref, P_top, P_bot, ppoly_T, ppoly_S ) result( z )
  type(neutral_diffusion_CS),intent(in) :: CS        !< Control structure with parameters for this module
  real,                      intent(in) :: z0        !< Lower bound of position, also serves as the
                                                     !! initial guess [nondim]
  real,                      intent(in) :: T_ref     !< Temperature at the searched from interface [C ~> degC]
  real,                      intent(in) :: S_ref     !< Salinity at the searched from interface [S ~> ppt]
  real,                      intent(in) :: P_ref     !< Pressure at the searched from interface [R L2 T-2 ~> Pa]
  real,                      intent(in) :: P_top     !< Pressure at top of layer being searched [R L2 T-2 ~> Pa]
  real,                      intent(in) :: P_bot     !< Pressure at bottom of layer being searched [R L2 T-2 ~> Pa]
  real, dimension(:),        intent(in) :: ppoly_T   !< Coefficients of the polynomial reconstruction of T within
                                                     !! the layer to be searched [C ~> degC]
  real, dimension(:),        intent(in) :: ppoly_S   !< Coefficients of the polynomial reconstruction of T within
                                                     !! the layer to be searched [S ~> ppt]
  real                                  :: z         !< Position where drho = 0 [nondim]
  ! Local variables


end function find_neutral_pos_full
module subroutine calc_delta_rho_and_derivs(CS, T1, S1, p1_in, T2, S2, p2_in, drho, &
                                     drdt1_out, drds1_out, drdt2_out, drds2_out )
  type(neutral_diffusion_CS)    :: CS        !< Neutral diffusion control structure
  real,           intent(in   ) :: T1        !< Temperature at point 1 [C ~> degC]
  real,           intent(in   ) :: S1        !< Salinity at point 1 [S ~> ppt]
  real,           intent(in   ) :: p1_in     !< Pressure at point 1 [R L2 T-2 ~> Pa]
  real,           intent(in   ) :: T2        !< Temperature at point 2 [C ~> degC]
  real,           intent(in   ) :: S2        !< Salinity at point 2 [S ~> ppt]
  real,           intent(in   ) :: p2_in     !< Pressure at point 2 [R L2 T-2 ~> Pa]
  real,           intent(  out) :: drho      !< Difference in density between the two points [R ~> kg m-3]
  real, optional, intent(  out) :: dRdT1_out !< drho_dt at point 1 [R C-1 ~> kg m-3 degC-1]
  real, optional, intent(  out) :: dRdS1_out !< drho_ds at point 1 [R S-1 ~> kg m-3 ppt-1]
  real, optional, intent(  out) :: dRdT2_out !< drho_dt at point 2 [R C-1 ~> kg m-3 degC-1]
  real, optional, intent(  out) :: dRdS2_out !< drho_ds at point 2 [R S-1 ~> kg m-3 ppt-1]
  ! Local variables

  ! Use the same reference pressure or the in-situ pressure
end subroutine calc_delta_rho_and_derivs
module function delta_rho_from_derivs( T1, S1, P1, dRdT1, dRdS1, &
                                T2, S2, P2, dRdT2, dRdS2  ) result (drho)
  real :: T1    !< Temperature at point 1 [C ~> degC]
  real :: S1    !< Salinity at point 1 [S ~> ppt]
  real :: P1    !< Pressure at point 1 [R L2 T-2 ~> Pa]
  real :: dRdT1 !< The partial derivative of density with temperature at point 1 [R C-1 ~> kg m-3 degC-1]
  real :: dRdS1 !< The partial derivative of density with salinity at point 1 [R S-1 ~> kg m-3 ppt-1]
  real :: T2    !< Temperature at point 2 [C ~> degC]
  real :: S2    !< Salinity at point 2 [S ~> ppt]
  real :: P2    !< Pressure at point 2 [R L2 T-2 ~> Pa]
  real :: dRdT2 !< The partial derivative of density with temperature at point 2 [R C-1 ~> kg m-3 degC-1]
  real :: dRdS2 !< The partial derivative of density with salinity at point 2 [R S-1 ~> kg m-3 ppt-1]
  ! Local variables
  real :: drho  ! The density difference [R ~> kg m-3]

end function delta_rho_from_derivs
module function absolute_position(n,ns,Pint,Karr,NParr,k_surface)
  integer, intent(in) :: n            !< Number of levels
  integer, intent(in) :: ns           !< Number of neutral surfaces
  real,    intent(in) :: Pint(n+1)    !< Position of interfaces [R L2 T-2 ~> Pa] or other units
  integer, intent(in) :: Karr(ns)     !< Index of interface above position
  real,    intent(in) :: NParr(ns)    !< Non-dimensional position within layer Karr(:) [nondim]
  integer, intent(in) :: k_surface    !< k-interface to query
  real                :: absolute_position !< The absolute position of a location [R L2 T-2 ~> Pa]
                                      !! or other units following Pint
  ! Local variables

end function absolute_position
module function absolute_positions(n,ns,Pint,Karr,NParr)
  integer, intent(in) :: n         !< Number of levels
  integer, intent(in) :: ns        !< Number of neutral surfaces
  real,    intent(in) :: Pint(n+1) !< Position of interface [R L2 T-2 ~> Pa] or other units
  integer, intent(in) :: Karr(ns)  !< Indexes of interfaces about positions
  real,    intent(in) :: NParr(ns) !< Non-dimensional positions within layers Karr(:) [nondim]

  real,  dimension(ns) :: absolute_positions !< Absolute positions [R L2 T-2 ~> Pa]
                                   !! or other units following Pint

  ! Local variables

end function absolute_positions
module subroutine neutral_surface_flux(nk, nsurf, deg, hl, hr, Tl, Tr, PiL, PiR, KoL, KoR, &
                                hEff, Flx, continuous, h_neglect, remap_CS, h_neglect_edge, &
                                coeff_l, coeff_r)
  integer,                      intent(in)    :: nk    !< Number of levels
  integer,                      intent(in)    :: nsurf !< Number of neutral surfaces
  integer,                      intent(in)    :: deg   !< Degree of polynomial reconstructions
  real, dimension(nk),          intent(in)    :: hl    !< Left-column layer thickness [H ~> m or kg m-2]
  real, dimension(nk),          intent(in)    :: hr    !< Right-column layer thickness [H ~> m or kg m-2]
  real, dimension(nk),          intent(in)    :: Tl    !< Left-column layer tracer in arbitrary concentration
                                                       !! units (e.g. [C ~> degC] for temperature)
  real, dimension(nk),          intent(in)    :: Tr    !< Right-column layer tracer in arbitrary concentration
                                                       !! units (e.g. [C ~> degC] for temperature)
  real, dimension(nsurf),       intent(in)    :: PiL   !< Fractional position of neutral surface
                                                       !! within layer KoL of left column [nondim]
  real, dimension(nsurf),       intent(in)    :: PiR   !< Fractional position of neutral surface
                                                       !! within layer KoR of right column [nondim]
  integer, dimension(nsurf),    intent(in)    :: KoL   !< Index of first left interface above neutral surface
  integer, dimension(nsurf),    intent(in)    :: KoR   !< Index of first right interface above neutral surface
  real, dimension(nsurf-1),     intent(in)    :: hEff  !< Effective thickness between two neutral
                                                       !! surfaces [H ~> m or kg m-2]
  real, dimension(nsurf-1),     intent(inout) :: Flx   !< Flux of tracer between pairs of neutral layers
                                                       !! in units  (conc H or conc H L2) that depend on
                                                       !! the presence and units of coeff_l and coeff_r.
                                                       !! If the tracer is temperature, this could have
                                                       !! units of [C H ~> degC m or degC kg m-2] or
                                                       !! [C H L2 ~> degC m3 or degC kg] if coeff_l has
                                                       !! units of [L2 ~> m2]
  logical,                      intent(in)    :: continuous !< True if using continuous reconstruction
  real,                         intent(in)    :: h_neglect !< A negligibly small width for the purpose
                                                       !! of cell reconstructions [H ~> m or kg m-2]
  type(remapping_CS), optional, intent(in)    :: remap_CS !< Remapping control structure used
                                                       !! to create sublayers
  real,               optional, intent(in)    :: h_neglect_edge !< A negligibly small width used for edge value
                                                       !! calculations if continuous is false [H ~> m or kg m-2]
  real, dimension(nk+1), optional, intent(in) :: coeff_l !< Left-column diffusivity  [L2 ~> m2] or [nondim]
  real, dimension(nk+1), optional, intent(in) :: coeff_r !< Right-column diffusivity [L2 ~> m2] or [nondim]

  ! Local variables
                                  ! columns in arbitrary concentration units (e.g. [C ~> degC] for temperature).
                                  ! columns in arbitrary concentration units (e.g. [C ~> degC] for temperature).
                        ! at various positions in the right column in arbitrary
                        ! concentration units (e.g. [C ~> degC] for temperature).
                        ! at various positions in the left column in arbitrary
                        ! concentration units (e.g. [C ~> degC] for temperature).
                        ! over various portions of the right and left columns in arbitrary
                        ! concentration units (e.g. [C ~> degC] for temperature).
                        ! at various positions between the right and left columns in arbitrary
                        ! concentration units (e.g. [C ~> degC] for temperature).
                   ! absent or in units copied from coeff_l and coeff_r [L2 ~> m2] or [nondim]
                              !! units (e.g. [C ~> degC] for temperature)
                              !! units (e.g. [C ~> degC] for temperature)
                              !! units (e.g. [C ~> degC] for temperature)
                              !! units (e.g. [C ~> degC] for temperature)
                              !! units (e.g. [C ~> degC] for temperature)
                              !! units (e.g. [C ~> degC] for temperature)
  ! Discontinuous reconstruction
                              !! units (e.g. [C ~> degC] for temperature)
                              !! units (e.g. [C ~> degC] for temperature)
                              ! sub-gridscale tracer concentrations in the left column, in arbitrary
                              ! concentration units (e.g. [C ~> degC] for temperature)
                              ! sub-gridscale tracer concentrations in the right column, in arbitrary
                              ! concentration units (e.g. [C ~> degC] for temperature)
                              ! gradient, which for temperature would be [C H-1 ~> degC m-1 or degC m2 kg-1].
                              ! gradient, which for temperature would be [C H-1 ~> degC m-1 or degC m2 kg-1].

end subroutine neutral_surface_flux
module subroutine neutral_surface_T_eval(nk, ns, k_sub, Ks, Ps, T_mean, T_int, deg, iMethod, T_poly, &
                                  T_top, T_bot, T_sub, T_top_int, T_bot_int, T_layer)
  integer,                   intent(in   ) :: nk        !< Number of cell averages
  integer,                   intent(in   ) :: ns        !< Number of neutral surfaces
  integer,                   intent(in   ) :: k_sub     !< Index of current neutral layer
  integer, dimension(ns),    intent(in   ) :: Ks        !< List of the layers associated with each neutral surface
  real, dimension(ns),       intent(in   ) :: Ps        !< List of the positions within a layer of each surface [nondim]
  real, dimension(nk),       intent(in   ) :: T_mean    !< Layer average of tracer in arbitrary concentration
                                                        !! units (e.g. [C ~> degC] for temperature)
  real, dimension(nk,2),     intent(in   ) :: T_int     !< Layer interface values of tracer from reconstruction
                                                        !! in concentration units (e.g. [C ~> degC] for temperature)
  integer,                   intent(in   ) :: deg       !< Degree of reconstruction polynomial (e.g. 1 is linear)
  integer,                   intent(in   ) :: iMethod   !< Method of integration to use
  real, dimension(nk,deg+1), intent(in   ) :: T_poly    !< Coefficients of polynomial reconstructions in arbitrary
                                                        !! concentration units (e.g. [C ~> degC] for temperature)
  real,                      intent(  out) :: T_top     !< Tracer value at top (across discontinuity if necessary) in
                                                        !! concentration units (e.g. [C ~> degC] for temperature)
  real,                      intent(  out) :: T_bot     !< Tracer value at bottom (across discontinuity if necessary)
                                                        !! in concentration units (e.g. [C ~> degC] for temperature)
  real,                      intent(  out) :: T_sub     !< Average of the tracer value over the sublayer in arbitrary
                                                        !! concentration units (e.g. [C ~> degC] for temperature)
  real,                      intent(  out) :: T_top_int !< Tracer value at the top interface of a neutral layer in
                                                        !! concentration units (e.g. [C ~> degC] for temperature)
  real,                      intent(  out) :: T_bot_int !< Tracer value at the bottom interface of a neutral layer in
                                                        !! concentration units (e.g. [C ~> degC] for temperature)
  real,                      intent(  out) :: T_layer   !< Cell-average tracer concentration in a layer that
                                                        !! the reconstruction belongs to in concentration
                                                        !! units (e.g. [C ~> degC] for temperature)


end subroutine neutral_surface_T_eval
module subroutine ppm_left_right_edge_values(nk, Tl, Ti, aL, aR)
  integer,                    intent(in)    :: nk !< Number of levels
  real, dimension(nk),        intent(in)    :: Tl !< Layer tracer (conc, e.g. degC) in arbitrary units [A ~> a]
  real, dimension(nk+1),      intent(in)    :: Ti !< Interface tracer (conc, e.g. degC) in arbitrary units [A ~> a]
  real, dimension(nk),        intent(inout) :: aL !< Left edge value of tracer (conc, e.g. degC)
                                                  !! in the same arbitrary units as Tl and Ti [A ~> a]
  real, dimension(nk),        intent(inout) :: aR !< Right edge value of tracer (conc, e.g. degC)
                                                  !! in the same arbitrary units as Tl and Ti [A ~> a]

  ! Setup reconstruction edge values
end subroutine ppm_left_right_edge_values
logical module function neutral_diffusion_unit_tests(verbose)
  logical, intent(in) :: verbose !< If true, write results to stdout

end function neutral_diffusion_unit_tests
logical module function ndiff_unit_tests_continuous(verbose)
  logical, intent(in) :: verbose !< If true, write results to stdout
  ! Local variables

end function ndiff_unit_tests_continuous
logical module function ndiff_unit_tests_discontinuous(verbose)
  logical, intent(in) :: verbose !< It true, write results to stdout
  ! Local variables
                                           ! arbitrary units [arbitrary]
                                           ! left and right columns
                                           ! of the left column or KoR of the right column [nondim]
                                           ! in the same units as hl and hr [arbitrary]

end function ndiff_unit_tests_discontinuous
logical module function test_fv_diff(verbose, hkm1, hk, hkp1, Skm1, Sk, Skp1, Ptrue, title)
  logical,          intent(in) :: verbose !< If true, write results to stdout
  real,             intent(in) :: hkm1  !< Left cell width [nondim]
  real,             intent(in) :: hk    !< Center cell width [nondim]
  real,             intent(in) :: hkp1  !< Right cell width [nondim]
  real,             intent(in) :: Skm1  !< Left cell average value in arbitrary units [arbitrary]
  real,             intent(in) :: Sk    !< Center cell average value in arbitrary units [arbitrary]
  real,             intent(in) :: Skp1  !< Right cell average value in arbitrary units [arbitrary]
  real,             intent(in) :: Ptrue !< True answer in arbitrary units [arbitrary]
  character(len=*), intent(in) :: title !< Title for messages

  ! Local variables

end function test_fv_diff
logical module function test_fvlsq_slope(verbose, hkm1, hk, hkp1, Skm1, Sk, Skp1, Ptrue, title)
  logical,          intent(in) :: verbose !< If true, write results to stdout
  real,             intent(in) :: hkm1  !< Left cell width in arbitrary units [B ~> b]
  real,             intent(in) :: hk    !< Center cell width in arbitrary units [B ~> b]
  real,             intent(in) :: hkp1  !< Right cell width in arbitrary units [B ~> b]
  real,             intent(in) :: Skm1  !< Left cell average value in arbitrary units [A ~> a]
  real,             intent(in) :: Sk    !< Center cell average value in arbitrary units [A ~> a]
  real,             intent(in) :: Skp1  !< Right cell average value in arbitrary units [A ~> a]
  real,             intent(in) :: Ptrue !< True answer in arbitrary units [A B-1 ~> a b-1]
  character(len=*), intent(in) :: title !< Title for messages

  ! Local variables

end function test_fvlsq_slope
logical module function test_ifndp(verbose, rhoNeg, Pneg, rhoPos, Ppos, Ptrue, title)
  logical,          intent(in) :: verbose !< If true, write results to stdout
  real,             intent(in) :: rhoNeg !< Lighter density [R ~> kg m-3]
  real,             intent(in) :: Pneg   !< Interface position of lighter density [nondim]
  real,             intent(in) :: rhoPos !< Heavier density [R ~> kg m-3]
  real,             intent(in) :: Ppos   !< Interface position of heavier density [nondim]
  real,             intent(in) :: Ptrue  !< True answer [nondim]
  character(len=*), intent(in) :: title  !< Title for messages

  ! Local variables

end function test_ifndp
logical module function test_data1d(verbose, nk, Po, Ptrue, title)
  logical,             intent(in) :: verbose !< If true, write results to stdout
  integer,             intent(in) :: nk    !< Number of layers
  real, dimension(nk), intent(in) :: Po    !< Calculated answer [arbitrary]
  real, dimension(nk), intent(in) :: Ptrue !< True answer [arbitrary]
  character(len=*),    intent(in) :: title !< Title for messages

  ! Local variables

end function test_data1d
logical module function test_data1di(verbose, nk, Po, Ptrue, title)
  logical,                intent(in) :: verbose !< If true, write results to stdout
  integer,                intent(in) :: nk    !< Number of layers
  integer, dimension(nk), intent(in) :: Po    !< Calculated answer [arbitrary]
  integer, dimension(nk), intent(in) :: Ptrue !< True answer [arbitrary]
  character(len=*),       intent(in) :: title !< Title for messages

  ! Local variables

end function test_data1di
logical module function test_nsp(verbose, ns, KoL, KoR, pL, pR, hEff, KoL0, KoR0, pL0, pR0, hEff0, title)
  logical,                intent(in) :: verbose !< If true, write results to stdout
  integer,                intent(in) :: ns    !< Number of surfaces
  integer, dimension(ns), intent(in) :: KoL   !< Index of first left interface above neutral surface
  integer, dimension(ns), intent(in) :: KoR   !< Index of first right interface above neutral surface
  real, dimension(ns),    intent(in) :: pL    !< Fractional position of neutral surface within layer
                                              !! KoL of left column [nondim]
  real, dimension(ns),    intent(in) :: pR    !< Fractional position of neutral surface within layer
                                              !! KoR of right column [nondim]
  real, dimension(ns-1),  intent(in) :: hEff  !< Effective thickness between two neutral surfaces [R L2 T-2 ~> Pa]
  integer, dimension(ns), intent(in) :: KoL0  !< Correct value for KoL
  integer, dimension(ns), intent(in) :: KoR0  !< Correct value for KoR
  real, dimension(ns),    intent(in) :: pL0   !< Correct value for pL [nondim]
  real, dimension(ns),    intent(in) :: pR0   !< Correct value for pR [nondim]
  real, dimension(ns-1),  intent(in) :: hEff0 !< Correct value for hEff [R L2 T-2 ~> Pa]
  character(len=*),       intent(in) :: title !< Title for messages

  ! Local variables

end function test_nsp
logical module function compare_nsp_row(KoL, KoR, pL, pR, KoL0, KoR0, pL0, pR0)
  integer,  intent(in) :: KoL   !< Index of first left interface above neutral surface
  integer,  intent(in) :: KoR   !< Index of first right interface above neutral surface
  real,     intent(in) :: pL    !< Fractional position of neutral surface within layer KoL of left column [nondim]
  real,     intent(in) :: pR    !< Fractional position of neutral surface within layer KoR of right column [nondim]
  integer,  intent(in) :: KoL0  !< Correct value for KoL
  integer,  intent(in) :: KoR0  !< Correct value for KoR
  real,     intent(in) :: pL0   !< Correct value for pL [nondim]
  real,     intent(in) :: pR0   !< Correct value for pR [nondim]

end function compare_nsp_row
logical module function test_rnp(expected_pos, test_pos, title)
  real,             intent(in) :: expected_pos !< The expected position [arbitrary]
  real,             intent(in) :: test_pos !< The position returned by the code [arbitrary]
  character(len=*), intent(in) :: title    !< A label for this test
  ! Local variables

end function test_rnp
module subroutine neutral_diffusion_end(CS)
  type(neutral_diffusion_CS), pointer :: CS  !< Neutral diffusion control structure

end subroutine neutral_diffusion_end
  end interface

end module MOM_neutral_diffusion
