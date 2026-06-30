! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> Diapycnal mixing and advection in isopycnal mode
module MOM_entrain_diffusive

use MOM_diag_mediator, only : post_data, register_diag_field, safe_alloc_ptr
use MOM_diag_mediator, only : diag_ctrl, time_type
use MOM_EOS,           only : calculate_density, calculate_density_derivs
use MOM_EOS,           only : calculate_specific_vol_derivs, EOS_domain
use MOM_error_handler, only : MOM_error, is_root_pe, FATAL, WARNING, NOTE
use MOM_file_parser,   only : get_param, log_version, param_file_type
use MOM_forcing_type,  only : forcing
use MOM_grid,          only : ocean_grid_type
use MOM_unit_scaling,  only : unit_scale_type
use MOM_variables,     only : thermo_var_ptrs
use MOM_verticalGrid,  only : verticalGrid_type

implicit none ; private

#include <MOM_memory.h>

public entrainment_diffusive, entrain_diffusive_init

! A note on unit descriptions in comments: MOM6 uses units that can be rescaled for dimensional
! consistency testing. These are noted in comments with units like Z, H, L, and T, along with
! their mks counterparts with notation like "a velocity [Z T-1 ~> m s-1]".  If the units
! vary with the Boussinesq approximation, the Boussinesq variant is given first.

!> The control structure holding parametes for the MOM_entrain_diffusive module
type, public :: entrain_diffusive_CS ; private
  logical :: initialized = .false. !< True if this control structure has been initialized.
  logical :: bulkmixedlayer  !< If true, a refined bulk mixed layer is used with
                             !! GV%nk_rho_varies variable density mixed & buffer layers.
  integer :: max_ent_it      !< The maximum number of iterations that may be used to
                             !! calculate the diapycnal entrainment.
  real    :: Tolerance_Ent   !< The tolerance with which to solve for entrainment values
                             !! [H ~> m or kg m-2].
  real    :: max_Ent         !< A large ceiling on the maximum permitted amount of entrainment
                             !! across each interface between the mixed and buffer layers within
                             !! a timestep [H ~> m or kg m-2].
  real    :: Rho_sig_off     !< The offset between potential density and a sigma value [R ~> kg m-3]
  type(diag_ctrl), pointer :: diag => NULL() !< A structure that is used to
                             !! regulate the timing of diagnostic output.
  integer :: id_Kd = -1      !< Diagnostic ID for diffusivity
  integer :: id_diff_work = -1 !< Diagnostic ID for mixing work
end type entrain_diffusive_CS


  interface
module subroutine entrainment_diffusive(h, tv, fluxes, dt, G, GV, US, CS, ea, eb, &
                                 kb_out, Kd_Lay, Kd_int)
  type(ocean_grid_type),      intent(in)  :: G  !< The ocean's grid structure.
  type(verticalGrid_type),    intent(in)  :: GV !< The ocean's vertical grid structure.
  type(unit_scale_type),      intent(in)  :: US !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),  &
                              intent(in)  :: h  !< Layer thicknesses [H ~> m or kg m-2].
  type(thermo_var_ptrs),      intent(in)  :: tv !< A structure containing pointers to any available
                                                !! thermodynamic fields. Absent fields have NULL
                                                !! ptrs.
  type(forcing),              intent(in)  :: fluxes !< A structure of surface fluxes that may
                                                !! be used.
  real,                       intent(in)  :: dt !< The time increment [T ~> s].
  type(entrain_diffusive_CS), intent(in)  :: CS !< The control structure returned by a previous
                                                !! call to entrain_diffusive_init.
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),  &
                              intent(out) :: ea !< The amount of fluid entrained from the layer
                                                !! above within this time step [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),  &
                              intent(out) :: eb !< The amount of fluid entrained from the layer
                                                !! below within this time step [H ~> m or kg m-2].
  integer, dimension(SZI_(G),SZJ_(G)),        &
                            intent(inout) :: kb_out !< The index of the lightest layer denser than
                                                !! the buffer layer.
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),  &
                              intent(in)  :: Kd_Lay !< The diapycnal diffusivity of layers
                                                !! [H Z T-1 ~> m2 s-1 or kg m-1 s-1].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)+1), &
                              intent(in)  :: Kd_int !< The diapycnal diffusivity of interfaces
                                                !! [H Z T-1 ~> m2 s-1 or kg m-1 s-1].

!   This subroutine calculates ea and eb, the rates at which a layer entrains
! from the layers above and below.  The entrainment rates are proportional to
! the buoyancy flux in a layer and inversely proportional to the density
! differences between layers.  The scheme that is used here is described in
! detail in Hallberg, Mon. Wea. Rev. 2000.

end subroutine entrainment_diffusive
module subroutine F_to_ent(F, h, kb, kmb, j, G, GV, CS, dsp1_ds, eakb, Ent_bl, ea, eb)
  type(ocean_grid_type),            intent(in)    :: G    !< The ocean's grid structure
  type(verticalGrid_type),          intent(in)    :: GV   !< The ocean's vertical grid structure
  real, dimension(SZI_(G),SZK_(GV)), intent(in)   :: F    !< The density flux through a layer within
                                                          !! a time step divided by the density
                                                          !! difference across the interface below
                                                          !! the layer [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                                    intent(in)    :: h    !< Layer thicknesses [H ~> m or kg m-2]
  integer, dimension(SZI_(G)),      intent(in)    :: kb   !< The index of the lightest layer denser than
                                                          !! the deepest buffer layer.
  integer,                          intent(in)    :: kmb  !< The number of mixed and buffer layers.
  integer,                          intent(in)    :: j    !< The meridional index upon which to work.
  type(entrain_diffusive_CS),       intent(in)    :: CS   !< This module's control structure.
  real, dimension(SZI_(G),SZK_(GV)), intent(in)   :: dsp1_ds !< The ratio of coordinate variable
                                                          !! differences across the interfaces below
                                                          !! a layer over the difference across the
                                                          !! interface above the layer [nondim].
  real, dimension(SZI_(G)),         intent(in)    :: eakb !< The entrainment from above by the layer
                                                          !! below the buffer layer [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZK_(GV)), intent(in)   :: Ent_bl !< The average entrainment upward and
                                                          !! downward across each interface around
                                                          !! the buffer layers [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                                    intent(inout) :: ea   !< The amount of fluid entrained from the layer
                                                          !! above within this time step [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                                    intent(inout) :: eb   !< The amount of fluid entrained from the layer
                                                          !! below within this time step [H ~> m or kg m-2].

                    ! after exchange with the layer below [H ~> m or kg m-2].

end subroutine F_to_ent
module subroutine set_Ent_bl(h, dtKd_int, tv, kb, kmb, do_i, G, GV, US, CS, j, Ent_bl, Sref, h_bl)
  type(ocean_grid_type),            intent(in)    :: G    !< The ocean's grid structure.
  type(verticalGrid_type),          intent(in)    :: GV   !< The ocean's vertical grid structure.
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                                    intent(in)    :: h    !< Layer thicknesses [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZK_(GV)+1), &
                                    intent(in)    :: dtKd_int !< The diapycnal diffusivity across
                                                          !! each interface times the time step
                                                          !! [H2 ~> m2 or kg2 m-4].
  type(thermo_var_ptrs),            intent(in)    :: tv   !< A structure containing pointers to any
                                                          !! available thermodynamic fields. Absent
                                                          !! fields have NULL ptrs.
  integer, dimension(SZI_(G)),      intent(inout) :: kb   !< The index of the lightest layer denser
                                                          !! than the buffer layer or 1 if there is
                                                          !! no buffer layer.
  integer,                          intent(in)    :: kmb  !< The number of mixed and buffer layers.
  logical, dimension(SZI_(G)),      intent(in)    :: do_i !< A logical variable indicating which
                                                          !! i-points to work on.
  type(unit_scale_type),            intent(in)    :: US   !< A dimensional unit scaling type
  type(entrain_diffusive_CS),       intent(in)    :: CS   !< This module's control structure.
  integer,                          intent(in)    :: j    !< The meridional index upon which to work.
  real, dimension(SZI_(G),SZK_(GV)+1), &
                                    intent(out)   :: Ent_bl !< The average entrainment upward and
                                                          !! downward across each interface around
                                                          !! the buffer layers [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZK_(GV)), intent(out)  :: Sref !< The coordinate potential density minus
                                                          !! 1000 for each layer [R ~> kg m-3].
  real, dimension(SZI_(G),SZK_(GV)), intent(out)  :: h_bl !< The thickness of each layer [H ~> m or kg m-2].

!   This subroutine sets the average entrainment across each of the interfaces
! between buffer layers within a timestep. It also causes thin and relatively
! light interior layers to be entrained by the deepest buffer layer.
!   Also find the initial coordinate potential densities (Sref) of each layer.
! Does there need to be limiting when the layers below are all thin?

  ! Local variables
end subroutine set_Ent_bl
module subroutine determine_dSkb(h_bl, Sref, Ent_bl, E_kb, is, ie, kmb, G, GV, limit, &
                          dSkb, ddSkb_dE, dSlay, ddSlay_dE, dS_anom_lim, do_i_in)
  type(ocean_grid_type),              intent(in)    :: G      !< The ocean's grid structure.
  type(verticalGrid_type),            intent(in)    :: GV     !< The ocean's vertical grid
                                                              !! structure.
  real, dimension(SZI_(G),SZK_(GV)),  intent(in)    :: h_bl   !< Layer thickness [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZK_(GV)),  intent(in)    :: Sref   !< Reference potential density [R ~> kg m-3]
  real, dimension(SZI_(G),SZK_(GV)),  intent(in)    :: Ent_bl !< The average entrainment upward and
                                                              !! downward across each interface
                                                              !! around the buffer layers [H ~> m or kg m-2].
  real, dimension(SZI_(G)),           intent(in)    :: E_kb   !< The entrainment by the top interior
                                                              !! layer [H ~> m or kg m-2].
  integer,                            intent(in)    :: is     !< The start of the i-index range to work on.
  integer,                            intent(in)    :: ie     !< The end of the i-index range to work on.
  integer,                            intent(in)    :: kmb    !< The number of mixed and buffer layers.
  logical,                            intent(in)    :: limit  !< If true, limit dSkb and dSlay to
                                                              !! avoid negative values.
  real, dimension(SZI_(G)),           intent(inout) :: dSkb   !< The limited potential density
                                                              !! difference across the interface
                                                              !! between the bottommost buffer layer
                                                              !! and the topmost interior layer. [R ~> kg m-3]
                                                              !! dSkb > 0.
  real, dimension(SZI_(G)), optional, intent(inout) :: ddSkb_dE !< The partial derivative of dSkb
                                                              !! with E [R H-1 ~> kg m-4 or m-1].
  real, dimension(SZI_(G)), optional, intent(inout) :: dSlay  !< The limited potential density
                                                              !! difference across the topmost
                                                              !! interior layer. 0 < dSkb [R ~> kg m-3]
  real, dimension(SZI_(G)), optional, intent(inout) :: ddSlay_dE !< The partial derivative of dSlay
                                                              !! with E [R H-1 ~> kg m-4 or m-1].
  real, dimension(SZI_(G)), optional, intent(inout) :: dS_anom_lim !< A limiting value to use for
                                                              !! the density anomalies below the
                                                              !! buffer layer [R ~> kg m-3].
  logical, dimension(SZI_(G)), optional, intent(in) :: do_i_in !< If present, determines which
                                                              !! columns are worked on.

! Note that dSkb, ddSkb_dE, dSlay, ddSlay_dE, and dS_anom_lim are declared
! intent inout  because they should not change where do_i_in is false.

!   This subroutine determines the reference density difference between the
! bottommost buffer layer and the first interior after the mixing between mixed
! and buffer layers and mixing with the layer below. Within the mixed and buffer
! layers, entrainment from the layer above is increased when it is necessary to
! keep the layers from developing a negative thickness; otherwise it equals
! Ent_bl.  At each interface, the upward and downward fluxes average out to
! Ent_bl, unless entrainment by the layer below is larger than twice Ent_bl.
!   The density difference across the first interior layer may also be returned.
! It could also be limited to avoid negative values or values that greatly
! exceed the density differences across an interface.
!   Additionally, the partial derivatives of dSkb and dSlay with E_kb could
! also be returned.

  ! Local variables
                    ! the buffer layers with the new density of the bottommost buffer layer [nondim]
                    ! after exchange with the layer below [H ~> m or kg m-2].
                    ! in roundoff and can be neglected [H ~> m or kg m-2].
                    ! added to ensure positive definiteness [H ~> m or kg m-2].
                       ! layer with the density difference across the interface above it [nondim]

end subroutine determine_dSkb
module subroutine F_kb_to_ea_kb(h_bl, Sref, Ent_bl, I_dSkbp1, F_kb, kmb, i, &
                         G, GV, CS, ea_kb, tol_in)
  type(ocean_grid_type),    intent(in)    :: G    !< The ocean's grid structure
  type(verticalGrid_type),  intent(in)    :: GV   !< The ocean's vertical grid structure
  real, dimension(SZI_(G),SZK_(GV)), &
                            intent(in)    :: h_bl !< Layer thickness, with the top interior
                                                  !! layer at k-index kmb+1 [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZK_(GV)), &
                            intent(in)    :: Sref !< The coordinate reference potential density,
                                                  !! with the value of the topmost interior layer
                                                  !! at index kmb+1 [R ~> kg m-3].
  real, dimension(SZI_(G),SZK_(GV)), &
                            intent(in)    :: Ent_bl !< The average entrainment upward and downward
                                                  !! across each interface around the buffer layers,
                                                  !! [H ~> m or kg m-2].
  real, dimension(SZI_(G)), intent(in)    :: I_dSkbp1 !< The inverse of the difference in reference
                                                  !! potential density across the base of the
                                                  !! uppermost interior layer [R-1 ~> m3 kg-1].
  real, dimension(SZI_(G)), intent(in)    :: F_kb !< The entrainment from below by the
                                                  !! uppermost interior layer [H ~> m or kg m-2]
  integer,                  intent(in)    :: kmb  !< The number of mixed and buffer layers.
  integer,                  intent(in)    :: i    !< The i-index to work on
  type(entrain_diffusive_CS), intent(in)  :: CS   !< This module's control structure.
  real, dimension(SZI_(G)), intent(inout) :: ea_kb !< The entrainment from above by the layer below
                                                  !! the buffer layer (i.e. layer kb) [H ~> m or kg m-2].
  real,           optional, intent(in)    :: tol_in !< A tolerance for the iterative determination
                                                  !! of the entrainment [H ~> m or kg m-2].

                         !  between the bottommost buffer layer and the topmost interior layer [R ~> kg m-3]
                         ! differences) found in the range min_ent < ent < max_ent [H ~> m or kg m-2].

end subroutine F_kb_to_ea_kb
module subroutine determine_Ea_kb(h_bl, dtKd_kb, Sref, I_dSkbp1, Ent_bl, ea_kbp1, &
                           min_eakb, max_eakb, kmb, is, ie, do_i, G, GV, CS, Ent, &
                           error, err_min_eakb0, err_max_eakb0, F_kb, dFdfm_kb)
  type(ocean_grid_type),            intent(in)  :: G        !< The ocean's grid structure.
  type(verticalGrid_type),          intent(in)  :: GV       !< The ocean's vertical grid structure.
  real, dimension(SZI_(G),SZK_(GV)), intent(in) :: h_bl     !< Layer thickness, with the top interior
                                                            !! layer at k-index kmb+1 [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZK_(GV)), intent(in) :: Sref     !< The coordinate reference potential
                                                            !! density, with the value of the
                                                            !! topmost interior layer at layer
                                                            !! kmb+1 [R ~> kg m-3].
  real, dimension(SZI_(G),SZK_(GV)), intent(in) :: Ent_bl   !< The average entrainment upward and
                                                            !! downward across each interface around
                                                            !! the buffer layers [H ~> m or kg m-2].
  real, dimension(SZI_(G)),         intent(in)  :: I_dSkbp1 !< The inverse of the difference in
                                                            !! reference potential density across
                                                            !! the base of the uppermost interior
                                                            !! layer [R-1 ~> m3 kg-1].
  real, dimension(SZI_(G)),         intent(in)  :: dtKd_kb  !< The diapycnal diffusivity in the top
                                                            !! interior layer times the time step
                                                            !! [H2 ~> m2 or kg2 m-4].
  real, dimension(SZI_(G)),         intent(in)  :: ea_kbp1  !< The entrainment from above by layer
                                                            !! kb+1 [H ~> m or kg m-2].
  real, dimension(SZI_(G)),         intent(in)  :: min_eakb !< The minimum permissible rate of
                                                            !! entrainment [H ~> m or kg m-2].
  real, dimension(SZI_(G)),         intent(in)  :: max_eakb !< The maximum permissible rate of
                                                            !! entrainment [H ~> m or kg m-2].
  integer,                          intent(in)  :: kmb      !< The number of mixed and buffer layers.
  integer,                          intent(in)  :: is       !< The start of the i-index range to work on.
  integer,                          intent(in)  :: ie       !< The end of the i-index range to work on.
  logical, dimension(SZI_(G)),      intent(in)  :: do_i     !< A logical variable indicating which
                                                            !! i-points to work on.
  type(entrain_diffusive_CS),       intent(in)  :: CS       !< This module's control structure.
  real, dimension(SZI_(G)),         intent(inout) :: Ent    !< The entrainment rate of the uppermost
                                                            !! interior layer [H ~> m or kg m-2].
                                                            !! The input value is the first guess.
  real, dimension(SZI_(G)), optional, intent(out) :: error  !< The error (locally defined in this
                                                            !! routine) associated with the returned
                                                            !! solution [H2 ~> m2 or kg2 m-4]
  real, dimension(SZI_(G)), optional, intent(in)  :: err_min_eakb0 !< The errors (locally defined)
                                                            !! associated with min_eakb when ea_kbp1 = 0,
                                                            !! returned from a previous call to this
                                                            !! subroutine [H2 ~> m2 or kg2 m-4].
  real, dimension(SZI_(G)), optional, intent(in)  :: err_max_eakb0 !< The errors (locally defined)
                                                            !! associated with min_eakb when ea_kbp1 = 0,
                                                            !! returned from a previous call to this
                                                            !! subroutine [H2 ~> m2 or kg2 m-4].
  real, dimension(SZI_(G)), optional, intent(out) :: F_kb   !< The entrainment from below by the
                                                            !! uppermost interior layer
                                                            !! corresponding to the returned
                                                            !! value of Ent [H ~> m or kg m-2].
  real, dimension(SZI_(G)), optional, intent(out) :: dFdfm_kb !< The partial derivative of F_kb with
                                                            !! ea_kbp1 [nondim].

!  This subroutine determines the entrainment from above by the top interior
! layer (labeled kb elsewhere) given an entrainment by the layer below it,
! constrained to be within the provided bounds.

  ! Local variables
                            ! ensure that it is positive [R ~> kg m-3].
end subroutine determine_Ea_kb
module subroutine find_maxF_kb(h_bl, Sref, Ent_bl, I_dSkbp1, min_ent_in, max_ent_in, &
                        kmb, is, ie, G, GV, CS, maxF, ent_maxF, do_i_in, &
                        F_lim_maxent, F_thresh)
  type(ocean_grid_type),      intent(in)  :: G        !< The ocean's grid structure.
  type(verticalGrid_type),    intent(in)  :: GV       !< The ocean's vertical grid structure.
  real, dimension(SZI_(G),SZK_(GV)), &
                              intent(in)  :: h_bl     !< Layer thickness [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZK_(GV)), &
                              intent(in)  :: Sref     !< Reference potential density [R ~> kg m-3].
  real, dimension(SZI_(G),SZK_(GV)), &
                              intent(in)  :: Ent_bl   !< The average entrainment upward and
                                                      !! downward across each interface around
                                                      !! the buffer layers [H ~> m or kg m-2].
  real, dimension(SZI_(G)),   intent(in)  :: I_dSkbp1 !< The inverse of the difference in
                                                      !! reference potential density across the
                                                      !! base of the uppermost interior layer
                                                      !! [R-1 ~> m3 kg-1].
  real, dimension(SZI_(G)),   intent(in)  :: min_ent_in !< The minimum value of ent to search,
                                                      !! [H ~> m or kg m-2].
  real, dimension(SZI_(G)),   intent(in)  :: max_ent_in !< The maximum value of ent to search,
                                                      !! [H ~> m or kg m-2].
  integer,                    intent(in)  :: kmb      !< The number of mixed and buffer layers.
  integer,                    intent(in)  :: is       !< The start of the i-index range to work on.
  integer,                    intent(in)  :: ie       !< The end of the i-index range to work on.
  type(entrain_diffusive_CS), intent(in)  :: CS       !< This module's control structure.
  real, dimension(SZI_(G)),   intent(out) :: maxF     !< The maximum value of F
                                                      !! = ent*ds_kb*I_dSkbp1 found in the range
                                                      !! min_ent < ent < max_ent [H ~> m or kg m-2].
  real, dimension(SZI_(G)), &
                    optional, intent(out) :: ent_maxF !< The value of ent at that maximum [H ~> m or kg m-2].
  logical, dimension(SZI_(G)), &
                    optional, intent(in)  :: do_i_in  !< A logical array indicating which columns
                                                      !! to work on.
  real, dimension(SZI_(G)), &
                    optional, intent(out) :: F_lim_maxent !< If present, do not apply the limit in
                                                      !! finding the maximum value, but return the
                                                      !! limited value at ent=max_ent_in in this
                                                      !! array [H ~> m or kg m-2].
  real, dimension(SZI_(G)), &
                    optional, intent(in)  :: F_thresh !< If F_thresh is present, return the first value
                                                      !! found that has F > F_thresh [H ~> m or kg m-2], or
                                                      !! the maximum root if it is absent.

! Maximize F = ent*ds_kb*I_dSkbp1 in the range min_ent < ent < max_ent.
! ds_kb may itself be limited to positive values in determine_dSkb, which gives
! the prospect of two local maxima in the range - one at max_ent_in with that
! minimum value of ds_kb, and the other due to the unlimited (potentially
! negative) value.  It is faster to find the true maximum by first finding the
! unlimited maximum and comparing it to the limited value at max_ent_in.
end subroutine find_maxF_kb
module subroutine entrain_diffusive_init(Time, G, GV, US, param_file, diag, CS, just_read_params)
  type(time_type),         intent(in)    :: Time !< The current model time.
  type(ocean_grid_type),   intent(in)    :: G    !< The ocean's grid structure.
  type(verticalGrid_type), intent(in)    :: GV   !< The ocean's vertical grid structure.
  type(unit_scale_type),   intent(in)    :: US   !< A dimensional unit scaling type
  type(param_file_type),   intent(in)    :: param_file !< A structure to parse for run-time
                                                 !! parameters.
  type(diag_ctrl), target, intent(inout) :: diag !< A structure that is used to regulate diagnostic
                                                 !! output.
  type(entrain_diffusive_CS), intent(inout) :: CS !< Entrainment diffusion control structure
  logical,                 intent(in)    :: just_read_params !< If true, this call will only read
                                                 !! and log parameters without registering
                                                 !! any diagnostics

  ! Local variables
  ! This include declares and sets the variable "version".

end subroutine entrain_diffusive_init
  end interface

end module MOM_entrain_diffusive
