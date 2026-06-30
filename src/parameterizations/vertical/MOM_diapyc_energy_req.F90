! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> Calculates the energy requirements of mixing.
module MOM_diapyc_energy_req

!! \author By Robert Hallberg, May 2015

use MOM_diag_mediator, only : diag_ctrl, Time_type, post_data, register_diag_field
use MOM_EOS,           only : calculate_specific_vol_derivs, calculate_density
use MOM_error_handler, only : MOM_error, FATAL, WARNING, MOM_mesg, is_root_pe
use MOM_file_parser,   only : get_param, log_version, param_file_type
use MOM_grid,          only : ocean_grid_type
use MOM_interface_heights, only : thickness_to_dz
use MOM_unit_scaling,  only : unit_scale_type
use MOM_variables,     only : thermo_var_ptrs
use MOM_verticalGrid,  only : verticalGrid_type

implicit none ; private

public diapyc_energy_req_init, diapyc_energy_req_calc, diapyc_energy_req_test, diapyc_energy_req_end

! A note on unit descriptions in comments: MOM6 uses units that can be rescaled for dimensional
! consistency testing. These are noted in comments with units like Z, H, L, and T, along with
! their mks counterparts with notation like "a velocity [Z T-1 ~> m s-1]".  If the units
! vary with the Boussinesq approximation, the Boussinesq variant is given first.

!> This control structure holds parameters for the MOM_diapyc_energy_req module
type, public :: diapyc_energy_req_CS ; private
  logical :: initialized = .false. !< A variable that is here because empty
                               !! structures are not permitted by some compilers.
  real :: test_Kh_scaling      !< A scaling factor for the diapycnal diffusivity [nondim]
  real :: ColHt_scaling        !< A scaling factor for the column height change correction term [nondim]
  real :: VonKar               !< The von Karman coefficient as used in this module [nondim]
  logical :: use_test_Kh_profile !< If true, use the internal test diffusivity profile in place of
                               !! any that might be passed in as an argument.
  type(diag_ctrl), pointer :: diag => NULL() !< A structure that is used to
                               !! regulate the timing of diagnostic output.

  !>@{ Diagnostic IDs
  integer :: id_ERt=-1, id_ERb=-1, id_ERc=-1, id_ERh=-1, id_Kddt=-1, id_Kd=-1
  integer :: id_CHCt=-1, id_CHCb=-1, id_CHCc=-1, id_CHCh=-1
  integer :: id_T0=-1, id_Tf=-1, id_S0=-1, id_Sf=-1, id_N2_0=-1, id_N2_f=-1
  integer :: id_h=-1, id_zInt=-1
  !>@}
end type diapyc_energy_req_CS


  interface
module subroutine diapyc_energy_req_test(h_3d, dt, tv, G, GV, US, CS, Kd_int)
  type(ocean_grid_type),          intent(in)    :: G    !< The ocean's grid structure.
  type(verticalGrid_type),        intent(in)    :: GV   !< The ocean's vertical grid structure.
  type(unit_scale_type),          intent(in)    :: US   !< A dimensional unit scaling type
  real, dimension(G%isd:G%ied,G%jsd:G%jed,GV%ke), &
                                  intent(in)    :: h_3d !< Layer thickness before entrainment [H ~> m or kg m-2].
  type(thermo_var_ptrs),          intent(inout) :: tv   !< A structure containing pointers to any
                                                        !! available thermodynamic fields.
                                                        !! Absent fields have NULL ptrs.
  real,                           intent(in)    :: dt   !< The amount of time covered by this call [T ~> s].
  type(diapyc_energy_req_CS),     pointer       :: CS   !< This module's control structure.
  real, dimension(G%isd:G%ied,G%jsd:G%jed,GV%ke+1), &
                        optional, intent(in)    :: Kd_int !< Interface diffusivities [H Z T-1 ~> m2 s-1 or kg m-1 s-1]

  ! Local variables
                 ! over the layer thicknesses [H Z-1 ~> nondim or kg m-3]
end subroutine diapyc_energy_req_test
module subroutine diapyc_energy_req_calc(h_in, dz_in, T_in, S_in, Kd, energy_Kd, dt, tv, &
                                  G, GV, US, may_print, CS)
  type(ocean_grid_type),    intent(in)    :: G    !< The ocean's grid structure.
  type(verticalGrid_type),  intent(in)    :: GV   !< The ocean's vertical grid structure.
  type(unit_scale_type),    intent(in)    :: US   !< A dimensional unit scaling type
  real, dimension(GV%ke),   intent(in)    :: h_in !< Layer thickness before entrainment,
                                                  !! [H ~> m or kg m-2]
  real, dimension(GV%ke),   intent(in)    :: dz_in !< Vertical distance across layers before
                                                  !! entrainment [Z ~> m]
  real, dimension(GV%ke),   intent(in)    :: T_in !< The layer temperatures [C ~> degC].
  real, dimension(GV%ke),   intent(in)    :: S_in !< The layer salinities [S ~> ppt].
  real, dimension(GV%ke+1), intent(in)    :: Kd   !< The interfaces diapycnal diffusivities
                                                  !! [H Z T-1 ~> m2 s-1 or kg m-1 s-1].
  real,                     intent(in)    :: dt   !< The amount of time covered by this call [T ~> s].
  real,                     intent(out)   :: energy_Kd !< The column-integrated rate of energy
                                                  !! consumption by diapycnal diffusion [R Z L2 T-3 ~> W m-2].
  type(thermo_var_ptrs),    intent(inout) :: tv   !< A structure containing pointers to any
                                                  !! available thermodynamic fields.
                                                  !! Absent fields have NULL ptrs.
  logical,        optional, intent(in)    :: may_print !< If present and true, write out diagnostics
                                                  !! of energy use.
  type(diapyc_energy_req_CS), &
                  optional, pointer       :: CS   !< This module's control structure.

!   This subroutine uses a substantially refactored tridiagonal equation for
! diapycnal mixing of temperature and salinity to estimate the potential energy
! change due to diapycnal mixing in a column of water.  It does this estimate
! 4 different ways, all of which should be equivalent, but reports only one.
! The various estimates are taken because they will later be used as templates
! for other bits of code.

end subroutine diapyc_energy_req_calc
module subroutine find_PE_chg(Kddt_h0, dKddt_h, hp_a, hp_b, Th_a, Sh_a, Th_b, Sh_b, &
                       dT_to_dPE_a, dS_to_dPE_a, dT_to_dPE_b, dS_to_dPE_b, &
                       pres_Z, dT_to_dColHt_a, dS_to_dColHt_a, dT_to_dColHt_b, dS_to_dColHt_b, &
                       PE_chg, dPEc_dKd, dPE_max, dPEc_dKd_0, PE_ColHt_cor)
  real, intent(in)  :: Kddt_h0  !< The previously used diffusivity at an interface times
                                !! the time step and  divided by the average of the
                                !! thicknesses around the interface [H ~> m or kg m-2].
  real, intent(in)  :: dKddt_h  !< The trial change in the diffusivity at an interface times
                                !! the time step and  divided by the average of the
                                !! thicknesses around the interface [H ~> m or kg m-2].
  real, intent(in)  :: hp_a     !< The effective pivot thickness of the layer above the
                                !! interface, given by h_k plus a term that
                                !! is a fraction (determined from the tridiagonal solver) of
                                !! Kddt_h for the interface above [H ~> m or kg m-2].
  real, intent(in)  :: hp_b     !< The effective pivot thickness of the layer below the
                                !! interface, given by h_k plus a term that
                                !! is a fraction (determined from the tridiagonal solver) of
                                !! Kddt_h for the interface above [H ~> m or kg m-2].
  real, intent(in)  :: Th_a     !< An effective temperature times a thickness in the layer
                                !! above, including implicit mixing effects with other
                                !! yet higher layers [C H ~> degC m or degC kg m-2].
  real, intent(in)  :: Sh_a     !< An effective salinity times a thickness in the layer
                                !! above, including implicit mixing effects with other
                                !! yet higher layers [S H ~> ppt m or ppt kg m-2].
  real, intent(in)  :: Th_b     !< An effective temperature times a thickness in the layer
                                !! below, including implicit mixing effects with other
                                !! yet lower layers [C H ~> degC m or degC kg m-2].
  real, intent(in)  :: Sh_b     !< An effective salinity times a thickness in the layer
                                !! below, including implicit mixing effects with other
                                !! yet lower layers [S H ~> ppt m or ppt kg m-2].
  real, intent(in)  :: dT_to_dPE_a !< A factor (pres_lay*mass_lay*dSpec_vol/dT) relating
                                !! a layer's temperature change to the change in column potential
                                !! energy, including all implicit diffusive changes in the
                                !! temperatures of all the layers above [R Z L2 T-2 C-1 ~> J m-2 degC-1].
  real, intent(in)  :: dS_to_dPE_a !< A factor (pres_lay*mass_lay*dSpec_vol/dS) relating
                                !! a layer's salinity change to the change in column potential
                                !! energy, including all implicit diffusive changes in the
                                !! salinities of all the layers above [R Z L2 T-2 S-1 ~> J m-2 ppt-1].
  real, intent(in)  :: dT_to_dPE_b !< A factor (pres_lay*mass_lay*dSpec_vol/dT) relating
                                !! a layer's temperature change to the change in column potential
                                !! energy, including all implicit diffusive changes in the
                                !! temperatures of all the layers below [R Z L2 T-2 C-1 ~> J m-2 degC-1].
  real, intent(in)  :: dS_to_dPE_b !< A factor (pres_lay*mass_lay*dSpec_vol/dS) relating
                                !! a layer's salinity change to the change in column potential
                                !! energy, including all implicit diffusive changes in the
                                !! salinities of all the layers below [R Z L2 T-2 S-1 ~> J m-2 ppt-1].
  real, intent(in)  :: pres_Z   !< The hydrostatic interface pressure, which relates
                                !! the changes in column thickness to the energy that is radiated
                                !! as gravity waves and unavailable to drive mixing [R L2 T-2 ~> J m-3].
  real, intent(in)  :: dT_to_dColHt_a !< A factor (mass_lay*dSColHtc_vol/dT) relating
                                !! a layer's temperature change to the change in column
                                !! height, including all implicit diffusive changes
                                !! in the temperatures of all the layers above [Z C-1 ~> m degC-1].
  real, intent(in)  :: dS_to_dColHt_a !< A factor (mass_lay*dSColHtc_vol/dS) relating
                                !! a layer's salinity change to the change in column
                                !! height, including all implicit diffusive changes
                                !! in the salinities of all the layers above [Z S-1 ~> m ppt-1].
  real, intent(in)  :: dT_to_dColHt_b !< A factor (mass_lay*dSColHtc_vol/dT) relating
                                !! a layer's temperature change to the change in column
                                !! height, including all implicit diffusive changes
                                !! in the temperatures of all the layers below [Z C-1 ~> m degC-1].
  real, intent(in)  :: dS_to_dColHt_b !< A factor (mass_lay*dSColHtc_vol/dS) relating
                                !! a layer's salinity change to the change in column
                                !! height, including all implicit diffusive changes
                                !! in the salinities of all the layers below [Z S-1 ~> m ppt-1].

  real, intent(out) :: PE_chg   !< The change in column potential energy from applying
                                !! Kddt_h at the present interface [R Z L2 T-2 ~> J m-2].
  real, optional, intent(out) :: dPEc_dKd !< The partial derivative of PE_chg with Kddt_h,
                                          !! [R Z L2 T-2 H-1 ~> J m-3 or J kg-1].
  real, optional, intent(out) :: dPE_max  !< The maximum change in column potential energy that could
                                          !! be realized by applying a huge value of Kddt_h at the
                                          !! present interface [R Z L2 T-2 ~> J m-2].
  real, optional, intent(out) :: dPEc_dKd_0 !< The partial derivative of PE_chg with Kddt_h in the
                                            !! limit where Kddt_h = 0 [R Z L2 T-2 H-1 ~> J m-3 or J kg-1].
  real, optional, intent(out) :: PE_ColHt_cor  !< The correction to PE_chg that is made due to a net
                                            !! change in the column height [R Z L2 T-2 ~> J m-2].

  ! Local variables
                   ! for the potential energy changes [H3 R Z L2 T-2 ~> J m or J kg3 m-8].
                     ! for the column height changes [H3 Z ~> m4 or kg3 m-5].

  !   The expression for the change in potential energy used here is derived
  ! from the expression for the final estimates of the changes in temperature
  ! and salinities, and then extensively manipulated to get it into its most
  ! succinct form. The derivation is not necessarily obvious, but it demonstrably
  ! works by comparison with separate calculations of the energy changes after
  ! the tridiagonal solver for the final changes in temperature and salinity are
  ! applied.

end subroutine find_PE_chg
module subroutine diapyc_energy_req_init(Time, G, GV, US, param_file, diag, CS)
  type(time_type),            intent(in)    :: Time        !< model time
  type(ocean_grid_type),      intent(in)    :: G           !< model grid structure
  type(verticalGrid_type),    intent(in)    :: GV          !< ocean vertical grid structure
  type(unit_scale_type),      intent(in)    :: US          !< A dimensional unit scaling type
  type(param_file_type),      intent(in)    :: param_file  !< file to parse for parameter values
  type(diag_ctrl),    target, intent(inout) :: diag        !< structure to regulate diagnostic output
  type(diapyc_energy_req_CS), pointer       :: CS          !< module control structure

! This include declares and sets the variable "version".

end subroutine diapyc_energy_req_init
module subroutine diapyc_energy_req_end(CS)
  type(diapyc_energy_req_CS), pointer :: CS !< Diapycnal energy requirement control structure that
                                            !! will be deallocated in this subroutine.
end subroutine diapyc_energy_req_end
  end interface

end module MOM_diapyc_energy_req
