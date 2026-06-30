! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> Does full convective adjustment of unstable regions via a strong diffusivity.
module MOM_full_convection

use MOM_grid,              only : ocean_grid_type
use MOM_interface_heights, only : thickness_to_dz
use MOM_unit_scaling,      only : unit_scale_type
use MOM_variables,         only : thermo_var_ptrs
use MOM_verticalGrid,      only : verticalGrid_type
use MOM_EOS,               only : calculate_density_derivs, EOS_domain

implicit none ; private

#include <MOM_memory.h>

public full_convection


  interface
module subroutine full_convection(G, GV, US, h, tv, T_adj, S_adj, p_surf, Kddt_smooth, halo)
  type(ocean_grid_type),   intent(in)    :: G     !< The ocean's grid structure
  type(verticalGrid_type), intent(in)    :: GV    !< The ocean's vertical grid structure
  type(unit_scale_type),   intent(in)    :: US    !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in)    :: h     !< Layer thicknesses [H ~> m or kg m-2]
  type(thermo_var_ptrs),   intent(in)    :: tv    !< A structure pointing to various
                                                  !! thermodynamic variables
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(out)   :: T_adj !< Adjusted potential temperature [C ~> degC].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(out)   :: S_adj !< Adjusted salinity [S ~> ppt].
  real, dimension(:,:),    pointer       :: p_surf !< The pressure at the ocean surface [R L2 T-2 ~> Pa] (or NULL).
  real,                    intent(in)    :: Kddt_smooth  !< A smoothing vertical diffusivity
                                                  !! times a timestep [H Z ~> m2 or kg m-1].
  integer,                 intent(in)    :: halo  !< Halo width over which to compute

  ! Local variables
                        ! in roundoff and can be neglected [H ~> m or kg m-2].
! logical :: use_EOS    ! If true, density is calculated from T & S using an equation of state.
            ! This array is discretized on tracer cells, but contains an extra
            ! layer at the top for algorithmic convenience.
end subroutine full_convection
module function is_unstable(dRho_dT, dRho_dS, h_a, h_b, mix_A, mix_B, T_a, T_b, S_a, S_b, &
                     Te_aa, Te_bb, Se_aa, Se_bb, d_A, d_B)
  real, intent(in) :: dRho_dT !< The derivative of in situ density with temperature [R C-1 ~> kg m-3 degC-1]
  real, intent(in) :: dRho_dS !< The derivative of in situ density with salinity [R S-1 ~> kg m-3 ppt-1]
  real, intent(in) :: h_a     !< The thickness of the layer above [H ~> m or kg m-2]
  real, intent(in) :: h_b     !< The thickness of the layer below [H ~> m or kg m-2]
  real, intent(in) :: mix_A   !< The time integrated mixing rate of the interface above [H ~> m or kg m-2]
  real, intent(in) :: mix_B   !< The time integrated mixing rate of the interface below [H ~> m or kg m-2]
  real, intent(in) :: T_a     !< The initial temperature of the layer above [C ~> degC]
  real, intent(in) :: T_b     !< The initial temperature of the layer below [C ~> degC]
  real, intent(in) :: S_a     !< The initial salinity of the layer below [S ~> ppt]
  real, intent(in) :: S_b     !< The initial salinity of the layer below [S ~> ppt]
  real, intent(in) :: Te_aa   !< The estimated temperature two layers above rescaled by d_A [C ~> degC]
  real, intent(in) :: Te_bb   !< The estimated temperature two layers below rescaled by d_B [C ~> degC]
  real, intent(in) :: Se_aa   !< The estimated salinity two layers above rescaled by d_A [S ~> ppt]
  real, intent(in) :: Se_bb   !< The estimated salinity two layers below rescaled by d_B [S ~> ppt]
  real, intent(in) :: d_A     !< The rescaling dependency across the interface above [nondim]
  real, intent(in) :: d_B     !< The rescaling dependency across the interface below [nondim]
  logical :: is_unstable !< The return value, true if the profile is statically unstable
                         !! around the interface in question.

  ! These expressions for the local stability are long, but they have been carefully
  ! grouped for accuracy even when the mixing rates are huge or tiny, and common
  ! positive definite factors that would appear in the final expression for the
  ! locally referenced potential density difference across an interface have been omitted.
end function is_unstable
module subroutine smoothed_dRdT_dRdS(h, dz, tv, Kddt, dR_dT, dR_dS, G, GV, US, j, p_surf, halo)
  type(ocean_grid_type),   intent(in)  :: G    !< The ocean's grid structure
  type(verticalGrid_type), intent(in)  :: GV   !< The ocean's vertical grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in)  :: h    !< Layer thicknesses [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZK_(GV)), &
                           intent(in)  :: dz   !< Height change across layers [Z ~> m]
  type(thermo_var_ptrs),   intent(in)  :: tv   !< A structure pointing to various
                                               !! thermodynamic variables
  real,                    intent(in)  :: Kddt !< A diffusivity times a time increment [H Z ~> m2 or kg m-1].
  real, dimension(SZI_(G),SZK_(GV)+1), &
                           intent(out) :: dR_dT !< Derivative of locally referenced
                                               !! potential density with temperature [R C-1 ~> kg m-3 degC-1]
  real, dimension(SZI_(G),SZK_(GV)+1), &
                           intent(out) :: dR_dS !< Derivative of locally referenced
                                               !! potential density with salinity [R S-1 ~> kg m-3 ppt-1]
  type(unit_scale_type),   intent(in)  :: US   !< A dimensional unit scaling type
  integer,                 intent(in)  :: j    !< The j-point to work on.
  real, dimension(:,:),    pointer     :: p_surf !< The pressure at the ocean surface [R L2 T-2 ~> Pa].
  integer,                 intent(in)  :: halo !< Halo width over which to compute

  ! Local variables
                                   ! between layers within in a timestep [H ~> m or kg m-2].
                                   ! [H ~> m or kg m-2].

end subroutine smoothed_dRdT_dRdS
  end interface

end module MOM_full_convection
