! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> A simple linear equation of state for sea water with constant coefficients
module MOM_EOS_linear

use MOM_EOS_base_type, only : EOS_base
use MOM_hor_index, only : hor_index_type

implicit none ; private

public linear_EOS
public int_density_dz_linear
public int_spec_vol_dp_linear
public avg_spec_vol_linear

!> The EOS_base implementation of a linear equation of state
type, extends (EOS_base) :: linear_EOS

  real :: Rho_T0_S0 !< The density at T=0, S=0 and p=0 [kg m-3].
  real :: dRho_dT   !< The derivative of density with temperature [kg m-3 degC-1].
  real :: dRho_dS   !< The derivative of density with salinity [kg m-3 ppt-1].
  real :: dRho_dp   !< The derivative of density with pressure [s2 m-2].

contains
  !> Implementation of the in-situ density as an elemental function [kg m-3]
  procedure :: density_elem => density_elem_linear
  !> Implementation of the in-situ density anomaly as an elemental function [kg m-3]
  procedure :: density_anomaly_elem => density_anomaly_elem_linear
  !> Implementation of the in-situ specific volume as an elemental function [m3 kg-1]
  procedure :: spec_vol_elem => spec_vol_elem_linear
  !> Implementation of the in-situ specific volume anomaly as an elemental function [m3 kg-1]
  procedure :: spec_vol_anomaly_elem => spec_vol_anomaly_elem_linear
  !> Implementation of the calculation of derivatives of density
  procedure :: calculate_density_derivs_elem => calculate_density_derivs_elem_linear
  !> Implementation of the calculation of second derivatives of density
  procedure :: calculate_density_second_derivs_elem => calculate_density_second_derivs_elem_linear
  !> Implementation of the calculation of derivatives of specific volume
  procedure :: calculate_specvol_derivs_elem => calculate_specvol_derivs_elem_linear
  !> Implementation of the calculation of compressibility
  procedure :: calculate_compress_elem => calculate_compress_elem_linear
  !> Implementation of the range query function
  procedure :: EOS_fit_range => EOS_fit_range_linear

  !> Instance specific function to set internal parameters
  procedure :: set_params_linear => set_params_linear

  !> Local implementation of generic calculate_density_array for efficiency
  procedure :: calculate_density_array => calculate_density_array_linear
  !> Local implementation of generic calculate_spec_vol_array for efficiency
  procedure :: calculate_spec_vol_array => calculate_spec_vol_array_linear

end type linear_EOS


  interface
real elemental module function density_elem_linear(this, T, S, pressure)
  class(linear_EOS), intent(in) :: this     !< This EOS
  real, intent(in) :: T        !< Potential temperature relative to the surface [degC]
  real, intent(in) :: S        !< Salinity [ppt]
  real, intent(in) :: pressure !< Pressure [Pa]

end function density_elem_linear
real elemental module function density_anomaly_elem_linear(this, T, S, pressure, rho_ref)
  class(linear_EOS), intent(in) :: this     !< This EOS
  real, intent(in) :: T        !< Potential temperature relative to the surface [degC]
  real, intent(in) :: S        !< Salinity [ppt]
  real, intent(in) :: pressure !< Pressure [Pa]
  real, intent(in) :: rho_ref  !< A reference density [kg m-3]

end function density_anomaly_elem_linear
real elemental module function spec_vol_elem_linear(this, T, S, pressure)
  class(linear_EOS), intent(in) :: this     !< This EOS
  real,              intent(in) :: T        !< Potential temperature relative to the surface [degC].
  real,              intent(in) :: S        !< Salinity [ppt].
  real,              intent(in) :: pressure !< Pressure [Pa].

end function spec_vol_elem_linear
real elemental module function spec_vol_anomaly_elem_linear(this, T, S, pressure, spv_ref)
  class(linear_EOS), intent(in) :: this     !< This EOS
  real,              intent(in) :: T        !< Potential temperature relative to the surface [degC].
  real,              intent(in) :: S        !< Salinity [ppt].
  real,              intent(in) :: pressure !< Pressure [Pa].
  real,              intent(in) :: spv_ref  !< A reference specific volume [m3 kg-1].

end function spec_vol_anomaly_elem_linear
elemental module subroutine calculate_density_derivs_elem_linear(this, T, S, pressure, dRho_dT, dRho_dS)
  class(linear_EOS),    intent(in)   :: this     !< This EOS
  real,    intent(in)  :: T        !< Potential temperature relative to the surface [degC].
  real,    intent(in)  :: S        !< Salinity [ppt].
  real,    intent(in)  :: pressure !< Pressure [Pa].
  real,    intent(out) :: drho_dT  !< The partial derivative of density with
                                   !! potential temperature [kg m-3 degC-1].
  real,    intent(out) :: drho_dS  !< The partial derivative of density with
                                   !! salinity [kg m-3 ppt-1].

end subroutine calculate_density_derivs_elem_linear
elemental module subroutine calculate_density_second_derivs_elem_linear(this, T, S, pressure, &
                                  drho_dS_dS, drho_dS_dT, drho_dT_dT, drho_dS_dP, drho_dT_dP)
  class(linear_EOS), intent(in) :: this !< This EOS
  real, intent(in)    :: T           !< Potential temperature relative to the surface [degC].
  real, intent(in)    :: S           !< Salinity [ppt].
  real, intent(in)    :: pressure    !< pressure [Pa].
  real, intent(inout) :: drho_dS_dS  !< The second derivative of density with
                                     !! salinity [kg m-3 ppt-2].
  real, intent(inout) :: drho_dS_dT  !< The second derivative of density with
                                     !! temperature and salinity [kg m-3 ppt-1 degC-1].
  real, intent(inout) :: drho_dT_dT  !< The second derivative of density with
                                     !! temperature [kg m-3 degC-2].
  real, intent(inout) :: drho_dS_dP  !< The second derivative of density with
                                     !! salinity and pressure [kg m-3 ppt-1 Pa-1].
  real, intent(inout) :: drho_dT_dP  !< The second derivative of density with
                                     !! temperature and pressure [kg m-3 degC-1 Pa-1].

end subroutine calculate_density_second_derivs_elem_linear
elemental module subroutine calculate_specvol_derivs_elem_linear(this, T, S, pressure, dSV_dT, dSV_dS)
  class(linear_EOS),  intent(in)    :: this     !< This EOS
  real,               intent(in)    :: T        !< Potential temperature [degC]
  real,               intent(in)    :: S        !< Salinity [ppt]
  real,               intent(in)    :: pressure !< pressure [Pa]
  real,               intent(inout) :: dSV_dS   !< The partial derivative of specific volume with
                                                !! salinity [m3 kg-1 ppt-1]
  real,               intent(inout) :: dSV_dT   !< The partial derivative of specific volume with
                                                !! potential temperature [m3 kg-1 degC-1]
  ! Local variables

  ! Sv = 1.0 / (Rho_T0_S0 + dRho_dT*T + dRho_dS*S)
end subroutine calculate_specvol_derivs_elem_linear
elemental module subroutine calculate_compress_elem_linear(this, T, S, pressure, rho, drho_dp)
  class(linear_EOS), intent(in)  :: this      !< This EOS
  real,              intent(in)  :: T         !< Potential temperature relative to the surface [degC].
  real,              intent(in)  :: S         !< Salinity [ppt].
  real,              intent(in)  :: pressure  !< pressure [Pa].
  real,              intent(out) :: rho       !< In situ density [kg m-3].
  real,              intent(out) :: drho_dp   !< The partial derivative of density with pressure
                                              !! (also the inverse of the square of sound speed)
                                              !! [s2 m-2].

end subroutine calculate_compress_elem_linear
module subroutine avg_spec_vol_linear(T, S, p_t, dp, SpV_avg, start, npts, Rho_T0_S0, dRho_dT, dRho_dS, dRho_dp)
  real, dimension(:), intent(in)    :: T         !< Potential temperature [degC]
  real, dimension(:), intent(in)    :: S         !< Salinity [ppt]
  real, dimension(:), intent(in)    :: p_t       !< Pressure at the top of the layer [Pa]
  real, dimension(:), intent(in)    :: dp        !< Pressure change in the layer [Pa]
  real, dimension(:), intent(inout) :: SpV_avg   !< The vertical average specific volume
                                                 !! in the layer [m3 kg-1]
  integer,            intent(in)    :: start     !< the starting point in the arrays.
  integer,            intent(in)    :: npts      !< the number of values to calculate.
  real,               intent(in)    :: Rho_T0_S0 !< The density at T=0, S=0 [kg m-3]
  real,               intent(in)    :: dRho_dT   !< The derivative of density with temperature
                                                 !! [kg m-3 degC-1]
  real,               intent(in)    :: dRho_dS   !< The derivative of density with salinity
                                                 !! [kg m-3 ppt-1]
  real,               intent(in)    :: dRho_dp   !< The derivative of density with pressure
                                                 !! [s2 m-2]
  ! Local variables

end subroutine avg_spec_vol_linear
module subroutine EoS_fit_range_linear(this, T_min, T_max, S_min, S_max, p_min, p_max)
  class(linear_EOS), intent(in) :: this !< This EOS
  real, optional, intent(out) :: T_min !< The minimum potential temperature over which this EoS is fitted [degC]
  real, optional, intent(out) :: T_max !< The maximum potential temperature over which this EoS is fitted [degC]
  real, optional, intent(out) :: S_min !< The minimum salinity over which this EoS is fitted [ppt]
  real, optional, intent(out) :: S_max !< The maximum salinity over which this EoS is fitted [ppt]
  real, optional, intent(out) :: p_min !< The minimum pressure over which this EoS is fitted [Pa]
  real, optional, intent(out) :: p_max !< The maximum pressure over which this EoS is fitted [Pa]

end subroutine EoS_fit_range_linear
module subroutine set_params_linear(this, Rho_T0_S0, dRho_dT, dRho_dS, dRho_dp)
  class(linear_EOS), intent(inout) :: this !< This EOS
  real, optional,    intent(in)    :: Rho_T0_S0 !< The density at T=0, S=0 [kg m-3]
  real, optional,    intent(in)    :: dRho_dT   !< The derivative of density with temperature,
                                                !! [kg m-3 degC-1]
  real, optional,    intent(in)    :: dRho_dS   !< The derivative of density with salinity,
                                                !! in [kg m-3 ppt-1]
  real, optional,    intent(in)    :: dRho_dp   !< The derivative of density with pressure,
                                                !! in [s2 m-2]

end subroutine set_params_linear
module subroutine int_density_dz_linear(T, S, z_t, z_b, rho_ref, rho_0, G_e, HI, &
                        Rho_T0_S0, dRho_dT, dRho_dS, dRho_dp, dpa, intz_dpa, intx_dpa, inty_dpa, &
                        bathyT, SSH, dz_neglect, MassWghtInterp, Z_0p)
  type(hor_index_type), intent(in)  :: HI        !< The horizontal index type for the arrays.
  real, dimension(HI%isd:HI%ied,HI%jsd:HI%jed), &
                        intent(in)  :: T         !< Potential temperature relative to the surface
                                                 !! [C ~> degC].
  real, dimension(HI%isd:HI%ied,HI%jsd:HI%jed), &
                        intent(in)  :: S         !< Salinity [S ~> ppt].
  real, dimension(HI%isd:HI%ied,HI%jsd:HI%jed), &
                        intent(in)  :: z_t       !< Height at the top of the layer in depth units [Z ~> m].
  real, dimension(HI%isd:HI%ied,HI%jsd:HI%jed), &
                        intent(in)  :: z_b       !< Height at the top of the layer [Z ~> m].
  real,                 intent(in)  :: rho_ref   !< A mean density [R ~> kg m-3], that
                                                 !! is subtracted out to reduce the magnitude of
                                                 !! each of the integrals.
  real,                 intent(in)  :: rho_0     !< A density [R ~> kg m-3], used to calculate
                                                 !! the pressure (as p~=-z*rho_0*G_e) used in
                                                 !! the equation of state.
  real,                 intent(in)  :: G_e       !< The Earth's gravitational acceleration
                                                 !! [L2 Z-1 T-2 ~> m s-2]
  real,                 intent(in)  :: Rho_T0_S0 !< The density at T=0, S=0 [R ~> kg m-3]
  real,                 intent(in)  :: dRho_dT   !< The derivative of density with temperature,
                                                 !! [R C-1 ~> kg m-3 degC-1]
  real,                 intent(in)  :: dRho_dS   !< The derivative of density with salinity,
                                                 !! in [R S-1 ~> kg m-3 ppt-1]
  real,                 intent(in)  :: dRho_dp   !< The derivative of density with pressure,
                                                 !! in [L-2 T2 ~> m-2 s2]
  real, dimension(HI%isd:HI%ied,HI%jsd:HI%jed), &
                        intent(out) :: dpa       !< The change in the pressure anomaly across the
                                                 !! layer [R L2 T-2 ~> Pa]
  real, dimension(HI%isd:HI%ied,HI%jsd:HI%jed), &
              optional, intent(out) :: intz_dpa  !< The integral through the thickness of the layer
                                                 !! of the pressure anomaly relative to the anomaly
                                                 !! at the top of the layer [R L2 Z T-2 ~> Pa m]
  real, dimension(HI%IsdB:HI%IedB,HI%jsd:HI%jed),  &
              optional, intent(out) :: intx_dpa  !< The integral in x of the difference between the
                                                 !! pressure anomaly at the top and bottom of the
                                                 !! layer divided by the x grid spacing [R L2 T-2 ~> Pa]
  real, dimension(HI%isd:HI%ied,HI%JsdB:HI%JedB),  &
              optional, intent(out) :: inty_dpa  !< The integral in y of the difference between the
                                                 !! pressure anomaly at the top and bottom of the
                                                 !! layer divided by the y grid spacing [R L2 T-2 ~> Pa]
  real, dimension(HI%isd:HI%ied,HI%jsd:HI%jed), &
              optional, intent(in)  :: bathyT    !< The depth of the bathymetry [Z ~> m].
  real, dimension(HI%isd:HI%ied,HI%jsd:HI%jed), &
              optional, intent(in)  :: SSH       !< The sea surface height [Z ~> m]
  real,       optional, intent(in)  :: dz_neglect !< A miniscule thickness change [Z ~> m].
  integer,    optional, intent(in)  :: MassWghtInterp !< A flag indicating whether and how to use
                                                 !! mass weighting to interpolate T/S in integrals
  real, dimension(HI%isd:HI%ied,HI%jsd:HI%jed), &
              optional, intent(in)  :: Z_0p      !< The height at which the pressure is 0 [Z ~> m]

  ! Local variables
                     ! 5 sub-column locations [R L2 T-2 ~> Pa]

  ! These array bounds work for the indexing convention of the input arrays, but
  ! on the computational domain defined for the output arrays.
end subroutine int_density_dz_linear
module subroutine int_spec_vol_dp_linear(T, S, p_t, p_b, alpha_ref, HI, Rho_T0_S0, &
               dRho_dT, dRho_dS, dRho_dp, dza, intp_dza, intx_dza, inty_dza, halo_size, &
               bathyP, P_surf, dP_neglect, MassWghtInterp)
  type(hor_index_type), intent(in)  :: HI        !< The ocean's horizontal index type.
  real, dimension(HI%isd:HI%ied,HI%jsd:HI%jed),  &
                        intent(in)  :: T         !< Potential temperature relative to the surface
                                                 !! [C ~> degC].
  real, dimension(HI%isd:HI%ied,HI%jsd:HI%jed),  &
                        intent(in)  :: S         !< Salinity [S ~> ppt].
  real, dimension(HI%isd:HI%ied,HI%jsd:HI%jed),  &
                        intent(in)  :: p_t       !< Pressure at the top of the layer [R L2 T-2 ~> Pa]
  real, dimension(HI%isd:HI%ied,HI%jsd:HI%jed),  &
                        intent(in)  :: p_b       !< Pressure at the top of the layer [R L2 T-2 ~> Pa]
  real,                 intent(in)  :: alpha_ref   !< A mean specific volume that is subtracted out
                            !! to reduce the magnitude of each of the integrals [R-1 ~> m3 kg-1].
                            !! The calculation is mathematically identical with different values of
                            !! alpha_ref, but this reduces the effects of roundoff.
  real,                 intent(in)  :: Rho_T0_S0 !< The density at T=0, S=0 [R ~> kg m-3]
  real,                 intent(in)  :: dRho_dT   !< The derivative of density with temperature
                                                 !! [R C-1 ~> kg m-3 degC-1]
  real,                 intent(in)  :: dRho_dS   !< The derivative of density with salinity,
                                                 !! in [R S-1 ~> kg m-3 ppt-1]
  real,                 intent(in)  :: dRho_dp   !< The derivative of density with pressure,
                                                 !! in [L-2 T2 ~> m-2 s2]
  real, dimension(HI%isd:HI%ied,HI%jsd:HI%jed), &
                        intent(out) :: dza       !< The change in the geopotential anomaly across
                                                 !! the layer [L2 T-2 ~> m2 s-2]
  real, dimension(HI%isd:HI%ied,HI%jsd:HI%jed), &
              optional, intent(out) :: intp_dza  !< The integral in pressure through the layer of the
                                                 !! geopotential anomaly relative to the anomaly at the
                                                 !! bottom of the layer [R L4 T-4 ~> Pa m2 s-2]
  real, dimension(HI%IsdB:HI%IedB,HI%jsd:HI%jed), &
              optional, intent(out) :: intx_dza  !< The integral in x of the difference between the
                                                 !! geopotential anomaly at the top and bottom of
                                                 !! the layer divided by the x grid spacing
                                                 !! [L2 T-2 ~> m2 s-2]
  real, dimension(HI%isd:HI%ied,HI%JsdB:HI%JedB), &
              optional, intent(out) :: inty_dza  !< The integral in y of the difference between the
                                                 !! geopotential anomaly at the top and bottom of
                                                 !! the layer divided by the y grid spacing
                                                 !! [L2 T-2 ~> m2 s-2]
  integer,    optional, intent(in)  :: halo_size !< The width of halo points on which to calculate dza.
  real, dimension(HI%isd:HI%ied,HI%jsd:HI%jed), &
              optional, intent(in)  :: bathyP    !< The pressure at the bathymetry [R L2 T-2 ~> Pa]
  real, dimension(HI%isd:HI%ied,HI%jsd:HI%jed), &
              optional, intent(in)  :: P_surf    !< The pressure at the ocean surface [R L2 T-2 ~> Pa]
  real,       optional, intent(in)  :: dP_neglect !< A miniscule pressure change with
                                                 !! the same units as p_t [R L2 T-2 ~> Pa]
  integer,    optional, intent(in)  :: MassWghtInterp !< A flag indicating whether and how to use
                                                 !! mass weighting to interpolate T/S in integrals
  ! Local variables
                     ! 5 sub-column locations [L2 T-2 ~> m2 s-2]

end subroutine int_spec_vol_dp_linear
module subroutine calculate_density_array_linear(this, T, S, pressure, rho, start, npts, rho_ref)
  class(linear_EOS),  intent(in)  :: this     !< This EOS
  real, dimension(:), intent(in)  :: T        !< Potential temperature relative to the surface [degC]
  real, dimension(:), intent(in)  :: S        !< Salinity [ppt]
  real, dimension(:), intent(in)  :: pressure !< Pressure [Pa]
  real, dimension(:), intent(out) :: rho      !< In situ density [kg m-3]
  integer,            intent(in)  :: start    !< The starting index for calculations
  integer,            intent(in)  :: npts     !< The number of values to calculate
  real,     optional, intent(in)  :: rho_ref  !< A reference density [kg m-3]

  ! Local variables

end subroutine calculate_density_array_linear
module subroutine calculate_spec_vol_array_linear(this, T, S, pressure, specvol, start, npts, spv_ref)
  class(linear_EOS),  intent(in) :: this      !< This EOS
  real, dimension(:), intent(in)  :: T        !< Potential temperature relative to the surface [degC]
  real, dimension(:), intent(in)  :: S        !< Salinity [ppt]
  real, dimension(:), intent(in)  :: pressure !< Pressure [Pa]
  real, dimension(:), intent(out) :: specvol  !< In situ specific volume [m3 kg-1]
  integer,            intent(in)  :: start    !< The starting index for calculations
  integer,            intent(in)  :: npts     !< The number of values to calculate
  real,     optional, intent(in)  :: spv_ref  !< A reference specific volume [m3 kg-1]

  ! Local variables

end subroutine calculate_spec_vol_array_linear
  end interface

end module MOM_EOS_linear
