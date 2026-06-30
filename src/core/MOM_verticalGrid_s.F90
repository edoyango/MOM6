submodule (MOM_verticalGrid) MOM_verticalGrid_s
#include <MOM_memory.h>
  implicit none
contains
module procedure verticalGridInit
  integer :: nk, H_power
  real    :: H_rescale_factor ! The integer power of 2 by which thicknesses are rescaled [nondim]
  real    :: rho_Kv  ! The density used convert input kinematic viscosities into dynamic viscosities
# include "version_variable.h"
  character(len=16) :: mdl = 'MOM_verticalGrid'
  if (associated(GV)) call MOM_error(FATAL, &
     'verticalGridInit: called with an associated GV pointer.')
  allocate(GV)

  ! Read all relevant parameters and write them to the model log.
  call log_version(param_file, mdl, version, &
                   "Parameters providing information about the vertical grid.", &
                   log_to_all=.true., debugging=.true.)
  call get_param(param_file, mdl, "G_EARTH", GV%g_Earth, &
                 "The gravitational acceleration of the Earth.", &
                 units="m s-2", default=9.80, scale=US%Z_to_m*US%m_s_to_L_T**2)
  call get_param(param_file, mdl, "RHO_0", GV%Rho0, &
                 "The mean ocean density used with BOUSSINESQ true to "//&
                 "calculate accelerations and the mass for conservation "//&
                 "properties, or with BOUSSINESQ false to convert some "//&
                 "parameters from vertical units of m to kg m-2.", &
                 units="kg m-3", default=1035.0, scale=US%kg_m3_to_R)
  call get_param(param_file, mdl, "BOUSSINESQ", GV%Boussinesq, &
                 "If true, make the Boussinesq approximation.", default=.true.)
  call get_param(param_file, mdl, "SEMI_BOUSSINESQ", GV%semi_Boussinesq, &
                 "If true, do non-Boussinesq pressure force calculations and use mass-based "//&
                 "thicknesses, but use RHO_0 to convert layer thicknesses into certain "//&
                 "height changes.  This only applies if BOUSSINESQ is false.", &
                 default=.true., do_not_log=GV%Boussinesq)
  if (GV%Boussinesq) GV%semi_Boussinesq = .true.
  call get_param(param_file, mdl, "RHO_KV_CONVERT", Rho_Kv, &
                 "The density used to convert input vertical distances into thickesses in "//&
                 "non-BOUSSINESQ mode, and to convert kinematic viscosities into dynamic "//&
                 "viscosities and similarly for vertical diffusivities.  GV%m_to_H is set "//&
                 "using this value, whereas GV%Z_to_H is set using RHO_0.  The default is "//&
                 "RHO_0, but this can be set separately to demonstrate the independence of the "//&
                 "non-Boussinesq solutions of the value of RHO_0.", &
                 units="kg m-3", default=GV%Rho0*US%R_to_kg_m3, scale=US%kg_m3_to_R, &
                 do_not_log=GV%Boussinesq)
  call get_param(param_file, mdl, "ANGSTROM", GV%Angstrom_Z, &
                 "The minimum layer thickness, usually one-Angstrom.", &
                 units="m", default=1.0e-10, scale=US%m_to_Z)
  call get_param(param_file, mdl, "H_RESCALE_POWER", H_power, &
                 "An integer power of 2 that is used to rescale the model's "//&
                 "intenal units of thickness.  Valid values range from -300 to 300.", &
                 units="nondim", default=0, debuggingParam=.true.)
  if (abs(H_power) > 300) call MOM_error(FATAL, "verticalGridInit: "//&
                 "H_RESCALE_POWER is outside of the valid range of -300 to 300.")
  H_rescale_factor = 1.0
  if (H_power /= 0) H_rescale_factor = 2.0**H_power
  if (.not.GV%Boussinesq) then
    call get_param(param_file, mdl, "H_TO_KG_M2", GV%H_to_kg_m2,&
                 "A constant that translates thicknesses from the model's "//&
                 "internal units of thickness to kg m-2.", units="kg m-2 H-1", &
                 default=1.0)
    GV%H_to_kg_m2 = GV%H_to_kg_m2 * H_rescale_factor
  else
    call get_param(param_file, mdl, "H_TO_M", GV%H_to_m, &
                 "A constant that translates the model's internal "//&
                 "units of thickness into m.", units="m H-1", default=1.0)
    GV%H_to_m = GV%H_to_m * H_rescale_factor
  endif
  ! This is not used:  GV%mks_g_Earth = US%L_T_to_m_s**2*US%m_to_Z * GV%g_Earth
  GV%g_Earth_Z_T2 = US%L_to_Z**2 * GV%g_Earth  ! This would result from scale=US%m_to_Z*US%T_to_s**2.
#ifdef STATIC_MEMORY_
  ! Here NK_ is a macro, while nk is a variable.
  call get_param(param_file, mdl, "NK", nk, &
                 "The number of model layers.", units="nondim", &
                 default=NK_)
  if (nk /= NK_) call MOM_error(FATAL, "verticalGridInit: " // &
       "Mismatched number of layers NK_ between MOM_memory.h and param_file")

#else
  call get_param(param_file, mdl, "NK", nk, &
                 "The number of model layers.", units="nondim", fail_if_missing=.true.)
#endif
  GV%ke = nk

  if (GV%Boussinesq) then
    GV%H_to_kg_m2 = US%R_to_kg_m3*GV%Rho0 * GV%H_to_m
    GV%kg_m2_to_H = 1.0 / GV%H_to_kg_m2
    GV%m_to_H = 1.0 / GV%H_to_m
    GV%H_to_MKS = GV%H_to_m
    GV%m2_s_to_HZ_T = GV%m_to_H * US%m_to_Z * US%T_to_s

    GV%H_to_Z = GV%H_to_m * US%m_to_Z
    GV%Z_to_H = US%Z_to_m * GV%m_to_H
  else
    GV%kg_m2_to_H = 1.0 / GV%H_to_kg_m2
    !  GV%m_to_H = US%R_to_kg_m3*GV%Rho0 * GV%kg_m2_to_H
    GV%m_to_H = US%R_to_kg_m3*rho_Kv * GV%kg_m2_to_H
    GV%H_to_MKS = GV%H_to_kg_m2
    GV%m2_s_to_HZ_T = US%R_to_kg_m3*rho_Kv * GV%kg_m2_to_H * US%m_to_Z * US%T_to_s
    GV%H_to_m = 1.0 / GV%m_to_H

    GV%H_to_Z = US%m_to_Z * ( GV%H_to_kg_m2 / (US%R_to_kg_m3*GV%Rho0) )
    GV%Z_to_H = US%Z_to_m * ( US%R_to_kg_m3*GV%Rho0 * GV%kg_m2_to_H )
  endif

  GV%Angstrom_H = (US%Z_to_m * GV%m_to_H) * GV%Angstrom_Z
  GV%Angstrom_m = US%Z_to_m * GV%Angstrom_Z

  GV%H_subroundoff = 1e-20 * max(GV%Angstrom_H, GV%m_to_H*1e-17)
  GV%dZ_subroundoff = 1e-20 * max(GV%Angstrom_Z, US%m_to_Z*1e-17)

  GV%H_to_Pa = US%L_T_to_m_s**2*US%m_to_Z * GV%g_Earth * GV%H_to_kg_m2

  GV%H_to_RZ = GV%H_to_kg_m2 * US%kg_m3_to_R * US%m_to_Z
  GV%RZ_to_H = GV%kg_m2_to_H * US%R_to_kg_m3 * US%Z_to_m

  GV%HZ_T_to_m2_s = 1.0 / GV%m2_s_to_HZ_T
  GV%HZ_T_to_MKS = GV%H_to_MKS * US%Z_to_m * US%s_to_T

  ! Note based on the above that for both Boussinsq and non-Boussinesq cases that:
  !     GV%Rho0 = GV%Z_to_H * GV%H_to_RZ
  !     1.0/GV%Rho0 = GV%H_to_Z * GV%RZ_to_H
  ! This is exact for power-of-2 scaling of the units, regardless of the value of Rho0, but
  ! the first term on the right hand side is invertable in Boussinesq mode, but the second
  ! is invertable when non-Boussinesq.

  ! Log derivative values.
  call log_param(param_file, mdl, "M to THICKNESS", GV%m_to_H*H_rescale_factor, units="H m-1")
  call log_param(param_file, mdl, "M to THICKNESS rescaled by 2^-n", GV%m_to_H, units="2^n H m-1")
  call log_param(param_file, mdl, "THICKNESS to M rescaled by 2^n", GV%H_to_m, units="2^-n m H-1")

  allocate( GV%sInterface(nk+1) )
  allocate( GV%sLayer(nk) )
  allocate( GV%g_prime(nk+1), source=0.0 )
  allocate( GV%Rlay(nk), source=0.0 )

end procedure verticalGridInit
module procedure get_thickness_units
  if (GV%Boussinesq) then
    get_thickness_units = "m"
  else
    get_thickness_units = "kg m-2"
  endif
end procedure get_thickness_units
module procedure get_flux_units
  if (GV%Boussinesq) then
    get_flux_units = "m3 s-1"
  else
    get_flux_units = "kg s-1"
  endif
end procedure get_flux_units
module procedure get_tr_flux_units
  integer :: cnt
  cnt = 0
  if (present(tr_units)) cnt = cnt+1
  if (present(tr_vol_conc_units)) cnt = cnt+1
  if (present(tr_mass_conc_units)) cnt = cnt+1

  if (cnt == 0) call MOM_error(FATAL, "get_tr_flux_units: One of the three "//&
    "arguments tr_units, tr_vol_conc_units, or tr_mass_conc_units "//&
    "must be present.")
  if (cnt > 1) call MOM_error(FATAL, "get_tr_flux_units: Only one of "//&
    "tr_units, tr_vol_conc_units, and tr_mass_conc_units may be present.")
  if (present(tr_units)) then
    if (GV%Boussinesq) then
      get_tr_flux_units = trim(tr_units)//" m3 s-1"
    else
      get_tr_flux_units = trim(tr_units)//" kg s-1"
    endif
  endif
  if (present(tr_vol_conc_units)) then
    if (GV%Boussinesq) then
      get_tr_flux_units = trim(tr_vol_conc_units)//" s-1"
    else
      get_tr_flux_units = trim(tr_vol_conc_units)//" m-3 kg s-1"
    endif
  endif
  if (present(tr_mass_conc_units)) then
    if (GV%Boussinesq) then
      get_tr_flux_units = trim(tr_mass_conc_units)//" kg-1 m3 s-1"
    else
      get_tr_flux_units = trim(tr_mass_conc_units)//" s-1"
    endif
  endif

end procedure get_tr_flux_units
module procedure setVerticalGridAxes
  integer :: k, nk
  nk = GV%ke

  GV%zAxisLongName = 'Target Potential Density'
  GV%zAxisUnits = 'kg m-3'
  do k=1,nk ; GV%sLayer(k) = scale*Rlay(k) ; enddo
  if (nk > 1) then
    GV%sInterface(1) = scale * (1.5*Rlay(1) - 0.5*Rlay(2))
    do K=2,nk ; GV%sInterface(K) = scale * 0.5*( Rlay(k-1) + Rlay(k) ) ; enddo
    GV%sInterface(nk+1) = scale * (1.5*Rlay(nk) - 0.5*Rlay(nk-1))
  else
    GV%sInterface(1) = 0.0 ; GV%sInterface(nk+1) = 2.0*scale*Rlay(nk)
  endif

end procedure setVerticalGridAxes
module procedure verticalGridEnd
  deallocate( GV%g_prime, GV%Rlay )
  deallocate( GV%sInterface , GV%sLayer )
  deallocate( GV )

end procedure verticalGridEnd
end submodule MOM_verticalGrid_s
