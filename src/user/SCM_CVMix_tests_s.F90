submodule (SCM_CVMix_tests) SCM_CVMix_tests_s
#include <MOM_memory.h>
  implicit none
contains
module procedure SCM_CVMix_tests_TS_init
  real :: UpperLayerTempMLD !< Upper layer Temp MLD thickness [Z ~> m].
  real :: UpperLayerSaltMLD !< Upper layer Salt MLD thickness [Z ~> m].
  real :: UpperLayerTemp !< Upper layer temperature (SST if thickness 0) [C ~> degC]
  real :: UpperLayerSalt !< Upper layer salinity (SSS if thickness 0) [S ~> ppt]
  real :: LowerLayerTemp !< Temp at top of lower layer [C ~> degC]
  real :: LowerLayerSalt !< Salt at top of lower layer [S ~> ppt]
  real :: LowerLayerdTdz !< Temp gradient in lower layer [C Z-1 ~> degC m-1].
  real :: LowerLayerdSdz !< Salt gradient in lower layer [S Z-1 ~> ppt m-1].
  real :: LowerLayerMinTemp !< Minimum temperature in lower layer [C ~> degC]
  real :: zC, DZ, top, bottom ! Depths and thicknesses [Z ~> m].
  integer :: i, j, k, is, ie, js, je, isd, ied, jsd, jed, nz
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed

  if (.not.just_read) call log_version(param_file, mdl, version)
  call get_param(param_file, mdl, "SCM_TEMP_MLD", UpperLayerTempMLD, &
                 'Initial temp mixed layer depth', &
                 units='m', default=0.0, scale=US%m_to_Z, do_not_log=just_read)
  call get_param(param_file, mdl, "SCM_SALT_MLD", UpperLayerSaltMLD, &
                 'Initial salt mixed layer depth', &
                 units='m', default=0.0, scale=US%m_to_Z, do_not_log=just_read)
  call get_param(param_file, mdl, "SCM_L1_SALT", UpperLayerSalt, &
                 'Layer 2 surface salinity', units="ppt", default=35.0, scale=US%ppt_to_S, do_not_log=just_read)
  call get_param(param_file, mdl, "SCM_L1_TEMP", UpperLayerTemp, &
                 'Layer 1 surface temperature', units="degC", default=20.0, scale=US%degC_to_C, do_not_log=just_read)
  call get_param(param_file, mdl, "SCM_L2_SALT", LowerLayerSalt, &
                 'Layer 2 surface salinity', units="ppt", default=35.0, scale=US%ppt_to_S, do_not_log=just_read)
  call get_param(param_file, mdl, "SCM_L2_TEMP", LowerLayerTemp, &
                 'Layer 2 surface temperature', units="degC", default=20.0, scale=US%degC_to_C, do_not_log=just_read)
  call get_param(param_file, mdl, "SCM_L2_DTDZ", LowerLayerdTdZ,     &
                 'Initial temperature stratification in layer 2', &
                 units='C/m', default=0.0, scale=US%degC_to_C*US%Z_to_m, do_not_log=just_read)
  call get_param(param_file, mdl, "SCM_L2_DSDZ", LowerLayerdSdZ,  &
                 'Initial salinity stratification in layer 2', &
                 units='PPT/m', default=0.0, scale=US%ppt_to_S*US%Z_to_m, do_not_log=just_read)
  call get_param(param_file, mdl, "SCM_L2_MINTEMP",LowerLayerMinTemp, &
                 'Layer 2 minimum temperature', units="degC", default=4.0, scale=US%degC_to_C, do_not_log=just_read)

  if (just_read) return ! All run-time parameters have been read, so return.

  do j=js,je ; do i=is,ie
    top = 0. ! Reference to surface
    bottom = 0.
    do k=1,nz
      bottom = bottom - h(i,j,k)       ! Interface below layer [Z ~> m]
      zC = 0.5*( top + bottom )        ! Z of middle of layer [Z ~> m]
      DZ = min(0., zC + UpperLayerTempMLD)
      T(i,j,k) = max(LowerLayerMinTemp,LowerLayerTemp + LowerLayerdTdZ * DZ)
      DZ = min(0., zC + UpperLayerSaltMLD)
      S(i,j,k) = LowerLayerSalt + LowerLayerdSdZ * DZ
      top = bottom
    enddo ! k
  enddo ; enddo

end procedure SCM_CVMix_tests_TS_init
module procedure SCM_CVMix_tests_surface_forcing_init
# include "version_variable.h"
  type(unit_scale_type), pointer :: US => NULL() !< A dimensional unit scaling type
  US => G%US

  if (associated(CS)) then
    call MOM_error(FATAL, "SCM_CVMix_tests_surface_forcing_init called with an associated "// &
                          "control structure.")
    return
  endif
  allocate(CS)

  ! Read all relevant parameters and write them to the model log.
  call log_version(param_file, mdl, version, "")
  call get_param(param_file, mdl, "SCM_USE_WIND_STRESS", CS%UseWindStress, &
                 "Wind Stress switch used in the SCM CVMix surface forcing.", &
                 default=.false.)
  call get_param(param_file, mdl, "SCM_USE_HEAT_FLUX", CS%UseHeatFlux, &
                 "Heat flux switch used in the SCM CVMix test surface forcing.", &
                 default=.false.)
  call get_param(param_file, mdl, "SCM_USE_EVAPORATION", CS%UseEvaporation, &
                 "Evaporation switch used in the SCM CVMix test surface forcing.", &
                 default=.false.)
  call get_param(param_file, mdl, "SCM_USE_DIURNAL_SW", CS%UseDiurnalSW, &
                 "Diurnal sw radation switch used in the SCM CVMix test surface forcing.", &
                 default=.false.)
  if (CS%UseWindStress) then
    call get_param(param_file, mdl, "SCM_TAU_X", CS%tau_x, &
                 "Constant X-dir wind stress used in the SCM CVMix test surface forcing.", &
                 units='N/m2', scale=US%kg_m2s_to_RZ_T*US%m_s_to_L_T, fail_if_missing=.true.)
    call get_param(param_file, mdl, "SCM_TAU_Y", CS%tau_y, &
                 "Constant y-dir wind stress used in the SCM CVMix test surface forcing.", &
                 units='N/m2', scale=US%kg_m2s_to_RZ_T*US%m_s_to_L_T, fail_if_missing=.true.)
  endif
  if (CS%UseHeatFlux) then
    call get_param(param_file, mdl, "SCM_HEAT_FLUX", CS%surf_HF, &
                 "Constant surface heat flux used in the SCM CVMix test surface forcing.", &
                 units='m K/s', scale=US%m_to_Z*US%degC_to_C*US%T_to_s, fail_if_missing=.true.)
  endif
  if (CS%UseEvaporation) then
    call get_param(param_file, mdl, "SCM_EVAPORATION", CS%surf_evap, &
                 "Constant surface evaporation used in the SCM CVMix test surface forcing.", &
                 units='m/s', scale=US%m_to_Z*US%T_to_s, fail_if_missing=.true.)
  endif
  if (CS%UseDiurnalSW) then
    call get_param(param_file, mdl, "SCM_DIURNAL_SW_MAX", CS%Max_sw, &
                 "Maximum diurnal sw radiation used in the SCM CVMix test surface forcing.", &
                 units='m K/s', scale=US%m_to_Z*US%degC_to_C*US%T_to_s, fail_if_missing=.true.)
  endif
  call get_param(param_file, mdl, "RHO_0", CS%Rho0, &
                 "The mean ocean density used with BOUSSINESQ true to "//&
                 "calculate accelerations and the mass for conservation "//&
                 "properties, or with BOUSSINESQ false to convert some "//&
                 "parameters from vertical units of m to kg m-2.", &
                 units="kg m-3", default=1035.0, scale=US%kg_m3_to_R)
  call get_param(param_file, mdl, "RESTORE_FLUX_RHO", CS%rho_restore, &
                 "The density that is used to convert piston velocities into salt or heat fluxes.", &
                 units="kg m-3", default=CS%Rho0*US%R_to_kg_m3, scale=US%kg_m3_to_R)

end procedure SCM_CVMix_tests_surface_forcing_init
module procedure SCM_CVMix_tests_wind_forcing
  integer :: i, j, is, ie, js, je, Isq, Ieq, Jsq, Jeq
  integer :: isd, ied, jsd, jed, IsdB, IedB, JsdB, JedB
  real    :: mag_tau  ! The magnitude of the wind stress [R Z2 T-2 ~> Pa]
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec
  Isq = G%IscB ; Ieq = G%IecB ; Jsq = G%JscB ; Jeq = G%JecB
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  IsdB = G%IsdB ; IedB = G%IedB ; JsdB = G%JsdB ; JedB = G%JedB

  do j=js,je ; do I=Isq,Ieq
    forces%taux(I,j) = CS%tau_x
  enddo ; enddo
  do J=Jsq,Jeq ; do i=is,ie
    forces%tauy(i,J) = CS%tau_y
  enddo ; enddo
  call pass_vector(forces%taux, forces%tauy, G%Domain, To_All)

  mag_tau = US%L_to_Z * sqrt((CS%tau_x*CS%tau_x) + (CS%tau_y*CS%tau_y))
  if (associated(forces%ustar)) then ; do j=js,je ; do i=is,ie
    forces%ustar(i,j) = sqrt( mag_tau / CS%Rho0 )
  enddo ; enddo ; endif

  if (associated(forces%tau_mag)) then ; do j=js,je ; do i=is,ie
    forces%tau_mag(i,j) = mag_tau
  enddo ; enddo ; endif

end procedure SCM_CVMix_tests_wind_forcing
module procedure SCM_CVMix_tests_buoyancy_forcing
  integer :: i, j, is, ie, js, je, Isq, Ieq, Jsq, Jeq
  integer :: isd, ied, jsd, jed, IsdB, IedB, JsdB, JedB
  real :: PI  ! The ratio of the circumference of a circle to its diameter [nondim]
  PI = 4.0*atan(1.0)

  ! Bounds for loops and memory allocation
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec
  Isq = G%IscB ; Ieq = G%IecB ; Jsq = G%JscB ; Jeq = G%JecB
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  IsdB = G%IsdB ; IedB = G%IedB ; JsdB = G%JsdB ; JedB = G%JedB

  if (CS%UseHeatFlux) then
    ! Note CVMix test inputs give Heat flux in [Z C T-1 ~> m K s-1]
    ! therefore must convert to [Q R Z T-1 ~> W m-2] by multiplying
    ! by Rho0*Cp
    do J=Jsq,Jeq ; do i=is,ie
      fluxes%sens(i,J) = CS%surf_HF * CS%rho_restore * fluxes%C_p
    enddo ; enddo
  endif

  if (CS%UseEvaporation) then
    do J=Jsq,Jeq ; do i=is,ie
    ! Note CVMix test inputs give evaporation in [Z T-1 ~> m s-1]
    ! This therefore must be converted to mass flux in [R Z T-1 ~> kg m-2 s-1]
    ! by multiplying by density and some unit conversion factors.
      fluxes%evap(i,J) = CS%surf_evap * CS%rho_restore
    enddo ; enddo
  endif

  if (CS%UseDiurnalSW) then
    do J=Jsq,Jeq ; do i=is,ie
    ! Note CVMix test inputs give max sw rad in [Z C T-1 ~> m degC s-1]
    ! therefore must convert to [Q R Z T-1 ~> W m-2] by multiplying by Rho0*Cp
    ! Note diurnal cycle peaks at Noon.
      fluxes%sw(i,J) = CS%Max_sw *  max(0.0, cos(2*PI*(time_type_to_real(DAY)/86400.0 - 0.5))) * &
                       CS%rho_restore * fluxes%C_p
    enddo ; enddo
  endif

end procedure SCM_CVMix_tests_buoyancy_forcing
end submodule SCM_CVMix_tests_s
