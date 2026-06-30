submodule (MOM_unit_scaling) MOM_unit_scaling_s
  implicit none
contains
module procedure unit_scaling_init
  integer :: Z_power, L_power, T_power, R_power, Q_power, C_power, S_power
  real    :: Z_rescale_factor, L_rescale_factor, T_rescale_factor, R_rescale_factor, Q_rescale_factor
  real    :: C_rescale_factor, S_rescale_factor
# include "version_variable.h"
  character(len=16) :: mdl = "MOM_unit_scaling"
  if (associated(US)) call MOM_error(FATAL, &
     'unit_scaling_init: called with an associated US pointer.')
  allocate(US)

  ! Read all relevant parameters and write them to the model log.
  call log_version(param_file, mdl, version, &
               "Parameters for doing unit scaling of variables.", debugging=.true.)
  call get_param(param_file, mdl, "Z_RESCALE_POWER", Z_power, &
               "An integer power of 2 that is used to rescale the model's "//&
               "internal units of depths and heights.  Valid values range from -300 to 300.", &
               default=0, debuggingParam=.true.)
  call get_param(param_file, mdl, "L_RESCALE_POWER", L_power, &
               "An integer power of 2 that is used to rescale the model's "//&
               "internal units of lateral distances.  Valid values range from -300 to 300.", &
               default=0, debuggingParam=.true.)
  call get_param(param_file, mdl, "T_RESCALE_POWER", T_power, &
               "An integer power of 2 that is used to rescale the model's "//&
               "internal units of time.  Valid values range from -300 to 300.", &
               default=0, debuggingParam=.true.)
  call get_param(param_file, mdl, "R_RESCALE_POWER", R_power, &
               "An integer power of 2 that is used to rescale the model's "//&
               "internal units of density.  Valid values range from -300 to 300.", &
               default=0, debuggingParam=.true.)
  call get_param(param_file, mdl, "Q_RESCALE_POWER", Q_power, &
               "An integer power of 2 that is used to rescale the model's "//&
               "internal units of heat content.  Valid values range from -300 to 300.", &
               default=0, debuggingParam=.true.)
  call get_param(param_file, mdl, "C_RESCALE_POWER", C_power, &
               "An integer power of 2 that is used to rescale the model's "//&
               "internal units of temperature.  Valid values range from -300 to 300.", &
               default=0, debuggingParam=.true.)
  call get_param(param_file, mdl, "S_RESCALE_POWER", S_power, &
               "An integer power of 2 that is used to rescale the model's "//&
               "internal units of salinity.  Valid values range from -300 to 300.", &
               default=0, debuggingParam=.true.)

  if (abs(Z_power) > 300) call MOM_error(FATAL, "unit_scaling_init: "//&
                 "Z_RESCALE_POWER is outside of the valid range of -300 to 300.")
  if (abs(L_power) > 300) call MOM_error(FATAL, "unit_scaling_init: "//&
                 "L_RESCALE_POWER is outside of the valid range of -300 to 300.")
  if (abs(T_power) > 300) call MOM_error(FATAL, "unit_scaling_init: "//&
                 "T_RESCALE_POWER is outside of the valid range of -300 to 300.")
  if (abs(R_power) > 300) call MOM_error(FATAL, "unit_scaling_init: "//&
                 "R_RESCALE_POWER is outside of the valid range of -300 to 300.")
  if (abs(Q_power) > 300) call MOM_error(FATAL, "unit_scaling_init: "//&
                 "Q_RESCALE_POWER is outside of the valid range of -300 to 300.")
  if (abs(C_power) > 300) call MOM_error(FATAL, "unit_scaling_init: "//&
                 "C_RESCALE_POWER is outside of the valid range of -300 to 300.")
  if (abs(S_power) > 300) call MOM_error(FATAL, "unit_scaling_init: "//&
                 "S_RESCALE_POWER is outside of the valid range of -300 to 300.")

  Z_rescale_factor = 1.0
  if (Z_power /= 0) Z_rescale_factor = 2.0**Z_power
  US%Z_to_m = 1.0 * Z_rescale_factor
  US%m_to_Z = 1.0 / Z_rescale_factor

  L_rescale_factor = 1.0
  if (L_power /= 0) L_rescale_factor = 2.0**L_power
  US%L_to_m = 1.0 * L_rescale_factor
  US%m_to_L = 1.0 / L_rescale_factor

  T_rescale_factor = 1.0
  if (T_power /= 0) T_rescale_factor = 2.0**T_power
  US%T_to_s = 1.0 * T_rescale_factor
  US%s_to_T = 1.0 / T_rescale_factor

  R_rescale_factor = 1.0
  if (R_power /= 0) R_rescale_factor = 2.0**R_power
  US%R_to_kg_m3 = 1.0 * R_rescale_factor
  US%kg_m3_to_R = 1.0 / R_rescale_factor

  Q_Rescale_factor = 1.0
  if (Q_power /= 0) Q_Rescale_factor = 2.0**Q_power
  US%Q_to_J_kg = 1.0 * Q_Rescale_factor
  US%J_kg_to_Q = 1.0 / Q_Rescale_factor

  C_Rescale_factor = 1.0
  if (C_power /= 0) C_Rescale_factor = 2.0**C_power
  US%C_to_degC = 1.0 * C_Rescale_factor
  US%degC_to_C = 1.0 / C_Rescale_factor

  S_Rescale_factor = 1.0
  if (S_power /= 0) S_Rescale_factor = 2.0**S_power
  US%S_to_ppt = 1.0 * S_Rescale_factor
  US%ppt_to_S = 1.0 / S_Rescale_factor

  call set_unit_scaling_combos(US)
end procedure unit_scaling_init
module procedure unit_no_scaling_init
  if (associated(US)) call MOM_error(FATAL, &
     'unit_scaling_init: called with an associated US pointer.')
  allocate(US)

  US%Z_to_m = 1.0 ; US%m_to_Z = 1.0
  US%L_to_m = 1.0 ; US%m_to_L = 1.0
  US%T_to_s = 1.0 ; US%s_to_T = 1.0
  US%R_to_kg_m3 = 1.0 ; US%kg_m3_to_R = 1.0
  US%Q_to_J_kg = 1.0 ; US%J_kg_to_Q = 1.0
  US%C_to_degC = 1.0 ; US%degC_to_C = 1.0
  US%S_to_ppt = 1.0 ; US%ppt_to_S = 1.0

  call set_unit_scaling_combos(US)
end procedure unit_no_scaling_init
module procedure set_unit_scaling_combos
  US%Z_to_L = US%Z_to_m * US%m_to_L
  US%L_to_Z = US%L_to_m * US%m_to_Z
  ! Horizontal velocities:
  US%L_T_to_m_s = US%L_to_m * US%s_to_T
  US%m_s_to_L_T = US%m_to_L * US%T_to_s
  ! Horizontal accelerations:
  US%L_T2_to_m_s2 = US%L_to_m * US%s_to_T**2
    ! It does not look like US%m_s2_to_L_T2 would be used, so it does not exist.
  ! Vertical diffusivities and viscosities:
  US%Z2_T_to_m2_s = US%Z_to_m**2 * US%s_to_T
  US%m2_s_to_Z2_T = US%m_to_Z**2 * US%T_to_s
  ! Column mass loads:
  US%RZ_to_kg_m2  = US%R_to_kg_m3 * US%Z_to_m
    ! It does not seem like US%kg_m2_to_RZ would be used enough in MOM6 to justify its existence.
  ! Vertical mass fluxes:
  US%kg_m2s_to_RZ_T = US%kg_m3_to_R * US%m_to_Z * US%T_to_s
  US%RZ_T_to_kg_m2s = US%R_to_kg_m3 * US%Z_to_m * US%s_to_T
  ! Turbulent kinetic energy vertical fluxes:
  US%RZ3_T3_to_W_m2 = US%R_to_kg_m3 * US%Z_to_m**3 * US%s_to_T**3
  US%W_m2_to_RZ3_T3 = US%kg_m3_to_R * US%m_to_Z**3 * US%T_to_s**3
  ! Vertical heat fluxes:
  US%W_m2_to_QRZ_T = US%J_kg_to_Q * US%kg_m3_to_R * US%m_to_Z * US%T_to_s
  US%QRZ_T_to_W_m2 = US%Q_to_J_kg * US%R_to_kg_m3 * US%Z_to_m * US%s_to_T
  ! Pressures:
  US%RL2_T2_to_Pa = US%R_to_kg_m3 * US%L_T_to_m_s**2
  US%Pa_to_RL2_T2 = US%kg_m3_to_R * US%m_s_to_L_T**2
  ! Wind stresses:
  US%RLZ_T2_to_Pa = US%R_to_kg_m3 * US%L_T_to_m_s**2 * US%Z_to_L
  US%Pa_to_RLZ_T2 = US%kg_m3_to_R * US%m_s_to_L_T**2 * US%L_to_Z
  ! Masses:
  US%RZL2_to_kg = US%R_to_kg_m3 * US%Z_to_m * US%L_to_m**2

end procedure set_unit_scaling_combos
module procedure fix_restart_unit_scaling
  US%m_to_Z_restart = 1.0 ! US%m_to_Z
  US%m_to_L_restart = 1.0 ! US%m_to_L
  US%s_to_T_restart = 1.0 ! US%s_to_T
  US%kg_m3_to_R_restart = 1.0 ! US%kg_m3_to_R
  US%J_kg_to_Q_restart = 1.0 ! US%J_kg_to_Q

  if (present(unscaled)) then ; if (unscaled) then
    US%m_to_Z_restart = 1.0
    US%m_to_L_restart = 1.0
    US%s_to_T_restart = 1.0
    US%kg_m3_to_R_restart = 1.0
    US%J_kg_to_Q_restart = 1.0
  endif ; endif

end procedure fix_restart_unit_scaling
module procedure unit_scaling_end
  deallocate( US )

end procedure unit_scaling_end
end submodule MOM_unit_scaling_s
