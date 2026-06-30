submodule (MOM_PressureForce) MOM_PressureForce_s
#include <MOM_memory.h>
  implicit none
contains
module procedure PressureForce
  if (CS%Analytic_FV_PGF) then
    if (GV%Boussinesq) then
      call PressureForce_FV_Bouss(h, tv, PFu, PFv, G, GV, US, CS%PressureForce_FV, &
                                   ALE_CSp, ADp, p_atm, pbce, eta)
    else
      call PressureForce_FV_nonBouss(h, tv, PFu, PFv, G, GV, US, CS%PressureForce_FV, &
                                      ALE_CSp, ADp, p_atm, pbce, eta)
    endif
  else
    if (GV%Boussinesq) then
      call PressureForce_Mont_Bouss(h, tv, PFu, PFv, G, GV, US, CS%PressureForce_Mont, &
                                    p_atm, pbce, eta)
    else
      call PressureForce_Mont_nonBouss(h, tv, PFu, PFv, G, GV, US, CS%PressureForce_Mont, &
                                       p_atm, pbce, eta)
    endif
  endif

end procedure PressureForce
module procedure PressureForce_init
#include "version_variable.h"
  character(len=40)  :: mdl = "MOM_PressureForce" ! This module's name.
  call log_version(param_file, mdl, version, "")
  call get_param(param_file, mdl, "ANALYTIC_FV_PGF", CS%Analytic_FV_PGF, &
                 "If true the pressure gradient forces are calculated "//&
                 "with a finite volume form that analytically integrates "//&
                 "the equations of state in pressure to avoid any "//&
                 "possibility of numerical thermobaric instability, as "//&
                 "described in Adcroft et al., O. Mod. (2008).", default=.true.)

  if (CS%Analytic_FV_PGF) then
    call PressureForce_FV_init(Time, G, GV, US, param_file, diag, &
             CS%PressureForce_FV, ADp, SAL_CSp, tides_CSp)
  else
    call PressureForce_Mont_init(Time, G, GV, US, param_file, diag, &
             CS%PressureForce_Mont, SAL_CSp, tides_CSp)
  endif
end procedure PressureForce_init
end submodule MOM_PressureForce_s
