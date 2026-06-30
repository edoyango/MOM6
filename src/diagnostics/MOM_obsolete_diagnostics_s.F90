submodule (MOM_obsolete_diagnostics) MOM_obsolete_diagnostics_s
#include <MOM_memory.h>
  implicit none
contains
module procedure register_obsolete_diagnostics
#include "version_variable.h"
  character(len=40)  :: mdl = "MOM_obsolete_diagnostics"  !< This module's name.
  logical :: foundEntry, causeFatal
  integer :: errType
  call log_version(param_file, mdl, version)
  call get_param(param_file, mdl, "OBSOLETE_DIAGNOSTIC_IS_FATAL", causeFatal,              &
                 "If an obsolete diagnostic variable appears in the diag_table, "//        &
                 "cause a FATAL error rather than issue a WARNING.", default=.true.)

  foundEntry = .false.
  ! Each obsolete entry, with replacement name is available.
  if (diag_found(diag, 'Net_Heat', 'net_heat_surface or net_heat_coupler')) foundEntry = .true.
  if (diag_found(diag, 'PmE', 'PRCmE'))                                     foundEntry = .true.
  if (diag_found(diag, 'froz_precip', 'fprec'))                             foundEntry = .true.
  if (diag_found(diag, 'liq_precip', 'lprec'))                              foundEntry = .true.
  if (diag_found(diag, 'virt_precip', 'vprec'))                             foundEntry = .true.
  if (diag_found(diag, 'froz_runoff', 'frunoff'))                           foundEntry = .true.
  if (diag_found(diag, 'liq_runoff', 'lrunoff'))                            foundEntry = .true.
  if (diag_found(diag, 'calving_heat_content', 'heat_content_frunoff'))     foundEntry = .true.
  if (diag_found(diag, 'precip_heat_content', 'heat_content_lprec'))        foundEntry = .true.
  if (diag_found(diag, 'evap_heat_content', 'heat_content_massout'))        foundEntry = .true.
  if (diag_found(diag, 'runoff_heat_content', 'heat_content_lrunoff'))      foundEntry = .true.
  if (diag_found(diag, 'latent_fprec'))                                     foundEntry = .true.
  if (diag_found(diag, 'latent_calve'))                                     foundEntry = .true.
  if (diag_found(diag, 'heat_rest', 'heat_restore'))                        foundEntry = .true.
  if (diag_found(diag, 'KPP_dTdt', 'KPP_NLT_dTdt'))                         foundEntry = .true.
  if (diag_found(diag, 'KPP_dSdt', 'KPP_NLT_dSdt'))                         foundEntry = .true.

  if (causeFatal) then ; errType = FATAL
  else ; errType = WARNING ; endif
  if (foundEntry .and. is_root_pe()) &
    call MOM_error(errType, 'MOM_obsolete_diagnostics: Obsolete diagnostics found in diag_table.')

end procedure register_obsolete_diagnostics
module procedure diag_found
  diag_found = found_in_diagtable(diag, varName)

  if (diag_found .and. is_root_pe()) then
    if (present(newVarName)) then
      call MOM_error(WARNING, 'MOM_obsolete_params: '//'diag_table entry "'// &
          trim(varName)//'" found. Use ''"'//trim(newVarName)//'" instead.' )
    else
      call MOM_error(WARNING, 'MOM_obsolete_params: '//'diag_table entry "'// &
          trim(varName)//'" is obsolete.' )
    endif
  endif

end procedure diag_found
end submodule MOM_obsolete_diagnostics_s
