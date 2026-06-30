submodule (user_revise_forcing) user_revise_forcing_s
  implicit none
contains
module procedure user_alter_forcing
  return

end procedure user_alter_forcing
module procedure user_revise_forcing_init
# include "version_variable.h"
  character(len=40) :: mdl = "user_revise_forcing" !< This module's name.
  call log_version(param_file, mdl, version)

end procedure user_revise_forcing_init
end submodule user_revise_forcing_s
