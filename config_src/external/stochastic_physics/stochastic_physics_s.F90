submodule (stochastic_physics) stochastic_physics_s
  implicit none
contains
module procedure init_stochastic_physics_ocn
  iret=0
  if (pert_epbl_in) then
    call MOM_error(WARNING, 'init_stochastic_physics_ocn: pert_epbl needs to be false if using the stub')
    iret=-1
  endif
  if (do_sppt_in) then
    call MOM_error(WARNING, 'init_stochastic_physics_ocn: do_sppt needs to be false if using the stub')
    iret=-1
  endif
  if (do_skeb_in) then
    call MOM_error(WARNING, 'init_stochastic_physics_ocn: do_skeb needs to be false if using the stub')
    iret=-1
  endif

  ! This stub function does not actually do anything.
  return
end procedure init_stochastic_physics_ocn
module procedure run_stochastic_physics_ocn
  return
end procedure run_stochastic_physics_ocn
end submodule stochastic_physics_s
