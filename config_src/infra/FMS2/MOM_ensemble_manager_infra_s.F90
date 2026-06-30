submodule (MOM_ensemble_manager_infra) MOM_ensemble_manager_infra_s
  implicit none
contains
module procedure ensemble_manager_init
  if (present(ensemble_suffix)) then
    call fms2_io_set_filename_appendix(trim(ensemble_suffix))
  else
    call FMS_ensemble_manager_init()
  endif

end procedure ensemble_manager_init
module procedure ensemble_pelist_setup
  call FMS_ensemble_pelist_setup(concurrent, atmos_npes, ocean_npes, land_npes, ice_npes, &
         Atm_pelist, Ocean_pelist, Land_pelist, Ice_pelist)

end procedure ensemble_pelist_setup
module procedure get_ensemble_id
  get_ensemble_id = FMS_get_ensemble_id()

end procedure get_ensemble_id
module procedure get_ensemble_size
  get_ensemble_size = FMS_get_ensemble_size()

end procedure get_ensemble_size
module procedure get_ensemble_pelist
  call FMS_get_ensemble_pelist(pelist, name)

end procedure get_ensemble_pelist
module procedure get_ensemble_filter_pelist
  call FMS_get_Ensemble_filter_pelist(pelist, name)

end procedure get_ensemble_filter_pelist
end submodule MOM_ensemble_manager_infra_s
