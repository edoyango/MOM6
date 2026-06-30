submodule (atmos_ocean_fluxes_mod) atmos_ocean_fluxes_mod_s
  implicit none
contains
module procedure aof_set_coupler_flux
  coupler_index = -1

end procedure aof_set_coupler_flux
end submodule atmos_ocean_fluxes_mod_s
