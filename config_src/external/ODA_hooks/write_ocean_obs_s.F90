submodule (write_ocean_obs_mod) write_ocean_obs_mod_s
  implicit none
contains
module procedure open_profile_file
  open_profile_file=-1
end procedure open_profile_file
module procedure write_profile
  return
end procedure write_profile
module procedure close_profile_file
  return
end procedure close_profile_file
module procedure write_ocean_obs_init
  return
end procedure write_ocean_obs_init
end submodule write_ocean_obs_mod_s
