submodule (ocean_da_core_mod) ocean_da_core_mod_s
  implicit none
contains
  module procedure ocean_da_core_init
    Profiles=>NULL()
    return
  end procedure ocean_da_core_init
  module procedure get_profiles
    Profiles=>NULL()
    Current_Profiles=>NULL()

    return
  end procedure get_profiles
end submodule ocean_da_core_mod_s
