submodule (MOM_database_comms) MOM_database_comms_s
  implicit none
contains
module procedure database_comms_init
  call MOM_error(WARNING,"dbcomms_init was compiled using the dummy module. If this was\n"//&
                       "a mistake, please follow the instructions in:\n"//&
                       "MOM6/config_src/external/dbclient/README.md")
end procedure database_comms_init
end submodule MOM_database_comms_s
