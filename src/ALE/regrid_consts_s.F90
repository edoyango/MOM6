submodule (regrid_consts) regrid_consts_s
  implicit none
contains
module procedure coordinateMode
  select case ( uppercase(trim(string)) )
    case (trim(REGRIDDING_LAYER_STRING)); coordinateMode = REGRIDDING_LAYER
    case (trim(REGRIDDING_ZSTAR_STRING)); coordinateMode = REGRIDDING_ZSTAR
    case (trim(REGRIDDING_ZSTAR_STRING_OLD)); coordinateMode = REGRIDDING_ZSTAR
    case (trim(REGRIDDING_RHO_STRING));   coordinateMode = REGRIDDING_RHO
    case (trim(REGRIDDING_SIGMA_STRING)); coordinateMode = REGRIDDING_SIGMA
    case (trim(REGRIDDING_HYCOM1_STRING)); coordinateMode = REGRIDDING_HYCOM1
    case (trim(REGRIDDING_HYBGEN_STRING)); coordinateMode = REGRIDDING_HYBGEN
    case (trim(REGRIDDING_ARBITRARY_STRING)); coordinateMode = REGRIDDING_ARBITRARY
    case (trim(REGRIDDING_SIGMA_SHELF_ZSTAR_STRING)); coordinateMode = REGRIDDING_SIGMA_SHELF_ZSTAR
    case (trim(REGRIDDING_ADAPTIVE_STRING)); coordinateMode = REGRIDDING_ADAPTIVE
    case default ; call MOM_error(FATAL, "coordinateMode: "//&
       "Unrecognized choice of coordinate ("//trim(string)//").")
  end select
end procedure coordinateMode
module procedure coordinateUnitsI
  select case ( coordMode )
    case (REGRIDDING_LAYER); coordinateUnitsI = "kg m^-3"
    case (REGRIDDING_ZSTAR); coordinateUnitsI = "m"
    case (REGRIDDING_SIGMA_SHELF_ZSTAR); coordinateUnitsI = "m"
    case (REGRIDDING_RHO);   coordinateUnitsI = "kg m^-3"
    case (REGRIDDING_SIGMA); coordinateUnitsI = "Non-dimensional"
    case (REGRIDDING_HYCOM1); coordinateUnitsI = "m"
    case (REGRIDDING_HYBGEN); coordinateUnitsI = "m"
    case (REGRIDDING_ADAPTIVE); coordinateUnitsI = "m"
    case default ; call MOM_error(FATAL, "coordinateUnts: "//&
       "Unrecognized coordinate mode.")
  end select
end procedure coordinateUnitsI
module procedure coordinateUnitsS
  integer :: coordMode
  coordMode = coordinateMode(string)
  coordinateUnitsS = coordinateUnitsI(coordMode)
end procedure coordinateUnitsS
module procedure state_dependent_char
  state_dependent_char = state_dependent_int( coordinateMode(string) )

end procedure state_dependent_char
module procedure state_dependent_int
  select case ( mode )
    case (REGRIDDING_LAYER); state_dependent_int = .true.
    case (REGRIDDING_ZSTAR); state_dependent_int = .false.
    case (REGRIDDING_SIGMA_SHELF_ZSTAR); state_dependent_int = .false.
    case (REGRIDDING_RHO);   state_dependent_int = .true.
    case (REGRIDDING_SIGMA); state_dependent_int = .false.
    case (REGRIDDING_HYCOM1); state_dependent_int = .true.
    case (REGRIDDING_HYBGEN); state_dependent_int = .true.
    case (REGRIDDING_ADAPTIVE); state_dependent_int = .true.
    case default ; call MOM_error(FATAL, "state_dependent: "//&
       "Unrecognized choice of coordinate.")
  end select
end procedure state_dependent_int
end submodule regrid_consts_s
