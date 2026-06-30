submodule (MOM_ice_shelf_state) MOM_ice_shelf_state_s
  implicit none
contains
module procedure ice_shelf_state_init
  integer :: isd, ied, jsd, jed
  isd = G%isd ; jsd = G%jsd ; ied = G%ied ; jed = G%jed

  if (associated(ISS)) then
    call MOM_error(FATAL, "MOM_ice_shelf_state.F90, ice_shelf_state_init: "// &
                          "called with an associated ice_shelf_state pointer.")
    return
  endif
  allocate(ISS)

  allocate(ISS%mass_shelf(isd:ied,jsd:jed), source=0.0 )
  allocate(ISS%area_shelf_h(isd:ied,jsd:jed), source=0.0 )
  allocate(ISS%melt_mask(isd:ied,jsd:jed), source=1.0 )
  allocate(ISS%h_shelf(isd:ied,jsd:jed), source=0.0 )
  allocate(ISS%dhdt_shelf(isd:ied,jsd:jed), source=0.0 )
  allocate(ISS%hmask(isd:ied,jsd:jed), source=-2.0 )

  allocate(ISS%tflux_ocn(isd:ied,jsd:jed), source=0.0 )
  allocate(ISS%water_flux(isd:ied,jsd:jed), source=0.0 )
  allocate(ISS%salt_flux(isd:ied,jsd:jed), source=0.0 )
  allocate(ISS%tflux_shelf(isd:ied,jsd:jed), source=0.0 )
  allocate(ISS%tfreeze(isd:ied,jsd:jed), source=0.0 )

  allocate(ISS%frazil(isd:ied,jsd:jed), source=0.0 )
  allocate(ISS%calving(isd:ied,jsd:jed), source=0.0 )
  allocate(ISS%calving_hflx(isd:ied,jsd:jed), source=0.0 )
end procedure ice_shelf_state_init
module procedure ice_shelf_state_end
  if (.not.associated(ISS)) return

  deallocate(ISS%mass_shelf, ISS%area_shelf_h, ISS%h_shelf, ISS%dhdt_shelf, ISS%hmask)

  deallocate(ISS%tflux_ocn, ISS%water_flux, ISS%salt_flux, ISS%tflux_shelf)
  deallocate(ISS%tfreeze, ISS%frazil)

  deallocate(ISS%calving, ISS%calving_hflx)

  deallocate(ISS)

end procedure ice_shelf_state_end
end submodule MOM_ice_shelf_state_s
