submodule (MOM_variables) MOM_variables_s
#include <MOM_memory.h>
  implicit none
contains
module procedure allocate_surface_state
  logical :: use_temp, alloc_integ, use_melt_potential, alloc_iceshelves, alloc_frazil, alloc_fco2
  logical :: even_turns  ! True if turns is absent or even
  integer :: tr_field_i_mem(4), tr_field_j_mem(4)
  integer :: is, ie, js, je, isd, ied, jsd, jed
  integer :: isdB, iedB, jsdB, jedB
  is  = G%isc ; ie  = G%iec ; js  = G%jsc ; je  = G%jec
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  isdB = G%isdB ; iedB = G%iedB ; jsdB = G%jsdB ; jedB = G%jedB

  use_temp = .true. ; if (present(use_temperature)) use_temp = use_temperature
  alloc_integ = .true. ; if (present(do_integrals)) alloc_integ = do_integrals
  use_melt_potential = .false. ; if (present(use_meltpot)) use_melt_potential = use_meltpot
  alloc_iceshelves = .false. ; if (present(use_iceshelves)) alloc_iceshelves = use_iceshelves
  alloc_frazil = .true. ; if (present(omit_frazil)) alloc_frazil = .not.omit_frazil
  alloc_fco2 = .false. ; if (present(use_marbl_tracers)) alloc_fco2 = use_marbl_tracers

  if (sfc_state%arrays_allocated) return

  if (use_temp) then
    allocate(sfc_state%SST(isd:ied,jsd:jed), source=0.0)
    allocate(sfc_state%SSS(isd:ied,jsd:jed), source=0.0)
  else
    allocate(sfc_state%sfc_density(isd:ied,jsd:jed), source=0.0)
  endif
  if (use_temp .and. alloc_frazil) then
    allocate(sfc_state%frazil(isd:ied,jsd:jed), source=0.0)
  endif
  allocate(sfc_state%sea_lev(isd:ied,jsd:jed), source=0.0)
  allocate(sfc_state%Hml(isd:ied,jsd:jed), source=0.0)
  allocate(sfc_state%u(IsdB:IedB,jsd:jed), source=0.0)
  allocate(sfc_state%v(isd:ied,JsdB:JedB), source=0.0)

  if (use_melt_potential) then
    allocate(sfc_state%melt_potential(isd:ied,jsd:jed), source=0.0)
  endif

  if (alloc_integ) then
    ! Allocate structures for the vertically integrated ocean_mass, ocean_heat, and ocean_salt.
    allocate(sfc_state%ocean_mass(isd:ied,jsd:jed), source=0.0)
    if (use_temp) then
      allocate(sfc_state%ocean_heat(isd:ied,jsd:jed), source=0.0)
      allocate(sfc_state%ocean_salt(isd:ied,jsd:jed), source=0.0)
    endif
  endif

  if (alloc_iceshelves) then
    allocate(sfc_state%taux_shelf(IsdB:IedB,jsd:jed), source=0.0)
    allocate(sfc_state%tauy_shelf(isd:ied,JsdB:JedB), source=0.0)
  endif

  ! The data fields in the coupler_2d_bc_type are never rotated.
  even_turns = .true. ; if (present(turns)) even_turns = (modulo(turns, 2) == 0)
  if (even_turns) then
    tr_field_i_mem(1:4) = (/is,is,ie,ie/) ; tr_field_j_mem(1:4) = (/js,js,je,je/)
  else
    tr_field_i_mem(1:4) = (/js,js,je,je/) ; tr_field_j_mem(1:4) = (/is,is,ie,ie/)
  endif
  if (present(gas_fields_ocn)) then
    call coupler_type_spawn(gas_fields_ocn, sfc_state%tr_fields, &
                            tr_field_i_mem, tr_field_j_mem, as_needed=.true.)
  elseif (present(sfc_state_in)) then
    if (coupler_type_initialized(sfc_state_in%tr_fields)) then
      call coupler_type_spawn(sfc_state_in%tr_fields, sfc_state%tr_fields, &
                              tr_field_i_mem, tr_field_j_mem, as_needed=.true.)
    endif
  endif

  if (alloc_fco2) then
    allocate(sfc_state%fco2(isd:ied,jsd:jed), source=0.0)
  endif

  sfc_state%arrays_allocated = .true.

end procedure allocate_surface_state
module procedure deallocate_surface_state
  if (.not.sfc_state%arrays_allocated) return

  if (allocated(sfc_state%melt_potential)) deallocate(sfc_state%melt_potential)
  if (allocated(sfc_state%SST)) deallocate(sfc_state%SST)
  if (allocated(sfc_state%SSS)) deallocate(sfc_state%SSS)
  if (allocated(sfc_state%sfc_density)) deallocate(sfc_state%sfc_density)
  if (allocated(sfc_state%sea_lev)) deallocate(sfc_state%sea_lev)
  if (allocated(sfc_state%Hml)) deallocate(sfc_state%Hml)
  if (allocated(sfc_state%u)) deallocate(sfc_state%u)
  if (allocated(sfc_state%v)) deallocate(sfc_state%v)
  if (allocated(sfc_state%ocean_mass)) deallocate(sfc_state%ocean_mass)
  if (allocated(sfc_state%ocean_heat)) deallocate(sfc_state%ocean_heat)
  if (allocated(sfc_state%ocean_salt)) deallocate(sfc_state%ocean_salt)
  if (allocated(sfc_state%fco2)) deallocate(sfc_state%fco2)
  call coupler_type_destructor(sfc_state%tr_fields)

  sfc_state%arrays_allocated = .false.

end procedure deallocate_surface_state
module procedure rotate_surface_state
  logical :: use_temperature, do_integrals, use_melt_potential, use_iceshelves
  use_temperature = allocated(sfc_state_in%SST) &
      .and. allocated(sfc_state_in%SSS)
  use_melt_potential = allocated(sfc_state_in%melt_potential)
  do_integrals = allocated(sfc_state_in%ocean_mass)
  use_iceshelves = allocated(sfc_state_in%taux_shelf) &
      .and. allocated(sfc_state_in%tauy_shelf)

  if (.not. sfc_state%arrays_allocated) then
    call allocate_surface_state(sfc_state, G, use_temperature=use_temperature, &
            do_integrals=do_integrals, use_meltpot=use_melt_potential, &
            use_iceshelves=use_iceshelves, sfc_state_in=sfc_state_in, turns=turns)
  endif

  if (use_temperature) then
    call rotate_array(sfc_state_in%SST, turns, sfc_state%SST)
    call rotate_array(sfc_state_in%SSS, turns, sfc_state%SSS)
  else
    call rotate_array(sfc_state_in%sfc_density, turns, sfc_state%sfc_density)
  endif

  call rotate_array(sfc_state_in%Hml, turns, sfc_state%Hml)
  call rotate_vector(sfc_state_in%u, sfc_state_in%v, turns, &
      sfc_state%u, sfc_state%v)
  call rotate_array(sfc_state_in%sea_lev, turns, sfc_state%sea_lev)

  if (use_melt_potential) then
    call rotate_array(sfc_state_in%melt_potential, turns, sfc_state%melt_potential)
  endif

  if (do_integrals) then
    call rotate_array(sfc_state_in%ocean_mass, turns, sfc_state%ocean_mass)
    if (use_temperature) then
      call rotate_array(sfc_state_in%ocean_heat, turns, sfc_state%ocean_heat)
      call rotate_array(sfc_state_in%ocean_salt, turns, sfc_state%ocean_salt)
      call rotate_array(sfc_state_in%SSS, turns, sfc_state%SSS)
    endif
  endif

  if (use_iceshelves) then
    call rotate_vector(sfc_state_in%taux_shelf, sfc_state_in%tauy_shelf, turns, &
        sfc_state%taux_shelf, sfc_state%tauy_shelf)
  endif

  if (use_temperature .and. allocated(sfc_state_in%frazil)) &
    call rotate_array(sfc_state_in%frazil, turns, sfc_state%frazil)

  ! Scalar transfers
  sfc_state%T_is_conT = sfc_state_in%T_is_conT
  sfc_state%S_is_absS = sfc_state_in%S_is_absS

  ! NOTE: Tracer fields are handled by FMS, so are left unrotated.  Any
  ! reads/writes to tr_fields must be appropriately rotated.
  if (coupler_type_initialized(sfc_state_in%tr_fields)) then
    call coupler_type_copy_data(sfc_state_in%tr_fields, sfc_state%tr_fields)
  endif
end procedure rotate_surface_state
module procedure alloc_BT_cont_type
  integer :: isd, ied, jsd, jed, IsdB, IedB, JsdB, JedB, nz
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed ; nz = GV%ke
  IsdB = G%IsdB ; IedB = G%IedB ; JsdB = G%JsdB ; JedB = G%JedB

  if (associated(BT_cont)) call MOM_error(FATAL, &
    "alloc_BT_cont_type called with an associated BT_cont_type pointer.")

  allocate(BT_cont)
  allocate(BT_cont%FA_u_WW(IsdB:IedB,jsd:jed), source=0.0)
  allocate(BT_cont%FA_u_W0(IsdB:IedB,jsd:jed), source=0.0)
  allocate(BT_cont%FA_u_E0(IsdB:IedB,jsd:jed), source=0.0)
  allocate(BT_cont%FA_u_EE(IsdB:IedB,jsd:jed), source=0.0)
  allocate(BT_cont%uBT_WW(IsdB:IedB,jsd:jed), source=0.0)
  allocate(BT_cont%uBT_EE(IsdB:IedB,jsd:jed), source=0.0)

  allocate(BT_cont%FA_v_SS(isd:ied,JsdB:JedB), source=0.0)
  allocate(BT_cont%FA_v_S0(isd:ied,JsdB:JedB), source=0.0)
  allocate(BT_cont%FA_v_N0(isd:ied,JsdB:JedB), source=0.0)
  allocate(BT_cont%FA_v_NN(isd:ied,JsdB:JedB), source=0.0)
  allocate(BT_cont%vBT_SS(isd:ied,JsdB:JedB), source=0.0)
  allocate(BT_cont%vBT_NN(isd:ied,JsdB:JedB), source=0.0)

  if (present(alloc_faces)) then ; if (alloc_faces) then
    allocate(BT_cont%h_u(IsdB:IedB,jsd:jed,1:nz), source=0.0)
    allocate(BT_cont%h_v(isd:ied,JsdB:JedB,1:nz), source=0.0)
  endif ; endif

end procedure alloc_BT_cont_type
module procedure dealloc_BT_cont_type
  if (.not.associated(BT_cont)) return

  deallocate(BT_cont%FA_u_WW) ; deallocate(BT_cont%FA_u_W0)
  deallocate(BT_cont%FA_u_E0) ; deallocate(BT_cont%FA_u_EE)
  deallocate(BT_cont%uBT_WW)  ; deallocate(BT_cont%uBT_EE)

  deallocate(BT_cont%FA_v_SS) ; deallocate(BT_cont%FA_v_S0)
  deallocate(BT_cont%FA_v_N0) ; deallocate(BT_cont%FA_v_NN)
  deallocate(BT_cont%vBT_SS)  ; deallocate(BT_cont%vBT_NN)

  if (allocated(BT_cont%h_u)) deallocate(BT_cont%h_u)
  if (allocated(BT_cont%h_v)) deallocate(BT_cont%h_v)

  deallocate(BT_cont)

end procedure dealloc_BT_cont_type
module procedure MOM_thermovar_chksum
  if (associated(tv%T)) &
    call hchksum(tv%T, mesg//" tv%T", G%HI, unscale=US%C_to_degC)
  if (associated(tv%S)) &
    call hchksum(tv%S, mesg//" tv%S", G%HI, unscale=US%S_to_ppt)
  if (associated(tv%frazil)) &
    call hchksum(tv%frazil, mesg//" tv%frazil", G%HI, unscale=US%Q_to_J_kg*US%RZ_to_kg_m2)
  if (associated(tv%salt_deficit)) &
    call hchksum(tv%salt_deficit, mesg//" tv%salt_deficit", G%HI, unscale=US%RZ_to_kg_m2*US%S_to_ppt)
  if (associated(tv%TempxPmE)) &
    call hchksum(tv%TempxPmE, mesg//" tv%TempxPmE", G%HI, unscale=US%RZ_to_kg_m2*US%C_to_degC)
end procedure MOM_thermovar_chksum
end submodule MOM_variables_s
