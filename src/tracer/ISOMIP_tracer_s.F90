submodule (ISOMIP_tracer) ISOMIP_tracer_s
#include <MOM_memory.h>
  implicit none
contains
module procedure register_ISOMIP_tracer
  character(len=80)  :: name, longname
# include "version_variable.h"
  character(len=40)  :: mdl = "ISOMIP_tracer" ! This module's name.
  character(len=200) :: inputdir
  character(len=48)  :: flux_units ! The units for tracer fluxes, usually
  real, pointer :: tr_ptr(:,:,:) => NULL() ! The tracer concentration [conc]
  integer :: isd, ied, jsd, jed, nz, m
  isd = HI%isd ; ied = HI%ied ; jsd = HI%jsd ; jed = HI%jed ; nz = GV%ke

  if (associated(CS)) then
    call MOM_error(FATAL, "ISOMIP_register_tracer called with an "// &
                            "associated control structure.")
  endif
  allocate(CS)

  ! Read all relevant parameters and write them to the model log.
  call log_version(param_file, mdl, version, "")
  call get_param(param_file, mdl, "ISOMIP_TRACER_IC_FILE", CS%tracer_IC_file, &
                 "The name of a file from which to read the initial "//&
                 "conditions for the ISOMIP tracers, or blank to initialize "//&
                 "them internally.", default=" ")
  if (len_trim(CS%tracer_IC_file) >= 1) then
    call get_param(param_file, mdl, "INPUTDIR", inputdir, default=".")
    inputdir = slasher(inputdir)
    CS%tracer_IC_file = trim(inputdir)//trim(CS%tracer_IC_file)
    call log_param(param_file, mdl, "INPUTDIR/ISOMIP_TRACER_IC_FILE", &
                   CS%tracer_IC_file)
  endif
  call get_param(param_file, mdl, "SPONGE", CS%use_sponge, &
                 "If true, sponges may be applied anywhere in the domain. "//&
                 "The exact location and properties of those sponges are "//&
                 "specified from MOM_initialization.F90.", default=.false.)

  allocate(CS%tr(isd:ied,jsd:jed,nz,NTR), source=0.0)

  do m=1,NTR
    write(name,'("tr_D",I0)') m
    write(longname,'("Concentration of ISOMIP Tracer ",I2.2)') m
    CS%tr_desc(m) = var_desc(name, units="kg kg-1", longname=longname, caller=mdl)
    if (GV%Boussinesq) then ; flux_units = "kg kg-1 m3 s-1"
    else ; flux_units = "kg s-1" ; endif

    ! This is needed to force the compiler not to do a copy in the registration
    ! calls.  Curses on the designers and implementers of Fortran90.
    tr_ptr => CS%tr(:,:,:,m)
    ! Register the tracer for horizontal advection, diffusion, and restarts.
    call register_tracer(tr_ptr, tr_Reg, param_file, HI, GV, &
                         name=name, longname=longname, units="kg kg-1", &
                         registry_diags=.true., flux_units=flux_units, &
                         restart_CS=restart_CS)

    !   Set coupled_tracers to be true (hard-coded above) to provide the surface
    ! values to the coupler (if any).  This is meta-code and its arguments will
    ! currently (deliberately) give fatal errors if it is used.
    if (CS%coupled_tracers) &
      CS%ind_tr(m) = atmos_ocn_coupler_flux(trim(name)//'_flux', &
          flux_type=' ', implementation=' ', caller="register_ISOMIP_tracer")
  enddo

  CS%tr_Reg => tr_Reg
  register_ISOMIP_tracer = .true.
end procedure register_ISOMIP_tracer
module procedure initialize_ISOMIP_tracer
  character(len=16) :: name     ! A variable's name in a NetCDF file.
  real :: h_neglect         ! A thickness that is so small it is usually lost
  integer :: i, j, k, is, ie, js, je, isd, ied, jsd, jed, nz, m
  integer :: IsdB, IedB, JsdB, JedB
  if (.not.associated(CS)) return
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  IsdB = G%IsdB ; IedB = G%IedB ; JsdB = G%JsdB ; JedB = G%JedB
  h_neglect = GV%H_subroundoff

  CS%Time => day
  CS%diag => diag

  if (.not.restart) then
    if (len_trim(CS%tracer_IC_file) >= 1) then
      !  Read the tracer concentrations from a netcdf file.
      if (.not.file_exists(CS%tracer_IC_file, G%Domain)) &
        call MOM_error(FATAL, "ISOMIP_initialize_tracer: Unable to open "// &
                        CS%tracer_IC_file)
      do m=1,NTR
        call query_vardesc(CS%tr_desc(m), name, caller="initialize_ISOMIP_tracer")
        call MOM_read_data(CS%tracer_IC_file, trim(name), CS%tr(:,:,:,m), G%Domain)
      enddo
    else
      do m=1,NTR
        do k=1,nz ; do j=js,je ; do i=is,ie
          CS%tr(i,j,k,m) = 0.0
        enddo ; enddo ; enddo
      enddo
    endif
  endif ! restart

! the following does not work in layer mode yet
!!  if ( CS%use_sponge ) then
  !   If sponges are used, this example damps tracers in sponges in the
  ! northern half of the domain to 1 and tracers in the southern half
  ! to 0.  For any tracers that are not damped in the sponge, the call
  ! to set_up_sponge_field can simply be omitted.
!    if (.not.associated(ALE_sponge_CSp)) &
!      call MOM_error(FATAL, "ISOMIP_initialize_tracer: "// &
!        "The pointer to ALEsponge_CSp must be associated if SPONGE is defined.")

!    allocate(temp(G%isd:G%ied,G%jsd:G%jed,nz))

!    do j=js,je ; do i=is,ie
!      if (G%geoLonT(i,j) >= 790.0 .AND. G%geoLonT(i,j) <= 800.0) then
!        temp(i,j,:) = 1.0
!      else
!        temp(i,j,:) = 0.0
!      endif
!    enddo ; enddo

      !   do m=1,NTR
!    do m=1,1
      ! This is needed to force the compiler not to do a copy in the sponge
      ! calls.  Curses on the designers and implementers of Fortran90.
!      tr_ptr => CS%tr(:,:,:,m)
!      call set_up_ALE_sponge_field(temp, G, tr_ptr, ALE_sponge_CSp)
!    enddo
!    deallocate(temp)
!  endif

end procedure initialize_ISOMIP_tracer
module procedure ISOMIP_tracer_column_physics
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)) :: h_work ! Used so that h can be modified [H ~> m or kg m-2]
  real :: melt(SZI_(G),SZJ_(G)) ! melt water (positive for melting, negative for freezing) [R Z T-1 ~> kg m-2 s-1]
  real :: mmax                ! The global maximum melting rate [R Z T-1 ~> kg m-2 s-1]
  integer :: i, j, k, is, ie, js, je, nz, m
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke

  if (.not.associated(CS)) return

  melt(:,:) = fluxes%iceshelf_melt(:,:)

  ! max. melt
  mmax = MAXVAL(melt(is:ie,js:je))
  call max_across_PEs(mmax)
  ! write(mesg,*) 'max melt = ', mmax
  ! call MOM_mesg(mesg, 5)
  ! dye melt water (m=1), dye = 1 if melt=max(melt)
  do m=1,NTR
    do j=js,je ; do i=is,ie
      if (melt(i,j) > 0.0) then ! melting
        CS%tr(i,j,1:2,m) = melt(i,j)/mmax ! inject dye in the ML
      else ! freezing
        CS%tr(i,j,1:2,m) = 0.0
      endif
    enddo ; enddo
  enddo

  if (present(evap_CFL_limit) .and. present(minimum_forcing_depth)) then
    do m=1,NTR
      do k=1,nz ;do j=js,je ; do i=is,ie
        h_work(i,j,k) = h_old(i,j,k)
      enddo ; enddo ; enddo
      call applyTracerBoundaryFluxesInOut(G, GV, CS%tr(:,:,:,m), dt, fluxes, h_work, &
                                          evap_CFL_limit, minimum_forcing_depth)
      call tracer_vertdiff(h_work, ea, eb, dt, CS%tr(:,:,:,m), G, GV)
    enddo
  else
    do m=1,NTR
      call tracer_vertdiff(h_old, ea, eb, dt, CS%tr(:,:,:,m), G, GV)
    enddo
  endif

end procedure ISOMIP_tracer_column_physics
module procedure ISOMIP_tracer_surface_state
  integer :: m, is, ie, js, je, isd, ied, jsd, jed
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed

  if (.not.associated(CS)) return

  if (CS%coupled_tracers) then
    do m=1,ntr
      !   This call loads the surface values into the appropriate array in the
      ! coupler-type structure.
      call set_coupler_type_data(CS%tr(:,:,1,m), CS%ind_tr(m), sfc_state%tr_fields, &
                   idim=(/isd, is, ie, ied/), jdim=(/jsd, js, je, jed/), turns=G%HI%turns)
    enddo
  endif

end procedure ISOMIP_tracer_surface_state
module procedure ISOMIP_tracer_end
  if (associated(CS)) then
    if (associated(CS%tr)) deallocate(CS%tr)
    deallocate(CS)
  endif
end procedure ISOMIP_tracer_end
end submodule ISOMIP_tracer_s
