submodule (DOME_tracer) DOME_tracer_s
#include <MOM_memory.h>
  implicit none
contains
module procedure register_DOME_tracer
  character(len=80)  :: name, longname
# include "version_variable.h"
  character(len=40)  :: mdl = "DOME_tracer" ! This module's name.
  character(len=48) :: flux_units ! The units for tracer fluxes, usually
  character(len=200) :: inputdir
  real, pointer :: tr_ptr(:,:,:) => NULL() ! A pointer to one of the tracers, perhaps in [g kg-1]
  integer :: isd, ied, jsd, jed, nz, m
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed ; nz = GV%ke

  if (associated(CS)) then
    call MOM_error(FATAL, "DOME_register_tracer called with an "// &
                          "associated control structure.")
  endif
  allocate(CS)

  ! Read all relevant parameters and write them to the model log.
  call log_version(param_file, mdl, version, "")
  call get_param(param_file, mdl, "DOME_TRACER_IC_FILE", CS%tracer_IC_file, &
                 "The name of a file from which to read the initial "//&
                 "conditions for the DOME tracers, or blank to initialize "//&
                 "them internally.", default=" ")
  if (len_trim(CS%tracer_IC_file) >= 1) then
    call get_param(param_file, mdl, "INPUTDIR", inputdir, default=".")
    inputdir = slasher(inputdir)
    CS%tracer_IC_file = trim(inputdir)//trim(CS%tracer_IC_file)
    call log_param(param_file, mdl, "INPUTDIR/DOME_TRACER_IC_FILE", &
                   CS%tracer_IC_file)
  endif
  call get_param(param_file, mdl, "DOME_TRACER_STRIPE_WIDTH", CS%stripe_width, &
                 "The meridional width of the vertical stripes in the initial condition "//&
                 "for the DOME tracers.", units=G%y_ax_unit_short, default=50.0)
  call get_param(param_file, mdl, "DOME_TRACER_STRIPE_LAT", CS%stripe_s_lat, &
                 "The southern latitude of the first vertical stripe in the initial condition "//&
                 "for the DOME tracers.", units=G%y_ax_unit_short, default=350.0)
  call get_param(param_file, mdl, "DOME_TRACER_SHEET_SPACING", CS%sheet_spacing, &
                 "The vertical spacing between successive horizontal sheets of tracer in the initial "//&
                 "conditions for the DOME tracers, and twice the thickness of these tracer sheets.", &
                 units="m", default=600.0, scale=US%m_to_Z)
  call get_param(param_file, mdl, "SPONGE", CS%use_sponge, &
                 "If true, sponges may be applied anywhere in the domain. "//&
                 "The exact location and properties of those sponges are "//&
                 "specified from MOM_initialization.F90.", default=.false.)

  allocate(CS%tr(isd:ied,jsd:jed,nz,NTR), source=0.0)

  do m=1,NTR
    write(name,'("tr_D",I0)') m
    write(longname,'("Concentration of DOME Tracer ",I2.2)') m
    CS%tr_desc(m) = var_desc(name, units="kg kg-1", longname=longname, caller=mdl)
    if (GV%Boussinesq) then ; flux_units = "kg kg-1 m3 s-1"
    else ; flux_units = "kg s-1" ; endif

    ! This is needed to force the compiler not to do a copy in the registration
    ! calls.  Curses on the designers and implementers of Fortran90.
    tr_ptr => CS%tr(:,:,:,m)
    ! Register the tracer for horizontal advection, diffusion, and restarts.
    call register_tracer(tr_ptr, tr_Reg, param_file, G%HI, GV, &
                         name=name, longname=longname, units="kg kg-1", &
                         registry_diags=.true., restart_CS=restart_CS, &
                         flux_units=trim(flux_units), flux_scale=GV%H_to_MKS)

    !   Set coupled_tracers to be true (hard-coded above) to provide the surface
    ! values to the coupler (if any).  This is meta-code and its arguments will
    ! currently (deliberately) give fatal errors if it is used.
    if (CS%coupled_tracers) &
      CS%ind_tr(m) = atmos_ocn_coupler_flux(trim(name)//'_flux', &
          flux_type=' ', implementation=' ', caller="register_DOME_tracer")
  enddo

  CS%tr_Reg => tr_Reg
  register_DOME_tracer = .true.
end procedure register_DOME_tracer
module procedure initialize_DOME_tracer
  real, allocatable :: temp(:,:,:) ! Target values for the tracers in the sponges, perhaps in [g kg-1]
  character(len=16) :: name     ! A variable's name in a NetCDF file.
  real, pointer :: tr_ptr(:,:,:) => NULL() ! A pointer to one of the tracers, perhaps in [g kg-1]
  real :: dz(SZI_(G),SZK_(GV)) ! Height change across layers [Z ~> m]
  real :: tr_y   ! Initial zonally uniform tracer concentrations, perhaps in [g kg-1]
  real :: dz_neglect        ! A thickness that is so small it is usually lost
  real :: e(SZK_(GV)+1)     ! Interface heights relative to the sea surface (negative down) [Z ~> m]
  real :: e_top  ! Height of the top of the tracer band relative to the sea surface [Z ~> m]
  real :: e_bot  ! Height of the bottom of the tracer band relative to the sea surface [Z ~> m]
  real :: d_tr   ! A change in tracer concentrations, in tracer units, perhaps [g kg-1]
  integer :: i, j, k, is, ie, js, je, isd, ied, jsd, jed, nz, m
  if (.not.associated(CS)) return
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed

  dz_neglect = GV%dz_subroundoff

  CS%Time => day
  CS%diag => diag

  if (.not.restart) then
    if (len_trim(CS%tracer_IC_file) >= 1) then
      !  Read the tracer concentrations from a netcdf file.
      if (.not.file_exists(CS%tracer_IC_file, G%Domain)) &
        call MOM_error(FATAL, "DOME_initialize_tracer: Unable to open "// &
                        CS%tracer_IC_file)
      do m=1,NTR
        call query_vardesc(CS%tr_desc(m), name, caller="initialize_DOME_tracer")
        call MOM_read_data(CS%tracer_IC_file, trim(name), CS%tr(:,:,:,m), G%Domain)
      enddo
    else
      do m=1,NTR
        do k=1,nz ; do j=js,je ; do i=is,ie
          CS%tr(i,j,k,m) = 1.0e-20 ! This could just as well be 0.
        enddo ; enddo ; enddo
      enddo

!    This sets a stripe of tracer across the basin.
      do m=2,min(6,NTR) ; do j=js,je ; do i=is,ie
        tr_y = 0.0
        if ((G%geoLatT(i,j) > (CS%stripe_s_lat + CS%stripe_width*real(m-2))) .and. &
            (G%geoLatT(i,j) < (CS%stripe_s_lat + CS%stripe_width*real(m-1)))) &
          tr_y = 1.0
        do k=1,nz
!      This adds the stripes of tracer to every layer.
            CS%tr(i,j,k,m) = CS%tr(i,j,k,m) + tr_y
        enddo
      enddo ; enddo ; enddo

      if (NTR >= 7) then
        do j=js,je
          call thickness_to_dz(h, tv, dz, j, G, GV)
          do i=is,ie
            e(1) = 0.0
            do k=1,nz
              e(K+1) = e(K) - dz(i,k)
              do m=7,NTR
                e_top = -CS%sheet_spacing * (real(m-6))
                e_bot = -CS%sheet_spacing * (real(m-6) + 0.5)
                if (e_top < e(K)) then
                  if (e_top < e(K+1)) then ; d_tr = 0.0
                  elseif (e_bot < e(K+1)) then
                    d_tr = 1.0 * (e_top-e(K+1)) / (dz(i,k)+dz_neglect)
                  else ; d_tr = 1.0 * (e_top-e_bot) / (dz(i,k)+dz_neglect)
                  endif
                elseif (e_bot < e(K)) then
                  if (e_bot < e(K+1)) then ; d_tr = 1.0
                  else ; d_tr = 1.0 * (e(K)-e_bot) / (dz(i,k)+dz_neglect)
                  endif
                else
                  d_tr = 0.0
                endif
                if (dz(i,k) < 2.0*GV%Angstrom_Z) d_tr=0.0
                CS%tr(i,j,k,m) = CS%tr(i,j,k,m) + d_tr
              enddo
            enddo
          enddo
        enddo
      endif

    endif
  endif ! restart

  if ( CS%use_sponge ) then
!   If sponges are used, this example damps tracers in sponges in the
! northern half of the domain to 1 and tracers in the southern half
! to 0.  For any tracers that are not damped in the sponge, the call
! to set_up_sponge_field can simply be omitted.
    if (.not.associated(sponge_CSp)) &
      call MOM_error(FATAL, "DOME_initialize_tracer: "// &
        "The pointer to sponge_CSp must be associated if SPONGE is defined.")

    allocate(temp(G%isd:G%ied,G%jsd:G%jed,nz))
    do k=1,nz ; do j=js,je ; do i=is,ie
      if (G%geoLatT(i,j) > 700.0 .and. (k > nz/2)) then
        temp(i,j,k) = 1.0
      else
        temp(i,j,k) = 0.0
      endif
    enddo ; enddo ; enddo

!   do m=1,NTR
    do m=1,1
      ! This pointer is needed to force the compiler not to do a copy in the sponge calls.
      tr_ptr => CS%tr(:,:,:,m)
      call set_up_sponge_field(temp, tr_ptr, G, GV, nz, sponge_CSp)
    enddo
    deallocate(temp)
  endif

end procedure initialize_DOME_tracer
module procedure DOME_tracer_column_physics
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)) :: h_work ! Used so that h can be modified [H ~> m or kg m-2]
  integer :: i, j, k, is, ie, js, je, nz, m
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke

  if (.not.associated(CS)) return

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

end procedure DOME_tracer_column_physics
module procedure DOME_tracer_surface_state
  integer :: m, is, ie, js, je, isd, ied, jsd, jed
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed

  if (.not.associated(CS)) return

  if (CS%coupled_tracers) then
    do m=1,NTR
      !   This call loads the surface values into the appropriate array in the
      ! coupler-type structure.
      call set_coupler_type_data(CS%tr(:,:,1,m), CS%ind_tr(m), sfc_state%tr_fields, &
                   idim=(/isd, is, ie, ied/), jdim=(/jsd, js, je, jed/), turns=G%HI%turns)
    enddo
  endif

end procedure DOME_tracer_surface_state
module procedure DOME_tracer_end
  if (associated(CS)) then
    if (associated(CS%tr)) deallocate(CS%tr)
    deallocate(CS)
  endif
end procedure DOME_tracer_end
end submodule DOME_tracer_s
