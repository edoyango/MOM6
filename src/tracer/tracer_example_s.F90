submodule (USER_tracer_example) USER_tracer_example_s
#include <MOM_memory.h>
  implicit none
contains
module procedure USER_register_tracer_example
  character(len=80)  :: name, longname
# include "version_variable.h"
  character(len=40)  :: mdl = "tracer_example" ! This module's name.
  character(len=200) :: inputdir
  character(len=48) :: flux_units ! The units for tracer fluxes, usually
  real, pointer :: tr_ptr(:,:,:) => NULL() ! A pointer to one of the tracers, perhaps in [g kg-1]
  integer :: isd, ied, jsd, jed, nz, m
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed ; nz = GV%ke

  if (associated(CS)) then
    call MOM_error(FATAL, "USER_register_tracer_example called with an "// &
                          "associated control structure.")
  endif
  allocate(CS)

  ! Read all relevant parameters and write them to the model log.
  call log_version(param_file, mdl, version, "")
  call get_param(param_file, mdl, "TRACER_EXAMPLE_IC_FILE", CS%tracer_IC_file, &
                 "The name of a file from which to read the initial conditions for "//&
                 "the tracer_example tracers, or blank to initialize them internally.", &
                 default=" ")
  if (len_trim(CS%tracer_IC_file) >= 1) then
    call get_param(param_file, mdl, "INPUTDIR", inputdir, default=".")
    CS%tracer_IC_file = trim(slasher(inputdir))//trim(CS%tracer_IC_file)
    call log_param(param_file, mdl, "INPUTDIR/TRACER_EXAMPLE_IC_FILE", &
                   CS%tracer_IC_file)
  endif
  call get_param(param_file, mdl, "SPONGE", CS%use_sponge, &
                 "If true, sponges may be applied anywhere in the domain. "//&
                 "The exact location and properties of those sponges are "//&
                 "specified from MOM_initialization.F90.", default=.false.)
  call get_param(param_file, mdl, "TRACER_EXAMPLE_STRIPE_WIDTH", CS%stripe_width, &
                 "The Gaussian width of the stripe in the initial condition for the "//&
                 "tracer_example tracers.", units="m", default=1.0e5, scale=US%m_to_L)
  call get_param(param_file, mdl, "TRACER_EXAMPLE_STRIPE_LAT", CS%stripe_lat, &
                 "The central latitude of the stripe in the initial condition for the "//&
                 "tracer_example tracers.", units=G%y_ax_unit_short, default=40.0)

  allocate(CS%tr(isd:ied,jsd:jed,nz,NTR), source=0.0)

  do m=1,NTR
    write(name,'("tr",I0)') m
    write(longname,'("Concentration of Tracer ",I2.2)') m
    CS%tr_desc(m) = var_desc(name, units="kg kg-1", longname=longname, caller=mdl)

    ! This needs to be changed if the units of tracer are changed above.
    if (GV%Boussinesq) then ; flux_units = "kg kg-1 m3 s-1"
    else ; flux_units = "kg s-1" ; endif

    ! This pointer is needed to force the compiler not to do a copy in the registration calls.
    tr_ptr => CS%tr(:,:,:,m)
    ! Register the tracer for horizontal advection, diffusion, and restarts.
    call register_tracer(tr_ptr, tr_Reg, param_file, G%HI, GV, &
                         name=name, longname=longname, units="kg kg-1", &
                         registry_diags=.true., flux_units=flux_units, &
                         restart_CS=restart_CS)

    !   Set coupled_tracers to be true (hard-coded above) to provide the surface
    ! values to the coupler (if any).  This is meta-code and its arguments will
    ! currently (deliberately) give fatal errors if it is used.
    if (CS%coupled_tracers) &
      CS%ind_tr(m) = atmos_ocn_coupler_flux(trim(name)//'_flux', &
          flux_type=' ', implementation=' ', caller="USER_register_tracer_example")
  enddo

  CS%tr_Reg => tr_Reg
  USER_register_tracer_example = .true.
end procedure USER_register_tracer_example
module procedure USER_initialize_tracer
  real, allocatable :: temp(:,:,:) ! Target values for the tracers in the sponges, perhaps in [g kg-1]
  character(len=32) :: name     ! A variable's name in a NetCDF file.
  real, pointer :: tr_ptr(:,:,:) => NULL() ! A pointer to one of the tracers, perhaps in [g kg-1]
  real :: PI     ! 3.1415926... calculated as 4*atan(1) [nondim]
  real :: tr_y   ! Initial zonally uniform tracer concentrations, perhaps in [g kg-1]
  real :: dist2  ! The distance squared from a line [L2 ~> m2].
  integer :: i, j, k, is, ie, js, je, isd, ied, jsd, jed, nz, m
  integer :: IsdB, IedB, JsdB, JedB, lntr
  if (.not.associated(CS)) return
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  IsdB = G%IsdB ; IedB = G%IedB ; JsdB = G%JsdB ; JedB = G%JedB

  lntr = NTR ! Avoids compile-time warning when NTR<2
  CS%Time => day
  CS%diag => diag

  if (.not.restart) then
    if (len_trim(CS%tracer_IC_file) >= 1) then
!  Read the tracer concentrations from a netcdf file.
      if (.not.file_exists(CS%tracer_IC_file, G%Domain)) &
        call MOM_error(FATAL, "USER_initialize_tracer: Unable to open "// &
                        CS%tracer_IC_file)
      do m=1,NTR
        call query_vardesc(CS%tr_desc(m), name, caller="USER_initialize_tracer")
        call MOM_read_data(CS%tracer_IC_file, trim(name), CS%tr(:,:,:,m), G%Domain)
      enddo
    else
      do m=1,NTR
        do k=1,nz ; do j=js,je ; do i=is,ie
          CS%tr(i,j,k,m) = 1.0e-20 ! This could just as well be 0.
        enddo ; enddo ; enddo
      enddo

!    This sets a stripe of tracer across the basin.
      PI = 4.0*atan(1.0)
      do j=js,je
        dist2 = (G%Rad_Earth_L * PI / 180.0)**2 * (G%geoLatT(i,j) - CS%stripe_lat)**2
        tr_y = 0.5 * exp( -dist2 / CS%stripe_width**2 )

        do k=1,nz ; do i=is,ie
!      This adds the stripes of tracer to every layer.
          CS%tr(i,j,k,1) = CS%tr(i,j,k,1) + tr_y
        enddo ; enddo
      enddo
    endif
  endif ! restart

  if ( CS%use_sponge ) then
!   If sponges are used, this example damps tracers in sponges in the
! northern half of the domain to 1 and tracers in the southern half
! to 0.  For any tracers that are not damped in the sponge, the call
! to set_up_sponge_field can simply be omitted.
    if (.not.associated(sponge_CSp)) &
      call MOM_error(FATAL, "USER_initialize_tracer: "// &
        "The pointer to sponge_CSp must be associated if SPONGE is defined.")

    allocate(temp(G%isd:G%ied,G%jsd:G%jed,nz))
    do k=1,nz ; do j=js,je ; do i=is,ie
      if ((G%geoLatT(i,j) > 0.5*G%len_lat + G%south_lat) .and. (k > nz/2)) then
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

  if (associated(OBC)) then
    call query_vardesc(CS%tr_desc(1), name, caller="USER_initialize_tracer")
    if (OBC%specified_v_BCs_exist_globally) then
      ! Steal from updated DOME in the fullness of time.
    else
      ! Steal from updated DOME in the fullness of time.
    endif
    ! All tracers but the first have 0 concentration in their inflows. As this
    ! is the default value, the following calls are unnecessary.
    !do m=2,lntr
    do m=2,ntr
      call query_vardesc(CS%tr_desc(m), name, caller="USER_initialize_tracer")
      ! Steal from updated DOME in the fullness of time.
    enddo
  endif

end procedure USER_initialize_tracer
module procedure tracer_column_physics
  real :: hold0(SZI_(G))       ! The original topmost layer thickness,
  real :: b1(SZI_(G))          ! b1 is a variable used by the tridiagonal solver [H ~> m or kg m-2].
  real :: c1(SZI_(G),SZK_(GV)) ! c1 is a variable used by the tridiagonal solver [nondim].
  real :: d1(SZI_(G))          ! d1=1-c1 is used by the tridiagonal solver [nondim].
  real :: h_neglect            ! A thickness that is so small it is usually lost
  real :: b_denom_1            ! The first term in the denominator of b1 [H ~> m or kg m-2].
  real :: diapyc_filt          ! A multiplicative filter that can be set to 0 to disable diapycnal
  real :: dye_up               ! The tracer concentration of upwelled water, perhaps in [g kg-1]?
  real :: dye_down             ! The tracer concentration of downwelled water, perhaps in [g kg-1]?
  integer :: i, j, k, is, ie, js, je, nz, m
  diapyc_filt = 1.0 ; dye_down = 0.0 ; dye_down = 0.0

  ! Uncomment the following line to dye downwelling.
!  diapyc_filt = 0.0 ; dye_down = 1.0
  ! Uncomment the following line to dye upwelling.
!  diapyc_filt = 0.0 ; dye_up = 1.0
  ! Uncomment the following line for tracer concentrations to be set
  ! to zero in any diapycnal motions.
!  diapyc_filt = 0.0 ; dye_down = 0.0 ; dye_down = 0.0

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke

  if (.not.associated(CS)) return
  h_neglect = GV%H_subroundoff

  do j=js,je
    do i=is,ie
!   The following line is appropriate for quantities like salinity
! that are left behind by evaporation, and any surface fluxes would
! be explicitly included in the flux structure.
      hold0(i) = h_old(i,j,1)
!   The following line is appropriate for quantities like temperature
! that can be assumed to have the same concentration in evaporation
! as they had in the water.  The explicit surface fluxes here would
! reflect differences in concentration from the ambient water, not
! the absolute fluxes.
  !   hold0(i) = h_old(i,j,1) + ea(i,j,1)
      b_denom_1 = h_old(i,j,1) + ea(i,j,1) + h_neglect
      b1(i) = 1.0 / (b_denom_1 + eb(i,j,1))
!       d1(i) = b_denom_1 * b1(i)
      d1(i) = diapyc_filt * (b_denom_1 * b1(i)) + (1.0 - diapyc_filt)
      do m=1,NTR
        CS%tr(i,j,1,m) = b1(i)*(hold0(i)*CS%tr(i,j,1,m) + dye_up*eb(i,j,1))
 !      Add any surface tracer fluxes to the preceding line.
      enddo
    enddo
    do k=2,nz ; do i=is,ie
      c1(i,k) = diapyc_filt * eb(i,j,k-1) * b1(i)
      b_denom_1 = h_old(i,j,k) + d1(i)*ea(i,j,k) + h_neglect
      b1(i) = 1.0 / (b_denom_1 + eb(i,j,k))
      d1(i) = diapyc_filt * (b_denom_1 * b1(i)) + (1.0 - diapyc_filt)
      do m=1,NTR
        CS%tr(i,j,k,m) = b1(i) * (h_old(i,j,k)*CS%tr(i,j,k,m) + &
                 ea(i,j,k)*(diapyc_filt*CS%tr(i,j,k-1,m) + dye_down) + &
                 eb(i,j,k)*dye_up)
      enddo
    enddo ; enddo
    do m=1,NTR ; do k=nz-1,1,-1 ; do i=is,ie
      CS%tr(i,j,k,m) = CS%tr(i,j,k,m) + c1(i,k+1)*CS%tr(i,j,k+1,m)
    enddo ; enddo ; enddo
  enddo

end procedure tracer_column_physics
module procedure USER_tracer_stock
  integer :: m
  USER_tracer_stock = 0
  if (.not.associated(CS)) return

  if (present(stock_index)) then ; if (stock_index > 0) then
    ! Check whether this stock is available from this routine.

    ! No stocks from this routine are being checked yet.  Return 0.
    return
  endif ; endif

  do m=1,NTR
    call query_vardesc(CS%tr_desc(m), name=names(m), units=units(m), caller="USER_tracer_stock")
    units(m) = trim(units(m))//" kg"
    stocks(m) = global_mass_int_EFP(h, G, GV, CS%tr(:,:,:,m), on_PE_only=.true.)
  enddo
  USER_tracer_stock = NTR

end procedure USER_tracer_stock
module procedure USER_tracer_surface_state
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

end procedure USER_tracer_surface_state
module procedure USER_tracer_example_end
  if (associated(CS)) then
    if (associated(CS%tr)) deallocate(CS%tr)
    deallocate(CS)
  endif
end procedure USER_tracer_example_end
end submodule USER_tracer_example_s
