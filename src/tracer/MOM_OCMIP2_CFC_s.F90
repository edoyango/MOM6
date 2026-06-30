submodule (MOM_OCMIP2_CFC) MOM_OCMIP2_CFC_s
#include <MOM_memory.h>
  implicit none
contains
module procedure register_OCMIP2_CFC
  character(len=40)  :: mdl = "MOM_OCMIP2_CFC" ! This module's name.
  character(len=200) :: inputdir ! The directory where NetCDF input files are.
# include "version_variable.h"
  real, dimension(:,:,:), pointer :: tr_ptr => NULL() ! A pointer to a CFC tracer [mol m-3]
  real :: a11_dflt(4), a12_dflt(4) ! Default values of the various coefficients
  real :: d11_dflt(4), d12_dflt(4) ! in the expressions for the solubility and
  real :: e11_dflt(3), e12_dflt(3) ! Schmidt numbers [various units by element].
  character(len=48) :: flux_units ! The units for tracer fluxes.
  integer :: isd, ied, jsd, jed, nz
  isd = HI%isd ; ied = HI%ied ; jsd = HI%jsd ; jed = HI%jed ; nz = GV%ke

  if (associated(CS)) then
    call MOM_error(FATAL, "register_OCMIP2_CFC called with an "// &
                          "associated control structure.")
  endif
  allocate(CS)

  ! This call sets default properties for the air-sea CFC fluxes and obtains the
  ! indicies for the CFC11 and CFC12 flux coupling.
  call flux_init_OCMIP2_CFC(CS, verbosity=3)
  if ((CS%ind_cfc_11_flux < 0) .or. (CS%ind_cfc_12_flux < 0)) then
    ! This is most likely to happen with the dummy version of atmos_ocn_coupler_flux
    ! used in ocean-only runs.
    call MOM_ERROR(WARNING, "CFCs are currently only set up to be run in " // &
                   " coupled model configurations, and will be disabled.")
    deallocate(CS)
    register_OCMIP2_CFC = .false.
    return
  endif

  ! Read all relevant parameters and write them to the model log.
  call log_version(param_file, mdl, version, "")
  call get_param(param_file, mdl, "CFC_IC_FILE", CS%IC_file, &
                 "The file in which the CFC initial values can be "//&
                 "found, or an empty string for internal initialization.", &
                 default=" ")
  if ((len_trim(CS%IC_file) > 0) .and. (scan(CS%IC_file,'/') == 0)) then
    ! Add the directory if CS%IC_file is not already a complete path.
    call get_param(param_file, mdl, "INPUTDIR", inputdir, default=".")
    CS%IC_file = trim(slasher(inputdir))//trim(CS%IC_file)
    call log_param(param_file, mdl, "INPUTDIR/CFC_IC_FILE", CS%IC_file)
  endif
  call get_param(param_file, mdl, "CFC_IC_FILE_IS_Z", CS%Z_IC_file, &
                 "If true, CFC_IC_FILE is in depth space, not layer space", &
                 default=.false.)
  call get_param(param_file, mdl, "TRACERS_MAY_REINIT", CS%tracers_may_reinit, &
                 "If true, tracers may go through the initialization code "//&
                 "if they are not found in the restart files.  Otherwise "//&
                 "it is a fatal error if tracers are not found in the "//&
                 "restart files of a restarted run.", default=.false.)

  !   The following vardesc types contain a package of metadata about each tracer,
  ! including, the name; units; longname; and grid information.
  CS%CFC11_name = "CFC11" ; CS%CFC12_name = "CFC12"
  CS%CFC11_desc = var_desc(CS%CFC11_name,"mol m-3","CFC-11 Concentration", caller=mdl)
  CS%CFC12_desc = var_desc(CS%CFC12_name,"mol m-3","CFC-12 Concentration", caller=mdl)
  if (GV%Boussinesq) then ; flux_units = "mol s-1"
  else ; flux_units = "mol m-3 kg s-1" ; endif

  allocate(CS%CFC11(isd:ied,jsd:jed,nz), source=0.0)
  allocate(CS%CFC12(isd:ied,jsd:jed,nz), source=0.0)

  ! This pointer assignment is needed to force the compiler not to do a copy in
  ! the registration calls.  Curses on the designers and implementers of F90.
  tr_ptr => CS%CFC11
  ! Register CFC11 for horizontal advection, diffusion, and restarts.
  call register_tracer(tr_ptr, tr_Reg, param_file, HI, GV, &
                       tr_desc=CS%CFC11_desc, registry_diags=.true., &
                       flux_units=flux_units, &
                       restart_CS=restart_CS, mandatory=.not.CS%tracers_may_reinit)
  ! Do the same for CFC12
  tr_ptr => CS%CFC12
  call register_tracer(tr_ptr, Tr_Reg, param_file, HI, GV, &
                       tr_desc=CS%CFC12_desc, registry_diags=.true., &
                       flux_units=flux_units, &
                       restart_CS=restart_CS, mandatory=.not.CS%tracers_may_reinit)

  ! Set and read the various empirical coefficients.

!-----------------------------------------------------------------------
! Default Schmidt number coefficients for CFC11 (_11) and CFC12 (_12) are given
! by Zheng et al (1998), JGR vol 103, C1.
!-----------------------------------------------------------------------
  a11_dflt(:) = (/ 3501.8, -210.31,  6.1851, -0.07513 /)
  a12_dflt(:) = (/ 3845.4, -228.95,  6.1908, -0.06743 /)
  call get_param(param_file, mdl, "CFC11_A1", CS%a1_11, &
                 "A coefficient in the Schmidt number of CFC11.", &
                 units="nondim", default=a11_dflt(1))
  call get_param(param_file, mdl, "CFC11_A2", CS%a2_11, &
                 "A coefficient in the Schmidt number of CFC11.", &
                 units="degC-1", default=a11_dflt(2))
  call get_param(param_file, mdl, "CFC11_A3", CS%a3_11, &
                 "A coefficient in the Schmidt number of CFC11.", &
                 units="degC-2", default=a11_dflt(3))
  call get_param(param_file, mdl, "CFC11_A4", CS%a4_11, &
                 "A coefficient in the Schmidt number of CFC11.", &
                 units="degC-3", default=a11_dflt(4))

  call get_param(param_file, mdl, "CFC12_A1", CS%a1_12, &
                 "A coefficient in the Schmidt number of CFC12.", &
                 units="nondim", default=a12_dflt(1))
  call get_param(param_file, mdl, "CFC12_A2", CS%a2_12, &
                 "A coefficient in the Schmidt number of CFC12.", &
                 units="degC-1", default=a12_dflt(2))
  call get_param(param_file, mdl, "CFC12_A3", CS%a3_12, &
                 "A coefficient in the Schmidt number of CFC12.", &
                 units="degC-2", default=a12_dflt(3))
  call get_param(param_file, mdl, "CFC12_A4", CS%a4_12, &
                 "A coefficient in the Schmidt number of CFC12.", &
                 units="degC-3", default=a12_dflt(4))

!-----------------------------------------------------------------------
! Solubility coefficients for alpha in mol/l/atm for CFC11 (_11) and CFC12 (_12)
! after Warner and Weiss (1985) DSR, vol 32.
!-----------------------------------------------------------------------
  d11_dflt(:) = (/ -229.9261, 319.6552, 119.4471, -1.39165 /)
  e11_dflt(:) = (/ -0.142382, 0.091459, -0.0157274 /)
  d12_dflt(:) = (/ -218.0971, 298.9702, 113.8049, -1.39165 /)
  e12_dflt(:) = (/ -0.143566, 0.091015, -0.0153924 /)

  call get_param(param_file, mdl, "CFC11_D1", CS%d1_11, &
                 "A coefficient in the solubility of CFC11.", &
                 units="none", default=d11_dflt(1))
  call get_param(param_file, mdl, "CFC11_D2", CS%d2_11, &
                 "A coefficient in the solubility of CFC11.", &
                 units="hK", default=d11_dflt(2))
  call get_param(param_file, mdl, "CFC11_D3", CS%d3_11, &
                 "A coefficient in the solubility of CFC11.", &
                 units="none", default=d11_dflt(3))
  call get_param(param_file, mdl, "CFC11_D4", CS%d4_11, &
                 "A coefficient in the solubility of CFC11.", &
                 units="hK-2", default=d11_dflt(4))
  call get_param(param_file, mdl, "CFC11_E1", CS%e1_11, &
                 "A coefficient in the solubility of CFC11.", &
                 units="PSU-1", default=e11_dflt(1))
  call get_param(param_file, mdl, "CFC11_E2", CS%e2_11, &
                 "A coefficient in the solubility of CFC11.", &
                 units="PSU-1 hK-1", default=e11_dflt(2))
  call get_param(param_file, mdl, "CFC11_E3", CS%e3_11, &
                 "A coefficient in the solubility of CFC11.", &
                 units="PSU-1 hK-2", default=e11_dflt(3))

  call get_param(param_file, mdl, "CFC12_D1", CS%d1_12, &
                 "A coefficient in the solubility of CFC12.", &
                 units="none", default=d12_dflt(1))
  call get_param(param_file, mdl, "CFC12_D2", CS%d2_12, &
                 "A coefficient in the solubility of CFC12.", &
                 units="hK", default=d12_dflt(2))
  call get_param(param_file, mdl, "CFC12_D3", CS%d3_12, &
                 "A coefficient in the solubility of CFC12.", &
                 units="none", default=d12_dflt(3))
  call get_param(param_file, mdl, "CFC12_D4", CS%d4_12, &
                 "A coefficient in the solubility of CFC12.", &
                 units="hK-2", default=d12_dflt(4))
  call get_param(param_file, mdl, "CFC12_E1", CS%e1_12, &
                 "A coefficient in the solubility of CFC12.", &
                 units="PSU-1", default=e12_dflt(1))
  call get_param(param_file, mdl, "CFC12_E2", CS%e2_12, &
                 "A coefficient in the solubility of CFC12.", &
                 units="PSU-1 hK-1", default=e12_dflt(2))
  call get_param(param_file, mdl, "CFC12_E3", CS%e3_12, &
                 "A coefficient in the solubility of CFC12.", &
                 units="PSU-1 hK-2", default=e12_dflt(3))

  CS%tr_Reg => tr_Reg
  CS%restart_CSp => restart_CS

  register_OCMIP2_CFC = .true.
end procedure register_OCMIP2_CFC
module procedure flux_init_OCMIP2_CFC
  character(len=128) :: default_ice_restart_file = 'ice_ocmip2_cfc.res.nc'
  character(len=128) :: default_ocean_restart_file = 'ocmip2_cfc.res.nc'
  integer :: ind_flux(2) ! Integer indices of the fluxes
  ind_flux(1) = atmos_ocn_coupler_flux('cfc_11_flux', &
       flux_type='air_sea_gas_flux', implementation='ocmip2', &
       param=(/ 9.36e-07, 9.7561e-06 /), &
       ice_restart_file=default_ice_restart_file, &
       ocean_restart_file=default_ocean_restart_file, &
       caller="register_OCMIP2_CFC", verbosity=verbosity)
  ind_flux(2) = atmos_ocn_coupler_flux('cfc_12_flux', &
       flux_type='air_sea_gas_flux', implementation='ocmip2', &
       param=(/ 9.36e-07, 9.7561e-06 /), &
       ice_restart_file=default_ice_restart_file, &
       ocean_restart_file=default_ocean_restart_file, &
       caller="register_OCMIP2_CFC", verbosity=verbosity)

  if (present(CS)) then ; if (associated(CS)) then
    CS%ind_cfc_11_flux = ind_flux(1)
    CS%ind_cfc_12_flux = ind_flux(2)
  endif ; endif

end procedure flux_init_OCMIP2_CFC
module procedure initialize_OCMIP2_CFC
  if (.not.associated(CS)) return

  CS%Time => day
  CS%diag => diag

  if (.not.restart .or. (CS%tracers_may_reinit .and. &
      .not.query_initialized(CS%CFC11, CS%CFC11_name, CS%restart_CSp))) then
    call init_tracer_CFC(h, CS%CFC11, CS%CFC11_name, CS%CFC11_land_val, &
                         CS%CFC11_IC_val, G, GV, US, CS)
    call set_initialized(CS%CFC11, CS%CFC11_name, CS%restart_CSp)
  endif

  if (.not.restart .or. (CS%tracers_may_reinit .and. &
      .not.query_initialized(CS%CFC12, CS%CFC12_name, CS%restart_CSp))) then
    call init_tracer_CFC(h, CS%CFC12, CS%CFC12_name, CS%CFC12_land_val, &
                         CS%CFC12_IC_val, G, GV, US, CS)
    call set_initialized(CS%CFC12, CS%CFC12_name, CS%restart_CSp)
  endif

  if (associated(OBC)) then
  ! Steal from updated DOME in the fullness of time.
  endif

end procedure initialize_OCMIP2_CFC
module procedure init_tracer_CFC
  logical :: OK
  integer :: i, j, k, is, ie, js, je, nz
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke

  if (len_trim(CS%IC_file) > 0) then
    !  Read the tracer concentrations from a netcdf file.
    if (.not.file_exists(CS%IC_file, G%Domain)) &
      call MOM_error(FATAL, "initialize_OCMIP2_CFC: Unable to open "//CS%IC_file)
    if (CS%Z_IC_file) then
      OK = tracer_Z_init(tr, h, CS%IC_file, name, G, GV, US)
      if (.not.OK) then
        OK = tracer_Z_init(tr, h, CS%IC_file, trim(name), G, GV, US)
        if (.not.OK) call MOM_error(FATAL,"initialize_OCMIP2_CFC: "//&
                "Unable to read "//trim(name)//" from "//&
                trim(CS%IC_file)//".")
      endif
    else
      call MOM_read_data(CS%IC_file, trim(name), tr, G%Domain)
    endif
  else
    do k=1,nz ; do j=js,je ; do i=is,ie
      if (G%mask2dT(i,j) < 0.5) then
        tr(i,j,k) = land_val
      else
        tr(i,j,k) = IC_val
      endif
    enddo ; enddo ; enddo
  endif

end procedure init_tracer_CFC
module procedure OCMIP2_CFC_column_physics
  real, dimension(SZI_(G),SZJ_(G)) :: &
    CFC11_flux, &    ! The fluxes of CFC11 and CFC12 into the ocean, in unscaled units of
    CFC12_flux       ! CFC concentrations times a vertical mass flux [mol R Z m-3 T-1 ~> mol kg m-3 s-1]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)) :: h_work ! Used so that h can be modified [H ~> m or kg m-2]
  integer :: i, j, k, is, ie, js, je, nz, idim(4), jdim(4)
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke
  idim(:) = (/G%isd, is, ie, G%ied/) ; jdim(:) = (/G%jsd, js, je, G%jed/)

  if (.not.associated(CS)) return

  ! These two calls unpack the fluxes from the input arrays.
  !   The -GV%Rho0 changes the sign convention of the flux and with the scaling factors changes
  ! the units of the flux from [conc m s-1] to [conc R Z T-1 ~> conc kg m-2 s-1].
  call extract_coupler_type_data(fluxes%tr_fluxes, CS%ind_cfc_11_flux, CFC11_flux, &
                                 scale_factor=-GV%Rho0*US%m_to_Z*US%T_to_s, idim=idim, jdim=jdim, turns=G%HI%turns)
  call extract_coupler_type_data(fluxes%tr_fluxes, CS%ind_cfc_12_flux, CFC12_flux, &
                                 scale_factor=-GV%Rho0*US%m_to_Z*US%T_to_s, idim=idim, jdim=jdim, turns=G%HI%turns)

  ! Use a tridiagonal solver to determine the concentrations after the
  ! surface source is applied and diapycnal advection and diffusion occurs.
  if (present(evap_CFL_limit) .and. present(minimum_forcing_depth)) then
    do k=1,nz ;do j=js,je ; do i=is,ie
      h_work(i,j,k) = h_old(i,j,k)
    enddo ; enddo ; enddo
    call applyTracerBoundaryFluxesInOut(G, GV, CS%CFC11, dt, fluxes, h_work, &
                                        evap_CFL_limit, minimum_forcing_depth)
    call tracer_vertdiff(h_work, ea, eb, dt, CS%CFC11, G, GV, sfc_flux=CFC11_flux)

    do k=1,nz ;do j=js,je ; do i=is,ie
      h_work(i,j,k) = h_old(i,j,k)
    enddo ; enddo ; enddo
    call applyTracerBoundaryFluxesInOut(G, GV, CS%CFC12, dt, fluxes, h_work, &
                                        evap_CFL_limit, minimum_forcing_depth)
    call tracer_vertdiff(h_work, ea, eb, dt, CS%CFC12, G, GV, sfc_flux=CFC12_flux)
  else
    call tracer_vertdiff(h_old, ea, eb, dt, CS%CFC11, G, GV, sfc_flux=CFC11_flux)
    call tracer_vertdiff(h_old, ea, eb, dt, CS%CFC12, G, GV, sfc_flux=CFC12_flux)
  endif

  ! Write out any desired diagnostics from tracer sources & sinks here.

end procedure OCMIP2_CFC_column_physics
module procedure OCMIP2_CFC_stock
  OCMIP2_CFC_stock = 0
  if (.not.associated(CS)) return

  if (present(stock_index)) then ; if (stock_index > 0) then
    ! Check whether this stock is available from this routine.

    ! No stocks from this routine are being checked yet.  Return 0.
    return
  endif ; endif

  call query_vardesc(CS%CFC11_desc, name=names(1), units=units(1), caller="OCMIP2_CFC_stock")
  call query_vardesc(CS%CFC12_desc, name=names(2), units=units(2), caller="OCMIP2_CFC_stock")
  units(1) = trim(units(1))//" kg" ; units(2) = trim(units(2))//" kg"

  stocks(1) = global_mass_int_EFP(h, G, GV, CS%CFC11, on_PE_only=.true.)
  stocks(2) = global_mass_int_EFP(h, G, GV, CS%CFC12, on_PE_only=.true.)

  OCMIP2_CFC_stock = 2

end procedure OCMIP2_CFC_stock
module procedure OCMIP2_CFC_surface_state
  real, dimension(SZI_(G),SZJ_(G)) :: &
    CFC11_Csurf, &  ! The CFC-11 surface concentrations times the Schmidt number term [mol m-3].
    CFC12_Csurf, &  ! The CFC-12 surface concentrations times the Schmidt number term [mol m-3].
    CFC11_alpha, &  ! The CFC-11 solubility [mol m-3 pptv-1].
    CFC12_alpha     ! The CFC-12 solubility [mol m-3 pptv-1].
  real :: ta        ! Absolute sea surface temperature [hectoKelvin] (Why use such bizzare units?)
  real :: sal       ! Surface salinity [PSU].
  real :: SST       ! Sea surface temperature [degC].
  real :: alpha_11  ! The solubility of CFC 11 [mol m-3 pptv-1].
  real :: alpha_12  ! The solubility of CFC 12 [mol m-3 pptv-1].
  real :: sc_11, sc_12 ! The Schmidt numbers of CFC 11 and CFC 12 [nondim].
  real :: sc_no_term   ! A term related to the Schmidt number [nondim].
  integer :: i, j, is, ie, js, je, idim(4), jdim(4)
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec
  idim(:) = (/G%isd, is, ie, G%ied/) ; jdim(:) = (/G%jsd, js, je, G%jed/)

  if (.not.associated(CS)) return

  do j=js,je ; do i=is,ie
    ta = max(0.01, (US%C_to_degC*sfc_state%SST(i,j) + 273.15) * 0.01) ! Why is this in hectoKelvin?
    sal = US%S_to_ppt*sfc_state%SSS(i,j) ; SST = US%C_to_degC*sfc_state%SST(i,j)
    !    Calculate solubilities using Warner and Weiss (1985) DSR, vol 32.
    ! The final result is in mol/cm3/pptv (1 part per trillion 1e-12)
    ! Use Bullister and Wisegavger for CCl4.
    ! The factor 1.e-09 converts from mol/(l * atm) to mol/(m3 * pptv).
    alpha_11 = exp(CS%d1_11 + CS%d2_11/ta + CS%d3_11*log(ta) + CS%d4_11*ta**2 +&
                   sal * ((CS%e3_11 * ta + CS%e2_11) * ta + CS%e1_11)) * &
               1.0e-09 * G%mask2dT(i,j)
    alpha_12 = exp(CS%d1_12 + CS%d2_12/ta + CS%d3_12*log(ta) + CS%d4_12*ta**2 +&
                   sal * ((CS%e3_12 * ta + CS%e2_12) * ta + CS%e1_12)) * &
               1.0e-09 * G%mask2dT(i,j)
    !   Calculate Schmidt numbers using coefficients given by
    ! Zheng et al (1998), JGR vol 103, C1.
    sc_11 = CS%a1_11 + SST * (CS%a2_11 + SST * (CS%a3_11 + SST * CS%a4_11)) * &
            G%mask2dT(i,j)
    sc_12 = CS%a1_12 + SST * (CS%a2_12 + SST * (CS%a3_12 + SST * CS%a4_12)) * &
            G%mask2dT(i,j)
    ! The abs here is to avoid NaNs. The model should be failing at this point.
    sc_no_term = sqrt(660.0 / (abs(sc_11) + 1.0e-30))
    CFC11_alpha(i,j) = alpha_11 * sc_no_term
    CFC11_Csurf(i,j) = CS%CFC11(i,j,1) * sc_no_term

    sc_no_term = sqrt(660.0 / (abs(sc_12) + 1.0e-30))
    CFC12_alpha(i,j) = alpha_12 * sc_no_term
    CFC12_Csurf(i,j) = CS%CFC12(i,j,1) * sc_no_term
  enddo ; enddo

  !   These calls load these values into the appropriate arrays in the
  ! coupler-type structure.
  call set_coupler_type_data(CFC11_alpha, CS%ind_cfc_11_flux, sfc_state%tr_fields, &
                             solubility=.true., idim=idim, jdim=jdim, turns=G%HI%turns)
  call set_coupler_type_data(CFC11_Csurf, CS%ind_cfc_11_flux, sfc_state%tr_fields, &
                             idim=idim, jdim=jdim, turns=G%HI%turns)
  call set_coupler_type_data(CFC12_alpha, CS%ind_cfc_12_flux, sfc_state%tr_fields, &
                             solubility=.true., idim=idim, jdim=jdim, turns=G%HI%turns)
  call set_coupler_type_data(CFC12_Csurf, CS%ind_cfc_12_flux, sfc_state%tr_fields, &
                             idim=idim, jdim=jdim, turns=G%HI%turns)

end procedure OCMIP2_CFC_surface_state
module procedure OCMIP2_CFC_end
  if (associated(CS)) then
    if (associated(CS%CFC11)) deallocate(CS%CFC11)
    if (associated(CS%CFC12)) deallocate(CS%CFC12)

    deallocate(CS)
  endif
end procedure OCMIP2_CFC_end
end submodule MOM_OCMIP2_CFC_s
