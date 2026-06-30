submodule (MOM_CFC_cap) MOM_CFC_cap_s
#include <MOM_memory.h>
  implicit none
contains
module procedure register_CFC_cap
  character(len=40)  :: mdl = "MOM_CFC_cap" ! This module's name.
# include "version_variable.h"
  character(len=200) :: inputdir ! The directory where NetCDF input files are.
  real, dimension(:,:,:), pointer :: tr_ptr => NULL() ! A pointer to a CFC tracer [mol kg-1]
  character(len=200) :: CFC_BC_file           ! filename with cfc11 and cfc12 data
  character(len=30)  :: CFC_BC_var_name       ! varname of field in CFC_BC_file
  character :: m2char
  integer :: isd, ied, jsd, jed, nz, m
  integer :: CFC_BC_data_year   ! specific year in CFC BC data calendar
  integer :: CFC_BC_model_year  ! model year corresponding to CFC_BC_data_year
  isd = HI%isd ; ied = HI%ied ; jsd = HI%jsd ; jed = HI%jed ; nz = GV%ke

  if (associated(CS)) then
    call MOM_error(FATAL, "register_CFC_cap called with an "// &
                          "associated control structure.")
  endif
  allocate(CS)

  ! Read all relevant parameters and write them to the model log.
  call log_version(param_file, mdl, version, "")
  call get_param(param_file, mdl, "DEBUG", CS%debug, &
                 "If true, write out verbose debugging data.", &
                 default=.false., debuggingParam=.true.)
  call get_param(param_file, mdl, "CFC_IC_FILE", CS%IC_file, &
                 "The file in which the CFC initial values can be "//&
                 "found, or an empty string for internal initialization.", &
                 default=" ")
  if ((len_trim(CS%IC_file) > 0) .and. (scan(CS%IC_file, '/') == 0)) then
    ! Add the directory if CS%IC_file is not already a complete path.
    call get_param(param_file, mdl, "INPUTDIR", inputdir, default=".")
    CS%IC_file = trim(slasher(inputdir))//trim(CS%IC_file)
    call log_param(param_file, mdl, "INPUTDIR/CFC_IC_FILE", CS%IC_file, &
                   "full path of CFC_IC_FILE")
  endif
  call get_param(param_file, mdl, "CFC_IC_FILE_IS_Z", CS%Z_IC_file, &
                 "If true, CFC_IC_FILE is in depth space, not layer space", &
                 default=.false.)
  call get_param(param_file, mdl, "TRACERS_MAY_REINIT", CS%tracers_may_reinit, &
                 "If true, tracers may go through the initialization code "//&
                 "if they are not found in the restart files.  Otherwise "//&
                 "it is a fatal error if tracers are not found in the "//&
                 "restart files of a restarted run.", default=.false.)
  do m=1,NTR
    write(m2char, "(I1)") m
    call get_param(param_file, mdl, "CFC1"//m2char//"_IC_VAL", CS%CFC_data(m)%IC_val, &
                   "Value that CFC_1"//m2char//" is set to when it is not read from a file.", &
                   units="mol kg-1", default=0.0)
  enddo

  ! the following params are not used in this module. Instead, they are used in
  ! the cap but are logged here to keep all the CFC cap params together.
  call get_param(param_file, mdl, "CFC_BC_FILE", CFC_BC_file, &
                 "The file in which the CFC-11 and CFC-12 atm concentrations can be "//&
                 "found (units must be parts per trillion).", default=" ")
  if (len_trim(CFC_BC_file) == 0) then
    call MOM_error(FATAL, "CFC_BC_FILE must be specified if USE_CFC_CAP=.true.")
  endif
  if (scan(CFC_BC_file, '/') == 0) then
    ! Add the directory if CFC_BC_file is not already a complete path.
    call get_param(param_file, mdl, "INPUTDIR", inputdir, default=".")
    CFC_BC_file = trim(slasher(inputdir))//trim(CFC_BC_file)
    call log_param(param_file, mdl, "INPUTDIR/CFC_BC_FILE", CFC_BC_file, &
                   "full path of CFC_BC_FILE")
  endif

  call get_param(param_file, mdl, "CFC_BC_DATA_YEAR", CFC_BC_data_year, &
                 "Specific year in CFC_BC_FILE data calendar", default=2000)
  call get_param(param_file, mdl, "CFC_BC_MODEL_YEAR", CFC_BC_model_year, &
                 "Model year corresponding to CFC_BC_MODEL_YEAR", default=2000)
  CS%CFC_BC_year_offset = CFC_BC_data_year - CFC_BC_model_year

  call get_param(param_file, mdl, "CFC11_NH_VARIABLE", CFC_BC_var_name, &
                 "Variable name of NH CFC-11 atm mole fraction in CFC_BC_FILE.", &
                 default="cfc11_nh")
  CS%cfc11_atm_nh_handle = init_external_field(CFC_BC_file, CFC_BC_var_name)

  call get_param(param_file, mdl, "CFC11_SH_VARIABLE", CFC_BC_var_name, &
                 "Variable name of SH CFC-11 atm mole fraction in CFC_BC_FILE.", &
                 default="cfc11_sh")
  CS%cfc11_atm_sh_handle = init_external_field(CFC_BC_file, CFC_BC_var_name)

  call get_param(param_file, mdl, "CFC12_NH_VARIABLE", CFC_BC_var_name, &
                 "Variable name of NH CFC-12 atm mole fraction in CFC_BC_FILE.", &
                 default="cfc12_nh")
  CS%cfc12_atm_nh_handle = init_external_field(CFC_BC_file, CFC_BC_var_name)

  call get_param(param_file, mdl, "CFC12_SH_VARIABLE", CFC_BC_var_name, &
                 "Variable name of SH CFC-12 atm mole fraction in CFC_BC_FILE.", &
                 default="cfc12_sh")
  CS%cfc12_atm_sh_handle = init_external_field(CFC_BC_file, CFC_BC_var_name)
!                                              domain=G%Domain%mpp_domain)

  ! The following vardesc types contain a package of metadata about each tracer,
  ! including, the name; units; longname; and grid information.
  do m=1,NTR
    write(m2char, "(I1)") m
    write(CS%CFC_data(m)%name, "(2A)") "CFC_1", m2char
    CS%CFC_data(m)%desc = var_desc(CS%CFC_data(m)%name, &
                                   "mol kg-1", &
                                   "Moles Per Unit Mass of CFC-1"//m2char//" in sea water", &
                                   caller=mdl)

    allocate(CS%CFC_data(m)%conc(isd:ied,jsd:jed,nz), source=0.0)
    allocate(CS%CFC_data(m)%sfc_flux(isd:ied,jsd:jed), source=0.0)

    ! This pointer assignment is needed to force the compiler not to do a copy in
    ! the registration calls.  Curses on the designers and implementers of F90.
    tr_ptr => CS%CFC_data(m)%conc
    ! Register CFC tracer for horizontal advection, diffusion, and restarts.
    call register_tracer(tr_ptr, tr_Reg, param_file, HI, GV, &
                        tr_desc=CS%CFC_data(m)%desc, registry_diags=.true., &
                        restart_CS=restart_CS, mandatory=.not.CS%tracers_may_reinit, &
                        Tr_out=CS%CFC_data(m)%tr_ptr)
  enddo

  CS%tr_Reg => tr_Reg
  CS%restart_CSp => restart_CS
  register_CFC_cap = .true.

end procedure register_CFC_cap
module procedure initialize_CFC_cap
  integer :: m
  character :: m2char
  if (.not.associated(CS)) return

  CS%Time => day
  CS%diag => diag

  do m=1,NTR
    if (.not.restart .or. (CS%tracers_may_reinit .and. &
        .not.query_initialized(CS%CFC_data(m)%conc, CS%CFC_data(m)%name, CS%restart_CSp))) then
      call init_tracer_CFC(h, CS%CFC_data(m)%conc, CS%CFC_data(m)%name, CS%CFC_data(m)%land_val, &
                          CS%CFC_data(m)%IC_val, G, GV, US, CS)
      call set_initialized(CS%CFC_data(m)%conc, CS%CFC_data(m)%name, CS%restart_CSp)
    endif

    ! cmor diagnostics
    ! units for cfc11_flux and cfc12_flux are [Conc R Z T-1 ~> mol m-2 s-1]
    ! CFC11 cmor conventions: http://clipc-services.ceda.ac.uk/dreq/u/42625c97b8fe75124a345962c4430982.html
    ! http://clipc-services.ceda.ac.uk/dreq/u/0940cbee6105037e4b7aa5579004f124.html
    ! CFC12 cmor conventions: http://clipc-services.ceda.ac.uk/dreq/u/3ab8e10027d7014f18f9391890369235.html
    ! http://clipc-services.ceda.ac.uk/dreq/u/e9e21426e4810d0bb2d3dddb24dbf4dc.html
    write(m2char, "(I1)") m
    CS%CFC_data(m)%id_cmor = register_diag_field('ocean_model', &
        'cfc1'//m2char, diag%axesTL, day, &
        'Mole Concentration of CFC1'//m2char//' in Sea Water', 'mol m-3', &
        conversion=GV%Rho0*US%R_to_kg_m3)

    CS%CFC_data(m)%id_sfc_flux = register_diag_field('ocean_model', &
        'cfc1'//m2char//'_flux', diag%axesT1, day, &
        'Gas exchange flux of CFC1'//m2char//' into the ocean ', &
        'mol m-2 s-1', conversion=US%RZ_T_to_kg_m2s, &
        cmor_field_name='fgcfc1'//m2char, &
        cmor_long_name='Surface Downward CFC1'//m2char//' Flux', &
        cmor_standard_name='surface_downward_cfc1'//m2char//'_flux')
  enddo


  if (associated(OBC)) then
  ! Steal from updated DOME in the fullness of time.
  ! GMM: TODO this must be coded if we intend to use this module in regional applications
  endif

end procedure initialize_CFC_cap
module procedure init_tracer_CFC
  logical :: OK
  integer :: i, j, k, is, ie, js, je, nz
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke

  if (len_trim(CS%IC_file) > 0) then
    !  Read the tracer concentrations from a netcdf file.
    if (.not.file_exists(CS%IC_file, G%Domain)) &
      call MOM_error(FATAL, "initialize_CFC_cap: Unable to open "//CS%IC_file)
    if (CS%Z_IC_file) then
      OK = tracer_Z_init(tr, h, CS%IC_file, name, G, GV, US)
      if (.not.OK) then
        OK = tracer_Z_init(tr, h, CS%IC_file, trim(name), G, GV, US)
        if (.not.OK) call MOM_error(FATAL,"initialize_CFC_cap: "//&
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
module procedure CFC_cap_column_physics
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)) :: h_work ! Used so that h can be modified [H ~> m or kg m-2]
  integer :: i, j, k, is, ie, js, je, nz, m
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke

  if (.not.associated(CS)) return

  ! Compute KPP nonlocal term if necessary
  if (present(KPP_CSp)) then
    if (associated(KPP_CSp) .and. present(nonLocalTrans)) then
      do m=1,NTR
        call KPP_NonLocalTransport(KPP_CSp, G, GV, h_old, nonLocalTrans, &
                                   CS%CFC_data(m)%sfc_flux(:,:), dt, CS%diag, &
                                   CS%CFC_data(m)%tr_ptr, CS%CFC_data(m)%conc(:,:,:), &
                                   flux_scale=GV%RZ_to_H)
      enddo
    endif
  endif

  ! Use a tridiagonal solver to determine the concentrations after the
  ! surface source is applied and diapycnal advection and diffusion occurs.
  if (present(evap_CFL_limit) .and. present(minimum_forcing_depth)) then
    do m=1,NTR
      do k=1,nz ;do j=js,je ; do i=is,ie
        h_work(i,j,k) = h_old(i,j,k)
      enddo ; enddo ; enddo
      call applyTracerBoundaryFluxesInOut(G, GV, CS%CFC_data(m)%conc, dt, fluxes, h_work, &
                                          evap_CFL_limit, minimum_forcing_depth)
      call tracer_vertdiff(h_work, ea, eb, dt, CS%CFC_data(m)%conc, G, GV, &
                           sfc_flux=CS%CFC_data(m)%sfc_flux)
    enddo
  else
    do m=1,NTR
      call tracer_vertdiff(h_old, ea, eb, dt, CS%CFC_data(m)%conc, G, GV, &
                           sfc_flux=CS%CFC_data(m)%sfc_flux)
    enddo
  endif

  ! If needed, write out any desired diagnostics from tracer sources & sinks here.
  do m=1,NTR
    if (CS%CFC_data(m)%id_cmor > 0) &
      call post_data(CS%CFC_data(m)%id_cmor, CS%CFC_data(m)%conc, CS%diag)

    if (CS%CFC_data(m)%id_sfc_flux > 0) &
      call post_data(CS%CFC_data(m)%id_sfc_flux, CS%CFC_data(m)%sfc_flux, CS%diag)
  enddo

end procedure CFC_cap_column_physics
module procedure CFC_cap_stock
  real :: stock_scale ! The dimensional scaling factor to convert stocks to kg [kg H-1 L-2 ~> kg m-3 or 1]
  real :: mass        ! The cell volume or mass [H L2 ~> m3 or kg]
  integer :: i, j, k, is, ie, js, je, nz, m
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke

  CFC_cap_stock = 0
  if (.not.associated(CS)) return

  if (present(stock_index)) then ; if (stock_index > 0) then
    ! Check whether this stock is available from this routine.

    ! No stocks from this routine are being checked yet.  Return 0.
    return
  endif ; endif

  do m=1,NTR
    call query_vardesc(CS%CFC_data(m)%desc, name=names(m), units=units(m), caller="CFC_cap_stock")
    units(m) = trim(units(m))//" kg"
    stocks(m) = global_mass_int_EFP(h, G, GV, CS%CFC_data(m)%conc, on_PE_only=.true.)
  enddo

  CFC_cap_stock = NTR

end procedure CFC_cap_stock
module procedure CFC_cap_set_forcing
  type(time_type) :: Time_external ! time value used in CFC_BC_file
  real, dimension(SZI_(G),SZJ_(G)) :: &
    kw_wo_sc_no_term, &  ! gas transfer velocity, without the Schmidt number term [Z T-1 ~> m s-1].
    kw, &                ! gas transfer velocity [Z T-1 ~> m s-1].
    cair, &              ! The surface gas concentration in equilibrium with the atmosphere
                         ! (saturation concentration) [mol kg-1].
    cfc11_atm, &         ! CFC11 atm mole fraction [pico mol/mol]
    cfc12_atm            ! CFC12 atm mole fraction [pico mol/mol]
  real :: cfc11_atm_nh   ! NH value for cfc11_atm [pico mol/mol]
  real :: cfc11_atm_sh   ! SH value for cfc11_atm [pico mol/mol]
  real :: cfc12_atm_nh   ! NH value for cfc12_atm [pico mol/mol]
  real :: cfc12_atm_sh   ! SH value for cfc12_atm [pico mol/mol]
  real :: ta             ! Absolute sea surface temperature [hectoKelvin]
  real :: sal            ! Surface salinity [PSU].
  real :: alpha_11       ! The solubility of CFC 11 [mol kg-1 atm-1].
  real :: alpha_12       ! The solubility of CFC 12 [mol kg-1 atm-1].
  real :: sc_11, sc_12   ! The Schmidt numbers of CFC 11 and CFC 12 [nondim].
  real :: kw_coeff       ! A coefficient used to compute the piston velocity [Z T-1 T2 L-2] = [Z T L-2 ~> s m-1]
  real, parameter :: pa_to_atm = 9.8692316931427e-6 ! factor for converting from Pa to atm [atm Pa-1].
  real :: press_to_atm   ! converts from model pressure units to atm [atm T2 R-1 L-2 ~> atm Pa-1]
  integer :: i, j, is, ie, js, je, m

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec

  ! Time_external = increment_date(day_start + day_interval/2, years=CS%CFC_BC_year_offset)
  Time_external = increment_date(day_start, years=CS%CFC_BC_year_offset)

  ! CFC11 atm mole fraction, convert from ppt (pico mol/mol) to mol/mol
  call time_interp_external(CS%cfc11_atm_nh_handle, Time_external, cfc11_atm_nh)
  cfc11_atm_nh = cfc11_atm_nh * 1.0e-12
  call time_interp_external(CS%cfc11_atm_sh_handle, Time_external, cfc11_atm_sh)
  cfc11_atm_sh = cfc11_atm_sh * 1.0e-12

  ! CFC12 atm mole fraction, convert from ppt (pico mol/mol) to mol/mol
  call time_interp_external(CS%cfc12_atm_nh_handle, Time_external, cfc12_atm_nh)
  cfc12_atm_nh = cfc12_atm_nh * 1.0e-12
  call time_interp_external(CS%cfc12_atm_sh_handle, Time_external, cfc12_atm_sh)
  cfc12_atm_sh = cfc12_atm_sh * 1.0e-12

  !---------------------------------------------------------------------
  !     Gas exchange/piston velocity parameter
  !---------------------------------------------------------------------
  ! From a = 0.251 cm/hr s^2/m^2 in Wannikhof 2014
  !        = 6.97e-7 [m/s s^2/m^2] [Z T-1 T2 L-2] = [Z T L-2 ~> s m-1]
  kw_coeff = (US%m_to_Z*US%s_to_T*US%L_to_m**2) * 6.97e-7

  ! set unit conversion factors
  press_to_atm = US%R_to_kg_m3*US%L_T_to_m_s**2 * pa_to_atm

  do j=js,je ; do i=is,ie
    if (G%geoLatT(i,j) < -10.0) then
      cfc11_atm(i,j) = cfc11_atm_sh
      cfc12_atm(i,j) = cfc12_atm_sh
    elseif (G%geoLatT(i,j) <= 10.0) then
      cfc11_atm(i,j) = cfc11_atm_sh + &
          (0.05 * G%geoLatT(i,j) + 0.5) * (cfc11_atm_nh - cfc11_atm_sh)
      cfc12_atm(i,j) = cfc12_atm_sh + &
          (0.05 * G%geoLatT(i,j) + 0.5) * (cfc12_atm_nh - cfc12_atm_sh)
    else
      cfc11_atm(i,j) = cfc11_atm_nh
      cfc12_atm(i,j) = cfc12_atm_nh
    endif
  enddo ; enddo

  do j=js,je ; do i=is,ie
    ! ta in hectoKelvin
    ta = max(0.01, (US%C_to_degC*sfc_state%SST(i,j) + 273.15) * 0.01)
    sal = US%S_to_ppt*sfc_state%SSS(i,j)

    ! Calculate solubilities
    call get_solubility(alpha_11, alpha_12, ta, sal , G%mask2dT(i,j))

    ! Calculate Schmidt numbers using coefficients given by
    ! Wanninkhof (2014); doi:10.4319/lom.2014.12.351.
    call comp_CFC_schmidt(US%C_to_degC*sfc_state%SST(i,j), sc_11, sc_12)

    kw_wo_sc_no_term(i,j) = kw_coeff * ((1.0 - fluxes%ice_fraction(i,j))*fluxes%u10_sqr(i,j))

    ! air concentrations and cfcs BC's fluxes
    ! CFC flux units: [mol kg-1 R Z T-1 ~> mol m-2 s-1]
    kw(i,j) = kw_wo_sc_no_term(i,j) * sqrt(660.0 / sc_11)
    cair(i,j) = press_to_atm * alpha_11 * cfc11_atm(i,j) * fluxes%p_surf_full(i,j)
    CS%CFC_data(1)%sfc_flux(i,j) = kw(i,j) * (cair(i,j) - CS%CFC_data(1)%conc(i,j,1)) * Rho0

    kw(i,j) = kw_wo_sc_no_term(i,j) * sqrt(660.0 / sc_12)
    cair(i,j) = press_to_atm * alpha_12 * cfc12_atm(i,j) * fluxes%p_surf_full(i,j)
    CS%CFC_data(2)%sfc_flux(i,j) = kw(i,j) * (cair(i,j) - CS%CFC_data(2)%conc(i,j,1)) * Rho0
  enddo ; enddo

  if (CS%debug) then
    do m=1,NTR
      call hchksum(CS%CFC_data(m)%sfc_flux, trim(CS%CFC_data(m)%name)//" sfc_flux", G%HI, &
                   unscale=US%RZ_T_to_kg_m2s)
    enddo
  endif

end procedure CFC_cap_set_forcing
module procedure get_solubility
  real, parameter :: d1_11 = -232.0411    ! [nondim]
  real, parameter :: d2_11 =  322.5546    ! [hectoKelvin-1]
  real, parameter :: d3_11 =  120.4956    ! [log(hectoKelvin)-1]
  real, parameter :: d4_11 =   -1.39165   ! [hectoKelvin-2]
  real, parameter :: e1_11 =   -0.146531  ! [PSU-1]
  real, parameter :: e2_11 =    0.093621  ! [PSU-1 hectoKelvin-1]
  real, parameter :: e3_11 =   -0.0160693 ! [PSU-2 hectoKelvin-2]
  real, parameter :: d1_12 = -220.2120    ! [nondim]
  real, parameter :: d2_12 =  301.8695    ! [hectoKelvin-1]
  real, parameter :: d3_12 =  114.8533    ! [log(hectoKelvin)-1]
  real, parameter :: d4_12 =   -1.39165   ! [hectoKelvin-2]
  real, parameter :: e1_12 =   -0.147718  ! [PSU-1]
  real, parameter :: e2_12 =    0.093175  ! [PSU-1 hectoKelvin-1]
  real, parameter :: e3_12 =   -0.0157340 ! [PSU-2 hectoKelvin-2]
  real :: factor ! introduce units to result [mol kg-1 atm-1]
  factor = 1.0
  alpha_11 = exp(d1_11 + d2_11/ta + d3_11*log(ta) + d4_11*ta**2 +&
                 sal * ((e3_11 * ta + e2_11) * ta + e1_11)) * &
             factor * mask
  alpha_12 = exp(d1_12 + d2_12/ta + d3_12*log(ta) + d4_12*ta**2 +&
                 sal * ((e3_12 * ta + e2_12) * ta + e1_12)) * &
             factor * mask

end procedure get_solubility
module procedure comp_CFC_schmidt
  real , parameter :: a_11 = 3579.2    ! CFC11 Schmidt number fit coefficient [nondim]
  real , parameter :: b_11 = -222.63   ! CFC11 Schmidt number fit coefficient [degC-1]
  real , parameter :: c_11 = 7.5749    ! CFC11 Schmidt number fit coefficient [degC-2]
  real , parameter :: d_11 = -0.14595  ! CFC11 Schmidt number fit coefficient [degC-3]
  real , parameter :: e_11 = 0.0011874 ! CFC11 Schmidt number fit coefficient [degC-4]
  real , parameter :: a_12 = 3828.1    ! CFC12 Schmidt number fit coefficient [nondim]
  real , parameter :: b_12 = -249.86   ! CFC12 Schmidt number fit coefficient [degC-1]
  real , parameter :: c_12 = 8.7603    ! CFC12 Schmidt number fit coefficient [degC-2]
  real , parameter :: d_12 = -0.1716   ! CFC12 Schmidt number fit coefficient [degC-3]
  real , parameter :: e_12 = 0.001408  ! CFC12 Schmidt number fit coefficient [degC-4]
  real             :: sst  ! A range-limited sea surface temperature [degC]
  sst = MAX(-2.0, MIN(40.0, sst_in))
  cfc11_sc = a_11 + sst * (b_11 + sst * (c_11 + sst * (d_11 + sst * e_11)))
  cfc12_sc = a_12 + sst * (b_12 + sst * (c_12 + sst * (d_12 + sst * e_12)))

end procedure comp_CFC_schmidt
module procedure CFC_cap_end
  integer :: m
  if (associated(CS)) then
    do m=1,NTR
      if (associated(CS%CFC_data(m)%conc)) deallocate(CS%CFC_data(m)%conc)
      if (associated(CS%CFC_data(m)%sfc_flux)) deallocate(CS%CFC_data(m)%sfc_flux)
    enddo

    deallocate(CS)
  endif
end procedure CFC_cap_end
module procedure CFC_cap_unit_tests
  real :: dummy1, dummy2 ! Test values of Schmidt numbers [nondim] or solubilities [mol kg-1 atm-1] for CFC11 and CFC12
  real :: ta  ! A test value of temperature [hectoKelvin]
  real :: sal ! A test value of salinity [ppt]
  character(len=120) :: test_name ! Title of the unit test
  CFC_cap_unit_tests = .false.
  write(stdout,*) '==== MOM_CFC_cap ======================='

  ! test comp_CFC_schmidt, Table 1 in Wanninkhof (2014); doi:10.4319/lom.2014.12.351
  test_name = 'Schmidt number calculation'
  call comp_CFC_schmidt(20.0, dummy1, dummy2)
  CFC_cap_unit_tests = CFC_cap_unit_tests .or. &
                       compare_values(verbose, test_name, dummy1, 1179.0, 0.5)
  CFC_cap_unit_tests = CFC_cap_unit_tests .or. &
                       compare_values(verbose, test_name, dummy2, 1188.0, 0.5)

  if (.not. CFC_cap_unit_tests) write(stdout,'(2x,a)') "Passed "//test_name

  test_name = 'Solubility function, SST = 1.0 C, and SSS = 10 psu'
  ta = max(0.01, (1.0 + 273.15) * 0.01) ; sal = 10.
  ! cfc1 = 3.238 10-2 mol kg-1 atm-1
  ! cfc2 = 7.943 10-3 mol kg-1 atm-1
  call get_solubility(dummy1, dummy2, ta, sal , 1.0)
  CFC_cap_unit_tests = CFC_cap_unit_tests .or. &
                       compare_values(verbose, test_name, dummy1, 3.238e-2, 5.0e-6)
  CFC_cap_unit_tests = CFC_cap_unit_tests .or. &
                       compare_values(verbose, test_name, dummy2, 7.943e-3, 5.0e-6)

  if (.not. CFC_cap_unit_tests) write(stdout,'(2x,a)')"Passed "//test_name

  test_name = 'Solubility function, SST = 20.0 C, and SSS = 35 psu'
  ta = max(0.01, (20.0 + 273.15) * 0.01) ; sal = 35.
  ! cfc1 = 0.881 10-2 mol kg-1 atm-1
  ! cfc2 = 2.446 10-3 mol kg-1 atm-1
  call get_solubility(dummy1, dummy2, ta, sal , 1.0)
  CFC_cap_unit_tests = CFC_cap_unit_tests .or. &
                       compare_values(verbose, test_name, dummy1, 8.8145e-3, 5.0e-8)
  CFC_cap_unit_tests = CFC_cap_unit_tests .or. &
                       compare_values(verbose, test_name, dummy2, 2.4462e-3, 5.0e-8)
  if (.not. CFC_cap_unit_tests) write(stdout,'(2x,a)')"Passed "//test_name

end procedure CFC_cap_unit_tests
module procedure compare_values
  real :: diff  ! Difference in values [A]
  diff = ans - calc

  compare_values = .false.
  if (diff > limit ) then
    compare_values = .true.
    write(stdout,*) "CFC_cap_unit_tests, UNIT TEST FAILED: ", test_name
    write(stdout,10) calc, ans
  elseif (verbose) then
    write(stdout,10) calc, ans
  endif

10 format("calc=",f22.16," ans",f22.16)
end procedure compare_values
end submodule MOM_CFC_cap_s
