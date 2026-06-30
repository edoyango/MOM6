submodule (MOM_self_attr_load) MOM_self_attr_load_s
#include <MOM_memory.h>
  implicit none
contains
module procedure calc_SAL
  real, dimension(SZI_(G),SZJ_(G)) :: bpa ! SSH or bottom pressure anomaly [Z ~> m] or [R L2 T-2 ~> Pa]
  integer :: n, m, l
  integer :: Isq, Ieq, Jsq, Jeq
  integer :: i, j
  call cpu_clock_begin(id_clock_SAL)

  Isq = G%IscB ; Ieq = G%IecB ; Jsq = G%JscB ; Jeq = G%JecB

  if (CS%use_bpa) then ; do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
    bpa(i,j) = eta(i,j) - CS%pbot_ref(i,j)
  enddo ; enddo ; else ; do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
    bpa(i,j) = eta(i,j)
  enddo ; enddo ; endif

  ! use the scalar approximation and/or iterative tidal SAL
  if (CS%use_sal_scalar .or. CS%use_tidal_sal_prev) then
    do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
      eta_sal(i,j) = CS%linear_scaling * bpa(i,j)
    enddo ; enddo

  ! use the spherical harmonics method
  elseif (CS%use_sal_sht) then
    call spherical_harmonics_forward(G, CS%sht, bpa, CS%Snm_Re, CS%Snm_Im, CS%sal_sht_Nd, tmp_scale=tmp_scale)

    ! Multiply scaling factors to each mode
    do m = 0,CS%sal_sht_Nd
      l = order2index(m, CS%sal_sht_Nd)
      do n = m,CS%sal_sht_Nd
        CS%Snm_Re(l+n-m) = CS%Snm_Re(l+n-m) * CS%Love_scaling(l+n-m)
        CS%Snm_Im(l+n-m) = CS%Snm_Im(l+n-m) * CS%Love_scaling(l+n-m)
      enddo
    enddo

    call spherical_harmonics_inverse(G, CS%sht, CS%Snm_Re, CS%Snm_Im, eta_sal, CS%sal_sht_Nd)
    ! Halo was not calculated in spherical harmonic transforms.
    call pass_var(eta_sal, G%domain)

  else
    do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
      eta_sal(i,j) = 0.0
    enddo ; enddo
  endif

  call cpu_clock_end(id_clock_SAL)
end procedure calc_SAL
module procedure scalar_SAL_sensitivity
  deta_sal_deta = CS%eta_prop
end procedure scalar_SAL_sensitivity
module procedure calc_love_scaling
  real :: coef_rhoE ! A scaling coefficient of solid Earth density. coef_rhoE = rhoW / rhoE with USE_BPA=False
  real, dimension(:), allocatable :: HDat, LDat, KDat ! Love numbers converted in CF reference frames [nondim]
  real :: H1, L1, K1 ! Temporary variables to store degree 1 Love numbers [nondim]
  integer :: n_tot ! Size of the stored Love numbers [nondim]
  integer :: nlm  ! Maximum spherical harmonics degree [nondim]
  integer :: n, m, l
  n_tot = size(Love_Data, dim=2)
  nlm = CS%sal_sht_Nd

  if (nlm+1 > n_tot) call MOM_error(FATAL, "MOM_tidal_forcing " // &
    "calc_love_scaling: maximum spherical harmonics degree is larger than " // &
    "the size of the stored Love numbers in MOM_load_love_number.")

  allocate(HDat(nlm+1), LDat(nlm+1), KDat(nlm+1))
  HDat(:) = Love_Data(2,1:nlm+1) ; LDat(:) = Love_Data(3,1:nlm+1) ; KDat(:) = Love_Data(4,1:nlm+1)

  ! Convert reference frames from CM to CF
  if (nlm > 0) then
    H1 = HDat(2) ; L1 = LDat(2) ;  K1 = KDat(2)
    HDat(2) = ( 2.0 / 3.0) * (H1 - L1)
    LDat(2) = (-1.0 / 3.0) * (H1 - L1)
    KDat(2) = (-1.0 / 3.0) * H1 - (2.0 / 3.0) * L1 - 1.0
  endif

  if (CS%use_bpa) then
    coef_rhoE = 1.0 / (rhoE * grav) ! [Z T2 L-2 R-1 ~> m Pa-1]
  else
    coef_rhoE = rhoW / rhoE ! [nondim]
  endif

  do m=0,nlm ; do n=m,nlm
    l = order2index(m, nlm)
    ! Love_scaling has the same as coef_rhoE.
    CS%Love_scaling(l+n-m) = (3.0 / real(2*n+1)) * coef_rhoE * (1.0 + KDat(n+1) - HDat(n+1))
  enddo ; enddo
end procedure calc_love_scaling
module procedure SAL_init
# include "version_variable.h"
  character(len=40)  :: mdl = "MOM_self_attr_load" ! This module's name.
  integer :: lmax ! Total modes of the real spherical harmonics [nondim]
  real :: rhoE    ! The average density of Earth [R ~> kg m-3].
  character(len=20)  :: bpa_config ! String for reference bottom pressure config option
  real :: tmp(G%isd:G%ied, G%jsd:G%jed) ! Temporary field storing mass returned by find_col_mass
  logical :: restart_sim ! If true, this is a restart run
  character(len=200) :: filename, ref_pbot_file, inputdir ! Strings for file/path
  character(len=200) :: ref_pbot_varname                  ! Variable name in file
  type(MOM_infra_file) :: IO_handle   ! used to write ref_pbot file
  type(vardesc) :: vars(1)            ! used to write ref_pbot file
  type(MOM_field) :: fields(1)        ! used to write ref_pbot file
  logical :: calculate_sal, tides, use_tidal_sal_file
  integer :: default_answer_date, tides_answer_date ! Recover old answers with tides
  real :: sal_scalar_value ! Scaling SAL factors [nondim]
  integer :: isd, ied, jsd, jed
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed

  ! Read all relevant parameters and write them to the model log.
  call log_version(param_file, mdl, version, "")

  call get_param(param_file, '', "TIDES", tides, default=.false., do_not_log=.True.)
  call get_param(param_file, '', "CALCULATE_SAL", calculate_sal, default=tides, do_not_log=.True.)
  if (.not. calculate_sal) return

  call get_param(param_file, mdl, "SAL_USE_BPA", CS%use_bpa, &
                 "If true, use bottom pressure anomaly to calculate self-attraction and "// &
                 "loading (SAL). Otherwise sea surface height anomaly is used, which is "// &
                 "only accurate for uniform density fluid.", default=.False.)
  if (CS%use_bpa) then
    allocate(CS%pbot_ref(isd:ied, jsd:jed), source=0.0)
    call get_param(param_file, mdl, "SAL_REF_PBOT_CONFIG", bpa_config, default="file", &
                   do_not_log=.True.)
    restart_sim = .False. ; if (present(restart_CS)) restart_sim = (.not. is_new_run(restart_CS))
    if (restart_sim .and. (trim(lowercase(bpa_config))/='file')) then
      call MOM_error(WARNING, "SAL_init: 'file' is not used by SAL_PBOT_REF_CONFIG for a restart "//&
                     "run, SAL_PBOT_REF_CONFIG is reset to 'file'.")
      bpa_config = 'file'
    endif
    call get_param(param_file, mdl, "SAL_REF_PBOT_CONFIG", bpa_config, &
                  "A string that determines how the reference bottom pressure for SAL "//&
                  "is specified:\n"//&
                  "\t init - calculated by thickness, temperature and salinity from \n"//&
                  "\t        initialization and assuming surface pressure is zero.\n"//&
                  "\t        This option can only be used by new simulations.\n"//&
                  "\t file - read from the file specified by REF_PBOT_FILE.", &
                  default="file", do_not_read=.True.)
    call get_param(param_file, '', "INPUTDIR", inputdir, default=".", do_not_log=.True.)
    call get_param(param_file, mdl, "REF_PBOT_FILE", ref_pbot_file, &
                   "Reference bottom pressure file used by self-attraction and loading (SAL).", &
                   default="pbot.nc")
    call get_param(param_file, mdl, "REF_PBOT_VARNAME", ref_pbot_varname, &
                   "The name of the variable in REF_PBOT_FILE with reference bottom "//&
                   "pressure.  The variable should have the unit of Pa.", default="pbot")
    filename = trim(slasher(inputdir))//trim(ref_pbot_file)
    call log_param(param_file, mdl, "INPUTDIR/REF_PBOT_FILE", filename)
    select case (trim(lowercase(bpa_config)))
      case ("file")
        call MOM_read_data(filename, trim(ref_pbot_varname), CS%pbot_ref, G%Domain,&
                           scale=US%Pa_to_RL2_T2)
      case ("init")
        call find_col_mass(h, tv, G, GV, US, tmp, CS%pbot_ref)
        ! Write reference bottom pressure file
        vars(1) = var_desc(trim(ref_pbot_varname), units="Pa", &
                           longname="Reference bottom pressure", &
                           hor_grid='h', z_grid='1', t_grid='1')
        call create_MOM_file(IO_handle, trim(filename), vars, 1, fields, G=G)
        call MOM_write_field(IO_handle, fields(1), G%Domain, CS%pbot_ref, unscale=US%RL2_T2_to_Pa)
        call IO_handle%close()
      case default
        call MOM_error(FATAL, "SAL_init: Unsupported SAL_PBOT_REF_CONFIG option "//trim(bpa_config))
    end select
    call pass_var(CS%pbot_ref, G%Domain)
  endif

  call get_param(param_file, mdl, "DEFAULT_ANSWER_DATE", default_answer_date, &
                 "This sets the default value for the various _ANSWER_DATE parameters.", &
                 default=99991231, do_not_log=.True.) ! used to check SAL_USE_BPA
  call get_param(param_file, '', "TIDES_ANSWER_DATE", tides_answer_date, &
                 default=default_answer_date, do_not_log=.True.) ! used to check SAL_USE_BPA
  if (tides_answer_date<=20250131 .and. CS%use_bpa) &
    call MOM_error(FATAL, trim(mdl) // ", SAL_init: SAL_USE_BPA needs to be false to recover "//&
                   "tide answers before 20250131.")
  call get_param(param_file, '', "TIDAL_SAL_FROM_FILE", use_tidal_sal_file, default=.false., &
                 do_not_log=.True.) ! used to set default of SAL_SCALAR_APPROX
  call get_param(param_file, mdl, "SAL_SCALAR_APPROX", CS%use_sal_scalar, &
                 "If true, use the scalar approximation to calculate self-attraction and "//&
                 "loading.", default=tides .and. (.not. use_tidal_sal_file))
  if (CS%use_sal_scalar .and. CS%use_bpa) &
    call MOM_error(WARNING, trim(mdl) // ", SAL_init: Using bottom pressure anomaly for scalar "//&
                   "approximation SAL is unsubstantiated.")
  call get_param(param_file, mdl, "SAL_SCALAR_VALUE", sal_scalar_value, "The constant of "//&
                 "proportionality between self-attraction and loading (SAL) geopotential "//&
                 "anomaly and barotropic geopotential anomaly. This is only used if "//&
                 "SAL_SCALAR_APPROX is true or USE_PREVIOUS_TIDES is true.", default=0.0, &
                 units="m m-1", do_not_log=.not.(CS%use_sal_scalar .or. CS%use_tidal_sal_prev), &
                 old_name='TIDE_SAL_SCALAR_VALUE')
  call get_param(param_file, '', "USE_PREVIOUS_TIDES", CS%use_tidal_sal_prev, &
                 default=.false., do_not_log=.True.)
  call get_param(param_file, mdl, "SAL_HARMONICS", CS%use_sal_sht, &
                 "If true, use the online spherical harmonics method to calculate "//&
                 "self-attraction and loading.", default=.false.)
  call get_param(param_file, mdl, "SAL_HARMONICS_DEGREE", CS%sal_sht_Nd, &
                 "The maximum degree of the spherical harmonics transformation used for "// &
                 "calculating the self-attraction and loading term.", &
                 default=0, do_not_log=(.not. CS%use_sal_sht))
  call get_param(param_file, mdl, "RHO_SOLID_EARTH", rhoE, &
                 "The mean solid earth density.  This is used for calculating the "// &
                 "self-attraction and loading term.", units="kg m-3", &
                 default=5517.0, scale=US%kg_m3_to_R, do_not_log=(.not. CS%use_sal_sht))

  ! Set scaling coefficients for scalar approximation
  if (CS%use_sal_scalar .or. CS%use_tidal_sal_prev) then
    if (CS%use_sal_scalar .and. CS%use_tidal_sal_prev) then
      CS%eta_prop = 2.0 * sal_scalar_value
    else
      CS%eta_prop = sal_scalar_value
    endif
    if (CS%use_bpa) then
      CS%linear_scaling = CS%eta_prop / (GV%Rho0 * GV%g_Earth)
    else
      CS%linear_scaling = CS%eta_prop
    endif
  else
    CS%eta_prop = 0.0 ; CS%linear_scaling = 0.0
  endif

  ! Set scaling coefficients for spherical harmonics
  if (CS%use_sal_sht) then
    lmax = calc_lmax(CS%sal_sht_Nd)
    allocate(CS%Snm_Re(lmax), source=0.0)
    allocate(CS%Snm_Im(lmax), source=0.0)

    allocate(CS%Love_scaling(lmax), source=0.0)
    call calc_love_scaling(GV%Rho0, rhoE, GV%g_Earth, CS)

    allocate(CS%sht)
    call spherical_harmonics_init(G, param_file, CS%sht)
  endif

  id_clock_SAL = cpu_clock_id('(Ocean SAL)', grain=CLOCK_MODULE)

end procedure SAL_init
module procedure SAL_end
  if (allocated(CS%pbot_ref)) deallocate(CS%pbot_ref)

  if (CS%use_sal_sht) then
    if (allocated(CS%Love_scaling)) deallocate(CS%Love_scaling)
    if (allocated(CS%Snm_Re)) deallocate(CS%Snm_Re)
    if (allocated(CS%Snm_Im)) deallocate(CS%Snm_Im)
    call spherical_harmonics_end(CS%sht)
    deallocate(CS%sht)
  endif
end procedure SAL_end
end submodule MOM_self_attr_load_s
