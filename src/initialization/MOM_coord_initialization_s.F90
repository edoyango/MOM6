submodule (MOM_coord_initialization) MOM_coord_initialization_s
  implicit none
contains
module procedure MOM_initialize_coord
  character(len=200) :: config
  logical :: debug
# include "version_variable.h"
  integer :: nz
  nz = GV%ke

  call callTree_enter("MOM_initialize_coord(), MOM_coord_initialization.F90")
  call log_version(PF, mdl, version, "")
  call get_param(PF, mdl, "DEBUG", debug, default=.false.)

! Set-up the layer densities, GV%Rlay, and reduced gravities, GV%g_prime.
  call get_param(PF, mdl, "COORD_CONFIG", config, &
                 "This specifies how layers are to be defined: \n"//&
                 " \t ALE or none - used to avoid defining layers in ALE mode \n"//&
                 " \t file - read coordinate information from the file \n"//&
                 " \t\t specified by (COORD_FILE).\n"//&
                 " \t BFB - Custom coords for buoyancy-forced basin case \n"//&
                 " \t\t based on SST_S, T_BOT and DRHO_DT.\n"//&
                 " \t linear - linear based on interfaces not layers \n"//&
                 " \t layer_ref - linear based on layer densities \n"//&
                 " \t ts_ref - use reference temperature and salinity \n"//&
                 " \t ts_range - use range of temperature and salinity \n"//&
                 " \t\t (T_REF and S_REF) to determine surface density \n"//&
                 " \t\t and GINT calculate internal densities. \n"//&
                 " \t gprime - use reference density (RHO_0) for surface \n"//&
                 " \t\t density and GINT calculate internal densities. \n"//&
                 " \t ts_profile - use temperature and salinity profiles \n"//&
                 " \t\t (read from COORD_FILE) to set layer densities. \n"//&
                 " \t USER - call a user modified routine.", &
                 default="none")
  select case ( trim(config) )
    case ("gprime")
      call set_coord_from_gprime(GV%Rlay, GV%g_prime, GV, US, PF)
    case ("layer_ref")
      call set_coord_from_layer_density(GV%Rlay, GV%g_prime, GV, US, PF)
    case ("linear")
      call set_coord_linear(GV%Rlay, GV%g_prime, GV, US, PF)
    case ("ts_ref")
      call set_coord_from_TS_ref(GV%Rlay, GV%g_prime, GV, US, PF, tv%eqn_of_state, tv%P_Ref)
    case ("ts_profile")
      call set_coord_from_TS_profile(GV%Rlay, GV%g_prime, GV, US, PF, tv%eqn_of_state, tv%P_Ref)
    case ("ts_range")
      call set_coord_from_TS_range(GV%Rlay, GV%g_prime, GV, US, PF, tv%eqn_of_state, tv%P_Ref)
    case ("file")
      call set_coord_from_file(GV%Rlay, GV%g_prime, GV, US, PF)
    case ("USER")
      call user_set_coord(GV%Rlay, GV%g_prime, GV, US, PF)
    case ("BFB")
      call BFB_set_coord(GV%Rlay, GV%g_prime, GV, US, PF)
    case ("none", "ALE")
      call set_coord_to_none(GV%Rlay, GV%g_prime, GV, US, PF)
    case default ; call MOM_error(FATAL,"MOM_initialize_coord: "// &
      "Unrecognized coordinate setup"//trim(config))
  end select
  ! There are nz+1 values of g_prime because it is an interface field, but the value at the bottom
  ! should not matter.  This is here just to avoid having an uninitialized value in some output.
  GV%g_prime(nz+1) = 10.0*GV%g_Earth

  if (debug) call chksum(US%R_to_kg_m3*GV%Rlay(:), "MOM_initialize_coord: Rlay ", 1, nz)
  if (debug) call chksum(US%m_to_Z*US%L_to_m**2*US%s_to_T**2*GV%g_prime(:), "MOM_initialize_coord: g_prime ", 1, nz)
  call setVerticalGridAxes( GV%Rlay, GV, scale=US%R_to_kg_m3 )

  ! Copy the maximum depth across from the input argument
  GV%max_depth = max_depth

  call callTree_leave('MOM_initialize_coord()')

end procedure MOM_initialize_coord
module procedure set_coord_from_gprime
  real :: g_int   ! Reduced gravities across the internal interfaces [L2 Z-1 T-2 ~> m s-2].
  real :: g_fs    ! Reduced gravity across the free surface [L2 Z-1 T-2 ~> m s-2].
  real :: Rlay_Ref ! The target density of the surface layer [R ~> kg m-3].
  character(len=40)  :: mdl = "set_coord_from_gprime" ! This subroutine's name.
  integer :: k, nz
  nz = GV%ke

  call callTree_enter(trim(mdl)//"(), MOM_coord_initialization.F90")

  call get_param(param_file, mdl, "GFS" , g_fs, &
                 "The reduced gravity at the free surface.", units="m s-2", &
                 default=GV%g_Earth*US%L_T_to_m_s**2*US%m_to_Z, scale=US%m_s_to_L_T**2*US%Z_to_m)
  call get_param(param_file, mdl, "GINT", g_int, &
                 "The reduced gravity across internal interfaces.", &
                 units="m s-2", fail_if_missing=.true., scale=US%m_s_to_L_T**2*US%Z_to_m)
  call get_param(param_file, mdl, "LIGHTEST_DENSITY", Rlay_Ref, &
                 "The reference potential density used for layer 1.", &
                 units="kg m-3", default=US%R_to_kg_m3*GV%Rho0, scale=US%kg_m3_to_R)

  g_prime(1) = g_fs
  do k=2,nz ; g_prime(k) = g_int ; enddo
  Rlay(1) = Rlay_Ref
  if (GV%Boussinesq .or. GV%semi_Boussinesq) then
    do k=2,nz ; Rlay(k) = Rlay(k-1) + g_prime(k)*(GV%Rho0/GV%g_Earth) ; enddo
  else
    do k=2,nz
      Rlay(k) = Rlay(k-1) * ((GV%g_Earth + 0.5*g_prime(k)) / (GV%g_Earth - 0.5*g_prime(k)))
    enddo
  endif

  call callTree_leave(trim(mdl)//'()')

end procedure set_coord_from_gprime
module procedure set_coord_from_layer_density
  real :: g_fs    ! Reduced gravity across the free surface [L2 Z-1 T-2 ~> m s-2].
  real :: Rlay_Ref! The surface layer's target density [R ~> kg m-3].
  real :: RLay_range ! The range of densities [R ~> kg m-3].
  character(len=40)  :: mdl = "set_coord_from_layer_density" ! This subroutine's name.
  integer :: k, nz
  nz = GV%ke

  call callTree_enter(trim(mdl)//"(), MOM_coord_initialization.F90")

  call get_param(param_file, mdl, "GFS", g_fs, &
                 "The reduced gravity at the free surface.", units="m s-2", &
                 default=GV%g_Earth*US%L_T_to_m_s**2*US%m_to_Z, scale=US%m_s_to_L_T**2*US%Z_to_m)
  call get_param(param_file, mdl, "LIGHTEST_DENSITY", Rlay_Ref, &
                 "The reference potential density used for layer 1.", &
                 units="kg m-3", default=US%R_to_kg_m3*GV%Rho0, scale=US%kg_m3_to_R)
  call get_param(param_file, mdl, "DENSITY_RANGE", Rlay_range, &
                 "The range of reference potential densities in the layers.", &
                 units="kg m-3", default=2.0, scale=US%kg_m3_to_R)

  Rlay(1) = Rlay_Ref
  do k=2,nz
    Rlay(k) = Rlay(k-1) + RLay_range/(real(nz-1))
  enddo
!    These statements set the interface reduced gravities.           !
  g_prime(1) = g_fs
  if (GV%Boussinesq .or. GV%semi_Boussinesq) then
    do k=2,nz
      g_prime(k) = (GV%g_Earth/GV%Rho0) * (Rlay(k) - Rlay(k-1))
    enddo
  else
    do k=2,nz
      g_prime(k) = 2.0*GV%g_Earth * (Rlay(k) - Rlay(k-1)) / (Rlay(k) + Rlay(k-1))
    enddo
  endif

  call callTree_leave(trim(mdl)//'()')
end procedure set_coord_from_layer_density
module procedure set_coord_from_TS_ref
  real :: T_ref   ! Reference temperature [C ~> degC]
  real :: S_ref   ! Reference salinity [S ~> ppt]
  real :: g_int   ! Reduced gravities across the internal interfaces [L2 Z-1 T-2 ~> m s-2].
  real :: g_fs    ! Reduced gravity across the free surface [L2 Z-1 T-2 ~> m s-2].
  character(len=40)  :: mdl = "set_coord_from_TS_ref" ! This subroutine's name.
  integer :: k, nz
  nz = GV%ke

  call callTree_enter(trim(mdl)//"(), MOM_coord_initialization.F90")

  call get_param(param_file, mdl, "T_REF", T_ref, &
                 "The initial temperature of the lightest layer.", &
                 units="degC", scale=US%degC_to_C, fail_if_missing=.true.)
  call get_param(param_file, mdl, "S_REF", S_ref, &
                 "The initial salinities.", units="ppt", default=35.0, scale=US%ppt_to_S)
  call get_param(param_file, mdl, "GFS", g_fs, &
                 "The reduced gravity at the free surface.", units="m s-2", &
                 default=GV%g_Earth*US%L_T_to_m_s**2*US%m_to_Z, scale=US%m_s_to_L_T**2*US%Z_to_m)
  call get_param(param_file, mdl, "GINT", g_int, &
                 "The reduced gravity across internal interfaces.", &
                 units="m s-2", fail_if_missing=.true., scale=US%m_s_to_L_T**2*US%Z_to_m)

!    These statements set the interface reduced gravities.           !
  g_prime(1) = g_fs
  do k=2,nz ; g_prime(k) = g_int ; enddo

!    The uppermost layer's density is set here.  Subsequent layers'  !
!  densities are determined from this value and the g values.        !
!        T0 = 28.228 ; S0 = 34.5848 ; Pref = P_Ref
  call calculate_density(T_ref, S_ref, P_ref, Rlay(1), eqn_of_state)

!    These statements set the layer densities.                       !
  if (GV%Boussinesq .or. GV%semi_Boussinesq) then
    do k=2,nz ; Rlay(k) = Rlay(k-1) + g_prime(k)*(GV%Rho0/GV%g_Earth) ; enddo
  else
    do k=2,nz
      Rlay(k) = Rlay(k-1) * ((GV%g_Earth + 0.5*g_prime(k)) / (GV%g_Earth - 0.5*g_prime(k)))
    enddo
  endif

  call callTree_leave(trim(mdl)//'()')
end procedure set_coord_from_TS_ref
module procedure set_coord_from_TS_profile
  real, dimension(GV%ke) :: T0   ! A profile of temperatures [C ~> degC]
  real, dimension(GV%ke) :: S0   ! A profile of salinities [S ~> ppt]
  real, dimension(GV%ke) :: Pref ! A array of reference pressures [R L2 T-2 ~> Pa]
  real :: g_fs    ! Reduced gravity across the free surface [L2 Z-1 T-2 ~> m s-2].
  integer :: k, nz
  character(len=40)  :: mdl = "set_coord_from_TS_profile" ! This subroutine's name.
  character(len=200) :: filename, coord_file, inputdir ! Strings for file/path
  character(len=64)  :: temp_var, salt_var ! Temperature and salinity names in files
  nz = GV%ke

  call callTree_enter(trim(mdl)//"(), MOM_coord_initialization.F90")

  call get_param(param_file, mdl, "GFS", g_fs, &
                 "The reduced gravity at the free surface.", units="m s-2", &
                 default=GV%g_Earth*US%L_T_to_m_s**2*US%m_to_Z, scale=US%m_s_to_L_T**2*US%Z_to_m)
  call get_param(param_file, mdl, "COORD_FILE", coord_file, &
                 "The file from which the coordinate temperatures and salinities are read.", &
                 fail_if_missing=.true.)
  call get_param(param_file, mdl, "TEMP_COORD_VAR", temp_var, &
                 "The coordinate reference profile variable name for potential temperature.", &
                 default="PTEMP")
  call get_param(param_file, mdl, "SALT_COORD_VAR", salt_var, &
                 "The coordinate reference profile variable name for salinity.", &
                 default="SALT")

  call get_param(param_file,  mdl, "INPUTDIR", inputdir, default=".")
  filename = trim(slasher(inputdir))//trim(coord_file)
  call log_param(param_file, mdl, "INPUTDIR/COORD_FILE", filename)

  call MOM_read_data(filename, temp_var, T0(:), scale=US%degC_to_C)
  call MOM_read_data(filename, salt_var, S0(:), scale=US%ppt_to_S)

  if (.not.file_exists(filename)) call MOM_error(FATAL, &
      " set_coord_from_TS_profile: Unable to open " //trim(filename))
!    These statements set the interface reduced gravities.           !
  g_prime(1) = g_fs
  do k=1,nz ; Pref(k) = P_Ref ; enddo
  call calculate_density(T0, S0, Pref, Rlay, eqn_of_state, (/1,nz/) )
  if (GV%Boussinesq .or. GV%semi_Boussinesq) then
    do k=2,nz
      g_prime(k) = (GV%g_Earth/GV%Rho0) * (Rlay(k) - Rlay(k-1))
    enddo
  else
    do k=2,nz
      g_prime(k) = 2.0*GV%g_Earth * (Rlay(k) - Rlay(k-1)) / (Rlay(k) + Rlay(k-1))
    enddo
  endif

  call callTree_leave(trim(mdl)//'()')
end procedure set_coord_from_TS_profile
module procedure set_coord_from_TS_range
  real, dimension(GV%ke) :: T0   ! A profile of temperatures [C ~> degC]
  real, dimension(GV%ke) :: S0   ! A profile of salinities [S ~> ppt]
  real, dimension(GV%ke) :: Pref ! A array of reference pressures [R L2 T-2 ~> Pa]
  real :: S_Ref   ! Default salinity range parameters [S ~> ppt].
  real :: T_Ref   ! Default temperature range parameters [C ~> degC].
  real :: S_Light, S_Dense ! Salinity range parameters [S ~> ppt].
  real :: T_Light, T_Dense ! Temperature range parameters [C ~> degC].
  real :: res_rat ! The ratio of density space resolution in the denser part
  real :: g_fs    ! Reduced gravity across the free surface [L2 Z-1 T-2 ~> m s-2].
  real :: a1, frac_dense, k_frac  ! Nondimensional temporary variables [nondim]
  integer :: k, nz, k_light
  character(len=40)  :: mdl = "set_coord_from_TS_range" ! This subroutine's name.
  nz = GV%ke

  call callTree_enter(trim(mdl)//"(), MOM_coord_initialization.F90")

  call get_param(param_file, mdl, "T_REF", T_Ref, &
                 "The default initial temperatures.", &
                 units="degC", default=10.0, scale=US%degC_to_C)
  call get_param(param_file, mdl, "TS_RANGE_T_LIGHT", T_Light, &
                 "The initial temperature of the lightest layer when "//&
                 "COORD_CONFIG is set to ts_range.", &
                 units="degC", default=US%C_to_degC*T_Ref, scale=US%degC_to_C)
  call get_param(param_file, mdl, "TS_RANGE_T_DENSE", T_Dense, &
                 "The initial temperature of the densest layer when "//&
                 "COORD_CONFIG is set to ts_range.", &
                 units="degC", default=US%C_to_degC*T_Ref, scale=US%degC_to_C)

  call get_param(param_file, mdl, "S_REF", S_Ref, &
                 "The default initial salinities.", &
                 units="ppt", default=35.0, scale=US%ppt_to_S)
  call get_param(param_file, mdl, "TS_RANGE_S_LIGHT", S_Light, &
                 "The initial lightest salinities when COORD_CONFIG is set to ts_range.", &
                 units="ppt", default=US%S_to_ppt*S_Ref, scale=US%ppt_to_S)
  call get_param(param_file, mdl, "TS_RANGE_S_DENSE", S_Dense, &
                 "The initial densest salinities when COORD_CONFIG is set to ts_range.", &
                 units="ppt", default=US%S_to_ppt*S_Ref, scale=US%ppt_to_S)

  call get_param(param_file, mdl, "TS_RANGE_RESOLN_RATIO", res_rat, &
                 "The ratio of density space resolution in the densest "//&
                 "part of the range to that in the lightest part of the "//&
                 "range when COORD_CONFIG is set to ts_range. Values "//&
                 "greater than 1 increase the resolution of the denser water.",&
                 default=1.0, units="nondim")

  call get_param(param_file, mdl, "GFS", g_fs, &
                 "The reduced gravity at the free surface.", units="m s-2", &
                 default=GV%g_Earth*US%L_T_to_m_s**2*US%m_to_Z, scale=US%m_s_to_L_T**2*US%Z_to_m)

  if ((GV%nk_rho_varies > 0) .and. (nz < GV%nk_rho_varies+2)) &
    call MOM_error(FATAL, "set_coord_from_TS_range requires that NZ >= NKML+NKBL+2.")

  k_light = GV%nk_rho_varies + 1

  ! Set T0(k) to range from T_LIGHT to T_DENSE, and similarly for S0(k).
  T0(k_light) = T_Light ; S0(k_light) = S_Light
  a1 = 2.0 * res_rat / (1.0 + res_rat)
  do k=k_light+1,nz
    k_frac = real(k-k_light)/real(nz-k_light)
    frac_dense = a1 * k_frac + (1.0 - a1) * k_frac**2
    T0(k) = frac_dense * (T_Dense - T_Light) + T_Light
    S0(k) = frac_dense * (S_Dense - S_Light) + S_Light
  enddo

  g_prime(1) = g_fs
  do k=1,nz ; Pref(k) = P_Ref ; enddo
  call calculate_density(T0, S0, Pref, Rlay, eqn_of_state, (/k_light,nz/) )
  ! Extrapolate target densities for the variable density mixed and buffer layers.
  do k=k_light-1,1,-1
    Rlay(k) = 2.0*Rlay(k+1) - Rlay(k+2)
  enddo
  if (GV%Boussinesq .or. GV%semi_Boussinesq) then
    do k=2,nz
      g_prime(k) = (GV%g_Earth/GV%Rho0) * (Rlay(k) - Rlay(k-1))
    enddo
  else
    do k=2,nz
      g_prime(k) = 2.0*GV%g_Earth * (Rlay(k) - Rlay(k-1)) / (Rlay(k) + Rlay(k-1))
    enddo
  endif

  call callTree_leave(trim(mdl)//'()')
end procedure set_coord_from_TS_range
module procedure set_coord_from_file
  real :: g_fs    ! Reduced gravity across the free surface [L2 Z-1 T-2 ~> m s-2].
  integer :: k, nz
  character(len=40)  :: mdl = "set_coord_from_file" ! This subroutine's name.
  character(len=40)  :: coord_var
  character(len=200) :: filename,coord_file,inputdir ! Strings for file/path
  nz = GV%ke

  call callTree_enter(trim(mdl)//"(), MOM_coord_initialization.F90")

  call get_param(param_file, mdl, "GFS", g_fs, &
                 "The reduced gravity at the free surface.", units="m s-2", &
                 default=GV%g_Earth*US%L_T_to_m_s**2*US%m_to_Z, scale=US%m_s_to_L_T**2*US%Z_to_m)
  call get_param(param_file, mdl, "INPUTDIR", inputdir, default=".")
  inputdir = slasher(inputdir)
  call get_param(param_file, mdl, "COORD_FILE", coord_file, &
                 "The file from which the coordinate densities are read.", &
                 fail_if_missing=.true.)
  call get_param(param_file, mdl, "COORD_VAR", coord_var, &
                 "The variable in COORD_FILE that is to be used for the "//&
                 "coordinate densities.", default="Layer")
  filename = trim(inputdir)//trim(coord_file)
  call log_param(param_file, mdl, "INPUTDIR/COORD_FILE", filename)
  if (.not.file_exists(filename)) call MOM_error(FATAL, &
      " set_coord_from_file: Unable to open "//trim(filename))

  call MOM_read_data(filename, coord_var, Rlay, scale=US%kg_m3_to_R)
  g_prime(1) = g_fs
  if (GV%Boussinesq .or. GV%semi_Boussinesq) then
    do k=2,nz
      g_prime(k) = (GV%g_Earth/GV%Rho0) * (Rlay(k) - Rlay(k-1))
    enddo
  else
    do k=2,nz
      g_prime(k) = 2.0*GV%g_Earth * (Rlay(k) - Rlay(k-1)) / (Rlay(k) + Rlay(k-1))
    enddo
  endif
  do k=1,nz ; if (g_prime(k) <= 0.0) then
    call MOM_error(FATAL, "MOM_initialization set_coord_from_file: "//&
       "Zero or negative g_primes read from variable "//"Layer"//" in file "//&
       trim(filename))
  endif ; enddo

  call callTree_leave(trim(mdl)//'()')
end procedure set_coord_from_file
module procedure set_coord_linear
  character(len=40)  :: mdl = "set_coord_linear" ! This subroutine
  real :: Rlay_ref, Rlay_range ! A reference density and its range [R ~> kg m-3]
  real :: g_fs  ! The reduced gravity across the free surface [L2 Z-1 T-2 ~> m s-2]
  integer :: k, nz
  nz = GV%ke

  call callTree_enter(trim(mdl)//"(), MOM_coord_initialization.F90")

  call get_param(param_file, mdl, "LIGHTEST_DENSITY", Rlay_Ref, &
                 "The reference potential density used for the surface interface.", &
                 units="kg m-3", default=US%R_to_kg_m3*GV%Rho0, scale=US%kg_m3_to_R)
  call get_param(param_file, mdl, "DENSITY_RANGE", Rlay_range, &
                 "The range of reference potential densities across all interfaces.", &
                 units="kg m-3", default=2.0, scale=US%kg_m3_to_R)
  call get_param(param_file, mdl, "GFS", g_fs, &
                 "The reduced gravity at the free surface.", units="m s-2", &
                 default=GV%g_Earth*US%L_T_to_m_s**2*US%m_to_Z, scale=US%m_s_to_L_T**2*US%Z_to_m)

  ! This following sets the target layer densities such that the
  ! surface interface has density Rlay_ref and the bottom
  ! is Rlay_range larger
  do k=1,nz
    Rlay(k) = Rlay_Ref + RLay_range*((real(k)-0.5)/real(nz))
  enddo
  ! These statements set the interface reduced gravities.
  g_prime(1) = g_fs
  if (GV%Boussinesq .or. GV%semi_Boussinesq) then
    do k=2,nz
      g_prime(k) = (GV%g_Earth/GV%Rho0) * (Rlay(k) - Rlay(k-1))
    enddo
  else
    do k=2,nz
      g_prime(k) = 2.0*GV%g_Earth * (Rlay(k) - Rlay(k-1)) / (Rlay(k) + Rlay(k-1))
    enddo
  endif

  call callTree_leave(trim(mdl)//'()')
end procedure set_coord_linear
module procedure set_coord_to_none
  real :: g_fs    ! Reduced gravity across the free surface [L2 Z-1 T-2 ~> m s-2].
  real :: Rlay_Ref ! The target density of the surface layer [R ~> kg m-3].
  character(len=40)  :: mdl = "set_coord_to_none" ! This subroutine's name.
  integer :: k, nz
  nz = GV%ke

  call callTree_enter(trim(mdl)//"(), MOM_coord_initialization.F90")

  call get_param(param_file, mdl, "GFS" , g_fs, &
                 "The reduced gravity at the free surface.", units="m s-2", &
                 default=GV%g_Earth*US%L_T_to_m_s**2*US%m_to_Z, scale=US%m_s_to_L_T**2*US%Z_to_m)
  call get_param(param_file, mdl, "LIGHTEST_DENSITY", Rlay_Ref, &
                 "The reference potential density used for layer 1.", &
                 units="kg m-3", default=US%R_to_kg_m3*GV%Rho0, scale=US%kg_m3_to_R)

  g_prime(1) = g_fs
  do k=2,nz ; g_prime(k) = 0. ; enddo
  Rlay(1) = Rlay_Ref
  if (GV%Boussinesq .or. GV%semi_Boussinesq) then
    do k=2,nz ; Rlay(k) = Rlay(k-1) + g_prime(k)*(GV%Rho0/GV%g_Earth) ; enddo
  else
    do k=2,nz
      Rlay(k) = Rlay(k-1) * ((GV%g_Earth + 0.5*g_prime(k)) / (GV%g_Earth - 0.5*g_prime(k)))
    enddo
  endif

  call callTree_leave(trim(mdl)//'()')

end procedure set_coord_to_none
module procedure write_vertgrid_file
  character(len=240) :: filepath
  type(vardesc) :: vars(2)
  type(MOM_field) :: fields(2)
  type(MOM_netCDF_file) :: IO_handle ! The I/O handle of the fileset
  filepath = trim(directory) // trim("Vertical_coordinate.nc")

  vars(1) = var_desc("R","kilogram meter-3","Target Potential Density",'1','L','1')
  vars(2) = var_desc("g","meter second-2","Reduced gravity",'1','i','1')

  call create_MOM_file(IO_handle, trim(filepath), vars, 2, fields, &
      SINGLE_FILE, GV=GV)

  call MOM_write_field(IO_handle, fields(1), GV%Rlay, unscale=US%R_to_kg_m3)
  call MOM_write_field(IO_handle, fields(2), GV%g_prime, unscale=US%L_T_to_m_s**2*US%m_to_Z)

  call IO_handle%close()

end procedure write_vertgrid_file
end submodule MOM_coord_initialization_s
