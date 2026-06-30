submodule (MOM_stoch_eos) MOM_stoch_eos_s
#include <MOM_memory.h>
  implicit none
contains
module procedure MOM_stoch_eos_init
  integer :: i,j
  MOM_stoch_eos_init = .false.

  CS%seed = 0

  call get_param(param_file, "MOM_stoch_eos", "STOCH_EOS", CS%use_stoch_eos, &
                 "If true, stochastic perturbations are applied "//&
                 "to the EOS in the PGF.", default=.false.)
  call get_param(param_file, "MOM_stoch_eos", "STANLEY_COEFF", CS%stanley_coeff, &
                 "Coefficient correlating the temperature gradient "//&
                 "and SGS T variance.", units="nondim", default=-1.0)
  call get_param(param_file, "MOM_stoch_eos", "STANLEY_A", CS%stanley_a, &
                 "Coefficient a which scales chi in stochastic perturbation of the "//&
                 "SGS T variance.", units="nondim", default=1.0, &
                 do_not_log=((CS%stanley_coeff<0.0) .or. .not.CS%use_stoch_eos))
  call get_param(param_file, "MOM_stoch_eos", "KD_SMOOTH", CS%kappa_smooth, &
                 "A diapycnal diffusivity that is used to interpolate "//&
                 "more sensible values of T & S into thin layers.", &
                 units="m2 s-1", default=1.0e-6, scale=GV%m2_s_to_HZ_T, &
                 do_not_log=(CS%stanley_coeff<0.0))

  ! Don't run anything if STANLEY_COEFF < 0
  if (CS%stanley_coeff >= 0.0) then
    if (.not.allocated(CS%pattern)) call MOM_error(FATAL, &
        "MOM_stoch_eos_CS%pattern is not allocated when it should be, suggesting that "//&
        "stoch_EOS_register_restarts() has not been called before MOM_stoch_eos_init().")

    allocate(CS%phi(G%isd:G%ied,G%jsd:G%jed), source=0.0)
    allocate(CS%l2_inv(G%isd:G%ied,G%jsd:G%jed), source=0.0)
    allocate(CS%rgauss(G%isd:G%ied,G%jsd:G%jed), source=0.0)
    call get_param(param_file, "MOM_stoch_eos", "SEED_STOCH_EOS", CS%seed, &
                 "Specfied seed for random number sequence ", default=0)
    call random_2d_constructor(CS%rn_CS, G%HI, Time, CS%seed)
    call random_2d_norm(CS%rn_CS, G%HI, CS%rgauss)
    ! fill array with approximation of grid area needed for decorrelation time-scale calculation
    do j=G%jsc,G%jec
      do i=G%isc,G%iec
        CS%l2_inv(i,j) = 1.0 / ( (G%dxT(i,j)**2) + (G%dyT(i,j)**2) )
      enddo
    enddo

    if (.not.query_initialized(CS%pattern, "stoch_eos_pattern", restart_CS) .or. &
        is_new_run(restart_CS)) then
      do j=G%jsc,G%jec ; do i=G%isc,G%iec
        CS%pattern(i,j) = CS%amplitude*CS%rgauss(i,j)
      enddo ; enddo
    endif

    !register diagnostics
    CS%id_tvar_sgs = register_diag_field('ocean_model', 'tvar_sgs', diag%axesTL, Time, &
      'Parameterized SGS Temperature Variance ', 'None')
    if (CS%use_stoch_eos) then
      CS%id_stoch_eos = register_diag_field('ocean_model', 'stoch_eos', diag%axesT1, Time, &
        'random pattern for EOS', 'None')
      CS%id_stoch_phi = register_diag_field('ocean_model', 'stoch_phi', diag%axesT1, Time, &
        'phi for EOS', 'None')
    endif
  endif

  ! This module is only used if explicitly enabled or a positive correlation coefficient is set.
  MOM_stoch_eos_init = CS%use_stoch_eos .or. (CS%stanley_coeff >= 0.0)

end procedure MOM_stoch_eos_init
module procedure stoch_EOS_register_restarts
  call get_param(param_file, "MOM_stoch_eos", "STANLEY_COEFF", CS%stanley_coeff, &
                 "Coefficient correlating the temperature gradient "//&
                 "and SGS T variance.", units="nondim", default=-1.0, do_not_log=.true.)

  if (CS%stanley_coeff >= 0.0) then
    allocate(CS%pattern(HI%isd:HI%ied,HI%jsd:HI%jed), source=0.0)
    call register_restart_field(CS%pattern, "stoch_eos_pattern", .false., restart_CS, &
                                "Random pattern for stoch EOS", "nondim")
  endif

end procedure stoch_EOS_register_restarts
module procedure MOM_stoch_eos_run
  real    :: ubar, vbar ! Averaged velocities [L T-1 ~> m s-1]
  real    :: phi        ! A temporal correlation factor [nondim]
  integer :: i, j
  if (.not.CS%use_stoch_eos) return

  call random_2d_constructor(CS%rn_CS, G%HI, Time, CS%seed)
  call random_2d_norm(CS%rn_CS, G%HI, CS%rgauss)

  ! advance AR(1)
  do j=G%jsc,G%jec
    do i=G%isc,G%iec
      ubar = 0.5*(u(I,j,1)*G%mask2dCu(I,j)+u(I-1,j,1)*G%mask2dCu(I-1,j))
      vbar = 0.5*(v(i,J,1)*G%mask2dCv(i,J)+v(i,J-1,1)*G%mask2dCv(i,J-1))
      phi = exp(-delt*CS%tfac * sqrt(((ubar**2) + (vbar**2))*CS%l2_inv(i,j)))
      CS%pattern(i,j) = phi*CS%pattern(i,j) + CS%amplitude*sqrt(1-phi**2)*CS%rgauss(i,j)
      CS%phi(i,j) = phi
    enddo
  enddo

end procedure MOM_stoch_eos_run
module procedure post_stoch_EOS_diags
  if (CS%id_stoch_eos > 0) call post_data(CS%id_stoch_eos, CS%pattern, diag)
  if (CS%id_stoch_phi > 0) call post_data(CS%id_stoch_phi, CS%phi, diag)
  if (CS%id_tvar_sgs > 0) call post_data(CS%id_tvar_sgs, tv%varT, diag)

end procedure post_stoch_EOS_diags
module procedure MOM_calc_varT
  real, dimension(SZI_(G), SZJ_(G), SZK_(GV)) :: &
    T, &          !> The temperature (or density) [C ~> degC], with the values in
                  !! in massless layers filled vertically by diffusion.
    S             !> The filled salinity [S ~> ppt], with the values in
                  !! in massless layers filled vertically by diffusion.
  real :: hl(5)              !> Copy of local stencil of H [H ~> m]
  real :: dTdi2, dTdj2       !> Differences in T variance [C2 ~> degC2]
  integer :: i, j, k

  ! Nothing happens if a negative correlation coefficient is set.
  if (CS%stanley_coeff < 0.0) return

  ! This block does a thickness weighted variance calculation and helps control for
  ! extreme gradients along layers which are vanished against topography. It is
  ! still a poor approximation in the interior when coordinates are strongly tilted.
  if (.not. associated(tv%varT)) allocate(tv%varT(G%isd:G%ied, G%jsd:G%jed, GV%ke), source=0.0)
  call vert_fill_TS(h, tv%T, tv%S, CS%kappa_smooth*dt, T, S, G, GV, US, halo_here=1, larger_h_denom=.true.)

  do k=1,G%ke
    do j=G%jsc,G%jec
      do i=G%isc,G%iec
        hl(1) = h(i,j,k) * G%mask2dT(i,j)
        hl(2) = h(i-1,j,k) * G%mask2dCu(I-1,j)
        hl(3) = h(i+1,j,k) * G%mask2dCu(I,j)
        hl(4) = h(i,j-1,k) * G%mask2dCv(i,J-1)
        hl(5) = h(i,j+1,k) * G%mask2dCv(i,J)

        ! SGS variance in i-direction [C2 ~> degC2]
        dTdi2 = ( ( G%mask2dCu(I  ,j) * (G%IdxCu(I  ,j) * ( T(i+1,j,k) - T(i,j,k) )) &
                  + G%mask2dCu(I-1,j) * (G%IdxCu(I-1,j) * ( T(i,j,k) - T(i-1,j,k) )) &
                ) * G%dxT(i,j) * 0.5 )**2
        ! SGS variance in j-direction [C2 ~> degC2]
        dTdj2 = ( ( G%mask2dCv(i,J  ) * (G%IdyCv(i,J  ) * ( T(i,j+1,k) - T(i,j,k) )) &
                  + G%mask2dCv(i,J-1) * (G%IdyCv(i,J-1) * ( T(i,j,k) - T(i,j-1,k) )) &
                ) * G%dyT(i,j) * 0.5 )**2
        tv%varT(i,j,k) = CS%stanley_coeff * ( dTdi2 + dTdj2 )
        ! Turn off scheme near land
        tv%varT(i,j,k) = tv%varT(i,j,k) * (minval(hl) / (maxval(hl) + GV%H_subroundoff))
      enddo
    enddo
  enddo
  ! if stochastic, perturb
  if (CS%use_stoch_eos) then
    do k=1,G%ke
      do j=G%jsc,G%jec
        do i=G%isc,G%iec
          tv%varT(i,j,k) = exp(CS%stanley_a * CS%pattern(i,j)) * tv%varT(i,j,k)
        enddo
      enddo
    enddo
  endif
end procedure MOM_calc_varT
end submodule MOM_stoch_eos_s
