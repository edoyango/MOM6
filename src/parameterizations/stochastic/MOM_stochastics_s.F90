submodule (MOM_stochastics) MOM_stochastics_s
#include <MOM_memory.h>
  implicit none
contains
module procedure stochastics_init
  integer, allocatable :: pelist(:) ! list of pes for this instance of the ocean
  integer :: mom_comm          ! list of pes for this instance of the ocean
  integer :: num_procs         ! number of processors to pass to stochastic physics
  integer :: iret              ! return code from stochastic physics
  integer :: pe_zero           !  root pe
  integer :: nxT, nxB          ! number of x-points including halo
  integer :: nyT, nyB          ! number of y-points including halo
  integer :: i, j, k           ! loop indices
  real    :: tmp(grid%isdB:grid%iedB,grid%jsdB:grid%jedB) ! Used to construct tapers
  integer :: taper_width       ! Width (in cells) of the taper that brings the stochastic velocity
# include "version_variable.h"
  character(len=40)  :: mdl = "ocean_stochastics_init"  ! This module's name.
  call callTree_enter("stochastic_init(), MOM_stochastics.F90")
  if (associated(CS)) then
    call MOM_error(WARNING, "MOM_stochastics_init called with an "// &
                            "associated control structure.")
    return
  else ; allocate(CS) ; endif

  CS%Time => Time
  CS%diag => diag

  ! Read all relevant parameters and write them to the model log.
  call log_version(param_file, mdl, version, "")

  ! get number of processors and PE list for stochastic physics initialization
  call get_param(param_file, mdl, "DO_SPPT", CS%do_sppt, &
                 "If true, then stochastically perturb the thermodynamic "//&
                 "tendencies of T,S, and h.  Amplitude and correlations are "//&
                 "controlled by the nam_stoch namelist in the UFS model only.", &
                 default=.false.)
  call get_param(param_file, mdl, "DO_SKEB", CS%do_skeb, &
                 "If true, then stochastically perturb the currents "//&
                 "using the stochastic kinetic energy backscatter scheme.",&
                 default=.false.)
  call get_param(param_file, mdl, "SKEB_NPASS", CS%skeb_npass, &
                 "number of passes of a 9-point smoother of the "//&
                 "dissipation estimate.", default=3, do_not_log=.not.CS%do_skeb)
  call get_param(param_file, mdl, "SKEB_TAPER_WIDTH", taper_width, &
                 "number of cells over which the stochastic velocity increment "//&
                 "is tapered to zero.", default=4, do_not_log=.not.CS%do_skeb)
  call get_param(param_file, mdl, "SKEB_USE_GM", CS%skeb_use_gm, &
                 "If true, adds GM work rate to the SKEBS amplitude.", &
                 default=.false., do_not_log=.not.CS%do_skeb)
  if ((.not. CS%do_skeb) .and. (CS%skeb_use_gm)) call MOM_error(FATAL, "If SKEB_USE_GM is True "//&
                 "then DO_SKEB must also be True.")
  call get_param(param_file, mdl, "SKEB_GM_COEF", CS%skeb_gm_coef, &
               "Fraction of GM work that is added to backscatter rate.", &
               units="nondim", default=0.0, do_not_log=.not.CS%skeb_use_gm)
  call get_param(param_file, mdl, "SKEB_USE_FRICT", CS%skeb_use_frict, &
                 "If true, adds horizontal friction dissipation rate "//&
                 "to the SKEBS amplitude.", default=.false., do_not_log=.not.CS%do_skeb)
  if ((.not. CS%do_skeb) .and. (CS%skeb_use_frict)) call MOM_error(FATAL, "If SKEB_USE_FRICT is "//&
                 "True then DO_SKEB must also be True.")
  call get_param(param_file, mdl, "SKEB_FRICT_COEF", CS%skeb_frict_coef, &
               "Fraction of horizontal friction work that is added to backscatter rate.", &
               units="nondim", default=0.0, do_not_log=.not.CS%skeb_use_frict)
  call get_param(param_file, mdl, "PERT_EPBL", CS%pert_epbl, &
                 "If true, then stochastically perturb the kinetic energy "//&
                 "production and dissipation terms.  Amplitude and correlations are "//&
                 "controlled by the nam_stoch namelist in the UFS model only.", &
                 default=.false.)

  if (CS%do_sppt .OR. CS%pert_epbl .OR. CS%do_skeb) then
    num_procs = num_PEs()
    allocate(pelist(num_procs))
    call Get_PElist(pelist,commID = mom_comm)
    pe_zero = root_PE()
    nxT = grid%ied - grid%isd + 1
    nyT = grid%jed - grid%jsd + 1
    nxB = grid%iedB - grid%isdB + 1
    nyB = grid%jedB - grid%jsdB + 1
    call init_stochastic_physics_ocn(dt, grid%geoLonT, grid%geoLatT, nxT, nyT, GV%ke, &
                                     grid%geoLonBu, grid%geoLatBu, nxB, nyB, &
                                     CS%pert_epbl, CS%do_sppt, CS%do_skeb, pe_zero, mom_comm, iret)
    if (iret/=0)  then
      call MOM_error(FATAL, "call to init_stochastic_physics_ocn failed")
      return
    endif

    if (CS%do_sppt) allocate(CS%sppt_wts(grid%isd:grid%ied,grid%jsd:grid%jed))
    if (CS%do_skeb) allocate(CS%skeb_wts(grid%isdB:grid%iedB,grid%jsdB:grid%jedB))
    if (CS%do_skeb) allocate(CS%skeb_diss(grid%isd:grid%ied,grid%jsd:grid%jed,GV%ke), source=0.)
    if (CS%pert_epbl) then
      allocate(CS%epbl1_wts(grid%isd:grid%ied,grid%jsd:grid%jed))
      allocate(CS%epbl2_wts(grid%isd:grid%ied,grid%jsd:grid%jed))
    endif
  endif

  CS%id_sppt_wts = register_diag_field('ocean_model', 'sppt_pattern', CS%diag%axesT1, Time, &
       'random pattern for sppt', 'None')
  CS%id_skeb_wts = register_diag_field('ocean_model', 'skeb_pattern', CS%diag%axesB1, Time, &
       'random pattern for skeb', 'None')
  CS%id_epbl1_wts = register_diag_field('ocean_model', 'epbl1_wts', CS%diag%axesT1, Time, &
      'random pattern for KE generation', 'None')
  CS%id_epbl2_wts = register_diag_field('ocean_model', 'epbl2_wts', CS%diag%axesT1, Time, &
      'random pattern for KE dissipation', 'None')
  CS%id_skebu = register_diag_field('ocean_model', 'skebu', CS%diag%axesCuL, Time, &
       'zonal current perts', 'None')
  CS%id_skebv = register_diag_field('ocean_model', 'skebv', CS%diag%axesCvL, Time, &
       'zonal current perts', 'None')
  CS%id_diss = register_diag_field('ocean_model', 'skeb_amp', CS%diag%axesTL, Time, &
       'SKEB amplitude', 'm s-1')
  CS%id_psi  = register_diag_field('ocean_model', 'psi', CS%diag%axesBL, Time, &
       'stream function', 'None')
  CS%id_skeb_taperu = register_static_field('ocean_model', 'skeb_taper_u', CS%diag%axesCu1, &
       'SKEB taper u', 'None', interp_method='none')
  CS%id_skeb_taperv = register_static_field('ocean_model', 'skeb_taper_v', CS%diag%axesCv1, &
       'SKEB taper v', 'None', interp_method='none')

  ! Initialize the "taper" fields. These fields multiply the components of the stochastic
  ! velocity increment in such a way as to smoothly taper them to zero at land boundaries.
  if ((CS%do_skeb) .or. (CS%id_skeb_taperu > 0) .or. (CS%id_skeb_taperv > 0)) then
    allocate(CS%taperCu(grid%IsdB:grid%IedB,grid%jsd:grid%jed))
    allocate(CS%taperCv(grid%isd:grid%ied,grid%JsdB:grid%JedB))
    ! Initialize taper from land mask
    do j=grid%jsd,grid%jed ; do I=grid%isdB,grid%iedB
      CS%taperCu(I,j) = grid%mask2dCu(I,j)
    enddo ; enddo
    do J=grid%jsdB,grid%jedB ; do i=grid%isd,grid%ied
      CS%taperCv(i,J) = grid%mask2dCv(i,J)
    enddo ; enddo
    ! Extend taper land
    do k=1,(taper_width / 2)
      do j=grid%jsc-1,grid%jec+1 ; do I=grid%iscB-1,grid%iecB+1
        tmp(I,j) = minval(CS%taperCu(I-1:I+1,j-1:j+1))
      enddo ; enddo
      do j=grid%jsc,grid%jec ; do I=grid%iscB,grid%iecB
        CS%taperCu(I,j) = minval(tmp(I-1:I+1,j-1:j+1))
      enddo ; enddo
      do J=grid%jscB-1,grid%jecB+1 ; do i=grid%isc-1,grid%iec+1
        tmp(i,J) = minval(CS%taperCv(i-1:i+1,J-1:J+1))
      enddo ; enddo
      do J=grid%jscB,grid%jecB ; do i=grid%isc,grid%iec
        CS%taperCv(i,J) = minval(tmp(i-1:i+1,J-1:J+1))
      enddo ; enddo
      ! Update halo
      call pass_vector(CS%taperCu, CS%taperCv, grid%Domain, SCALAR_PAIR)
    enddo
    ! Smooth tapers. Each call smooths twice.
    do k=1,(taper_width - (taper_width/2))
      call smooth_x9_uv(grid, CS%taperCu, CS%taperCv, zero_land=.true.)
      call pass_vector(CS%taperCu, CS%taperCv, grid%Domain, SCALAR_PAIR)
    enddo
  endif

  !call uvchksum("SKEB taper [uv]", CS%taperCu, CS%taperCv, grid%HI)

  if (CS%id_skeb_taperu > 0) call post_data(CS%id_skeb_taperu, CS%taperCu, CS%diag, .true.)
  if (CS%id_skeb_taperv > 0) call post_data(CS%id_skeb_taperv, CS%taperCv, CS%diag, .true.)

  if (CS%do_sppt .OR. CS%pert_epbl .OR. CS%do_skeb) &
    call MOM_mesg('            === COMPLETED MOM STOCHASTIC INITIALIZATION =====')

  call callTree_leave("stochastic_init(), MOM_stochastics.F90")

end procedure stochastics_init
module procedure update_stochastics
  call callTree_enter("update_stochastics(), MOM_stochastics.F90")

! update stochastic physics patterns before running next time-step
  call run_stochastic_physics_ocn(CS%sppt_wts,CS%skeb_wts,CS%epbl1_wts,CS%epbl2_wts)

  call callTree_leave("update_stochastics(), MOM_stochastics.F90")

end procedure update_stochastics
module procedure apply_skeb
  real, dimension(SZIB_(grid),SZJB_(grid),SZK_(GV)) :: psi         !< Streamfunction for stochastic velocity increments
  real, dimension(SZIB_(grid),SZJ_(grid) ,SZK_(GV)) :: ustar       !< Stochastic u velocity increment [L T-1 ~> m s-1]
  real, dimension(SZI_(grid) ,SZJB_(grid),SZK_(GV)) :: vstar       !< Stochastic v velocity increment [L T-1 ~> m s-1]
  real, dimension(SZI_(grid),SZJ_(grid))            :: diss_tmp    !< Temporary array used in smoothing skeb_diss
  real, dimension(3,3) :: local_weights                            !< 3x3 stencil weights used in smoothing skeb_diss
  real    :: shr,ten,tot,kh
  integer :: i,j,k,iter
  integer, dimension(2) :: EOSdom ! The i-computational domain for the equation of state
  call callTree_enter("apply_skeb(), MOM_stochastics.F90")

  if ((.not. CS%skeb_use_gm) .and. (.not. CS%skeb_use_frict)) then
    ! fill in halos with zeros
    do k=1,GV%ke
      do j=grid%jsd,grid%jed ; do i=grid%isd,grid%ied
        CS%skeb_diss(i,j,k) = 0.0
      enddo ; enddo
    enddo

    !kh needs to be scaled

    kh=1!(120*111)**2
    do k=1,GV%ke
      do j=grid%jsc,grid%jec ; do i=grid%isc,grid%iec
        ! Shear
        shr = (vc(i,J,k)-vc(i-1,J,k))*grid%mask2dCv(i,J)*grid%mask2dCv(i-1,J)*grid%IdxCv(i,J)+&
              (uc(I,j,k)-uc(I,j-1,k))*grid%mask2dCu(I,j)*grid%mask2dCu(I,j-1)*grid%IdyCu(I,j)
        ! Tension
        ten = (vc(i,J,k)-vc(i-1,J,k))*grid%mask2dCv(i,J)*grid%mask2dCv(i-1,J)*grid%IdyCv(i,J)+&
              (uc(I,j,k)-uc(I,j-1,k))*grid%mask2dCu(I,j)*grid%mask2dCu(I,j-1)*grid%IdxCu(I,j)

        tot = sqrt( shr**2 + ten**2 ) * grid%mask2dT(i,j)
        CS%skeb_diss(i,j,k) = tot**3 * kh * grid%areaT(i,j)!!**2
      enddo ; enddo
    enddo
  endif ! Sets CS%skeb_diss without GM or FrictWork

  ! smooth dissipation skeb_npass times
  do iter=1,CS%skeb_npass
    if (mod(iter,2) == 1) call pass_var(CS%skeb_diss, grid%domain)
    do k=1,GV%ke
      do j=grid%jsc-1,grid%jec+1 ; do i=grid%isc-1,grid%iec+1
        ! This does not preserve rotational symmetry
        local_weights = grid%mask2dT(i-1:i+1,j-1:j+1)*grid%areaT(i-1:i+1,j-1:j+1)
        diss_tmp(i,j) = sum(local_weights*CS%skeb_diss(i-1:i+1,j-1:j+1,k)) / &
                       (sum(local_weights) + 1.E-16)
      enddo ; enddo
      do j=grid%jsc-1,grid%jec+1 ; do i=grid%isc-1,grid%iec+1
        if (grid%mask2dT(i,j)==0.) cycle
        CS%skeb_diss(i,j,k) = diss_tmp(i,j)
      enddo ; enddo
    enddo
  enddo
  call pass_var(CS%skeb_diss, grid%domain)

  ! call hchksum(CS%skeb_diss, "SKEB DISS", grid%HI, haloshift=2)
  ! call qchksum(CS%skeb_wts, "SKEB WTS", grid%HI, haloshift=1)

  do k=1,GV%ke
    do J=grid%jscB-1,grid%jecB ; do I=grid%iscB-1,grid%iecB
      psi(I,J,k) = sqrt(0.25 * dt * max((CS%skeb_diss(i  ,j  ,k) + CS%skeb_diss(i+1,j+1,k)) + &
                                        (CS%skeb_diss(i  ,j+1,k) + CS%skeb_diss(i+1,j  ,k)), 0.) ) &
                                  * CS%skeb_wts(I,J)
    enddo ; enddo
  enddo
  !call qchksum(psi,"SKEB PSI", grid%HI, haloshift=1)
  !call pass_var(psi, grid%domain, position=CORNER)
  do k=1,GV%ke
    do j=grid%jsc,grid%jec ; do I=grid%iscB,grid%iecB
      ustar(I,j,k) = - (psi(I,J,k) - psi(I,J-1,k)) * CS%taperCu(I,j) * grid%IdyCu(I,j)
      uc(I,j,k) = uc(I,j,k) + ustar(I,j,k)
    enddo ; enddo
    do J=grid%jscB,grid%jecB ; do i=grid%isc,grid%iec
      vstar(i,J,k) =   (psi(I,J,k) - psi(I-1,J,k)) * CS%taperCv(i,J) * grid%IdxCv(i,J)
      vc(i,J,k) = vc(i,J,k) + vstar(i,J,k)
    enddo ; enddo
  enddo

  !call uvchksum("SKEB increment [uv]", ustar, vstar, grid%HI)

  call enable_averages(dt, Time_end, CS%diag)
  if (CS%id_diss > 0) then
     call post_data(CS%id_diss, sqrt(dt * max(CS%skeb_diss(:,:,:), 0.)), CS%diag)
  endif
  if (CS%id_skeb_wts > 0) then
     call post_data(CS%id_skeb_wts, CS%skeb_wts, CS%diag)
  endif
  if (CS%id_skebu > 0) then
     call post_data(CS%id_skebu, ustar(:,:,:), CS%diag)
  endif
  if (CS%id_skebv > 0) then
     call post_data(CS%id_skebv, vstar(:,:,:), CS%diag)
  endif
  if (CS%id_psi > 0) then
     call post_data(CS%id_psi, psi(:,:,:), CS%diag)
  endif
  call disable_averaging(CS%diag)
  CS%skeb_diss(:,:,:) = 0.0 ! Must zero before next time step.

  call callTree_leave("apply_skeb(), MOM_stochastics.F90")

end procedure apply_skeb
module procedure smooth_x9_uv
  real :: fu_prev(SZIB_(G),SZJ_(G))  ! The value of the u-point field at the previous iteration [arbitrary]
  real :: fv_prev(SZI_(G),SZJB_(G))  ! The value of the v-point field at the previous iteration [arbitrary]
  real :: Iwts             ! The inverse of the sum of the weights [nondim]
  logical :: zero_land_val ! The value of the zero_land optional argument or .true. if it is absent.
  integer :: i, j, s, is, ie, js, je, Isq, Ieq, Jsq, Jeq
  is  = G%isc  ; ie  = G%iec  ; js  = G%jsc  ; je  = G%jec
  Isq = G%IscB ; Ieq = G%IecB ; Jsq = G%JscB ; Jeq = G%JecB

  zero_land_val = .true. ; if (present(zero_land)) zero_land_val = zero_land

  do s=1,0,-1
    fu_prev(:,:) = field_u(:,:)
    ! apply smoothing on field_u using rotationally symmetric expressions.
    do j=js-s,je+s ; do I=Isq-s,Ieq+s ; if (G%mask2dCu(I,j) > 0.0) then
      Iwts = 0.0625
      if (.not. zero_land_val) &
        Iwts = 1.0 / ( (4.0*G%mask2dCu(I,j) + &
                        ( 2.0*((G%mask2dCu(I-1,j) + G%mask2dCu(I+1,j)) + &
                               (G%mask2dCu(I,j-1) + G%mask2dCu(I,j+1))) + &
                         ((G%mask2dCu(I-1,j-1) + G%mask2dCu(I+1,j+1)) + &
                          (G%mask2dCu(I-1,j+1) + G%mask2dCu(I+1,j-1))) ) ) + 1.0e-16 )
      field_u(I,j) = Iwts * ( 4.0*G%mask2dCu(I,j) * fu_prev(I,j) &
                            + (2.0*((G%mask2dCu(I-1,j) * fu_prev(I-1,j) + G%mask2dCu(I+1,j) * fu_prev(I+1,j)) + &
                                    (G%mask2dCu(I,j-1) * fu_prev(I,j-1) + G%mask2dCu(I,j+1) * fu_prev(I,j+1))) &
                              + ((G%mask2dCu(I-1,j-1) * fu_prev(I-1,j-1) + G%mask2dCu(I+1,j+1) * fu_prev(I+1,j+1)) + &
                                 (G%mask2dCu(I-1,j+1) * fu_prev(I-1,j+1) + G%mask2dCu(I+1,j-1) * fu_prev(I-1,j-1))) ))
    endif ; enddo ; enddo

    fv_prev(:,:) = field_v(:,:)
    ! apply smoothing on field_v using rotationally symmetric expressions.
    do J=Jsq-s,Jeq+s ; do i=is-s,ie+s ; if (G%mask2dCv(i,J) > 0.0) then
      Iwts = 0.0625
      if (.not. zero_land_val) &
        Iwts = 1.0 / ( (4.0*G%mask2dCv(i,J) + &
                        ( 2.0*((G%mask2dCv(i-1,J) + G%mask2dCv(i+1,J)) + &
                               (G%mask2dCv(i,J-1) + G%mask2dCv(i,J+1))) + &
                         ((G%mask2dCv(i-1,J-1) + G%mask2dCv(i+1,J+1)) + &
                          (G%mask2dCv(i-1,J+1) + G%mask2dCv(i+1,J-1))) ) ) + 1.0e-16 )
      field_v(i,J) = Iwts * ( 4.0*G%mask2dCv(i,J) * fv_prev(i,J) &
                            + (2.0*((G%mask2dCv(i-1,J) * fv_prev(i-1,J) + G%mask2dCv(i+1,J) * fv_prev(i+1,J)) + &
                                    (G%mask2dCv(i,J-1) * fv_prev(i,J-1) + G%mask2dCv(i,J+1) * fv_prev(i,J+1))) &
                              + ((G%mask2dCv(i-1,J-1) * fv_prev(i-1,J-1) + G%mask2dCv(i+1,J+1) * fv_prev(i+1,J+1)) + &
                                 (G%mask2dCv(i-1,J+1) * fv_prev(i-1,J+1) + G%mask2dCv(i+1,J-1) * fv_prev(i-1,J-1))) ))
    endif ; enddo ; enddo
  enddo

end procedure smooth_x9_uv
end submodule MOM_stochastics_s
