submodule (MOM_Zanna_Bolton) MOM_Zanna_Bolton_s
#include <MOM_memory.h>
  implicit none
contains
module procedure ZB2020_init
  real :: subroundoff_Cor     ! A negligible parameter which avoids division by zero
  integer :: is, ie, js, je, Isq, Ieq, Jsq, Jeq
  integer :: i, j
#include "version_variable.h"
  character(len=40)  :: mdl = "MOM_Zanna_Bolton" ! This module's name.
  is  = G%isc  ; ie  = G%iec  ; js  = G%jsc  ; je  = G%jec
  Isq = G%IscB ; Ieq = G%IecB ; Jsq = G%JscB ; Jeq = G%JecB

  call log_version(param_file, mdl, version, "")

  call get_param(param_file, mdl, "USE_ZB2020", use_ZB2020, &
                 "If true, turns on Zanna-Bolton-2020 (ZB) " //&
                 "subgrid momentum parameterization of mesoscale eddies.", default=.false.)
  if (.not. use_ZB2020) return

  call get_param(param_file, mdl, "ZB2020_USE_ANN", CS%use_ann, &
                 "ANN inference of momentum fluxes", default=.false.)

  call get_param(param_file, mdl, "ZB2020_ANN_STENCIL_SIZE", CS%stencil_size, &
                 "ANN stencil size", default=3)

  call get_param(param_file, mdl, "ZB2020_ANN_FILE_TALL", CS%ann_file_Tall, &
                 "ANN parameters for prediction of Txy, Txx and Tyy netcdf input", &
                 default="INPUT/EXP1/Tall.nc")

  call get_param(param_file, mdl, "ZB_SCALING", CS%amplitude, &
                 "The nondimensional scaling factor in ZB model, " //&
                 "typically 0.5-2.5", units="nondim", default=0.5)

  call get_param(param_file, mdl, "ZB_TRACE_MODE", CS%ZB_type, &
                 "Select how to compute the trace part of ZB model:\n" //&
                 "\t 0 - both deviatoric and trace components are computed\n" //&
                 "\t 1 - only deviatoric component is computed\n" //&
                 "\t 2 - only trace component is computed", default=0)

  call get_param(param_file, mdl, "ZB_SCHEME", CS%ZB_cons, &
                 "Select a discretization scheme for ZB model:\n" //&
                 "\t 0 - non-conservative scheme\n" //&
                 "\t 1 - conservative scheme for deviatoric component", default=1)

  call get_param(param_file, mdl, "VG_SHARP_PASS", CS%HPF_iter, &
                "Number of sharpening passes for the Velocity Gradient (VG) components " //&
                "in ZB model.", default=0)

  call get_param(param_file, mdl, "STRESS_SMOOTH_PASS", CS%Stress_iter, &
                 "Number of smoothing passes for the Stress tensor components " //&
                 "in ZB model.", default=0)

  call get_param(param_file, mdl, "ZB_KLOWER_R_DISS", CS%Klower_R_diss, &
                 "Attenuation of " //&
                 "the ZB parameterization in the regions of " //&
                 "geostrophically-unbalanced flows (Klower 2018, Juricke2020,2019). " //&
                 "Subgrid stress is multiplied by 1/(1+(shear/(f*R_diss))):\n" //&
                 "\t R_diss=-1. - attenuation is not used\n\t R_diss= 1. - typical value", &
                 units="nondim", default=-1.)

  call get_param(param_file, mdl, "ZB_KLOWER_SHEAR", CS%Klower_shear, &
                 "Type of expression for shear in Klower formula:\n" //&
                 "\t 0: sqrt(sh_xx**2 + sh_xy**2)\n" //&
                 "\t 1: sqrt(sh_xx**2 + sh_xy**2 + vort_xy**2)", &
                 default=1, do_not_log=.not.CS%Klower_R_diss>0)

  call get_param(param_file, mdl, "ZB_MARCHING_HALO", CS%Marching_halo, &
                 "The number of filter iterations per single MPI " //&
                 "exchange", default=4, do_not_log=(CS%Stress_iter==0).and.(CS%HPF_iter==0))

  ! Register fields for output from this module.
  CS%diag => diag

  CS%id_ZB2020u = register_diag_field('ocean_model', 'ZB2020u', diag%axesCuL, Time, &
      'Zonal Acceleration from Zanna-Bolton 2020', 'm s-2', conversion=US%L_T2_to_m_s2)
  CS%id_ZB2020v = register_diag_field('ocean_model', 'ZB2020v', diag%axesCvL, Time, &
      'Meridional Acceleration from Zanna-Bolton 2020', 'm s-2', conversion=US%L_T2_to_m_s2)
  CS%id_KE_ZB2020 = register_diag_field('ocean_model', 'KE_ZB2020', diag%axesTL, Time, &
      'Kinetic Energy Source from Horizontal Viscosity', &
      'm3 s-3', conversion=GV%H_to_m*(US%L_T_to_m_s**2)*US%s_to_T)

  CS%id_Txx = register_diag_field('ocean_model', 'Txx', diag%axesTL, Time, &
      'Diagonal term (Txx) in the ZB stress tensor', 'm2 s-2', conversion=US%L_T_to_m_s**2)

  CS%id_Tyy = register_diag_field('ocean_model', 'Tyy', diag%axesTL, Time, &
      'Diagonal term (Tyy) in the ZB stress tensor', 'm2 s-2', conversion=US%L_T_to_m_s**2)

  CS%id_Txy = register_diag_field('ocean_model', 'Txy', diag%axesBL, Time, &
      'Off-diagonal term (Txy) in the ZB stress tensor', 'm2 s-2', conversion=US%L_T_to_m_s**2)

  if (CS%Klower_R_diss > 0) then
    CS%id_cdiss = register_diag_field('ocean_model', 'c_diss', diag%axesTL, Time, &
        'Klower (2018) attenuation coefficient', 'nondim')
  endif

  ! Clock IDs
  ! Only module is measured with syncronization. While smaller
  ! parts are measured without - because these are nested clocks.
  CS%id_clock_module = cpu_clock_id('(Ocean Zanna-Bolton-2020)', grain=CLOCK_MODULE)
  CS%id_clock_copy = cpu_clock_id('(ZB2020 copy fields)', grain=CLOCK_ROUTINE, sync=.false.)
  CS%id_clock_cdiss = cpu_clock_id('(ZB2020 compute c_diss)', grain=CLOCK_ROUTINE, sync=.false.)
  CS%id_clock_stress = cpu_clock_id('(ZB2020 compute stress)', grain=CLOCK_ROUTINE, sync=.false.)
  CS%id_clock_stress_ANN = cpu_clock_id('(ZB2020 compute stress ANN)', grain=CLOCK_ROUTINE, sync=.false.)
  CS%id_clock_divergence = cpu_clock_id('(ZB2020 compute divergence)', grain=CLOCK_ROUTINE, sync=.false.)
  CS%id_clock_mpi = cpu_clock_id('(ZB2020 filter MPI exchanges)', grain=CLOCK_ROUTINE, sync=.false.)
  CS%id_clock_filter = cpu_clock_id('(ZB2020 filter no MPI)', grain=CLOCK_ROUTINE, sync=.false.)
  CS%id_clock_post = cpu_clock_id('(ZB2020 post data)', grain=CLOCK_ROUTINE, sync=.false.)
  CS%id_clock_source = cpu_clock_id('(ZB2020 compute energy source)', grain=CLOCK_ROUTINE, sync=.false.)

  CS%subroundoff_shear = 1e-30 * US%T_to_s
  if (CS%use_ann) then
    call ANN_init(CS%ann_Tall, CS%ann_file_Tall)
  endif

  ! Allocate memory
  ! We set the stress tensor and velocity gradient tensor to zero
  ! with full halo because they potentially may be filtered
  ! with marching halo algorithm
  allocate(CS%sh_xx(SZI_(G),SZJ_(G),SZK_(GV)), source=0.)
  allocate(CS%sh_xy(SZIB_(G),SZJB_(G),SZK_(GV)), source=0.)
  allocate(CS%vort_xy(SZIB_(G),SZJB_(G),SZK_(GV)), source=0.)
  allocate(CS%hq(SZIB_(G),SZJB_(G),SZK_(GV)))

  allocate(CS%Txx(SZI_(G),SZJ_(G),SZK_(GV)), source=0.)
  allocate(CS%Tyy(SZI_(G),SZJ_(G),SZK_(GV)), source=0.)
  allocate(CS%Txy(SZIB_(G),SZJB_(G),SZK_(GV)), source=0.)
  allocate(CS%kappa_h(SZI_(G),SZJ_(G)))
  allocate(CS%kappa_q(SZIB_(G),SZJB_(G)))

  ! Precomputing the scaling coefficient
  ! Mask is included to automatically satisfy B.C.
  do j=js-2,je+2 ; do i=is-2,ie+2
    CS%kappa_h(i,j) = -CS%amplitude * G%areaT(i,j) * G%mask2dT(i,j)
  enddo ; enddo

  do J=Jsq-2,Jeq+2 ; do I=Isq-2,Ieq+2
    CS%kappa_q(I,J) = -CS%amplitude * G%areaBu(I,J) * G%mask2dBu(I,J)
  enddo ; enddo

  if (CS%Klower_R_diss > 0) then
    allocate(CS%ICoriolis_h(SZI_(G),SZJ_(G)))
    allocate(CS%c_diss(SZI_(G),SZJ_(G),SZK_(GV)))

    subroundoff_Cor = 1e-30 * US%T_to_s
    ! Precomputing 1/(f * R_diss)
    do j=js-1,je+1 ; do i=is-1,ie+1
      CS%ICoriolis_h(i,j) = 1. / ((abs(0.25 * ((G%CoriolisBu(I,J) + G%CoriolisBu(I-1,J-1)) &
                          + (G%CoriolisBu(I-1,J) + G%CoriolisBu(I,J-1)))) + subroundoff_Cor) &
                          * CS%Klower_R_diss)
    enddo ; enddo
  endif

  if (CS%Stress_iter > 0 .or. CS%HPF_iter > 0) then
    ! Include 1/16. factor to the mask for filter implementation
    allocate(CS%maskw_h(SZI_(G),SZJ_(G))) ; CS%maskw_h(:,:) = G%mask2dT(:,:) * 0.0625
    allocate(CS%maskw_q(SZIB_(G),SZJB_(G))) ; CS%maskw_q(:,:) = G%mask2dBu(:,:) * 0.0625
  endif

  ! Initialize MPI group passes
  if (CS%Stress_iter > 0) then
    ! reduce size of halo exchange accordingly to
    ! Marching halo, number of iterations and the array size
    ! But let exchange width be at least 1
    CS%Stress_halo = max(min(CS%Marching_halo, CS%Stress_iter, &
                             G%Domain%nihalo, G%Domain%njhalo), 1)

    call create_group_pass(CS%pass_Tq, CS%Txy, G%Domain, halo=CS%Stress_halo, &
      position=CORNER)
    call create_group_pass(CS%pass_Th, CS%Txx, G%Domain, halo=CS%Stress_halo)
    call create_group_pass(CS%pass_Th, CS%Tyy, G%Domain, halo=CS%Stress_halo)
  endif

  if (CS%HPF_iter > 0) then
    ! The minimum halo size is 2 because it is requirement for the
    ! outputs of function filter_velocity_gradients
    CS%HPF_halo = max(min(CS%Marching_halo, CS%HPF_iter, &
                          G%Domain%nihalo, G%Domain%njhalo), 2)

    call create_group_pass(CS%pass_xx, CS%sh_xx, G%Domain, halo=CS%HPF_halo)
    call create_group_pass(CS%pass_xy, CS%sh_xy, G%Domain, halo=CS%HPF_halo, &
      position=CORNER)
    call create_group_pass(CS%pass_xy, CS%vort_xy, G%Domain, halo=CS%HPF_halo, &
      position=CORNER)
  endif

end procedure ZB2020_init
module procedure ZB2020_end
  deallocate(CS%sh_xx)
  deallocate(CS%sh_xy)
  deallocate(CS%vort_xy)
  deallocate(CS%hq)

  deallocate(CS%Txx)
  deallocate(CS%Tyy)
  deallocate(CS%Txy)
  deallocate(CS%kappa_h)
  deallocate(CS%kappa_q)

  if (CS%Klower_R_diss > 0) then
    deallocate(CS%ICoriolis_h)
    deallocate(CS%c_diss)
  endif

  if (CS%Stress_iter > 0 .or. CS%HPF_iter > 0) then
    deallocate(CS%maskw_h)
    deallocate(CS%maskw_q)
  endif

  if (CS%use_ann) then
    call ANN_end(CS%ann_Tall)
  endif

end procedure ZB2020_end
module procedure ZB2020_copy_gradient_and_thickness
  integer :: is, ie, js, je, Isq, Ieq, Jsq, Jeq
  integer :: i, j
  call cpu_clock_begin(CS%id_clock_copy)

  is  = G%isc  ; ie  = G%iec  ; js  = G%jsc  ; je  = G%jec
  Isq = G%IscB ; Ieq = G%IecB ; Jsq = G%JscB ; Jeq = G%JecB

  do J=js-1,Jeq ; do I=is-1,Ieq
    CS%hq(I,J,k) = hq(I,J)
  enddo ; enddo

  ! No physical B.C. is required for
  ! sh_xx in ZB2020. However, filtering
  ! may require BC
  do j=Jsq-1,je+2 ; do i=Isq-1,ie+2
    CS%sh_xx(i,j,k) = sh_xx(i,j) * G%mask2dT(i,j)
  enddo ; enddo

  ! We multiply by mask to remove
  ! implicit dependence on CS%no_slip
  ! flag in hor_visc module
  do J=js-2,Jeq+1 ; do I=is-2,Ieq+1
    CS%sh_xy(I,J,k) = sh_xy(I,J) * G%mask2dBu(I,J)
  enddo ; enddo

  do J=js-2,Jeq+1 ; do I=is-2,Ieq+1
    CS%vort_xy(I,J,k) = vort_xy(I,J) * G%mask2dBu(I,J)
  enddo ; enddo

  call cpu_clock_end(CS%id_clock_copy)

end procedure ZB2020_copy_gradient_and_thickness
module procedure ZB2020_lateral_stress
  call cpu_clock_begin(CS%id_clock_module)

  ! Compute attenuation if specified
  call compute_c_diss(G, GV, CS)

  ! Sharpen velocity gradients if specified
  call filter_velocity_gradients(G, GV, CS)

  ! Compute the stress tensor given the
  ! (optionally sharpened) velocity gradients
  if (CS%use_ann) then
    call compute_stress_ANN_collocated(G, GV, CS)
  else
    call compute_stress(G, GV, CS)
  endif

  ! Smooth the stress tensor if specified
  call filter_stress(G, GV, CS)

  ! Update the acceleration due to eddy viscosity (diffu, diffv)
  ! with the ZB2020 lateral parameterization
  call compute_stress_divergence(u, v, h, diffu, diffv,    &
                                 dx2h, dy2h, dx2q, dy2q, &
                                 G, GV, CS)

  call cpu_clock_begin(CS%id_clock_post)
  if (CS%id_Txx>0)       call post_data(CS%id_Txx, CS%Txx, CS%diag)
  if (CS%id_Tyy>0)       call post_data(CS%id_Tyy, CS%Tyy, CS%diag)
  if (CS%id_Txy>0)       call post_data(CS%id_Txy, CS%Txy, CS%diag)

  if (CS%id_cdiss>0)     call post_data(CS%id_cdiss, CS%c_diss, CS%diag)
  call cpu_clock_end(CS%id_clock_post)

  call cpu_clock_end(CS%id_clock_module)

end procedure ZB2020_lateral_stress
module procedure compute_c_diss
  integer :: is, ie, js, je, Isq, Ieq, Jsq, Jeq, nz
  integer :: i, j, k
  real :: shear ! Shear in Klower2018 formula at h points [T-1 ~> s-1]
  if (.not. CS%Klower_R_diss > 0) &
    return

  call cpu_clock_begin(CS%id_clock_cdiss)

  is  = G%isc  ; ie  = G%iec  ; js  = G%jsc  ; je  = G%jec ; nz = GV%ke
  Isq = G%IscB ; Ieq = G%IecB ; Jsq = G%JscB ; Jeq = G%JecB

  do k=1,nz

    ! sqrt(sh_xx**2 + sh_xy**2)
    if (CS%Klower_shear == 0) then
      do j=js-1,je+1 ; do i=is-1,ie+1
        shear = sqrt(CS%sh_xx(i,j,k)**2 + 0.25 * (          &
                     ((CS%sh_xy(I-1,J-1,k)**2) + (CS%sh_xy(I,J  ,k)**2)) &
                   + ((CS%sh_xy(I-1,J  ,k)**2) + (CS%sh_xy(I,J-1,k)**2)) &
                    ))
        CS%c_diss(i,j,k) = 1. / (1. + shear * CS%ICoriolis_h(i,j))
      enddo ; enddo

    ! sqrt(sh_xx**2 + sh_xy**2 + vort_xy**2)
    elseif (CS%Klower_shear == 1) then
      do j=js-1,je+1 ; do i=is-1,ie+1
        shear = sqrt(CS%sh_xx(i,j,k)**2 + 0.25 * (             &
                     ((CS%sh_xy(I-1,J-1,k)**2 + CS%vort_xy(I-1,J-1,k)**2) &
                   +  (CS%sh_xy(I,J,k)**2     + CS%vort_xy(I,J,k)**2))    &
                   + ((CS%sh_xy(I-1,J,k)**2   + CS%vort_xy(I-1,J,k)**2)   &
                   +  (CS%sh_xy(I,J-1,k)**2   + CS%vort_xy(I,J-1,k)**2))  &
                    ))
        CS%c_diss(i,j,k) = 1. / (1. + shear * CS%ICoriolis_h(i,j))
      enddo ; enddo
    endif

  enddo ! end of k loop

  call cpu_clock_end(CS%id_clock_cdiss)

end procedure compute_c_diss
module procedure compute_stress
  real :: &
    vort_xy_h, &  ! Vorticity interpolated to h point [T-1 ~> s-1]
    sh_xy_h       ! Shearing strain interpolated to h point [T-1 ~> s-1]
  real :: &
    sh_xx_q       ! Horizontal tension interpolated to q point [T-1 ~> s-1]
  real :: sum_sq  ! 1/2*(vort_xy^2 + sh_xy^2 + sh_xx^2) in h point [T-2 ~> s-2]
  real :: vort_sh ! vort_xy*sh_xy in h point [T-2 ~> s-2]
  integer :: is, ie, js, je, Isq, Ieq, Jsq, Jeq, nz
  integer :: i, j, k
  logical :: sum_sq_flag ! Flag to compute trace
  logical :: vort_sh_scheme_0, vort_sh_scheme_1 ! Flags to compute diagonal trace-free part
  call cpu_clock_begin(CS%id_clock_stress)

  is  = G%isc  ; ie  = G%iec  ; js  = G%jsc  ; je  = G%jec ; nz = GV%ke
  Isq = G%IscB ; Ieq = G%IecB ; Jsq = G%JscB ; Jeq = G%JecB

  sum_sq = 0.
  vort_sh = 0.

  sum_sq_flag = CS%ZB_type /= 1
  vort_sh_scheme_0 = CS%ZB_type /= 2 .and. CS%ZB_cons == 0
  vort_sh_scheme_1 = CS%ZB_type /= 2 .and. CS%ZB_cons == 1

  do k=1,nz

    ! compute Txx, Tyy tensor
    do j=js-1,je+1 ; do i=is-1,ie+1
      ! It is assumed that B.C. is applied to sh_xy and vort_xy
      sh_xy_h = 0.25 * ( (CS%sh_xy(I-1,J-1,k) + CS%sh_xy(I,J,k)) &
                       + (CS%sh_xy(I-1,J,k) + CS%sh_xy(I,J-1,k)) )

      vort_xy_h = 0.25 * ( (CS%vort_xy(I-1,J-1,k) + CS%vort_xy(I,J,k)) &
                         + (CS%vort_xy(I-1,J,k) + CS%vort_xy(I,J-1,k)) )

      if (sum_sq_flag) then
        sum_sq = 0.5 *                          &
          ((vort_xy_h * vort_xy_h               &
           + sh_xy_h * sh_xy_h)                 &
           + CS%sh_xx(i,j,k) * CS%sh_xx(i,j,k)  &
            )
      endif

      if (vort_sh_scheme_0) &
        vort_sh = vort_xy_h * sh_xy_h

      if (vort_sh_scheme_1) then
        ! It is assumed that B.C. is applied to sh_xy and vort_xy
        vort_sh = 0.25 * (                                                      &
          (((G%areaBu(I-1,J-1) * CS%vort_xy(I-1,J-1,k)) * CS%sh_xy(I-1,J-1,k))  + &
           ((G%areaBu(I  ,J  ) * CS%vort_xy(I  ,J  ,k)) * CS%sh_xy(I  ,J  ,k))) + &
          (((G%areaBu(I-1,J  ) * CS%vort_xy(I-1,J  ,k)) * CS%sh_xy(I-1,J  ,k))  + &
           ((G%areaBu(I  ,J-1) * CS%vort_xy(I  ,J-1,k)) * CS%sh_xy(I  ,J-1,k)))   &
          ) * G%IareaT(i,j)
      endif

      ! B.C. is already applied in kappa_h
      CS%Txx(i,j,k) = CS%kappa_h(i,j) * (- vort_sh + sum_sq)
      CS%Tyy(i,j,k) = CS%kappa_h(i,j) * (+ vort_sh + sum_sq)

    enddo ; enddo

    ! Here we assume that Txy is initialized to zero
    if (CS%ZB_type /= 2) then
      do J=Jsq-1,Jeq+1 ; do I=Isq-1,Ieq+1
        sh_xx_q = 0.25 * ( (CS%sh_xx(i+1,j+1,k) + CS%sh_xx(i,j,k)) &
                         + (CS%sh_xx(i+1,j,k) + CS%sh_xx(i,j+1,k)))
        ! B.C. is already applied in kappa_q
        CS%Txy(I,J,k) = CS%kappa_q(I,J) * (CS%vort_xy(I,J,k) * sh_xx_q)

      enddo ; enddo
    endif

  enddo ! end of k loop

  call cpu_clock_end(CS%id_clock_stress)

end procedure compute_stress
module procedure compute_stress_ANN_collocated
  integer :: is, ie, js, je, Isq, Ieq, Jsq, Jeq, nz
  integer :: i, j, k, m
  integer :: ii, jj
  integer :: nij
  real, allocatable :: x(:,:)        ! Vector of non-dimensional input features
  real, allocatable :: y(:,:)        ! Vector of nondimensional
  real :: yy(3)                      ! Vector of dimensional
  real :: tmp                        ! Temporal value of squared norm [T-2 ~> s-2]
  integer :: offset                  ! Half the stencil size. Used for selection
  integer :: stencil_points          ! The number of points after flattening
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)) :: &
        sh_xy_h,   & ! sh_xy interpolated to the center [T-1 ~> s-1]
        vort_xy_h, & ! vort_xy interpolated to the center [T-1 ~> s-1]
        norm_h       ! Norm of input feautres in center points [T-1 ~> s-1]
  real, dimension(SZI_(G),SZJ_(G)) :: &
        sqr_h, & ! Squared norm of velocity gradients in center points [T-2 ~> s-2]
        Txy      ! Predicted Txy in center points                      [T-1 ~> s-1]
  call cpu_clock_begin(CS%id_clock_stress_ANN)

  is  = G%isc  ; ie  = G%iec  ; js  = G%jsc  ; je  = G%jec ; nz = GV%ke
  Isq = G%IscB ; Ieq = G%IecB ; Jsq = G%JscB ; Jeq = G%JecB

  ! Number of horizontal grid points in ANN inference loop below
  nij = (ie - is + 5) * (je - js + 5)
  allocate(x(nij, 3 * CS%stencil_size**2))
  allocate(y(nij, 3))

  sh_xy_h = 0.
  vort_xy_h = 0.
  norm_h = 0.

  call pass_var(CS%sh_xy, G%Domain, clock=CS%id_clock_mpi, position=CORNER)
  call pass_var(CS%sh_xx, G%Domain, clock=CS%id_clock_mpi)
  call pass_var(CS%vort_xy, G%Domain, clock=CS%id_clock_mpi, position=CORNER)

  offset = (CS%stencil_size-1)/2
  stencil_points = CS%stencil_size**2

  ! Interpolate input features
  do k=1,nz
    do j=js-2,je+2 ; do i=is-2,ie+2
      ! It is assumed that B.C. is applied to sh_xy and vort_xy
      sh_xy_h(i,j,k) = 0.25 * ( (CS%sh_xy(I-1,J-1,k) + CS%sh_xy(I,J,k)) &
                              + (CS%sh_xy(I-1,J,k) + CS%sh_xy(I,J-1,k)) )

      vort_xy_h(i,j,k) = 0.25 * ( (CS%vort_xy(I-1,J-1,k) + CS%vort_xy(I,J,k)) &
                                + (CS%vort_xy(I-1,J,k) + CS%vort_xy(I,J-1,k)) )

      sqr_h(i,j) = (((CS%sh_xx(i,j,k)**2) + (sh_xy_h(i,j,k)**2)) + (vort_xy_h(i,j,k)**2)) * G%mask2dT(i,j)
    enddo ; enddo

    do j=js,je ; do i=is,ie
      tmp = 0.0
      do jj=j-offset,j+offset ; do ii=i-offset,i+offset
        tmp = tmp + sqr_h(ii,jj)
      enddo ; enddo
      norm_h(i,j,k) = sqrt(tmp)
    enddo ; enddo
  enddo

  call pass_var(sh_xy_h, G%Domain, clock=CS%id_clock_mpi)
  call pass_var(vort_xy_h, G%Domain, clock=CS%id_clock_mpi)
  call pass_var(norm_h, G%Domain, clock=CS%id_clock_mpi)

  do k=1,nz
    m = 0
    do j=js-2,je+2 ; do i=is-2,ie+2
      m = m + 1
      x(m,1:stencil_points) =                                                            &
                        RESHAPE(sh_xy_h(i-offset:i+offset,                               &
                                        j-offset:j+offset,k), (/stencil_points/))
      x(m,stencil_points+1:2*stencil_points) =                                           &
                        RESHAPE(CS%sh_xx(i-offset:i+offset,                              &
                                         j-offset:j+offset,k), (/stencil_points/))
      x(m,2*stencil_points+1:3*stencil_points) =                                         &
                        RESHAPE(vort_xy_h(i-offset:i+offset,                             &
                                          j-offset:j+offset,k), (/stencil_points/))

      x(m,:) = x(m,:) / (norm_h(i,j,k) + CS%subroundoff_shear)
    enddo ; enddo

    call ANN_apply_array_sio(nij, x, y, CS%ann_Tall)

    m = 0
    do j=js-2,je+2 ; do i=is-2,ie+2
      m = m+1
      yy(:) = y(m, :) * norm_h(i,j,k) * norm_h(i,j,k) * CS%kappa_h(i,j)

      Txy(i,j)      = yy(1)
      CS%Txx(i,j,k) = yy(2)
      CS%Tyy(i,j,k) = yy(3)
    enddo ; enddo

    do J=Jsq-1,Jeq+1 ; do I=Isq-1,Ieq+1
      CS%Txy(I,J,k) = 0.25 * ( (Txy(i+1,j+1) + Txy(i,j)) &
                             + (Txy(i+1,j)   + Txy(i,j+1))) * G%mask2dBu(I,J)
    enddo ; enddo

  enddo ! end of k loop

  call pass_var(CS%Txy, G%Domain, clock=CS%id_clock_mpi, position=CORNER)
  call pass_var(CS%Txx, G%Domain, clock=CS%id_clock_mpi)
  call pass_var(CS%Tyy, G%Domain, clock=CS%id_clock_mpi)

  deallocate(x)
  deallocate(y)

  call cpu_clock_end(CS%id_clock_stress_ANN)

end procedure compute_stress_ANN_collocated
module procedure compute_stress_divergence
  real, dimension(SZI_(G),SZJ_(G)) :: &
        Mxx, & ! Subgrid stress Txx multiplied by thickness and dy^2 [H L4 T-2 ~> m5 s-2]
        Myy    ! Subgrid stress Tyy multiplied by thickness and dx^2 [H L4 T-2 ~> m5 s-2]
  real, dimension(SZIB_(G),SZJB_(G)) :: &
        Mxy    ! Subgrid stress Txy multiplied by thickness [H L2 T-2 ~> m3 s-2]
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)) :: &
        ZB2020u           !< Zonal acceleration due to convergence of
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)) :: &
        ZB2020v           !< Meridional acceleration due to convergence
  real :: h_u ! Thickness interpolated to u points [H ~> m or kg m-2].
  real :: h_v ! Thickness interpolated to v points [H ~> m or kg m-2].
  real :: fx  ! Zonal acceleration      [L T-2 ~> m s-2]
  real :: fy  ! Meridional acceleration [L T-2 ~> m s-2]
  real :: h_neglect    ! Thickness so small it can be lost in
  integer :: is, ie, js, je, Isq, Ieq, Jsq, Jeq, nz
  integer :: i, j, k
  logical :: save_ZB2020u, save_ZB2020v ! Save the acceleration due to ZB2020 model
  call cpu_clock_begin(CS%id_clock_divergence)

  save_ZB2020u = (CS%id_ZB2020u > 0) .or. (CS%id_KE_ZB2020 > 0)
  save_ZB2020v = (CS%id_ZB2020v > 0) .or. (CS%id_KE_ZB2020 > 0)

  is  = G%isc  ; ie  = G%iec  ; js  = G%jsc  ; je  = G%jec ; nz = GV%ke
  Isq = G%IscB ; Ieq = G%IecB ; Jsq = G%JscB ; Jeq = G%JecB

  h_neglect  = GV%H_subroundoff

  do k=1,nz
    if (CS%Klower_R_diss > 0) then
      do J=js-1,Jeq ; do I=is-1,Ieq
          Mxy(I,J) = (CS%Txy(I,J,k) *                                         &
                      (0.25 * ( (CS%c_diss(i,j  ,k) + CS%c_diss(i+1,j+1,k))   &
                              + (CS%c_diss(i,j+1,k) + CS%c_diss(i+1,j  ,k)))  &
                      )                                                       &
                     ) * CS%hq(I,J,k)
      enddo ; enddo
    else
      do J=js-1,Jeq ; do I=is-1,Ieq
        Mxy(I,J) = CS%Txy(I,J,k) * CS%hq(I,J,k)
      enddo ; enddo
    endif

    if (CS%Klower_R_diss > 0) then
      do j=js-1,je+1 ; do i=is-1,ie+1
        Mxx(i,j) = ((CS%Txx(i,j,k) * CS%c_diss(i,j,k)) * h(i,j,k)) * dy2h(i,j)
        Myy(i,j) = ((CS%Tyy(i,j,k) * CS%c_diss(i,j,k)) * h(i,j,k)) * dx2h(i,j)
      enddo ; enddo
    else
      do j=js-1,je+1 ; do i=is-1,ie+1
        Mxx(i,j) = ((CS%Txx(i,j,k)) * h(i,j,k)) * dy2h(i,j)
        Myy(i,j) = ((CS%Tyy(i,j,k)) * h(i,j,k)) * dx2h(i,j)
      enddo ; enddo
    endif

    ! Evaluate du/dt=1/h x.Div(h T) (Line 1495 of MOM_hor_visc.F90)
    do j=js,je ; do I=Isq,Ieq
      h_u = 0.5 * (G%mask2dT(i,j)*h(i,j,k) + G%mask2dT(i+1,j)*h(i+1,j,k)) + h_neglect
      fx =  ((G%IdyCu(I,j)*(Mxx(i+1,j) - Mxx(i,j)) + &
              G%IdxCu(I,j)*((dx2q(I,J)*Mxy(I,J)) - (dx2q(I,J-1)*Mxy(I,J-1)))) * &
              G%IareaCu(I,j)) / h_u
      diffu(I,j,k) = diffu(I,j,k) + fx
      if (save_ZB2020u) &
        ZB2020u(I,j,k) = fx
    enddo ; enddo

    ! Evaluate dv/dt=1/h y.Div(h T) (Line 1517 of MOM_hor_visc.F90)
    do J=Jsq,Jeq ; do i=is,ie
      h_v = 0.5 * (G%mask2dT(i,j)*h(i,j,k) + G%mask2dT(i,j+1)*h(i,j+1,k)) + h_neglect
      fy =  ((G%IdxCv(i,J)*(Myy(i,j+1) - Myy(i,j)) + &
              G%IdyCv(i,J)*((dy2q(I,J)*Mxy(I,J)) - (dy2q(I-1,J)*Mxy(I-1,J)))) * &
              G%IareaCv(i,J)) / h_v
      diffv(i,J,k) = diffv(i,J,k) + fy
      if (save_ZB2020v) &
        ZB2020v(i,J,k) = fy
    enddo ; enddo

  enddo ! end of k loop

  call cpu_clock_end(CS%id_clock_divergence)

  call cpu_clock_begin(CS%id_clock_post)
  if (CS%id_ZB2020u>0)   call post_data(CS%id_ZB2020u, ZB2020u, CS%diag)
  if (CS%id_ZB2020v>0)   call post_data(CS%id_ZB2020v, ZB2020v, CS%diag)
  call cpu_clock_end(CS%id_clock_post)

  call compute_energy_source(u, v, h, ZB2020u, ZB2020v, G, GV, CS)

end procedure compute_stress_divergence
module procedure filter_velocity_gradients
  real, dimension(SZI_(G), SZJ_(G), SZK_(GV)) :: &
        sh_xx          ! Copy of CS%sh_xx [T-1 ~> s-1]
  real, dimension(SZIB_(G),SZJB_(G),SZK_(GV)) :: &
        sh_xy, vort_xy ! Copy of CS%sh_xy and CS%vort_xy [T-1 ~> s-1]
  integer :: xx_halo, xy_halo, vort_halo ! currently available halo for gradient components
  integer :: xx_iter, xy_iter, vort_iter ! remaining number of iterations
  integer :: niter                       ! required number of iterations
  integer :: is, ie, js, je, Isq, Ieq, Jsq, Jeq, nz
  integer :: i, j, k
  niter = CS%HPF_iter

  if (niter == 0) return

  is  = G%isc  ; ie  = G%iec  ; js  = G%jsc  ; je  = G%jec ; nz = GV%ke
  Isq = G%IscB ; Ieq = G%IecB ; Jsq = G%JscB ; Jeq = G%JecB

  if (.not. G%symmetric) &
    call do_group_pass(CS%pass_xx, G%Domain, &
      clock=CS%id_clock_mpi)

  ! This is just copy of the array
  call cpu_clock_begin(CS%id_clock_filter)
  do k=1,nz
    ! Halo of size 2 is valid
    do j=js-2,je+2 ; do i=is-2,ie+2
      sh_xx(i,j,k) = CS%sh_xx(i,j,k)
    enddo ; enddo
    ! Only halo of size 1 is valid
    do J=Jsq-1,Jeq+1 ; do I=Isq-1,Ieq+1
      sh_xy(I,J,k) = CS%sh_xy(I,J,k)
      vort_xy(I,J,k) = CS%vort_xy(I,J,k)
    enddo ; enddo
  enddo
  call cpu_clock_end(CS%id_clock_filter)

  xx_halo = 2 ; xy_halo = 1 ; vort_halo = 1
  xx_iter = niter ; xy_iter = niter ; vort_iter = niter

  do while &
    (xx_iter >  0 .or. xy_iter >  0 .or. & ! filter iterations remain to be done
    xx_halo < 2 .or. xy_halo < 1)        ! there is no halo for VG tensor

    ! ---------- filtering sh_xx ---------
    if (xx_halo < 2) then
      call complete_group_pass(CS%pass_xx, G%Domain, clock=CS%id_clock_mpi)
      xx_halo = CS%HPF_halo
    endif

    call filter_hq(G, GV, CS, xx_halo, xx_iter, h=CS%sh_xx)

    if (xx_halo < 2) &
      call start_group_pass(CS%pass_xx, G%Domain, clock=CS%id_clock_mpi)

    ! ------ filtering sh_xy, vort_xy ----
    if (xy_halo < 1) then
      call complete_group_pass(CS%pass_xy, G%Domain, clock=CS%id_clock_mpi)
      xy_halo = CS%HPF_halo ; vort_halo = CS%HPF_halo
    endif

    call filter_hq(G, GV, CS, xy_halo, xy_iter, q=CS%sh_xy)
    call filter_hq(G, GV, CS, vort_halo, vort_iter, q=CS%vort_xy)

    if (xy_halo < 1) &
      call start_group_pass(CS%pass_xy, G%Domain, clock=CS%id_clock_mpi)

  enddo

  ! We implement sharpening by computing residual
  ! B.C. are already applied to all fields
  call cpu_clock_begin(CS%id_clock_filter)
  do k=1,nz
    do j=js-2,je+2 ; do i=is-2,ie+2
      CS%sh_xx(i,j,k) = sh_xx(i,j,k) - CS%sh_xx(i,j,k)
    enddo ; enddo
    do J=Jsq-1,Jeq+1 ; do I=Isq-1,Ieq+1
      CS%sh_xy(I,J,k) = sh_xy(I,J,k) - CS%sh_xy(I,J,k)
      CS%vort_xy(I,J,k) = vort_xy(I,J,k) - CS%vort_xy(I,J,k)
    enddo ; enddo
  enddo
  call cpu_clock_end(CS%id_clock_filter)

  if (.not. G%symmetric) &
    call do_group_pass(CS%pass_xy, G%Domain, &
      clock=CS%id_clock_mpi)

end procedure filter_velocity_gradients
module procedure filter_stress
  integer :: Txx_halo, Tyy_halo, Txy_halo ! currently available halo for stress components
  integer :: Txx_iter, Tyy_iter, Txy_iter ! remaining number of iterations
  integer :: niter                        ! required number of iterations
  niter = CS%Stress_iter

  if (niter == 0) return

  Txx_halo = 1 ; Tyy_halo = 1 ; Txy_halo = 1 ; ! these are required halo for Txx, Tyy, Txy
  Txx_iter = niter ; Tyy_iter = niter ; Txy_iter = niter

  do while &
      (Txx_iter >  0 .or. Txy_iter >  0 .or. & ! filter iterations remain to be done
       Txx_halo < 1 .or. Txy_halo < 1)         ! there is no halo for Txx or Txy

    ! ---------- filtering Txy -----------
    if (Txy_halo < 1) then
      call complete_group_pass(CS%pass_Tq, G%Domain, clock=CS%id_clock_mpi)
      Txy_halo = CS%Stress_halo
    endif

    call filter_hq(G, GV, CS, Txy_halo, Txy_iter, q=CS%Txy)

    if (Txy_halo < 1) &
       call start_group_pass(CS%pass_Tq, G%Domain, clock=CS%id_clock_mpi)

    ! ------- filtering Txx, Tyy ---------
    if (Txx_halo < 1) then
      call complete_group_pass(CS%pass_Th, G%Domain, clock=CS%id_clock_mpi)
      Txx_halo = CS%Stress_halo ; Tyy_halo = CS%Stress_halo
    endif

    call filter_hq(G, GV, CS, Txx_halo, Txx_iter, h=CS%Txx)
    call filter_hq(G, GV, CS, Tyy_halo, Tyy_iter, h=CS%Tyy)

    if (Txx_halo < 1) &
      call start_group_pass(CS%pass_Th, G%Domain, clock=CS%id_clock_mpi)

  enddo

end procedure filter_stress
module procedure filter_hq
  logical :: direction ! The direction of the first 1D filter
  direction = (MOD(G%first_direction,2) == 0)

  call cpu_clock_begin(CS%id_clock_filter)

  if (present(h)) then
    call filter_3D(h, CS%maskw_h,                  &
              G%isd, G%ied, G%jsd, G%jed,          &
              G%isc, G%iec, G%jsc, G%jec, GV%ke,   &
              current_halo, remaining_iterations,  &
              direction)
  endif

  if (present(q)) then
    call filter_3D(q, CS%maskw_q,                  &
            G%IsdB, G%IedB, G%JsdB, G%JedB,        &
            G%IscB, G%IecB, G%JscB, G%JecB, GV%ke, &
            current_halo, remaining_iterations,    &
            direction)
  endif

  call cpu_clock_end(CS%id_clock_filter)
end procedure filter_hq
module procedure filter_3D
  real, parameter :: weight = 2. ! Filter weight [nondim]
  integer :: i, j, k, iter, niter, halo
  real :: tmp(isd:ied, jsd:jed) ! Array with temporary results [arbitrary]
  niter = min(current_halo, remaining_iterations)
  if (niter == 0) return ! nothing to do

  ! Update remaining iterations
  remaining_iterations = remaining_iterations - niter
  ! Update halo information
  current_halo = current_halo - niter

  do k=1,Nz
    halo = niter-1 + &
      current_halo ! Save as many halo points as possible
    do iter=1,niter

      if (direction) then
        do j = js-halo, je+halo ; do i = is-halo-1, ie+halo+1
          tmp(i,j) = weight * x(i,j,k) + (x(i,j-1,k) + x(i,j+1,k))
        enddo ; enddo

        do j = js-halo, je+halo ; do i = is-halo, ie+halo
          x(i,j,k) = (weight * tmp(i,j) + (tmp(i-1,j) + tmp(i+1,j))) * maskw(i,j)
        enddo ; enddo
      else
        do j = js-halo-1, je+halo+1 ; do i = is-halo, ie+halo
          tmp(i,j) = weight * x(i,j,k) + (x(i-1,j,k) + x(i+1,j,k))
        enddo ; enddo

        do j = js-halo, je+halo ; do i = is-halo, ie+halo
          x(i,j,k) = (weight * tmp(i,j) + (tmp(i,j-1) + tmp(i,j+1))) * maskw(i,j)
        enddo ; enddo
      endif

      halo = halo - 1
    enddo
  enddo

end procedure filter_3D
module procedure compute_energy_source
  real :: KE_term(SZI_(G),SZJ_(G),SZK_(GV)) ! A term in the kinetic energy budget
  real :: KE_u(SZIB_(G),SZJ_(G))            ! The area integral of a KE term in a layer at u-points
  real :: KE_v(SZI_(G),SZJB_(G))            ! The area integral of a KE term in a layer at v-points
  real :: uh                                ! Transport through zonal faces = u*h*dy,
  real :: vh                                ! Transport through meridional faces = v*h*dx,
  type(group_pass_type) :: pass_KE_uv       ! A handle used for group halo passes
  integer :: is, ie, js, je, Isq, Ieq, Jsq, Jeq, nz
  integer :: i, j, k
  if (CS%id_KE_ZB2020 > 0) then
    call cpu_clock_begin(CS%id_clock_source)
    call create_group_pass(pass_KE_uv, KE_u, KE_v, G%Domain, To_North+To_East)

    is  = G%isc  ; ie  = G%iec  ; js  = G%jsc  ; je  = G%jec ; nz = GV%ke
    Isq = G%IscB ; Ieq = G%IecB ; Jsq = G%JscB ; Jeq = G%JecB

    KE_term(:,:,:) = 0.
    ! Calculate the KE source from Zanna-Bolton2020 [H L2 T-3 ~> m3 s-3].
    do k=1,nz
      KE_u(:,:) = 0.
      KE_v(:,:) = 0.
      do j=js,je ; do I=Isq,Ieq
        uh = u(I,j,k) * 0.5 * (G%mask2dT(i,j)*h(i,j,k) + G%mask2dT(i+1,j)*h(i+1,j,k)) * &
          G%dyCu(I,j)
        KE_u(I,j) = uh * G%dxCu(I,j) * fx(I,j,k)
      enddo ; enddo
      do J=Jsq,Jeq ; do i=is,ie
        vh = v(i,J,k) * 0.5 * (G%mask2dT(i,j)*h(i,j,k) + G%mask2dT(i,j+1)*h(i,j+1,k)) * &
          G%dxCv(i,J)
        KE_v(i,J) = vh * G%dyCv(i,J) * fy(i,J,k)
      enddo ; enddo
      call do_group_pass(pass_KE_uv, G%domain, clock=CS%id_clock_mpi)
      do j=js,je ; do i=is,ie
        KE_term(i,j,k) = 0.5 * G%IareaT(i,j) &
            * ((KE_u(I,j) + KE_u(I-1,j)) + (KE_v(i,J) + KE_v(i,J-1)))
      enddo ; enddo
    enddo

    call cpu_clock_end(CS%id_clock_source)

    call cpu_clock_begin(CS%id_clock_post)
    call post_data(CS%id_KE_ZB2020, KE_term, CS%diag)
    call cpu_clock_end(CS%id_clock_post)
  endif

end procedure compute_energy_source
end submodule MOM_Zanna_Bolton_s
