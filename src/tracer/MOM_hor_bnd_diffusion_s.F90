submodule (MOM_hor_bnd_diffusion) MOM_hor_bnd_diffusion_s
#include <MOM_memory.h>
  implicit none
contains
module procedure hor_bnd_diffusion_init
  character(len=80)  :: string ! Temporary strings
  logical :: boundary_extrap   ! controls if boundary extrapolation is used in the HBD code
  logical :: om4_remap_via_sub_cells ! Use the OM4-era ramap_via_sub_cells for HBD
  logical :: debug             !< If true, write verbose checksums for debugging purposes
  if (ASSOCIATED(CS)) then
    call MOM_error(FATAL, "hor_bnd_diffusion_init called with associated control structure.")
    return
  endif

  ! Log this module and master switch for turning it on/off
  call get_param(param_file, mdl, "USE_HORIZONTAL_BOUNDARY_DIFFUSION", hor_bnd_diffusion_init, &
                 default=.false., do_not_log=.true.)
  call log_version(param_file, mdl, version, &
           "This module implements horizontal diffusion of tracers near boundaries", &
           all_default=.not.hor_bnd_diffusion_init)
  call get_param(param_file, mdl, "USE_HORIZONTAL_BOUNDARY_DIFFUSION", hor_bnd_diffusion_init, &
                 "If true, enables the horizontal boundary tracer's diffusion module.", &
                 default=.false.)
  if (.not. hor_bnd_diffusion_init) return

  allocate(CS)
  CS%diag => diag
  CS%H_subroundoff = GV%H_subroundoff
  call extract_diabatic_member(diabatic_CSp, KPP_CSp=CS%KPP_CSp)
  call extract_diabatic_member(diabatic_CSp, energetic_PBL_CSp=CS%energetic_PBL_CSp)

  ! max. number of vertical layers
  CS%hbd_nk = 2 + (GV%ke*2)
  ! allocate the hbd grids and k_max
  allocate(CS%hbd_grd_u(SZIB_(G),SZJ_(G),CS%hbd_nk), source=0.0)
  allocate(CS%hbd_grd_v(SZI_(G),SZJB_(G),CS%hbd_nk), source=0.0)
  allocate(CS%hbd_u_kmax(SZIB_(G),SZJ_(G)), source=0)
  allocate(CS%hbd_v_kmax(SZI_(G),SZJB_(G)), source=0)

  CS%surface_boundary_scheme = -1
  if ( .not. ASSOCIATED(CS%energetic_PBL_CSp) .and. .not. ASSOCIATED(CS%KPP_CSp) ) then
    call MOM_error(FATAL,"Horizontal boundary diffusion is true, but no valid boundary layer scheme was found")
  endif

  ! Read all relevant parameters and write them to the model log.
  call get_param(param_file, mdl, "HBD_LINEAR_TRANSITION", CS%linear, &
                 "If True, apply a linear transition at the base/top of the boundary. \n"//&
                 "The flux will be fully applied at k=k_min and zero at k=k_max.", default=.false.)
  call get_param(param_file, mdl, "APPLY_LIMITER", CS%limiter, &
                   "If True, apply a flux limiter in the native grid.", default=.true.)
  call get_param(param_file, mdl, "APPLY_LIMITER_REMAP", CS%limiter_remap, &
                   "If True, apply a flux limiter in the remapped grid.", default=.false.)
  call get_param(param_file, mdl, "HBD_BOUNDARY_EXTRAP", boundary_extrap, &
                 "Use boundary extrapolation in HBD code", &
                 default=.false.)
  call get_param(param_file, mdl, "HBD_REMAPPING_SCHEME", string, &
                 "This sets the reconstruction scheme used "//&
                 "for vertical remapping for all variables. "//&
                 "It can be one of the following schemes: "//&
                 trim(remappingSchemesDoc), default=remappingDefaultScheme)
  call get_param(param_file, mdl, "REMAPPING_USE_OM4_SUBCELLS", om4_remap_via_sub_cells, &
                 do_not_log=.true., default=.true.)

  call get_param(param_file, mdl, "HBD_REMAPPING_USE_OM4_SUBCELLS", om4_remap_via_sub_cells, &
                 "If true, use the OM4 remapping-via-subcells algorithm for horizontal boundary diffusion. "//&
                 "See REMAPPING_USE_OM4_SUBCELLS for details. "//&
                 "We recommend setting this option to false.", default=om4_remap_via_sub_cells)

  ! GMM, TODO: add HBD params to control optional arguments in initialize_remapping.
  call initialize_remapping( CS%remap_CS, string, boundary_extrapolation=boundary_extrap, &
                             om4_remap_via_sub_cells=om4_remap_via_sub_cells, &
                             check_reconstruction=.false., check_remapping=.false., &
                             h_neglect=CS%H_subroundoff, h_neglect_edge=CS%H_subroundoff)
  call extract_member_remapping_CS(CS%remap_CS, degree=CS%deg)
  call get_param(param_file, mdl, "DEBUG", debug, &
                 default=.false., debuggingParam=.true., do_not_log=.true.)
  call get_param(param_file, mdl, "HBD_DEBUG", CS%debug, &
                 "If true, write out verbose debugging data in the HBD module.", &
                 default=debug, debuggingParam=.true.)

  id_clock_hbd = cpu_clock_id('(Ocean HBD)', grain=CLOCK_MODULE)

end procedure hor_bnd_diffusion_init
module procedure hor_bnd_diffusion
  real, dimension(SZI_(G),SZJ_(G))           :: hbl         !< Boundary layer depth [H ~> m or kg m-2]
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)) :: uFlx        !< Zonal flux of tracer [conc H L2 ~> conc m3 or conc kg]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)) :: vFlx        !< Meridional flux of tracer
  real, dimension(SZIB_(G),SZJ_(G))          :: uwork_2d    !< Layer summed u-flux transport
  real, dimension(SZI_(G),SZJB_(G))          :: vwork_2d    !< Layer summed v-flux transport
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV))  :: tendency    !< tendency array for diagnostics at first in
  real, dimension(SZI_(G),SZJ_(G))           :: tendency_2d !< depth integrated content tendency for diagnostics in
  type(tracer_type), pointer                 :: tracer => NULL() !< Pointer to the current tracer [conc]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV))  :: tracer_old  !< local copy of the initial tracer concentration,
  real :: tracer_int_prev !< Globally integrated tracer before HBD is applied, in mks units [conc kg]
  real :: tracer_int_end  !< Integrated tracer after HBD is applied, in mks units [conc kg]
  real    :: Idt          !< inverse of the time step [T-1 ~> s-1]
  character(len=256) :: mesg !< Message for error messages.
  integer :: i, j, k, m   !< indices to loop over
  call cpu_clock_begin(id_clock_hbd)
  Idt = 1./dt

  if (associated(visc%h_ML)) then
    hbl(:,:) = visc%h_ML(:,:)
  else
    call MOM_error(FATAL, "hor_bnd_diffusion requires that visc%h_ML is associated.")
  endif
  ! This halo update is probably not necessary because visc%h_ML has valid halo data.
  call pass_var(hbl, G%Domain, halo=1)

  ! build HBD grid
  call hbd_grid(SURFACE, G, GV, hbl, h, CS)

  do m = 1,Reg%ntr
    ! current tracer
    tracer => Reg%tr(m)

    if (CS%debug) then
      call hchksum(tracer%t, "before HBD "//tracer%name, G%HI, unscale=tracer%conc_scale)
    endif

    ! for diagnostics
    if (tracer%id_hbdxy_conc > 0 .or. tracer%id_hbdxy_cont > 0 .or. tracer%id_hbdxy_cont_2d > 0 .or. CS%debug) then
      tendency(:,:,:) = 0.0
      tracer_old(:,:,:) = tracer%t(:,:,:)
    endif

    ! Diffusive fluxes in the i- and j-direction
    uFlx(:,:,:) = 0.
    vFlx(:,:,:) = 0.

    ! HBD layer by layer
    do j=G%jsc,G%jec
      do i=G%isc-1,G%iec
        if (G%mask2dCu(I,j)>0.) then
           call fluxes_layer_method(SURFACE, GV%ke, hbl(I,j), hbl(I+1,j),  &
            h(I,j,:), h(I+1,j,:), tracer%t(I,j,:), tracer%t(I+1,j,:), &
            Coef_x(I,j,:), uFlx(I,j,:), G%areaT(I,j), G%areaT(I+1,j), CS%hbd_u_kmax(I,j), &
            CS%hbd_grd_u(I,j,:), CS)
        endif
      enddo
    enddo
    do J=G%jsc-1,G%jec
      do i=G%isc,G%iec
        if (G%mask2dCv(i,J)>0.) then
          call fluxes_layer_method(SURFACE, GV%ke, hbl(i,J), hbl(i,J+1),  &
            h(i,J,:), h(i,J+1,:), tracer%t(i,J,:), tracer%t(i,J+1,:), &
            Coef_y(i,J,:), vFlx(i,J,:), G%areaT(i,J), G%areaT(i,J+1), CS%hbd_v_kmax(i,J), &
            CS%hbd_grd_v(i,J,:), CS)
        endif
      enddo
    enddo

    ! Update the tracer fluxes
    do k=1,GV%ke ; do j=G%jsc,G%jec ; do i=G%isc,G%iec
      if (G%mask2dT(i,j)>0.) then
        tracer%t(i,j,k) = tracer%t(i,j,k) + (( (uFlx(I-1,j,k)-uFlx(I,j,k)) ) + ( (vFlx(i,J-1,k)-vFlx(i,J,k) ) ))* &
                          G%IareaT(i,j) / ( h(i,j,k) + GV%H_subroundoff )

        if (tracer%id_hbdxy_conc > 0  .or. tracer%id_hbdxy_cont > 0 .or. tracer%id_hbdxy_cont_2d > 0 ) then
          tendency(i,j,k) = ((uFlx(I-1,j,k)-uFlx(I,j,k)) + (vFlx(i,J-1,k)-vFlx(i,J,k)))  * &
                            G%IareaT(i,j) * Idt
        endif
      endif
    enddo ; enddo ; enddo

    ! Do user controlled underflow of the tracer concentrations.
    if (tracer%conc_underflow > 0.0) then
      do k=1,GV%ke ; do j=G%jsc,G%jec ; do i=G%isc,G%iec
        if (abs(tracer%t(i,j,k)) < tracer%conc_underflow) tracer%t(i,j,k) = 0.0
      enddo ; enddo ; enddo
    endif

    if (CS%debug) then
      call hchksum(tracer%t, "after HBD "//tracer%name, G%HI, unscale=tracer%conc_scale)
      ! tracer (native grid) integrated tracer amounts before and after HBD
      tracer_int_prev = global_mass_integral(h, G, GV, tracer_old, scale=tracer%conc_scale)
      tracer_int_end = global_mass_integral(h, G, GV, tracer%t, scale=tracer%conc_scale)
      write(mesg,*) 'Total '//tracer%name//' before/after HBD:', tracer_int_prev, tracer_int_end
      call MOM_mesg(mesg)
    endif

    ! Post the tracer diagnostics
    if (tracer%id_hbd_dfx>0)      call post_data(tracer%id_hbd_dfx, uFlx(:,:,:)*Idt, CS%diag)
    if (tracer%id_hbd_dfy>0)      call post_data(tracer%id_hbd_dfy, vFlx(:,:,:)*Idt, CS%diag)
    if (tracer%id_hbd_dfx_2d>0) then
      uwork_2d(:,:) = 0.
      do k=1,GV%ke ; do j=G%jsc,G%jec ; do I=G%isc-1,G%iec
        uwork_2d(I,j) = uwork_2d(I,j) + (uFlx(I,j,k) * Idt)
      enddo ; enddo ; enddo
      call post_data(tracer%id_hbd_dfx_2d, uwork_2d, CS%diag)
    endif

    if (tracer%id_hbd_dfy_2d>0) then
      vwork_2d(:,:) = 0.
      do k=1,GV%ke ; do J=G%jsc-1,G%jec ; do i=G%isc,G%iec
        vwork_2d(i,J) = vwork_2d(i,J) + (vFlx(i,J,k) * Idt)
      enddo ; enddo ; enddo
      call post_data(tracer%id_hbd_dfy_2d, vwork_2d, CS%diag)
    endif

    ! post tendency of tracer content
    if (tracer%id_hbdxy_cont > 0) then
      call post_data(tracer%id_hbdxy_cont, tendency, CS%diag)
    endif

    ! post depth summed tendency for tracer content
    if (tracer%id_hbdxy_cont_2d > 0) then
      tendency_2d(:,:) = 0.
      do j=G%jsc,G%jec ; do i=G%isc,G%iec
        do k=1,GV%ke
          tendency_2d(i,j) = tendency_2d(i,j) + tendency(i,j,k)
        enddo
      enddo ; enddo
      call post_data(tracer%id_hbdxy_cont_2d, tendency_2d, CS%diag)
    endif

    ! post tendency of tracer concentration; this step must be
    ! done after posting tracer content tendency, since we alter
    ! the tendency array and its units.
    if (tracer%id_hbdxy_conc > 0) then
      do k=1,GV%ke ; do j=G%jsc,G%jec ; do i=G%isc,G%iec
        tendency(i,j,k) =  tendency(i,j,k) / ( h(i,j,k) + CS%H_subroundoff )
      enddo ; enddo ; enddo
      call post_data(tracer%id_hbdxy_conc, tendency, CS%diag)
    endif

  enddo

  call cpu_clock_end(id_clock_hbd)

end procedure hor_bnd_diffusion
module procedure hbd_grid
  real, allocatable :: dz_top(:) !< temporary HBD grid given by merge_interfaces           [H ~> m or kg m-2]
  integer :: nk, i, j, k         !< number of layers in the HBD grid, and integers used in do-loops
  CS%hbd_grd_u(:,:,:) = 0.0
  CS%hbd_grd_v(:,:,:) = 0.0
  CS%hbd_u_kmax(:,:)  = 0
  CS%hbd_v_kmax(:,:)  = 0

  do j=G%jsc,G%jec
    do I=G%isc-1,G%iec
      if (G%mask2dCu(I,j)>0.) then
        call merge_interfaces(GV%ke, h(I,j,:), h(I+1,j,:), hbl(I,j), hbl(I+1,j), &
                              CS%H_subroundoff, dz_top)
        nk = SIZE(dz_top)
        if (nk > CS%hbd_nk) then
          write(*,*)'nk, CS%hbd_nk', nk, CS%hbd_nk
          call MOM_error(FATAL,"Houston, we've had a problem in hbd_grid, u-points (nk cannot be > CS%hbd_nk)")
        endif

        CS%hbd_u_kmax(I,j) = nk

        ! set the HBD grid to dz_top
        do k=1,nk
          CS%hbd_grd_u(I,j,k) = dz_top(k)
        enddo
        deallocate(dz_top)
      endif
    enddo
  enddo

  do J=G%jsc-1,G%jec
    do i=G%isc,G%iec
      if (G%mask2dCv(i,J)>0.) then
        call merge_interfaces(GV%ke, h(i,J,:), h(i,J+1,:), hbl(i,J), hbl(i,J+1), &
                              CS%H_subroundoff, dz_top)

        nk = SIZE(dz_top)
        if (nk > CS%hbd_nk) then
          write(*,*)'nk, CS%hbd_nk', nk, CS%hbd_nk
          call MOM_error(FATAL,"Houston, we've had a problem in hbd_grid, v-points (nk cannot be > CS%hbd_nk)")
        endif

        CS%hbd_v_kmax(i,J) = nk

        ! set the HBD grid to dz_top
        do k=1,nk
          CS%hbd_grd_v(i,J,k) = dz_top(k)
        enddo
        deallocate(dz_top)
      endif
    enddo
  enddo

end procedure hbd_grid
module procedure harmonic_mean
  if (h1 + h2 == 0.) then
    harmonic_mean = 0.
  else
    harmonic_mean = 2.*(h1*h2)/(h1+h2)
  endif
end procedure harmonic_mean
module procedure find_minimum
  real :: minimum ! Minimum value in the same units as x [arbitrary]
  integer :: location
  integer :: i
  minimum  = x(s)   ! assume the first is the min
  location = s      ! record its position
  do i = s+1, e     ! start with next elements
    if (x(i) < minimum) then !   if x(i) less than the min?
      minimum  = x(i)   !      Yes, a new minimum found
      location = i                !      record its position
    endif
  enddo
  find_minimum = location          ! return the position
end procedure find_minimum
module procedure swap
  real :: tmp ! A temporary copy of a [arbitrary]
  tmp = a
  a = b
  b = tmp
end procedure swap
module procedure sort
  integer :: i, location
  do i = 1, n-1
    location = find_minimum(x, i, n) ! find min from this to last
    call swap(x(i), x(location))     ! swap this and the minimum
  enddo
end procedure sort
module procedure unique
  real, dimension(n) :: tmp ! The list of unique values [arbitrary]
  integer :: i, j, ii
  real :: min_val, max_val ! The minimum and maximum values in the list [arbitrary]
  logical :: limit
  limit = .false.
  if (present(val_max)) then
    limit = .true.
    if (val_max > MAXVAL(val)) then
      if (is_root_pe()) write(*,*)'val_max, MAXVAL(val)',val_max, MAXVAL(val)
      call MOM_error(FATAL,"Houston, we've had a problem in unique (val_max cannot be > MAXVAL(val))")
    endif
  endif

  tmp(:) = 0.
  min_val = MINVAL(val)-1
  max_val = MAXVAL(val)
  i = 0
  do while (min_val<max_val)
    i = i+1
    min_val = MINVAL(val, mask=val>min_val)
    tmp(i) = min_val
  enddo
  ii = i
  if (limit) then
    do j=1,ii
      if (tmp(j) <= val_max) i = j
    enddo
  endif
  allocate(val_unique(i), source=tmp(1:i))
end procedure unique
module procedure merge_interfaces
  integer                         :: n           !< Number of layers in eta_all
  real, dimension(nk+1)           :: eta_L, eta_R!< Interfaces in the left and right columns [H ~> m or kg m-2]
  real, dimension(:), allocatable :: eta_all     !< Combined list of interfaces in the left and right columns
  real, dimension(:), allocatable :: eta_unique  !< Combined list of unique interfaces (eta_L, eta_R), possibly
  real                            :: min_depth   !< Minimum depth [H ~> m or kg m-2]
  real                            :: max_depth   !< Maximum depth [H ~> m or kg m-2]
  real                            :: max_bld     !< Deepest BLD [H ~> m or kg m-2]
  integer :: k, kk, nk1                          !< loop indices (k and kk) and array size (nk1)
  n = (2*nk)+3
  allocate(eta_all(n))
  ! compute and merge interfaces
  eta_L(:) = 0.0 ; eta_R(:) = 0.0 ; eta_all(:) = 0.0
  kk = 0
  do k=2,nk+1
    eta_L(k) = eta_L(k-1) + h_L(k-1)
    eta_R(k) = eta_R(k-1) + h_R(k-1)
    kk = kk + 2
    eta_all(kk)   = eta_L(k)
    eta_all(kk+1) = eta_R(k)
  enddo

  ! add hbl_L and hbl_R into eta_all
  eta_all(kk+2) = hbl_L
  eta_all(kk+3) = hbl_R

  ! find maximum depth
  min_depth = MIN(MAXVAL(eta_L), MAXVAL(eta_R))
  max_bld = MAX(hbl_L, hbl_R)
  max_depth = MIN(min_depth, max_bld)

  ! sort eta_all
  call sort(eta_all, n)
  ! remove duplicates from eta_all and sets maximum depth
  call unique(eta_all, n, eta_unique, max_depth)

  nk1 = SIZE(eta_unique)
  allocate(h(nk1-1))
  do k=1,nk1-1
    h(k) = (eta_unique(k+1) - eta_unique(k)) + H_subroundoff
  enddo
end procedure merge_interfaces
module procedure flux_limiter
  real :: F_max !< maximum flux allowed [conc H L2 ~> conc m3 or conc kg]
  F_max = -0.2 * ((area_R*(phi_R*h_R))-(area_L*(phi_L*h_L)))

  if ( SIGN(1.,F_layer) == SIGN(1., F_max)) then
    ! Apply flux limiter calculated above
    if (F_max >= 0.) then
      F_layer = MIN(F_layer,F_max)
    else
      F_layer = MAX(F_layer,F_max)
    endif
  else
    F_layer = 0.0
  endif
end procedure flux_limiter
module procedure boundary_k_range
  real :: htot ! Summed thickness [H ~> m or kg m-2]
  integer :: k
  if ( boundary == SURFACE ) then
    k_top = 1
    zeta_top = 0.
    htot = 0.
    k_bot = 1
    zeta_bot = 0.
    if (hbl == 0.) return
    if (hbl >= SUM(h(:))) then
      k_bot = nk
      zeta_bot = 1.
      return
    endif
    do k=1,nk
      htot = htot + h(k)
      if ( htot >= hbl) then
        k_bot = k
        zeta_bot = 1 - (htot - hbl)/h(k)
        return
      endif
    enddo

  ! Bottom boundary layer
  elseif ( boundary == BOTTOM ) then
    k_top = nk
    zeta_top = 1.
    k_bot = nk
    zeta_bot = 0.
    htot = 0.
    if (hbl == 0.) return
    if (hbl >= SUM(h(:))) then
      k_top = 1
      zeta_top = 1.
      return
    endif
    do k=nk,1,-1
      htot = htot + h(k)
      if (htot >= hbl) then
        k_top = k
        zeta_top = 1 - (htot - hbl)/h(k)
        return
      endif
    enddo
  else
    call MOM_error(FATAL,"Houston, we've had a problem in boundary_k_range")
  endif

end procedure boundary_k_range
module procedure fluxes_layer_method
  real, allocatable :: phi_L_z(:)    !< Tracer values in the ztop grid (left)                           [conc]
  real, allocatable :: phi_R_z(:)    !< Tracer values in the ztop grid (right)                          [conc]
  real, allocatable :: F_layer_z(:)  !< Diffusive flux at U/V-point in the ztop grid    [H L2 conc ~> m3 conc]
  real, allocatable :: khtr_ul_z(:)  !< khtr_u at layer centers in the ztop grid        [H L2 conc ~> m3 conc]
  real, dimension(ke) :: h_vel       !< Thicknesses at u- and v-points in the native grid
  real, dimension(ke) :: khtr_ul     !< khtr_u at the vertical layer of the native grid             [L2 ~> m2]
  real    :: htot                    !< Total column thickness                              [H ~> m or kg m-2]
  integer :: k                       !< Index used in the vertical direction
  integer :: k_bot_min               !< Minimum k-index for the bottom
  integer :: k_bot_max               !< Maximum k-index for the bottom
  integer :: k_bot_diff              !< Difference between bottom left and right k-indices
  integer :: k_top_L, k_bot_L        !< k-indices left native grid
  integer :: k_top_R, k_bot_R        !< k-indices right native grid
  real    :: zeta_top_L, zeta_top_R  !< distance from the top of a layer to the boundary
  real    :: zeta_bot_L, zeta_bot_R  !< distance from the bottom of a layer to the boundary
  real    :: wgt                     !< weight to be used in the linear transition to the interior    [nondim]
  real    :: a                       !< coefficient used in the linear transition to the interior     [nondim]
  real    :: tmp1, tmp2              !< dummy variables                                     [H ~> m or kg m-2]
  real    :: htot_max                !< depth below which no fluxes should be applied       [H ~> m or kg m-2]
  F_layer(:) = 0.0
  khtr_ul(:) = 0.0
  if (hbl_L == 0. .or. hbl_R == 0.) then
    return
  endif

  ! allocate arrays
  allocate(phi_L_z(nk), source=0.0)
  allocate(phi_R_z(nk), source=0.0)
  allocate(F_layer_z(nk), source=0.0)
  allocate(khtr_ul_z(nk), source=0.0)

  ! remap tracer to dz_top
  call remapping_core_h(CS%remap_cs, ke, h_L(:), phi_L(:), nk, dz_top(:), phi_L_z(:))
  call remapping_core_h(CS%remap_cs, ke, h_R(:), phi_R(:), nk, dz_top(:), phi_R_z(:))

  ! thicknesses at velocity points & khtr_u at layer centers
  do k = 1,ke
    h_vel(k)   = harmonic_mean(h_L(k), h_R(k))
    ! GMM, writing 0.5 * (A(k) + A(k+1)) as A(k) + 0.5 * (A(k+1) - A(k)) to recover
    ! answers with depth-independent khtr
    khtr_ul(k) = khtr_u(k) + 0.5 * (khtr_u(k+1) - khtr_u(k))
  enddo

  ! remap khtr_ul to khtr_ul_z
  call remapping_core_h(CS%remap_cs, ke, h_vel(:), khtr_ul(:), nk, dz_top(:), khtr_ul_z(:))

  ! Calculate vertical indices containing the boundary layer in dz_top
  call boundary_k_range(boundary, nk, dz_top, hbl_L, k_top_L, zeta_top_L, k_bot_L, zeta_bot_L)
  call boundary_k_range(boundary, nk, dz_top, hbl_R, k_top_R, zeta_top_R, k_bot_R, zeta_bot_R)

  if (boundary == SURFACE) then
    k_bot_min = MIN(k_bot_L, k_bot_R)
    k_bot_max = MAX(k_bot_L, k_bot_R)
    k_bot_diff = (k_bot_max - k_bot_min)

    ! tracer flux where the minimum BLD intersects layer
    if ((CS%linear) .and. (k_bot_diff > 1)) then
      ! apply linear decay at the base of hbl
      do k = k_bot_min,1,-1
        F_layer_z(k) = -(dz_top(k) * khtr_ul_z(k)) * (phi_R_z(k) - phi_L_z(k))
        if (CS%limiter_remap) call flux_limiter(F_layer_z(k), area_L, area_R, phi_L_z(k), &
                                          phi_R_z(k), dz_top(k), dz_top(k))
      enddo
      htot = 0.0
      do k = k_bot_min+1,k_bot_max, 1
        htot = htot + dz_top(k)
      enddo

      a = -1.0/htot
      htot = 0.
      do k = k_bot_min+1,k_bot_max, 1
        wgt = (a*(htot + (dz_top(k) * 0.5))) + 1.0
        F_layer_z(k) = -(dz_top(k) * khtr_ul_z(k)) * (phi_R_z(k) - phi_L_z(k)) * wgt
        htot = htot + dz_top(k)
        if (CS%limiter_remap) call flux_limiter(F_layer_z(k), area_L, area_R, phi_L_z(k), &
                                          phi_R_z(k), dz_top(k), dz_top(k))
      enddo
    else
      do k = k_bot_min,1,-1
        F_layer_z(k) = -(dz_top(k) * khtr_ul_z(k)) * (phi_R_z(k) - phi_L_z(k))
        if (CS%limiter_remap) call flux_limiter(F_layer_z(k), area_L, area_R, phi_L_z(k), &
                                          phi_R_z(k), dz_top(k), dz_top(k))
      enddo
    endif
  endif

  !GMM, TODO: boundary == BOTTOM

  ! remap flux to h_vel (native grid)
  call reintegrate_column(nk, dz_top(:), F_layer_z(:), ke, h_vel(:), F_layer(:))

  ! used to avoid fluxes below hbl
  if (CS%linear) then
    htot_max = MAX(hbl_L, hbl_R)
  else
    htot_max = MIN(hbl_L, hbl_R)
  endif

  tmp1 = 0.0 ; tmp2 = 0.0
  do k = 1,ke
    ! apply flux_limiter
    if (CS%limiter .and. F_layer(k) /= 0.) then
       call flux_limiter(F_layer(k), area_L, area_R, phi_L(k), phi_R(k), h_L(k), h_R(k))
    endif

    ! if tracer point is below htot_max, set flux to zero
    if (MAX(tmp1+(h_L(k)*0.5), tmp2+(h_R(k)*0.5)) > htot_max) then
      F_layer(k) = 0.
    endif

    tmp1 = tmp1 + h_L(k)
    tmp2 = tmp2 + h_R(k)
  enddo

  ! deallocated arrays
  deallocate(phi_L_z)
  deallocate(phi_R_z)
  deallocate(F_layer_z)
  deallocate(khtr_ul_z)

end procedure fluxes_layer_method
module procedure near_boundary_unit_tests
  integer, parameter    :: nk = 2               ! Number of layers
  real, dimension(nk+1) :: eta1                 ! Updated interfaces with one extra value           [m]
  real, dimension(:), allocatable :: h1         ! Updated list of layer thicknesses or other field  [m] or [arbitrary]
  real, dimension(nk)   :: phi_L, phi_R         ! Tracer values (left and right column)             [conc]
  real, dimension(nk)   :: h_L, h_R             ! Layer thickness (left and right)                  [m]
  real, dimension(nk+1) :: khtr_u               ! Horizontal diffusivities at U-point and interfaces[m2 s-1]
  real                  :: hbl_L, hbl_R         ! Depth of the boundary layer (left and right)      [m]
  real, dimension(nk)   :: F_layer              ! Diffusive flux within each layer at U-point       [conc m3 s-1]
  character(len=120)    :: test_name            ! Title of the unit test
  integer               :: k_top                ! Index of cell containing top of boundary
  real                  :: zeta_top             ! Fractional position in the cell of the top        [nondim]
  integer               :: k_bot                ! Index of cell containing bottom of boundary
  real                  :: zeta_bot             ! Fractional position in the cell of the bottom     [nondim]
  type(hbd_CS), pointer :: CS
  allocate(CS)
  ! fill required fields in CS
  CS%linear=.false.
  CS%H_subroundoff = 1.0E-20
  CS%debug=.false.
  CS%limiter=.false.
  CS%limiter_remap=.false.
  CS%hbd_nk = 2 + (2*2)
  call initialize_remapping( CS%remap_CS, 'PLM', boundary_extrapolation=.true., &
                             om4_remap_via_sub_cells=.true., & ! ### see fail below when using fixed remapping alg.
                             check_reconstruction=.true., check_remapping=.true., &
                             h_neglect=CS%H_subroundoff, h_neglect_edge=CS%H_subroundoff)
  call extract_member_remapping_CS(CS%remap_CS, degree=CS%deg)
  allocate(CS%hbd_grd_u(1,1,CS%hbd_nk), source=0.0)
  allocate(CS%hbd_u_kmax(1,1), source=0)
  near_boundary_unit_tests = .false.
  write(stdout,*) '==== MOM_hor_bnd_diffusion ======================='

  ! Unit tests for boundary_k_range
  test_name = 'Surface boundary spans the entire top cell'
  h_L = (/5.,5./)
  call boundary_k_range(SURFACE, nk, h_L, 5., k_top, zeta_top, k_bot, zeta_bot)
  near_boundary_unit_tests = near_boundary_unit_tests .or. &
                             test_boundary_k_range(k_top, zeta_top, k_bot, zeta_bot, 1, 0., 1, 1., test_name, verbose)

  test_name = 'Surface boundary spans the entire column'
  h_L = (/5.,5./)
  call boundary_k_range(SURFACE, nk, h_L, 10., k_top, zeta_top, k_bot, zeta_bot)
  near_boundary_unit_tests = near_boundary_unit_tests .or. &
                             test_boundary_k_range(k_top, zeta_top, k_bot, zeta_bot, 1, 0., 2, 1., test_name, verbose)

  test_name = 'Bottom boundary spans the entire bottom cell'
  h_L = (/5.,5./)
  call boundary_k_range(BOTTOM, nk, h_L, 5., k_top, zeta_top, k_bot, zeta_bot)
  near_boundary_unit_tests = near_boundary_unit_tests .or. &
                             test_boundary_k_range(k_top, zeta_top, k_bot, zeta_bot, 2, 1., 2, 0., test_name, verbose)

  test_name = 'Bottom boundary spans the entire column'
  h_L = (/5.,5./)
  call boundary_k_range(BOTTOM, nk, h_L, 10., k_top, zeta_top, k_bot, zeta_bot)
  near_boundary_unit_tests = near_boundary_unit_tests .or. &
                             test_boundary_k_range(k_top, zeta_top, k_bot, zeta_bot, 1, 1., 2, 0., test_name, verbose)

  test_name = 'Surface boundary intersects second layer'
  h_L = (/10.,10./)
  call boundary_k_range(SURFACE, nk, h_L, 17.5, k_top, zeta_top, k_bot, zeta_bot)
  near_boundary_unit_tests = near_boundary_unit_tests .or. &
                             test_boundary_k_range(k_top, zeta_top, k_bot, zeta_bot, 1, 0., 2, 0.75, test_name, verbose)

  test_name = 'Surface boundary intersects first layer'
  h_L = (/10.,10./)
  call boundary_k_range(SURFACE, nk, h_L, 2.5, k_top, zeta_top, k_bot, zeta_bot)
  near_boundary_unit_tests = near_boundary_unit_tests .or. &
                             test_boundary_k_range(k_top, zeta_top, k_bot, zeta_bot, 1, 0., 1, 0.25, test_name, verbose)

  test_name = 'Surface boundary is deeper than column thickness'
  h_L = (/10.,10./)
  call boundary_k_range(SURFACE, nk, h_L, 21.0, k_top, zeta_top, k_bot, zeta_bot)
  near_boundary_unit_tests = near_boundary_unit_tests .or. &
                             test_boundary_k_range(k_top, zeta_top, k_bot, zeta_bot, 1, 0., 2, 1., test_name, verbose)

  test_name = 'Bottom boundary intersects first layer'
  h_L = (/10.,10./)
  call boundary_k_range(BOTTOM, nk, h_L, 17.5, k_top, zeta_top, k_bot, zeta_bot)
  near_boundary_unit_tests = near_boundary_unit_tests .or. &
                             test_boundary_k_range(k_top, zeta_top, k_bot, zeta_bot, 1, 0.75, 2, 0., test_name, verbose)

  test_name = 'Bottom boundary intersects second layer'
  h_L = (/10.,10./)
  call boundary_k_range(BOTTOM, nk, h_L, 2.5, k_top, zeta_top, k_bot, zeta_bot)
  near_boundary_unit_tests = near_boundary_unit_tests .or. &
                             test_boundary_k_range(k_top, zeta_top, k_bot, zeta_bot, 2, 0.25, 2, 0., test_name, verbose)

  if (.not. near_boundary_unit_tests) write(stdout,*) 'Passed boundary_k_range'

  ! unit tests for sorting array and finding unique values
  test_name = 'Sorting array'
  eta1 = (/1., 0., 0.1/)
  call sort(eta1, nk+1)
  near_boundary_unit_tests = near_boundary_unit_tests .or. &
                             test_layer_fluxes( verbose, nk+1, test_name, eta1, (/0., 0.1, 1./) )

  test_name = 'Unique values'
  call unique((/0., 1., 1., 2./), nk+2, h1)
  near_boundary_unit_tests = near_boundary_unit_tests .or. &
                             test_layer_fluxes( verbose, nk+1, test_name, h1, (/0., 1., 2./) )
  deallocate(h1)

  test_name = 'Unique values with maximum depth'
  call unique((/0., 1., 1., 2., 3./), nk+3, h1, 2.)
  near_boundary_unit_tests = near_boundary_unit_tests .or. &
                             test_layer_fluxes( verbose, nk+1, test_name, h1, (/0., 1., 2./) )
  deallocate(h1)

  if (.not. near_boundary_unit_tests) write(stdout,*) 'Passed sort and unique'

  ! unit tests for merge_interfaces
  test_name = 'h_L = h_R and BLD_L = BLD_R'
  call merge_interfaces(nk, (/1., 2./), (/1., 2./), 1.5, 1.5, CS%H_subroundoff, h1)
  near_boundary_unit_tests = near_boundary_unit_tests .or. &
                             test_layer_fluxes( verbose, nk, test_name, h1, (/1., 0.5/) )
  deallocate(h1)

  test_name = 'h_L = h_R and BLD_L /= BLD_R'
  call merge_interfaces(nk, (/1., 2./), (/1., 2./), 0.5, 1.5, CS%H_subroundoff, h1)
  near_boundary_unit_tests = near_boundary_unit_tests .or. &
                             test_layer_fluxes( verbose, nk+1, test_name, h1, (/0.5, 0.5, 0.5/) )
  deallocate(h1)

  test_name = 'h_L /= h_R and BLD_L = BLD_R'
  call merge_interfaces(nk, (/1., 3./), (/2., 2./), 1.5, 1.5, CS%H_subroundoff, h1)
  near_boundary_unit_tests = near_boundary_unit_tests .or. &
                             test_layer_fluxes( verbose, nk, test_name, h1, (/1., 0.5/) )
  deallocate(h1)

  test_name = 'h_L /= h_R and BLD_L /= BLD_R'
  call merge_interfaces(nk, (/1., 3./), (/2., 2./), 0.5, 1.5, CS%H_subroundoff, h1)
  near_boundary_unit_tests = near_boundary_unit_tests .or. &
                             test_layer_fluxes( verbose, nk+1, test_name, h1, (/0.5, 0.5, 0.5/) )
  deallocate(h1)

  test_name = 'Left deeper than right, h_L /= h_R and BLD_L /= BLD_R'
  call merge_interfaces(nk, (/2., 3./), (/2., 2./), 1.0, 2.0, CS%H_subroundoff, h1)
  near_boundary_unit_tests = near_boundary_unit_tests .or. &
                             test_layer_fluxes( verbose, nk, test_name, h1, (/1., 1./) )
  deallocate(h1)

  test_name = 'Left has zero thickness, h_L /= h_R and BLD_L = BLD_R'
  call merge_interfaces(nk, (/4., 0./), (/2., 2./), 2.0, 2.0, CS%H_subroundoff, h1)
  near_boundary_unit_tests = near_boundary_unit_tests .or. &
                             test_layer_fluxes( verbose, nk-1, test_name, h1, (/2./) )
  deallocate(h1)

  test_name = 'Left has zero thickness, h_L /= h_R and BLD_L /= BLD_R'
  call merge_interfaces(nk, (/4., 0./), (/2., 2./), 1.0, 2.0, CS%H_subroundoff, h1)
  near_boundary_unit_tests = near_boundary_unit_tests .or. &
                             test_layer_fluxes( verbose, nk, test_name, h1, (/1., 1./) )
  deallocate(h1)

  test_name = 'Right has zero thickness, h_L /= h_R and BLD_L = BLD_R'
  call merge_interfaces(nk, (/2., 2./), (/0., 4./), 2.0, 2.0, CS%H_subroundoff, h1)
  near_boundary_unit_tests = near_boundary_unit_tests .or. &
                             test_layer_fluxes( verbose, nk-1, test_name, h1, (/2./) )
  deallocate(h1)

  test_name = 'Right has zero thickness, h_L /= h_R and BLD_L /= BLD_R'
  call merge_interfaces(nk, (/2., 2./), (/0., 4./), 1.0, 2.0, CS%H_subroundoff, h1)
  near_boundary_unit_tests = near_boundary_unit_tests .or. &
                             test_layer_fluxes( verbose, nk, test_name, h1, (/1., 1./) )
  deallocate(h1)

  test_name = 'Right deeper than left, h_L /= h_R and BLD_L = BLD_R'
  call merge_interfaces(nk+1, (/2., 2., 0./), (/2., 2., 1./), 4., 4., CS%H_subroundoff, h1)
  near_boundary_unit_tests = near_boundary_unit_tests .or. &
                             test_layer_fluxes( verbose, nk, test_name, h1, (/2., 2./) )
  deallocate(h1)

  test_name = 'Right and left small values at bottom, h_L /= h_R and BLD_L = BLD_R'
  call merge_interfaces(nk+2, (/2., 2., 1., 1./), (/1., 1., .5, .5/), 3., 3., CS%H_subroundoff, h1)
  near_boundary_unit_tests = near_boundary_unit_tests .or. &
                             test_layer_fluxes( verbose, nk+2, test_name, h1, (/1., 1., .5, .5/) )
  deallocate(h1)

  if (.not. near_boundary_unit_tests) write(stdout,*) 'Passed merge interfaces'

  ! All cases in this section have hbl which are equal to the column thicknesses
  test_name = 'Equal hbl and same layer thicknesses (gradient from right to left)'
  hbl_L = 2. ; hbl_R = 2.
  h_L = (/2.,2./) ; h_R = (/2.,2./)
  phi_L = (/0.,0./) ; phi_R = (/1.,1./)
  khtr_u = (/1.,1.,1./)
  call hbd_grid_test(SURFACE, hbl_L, hbl_R, h_L, h_R, CS)
  call fluxes_layer_method(SURFACE, nk, hbl_L, hbl_R, h_L, h_R, phi_L, phi_R, &
                           khtr_u, F_layer, 1., 1., CS%hbd_u_kmax(1,1), CS%hbd_grd_u(1,1,:), CS)
  near_boundary_unit_tests = near_boundary_unit_tests .or. &
                             test_layer_fluxes( verbose, nk, test_name, F_layer, (/-2.0,0.0/) )

  test_name = 'Equal hbl and same layer thicknesses (gradient from left to right)'
  hbl_L = 2. ; hbl_R = 2.
  h_L = (/2.,2./) ; h_R = (/2.,2./)
  phi_L = (/2.,1./) ; phi_R = (/1.,1./)
  khtr_u = (/0.5,0.5,0.5/)
  call hbd_grid_test(SURFACE, hbl_L, hbl_R, h_L, h_R, CS)
  call fluxes_layer_method(SURFACE, nk, hbl_L, hbl_R, h_L, h_R, phi_L, phi_R, &
                           khtr_u, F_layer, 1., 1., CS%hbd_u_kmax(1,1), CS%hbd_grd_u(1,1,:), CS)
  near_boundary_unit_tests = near_boundary_unit_tests .or. &
                             test_layer_fluxes( verbose, nk, test_name, F_layer, (/1.0,0.0/) )

  test_name = 'hbl < column thickness, hbl same, linear profile right, khtr=2'
  hbl_L = 2 ; hbl_R = 2
  h_L = (/1.,2./) ; h_R = (/1.,2./)
  phi_L = (/0.,0./) ; phi_R = (/0.5,2./)
  khtr_u = (/2.,2.,2./)
  call hbd_grid_test(SURFACE, hbl_L, hbl_R, h_L, h_R, CS)
  call fluxes_layer_method(SURFACE, nk, hbl_L, hbl_R, h_L, h_R, phi_L, phi_R, &
                           khtr_u, F_layer, 1., 1., CS%hbd_u_kmax(1,1), CS%hbd_grd_u(1,1,:), CS)
 ! ### This test fails when om4_remap_via_sub_cells=.false.
  near_boundary_unit_tests = near_boundary_unit_tests .or. &
                             test_layer_fluxes( verbose, nk, test_name, F_layer, (/-1.0,-4.0/) )

  test_name = 'Different hbl and different column thicknesses (zero gradient)'
  hbl_L = 12 ; hbl_R = 20
  h_L = (/6.,6./) ; h_R = (/10.,10./)
  phi_L = (/1.,1./) ; phi_R = (/1.,1./)
  khtr_u = (/1.,1.,1./)
  call hbd_grid_test(SURFACE, hbl_L, hbl_R, h_L, h_R, CS)
  call fluxes_layer_method(SURFACE, nk, hbl_L, hbl_R, h_L, h_R, phi_L, phi_R, &
                           khtr_u, F_layer, 1., 1., CS%hbd_u_kmax(1,1), CS%hbd_grd_u(1,1,:), CS)
  near_boundary_unit_tests = near_boundary_unit_tests .or. &
                             test_layer_fluxes( verbose, nk, test_name, F_layer, (/0.,0./) )

  test_name = 'Different hbl and different column thicknesses (gradient from left to right)'

  hbl_L = 15 ; hbl_R = 10.
  h_L = (/10.,5./) ; h_R = (/10.,0./)
  phi_L = (/1.,1./) ; phi_R = (/0.,0./)
  khtr_u = (/1.,1.,1./)
  call hbd_grid_test(SURFACE, hbl_L, hbl_R, h_L, h_R, CS)
  call fluxes_layer_method(SURFACE, nk, hbl_L, hbl_R, h_L, h_R, phi_L, phi_R, &
                           khtr_u, F_layer, 1., 1., CS%hbd_u_kmax(1,1), CS%hbd_grd_u(1,1,:), CS)
  near_boundary_unit_tests = near_boundary_unit_tests .or. &
                             test_layer_fluxes( verbose, nk, test_name, F_layer, (/10.,0.0/) )

  if (.not. near_boundary_unit_tests) write(stdout,*) 'Passed fluxes_layer_method'

end procedure near_boundary_unit_tests
module procedure test_layer_fluxes
  integer :: k
  test_layer_fluxes = .false.
  do k=1,nk
    if ( F_calc(k) /= F_ans(k) ) then
      test_layer_fluxes = .true.
      write(stdout,*) "MOM_hor_bnd_diffusion, UNIT TEST FAILED: ", test_name
      write(stdout,10) k, F_calc(k), F_ans(k)
    elseif (verbose) then
      write(stdout,10) k, F_calc(k), F_ans(k)
    endif
  enddo

10 format("Layer=",i3," F_calc=",f20.16," F_ans",f20.16)
end procedure test_layer_fluxes
module procedure test_boundary_k_range
  test_boundary_k_range = k_top /= k_top_ans
  test_boundary_k_range = test_boundary_k_range .or. (zeta_top /= zeta_top_ans)
  test_boundary_k_range = test_boundary_k_range .or. (k_bot /= k_bot_ans)
  test_boundary_k_range = test_boundary_k_range .or. (zeta_bot /= zeta_bot_ans)

  if (test_boundary_k_range) write(stdout,*) "UNIT TEST FAILED: ", test_name
  if (test_boundary_k_range .or. verbose) then
    write(stdout,20) "k_top", k_top, "k_top_ans", k_top_ans
    write(stdout,20) "k_bot", k_bot, "k_bot_ans", k_bot_ans
    write(stdout,30) "zeta_top", zeta_top, "zeta_top_ans", zeta_top_ans
    write(stdout,30) "zeta_bot", zeta_bot, "zeta_bot_ans", zeta_bot_ans
  endif

  20 format(A,"=",i3,1X,A,"=",i3)
  30 format(A,"=",f20.16,1X,A,"=",f20.16)


end procedure test_boundary_k_range
module procedure hbd_grid_test
  real, allocatable :: dz_top(:)     !< temporary HBD grid given by merge_interfaces       [H ~> m or kg m-2]
  integer           :: nk, k         !< number of layers in the HBD grid, and integers used in do-loops
  CS%hbd_grd_u(1,1,:) = 0.0
  CS%hbd_u_kmax(1,1)  = 0

  call merge_interfaces(2, h_L, h_R, hbl_L, hbl_R, CS%H_subroundoff, dz_top)
  nk = SIZE(dz_top)
  if (nk > CS%hbd_nk) then
    write(*,*)'nk, CS%hbd_nk', nk, CS%hbd_nk
    call MOM_error(FATAL,"Houston, we've had a problem in hbd_grid_test, (nk cannot be > CS%hbd_nk)")
  endif

  CS%hbd_u_kmax(1,1) = nk

  ! set the HBD grid to dz_top
  do k=1,nk
    CS%hbd_grd_u(1,1,k) = dz_top(k)
  enddo
  deallocate(dz_top)

end procedure hbd_grid_test
module procedure hor_bnd_diffusion_end
  if (associated(CS)) deallocate(CS)

end procedure hor_bnd_diffusion_end
end submodule MOM_hor_bnd_diffusion_s
