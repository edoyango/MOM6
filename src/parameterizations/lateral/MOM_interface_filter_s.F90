submodule (MOM_interface_filter) MOM_interface_filter_s
#include <MOM_memory.h>
  implicit none
contains
module procedure interface_filter
  real :: e(SZI_(G),SZJ_(G),SZK_(GV)+1)  ! Heights of interfaces, relative to mean
  real :: de_smooth(SZI_(G),SZJ_(G),SZK_(GV)+1) ! Change in the heights of interfaces after one pass
  real :: uhD(SZIB_(G),SZJ_(G),SZK_(GV)) ! Smoothing u*h fluxes within a timestep [L2 H ~> m3 or kg]
  real :: vhD(SZI_(G),SZJB_(G),SZK_(GV)) ! Smoothing v*h fluxes within a timestep [L2 H ~> m3 or kg]
  real, dimension(SZIB_(G),SZJ_(G)) :: &
    Lsm2_u        ! Interface height squared smoothing lengths per timestep at u-points [L2 ~> m2]
  real, dimension(SZI_(G),SZJB_(G)) :: &
    Lsm2_v        ! Interface height squared smoothing lengths per timestep at v-points [L2 ~> m2]
  real :: diag_sfn_x(SZIB_(G),SZJ_(G),SZK_(GV)+1)       ! Diagnostic of the x-face streamfunction
  real :: diag_sfn_y(SZI_(G),SZJB_(G),SZK_(GV)+1)       ! Diagnostic of the y-face streamfunction
  real :: filter_strength ! The amount of filtering within a each iteration [nondim]
  real :: Idt       ! The inverse of the timestep [T-1 ~> s-1]
  real :: h_neglect ! A thickness that is so small it is usually lost
  integer :: itt, filter_itts  ! The number of iterations of the filter, set as 1/2 the power.
  integer :: i, j, k, is, ie, js, je, nz, hs
  if (.not. CS%initialized) call MOM_error(FATAL, "MOM_interface_filter: "//&
         "Module must be initialized before it is used.")

  if ((.not.CS%interface_filter) .or. (CS%filter_rate <= 0.0) .or. (CS%filter_order < 2)) return

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke
  h_neglect = GV%H_subroundoff

  filter_itts = CS%filter_order / 2
  Idt = 1.0 / dt

  if (filter_itts > min(G%isc-G%isd, G%jsc-G%jsd)) call MOM_error(FATAL, &
    "interface_filter: The halos are not wide enough to accommodate the filter "//&
    "order specified by INTERFACE_FILTER_ORDER.")

  ! Calculates interface heights, e, in [Z ~> m].
  call find_eta(h, tv, G, GV, US, e, halo_size=filter_itts)

  ! Set the smoothing length scales to apply at each iteration.
  if (filter_itts == 1) then
    filter_strength = min(CS%filter_rate*dt, CS%max_smoothing_CFL)
  elseif (filter_itts == 2) then
    filter_strength = min(sqrt(CS%filter_rate*dt), CS%max_smoothing_CFL)
  else
    filter_strength = min((CS%filter_rate*dt)**(1.0/filter_itts), CS%max_smoothing_CFL)
  endif
  hs = filter_itts-1
  if (CS%isotropic_filter) then
    !$OMP parallel do default(shared)
    do j=js-hs,je+hs ; do I=is-(hs+1),ie+hs
      Lsm2_u(I,j) = (0.25*filter_strength) / ((G%IdxCu(I,j)**2) + (G%IdyCu(I,j)**2))
    enddo ; enddo
    !$OMP parallel do default(shared)
    do J=js-(hs+1),je+hs ; do i=is-hs,ie+hs
      Lsm2_v(i,J) = (0.25*filter_strength) / ((G%IdxCv(i,J)**2) + (G%IdyCv(i,J)**2))
    enddo ; enddo
  else
    !$OMP parallel do default(shared)
    do j=js-hs,je+hs ; do I=is-(hs+1),ie+hs
      Lsm2_u(I,j) = (0.125*filter_strength) *  (min(G%areaT(i,j), G%areaT(i+1,j)) * G%IdyCu(I,j))**2
    enddo ; enddo
    !$OMP parallel do default(shared)
    do J=js-(hs+1),je+hs ; do i=is-hs,ie+hs
      Lsm2_v(i,J) = (0.125*filter_strength) *  (min(G%areaT(i,j), G%areaT(i,j+1)) * G%IdxCv(i,J))**2
    enddo ; enddo
  endif

  if (CS%debug) then
    call uvchksum("Kh_[uv]", Lsm2_u, Lsm2_v, G%HI, haloshift=hs, &
                  unscale=US%L_to_m**2, scalar_pair=.true.)
    call hchksum(h, "interface_filter_1 h", G%HI, haloshift=hs+1, unscale=GV%H_to_m)
    call hchksum(e, "interface_filter_1 e", G%HI, haloshift=hs+1, unscale=US%Z_to_m)
  endif

  ! Calculate uhD, vhD from h, e, Lsm2_u, Lsm2_v
  call filter_interface(h, e, Lsm2_u, Lsm2_v, uhD, vhD, tv, G, GV, US, halo_size=filter_itts-1)


  do itt=2,filter_itts
    hs = (filter_itts - itt) + 1  ! Set the halo size to work on.
    !$OMP parallel do default(shared)
    do j=js-hs,je+hs
      do i=is-hs,ie+hs ; de_smooth(i,j,nz+1) = 0.0 ; enddo

      if (allocated(tv%SpV_avg)) then
        ! This is the fully non-Boussinesq version.
        do k=nz,1,-1 ; do i=is-hs,ie+hs
          de_smooth(i,j,K) = de_smooth(i,j,K+1) + (GV%H_to_RZ * tv%SpV_avg(i,j,k)) * G%IareaT(i,j) * &
              ((uhD(I,j,k) - uhD(I-1,j,k)) + (vhD(i,J,k) - vhD(i,J-1,k)))
        enddo ; enddo
      else
        do k=nz,1,-1 ; do i=is-hs,ie+hs
          de_smooth(i,j,K) = de_smooth(i,j,K+1) + GV%H_to_Z * G%IareaT(i,j) * &
              ((uhD(I,j,k) - uhD(I-1,j,k)) + (vhD(i,J,k) - vhD(i,J-1,k)))
        enddo ; enddo
      endif
    enddo

    ! Calculate uhD, vhD from h, de_smooth, Lsm2_u, Lsm2_v
    call filter_interface(h, de_smooth, Lsm2_u, Lsm2_v, uhD, vhD, tv, G, GV, US, halo_size=filter_itts-itt)
  enddo

  ! Offer diagnostic fields for averaging. This must occur before updating the layer thicknesses
  ! so that the diagnostics can be remapped properly to other diagnostic vertical coordinates.
  if (query_averaging_enabled(CS%diag)) then
    if (CS%id_sfn_x > 0) then
      diag_sfn_x(:,:,1) = 0.0 ; diag_sfn_x(:,:,nz+1) = 0.0
      do K=nz,2,-1 ; do j=js,je ; do I=is-1,ie
        if (CS%id_sfn_x>0) diag_sfn_x(I,j,K) = diag_sfn_x(I,j,K+1) + uhD(I,j,k)
      enddo ; enddo ; enddo
      call post_data(CS%id_sfn_x, diag_sfn_x, CS%diag)
    endif
    if (CS%id_sfn_y > 0) then
      diag_sfn_y(:,:,1) = 0.0 ; diag_sfn_y(:,:,nz+1) = 0.0
      do K=nz,2,-1 ; do J=js-1,je ; do i=is,ie
        diag_sfn_y(i,J,K) = diag_sfn_y(i,J,K+1) + vhD(i,J,k)
      enddo ; enddo ; enddo
      call post_data(CS%id_sfn_y, diag_sfn_y, CS%diag)
    endif
    if (CS%id_uh_sm > 0) call post_data(CS%id_uh_sm, Idt*uhD(:,:,:), CS%diag)
    if (CS%id_vh_sm > 0) call post_data(CS%id_vh_sm, Idt*vhD(:,:,:), CS%diag)
    if (CS%id_L2_u > 0) call post_data(CS%id_L2_u, Lsm2_u, CS%diag)
    if (CS%id_L2_v > 0) call post_data(CS%id_L2_v, Lsm2_v, CS%diag)
  endif

  ! Update the layer thicknesses, and store the transports that will be needed for the tracers.
  !$OMP parallel do default(shared)
  do k=1,nz
    do j=js,je ; do I=is-1,ie
      uhtr(I,j,k) = uhtr(I,j,k) + uhD(I,j,k)
    enddo ; enddo
    do J=js-1,je ; do i=is,ie
      vhtr(i,J,k) = vhtr(i,J,k) + vhD(i,J,k)
    enddo ; enddo
    do j=js,je ; do i=is,ie
      h(i,j,k) = h(i,j,k) - G%IareaT(i,j) * &
          ((uhD(I,j,k) - uhD(I-1,j,k)) + (vhD(i,J,k) - vhD(i,J-1,k)))
      if (h(i,j,k) < GV%Angstrom_H) h(i,j,k) = GV%Angstrom_H
    enddo ; enddo

    ! Store the transports associated with the smoothing if they are needed for diagnostics.
    if (associated(CDp%uh_smooth)) then ; do j=js,je ; do I=is-1,ie
      CDp%uh_smooth(I,j,k) = uhD(I,j,k)*Idt
    enddo ; enddo ; endif
    if (associated(CDp%vh_smooth)) then ; do J=js-1,je ; do i=is,ie
      CDp%vh_smooth(i,J,k) = vhD(i,J,k)*Idt
    enddo ; enddo ; endif

  enddo

  if (CS%debug) then
    call uvchksum("interface_filter [uv]hD", uhD, vhD, &
                  G%HI, haloshift=0, unscale=GV%H_to_m*US%L_to_m**2)
    call uvchksum("interface_filter [uv]htr", uhtr, vhtr, &
                  G%HI, haloshift=0, unscale=US%L_to_m**2*GV%H_to_m)
    call hchksum(h, "interface_filter h", G%HI, haloshift=0, unscale=GV%H_to_m)
  endif

end procedure interface_filter
module procedure filter_interface
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)) :: &
    h_avail       ! The mass available for diffusion out of each face [H L2 ~> m3 or kg].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)+1) :: &
    h_avail_rsum  ! The running sum of h_avail above an interface [H L2 ~> m3 or kg].
  real :: uhtot(SZIB_(G),SZJ_(G))  ! The vertical sum of uhD [H L2 ~> m3 or kg].
  real :: vhtot(SZI_(G),SZJB_(G))  ! The vertical sum of vhD [H L2 ~> m3 or kg].
  real :: Slope         ! The slope of density surfaces, calculated in a way that is always
  real :: Sfn_est       ! A preliminary estimate (before limiting) of the overturning
  real :: Sfn           ! The overturning streamfunction [H L2 ~> m3 or kg].
  real :: Rho_avg       ! The in situ density averaged to an interface [R ~> kg m-3]
  real :: h_neglect     ! A thickness that is so small it is usually lost
  real :: hn_2          ! Half of h_neglect [H ~> m or kg m-2].
  integer :: i, j, k, is, ie, js, je, nz, hs
  hs = 0 ; if (present(halo_size)) hs = halo_size
  is = G%isc-hs ; ie = G%iec+hs ; js = G%jsc-hs ; je = G%jec+hs ; nz = GV%ke

  h_neglect = GV%H_subroundoff ; hn_2 = 0.5*h_neglect

  ! Find the maximum and minimum permitted streamfunction.
  !$OMP parallel do default(shared)
  do j=js-1,je+1
    do i=is-1,ie+1
      h_avail_rsum(i,j,1) = 0.0
      h_avail(i,j,1) = max(0.25*G%areaT(i,j)*(h(i,j,1)-GV%Angstrom_H),0.0)
      h_avail_rsum(i,j,2) = h_avail(i,j,1)
    enddo
    do k=2,nz ; do i=is-1,ie+1
      h_avail(i,j,k) = max(0.25*G%areaT(i,j)*(h(i,j,k)-GV%Angstrom_H),0.0)
      h_avail_rsum(i,j,k+1) = h_avail_rsum(i,j,k) + h_avail(i,j,k)
    enddo ; enddo
  enddo

  !$OMP parallel do default(shared) private(Slope,Sfn_est,Sfn)
  do j=js,je
    do I=is-1,ie ; uhtot(I,j) = 0.0 ; enddo
    do K=nz,2,-1
      do I=is-1,ie
        Slope = (e(i,j,K)-e(i+1,j,K)) * G%IdxCu_OBCmask(I,j)

        if (allocated(tv%SpV_avg)) then
          ! This is the fully non-Boussinesq version.
          Rho_avg = ( ((h(i,j,k) + h(i,j,k-1)) + (h(i+1,j,k) + h(i+1,j,k-1))) + 4.0*hn_2 ) / &
                ( ((h(i,j,k)+hn_2) * tv%SpV_avg(i,j,k)   + (h(i,j,k-1)+hn_2) * tv%SpV_avg(i,j,k-1)) + &
                  ((h(i+1,j,k)+hn_2)*tv%SpV_avg(i+1,j,k) + (h(i+1,j,k-1)+hn_2)*tv%SpV_avg(i+1,j,k-1)) )
          Sfn_est = (Lsm2_u(I,j)*G%dy_Cu(I,j)) * (GV%RZ_to_H * Slope) * Rho_avg
        else
          Sfn_est = (Lsm2_u(I,j)*G%dy_Cu(I,j)) * (GV%Z_to_H * Slope)
        endif

        ! Make sure that there is enough mass above to allow the streamfunction
        ! to satisfy the boundary condition of 0 at the surface.
        Sfn = min(max(Sfn_est, -h_avail_rsum(i,j,K)), h_avail_rsum(i+1,j,K))

        ! The actual transport is limited by the mass available in the two
        ! neighboring grid cells.
        uhD(I,j,k) = max(min((Sfn - uhtot(I,j)), h_avail(i,j,k)), &
                         -h_avail(i+1,j,k))

        ! sfn_x(I,j,K) = max(min(Sfn, uhtot(I,j)+h_avail(i,j,k)), &
        !                    uhtot(I,j)-h_avail(i+1,j,K))

        uhtot(I,j) = uhtot(I,j) + uhD(I,j,k)

      enddo
    enddo ! end of k-loop

    ! In layer 1, enforce the boundary conditions that Sfn(z=0) = 0.0
    do I=is-1,ie ; uhD(I,j,1) = -uhtot(I,j) ; enddo
  enddo ! end of j-loop

  ! Calculate the meridional fluxes and gradients.

  !$OMP parallel do default(shared) private(Slope,Sfn_est,Sfn)
  do J=js-1,je
    do i=is,ie ; vhtot(i,J) = 0.0 ; enddo
    do K=nz,2,-1
      do i=is,ie
        Slope = (e(i,j,K)-e(i,j+1,K)) * G%IdyCv_OBCmask(i,J)

        if (allocated(tv%SpV_avg)) then
          ! This is the fully non-Boussinesq version.
          Rho_avg = ( ((h(i,j,k) + h(i,j,k-1)) + (h(i,j+1,k) + h(i,j+1,k-1))) + 4.0*hn_2 ) / &
                ( ((h(i,j,k)+hn_2) * tv%SpV_avg(i,j,k)   + (h(i,j,k-1)+hn_2) * tv%SpV_avg(i,j,k-1)) + &
                  ((h(i,j+1,k)+hn_2)*tv%SpV_avg(i,j+1,k) + (h(i,j+1,k-1)+hn_2)*tv%SpV_avg(i,j+1,k-1)) )
          Sfn_est = (Lsm2_v(i,J)*G%dx_Cv(i,J)) * (GV%RZ_to_H * Slope) * Rho_avg
        else
          Sfn_est = (Lsm2_v(i,J)*G%dx_Cv(i,J)) * (GV%Z_to_H * Slope)
        endif

        ! Make sure that there is enough mass above to allow the streamfunction
        ! to satisfy the boundary condition of 0 at the surface.
        Sfn = min(max(Sfn_est, -h_avail_rsum(i,j,K)), h_avail_rsum(i,j+1,K))

        ! The actual transport is limited by the mass available in the two neighboring grid cells.
        vhD(i,J,k) = max(min((Sfn - vhtot(i,J)), h_avail(i,j,k)), -h_avail(i,j+1,k))

        ! sfn_y(i,J,K) = max(min(Sfn, vhtot(i,J)+h_avail(i,j,k)), &
        !                    vhtot(i,J)-h_avail(i,j+1,k))

        vhtot(i,J) = vhtot(i,J)  + vhD(i,J,k)

      enddo
    enddo ! end of k-loop
    ! In layer 1, enforce the boundary conditions that Sfn(z=0) = 0.0
    do i=is,ie ; vhD(i,J,1) = -vhtot(i,J) ; enddo
  enddo ! end of j-loop

end procedure filter_interface
module procedure interface_filter_init
  character(len=40)  :: mdl = "MOM_interface_filter" ! This module's name.
# include "version_variable.h"
  real :: interface_filter_time   ! The grid-scale interface height filtering timescale [T ~> s]
  CS%initialized = .true.
  CS%diag => diag

  ! Read all relevant parameters and write them to the model log.
  call log_version(param_file, mdl, version, "")
  call get_param(param_file, mdl, "INTERFACE_FILTER_TIME", interface_filter_time, &
                 "If positive, interface heights are subjected to a grid-scale "//&
                 "dependent biharmonic filter, using a rate based on this timescale.", &
                 default=0.0, units="s", scale=US%s_to_T)
  CS%filter_rate = 0.0
  if (interface_filter_time > 0.0) CS%filter_rate = 1.0 / interface_filter_time
  CS%interface_filter  = (interface_filter_time > 0.0)
  call get_param(param_file, mdl, "INTERFACE_FILTER_MAX_CFL", CS%max_smoothing_CFL, &
                 "The maximum value of the local CFL ratio that "//&
                 "is permitted for the interface height smoothing. 1.0 is the "//&
                 "marginally unstable value.", units="nondimensional", default=0.8)
  if (CS%max_smoothing_CFL < 0.0) CS%max_smoothing_CFL = 0.0

  call get_param(param_file, mdl, "INTERFACE_FILTER_ORDER", CS%filter_order, &
                 "The even power of the interface height smoothing.  "//&
                 "At present valid values are 0, 2, 4 or 6.", default=4)
  if (CS%filter_order == 0) then
    CS%filter_rate = 0.0
  elseif ((CS%filter_order /= 2) .and. (CS%filter_order /= 4) .and. (CS%filter_order /= 6)) then
    call MOM_error(FATAL, "Unsupported value of INTERFACE_FILTER_ORDER specified.  "//&
                          "Only 0, 2, 4 or 6 are supported.")
  endif
  call get_param(param_file, mdl, "INTERFACE_FILTER_ISOTROPIC", CS%isotropic_filter, &
                 "If true, use the same filtering lengthscales in both directions; "//&
                 "otherwise use filtering lengthscales in each direction that scale "//&
                 "with the grid spacing in that direction.", default=.true.)

  call get_param(param_file, mdl, "DEBUG", CS%debug, &
                 "If true, write out verbose debugging data.", &
                 default=.false., debuggingParam=.true.)

  if (CS%filter_order > 0) then
    CS%id_uh_sm = register_diag_field('ocean_model', 'uh_smooth', diag%axesCuL, Time, &
           'Interface Smoothing Zonal Thickness Flux', &
           'kg s-1', conversion=GV%H_to_kg_m2*US%L_to_m**2*US%s_to_T, &
           y_cell_method='sum', v_extensive=.true.)
    CS%id_vh_sm = register_diag_field('ocean_model', 'vh_smooth', diag%axesCvL, Time, &
           'Interface Smoothing Meridional Thickness Flux', &
           'kg s-1', conversion=GV%H_to_kg_m2*US%L_to_m**2*US%s_to_T, &
           x_cell_method='sum', v_extensive=.true.)

    CS%id_L2_u = register_diag_field('ocean_model', 'Lsmooth2_u', diag%axesCu1, Time, &
           'Interface height smoothing length-scale squared at U-points', &
           'm2', conversion=US%L_to_m**2)
    CS%id_L2_v = register_diag_field('ocean_model', 'Lsmooth2_u', diag%axesCv1, Time, &
           'Interface height smoothing length-scale squared at V-points', &
           'm2', conversion=US%L_to_m**2)

    CS%id_sfn_x =  register_diag_field('ocean_model', 'Smooth_sfn_x', diag%axesCui, Time, &
           'Interface smoothing Zonal Overturning Streamfunction', &
           'm3 s-1', conversion=GV%H_to_m*US%L_to_m**2*US%s_to_T)
    CS%id_sfn_y =  register_diag_field('ocean_model', 'Smooth_sfn_y', diag%axesCvi, Time, &
           'Interface smoothing Meridional Overturning Streamfunction', &
           'm3 s-1', conversion=GV%H_to_m*US%L_to_m**2*US%s_to_T)
  endif

end procedure interface_filter_init
module procedure interface_filter_end
  if (associated(CDp%uh_smooth)) deallocate(CDp%uh_smooth)
  if (associated(CDp%vh_smooth)) deallocate(CDp%vh_smooth)

end procedure interface_filter_end
end submodule MOM_interface_filter_s
