submodule (MOM_interface_heights) MOM_interface_heights_s
#include <MOM_memory.h>
  implicit none
contains
module procedure find_dz_for_eta
  real :: p(SZI_(G),SZJ_(G),SZK_(GV)+1)   ! Hydrostatic pressure at each interface [R L2 T-2 ~> Pa]
  real :: dz_geo(SZI_(G),SZJ_(G)) ! The change in geopotential height across a layer [L2 T-2 ~> m2 s-2]
  real :: SpV_lay_conv(SZK_(GV))  ! The prescribed layer specific volume times a conversion factor from
  real :: I_gEarth                ! The inverse of the gravitational acceleration times the
  integer :: i, j, k, isv, iev, jsv, jev, nz, halo
  halo = 0 ; if (present(halo_size)) halo = max(0,halo_size)

  isv = G%isc-halo ; iev = G%iec+halo ; jsv = G%jsc-halo ; jev = G%jec+halo
  nz = GV%ke

  if ((isv<G%isd) .or. (iev>G%ied) .or. (jsv<G%jsd) .or. (jev>G%jed)) &
    call MOM_error(FATAL,"find_dz_for_eta called with an overly large halo_size.")

  if (GV%Boussinesq) then
    do k=1,nz ; do j=jsv,jev ; do i=isv,iev
      dz_lay(i,j,K) = h(i,j,k)*GV%H_to_Z
    enddo ; enddo ; enddo
  elseif (associated(tv%eqn_of_state)) then
    I_gEarth = 1.0 / GV%g_Earth
    !$OMP parallel do default(shared)
    do j=jsv,jev
      if (associated(tv%p_surf)) then
        do i=isv,iev ; p(i,j,1) = tv%p_surf(i,j) ; enddo
      else
        do i=isv,iev ; p(i,j,1) = 0.0 ; enddo
      endif
      do k=1,nz ; do i=isv,iev
        p(i,j,K+1) = p(i,j,K) + GV%g_Earth*GV%H_to_RZ*h(i,j,k)
      enddo ; enddo
    enddo
    !$OMP parallel do default(shared) private(dz_geo)
    do k=1,nz
      call int_specific_vol_dp(tv%T(:,:,k), tv%S(:,:,k), p(:,:,K), p(:,:,K+1), &
                               0.0, G%HI, tv%eqn_of_state, US, dz_geo, halo_size=halo)
      do j=jsv,jev ; do i=isv,iev
        dz_lay(i,j,K) = I_gEarth * dz_geo(i,j)
      enddo ; enddo
    enddo
  else ! non-Boussinesq but with no equation of state
    do k=1,nz ; do j=jsv,jev ; do i=isv,iev
      dz_lay(i,j,K) = GV%H_to_RZ*h(i,j,k) / GV%Rlay(k)
    enddo ; enddo ; enddo
    ! This would be faster but could change answers.
    ! do k=1,nz ; SpV_lay_conv(k) = GV%H_to_RZ / GV%Rlay(k) ; enddo
    ! do k=1,nz ; do j=jsv,jev ; do i=isv,iev
    !   dz_lay(i,j,K) = h(i,j,k) * SpV_lay_conv(k)
    ! enddo ; enddo ; enddo
  endif

  ! To find eta, do the following:
  ! do j=jsv,jev ; do i=isv,iev ; eta(i,j,nz+1) = -(G%bathyT(i,j) + dZ_ref) ; enddo ; enddo
  ! do k=nz,1,-1 ; do j=jsv,jev ; do i=isv,iev
  !   eta(i,j,K) = eta(i,j,K+1) + dz_lay(i,j,K)
  ! enddo ; enddo ; enddo

end procedure find_dz_for_eta
module procedure find_eta_3d
  real :: dz_lay(SZI_(G),SZJ_(G),SZK_(GV)) ! The change in height across a layer [Z ~> m]
  real :: dilate(SZI_(G))                 ! A non-dimensional dilation factor [nondim]
  real :: htot(SZI_(G))                   ! total thickness [H ~> m or kg m-2]
  real :: dZ_ref    ! The difference in the reference height between G%bathyT and eta [Z ~> m].
  integer :: i, j, k, isv, iev, jsv, jev, nz, halo
  halo = 0 ; if (present(halo_size)) halo = max(0,halo_size)

  isv = G%isc-halo ; iev = G%iec+halo ; jsv = G%jsc-halo ; jev = G%jec+halo
  nz = GV%ke

  if ((isv<G%isd) .or. (iev>G%ied) .or. (jsv<G%jsd) .or. (jev>G%jed)) &
    call MOM_error(FATAL,"find_eta called with an overly large halo_size.")

  dZ_ref = 0.0 ; if (present(dZref)) dZ_ref = dZref

  if (GV%Boussinesq) then
    !$OMP parallel default(shared) private(dilate,htot)
    !$OMP do
    do j=jsv,jev ; do i=isv,iev ; eta(i,j,nz+1) = -(G%bathyT(i,j) + dZ_ref) ; enddo ; enddo
    !$OMP do
    do j=jsv,jev ; do k=nz,1,-1 ; do i=isv,iev
      eta(i,j,K) = eta(i,j,K+1) + h(i,j,k)*GV%H_to_Z
    enddo ; enddo ; enddo
    if (present(eta_bt)) then
      ! Dilate the water column to agree with the free surface height
      ! that is used for the dynamics.
      !$OMP do
      do j=jsv,jev    !$OMP parallel do default(shared)

        do i=isv,iev
          dilate(i) = (eta_bt(i,j)*GV%H_to_Z + G%bathyT(i,j)) / &
                      (eta(i,j,1) + (G%bathyT(i,j) + dZ_ref))
        enddo
        do k=1,nz ; do i=isv,iev
          eta(i,j,K) = dilate(i) * (eta(i,j,K) + (G%bathyT(i,j) + dZ_ref)) - &
                       (G%bathyT(i,j) + dZ_ref)
        enddo ; enddo
      enddo
    endif
    !$OMP end parallel
  else
    call find_dz_for_eta(h, tv, G, GV, US, dz_lay, halo_size)
    !$OMP parallel default(shared) private(dilate,htot)
    !$OMP do
    do j=jsv,jev
      do i=isv,iev ; eta(i,j,nz+1) = -(G%bathyT(i,j) + dZ_ref) ; enddo
      do k=nz,1,-1 ; do i=isv,iev
        eta(i,j,K) = eta(i,j,K+1) + dz_lay(i,j,k)
      enddo ; enddo
    enddo

    if (present(eta_bt)) then
      ! Dilate the water column to agree with the free surface height
      ! from the time-averaged barotropic solution.
      !$OMP do
      do j=jsv,jev
        do i=isv,iev ; htot(i) = GV%H_subroundoff ; enddo
        do k=1,nz ; do i=isv,iev ; htot(i) = htot(i) + h(i,j,k) ; enddo ; enddo
        do i=isv,iev ; dilate(i) = eta_bt(i,j) / htot(i) ; enddo
        do k=1,nz ; do i=isv,iev
          eta(i,j,K) = dilate(i) * (eta(i,j,K) + (G%bathyT(i,j) + dZ_ref)) - &
                       (G%bathyT(i,j) + dZ_ref)
        enddo ; enddo
      enddo
    endif
    !$OMP end parallel
  endif

end procedure find_eta_3d
module procedure find_eta_2d
  real :: dz_lay(SZI_(G),SZJ_(G),SZK_(GV)) ! The change in height across a layer [Z ~> m]
  real :: htot(SZI_(G))  ! The sum of all layers' thicknesses [H ~> m or kg m-2].
  real :: dZ_ref    ! The difference in the reference height between G%bathyT and eta [Z ~> m].
  integer :: i, j, k, is, ie, js, je, nz, halo
  halo = 0 ; if (present(halo_size)) halo = max(0,halo_size)
  is = G%isc-halo ; ie = G%iec+halo ; js = G%jsc-halo ; je = G%jec+halo
  nz = GV%ke

  dZ_ref = 0.0 ; if (present(dZref)) dZ_ref = dZref

  if (GV%Boussinesq) then
    if (present(eta_bt)) then
      !$OMP parallel do default(shared)
      do j=js,je ; do i=is,ie
        eta(i,j) = GV%H_to_Z*eta_bt(i,j) - dZ_ref
      enddo ; enddo
    else
      !$OMP parallel do default(shared)
      do j=js,je
        do i=is,ie ; eta(i,j) = -(G%bathyT(i,j) + dZ_ref) ; enddo
        do k=1,nz ; do i=is,ie
          eta(i,j) = eta(i,j) + h(i,j,k)*GV%H_to_Z
        enddo ; enddo
      enddo
    endif
  else
    call find_dz_for_eta(h, tv, G, GV, US, dz_lay, halo_size)
    !$OMP parallel default(shared) private(htot)
    !$OMP do
    do j=js,je
      do i=is,ie ; eta(i,j) = -(G%bathyT(i,j) + dZ_ref) ; enddo
      do k=1,nz ; do i=is,ie
        eta(i,j) = eta(i,j) + dz_lay(i,j,k)
      enddo ; enddo
    enddo
    if (present(eta_bt)) then
      !   Dilate the water column to agree with the time-averaged column
      ! mass from the barotropic solution.
      !$OMP do
      do j=js,je
        do i=is,ie ; htot(i) = GV%H_subroundoff ; enddo
        do k=1,nz ; do i=is,ie ; htot(i) = htot(i) + h(i,j,k) ; enddo ; enddo
        do i=is,ie
          eta(i,j) = (eta_bt(i,j) / htot(i)) * (eta(i,j) + (G%bathyT(i,j) + dZ_ref)) - &
                     (G%bathyT(i,j) + dZ_ref)
        enddo
      enddo
    endif
    !$OMP end parallel
  endif

end procedure find_eta_2d
module procedure calc_derived_thermo
  real, dimension(SZI_(G),SZJ_(G)) :: p_t  ! Hydrostatic pressure atop a layer [R L2 T-2 ~> Pa]
  real, dimension(SZI_(G),SZJ_(G)) :: dp   ! Pressure change across a layer [R L2 T-2 ~> Pa]
  real, dimension(SZK_(GV)) :: SpV_lay     ! The specific volume of each layer when no equation of
  logical :: do_debug  ! If true, write checksums for debugging.
  integer :: i, j, k, is, ie, js, je, halos, nz
  do_debug = .false. ; if (present(debug)) do_debug = debug
  halos = 0 ; if (present(halo)) halos = max(0,halo)
  is = G%isc-halos ; ie = G%iec+halos ; js = G%jsc-halos ; je = G%jec+halos ; nz = GV%ke

  if (allocated(tv%Spv_avg) .and. associated(tv%eqn_of_state)) then
    if (associated(tv%p_surf)) then
      do j=js,je ; do i=is,ie ; p_t(i,j) = tv%p_surf(i,j) ; enddo ; enddo
    else
      do j=js,je ; do i=is,ie ; p_t(i,j) = 0.0 ; enddo ; enddo
    endif
    do k=1,nz
      do j=js,je ; do i=is,ie
        dp(i,j) = GV%g_Earth*GV%H_to_RZ*h(i,j,k)
      enddo ; enddo
      call avg_specific_vol(tv%T(:,:,k), tv%S(:,:,k), p_t, dp, G%HI, tv%eqn_of_state, tv%SpV_avg(:,:,k), halo)
      if (k<nz) then ; do j=js,je ; do i=is,ie
        p_t(i,j) = p_t(i,j) + dp(i,j)
      enddo ; enddo ; endif
    enddo
    tv%valid_SpV_halo = halos

    if (do_debug) then
      call hchksum(h, "derived_thermo h", G%HI, haloshift=halos, unscale=GV%H_to_MKS)
      if (associated(tv%p_surf)) call hchksum(tv%p_surf, "derived_thermo p_surf", G%HI, &
                                              haloshift=halos, unscale=US%RL2_T2_to_Pa)
      call hchksum(tv%T, "derived_thermo T", G%HI, haloshift=halos, unscale=US%C_to_degC)
      call hchksum(tv%S, "derived_thermo S", G%HI, haloshift=halos, unscale=US%S_to_ppt)
    endif
  elseif (allocated(tv%Spv_avg)) then
    do k=1,nz ; SpV_lay(k) = 1.0 / GV%Rlay(k) ; enddo
    do k=1,nz ; do j=js,je ; do i=is,ie
      tv%SpV_avg(i,j,k) = SpV_lay(k)
    enddo ; enddo ; enddo
    tv%valid_SpV_halo = halos
  endif

end procedure calc_derived_thermo
module procedure find_col_avg_SpV
  real :: h_tot(SZI_(G))        ! Sum of the layer thicknesses [H ~> m or kg m-2]
  real :: SpV_x_h_tot(SZI_(G))  ! Vertical sum of the layer average specific volume times
  real :: I_rho                 ! The inverse of the Boussiensq reference density [R-1 ~> m3 kg-1]
  real :: SpV_lay(SZK_(GV))     ! The inverse of the layer target potential densities [R-1 ~> m3 kg-1]
  character(len=128) :: mesg    ! A string for error messages
  integer :: i, j, k, is, ie, js, je, nz, halo
  halo = 0 ; if (present(halo_size)) halo = max(0,halo_size)

  is = G%isc-halo ; ie = G%iec+halo ; js = G%jsc-halo ; je = G%jec+halo
  nz = GV%ke

  if (GV%Boussinesq) then
    I_rho = 1.0 / GV%Rho0
    do j=js,je ; do i=is,ie
      SpV_avg(i,j) = I_rho
    enddo ; enddo
  elseif (.not.allocated(tv%SpV_avg)) then
    do k=1,nz ; Spv_lay(k) = 1.0 / GV%Rlay(k) ; enddo
    do j=js,je
      do i=is,ie ; SpV_x_h_tot(i) = 0.0 ; h_tot(i) = 0.0 ; enddo
      do k=1,nz ; do i=is,ie
        h_tot(i) = h_tot(i) + max(h(i,j,k), GV%H_subroundoff)
        SpV_x_h_tot(i) = SpV_x_h_tot(i) + Spv_lay(k)*max(h(i,j,k), GV%H_subroundoff)
      enddo ; enddo
      do i=is,ie ; SpV_avg(i,j) = SpV_x_h_tot(i) / h_tot(i) ; enddo
    enddo
  else
    ! Check that SpV_avg has been set.
    if ((allocated(tv%SpV_avg)) .and. (tv%valid_SpV_halo < halo)) then
      if (tv%valid_SpV_halo < 0) then
        mesg = "invalid values of SpV_avg."
      else
        write(mesg, '("insufficiently large SpV_avg halos of width ", i2, " but ", i2," is needed.")') &
                     tv%valid_SpV_halo, halo
      endif
      call MOM_error(FATAL, "find_col_avg_SpV called in fully non-Boussinesq mode with "//trim(mesg))
    endif

    do j=js,je
      do i=is,ie ; SpV_x_h_tot(i) = 0.0 ; h_tot(i) = 0.0 ; enddo
      do k=1,nz ; do i=is,ie
        h_tot(i) = h_tot(i) + max(h(i,j,k), GV%H_subroundoff)
        SpV_x_h_tot(i) = SpV_x_h_tot(i) + tv%SpV_avg(i,j,k)*max(h(i,j,k), GV%H_subroundoff)
      enddo ; enddo
      do i=is,ie ; SpV_avg(i,j) = SpV_x_h_tot(i) / h_tot(i) ; enddo
    enddo
  endif

end procedure find_col_avg_SpV
module procedure find_col_mass
  real :: I_gEarth ! The inverse of GV%g_Earth [T2 Z L-2 ~> s2 m-1]
  real, dimension(SZI_(G),SZJ_(G)) :: &
    z_top, & ! Height of the top of a layer [Z ~> m].
    z_bot, & ! Height of the bottom of a layer [Z ~> m].
    dp       ! Change in hydrostatic pressure across a layer [R L2 T-2 ~> Pa].
  integer :: i, j, k, is, ie, js, je, isq, ieq, jsq, jeq, nz
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec
  isq = G%iscB ; ieq = G%iecB ; jsq = G%jscB ; jeq = G%jecB
  nz = GV%ke

  do j=js,je ; do i=is,ie ; mass(i,j) = 0.0 ; enddo ; enddo
  if (GV%Boussinesq) then
    if (associated(tv%eqn_of_state)) then
      I_gEarth = 1.0 / GV%g_Earth
      do j=jsq,jeq+1 ; do i=isq,ieq+1 ; z_bot(i,j) = 0.0 ; enddo ; enddo
      do k=1,nz
        ! NOTE: int_density_z expects z_top and z_bot values from [ij]sq to [ij]eq+1
        do j=jsq,jeq+1 ; do i=isq,ieq+1
          z_top(i,j) = z_bot(i,j)
          z_bot(i,j) = z_top(i,j) - GV%H_to_Z * h(i,j,k)
        enddo ; enddo
        call int_density_dz(tv%T(:,:,k), tv%S(:,:,k), z_top, z_bot, 0.0, GV%Rho0, GV%g_Earth, &
                            G%HI, tv%eqn_of_state, US, dp)
        do j=js,je ; do i=is,ie
          mass(i,j) = mass(i,j) + dp(i,j) * I_gEarth
        enddo ; enddo
      enddo
    else
      do k=1,nz ; do j=js,je ; do i=is,ie
        mass(i,j) = mass(i,j) + (GV%H_to_Z * GV%Rlay(k)) * h(i,j,k)
      enddo ; enddo ; enddo
    endif
  else
    do k=1,nz ; do j=js,je ; do i=is,ie
      mass(i,j) = mass(i,j) + GV%H_to_RZ * h(i,j,k)
    enddo ; enddo ; enddo
  endif

  if (present(p_bot)) then
    do j=js,je ; do i=is,ie
      p_bot(i,j) = GV%g_Earth * mass(i,j)
    enddo ; enddo
    if (present(p_surf) .and. associated(p_surf)) then ; do j=js,je ; do i=is,ie
      p_bot(i,j) = p_bot(i,j) + p_surf(i,j)
    enddo ; enddo ; endif
  endif

end procedure find_col_mass
module procedure find_rho_bottom
  real :: hb(SZI_(G))         ! Running sum of the thickness in the bottom boundary layer [H ~> m or kg m-2]
  real :: SpV_h_bot(SZI_(G))  ! Running sum of the specific volume times thickness in the bottom
  real :: dz_bbl_rem(SZI_(G)) ! Vertical extent of the boundary layer that has yet to be accounted
  real :: h_bbl_frac(SZI_(G)) ! Thickness of the fractional layer that makes up the top of the
  real :: T_bbl(SZI_(G))      ! Temperature of the fractional layer that makes up the top of the
  real :: S_bbl(SZI_(G))      ! Salinity of the fractional layer that makes up the top of the
  real :: P_bbl(SZI_(G))      ! Pressure the top of the boundary layer [R L2 T-2 ~> Pa]
  real :: dp(SZI_(G))         ! Pressure change across the fractional layer that makes up the top
  real :: SpV_bbl(SZI_(G))    ! In situ specific volume of the fractional layer that makes up the
  real :: frac_in             ! The fraction of a layer that is within the bottom boundary layer [nondim]
  logical :: do_i(SZI_(G)), do_any
  logical :: use_EOS
  integer, dimension(2) :: EOSdom ! The i-computational domain for the equation of state
  integer :: i, k, is, ie, nz
  is = G%isc ; ie = G%iec ; nz = GV%ke

  use_EOS = associated(tv%T) .and. associated(tv%S) .and. associated(tv%eqn_of_state)

  if (GV%Boussinesq .or. GV%semi_Boussinesq .or. .not.allocated(tv%SpV_avg)) then
    do i=is,ie
      rho_bot(i) = GV%Rho0
    enddo

    ! Obtain bottom boundary layer thickness and index of top layer
    do i=is,ie
      hb(i) = 0.0 ; h_bot(i) = 0.0 ; k_bot(i) = nz
      dz_bbl_rem(i) = G%mask2dT(i,j) * max(0.0, dz_avg(i))
      do_i(i) = .true.
      if (G%mask2dT(i,j) <= 0.0) then
        h_bbl_frac(i) = 0.0
        do_i(i) = .false.
      endif
    enddo

    do k=nz,1,-1
      do_any = .false.
      do i=is,ie ; if (do_i(i)) then
        if (dz(i,k) < dz_bbl_rem(i)) then
          ! This layer is fully within the averaging depth.
          dz_bbl_rem(i) = dz_bbl_rem(i) - dz(i,k)
          hb(i) = hb(i) + h(i,j,k)
          k_bot(i) = k
          do_any = .true.
        else
          if (dz(i,k) > 0.0) then
            frac_in = dz_bbl_rem(i) / dz(i,k)
            if (frac_in >= 0.5) k_bot(i) = k ! update bbl top index if >= 50% of layer
          else
            frac_in = 0.0
          endif
          h_bbl_frac(i) = frac_in * h(i,j,k)
          dz_bbl_rem(i) = 0.0
          do_i(i) = .false.
        endif
      endif ; enddo
      if (.not.do_any) exit
    enddo
    do i=is,ie ; if (do_i(i)) then
      ! The nominal bottom boundary layer is thicker than the water column, but layer 1 is
      ! already included in the averages.  These values are set so that the call to find
      ! the layer-average specific volume will behave sensibly.
      h_bbl_frac(i) = 0.0
    endif ; enddo

    do i=is,ie
      if (hb(i) + h_bbl_frac(i) < GV%H_subroundoff) h_bbl_frac(i) = GV%H_subroundoff
      h_bot(i) = hb(i) + h_bbl_frac(i)
    enddo

  else
    ! Check that SpV_avg has been set.
    if (tv%valid_SpV_halo < 0) call MOM_error(FATAL, &
        "find_rho_bottom called in fully non-Boussinesq mode with invalid values of SpV_avg.")

    ! Set the bottom density to the inverse of the in situ specific volume averaged over the
    ! specified distance, with care taken to avoid having compressibility lead to an imprint
    ! of the layer thicknesses on this density.
    do i=is,ie
      hb(i) = 0.0 ; SpV_h_bot(i) = 0.0 ; h_bot(i) = 0.0 ; k_bot(i) = nz
      dz_bbl_rem(i) = G%mask2dT(i,j) * max(0.0, dz_avg(i))
      do_i(i) = .true.
      if (G%mask2dT(i,j) <= 0.0) then
        ! Set acceptable values for calling the equation of state over land.
        T_bbl(i) = 0.0 ; S_bbl(i) = 0.0 ; dp(i) = 0.0 ; P_bbl(i) = 0.0
        SpV_bbl(i) = 1.0 ! This value is arbitrary, provided it is non-zero.
        h_bbl_frac(i) = 0.0
        do_i(i) = .false.
      endif
    enddo

    do k=nz,1,-1
      do_any = .false.
      do i=is,ie ; if (do_i(i)) then
        if (dz(i,k) < dz_bbl_rem(i)) then
          ! This layer is fully within the averaging depth.
          SpV_h_bot(i) = SpV_h_bot(i) + h(i,j,k) * tv%SpV_avg(i,j,k)
          dz_bbl_rem(i) = dz_bbl_rem(i) - dz(i,k)
          hb(i) = hb(i) + h(i,j,k)
          k_bot(i) = k
          do_any = .true.
        else
          if (dz(i,k) > 0.0) then
            frac_in = dz_bbl_rem(i) / dz(i,k)
            if (frac_in >= 0.5) k_bot(i) = k ! update bbl top index if >= 50% of layer
          else
            frac_in = 0.0
          endif
          if (use_EOS) then
            ! Store the properties of this layer to determine the average
            ! specific volume of the portion that is within the BBL.
            T_bbl(i) = tv%T(i,j,k) ; S_bbl(i) = tv%S(i,j,k)
            dp(i) = frac_in * (GV%g_Earth*GV%H_to_RZ * h(i,j,k))
            P_bbl(i) = pres_int(i,K) + (1.0-frac_in) * (GV%g_Earth*GV%H_to_RZ * h(i,j,k))
          else
            SpV_bbl(i) = tv%SpV_avg(i,j,k)
          endif
          h_bbl_frac(i) = frac_in * h(i,j,k)
          dz_bbl_rem(i) = 0.0
          do_i(i) = .false.
        endif
      endif ; enddo
      if (.not.do_any) exit
    enddo
    do i=is,ie ; if (do_i(i)) then
      ! The nominal bottom boundary layer is thicker than the water column, but layer 1 is
      ! already included in the averages.  These values are set so that the call to find
      ! the layer-average specific volume will behave sensibly.
      if (use_EOS) then
        T_bbl(i) = tv%T(i,j,1) ; S_bbl(i) = tv%S(i,j,1)
        dp(i) = 0.0
        P_bbl(i) = pres_int(i,1)
      else
        SpV_bbl(i) = tv%SpV_avg(i,j,1)
      endif
      h_bbl_frac(i) = 0.0
    endif ; enddo

    if (use_EOS) then
      ! Find the average specific volume of the fractional layer atop the BBL.
      EOSdom(:) = EOS_domain(G%HI)
      call average_specific_vol(T_bbl, S_bbl, P_bbl, dp, SpV_bbl, tv%eqn_of_state, EOSdom)
    endif

    do i=is,ie
      if (hb(i) + h_bbl_frac(i) < GV%H_subroundoff) h_bbl_frac(i) = GV%H_subroundoff
      rho_bot(i) = G%mask2dT(i,j) * (hb(i) + h_bbl_frac(i)) / (SpV_h_bot(i) + h_bbl_frac(i)*SpV_bbl(i))
      h_bot(i) = hb(i) + h_bbl_frac(i)
    enddo
  endif

end procedure find_rho_bottom
module procedure dz_to_thickness_tv
  integer :: i, j, k, is, ie, js, je, halo, nz
  halo = 0 ; if (present(halo_size)) halo = max(0,halo_size)
  is = G%isc-halo ; ie = G%iec+halo ; js = G%jsc-halo ; je = G%jec+halo ; nz = GV%ke

  if (GV%Boussinesq) then
    do k=1,nz ; do j=js,je ; do i=is,ie
      h(i,j,k) = GV%Z_to_H * dz(i,j,k)
    enddo ; enddo ; enddo
  else
    if (associated(tv%eqn_of_state)) then
      if (associated(tv%p_surf)) then
        call dz_to_thickness_EOS(dz, tv%T, tv%S, tv%eqn_of_state, h, G, GV, US, halo, tv%p_surf)
      else
        call dz_to_thickness_EOS(dz, tv%T, tv%S, tv%eqn_of_state, h, G, GV, US, halo)
      endif
    else
      do k=1,nz ; do j=js,je ; do i=is,ie
        h(i,j,k) = (GV%RZ_to_H * GV%Rlay(k)) * dz(i,j,k)
      enddo ; enddo ; enddo
    endif
  endif

end procedure dz_to_thickness_tv
module procedure dz_to_thickness_EOS
  real, dimension(SZI_(G),SZJ_(G)) :: &
    p_top, p_bot                  ! Pressure at the interfaces above and below a layer [R L2 T-2 ~> Pa]
  real :: dp(SZI_(G),SZJ_(G))     ! Pressure change across a layer [R L2 T-2 ~> Pa]
  real :: dz_geo(SZI_(G),SZJ_(G)) ! The change in geopotential height across a layer [L2 T-2 ~> m2 s-2]
  real :: rho(SZI_(G))            ! The in situ density [R ~> kg m-3]
  real :: dp_adj                  ! The amount by which to change the bottom pressure in an
  real :: I_gEarth                ! Unit conversion factors divided by the gravitational
  logical :: do_more(SZI_(G),SZJ_(G)) ! If true, additional iterations would be beneficial.
  logical :: do_any               ! True if there are points in this layer that need more itertions.
  integer, dimension(2) :: EOSdom ! The i-computational domain for the equation of state
  integer :: i, j, k, is, ie, js, je, halo, nz
  integer :: itt, max_itt
  halo = 0 ; if (present(halo_size)) halo = max(0,halo_size)
  is = G%isc-halo ; ie = G%iec+halo ; js = G%jsc-halo ; je = G%jec+halo ; nz = GV%ke
  max_itt = 10

  if (GV%Boussinesq) then
    do k=1,nz ; do j=js,je ; do i=is,ie
      h(i,j,k) = GV%Z_to_H * dz(i,j,k)
    enddo ; enddo ; enddo
  else
    I_gEarth = GV%RZ_to_H / GV%g_Earth

    if (present(p_surf)) then
      do j=js,je ; do i=is,ie
        p_bot(i,j) = 0.0 ; p_top(i,j) = p_surf(i,j)
      enddo ; enddo
    else
      do j=js,je ; do i=is,ie
        p_bot(i,j) = 0.0 ; p_top(i,j) = 0.0
      enddo ; enddo
    endif
    EOSdom(:) = EOS_domain(G%HI)

    ! The iterative approach here is inherited from very old code that was in the
    ! MOM_state_initialization module.  It does converge, but it is very inefficient and
    ! should be revised, although doing so would change answers in non-Boussinesq mode.
    do k=1,nz
      do j=js,je
        do i=is,ie ; p_top(i,j) = p_bot(i,j) ; enddo
        call calculate_density(Temp(:,j,k), Saln(:,j,k), p_top(:,j), rho, &
                               EoS, EOSdom)
        ! The following two expressions are mathematically equivalent.
        if (GV%semi_Boussinesq) then
          do i=is,ie
            p_bot(i,j) = p_top(i,j) + (GV%g_Earth*GV%H_to_Z) * ((GV%Z_to_H*dz(i,j,k)) * rho(i))
            dp(i,j) = (GV%g_Earth*GV%H_to_Z) * ((GV%Z_to_H*dz(i,j,k)) * rho(i))
          enddo
        else
          do i=is,ie
            p_bot(i,j) = p_top(i,j) + rho(i) * (GV%g_Earth * dz(i,j,k))
            dp(i,j) = rho(i) * (GV%g_Earth * dz(i,j,k))
          enddo
        endif
      enddo

      do_more(:,:) = .true.
      do itt=1,max_itt
        do_any = .false.
        call int_specific_vol_dp(Temp(:,:,k), Saln(:,:,k), p_top, p_bot, 0.0, G%HI, EoS, US, dz_geo)
        if (itt < max_itt) then ; do j=js,je
          call calculate_density(Temp(:,j,k), Saln(:,j,k), p_bot(:,j), rho, EoS, EOSdom)
          ! Use Newton's method to correct the bottom value.
          ! The hydrostatic equation is sufficiently linear that no bounds-checking is needed.
          if (GV%semi_Boussinesq) then
            do i=is,ie
              dp_adj = rho(i) * ((GV%g_Earth*GV%H_to_Z)*(GV%Z_to_H*dz(i,j,k)) - dz_geo(i,j))
              p_bot(i,j) = p_bot(i,j) + dp_adj
              dp(i,j) = dp(i,j) + dp_adj
            enddo
            do_any = .true. ! To avoid changing answers, always use the maximum number of itertions.
          else
            do i=is,ie ; if (do_more(i,j)) then
              dp_adj = rho(i) * (GV%g_Earth*dz(i,j,k) - dz_geo(i,j))
              p_bot(i,j) = p_bot(i,j) + dp_adj
              dp(i,j) = dp(i,j) + dp_adj
              ! Check for convergence to roundoff.
              do_more(i,j) = (abs(dp_adj) > 1.0e-15*dp(i,j))
              if (do_more(i,j)) do_any = .true.
            endif ; enddo
          endif
        enddo ; endif
        if (.not.do_any) exit
      enddo

      if (GV%semi_Boussinesq) then
        do j=js,je ; do i=is,ie
          h(i,j,k) = (p_bot(i,j) - p_top(i,j)) * I_gEarth
        enddo ; enddo
      else
        do j=js,je ; do i=is,ie
          h(i,j,k) = dp(i,j) * I_gEarth
        enddo ; enddo
      endif
    enddo
  endif

end procedure dz_to_thickness_EOS
module procedure dz_to_thickness_simple
  logical :: layered  ! If true and the model is non-Boussinesq, do calculations appropriate for use
  integer :: i, j, k, is, ie, js, je, halo, nz
  halo = 0 ; if (present(halo_size)) halo = max(0,halo_size)
  layered = .false. ; if (present(layer_mode)) layered = layer_mode
  is = G%isc-halo ; ie = G%iec+halo ; js = G%jsc-halo ; je = G%jec+halo ; nz = GV%ke

  if (GV%Boussinesq) then
    do k=1,nz ; do j=js,je ; do i=is,ie
      h(i,j,k) = GV%Z_to_H * dz(i,j,k)
    enddo ; enddo ; enddo
  elseif (layered) then
    do k=1,nz ; do j=js,je ; do i=is,ie
      h(i,j,k) = (GV%RZ_to_H * GV%Rlay(k)) * dz(i,j,k)
    enddo ; enddo ; enddo
  else
    do k=1,nz ; do j=js,je ; do i=is,ie
      h(i,j,k) = (US%Z_to_m * GV%m_to_H) * dz(i,j,k)
    enddo ; enddo ; enddo
  endif

end procedure dz_to_thickness_simple
module procedure thickness_to_dz_3d
  character(len=128) :: mesg    ! A string for error messages
  integer :: i, j, k, is, ie, js, je, halo, nz
  halo = 0 ; if (present(halo_size)) halo = max(0,halo_size)
  is = G%isc-halo ; ie = G%iec+halo ; js = G%jsc-halo ; je = G%jec+halo ; nz = GV%ke

  if ((.not.GV%Boussinesq) .and. allocated(tv%SpV_avg))  then
    if ((allocated(tv%SpV_avg)) .and. (tv%valid_SpV_halo < halo)) then
      if (tv%valid_SpV_halo < 0) then
        mesg = "invalid values of SpV_avg."
      else
        write(mesg, '("insufficiently large SpV_avg halos of width ", i2, " but ", i2," is needed.")') &
                     tv%valid_SpV_halo, halo
      endif
      call MOM_error(FATAL, "thickness_to_dz called in fully non-Boussinesq mode with "//trim(mesg))
    endif

    do k=1,nz ; do j=js,je ; do i=is,ie
      dz(i,j,k) = GV%H_to_RZ * h(i,j,k) * tv%SpV_avg(i,j,k)
    enddo ; enddo ; enddo
  else
    do k=1,nz ; do j=js,je ; do i=is,ie
      dz(i,j,k) = GV%H_to_Z * h(i,j,k)
    enddo ; enddo ; enddo
  endif

end procedure thickness_to_dz_3d
module procedure thickness_to_dz_jslice
  character(len=128) :: mesg    ! A string for error messages
  integer :: i, k, is, ie, halo, nz
  halo = 0 ; if (present(halo_size)) halo = max(0,halo_size)
  is = G%isc-halo ; ie = G%iec+halo ; nz = GV%ke

  if ((.not.GV%Boussinesq) .and. allocated(tv%SpV_avg))  then
    if ((allocated(tv%SpV_avg)) .and. (tv%valid_SpV_halo < halo)) then
      if (tv%valid_SpV_halo < 0) then
        mesg = "invalid values of SpV_avg."
      else
        write(mesg, '("insufficiently large SpV_avg halos of width ", i2, " but ", i2," is needed.")') &
                     tv%valid_SpV_halo, halo
      endif
      call MOM_error(FATAL, "thickness_to_dz called in fully non-Boussinesq mode with "//trim(mesg))
    endif

    do k=1,nz ; do i=is,ie
      dz(i,k) = GV%H_to_RZ * h(i,j,k) * tv%SpV_avg(i,j,k)
    enddo ; enddo
  else
    do k=1,nz ; do i=is,ie
      dz(i,k) = GV%H_to_Z * h(i,j,k)
    enddo ; enddo
  endif

end procedure thickness_to_dz_jslice
module procedure convert_MLD_to_ML_thickness
  real :: MLD_rem(SZI_(G)) ! The vertical extent of the MLD_in that has not yet been accounted for [Z ~> m]
  character(len=128) :: mesg    ! A string for error messages
  logical :: keep_going
  integer :: i, j, k, is, ie, js, je, nz, halos
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke

  halos = 0 ; if (present(halo)) halos = halo
  if (present(halo)) then
    is = G%isc-halo ; ie = G%iec+halo ; js = G%jsc-halo ; je = G%jec+halo
  endif

  if (GV%Boussinesq .or. (.not.allocated(tv%SpV_avg))) then
    do j=js,je ; do i=is,ie
      h_MLD(i,j) = GV%Z_to_H * MLD_in(i,j)
    enddo ; enddo
  else  ! The fully non-Boussinesq conversion between height in MLD_in and thickness.
    if ((allocated(tv%SpV_avg)) .and. (tv%valid_SpV_halo < halos)) then
      if (tv%valid_SpV_halo < 0) then
        mesg = "invalid values of SpV_avg."
      else
        write(mesg, '("insufficiently large SpV_avg halos of width ", i2, " but ", i2," is needed.")') &
                     tv%valid_SpV_halo, halos
      endif
      call MOM_error(FATAL, "convert_MLD_to_ML_thickness called in fully non-Boussinesq mode with "//trim(mesg))
    endif

    do j=js,je
      do i=is,ie ; MLD_rem(i) = MLD_in(i,j) ; h_MLD(i,j) = 0.0 ; enddo
      do k=1,nz
        keep_going = .false.
        do i=is,ie ; if (MLD_rem(i) > 0.0) then
          if (MLD_rem(i) > GV%H_to_RZ * h(i,j,k) * tv%SpV_avg(i,j,k)) then
            h_MLD(i,j) = h_MLD(i,j) + h(i,j,k)
            MLD_rem(i) = MLD_rem(i) - GV%H_to_RZ * h(i,j,k) * tv%SpV_avg(i,j,k)
            keep_going = .true.
          else
            h_MLD(i,j) = h_MLD(i,j) + GV%RZ_to_H * MLD_rem(i) / tv%SpV_avg(i,j,k)
            MLD_rem(i) = 0.0
          endif
        endif ; enddo
        if (.not.keep_going) exit
      enddo
    enddo
  endif

end procedure convert_MLD_to_ML_thickness
end submodule MOM_interface_heights_s
