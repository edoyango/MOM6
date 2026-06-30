submodule (MOM_sponge) MOM_sponge_s
#include <MOM_memory.h>
  implicit none
contains
module procedure initialize_sponge
# include "version_variable.h"
  character(len=40)  :: mdl = "MOM_sponge"  ! This module's name.
  logical :: use_sponge
  integer :: i, j, k, col, total_sponge_cols
  if (associated(CS)) then
    call MOM_error(WARNING, "initialize_sponge called with an associated "// &
                            "control structure.")
    return
  endif

! Set default, read and log parameters
  call log_version(param_file, mdl, version)
  call get_param(param_file, mdl, "SPONGE", use_sponge, &
                 "If true, sponges may be applied anywhere in the domain. "//&
                 "The exact location and properties of those sponges are "//&
                 "specified from MOM_initialization.F90.", default=.false.)

  if (.not.use_sponge) return
  allocate(CS)

  if (present(Iresttime_i_mean) .neqv. present(int_height_i_mean)) &
    call MOM_error(FATAL, "initialize_sponge:  The optional arguments \n"//&
           "Iresttime_i_mean and int_height_i_mean must both be present \n"//&
           "if either one is.")

  CS%do_i_mean_sponge = present(Iresttime_i_mean)

  CS%nz = GV%ke

  ! CS%bulkmixedlayer may be set later via a call to set_up_sponge_ML_density.
  CS%bulkmixedlayer = .false.

  CS%num_col = 0 ; CS%fldno = 0
  do j=G%jsc,G%jec ; do i=G%isc,G%iec
    if ((Iresttime(i,j) > 0.0) .and. (G%mask2dT(i,j) > 0.0)) &
      CS%num_col = CS%num_col + 1
  enddo ; enddo

  if (CS%num_col > 0) then

    allocate(CS%Iresttime_col(CS%num_col), source=0.0)
    allocate(CS%col_i(CS%num_col), source=0)
    allocate(CS%col_j(CS%num_col), source=0)

    col = 1
    do j=G%jsc,G%jec ; do i=G%isc,G%iec
      if ((Iresttime(i,j) > 0.0) .and. (G%mask2dT(i,j) > 0.0)) then
        CS%col_i(col) = i ; CS%col_j(col) = j
        CS%Iresttime_col(col) = Iresttime(i,j)
        col = col +1
      endif
    enddo ; enddo

    allocate(CS%Ref_eta(CS%nz+1,CS%num_col))
    do col=1,CS%num_col ; do K=1,CS%nz+1
      CS%Ref_eta(K,col) = int_height(CS%col_i(col),CS%col_j(col),K)
    enddo ; enddo

  endif

  if (CS%do_i_mean_sponge) then
    allocate(CS%Iresttime_im(G%jsd:G%jed), source=0.0)
    allocate(CS%Ref_eta_im(G%jsd:G%jed,GV%ke+1), source=0.0)

    do j=G%jsc,G%jec
      CS%Iresttime_im(j) = Iresttime_i_mean(j)
    enddo
    do K=1,CS%nz+1 ; do j=G%jsc,G%jec
      CS%Ref_eta_im(j,K) = int_height_i_mean(j,K)
    enddo ; enddo
  endif

  total_sponge_cols = CS%num_col
  call sum_across_PEs(total_sponge_cols)

  call log_param(param_file, mdl, "!Total sponge columns", total_sponge_cols, &
                 "The total number of columns where sponges are applied.")

end procedure initialize_sponge
module procedure init_sponge_diags
  if (.not.associated(CS)) return

  CS%diag => diag
  CS%id_w_sponge = register_diag_field('ocean_model', 'w_sponge', diag%axesTi, &
      Time, 'The diapycnal motion due to the sponges', 'm s-1', conversion=GV%H_to_m*US%s_to_T)

end procedure init_sponge_diags
module procedure set_up_sponge_field
  integer :: j, k, col
  character(len=256) :: mesg ! String for error messages
  if (.not.associated(CS)) return

  CS%fldno = CS%fldno + 1

  if (CS%fldno > MAX_FIELDS_) then
    write(mesg,'("Increase MAX_FIELDS_ to at least ",I0," in MOM_memory.h or decrease &
           &the number of fields to be damped in the call to &
           &initialize_sponge." )') CS%fldno
    call MOM_error(FATAL,"set_up_sponge_field: "//mesg)
  endif

  allocate(CS%Ref_val(CS%fldno)%p(CS%nz,CS%num_col), source=0.0)
  do col=1,CS%num_col
    do k=1,nlay
      CS%Ref_val(CS%fldno)%p(k,col) = sp_val(CS%col_i(col),CS%col_j(col),k)
    enddo
    do k=nlay+1,CS%nz
      CS%Ref_val(CS%fldno)%p(k,col) = 0.0
    enddo
  enddo

  CS%var(CS%fldno)%p => f_ptr

  if (nlay/=CS%nz) then
    write(mesg,'("Danger: Sponge reference fields require nz (",I0,") layers.&
        & A field with ",I0," layers was passed to set_up_sponge_field.")') &
          CS%nz, nlay
    if (is_root_pe()) call MOM_error(WARNING, "set_up_sponge_field: "//mesg)
  endif

  if (CS%do_i_mean_sponge) then
    if (.not.present(sp_val_i_mean)) call MOM_error(FATAL, &
      "set_up_sponge_field: sp_val_i_mean must be present with i-mean sponges.")

    allocate(CS%Ref_val_im(CS%fldno)%p(G%jsd:G%jed,CS%nz), source=0.0)
    do k=1,CS%nz ; do j=G%jsc,G%jec
      CS%Ref_val_im(CS%fldno)%p(j,k) = sp_val_i_mean(j,k)
    enddo ; enddo
  endif

end procedure set_up_sponge_field
module procedure set_up_sponge_ML_density
  integer :: j, col
  if (.not.associated(CS)) return

  if (associated(CS%Rcv_ml_ref)) then
    call MOM_error(FATAL, "set_up_sponge_ML_density appears to have been "//&
                           "called twice.")
  endif

  CS%bulkmixedlayer = .true.
  allocate(CS%Rcv_ml_ref(CS%num_col), source=0.0)
  do col=1,CS%num_col
    CS%Rcv_ml_ref(col) = sp_val(CS%col_i(col),CS%col_j(col))
  enddo

  if (CS%do_i_mean_sponge) then
    if (.not.present(sp_val_i_mean)) call MOM_error(FATAL, &
      "set_up_sponge_field: sp_val_i_mean must be present with i-mean sponges.")

    allocate(CS%Rcv_ml_ref_im(G%jsd:G%jed), source=0.0)
    do j=G%jsc,G%jec
      CS%Rcv_ml_ref_im(j) = sp_val_i_mean(j)
    enddo
  endif

end procedure set_up_sponge_ML_density
module procedure apply_sponge
  real, dimension(SZI_(G), SZJ_(G), SZK_(GV)+1) :: &
    w_int, &       ! Water moved upward across an interface within a timestep,
                   ! [H ~> m or kg m-2].
    e_D            ! Interface heights that are dilated to have a value of 0
                   ! at the surface [Z ~> m].
  real, dimension(SZI_(G), SZJ_(G)) :: &
    eta_anom, &    ! Anomalies in the interface height, relative to the i-mean
                   ! target value [Z ~> m].
    fld_anom       ! Anomalies in a tracer concentration, relative to the
                   ! i-mean target value [various]
  real, dimension(SZJ_(G), SZK_(GV)+1) :: &
    eta_mean_anom  ! The i-mean interface height anomalies [Z ~> m].
  real, allocatable, dimension(:,:,:) :: &
    fld_mean_anom  ! The i-mean tracer concentration anomalies [various]
  real, dimension(SZI_(G), SZK_(GV)+1) :: &
    h_above, &     ! The total thickness above an interface [H ~> m or kg m-2].
    h_below        ! The total thickness below an interface [H ~> m or kg m-2].
  real, dimension(SZI_(G)) :: &
    dilate         ! A nondimensional factor by which to dilate layers to
                   ! give 0 at the surface [nondim].

  real :: e(SZK_(GV)+1)  ! The interface heights [Z ~> m], usually negative.
  real :: dz_to_h(SZK_(GV)+1)  ! Factors used to convert interface height movement
                   ! to thickness fluxes [H Z-1 ~> nondim or kg m-3]
  real :: e0       ! The height of the free surface [Z ~> m].
  real :: e_str    ! A nondimensional amount by which the reference
                   ! profile must be stretched for the free surfaces
                   ! heights in the two profiles to agree [nondim].
  real :: w_mean   ! The vertical displacement of water moving upward through an
                   ! interface within 1 timestep [Z ~> m].
  real :: w        ! The thickness of water moving upward through an
                   ! interface within 1 timestep [H ~> m or kg m-2].
  real :: wm       ! wm is w if w is negative and 0 otherwise [H ~> m or kg m-2].
  real :: wb       ! w at the interface below a layer [H ~> m or kg m-2].
  real :: wpb      ! wpb is wb if wb is positive and 0 otherwise [H ~> m or kg m-2].
  real :: ea_k     ! Water entrained from above within a timestep [H ~> m or kg m-2]
  real :: eb_k     ! Water entrained from below within a timestep [H ~> m or kg m-2]
  real :: damp     ! The timestep times the local damping coefficient [nondim].
  real :: I1pdamp  ! I1pdamp is 1/(1 + damp). [nondim]
  real :: damp_1pdamp ! damp_1pdamp is damp/(1 + damp). [nondim]
  real :: Idt      ! The inverse of the timestep [T-1 ~> s-1]
  integer :: c, m, nkmb, i, j, k, is, ie, js, je, nz
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke

  if (.not.associated(CS)) return
  if (CS%bulkmixedlayer) nkmb = GV%nk_rho_varies
  if (CS%bulkmixedlayer .and. (.not.present(Rcv_ml))) &
    call MOM_error(FATAL, "Rml must be provided to apply_sponge when using "//&
                           "a bulk mixed layer.")

  if ((CS%id_w_sponge > 0) .or. CS%do_i_mean_sponge) then
    do k=1,nz+1 ; do j=js,je ; do i=is,ie
      w_int(i,j,K) = 0.0
    enddo ; enddo ; enddo
  endif

  if (CS%do_i_mean_sponge) then
    ! Apply forcing to restore the zonal-mean properties to prescribed values.

    if (CS%bulkmixedlayer) call MOM_error(FATAL, "apply_sponge is not yet set up to "//&
                  "work properly with i-mean sponges and a bulk mixed layer.")

    do j=js,je ; do i=is,ie ; e_D(i,j,nz+1) = -G%bathyT(i,j) ; enddo ; enddo
    if ((.not.GV%Boussinesq) .and. allocated(tv%SpV_avg)) then
      do k=nz,1,-1 ; do j=js,je ; do i=is,ie
        e_D(i,j,K) = e_D(i,j,K+1) + GV%H_to_RZ * h(i,j,k) * tv%SpV_avg(i,j,k)
      enddo ; enddo ; enddo
    else
      do k=nz,1,-1 ; do j=js,je ; do i=is,ie
        e_D(i,j,K) = e_D(i,j,K+1) + h(i,j,k)*GV%H_to_Z
      enddo ; enddo ; enddo
    endif
    do j=js,je
      do i=is,ie
        dilate(i) = (G%bathyT(i,j) + G%Z_ref) / (e_D(i,j,1) + G%bathyT(i,j))
      enddo
      do k=1,nz+1 ; do i=is,ie
        e_D(i,j,K) = dilate(i) * (e_D(i,j,K) + G%bathyT(i,j)) - (G%bathyT(i,j) + G%Z_ref)
      enddo ; enddo
    enddo

    do k=2,nz
      do j=js,je ; do i=is,ie
        eta_anom(i,j) = e_D(i,j,k) - CS%Ref_eta_im(j,k)
        if (CS%Ref_eta_im(j,K) < -(G%bathyT(i,j) + G%Z_ref)) eta_anom(i,j) = 0.0
      enddo ; enddo
      call global_i_mean(eta_anom(:,:), eta_mean_anom(:,K), G, tmp_scale=US%Z_to_m)
    enddo

    if (CS%fldno > 0) allocate(fld_mean_anom(G%isd:G%ied,nz,CS%fldno))
    do m=1,CS%fldno
      do j=js,je ; do i=is,ie
        fld_anom(i,j) = CS%var(m)%p(i,j,k) - CS%Ref_val_im(m)%p(j,k)
      enddo ; enddo
      call global_i_mean(fld_anom(:,:), fld_mean_anom(:,k,m), G, h(:,:,k))
    enddo

    do j=js,je ; if (CS%Iresttime_im(j) > 0.0) then
      damp = dt * CS%Iresttime_im(j) ; damp_1pdamp = damp / (1.0 + damp)

      do i=is,ie
        h_above(i,1) = 0.0 ; h_below(i,nz+1) = 0.0
      enddo
      do K=nz,1,-1 ; do i=is,ie
        h_below(i,K) = h_below(i,K+1) + max(h(i,j,k)-GV%Angstrom_H, 0.0)
      enddo ; enddo
      do K=2,nz+1 ; do i=is,ie
        h_above(i,K) = h_above(i,K-1) + max(h(i,j,k-1)-GV%Angstrom_H, 0.0)
      enddo ; enddo

      ! In both blocks below, w is positive for an upward (lightward) flux of mass,
      ! resulting in the downward movement of an interface.
      if ((.not.GV%Boussinesq) .and. allocated(tv%SpV_avg)) then
        do K=2,nz
          w_mean = damp_1pdamp * eta_mean_anom(j,K)
          do i=is,ie
            w = w_mean * 2.0*GV%RZ_to_H / (tv%SpV_avg(i,j,k-1) + tv%SpV_avg(i,j,k))
            if (w > 0.0) then
              w_int(i,j,K) = min(w, h_below(i,K))
              eb(i,j,k-1) = eb(i,j,k-1) + w_int(i,j,K)
            else
              w_int(i,j,K) = max(w, -h_above(i,K))
              ea(i,j,k) = ea(i,j,k) - w_int(i,j,K)
            endif
          enddo
        enddo
      else
        do K=2,nz
          w = damp_1pdamp * eta_mean_anom(j,K) * GV%Z_to_H
          if (w > 0.0) then
            do i=is,ie
              w_int(i,j,K) = min(w, h_below(i,K))
              eb(i,j,k-1) = eb(i,j,k-1) + w_int(i,j,K)
            enddo
          else
            do i=is,ie
              w_int(i,j,K) = max(w, -h_above(i,K))
              ea(i,j,k) = ea(i,j,k) - w_int(i,j,K)
            enddo
          endif
        enddo
      endif
      do k=1,nz ; do i=is,ie
        ea_k = max(0.0, -w_int(i,j,K))
        eb_k = max(0.0, w_int(i,j,K+1))
        do m=1,CS%fldno
          CS%var(m)%p(i,j,k) = (h(i,j,k)*CS%var(m)%p(i,j,k) + &
              CS%Ref_val_im(m)%p(j,k) * (ea_k + eb_k)) / &
                     (h(i,j,k) + (ea_k + eb_k)) - &
              damp_1pdamp * fld_mean_anom(j,k,m)
        enddo

        h(i,j,k) = max(h(i,j,k) + (w_int(i,j,K+1) - w_int(i,j,K)), &
                       min(h(i,j,k), GV%Angstrom_H))
      enddo ; enddo
    endif ; enddo

    if (CS%fldno > 0) deallocate(fld_mean_anom)

  endif

  do c=1,CS%num_col
    i = CS%col_i(c) ; j = CS%col_j(c)
    damp = dt * CS%Iresttime_col(c)

    e(1) = 0.0 ; e0 = 0.0
    if ((.not.GV%Boussinesq) .and. allocated(tv%SpV_avg)) then
      do K=1,nz
        e(K+1) = e(K) - GV%H_to_RZ * h(i,j,k) * tv%SpV_avg(i,j,k)
      enddo
      dz_to_h(1) = GV%RZ_to_H / tv%SpV_avg(i,j,1)
      do K=2,nz
        dz_to_h(K) = 2.0*GV%RZ_to_H / (tv%SpV_avg(i,j,k-1) + tv%SpV_avg(i,j,k))
      enddo
    else
      do K=1,nz
        e(K+1) = e(K) - h(i,j,k)*GV%H_to_Z
        dz_to_h(K) = GV%Z_to_H
      enddo
    endif
    e_str = e(nz+1) / CS%Ref_eta(nz+1,c)

    if ( CS%bulkmixedlayer ) then
      I1pdamp = 1.0 / (1.0 + damp)
      if (associated(CS%Rcv_ml_ref)) &
        Rcv_ml(i,j) = I1pdamp * (Rcv_ml(i,j) + CS%Rcv_ml_ref(c)*damp)
      do k=1,nkmb
        do m=1,CS%fldno
          CS%var(m)%p(i,j,k) = I1pdamp * &
              (CS%var(m)%p(i,j,k) + CS%Ref_val(m)%p(k,c)*damp)
        enddo
      enddo

      wpb = 0.0 ; wb = 0.0
      do k=nz,nkmb+1,-1
        if (GV%Rlay(k) > Rcv_ml(i,j)) then
          w = MIN((((e(K)-e0) - e_str*CS%Ref_eta(K,c)) * damp)*dz_to_h(K), &
                    ((wb + h(i,j,k)) - GV%Angstrom_H))
          wm = 0.5*(w-ABS(w))
          do m=1,CS%fldno
            CS%var(m)%p(i,j,k) = (h(i,j,k)*CS%var(m)%p(i,j,k) + &
                     CS%Ref_val(m)%p(k,c)*(damp*h(i,j,k) + (wpb - wm))) / &
                     (h(i,j,k)*(1.0 + damp) + (wpb - wm))
          enddo
        else
          do m=1,CS%fldno
            CS%var(m)%p(i,j,k) = I1pdamp * &
              (CS%var(m)%p(i,j,k) + CS%Ref_val(m)%p(k,c)*damp)
          enddo
          w = wb + (h(i,j,k) - GV%Angstrom_H)
          wm = 0.5*(w-ABS(w))
        endif
        eb(i,j,k) = eb(i,j,k) + wpb
        ea(i,j,k) = ea(i,j,k) - wm
        h(i,j,k)  = h(i,j,k)  + (wb - w)
        wb = w
        wpb = w - wm
      enddo

      if (wb < 0) then
        do k=nkmb,1,-1
          w = MIN((wb + (h(i,j,k) - GV%Angstrom_H)),0.0)
          h(i,j,k)  = h(i,j,k)  + (wb - w)
          ea(i,j,k) = ea(i,j,k) - w
          wb = w
        enddo
      else
        w = wb
        do k=GV%nkml,nkmb
          eb(i,j,k) = eb(i,j,k) + w
        enddo

        k = GV%nkml
        h(i,j,k) = h(i,j,k) + w
        do m=1,CS%fldno
          CS%var(m)%p(i,j,k) = (CS%var(m)%p(i,j,k)*h(i,j,k) + &
                                CS%Ref_val(m)%p(k,c)*w) / (h(i,j,k) + w)
        enddo
      endif

      do k=1,nkmb
        do m=1,CS%fldno
          CS%var(m)%p(i,j,k) = I1pdamp * &
              (CS%var(m)%p(i,j,k) + CS%Ref_val(m)%p(GV%nkml,c)*damp)
        enddo
      enddo

    else                                          ! not BULKMIXEDLAYER

      wpb = 0.0
      wb = 0.0
      do k=nz,1,-1
        w = MIN((((e(K)-e0) - e_str*CS%Ref_eta(K,c)) * damp)*dz_to_h(K), &
                  ((wb + h(i,j,k)) - GV%Angstrom_H))
        wm = 0.5*(w - ABS(w))
        do m=1,CS%fldno
          CS%var(m)%p(i,j,k) = (h(i,j,k)*CS%var(m)%p(i,j,k) + &
              CS%Ref_val(m)%p(k,c) * (damp*h(i,j,k) + (wpb - wm))) / &
                     (h(i,j,k)*(1.0 + damp) + (wpb - wm))
        enddo
        eb(i,j,k) = eb(i,j,k) + wpb
        ea(i,j,k) = ea(i,j,k) - wm
        h(i,j,k)  = h(i,j,k)  + (wb - w)
        wb = w
        wpb = w - wm
      enddo

    endif                                         ! end BULKMIXEDLAYER
  enddo ! end of c loop

  if (associated(CS%diag)) then ; if (query_averaging_enabled(CS%diag)) then
    if (CS%id_w_sponge > 0) then
      Idt = 1.0 / dt
      do k=1,nz+1 ; do j=js,je ; do i=is,ie
        w_int(i,j,K) = w_int(i,j,K) * Idt ! Scale values by clobbering array since it is local
      enddo ; enddo ; enddo
      call post_data(CS%id_w_sponge, w_int(:,:,:), CS%diag)
    endif
  endif ; endif

end procedure apply_sponge
module procedure sponge_end
  integer :: m
  if (.not.associated(CS)) return

  if (associated(CS%col_i)) deallocate(CS%col_i)
  if (associated(CS%col_j)) deallocate(CS%col_j)

  if (associated(CS%Iresttime_col)) deallocate(CS%Iresttime_col)
  if (associated(CS%Rcv_ml_ref)) deallocate(CS%Rcv_ml_ref)
  if (associated(CS%Ref_eta)) deallocate(CS%Ref_eta)

  if (associated(CS%Iresttime_im)) deallocate(CS%Iresttime_im)
  if (associated(CS%Rcv_ml_ref_im)) deallocate(CS%Rcv_ml_ref_im)
  if (associated(CS%Ref_eta_im)) deallocate(CS%Ref_eta_im)

  do m=1,CS%fldno
    if (associated(CS%Ref_val(CS%fldno)%p)) deallocate(CS%Ref_val(CS%fldno)%p)
    if (associated(CS%Ref_val_im(CS%fldno)%p)) &
      deallocate(CS%Ref_val_im(CS%fldno)%p)
  enddo

  deallocate(CS)

end procedure sponge_end
end submodule MOM_sponge_s
