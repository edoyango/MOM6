submodule (pseudo_salt_tracer) pseudo_salt_tracer_s
#include <MOM_memory.h>
  implicit none
contains
module procedure register_pseudo_salt_tracer
  character(len=40)  :: mdl = "pseudo_salt_tracer" ! This module's name.
  character(len=48)  :: var_name ! The variable's name.
# include "version_variable.h"
  real, pointer :: tr_ptr(:,:,:) => NULL() ! The tracer concentration [ppt]
  integer :: isd, ied, jsd, jed, nz
  isd = HI%isd ; ied = HI%ied ; jsd = HI%jsd ; jed = HI%jed ; nz = GV%ke

  if (associated(CS)) then
    call MOM_error(FATAL, "register_pseudo_salt_tracer called with an "// &
                          "associated control structure.")
  endif
  allocate(CS)

  ! Read all relevant parameters and write them to the model log.
  call log_version(param_file, mdl, version, "")

  allocate(CS%ps(isd:ied,jsd:jed,nz), source=0.0)

  CS%tr_desc = var_desc(trim("pseudo_salt"), "psu", &
                     "Pseudo salt passive tracer", caller=mdl)

  tr_ptr => CS%ps(:,:,:)
  call query_vardesc(CS%tr_desc, name=var_name, caller="register_pseudo_salt_tracer")
  ! Register the tracer for horizontal advection, diffusion, and restarts.
  call register_tracer(tr_ptr, tr_Reg, param_file, HI, GV, name="pseudo_salt", &
                       longname="Pseudo salt passive tracer", units="psu", &
                       registry_diags=.true., restart_CS=restart_CS, &
                       mandatory=.not.CS%pseudo_salt_may_reinit, Tr_out=CS%tr_ptr)

  CS%tr_Reg => tr_Reg
  CS%restart_CSp => restart_CS
  register_pseudo_salt_tracer = .true.

end procedure register_pseudo_salt_tracer
module procedure initialize_pseudo_salt_tracer
  character(len=16) :: name     ! A variable's name in a NetCDF file
  integer :: i, j, k, isd, ied, jsd, jed, nz
  if (.not.associated(CS)) return

  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed ; nz = GV%ke

  CS%Time => day
  CS%diag => diag
  name = "pseudo_salt"

  call query_vardesc(CS%tr_desc, name=name, caller="initialize_pseudo_salt_tracer")
  if ((.not.restart) .or. (.not.query_initialized(CS%ps, name, CS%restart_CSp))) then
    do k=1,nz ; do j=jsd,jed ; do i=isd,ied
      CS%ps(i,j,k) = US%S_to_ppt*tv%S(i,j,k)
    enddo ; enddo ; enddo
    call set_initialized(CS%ps, name, CS%restart_CSp)
  endif

  if (associated(OBC)) then
  ! Steal from updated DOME in the fullness of time.
  endif

  CS%id_psd = register_diag_field("ocean_model", "pseudo_salt_diff", CS%diag%axesTL, &
        day, "Difference between pseudo salt passive tracer and salt tracer", "psu")
  if (.not.allocated(CS%diff)) allocate(CS%diff(isd:ied,jsd:jed,nz), source=0.0)

end procedure initialize_pseudo_salt_tracer
module procedure pseudo_salt_tracer_column_physics
  real :: net_salt_rate(SZI_(G),SZJ_(G)) ! Net salt flux into the ocean
  real :: net_salt(SZI_(G),SZJ_(G)) ! Net salt flux into the ocean integrated over
  real :: htot(SZI_(G))       ! Total ocean depth [H ~> m or kg m-2]
  real :: FluxRescaleDepth    ! Minimum total ocean depth at which fluxes start to be scaled
  real :: Ih_limit            ! Inverse of FluxRescaleDepth or 0 for no limiting [H-1 ~> m-1 or m2 kg-1]
  real :: scale               ! Scale scales away fluxes if depth < FluxRescaleDepth [nondim]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)) :: h_work ! Used so that h can be modified [H ~> m or kg m-2]
  integer :: i, j, k, is, ie, js, je, nz
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke

  if (.not.associated(CS)) return
  if (.not.associated(CS%ps)) return

  if (debug) then
    call hchksum(tv%S,"salt pre pseudo-salt vertdiff", G%HI, unscale=US%S_to_ppt)
    call hchksum(CS%ps,"pseudo_salt pre pseudo-salt vertdiff", G%HI)
  endif

  FluxRescaleDepth = max( GV%Angstrom_H, 1.e-30*GV%m_to_H )
  Ih_limit  = 0.0 ; if (FluxRescaleDepth > 0.0) Ih_limit  = 1.0 / FluxRescaleDepth

  ! Compute KPP nonlocal term if necessary
  if (present(KPP_CSp)) then
    if (associated(KPP_CSp) .and. present(nonLocalTrans)) then
      ! Determine the salt flux, including limiting for small total ocean depths.
      net_salt_rate(:,:) = 0.0
      if (associated(fluxes%salt_flux)) then
        do j=js,je
          do i=is,ie ; htot(i) = h_old(i,j,1) ; enddo
          do k=2,nz ; do i=is,ie ; htot(i) = htot(i) + h_old(i,j,k) ; enddo ; enddo
          do i=is,ie
            scale = 1.0 ; if ((Ih_limit > 0.0) .and. (htot(i)*Ih_limit < 1.0)) scale = htot(i)*Ih_limit
            net_salt_rate(i,j) = (scale * (1000.0 * fluxes%salt_flux(i,j))) * GV%RZ_to_H
          enddo
        enddo
      endif
      call KPP_NonLocalTransport(KPP_CSp, G, GV, h_old, nonLocalTrans, net_salt_rate, &
                                 dt, CS%diag, CS%tr_ptr, CS%ps(:,:,:))
    endif
  endif

  ! This uses applyTracerBoundaryFluxesInOut, usually in ALE mode
  if (present(evap_CFL_limit) .and. present(minimum_forcing_depth)) then
    ! This option uses applyTracerBoundaryFluxesInOut, usually in ALE mode

    ! Determine the time-integrated salt flux, including limiting for small total ocean depths.
    net_Salt(:,:) = 0.0
    do j=js,je
      do i=is,ie ; htot(i) = h_old(i,j,1) ; enddo
      do k=2,nz ; do i=is,ie ; htot(i) = htot(i) + h_old(i,j,k) ; enddo ; enddo
      do i=is,ie
        scale = 1.0 ; if ((Ih_limit > 0.0) .and. (htot(i)*Ih_limit < 1.0)) scale = htot(i)*Ih_limit
        net_salt(i,j) = (scale * dt * (1000.0 * fluxes%salt_flux(i,j))) * GV%RZ_to_H
      enddo
    enddo

    do k=1,nz ; do j=js,je ; do i=is,ie
      h_work(i,j,k) = h_old(i,j,k)
    enddo ; enddo ; enddo
    call applyTracerBoundaryFluxesInOut(G, GV, CS%ps, dt, fluxes, h_work, evap_CFL_limit, &
                                        minimum_forcing_depth, out_flux_optional=net_salt)
    call tracer_vertdiff(h_work, ea, eb, dt, CS%ps, G, GV)
  else
    call tracer_vertdiff(h_old, ea, eb, dt, CS%ps, G, GV)
  endif

  if (debug) then
    call hchksum(tv%S, "salt post pseudo-salt vertdiff", G%HI, unscale=US%S_to_ppt)
    call hchksum(CS%ps, "pseudo_salt post pseudo-salt vertdiff", G%HI)
  endif

  if (allocated(CS%diff)) then
    do k=1,nz ; do j=js,je ; do i=is,ie
      CS%diff(i,j,k) = CS%ps(i,j,k) - US%S_to_ppt*tv%S(i,j,k)
    enddo ; enddo ; enddo
    if (CS%id_psd>0) call post_data(CS%id_psd, CS%diff, CS%diag)
  endif

end procedure pseudo_salt_tracer_column_physics
module procedure pseudo_salt_stock
  pseudo_salt_stock = 0
  if (.not.associated(CS)) return
  if (.not.allocated(CS%diff)) return

  if (present(stock_index)) then ; if (stock_index > 0) then
    ! Check whether this stock is available from this routine.

    ! No stocks from this routine are being checked yet.  Return 0.
    return
  endif ; endif

  call query_vardesc(CS%tr_desc, name=names(1), units=units(1), caller="pseudo_salt_stock")
  units(1) = trim(units(1))//" kg"
  stocks(1) = global_mass_int_EFP(h, G, GV, CS%diff, on_PE_only=.true.)

  pseudo_salt_stock = 1

end procedure pseudo_salt_stock
module procedure pseudo_salt_tracer_surface_state
  if (.not.associated(CS)) return

  ! By design, this tracer package does not return any surface states.

end procedure pseudo_salt_tracer_surface_state
module procedure pseudo_salt_tracer_end
  if (associated(CS)) then
    if (associated(CS%ps)) deallocate(CS%ps)
    if (allocated(CS%diff)) deallocate(CS%diff)
    deallocate(CS)
  endif
end procedure pseudo_salt_tracer_end
end submodule pseudo_salt_tracer_s
