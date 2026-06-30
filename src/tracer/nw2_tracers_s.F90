submodule (nw2_tracers) nw2_tracers_s
#include <MOM_memory.h>
  implicit none
contains
module procedure register_nw2_tracers
# include "version_variable.h"
  character(len=40)  :: mdl = "nw2_tracers" ! This module's name.
  character(len=8)  :: var_name ! The variable's name.
  real, pointer :: tr_ptr(:,:,:) => NULL() ! The tracer concentration [conc]
  integer :: isd, ied, jsd, jed, nz, m, ig
  integer :: n_groups ! Number of groups of three tracers (i.e. # tracers/3)
  real, allocatable, dimension(:) :: timescale_in_days ! Damping timescale [days]
  type(vardesc) :: tr_desc ! Descriptions and metadata for the tracers
  isd = HI%isd ; ied = HI%ied ; jsd = HI%jsd ; jed = HI%jed ; nz = GV%ke

  if (associated(CS)) then
    call MOM_error(FATAL, "register_nw2_tracer called with an "// &
                          "associated control structure.")
  endif
  allocate(CS)

  ! Read all relevant parameters and write them to the model log.
  call log_version(param_file, mdl, version, "")
  call get_param(param_file, mdl, "NW2_TRACER_GROUPS", n_groups, &
                 "The number of tracer groups where a group is of three tracers "//&
                 "initialized and restored to sin(x), y and z, respectively. Each "//&
                 "group is restored with an independent restoration rate.", &
                 default=3)
  allocate(timescale_in_days(n_groups))
  timescale_in_days = (/365., 730., 1460./)
  call get_param(param_file, mdl, "NW2_TRACER_RESTORE_TIMESCALE", timescale_in_days, &
                 "A list of timescales, one for each tracer group.", &
                 units="days")

  CS%ntr = 3 * n_groups
  allocate(CS%tr(isd:ied,jsd:jed,nz,CS%ntr), source=0.0)
  allocate(CS%restore_rate(CS%ntr))

  do m=1,CS%ntr
    write(var_name(1:8),'(a6,i2.2)') 'tracer',m
    tr_desc = var_desc(var_name, "1", "Ideal Tracer", caller=mdl)
    ! This is needed to force the compiler not to do a copy in the registration
    ! calls.  Curses on the designers and implementers of Fortran90.
    tr_ptr => CS%tr(:,:,:,m)
    ! Register the tracer for horizontal advection, diffusion, and restarts.
    call register_tracer(tr_ptr, tr_Reg, param_file, HI, GV, tr_desc=tr_desc, &
                         registry_diags=.true., restart_CS=restart_CS, mandatory=.false.)
    ig = int( (m+2)/3 ) ! maps (1,2,3)->1, (4,5,6)->2, ...
    CS%restore_rate(m) = 1.0 / ( timescale_in_days(ig) * 86400.0*US%s_to_T )
  enddo

  CS%tr_Reg => tr_Reg
  CS%restart_CSp => restart_CS
  register_nw2_tracers = .true.
end procedure register_nw2_tracers
module procedure initialize_nw2_tracers
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)+1) :: eta ! Interface heights [Z ~> m]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)) :: dz  ! Vertical extent of layers [Z ~> m]
  real :: rscl ! z* scaling factor [nondim]
  character(len=8)  :: var_name ! The variable's name.
  integer :: i, j, k, m
  if (.not.associated(CS)) return

  CS%Time => day
  CS%diag => diag

  ! Calculate z* interface positions
  call thickness_to_dz(h, tv, dz, G, GV, US)

  if (GV%Boussinesq) then
    ! First calculate interface positions in z-space (m)
    do j=G%jsc,G%jec ; do i=G%isc,G%iec
      eta(i,j,GV%ke+1) = - G%mask2dT(i,j) * G%bathyT(i,j)
    enddo ; enddo
    do k=GV%ke,1,-1 ; do j=G%jsc,G%jec ; do i=G%isc,G%iec
      eta(i,j,K) = eta(i,j,K+1) + G%mask2dT(i,j) * dz(i,j,k)
    enddo ; enddo ; enddo
    ! Re-calculate for interface positions in z*-space (m)
    do j=G%jsc,G%jec ; do i=G%isc,G%iec
      if (G%bathyT(i,j)>0.) then
        rscl = G%bathyT(i,j) / ( eta(i,j,1) + G%bathyT(i,j) )
        do K=GV%ke, 1, -1
          eta(i,j,K) = eta(i,j,K+1) + G%mask2dT(i,j) * dz(i,j,k) * rscl
        enddo
      endif
    enddo ; enddo
  else
    call MOM_error(FATAL, "NW2 tracers assume Boussinesq mode")
  endif

  do m=1,CS%ntr
    ! Initialize only if this is not a restart or we are using a restart
    ! in which the tracers were not present
    write(var_name(1:8),'(a6,i2.2)') 'tracer',m
    if ((.not.restart) .or. &
        (.not. query_initialized(CS%tr(:,:,:,m), var_name, CS%restart_CSp))) then
      do k=1,GV%ke ; do j=G%jsc,G%jec ; do i=G%isc,G%iec
          CS%tr(i,j,k,m) = nw2_tracer_dist(m, G, GV, eta, i, j, k)
      enddo ; enddo ; enddo
      call set_initialized(CS%tr(:,:,:,m), var_name, CS%restart_CSp)
    endif ! restart
  enddo ! Tracer loop

end procedure initialize_nw2_tracers
module procedure nw2_tracer_column_physics
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)) :: h_work ! Used so that h can be modified [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)+1) :: eta ! Interface heights [Z ~> m]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)) :: dz  ! Vertical extent of layers [Z ~> m]
  integer :: i, j, k, m
  real :: dt_x_rate ! dt * restoring rate [nondim]
  real :: rscl ! z* scaling factor [nondim]
  real :: target_value ! tracer target value for damping [conc]
  if (present(evap_CFL_limit) .and. present(minimum_forcing_depth)) then
    do m=1,CS%ntr
      do k=1,GV%ke ; do j=G%jsc,G%jec ; do i=G%isc,G%iec
        h_work(i,j,k) = h_old(i,j,k)
      enddo ; enddo ; enddo
      call applyTracerBoundaryFluxesInOut(G, GV, CS%tr(:,:,:,m), dt, fluxes, h_work, &
                                          evap_CFL_limit, minimum_forcing_depth)
      call tracer_vertdiff(h_work, ea, eb, dt, CS%tr(:,:,:,m), G, GV)
    enddo
  else
    do m=1,CS%ntr
      call tracer_vertdiff(h_old, ea, eb, dt, CS%tr(:,:,:,m), G, GV)
    enddo
  endif

  ! Calculate z* interface positions
  call thickness_to_dz(h_new, tv, dz, G, GV, US)

  if (GV%Boussinesq) then
    ! First calculate interface positions in z-space [Z ~> m]
    do j=G%jsc,G%jec ; do i=G%isc,G%iec
      eta(i,j,GV%ke+1) = - G%mask2dT(i,j) * G%bathyT(i,j)
    enddo ; enddo
    do k=GV%ke,1,-1 ; do j=G%jsc,G%jec ; do i=G%isc,G%iec
      eta(i,j,K) = eta(i,j,K+1) + G%mask2dT(i,j) * dz(i,j,k)
    enddo ; enddo ; enddo
    ! Re-calculate for interface positions in z*-space [Z ~> m]
    do j=G%jsc,G%jec ; do i=G%isc,G%iec
      if (G%bathyT(i,j)>0.) then
        rscl = G%bathyT(i,j) / ( eta(i,j,1) + G%bathyT(i,j) )
        do K=GV%ke, 1, -1
          eta(i,j,K) = eta(i,j,K+1) + G%mask2dT(i,j) * dz(i,j,k) * rscl
        enddo
      endif
    enddo ; enddo
  else
    call MOM_error(FATAL, "NW2 tracers assume Boussinesq mode")
  endif

  do m=1,CS%ntr
    dt_x_rate = dt * CS%restore_rate(m)
    !$OMP parallel do default(shared) private(target_value)
    do k=1,GV%ke ; do j=G%jsc,G%jec ; do i=G%isc,G%iec
      target_value = nw2_tracer_dist(m, G, GV, eta, i, j, k)
      CS%tr(i,j,k,m) = CS%tr(i,j,k,m) + G%mask2dT(i,j) * dt_x_rate * ( target_value - CS%tr(i,j,k,m) )
    enddo ; enddo ; enddo
  enddo

end procedure nw2_tracer_column_physics
module procedure nw2_tracer_dist
  real :: pi ! 3.1415... [nondim]
  real :: x, y, z ! non-dimensional relative positions [nondim]
  pi = 2.*acos(0.)
  x = ( G%geolonT(i,j) - G%west_lon ) / G%len_lon ! 0 ... 1
  y = -G%geolatT(i,j) / G%south_lat ! -1 ... 1
  z = - 0.5 * ( eta(i,j,K) + eta(i,j,K+1) ) / GV%max_depth ! 0 ... 1
  select case ( mod(m-1,3) )
  case (0) ! sin(2 pi x/L)
    nw2_tracer_dist = sin( 2.0 * pi * x )
  case (1) ! y/L
    nw2_tracer_dist = y
  case (2) ! -z/L
    nw2_tracer_dist = -z
  case default
    stop 'This should not happen. Died in nw2_tracer_dist()!'
  end select
  nw2_tracer_dist = nw2_tracer_dist * G%mask2dT(i,j)
end procedure nw2_tracer_dist
module procedure nw2_tracers_end
  if (associated(CS)) then
    if (associated(CS%tr)) deallocate(CS%tr)
    deallocate(CS)
  endif
end procedure nw2_tracers_end
end submodule nw2_tracers_s
