submodule (user_change_diffusivity) user_change_diffusivity_s
#include <MOM_memory.h>
  implicit none
contains
module procedure user_change_diff
  real :: Rcv(SZI_(G),SZK_(GV)) ! The coordinate density in layers [R ~> kg m-3].
  real :: p_ref(SZI_(G))       ! An array of tv%P_Ref pressures [R L2 T-2 ~> Pa].
  real :: rho_fn      ! The density dependence of the input function, 0-1 [nondim].
  real :: lat_fn      ! The latitude dependence of the input function, 0-1 [nondim].
  logical :: use_EOS  ! If true, density is calculated from T & S using an
  logical :: store_Kd_add  ! Save the added diffusivity as a diagnostic if true.
  integer, dimension(2) :: EOSdom ! The i-computational domain for the equation of state
  integer :: i, j, k, is, ie, js, je, nz
  integer :: isd, ied, jsd, jed
  character(len=200) :: mesg
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed

  if (.not.associated(CS)) call MOM_error(FATAL,"user_set_diffusivity: "//&
         "Module must be initialized before it is used.")

  if (.not.CS%initialized) call MOM_error(FATAL,"user_set_diffusivity: "//&
         "Module must be initialized before it is used.")

  use_EOS = associated(tv%eqn_of_state)
  if (.not.use_EOS) return
  store_Kd_add = .false.
  if (present(Kd_int_add)) store_Kd_add = associated(Kd_int_add)

  if (.not.range_OK(CS%lat_range)) then
    write(mesg, '(4(1pe15.6))') CS%lat_range(1:4)
    call MOM_error(FATAL, "user_set_diffusivity: bad latitude range: \n  "//&
                    trim(mesg))
  endif
  if (.not.range_OK(CS%rho_range)) then
    write(mesg, '(4(1pe15.6))') CS%rho_range(1:4)
    call MOM_error(FATAL, "user_set_diffusivity: bad density range: \n  "//&
                    trim(mesg))
  endif

  if (store_Kd_add) Kd_int_add(:,:,:) = 0.0

  do i=is,ie ; p_ref(i) = tv%P_Ref ; enddo
  EOSdom(:) = EOS_domain(G%HI)
  do j=js,je
    if (present(T_f) .and. present(S_f)) then
      do k=1,nz
        call calculate_density(T_f(:,j,k), S_f(:,j,k), p_ref, Rcv(:,k), tv%eqn_of_state, EOSdom)
      enddo
    else
      do k=1,nz
        call calculate_density(tv%T(:,j,k), tv%S(:,j,k), p_ref, Rcv(:,k), tv%eqn_of_state, EOSdom)
      enddo
    endif

    if (present(Kd_lay)) then
      do k=1,nz ; do i=is,ie
        if (CS%use_abs_lat) then
          lat_fn = val_weights(abs(G%geoLatT(i,j)), CS%lat_range)
        else
          lat_fn = val_weights(G%geoLatT(i,j), CS%lat_range)
        endif
        rho_fn = val_weights(Rcv(i,k), CS%rho_range)
        if (rho_fn * lat_fn > 0.0) &
          Kd_lay(i,j,k) = Kd_lay(i,j,k) + CS%Kd_add * rho_fn * lat_fn
      enddo ; enddo
    endif
    if (present(Kd_int)) then
      do K=2,nz ; do i=is,ie
        if (CS%use_abs_lat) then
          lat_fn = val_weights(abs(G%geoLatT(i,j)), CS%lat_range)
        else
          lat_fn = val_weights(G%geoLatT(i,j), CS%lat_range)
        endif
        rho_fn = val_weights( 0.5*(Rcv(i,k-1) + Rcv(i,k)), CS%rho_range)
        if (rho_fn * lat_fn > 0.0) then
          Kd_int(i,j,K) = Kd_int(i,j,K) + CS%Kd_add * rho_fn * lat_fn
          if (store_Kd_add) Kd_int_add(i,j,K) = CS%Kd_add * rho_fn * lat_fn
        endif
      enddo ; enddo
    endif
  enddo

end procedure user_change_diff
module procedure range_OK
  OK = ((range(1) <= range(2)) .and. (range(2) <= range(3)) .and. &
        (range(3) <= range(4)))

end procedure range_OK
module procedure val_weights
  real :: x   ! A nondimensional number between 0 and 1 [nondim].
  ans = 0.0
  if ((val > range(1)) .and. (val < range(4))) then
    if (val < range(2)) then
      ! x goes from 0 to 1; ans goes from 0 to 1, with 0 derivatives at the ends.
      x = (val - range(1)) / (range(2) - range(1))
      ans = x**2 * (3.0 - 2.0 * x)
    elseif (val > range(3)) then
      ! x goes from 0 to 1; ans goes from 0 to 1, with 0 derivatives at the ends.
      x = (range(4) - val) / (range(4) - range(3))
      ans = x**2 * (3.0 - 2.0 * x)
    else
      ans = 1.0
    endif
  endif

end procedure val_weights
module procedure user_change_diff_init
# include "version_variable.h"
  character(len=40)  :: mdl = "user_set_diffusivity"  ! This module's name.
  character(len=200) :: mesg
  if (associated(CS)) then
    call MOM_error(WARNING, "diabatic_entrain_init called with an associated "// &
                            "control structure.")
    return
  endif
  allocate(CS)

  CS%initialized = .true.
  CS%diag => diag

  ! Read all relevant parameters and write them to the model log.
  call log_version(param_file, mdl, version, "")
  call get_param(param_file, mdl, "USER_KD_ADD", CS%Kd_add, &
                 "A user-specified additional diffusivity over a range of "//&
                 "latitude and density.", default=0.0, units="m2 s-1", scale=GV%m2_s_to_HZ_T)
  if (CS%Kd_add /= 0.0) then
    call get_param(param_file, mdl, "USER_KD_ADD_LAT_RANGE", CS%lat_range(:), &
                 "Four successive values that define a range of latitudes "//&
                 "over which the user-specified extra diffusivity is "//&
                 "applied.  The four values specify the latitudes at "//&
                 "which the extra diffusivity starts to increase from 0, "//&
                 "hits its full value, starts to decrease again, and is "//&
                 "back to 0.", units="degrees_N", defaults=(/-1.0e9,-1.0e9,-1.0e9,-1.0e9/))
    call get_param(param_file, mdl, "USER_KD_ADD_RHO_RANGE", CS%rho_range(:), &
                 "Four successive values that define a range of potential "//&
                 "densities over which the user-given extra diffusivity "//&
                 "is applied.  The four values specify the density at "//&
                 "which the extra diffusivity starts to increase from 0, "//&
                 "hits its full value, starts to decrease again, and is "//&
                 "back to 0.", units="kg m-3", defaults=(/-1.0e9,-1.0e9,-1.0e9,-1.0e9/),&
                 scale=US%kg_m3_to_R)
    call get_param(param_file, mdl, "USER_KD_ADD_USE_ABS_LAT", CS%use_abs_lat, &
                 "If true, use the absolute value of latitude when "//&
                 "checking whether a point fits into range of latitudes.", &
                 default=.false.)
  endif

  if (.not.range_OK(CS%lat_range)) then
    write(mesg, '(4(1pe15.6))') CS%lat_range(1:4)
    call MOM_error(FATAL, "user_set_diffusivity: bad latitude range: \n  "//&
                    trim(mesg))
  endif
  if (.not.range_OK(CS%rho_range)) then
    write(mesg, '(4(1pe15.6))') CS%rho_range(1:4)
    call MOM_error(FATAL, "user_set_diffusivity: bad density range: \n  "//&
                    trim(mesg))
  endif

end procedure user_change_diff_init
module procedure user_change_diff_end
  if (associated(CS)) deallocate(CS)

end procedure user_change_diff_end
end submodule user_change_diffusivity_s
