submodule (Recon1d_type) Recon1d_type_s
  implicit none
contains
module procedure x
  real :: xl, xr, xo ! Left/right bounds and guess [nondim]
  real :: fl, fr ! Left right values [A]
  real :: slp ! Difference across cell or derivative wrt nondim x [A]
  real :: f_at_x ! Value at current x [A]
  integer :: iter
  x = 0.5 ! Fall back for special conditions
  fl = this%f(k, 0.)
  fr = this%f(k, 1.)
  slp = fr - fl
  if ( ( fl - t ) * ( t - fr ) > 0. ) then
    ! t is inside the range fl..fr
    xl = 0.
    xr = 1.
    xo = ( t - this%f(k, 0.) ) / slp ! First guess by regula falsi
    f_at_x = this%f(k, xo)
    do iter = 1,10
      slp = this%dfdx(k, xo)
      x = xo - ( f_at_x - t ) / slp ! Newton-Raphson step
      if ( x < xl ) x = 0.5 * ( xl + xo ) ! Replace with bi-section
      if ( x > xr ) x = 0.5 * ( xr + xo ) ! Replace with bi-section
      f_at_x = this%f(k, x)
      if ( abs(f_at_x - t) <= 0. .or. abs(x - xo) < this%x_tolerance ) return
      if ( f_at_x < t ) xl = x ! Replace left bound
      if ( f_at_x > t ) xr = x ! Replace right bound
      xo = x
    enddo
  elseif ( abs(slp) > 0. ) then
    slp = sign(1., slp)
    ! if t>u_mean & slp=1 then x=1
    ! if t<u_mean & slp=1 then x=0
    ! if t>u_mean & slp=-1 then x=0
    ! if t<u_mean & slp=-1 then x=1
    x = 0.5 + slp * sign(0.5, t - this%u_mean(k))
  else
    ! slp=0 so estimate "direction" from neighbors
    slp = this%f(min(k+1,this%n), 0.) - this%f(max(k-1,1), 1.)
    if ( abs(slp) > 0. ) slp = sign(1., slp)
    ! if t>u_mean & slp=1 then x=1
    ! if t<u_mean & slp=1 then x=0
    ! if t>u_mean & slp=-1 then x=0
    ! if t<u_mean & slp=-1 then x=1
    ! if t=u_mean then x=0.5
    ! if slp=0 then x=0.5
    if ( abs(t - this%u_mean(k)) > 0. ) x = 0.5 + slp * sign(0.5, t - this%u_mean(k))
  endif
end procedure x
module procedure remap_to_sub_grid
  integer :: i_sub ! Index of sub-cell
  integer :: i0 ! Index into h0(1:n0), source column
  integer :: i_max ! Used to record which sub-cell is the largest contribution of a source cell
  real :: dh_max ! Used to record which sub-cell is the largest contribution of a source cell [H]
  real :: xa, xb ! Non-dimensional position within a source cell (0..1) [nondim]
  real :: dh ! The width of the sub-cell [H]
  real :: duh ! The total amount of accumulated stuff (u*h) [A H]
  real :: dh0_eff ! Running sum of source cell thickness [H]
  integer :: i0_last_thick_cell, n0
  n0 = this%n

  i0_last_thick_cell = 0
  do i0 = 1, n0
!   ul = this%f(i0, 0.)
!   ur = this%f(i0, 1.)
!   u0_min(i0) = min(ul, ur)
!   u0_max(i0) = max(ul, ur)
    if (h0(i0)>0.) i0_last_thick_cell = i0
  enddo

  ! Loop over each sub-cell to calculate average/integral values within each sub-cell.
  ! Uses: h_sub, isub_src, h0_eff
  ! Sets: u_sub, uh_sub
  xa = 0.
  dh0_eff = 0.
  u02_err = 0.
  do i_sub = 1, n0+n1

    ! Sub-cell thickness from loop above
    dh = h_sub(i_sub)

    ! Source cell
    i0 = isub_src(i_sub)

    ! Evaluate average and integral for sub-cell i_sub.
    ! Integral is over distance dh but expressed in terms of non-dimensional
    ! positions with source cell from xa to xb  (0 <= xa <= xb <= 1).
    dh0_eff = dh0_eff + dh ! Cumulative thickness within the source cell
    if (h0(i0)>0.) then
      xb = dh0_eff / h0(i0) ! This expression yields xa <= xb <= 1.0
      xb = min(1., xb) ! This is only needed when the total target column is wider than the source column
      u_sub(i_sub) = this%average( i0, xa, xb )
    else ! Vanished cell
      xb = 1.
      u_sub(i_sub) = u0(i0)
    endif
!   u_sub(i_sub) = max( u_sub(i_sub), u0_min(i0) )
!   u_sub(i_sub) = min( u_sub(i_sub), u0_max(i0) )
    uh_sub(i_sub) = dh * u_sub(i_sub)

    if (isub_src(i_sub+1) /= i0) then
      ! If the next sub-cell is in a different source cell, reset the position counters
      dh0_eff = 0.
      xa = 0.
    else
      xa = xb ! Next integral will start at end of last
    endif

  enddo
  i_sub = n0+n1+1
  ! Sub-cell thickness from loop above
  dh = h_sub(i_sub)
  ! Source cell
  i0 = isub_src(i_sub)

  ! Evaluate average and integral for sub-cell i_sub.
  ! Integral is over distance dh but expressed in terms of non-dimensional
  ! positions with source cell from xa to xb  (0 <= xa <= xb <= 1).
  dh0_eff = dh0_eff + dh ! Cumulative thickness within the source cell
  if (h0(i0)>0.) then
    xb = dh0_eff / h0(i0) ! This expression yields xa <= xb <= 1.0
    xb = min(1., xb) ! This is only needed when the total target column is wider than the source column
    u_sub(i_sub) = this%average( i0, xa, xb )
  else ! Vanished cell
    xb = 1.
    u_sub(i_sub) = u0(i0)
  endif
! u_sub(i_sub) = max( u_sub(i_sub), u0_min(i0) )
! u_sub(i_sub) = min( u_sub(i_sub), u0_max(i0) )
  uh_sub(i_sub) = dh * u_sub(i_sub)

  ! Loop over each source cell substituting the integral/average for the thickest sub-cell (within
  ! the source cell) with the residual of the source cell integral minus the other sub-cell integrals
  ! aka a genius algorithm for accurate conservation when remapping from Robert Hallberg (\@Hallberg-NOAA).
  ! Uses: i0_last_thick_cell, isrc_max, h_sub, isrc_start, isrc_end, uh_sub, u0, h0
  ! Updates: uh_sub
  do i0 = 1, i0_last_thick_cell
    i_max = isrc_max(i0)
    dh_max = h_sub(i_max)
    if (dh_max > 0.) then
      ! duh will be the sum of sub-cell integrals within the source cell except for the thickest sub-cell.
      duh = 0.
      do i_sub = isrc_start(i0), isrc_end(i0)
        if (i_sub /= i_max) duh = duh + uh_sub(i_sub)
      enddo
      uh_sub(i_max) = u0(i0)*h0(i0) - duh
      u02_err = u02_err + max( abs(uh_sub(i_max)), abs(u0(i0)*h0(i0)), abs(duh) )
    endif
  enddo

  ! This should not generally be used
  if (this%check) then
    if ( this%check_reconstruction(h0, u0) ) stop 910 ! A debugger is required to understand why this failed
  endif

end procedure remap_to_sub_grid
module procedure a_set_debug
  this%debug = .true.

end procedure a_set_debug
end submodule Recon1d_type_s
