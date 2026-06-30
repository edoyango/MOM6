submodule (Recon1d_PLM_CWK) Recon1d_PLM_CWK_s
  implicit none
contains
module procedure reconstruct
  real :: slp ! The PLM slopes (difference across cell) [A]
  real :: sigma_l, sigma_c, sigma_r ! Left, central and right slope estimates as
  real :: u_min, u_max ! Minimum and maximum value across cell [A]
  real :: u_l, u_r, u_c ! Left, right, and center values [A]
  integer :: k, n
  n = this%n

  ! Loop over all cells
  do k = 1, n
    this%u_mean(k) = u(k)
  enddo

  ! Boundary cells use PCM
  this%ul(1) = u(1)
  this%ur(1) = u(1)

  ! Loop over interior cells
  do k = 2, n-1
    u_l = u(k-1)
    u_c = u(k)
    u_r = u(k+1)

    ! Side differences
    sigma_r = u_r - u_c
    sigma_l = u_c - u_l

    ! This is the second order slope given by equation 1.7 of
    ! Piecewise Parabolic Method, Colella and Woodward (1984),
    ! but for uniform resolution.
    sigma_c = 0.5 * ( u_r - u_l )

    ! Limit slope so that reconstructions are bounded by neighbors
    u_min = min( u_l, u_c, u_r )
    u_max = max( u_l, u_c, u_r )

    if ( (sigma_l * sigma_r) > 0.0 ) then
      ! This limits the slope so that the edge values are bounded by the two cell averages spanning the edge
      slp = sign( min( abs(sigma_c), 2.*min( u_c - u_min, u_max - u_c ) ), sigma_c )
    else
      ! Extrema in the mean values require a PCM reconstruction
      slp = 0.0
    endif

    ! Left edge
    u_min = min( u_c, u_l )
    u_max = max( u_c, u_l )
    u_l = u_c - 0.5 * slp
    this%ul(k) = max( min( u_l, u_max), u_min )

    ! Right edge
    u_min = min( u_c, u_r )
    u_max = max( u_c, u_r )
    u_r = u_c + 0.5 * slp
    this%ur(k) = max( min( u_r, u_max), u_min )
  enddo

  ! Boundary cells use PCM
  this%ul(n) = u(n)
  this%ur(n) = u(n)

end procedure reconstruct
end submodule Recon1d_PLM_CWK_s
