submodule (MOM_spherical_harmonics) MOM_spherical_harmonics_s
#include <MOM_memory.h>
  implicit none
contains
module procedure spherical_harmonics_forward
  integer :: Nmax ! Local copy of the maximum degree of the spherical harmonics
  integer :: Ltot ! Local copy of the number of spherical harmonics
  real, dimension(SZI_(G),SZJ_(G)) :: &
    pmn,   & ! Current associated Legendre polynomials of degree n and order m [nondim]
    pmnm1, & ! Associated Legendre polynomials of degree n-1 and order m [nondim]
    pmnm2    ! Associated Legendre polynomials of degree n-2 and order m [nondim]
  real, allocatable, dimension(:,:,:) :: &
    Snm_Re_raw, & ! Array of un-summed real spherical harmonics transform coefficients for
                  ! reproducing sums in the same arbitrary units as var, [a] or [A ~> a]
    Snm_Im_raw    ! Array of un-summed imaginary spherical harmonics transform coefficients for
                  ! reproducing sums in the same arbitrary units as var, [a] or [A ~> a]
  real :: sum_tot ! The total of all components output by the reproducing sum in the same
                  ! arbitrary units as var, [a] or [A ~> a]
  integer :: i, j, k
  integer :: is, ie, js, je, isd, ied, jsd, jed
  integer :: m, n, l

  if (.not.CS%initialized) call MOM_error(FATAL, "MOM_spherical_harmonics " // &
    "spherical_harmonics_forward: Module must be initialized before it is used.")

  if (id_clock_sht>0) call cpu_clock_begin(id_clock_sht)
  if (id_clock_sht_forward>0) call cpu_clock_begin(id_clock_sht_forward)

  Nmax = CS%ndegree ; if (present(Nd)) Nmax = Nd
  Ltot = calc_lmax(Nmax)

  is  = G%isc ; ie  = G%iec ; js  = G%jsc ; je  = G%jec
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed

  do j=jsd,jed ; do i=isd,ied
    pmn(i,j) = 0.0 ; pmnm1(i,j) = 0.0 ; pmnm2(i,j) = 0.0
  enddo ; enddo

  do l=1,Ltot ; Snm_Re(l) = 0.0 ; Snm_Im(l) = 0.0 ; enddo

  if (CS%reprod_sum) then
    allocate(Snm_Re_raw(is:ie, js:je, Ltot), source=0.0)
    allocate(Snm_Im_raw(is:ie, js:je, Ltot), source=0.0)
    do m=0,Nmax
      l = order2index(m, Nmax)

      do j=js,je ; do i=is,ie
        Snm_Re_raw(i,j,l) = var(i,j) * CS%Pmm(i,j,m+1) * CS%cos_lonT_wtd(i,j,m+1)
        Snm_Im_raw(i,j,l) = var(i,j) * CS%Pmm(i,j,m+1) * CS%sin_lonT_wtd(i,j,m+1)
        pmnm2(i,j) = 0.0
        pmnm1(i,j) = CS%Pmm(i,j,m+1)
      enddo ; enddo

      do n = m+1, Nmax ; do j=js,je ; do i=is,ie
        pmn(i,j) = &
          CS%a_recur(n+1,m+1) * CS%cos_clatT(i,j) * pmnm1(i,j) - CS%b_recur(n+1,m+1) * pmnm2(i,j)
        Snm_Re_raw(i,j,l+n-m) = var(i,j) * pmn(i,j) * CS%cos_lonT_wtd(i,j,m+1)
        Snm_Im_raw(i,j,l+n-m) = var(i,j) * pmn(i,j) * CS%sin_lonT_wtd(i,j,m+1)
        pmnm2(i,j) = pmnm1(i,j)
        pmnm1(i,j) = pmn(i,j)
      enddo ; enddo ; enddo
    enddo
  else
    do m=0,Nmax
      l = order2index(m, Nmax)

      do j=js,je ; do i=is,ie
        Snm_Re(l) = Snm_Re(l) + var(i,j) * CS%Pmm(i,j,m+1) * CS%cos_lonT_wtd(i,j,m+1)
        Snm_Im(l) = Snm_Im(l) + var(i,j) * CS%Pmm(i,j,m+1) * CS%sin_lonT_wtd(i,j,m+1)
        pmnm2(i,j) = 0.0
        pmnm1(i,j) = CS%Pmm(i,j,m+1)
      enddo ; enddo

      do n=m+1, Nmax ; do j=js,je ; do i=is,ie
        pmn(i,j) = &
          CS%a_recur(n+1,m+1) * CS%cos_clatT(i,j) * pmnm1(i,j) - CS%b_recur(n+1,m+1) * pmnm2(i,j)
        Snm_Re(l+n-m) = Snm_Re(l+n-m) + var(i,j) * pmn(i,j) * CS%cos_lonT_wtd(i,j,m+1)
        Snm_Im(l+n-m) = Snm_Im(l+n-m) + var(i,j) * pmn(i,j) * CS%sin_lonT_wtd(i,j,m+1)
        pmnm2(i,j) = pmnm1(i,j)
        pmnm1(i,j) = pmn(i,j)
      enddo ; enddo ; enddo
    enddo
  endif

  if (id_clock_sht_global_sum>0) call cpu_clock_begin(id_clock_sht_global_sum)

  if (CS%reprod_sum) then
    sum_tot = reproducing_sum(Snm_Re_raw(:,:,1:Ltot), sums=Snm_Re(1:Ltot), unscale=tmp_scale)
    sum_tot = reproducing_sum(Snm_Im_raw(:,:,1:Ltot), sums=Snm_Im(1:Ltot), unscale=tmp_scale)
    deallocate(Snm_Re_raw, Snm_Im_raw)
  else
    call sum_across_PEs(Snm_Re, Ltot)
    call sum_across_PEs(Snm_Im, Ltot)
  endif

  if (id_clock_sht_global_sum>0) call cpu_clock_end(id_clock_sht_global_sum)
  if (id_clock_sht_forward>0) call cpu_clock_end(id_clock_sht_forward)
  if (id_clock_sht>0) call cpu_clock_end(id_clock_sht)
end procedure spherical_harmonics_forward
module procedure spherical_harmonics_inverse
  integer :: Nmax ! Local copy of the maximum degree of the spherical harmonics [nondim]
  real    :: mFac ! A constant multiplier. mFac = 1 (if m==0) or 2 (if m>0) [nondim]
  real, dimension(SZI_(G),SZJ_(G)) :: &
    pmn,   & ! Current associated Legendre polynomials of degree n and order m [nondim]
    pmnm1, & ! Associated Legendre polynomials of degree n-1 and order m [nondim]
    pmnm2    ! Associated Legendre polynomials of degree n-2 and order m [nondim]
  integer :: i, j, k
  integer :: is, ie, js, je, isd, ied, jsd, jed
  integer :: m, n, l
  if (.not.CS%initialized) call MOM_error(FATAL, "MOM_spherical_harmonics " // &
    "spherical_harmonics_inverse: Module must be initialized before it is used.")

  if (id_clock_sht>0) call cpu_clock_begin(id_clock_sht)
  if (id_clock_sht_inverse>0) call cpu_clock_begin(id_clock_sht_inverse)

  Nmax = CS%ndegree ; if (present(Nd)) Nmax = Nd

  is  = G%isc ; ie  = G%iec ; js  = G%jsc ; je  = G%jec
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed

  do j=jsd,jed ; do i=isd,ied
    pmn(i,j) = 0.0 ; pmnm1(i,j) = 0.0 ; pmnm2(i,j) = 0.0
    var(i,j) = 0.0
  enddo ; enddo

  do m=0,Nmax
    mFac = sign(1.0, m-0.5)*0.5 + 1.5
    l = order2index(m, Nmax)

    do j=js,je ; do i=is,ie
      var(i,j) = var(i,j) &
        + mFac * CS%Pmm(i,j,m+1) * (  Snm_Re(l) * CS%cos_lonT(i,j,m+1) &
                                    + Snm_Im(l) * CS%sin_lonT(i,j,m+1))
      pmnm2(i,j) = 0.0
      pmnm1(i,j) = CS%Pmm(i,j,m+1)
    enddo ; enddo

    do n=m+1,Nmax ; do j=js,je ; do i=is,ie
      pmn(i,j) = &
        CS%a_recur(n+1,m+1) * CS%cos_clatT(i,j) * pmnm1(i,j) - CS%b_recur(n+1,m+1) * pmnm2(i,j)
      var(i,j) = var(i,j) &
        + mFac * pmn(i,j) * (  Snm_Re(l+n-m) * CS%cos_lonT(i,j,m+1) &
                             + Snm_Im(l+n-m) * CS%sin_lonT(i,j,m+1))
      pmnm2(i,j) = pmnm1(i,j)
      pmnm1(i,j) = pmn(i,j)
    enddo ; enddo ; enddo
  enddo

  if (id_clock_sht_inverse>0) call cpu_clock_end(id_clock_sht_inverse)
  if (id_clock_sht>0) call cpu_clock_end(id_clock_sht)
end procedure spherical_harmonics_inverse
module procedure spherical_harmonics_init
  real, parameter :: PI = 4.0*atan(1.0) ! 3.1415926... calculated as 4*atan(1) [nondim]
  real, parameter :: RADIAN = PI / 180.0 ! Degree to Radian constant [radian degree-1]
  real, dimension(SZI_(G),SZJ_(G)) :: sin_clatT ! sine of colatitude at the t-cells [nondim].
  real :: Pmm_coef ! = sqrt{ 1.0/(4.0*PI) * prod[(2k+1)/2k)] } [nondim].
  integer :: is, ie, js, je
  integer :: i, j, k
  integer :: m, n
  integer :: Nd_SAL ! Maximum degree for SAL
# include "version_variable.h"
  character(len=40) :: mdl = "MOM_spherical_harmonics" ! This module's name.
  if (CS%initialized) return
  CS%initialized = .True.

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec

  call log_version(param_file, mdl, version, "")
  call get_param(param_file, mdl, "SAL_HARMONICS_DEGREE", Nd_SAL, "", default=0, do_not_log=.true.)
  CS%ndegree = Nd_SAL
  CS%lmax = calc_lmax(CS%ndegree)
  call get_param(param_file, mdl, "SHT_REPRODUCING_SUM", CS%reprod_sum, &
                 "If true, use reproducing sums (invariant to PE layout) in inverse transform "// &
                 "of spherical harmonics. Otherwise use a simple sum of floating point numbers. ", &
                 default=.False.)

  ! Calculate recurrence relationship coefficients
  allocate(CS%a_recur(CS%ndegree+1, CS%ndegree+1), source=0.0)
  allocate(CS%b_recur(CS%ndegree+1, CS%ndegree+1), source=0.0)
  do m=0,CS%ndegree ; do n=m+1,CS%ndegree
    ! These expressione will give NaNs with 32-bit integers for n > 23170, but this is trapped elsewhere.
    CS%a_recur(n+1,m+1) = sqrt(real((2*n-1) * (2*n+1)) / real((n-m) * (n+m)))
    CS%b_recur(n+1,m+1) = sqrt((real(2*n+1) * real((n+m-1) * (n-m-1))) / (real((n-m) * (n+m)) * real(2*n-3)))
  enddo ; enddo

  ! Calculate complex exponential factors
  allocate(CS%cos_lonT_wtd(is:ie, js:je, CS%ndegree+1), source=0.0)
  allocate(CS%sin_lonT_wtd(is:ie, js:je, CS%ndegree+1), source=0.0)
  allocate(CS%cos_lonT(is:ie, js:je, CS%ndegree+1), source=0.0)
  allocate(CS%sin_lonT(is:ie, js:je, CS%ndegree+1), source=0.0)
  do m=0,CS%ndegree
    do j=js,je ; do i=is,ie
      CS%cos_lonT(i,j,m+1)     = cos(real(m) * (G%geolonT(i,j)*RADIAN))
      CS%sin_lonT(i,j,m+1)     = sin(real(m) * (G%geolonT(i,j)*RADIAN))
      CS%cos_lonT_wtd(i,j,m+1) = CS%cos_lonT(i,j,m+1) * G%areaT(i,j) / G%Rad_Earth_L**2
      CS%sin_lonT_wtd(i,j,m+1) = CS%sin_lonT(i,j,m+1) * G%areaT(i,j) / G%Rad_Earth_L**2
    enddo ; enddo
  enddo

  ! Calculate sine and cosine of colatitude
  allocate(CS%cos_clatT(is:ie, js:je), source=0.0)
  do j=js,je ; do i=is,ie
    CS%cos_clatT(i,j) = cos(0.5*PI - G%geolatT(i,j)*RADIAN)
    sin_clatT(i,j)    = sin(0.5*PI - G%geolatT(i,j)*RADIAN)
  enddo ; enddo

  ! Calculate the diagonal elements of the associated Legendre polynomials (n=m)
  allocate(CS%Pmm(is:ie,js:je,m+1), source=0.0)
  do m=0,CS%ndegree
    Pmm_coef = 1.0/(4.0*PI)
    do k=1,m ; Pmm_coef = Pmm_coef * (real(2*k+1) / real(2*k)) ; enddo
    Pmm_coef = sqrt(Pmm_coef)
    do j=js,je ; do i=is,ie
      CS%Pmm(i,j,m+1) = Pmm_coef * (sin_clatT(i,j)**m)
    enddo ; enddo
  enddo

  id_clock_sht = cpu_clock_id('(Ocean spherical harmonics)', grain=CLOCK_MODULE)
  id_clock_sht_forward = cpu_clock_id('(Ocean SHT forward)', grain=CLOCK_ROUTINE)
  id_clock_sht_inverse = cpu_clock_id('(Ocean SHT inverse)', grain=CLOCK_ROUTINE)
  id_clock_sht_global_sum = cpu_clock_id('(Ocean SHT global sum)', grain=CLOCK_LOOP)

end procedure spherical_harmonics_init
module procedure spherical_harmonics_end
  deallocate(CS%cos_clatT)
  deallocate(CS%Pmm)
  deallocate(CS%cos_lonT_wtd, CS%sin_lonT_wtd, CS%cos_lonT, CS%sin_lonT)
  deallocate(CS%a_recur, CS%b_recur)
end procedure spherical_harmonics_end
module procedure calc_lmax
  lmax = (Nd+2) * (Nd+1) / 2
end procedure calc_lmax
module procedure order2index
  l = ((Nd+1) + (Nd+1-(m-1)))*m/2 + 1
end procedure order2index
end submodule MOM_spherical_harmonics_s
