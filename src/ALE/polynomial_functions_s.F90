submodule (polynomial_functions) polynomial_functions_s
  implicit none
contains
module procedure evaluation_polynomial
  integer :: k
  real    :: f    ! value of polynomial at x in arbitrary units [A]
  f = 0.0
  do k = 1,ncoef
    f = f + coeff(k) * ( x**(k-1) )
  enddo

  evaluation_polynomial = f

end procedure evaluation_polynomial
module procedure first_derivative_polynomial
  integer                               :: k
  real                                  :: f    ! value of the derivative at x in [A H-1]
  f = 0.0
  do k = 2,ncoef
    f = f + REAL(k-1)*coeff(k) * ( x**(k-2) )
  enddo

  first_derivative_polynomial = f

end procedure first_derivative_polynomial
module procedure integration_polynomial
  integer :: k
  real    :: integral  ! The integral of the polynomial over the specified range in [A H]
  integral = 0.0

  do k = 1,npoly+1
    integral = integral + Coeff(k) * (xi1**k - xi0**k) / real(k)
  enddo
!
!One non-answer-changing way of unrolling the above is:
!  k=1
!  integral = integral + Coeff(k) * (xi1**k - xi0**k) / real(k)
!  if (npoly>=1) then
!    k=2
!    integral = integral + Coeff(k) * (xi1**k - xi0**k) / real(k)
!  endif
!  if (npoly>=2) then
!    k=3
!    integral = integral + Coeff(k) * (xi1**k - xi0**k) / real(k)
!  endif
!  if (npoly>=3) then
!    k=4
!    integral = integral + Coeff(k) * (xi1**k - xi0**k) / real(k)
!  endif
!  if (npoly>=4) then
!    k=5
!    integral = integral + Coeff(k) * (xi1**k - xi0**k) / real(k)
!  endif
!
  integration_polynomial = integral

end procedure integration_polynomial
end submodule polynomial_functions_s
