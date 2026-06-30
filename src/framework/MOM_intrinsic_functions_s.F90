submodule (MOM_intrinsic_functions) MOM_intrinsic_functions_s
  implicit none
contains
module procedure invcosh
#ifdef __INTEL_COMPILER
  invcosh = acosh(x)
#else
  invcosh = log(x+sqrt(x*x-1))
#endif

end procedure invcosh
module procedure cuberoot
  real :: asx ! The absolute value of x rescaled by an integer power of 8 to put it into
  real :: root_asx ! The cube root of asx [B]
  real :: ra_3 ! root_asx cubed [B3]
  real :: num ! The numerator of an expression for the evolving estimate of the cube root of asx
  real :: den ! The denominator of an expression for the evolving estimate of the cube root of asx
  real :: num_prev ! The numerator of an expression for the previous iteration of the evolving estimate
  real :: np_3 ! num_prev cubed  [B3 D3]
  real :: den_prev ! The denominator of an expression for the previous iteration of the evolving estimate of
  real :: dp_3 ! den_prev cubed  [C3]
  real :: r0  ! Initial value of the iterative solver. [B C]
  real :: r0_3 ! r0 cubed [B3 C3]
  integer :: itt
  integer(kind=int64) :: e_x, s_x
  if ((x >= 0.0) .eqv. (x <= 0.0)) then
    ! Return 0 for an input of 0, or NaN for a NaN input.
    root = x
  else
    call rescale_cbrt(x, asx, e_x, s_x)

    !   Iteratively determine root_asx = asx**1/3 using Halley's method and then Newton's method,
    ! noting that Halley's method onverges monotonically and needs no bounding.  Halley's method is
    ! slightly more complicated that Newton's method, but converges in a third fewer iterations.
    !   Keeping the estimates in a fractional form Root = num / den allows this calculation with
    ! no real divisions during the iterations before doing a single real division at the end,
    ! and it is therefore more computationally efficient.

    ! This first estimate gives the same magnitude of errors for 0.125 and 1.0 after two iterations.
    ! The first iteration is applied explicitly.
    r0 = 0.707106
    r0_3 = r0 * r0 * r0
    num = r0 * (r0_3 + 2.0 * asx)
    den = 2.0 * r0_3 + asx

    do itt=1,2
      ! Halley's method iterates estimates as Root = Root * (Root**3 + 2.*asx) / (2.*Root**3 + asx).
      num_prev = num ; den_prev = den

      ! Pre-compute these as integer powers, to avoid `pow()`-like intrinsics.
      np_3 = num_prev * num_prev * num_prev
      dp_3 = den_prev * den_prev * den_prev

      num = num_prev * (np_3 + 2.0 * asx * dp_3)
      den = den_prev * (2.0 * np_3 + asx * dp_3)
      ! Equivalent to:  root_asx = root_asx * (root_asx**3 + 2.*asx) / (2.*root_asx**3 + asx)
    enddo
    ! At this point the error in root_asx is better than 1 part in 3e14.
    root_asx = num / den

    ! One final iteration with Newton's method polishes up the root and gives a solution
    ! that is within the last bit of the true solution.
    ra_3 = root_asx * root_asx * root_asx
    root_asx = root_asx - (ra_3 - asx) / (3.0 * (root_asx * root_asx))

    root = descale(root_asx, e_x, s_x)
  endif
end procedure cuberoot
module procedure rescale_cbrt
  integer(kind=int64) :: xb
  integer(kind=int64) :: e_a
  integer(kind=int64) :: e_x
  xb = transfer(a, 1_int64)
  s_a = ibits(xb, signbit, 1)
  e_a = ibits(xb, expbit, explen) - bias

  ! The floating-point form of `a` with exponent `e` is
  !
  !   a = s * (1 + m) * 2**e
  !
  ! where (1+m) ∈ [1,2).  We want to split 2**e so that (1+m) is rescaled to
  ! the range [0.125, 1); that is, [2**-3, 2**0).
  !
  ! First decompose the exponent `e` into quotient-remainder form:
  !
  !   e = 3⌊e/3⌋ + modulo(e,3)
  !
  ! Since modulo(e,3) ∈ {0,1,2}, the second term of the following expression is
  ! in {-3,-2,-1}.
  !
  !   e = 3 * (⌊e/3⌋ + 1) + (modulo(e,3) - 3).
  !
  ! Here, (modulo(e,3) - 3) is in the range [2**-3, 1) and holds the
  ! floating-point exponent of `x`.
  !
  ! Fortran integer division is round-to-zero.  To convert to floor division,
  ! we use the sign() intrinsic to shift negative values downward.
  !
  !   ⌊e/3⌋ = (e + sign(1,e) - 1) / 3
  !
  ! ⌊e/3⌋ + 1 reduces to the form below.  This is what we call the integral
  ! cube-root of `a` in the description above.

  e_r = (e_a + sign(1_int64, e_a) + 2) / 3

  ! modulo() is not implemented on all systems, so compute the remainder as
  ! r = n - 3*q.

  e_x = e_a - e_r * 3

  ! Insert the new 11-bit exponent into xb and write to x and extend the
  ! bitcount to 12, so that the sign bit is zero and x is always positive.
  call mvbits(e_x + bias, 0, explen + 1, xb, fraclen)
  x = transfer(xb, 1.)
end procedure rescale_cbrt
module procedure descale
  integer(kind=int64) :: xb
  integer(kind=int64) :: e_x
  xb = transfer(x, 1_int64)
  e_x = ibits(xb, expbit, explen)
  call mvbits(e_a + e_x, 0, explen, xb, expbit)
  call mvbits(s_a, 0, 1, xb, signbit)
  a = transfer(xb, 1.)
end procedure descale
module procedure intrinsic_functions_unit_tests
  real :: testval  ! A test value for self-consistency testing [nondim]
  logical :: v
  integer :: n
  fail = .false.
  v = verbose
  write(stdout,*) '==== MOM_intrinsic_functions: intrinsic_functions_unit_tests ==='

  fail = fail .or. Test_cuberoot(v, 1.2345678901234e9)
  fail = fail .or. Test_cuberoot(v, -9.8765432109876e-21)
  fail = fail .or. Test_cuberoot(v, 64.0)
  fail = fail .or. Test_cuberoot(v, -0.5000000000001)
  fail = fail .or. Test_cuberoot(v, 0.0)
  fail = fail .or. Test_cuberoot(v, 1.0)
  fail = fail .or. Test_cuberoot(v, 0.125)
  fail = fail .or. Test_cuberoot(v, 0.965)
  fail = fail .or. Test_cuberoot(v, 1.0 - epsilon(1.0))
  fail = fail .or. Test_cuberoot(v, 1.0 - 0.5*epsilon(1.0))

  testval = 1.0e-99
  v = .false.
  do n=-160,160
    fail = fail .or. Test_cuberoot(v, testval)
    testval = (-2.908 * (1.414213562373 + 1.2345678901234e-5*n)) * testval
  enddo
end procedure intrinsic_functions_unit_tests
module procedure Test_cuberoot
  real :: diff ! The difference between val and the cube root of its cube [A].
  diff = val - cuberoot(val)**3
  Test_cuberoot = (abs(diff) > 2.0e-15*abs(val))

  if (Test_cuberoot) then
    write(stdout, '("For val = ",ES22.15,", (val - cuberoot(val**3))) = ",ES9.2," <-- FAIL")') val, diff
  elseif (verbose) then
    write(stdout, '("For val = ",ES22.15,", (val - cuberoot(val**3))) = ",ES9.2)') val, diff

  endif
end procedure Test_cuberoot
end submodule MOM_intrinsic_functions_s
