! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> A module with intrinsic functions that are used by MOM but are not supported
!! by some compilers.

#include "intrinsics/MOM_exp.h"

module MOM_intrinsic_functions

use, intrinsic :: iso_fortran_env, only : stdout => output_unit, stderr => error_unit
use, intrinsic :: iso_fortran_env, only : int32, int64
use, intrinsic :: ieee_arithmetic, only : ieee_rint
use MOM_exp_data_n128, only : ndiv, idiv_scale_lookup, idiv_residual_lookup

implicit none ; private

public :: invcosh, cuberoot, nth_root, exp_repro
public :: intrinsic_functions_unit_tests

! Floating point model, if bit layout from high to low is (sign, exp, frac)

integer, parameter :: bias = maxexponent(1.) - 1
  !< The double precision exponent offset
integer, parameter :: signbit = storage_size(1.) - 1
  !< Position of sign bit
integer, parameter :: explen = 1 + ceiling(log(real(bias))/log(2.))
  !< Bit size of exponent
integer, parameter :: expbit = signbit - explen
  !< Position of lowest exponent bit
integer, parameter :: fraclen = expbit
  !< Length of fractional part

! exp_repro() floating point model and constants.
! NOTE (E3): signbit/expbit above and their exp_repro equivalents below compute
! the same quantities via different formulas (storage_size()-based vs
! digits()-based); despite both evaluating to 63/52 for real64, Fortran
! forbids the duplicate top-level names, hence the "_exp" suffix here.

real, parameter :: real_mold = 0.
  !< Real mold for transfers and numerical format queries, for exp_repro()

integer, parameter :: int_kind &
    = merge(int64, int32, storage_size(real_mold) > storage_size(0_int32))
  !< Integer kind with the same storage size as default real
integer(kind=int_kind), parameter :: int_mold = 0
  !< Integer mold value

integer, parameter :: expbit_exp = digits(real_mold) - 1
  !< Position of lowest exponent bit, for exp_repro()
  ! NOTE: digits() includes the implicit leading digit
integer, parameter :: signbit_exp = storage_size(real_mold) - 1
  !< Position of sign bit, for exp_repro()
integer, parameter :: expwidth_exp = signbit_exp - expbit_exp
  !< Number of exponent bits, for exp_repro()
integer, parameter :: expbias_exp = maxexponent(real_mold) - 1
  !< Exponent bias, for exp_repro()

! IEEE 754 special values, for exp_repro()
integer(kind=int_kind), parameter :: pos_inf_bits &
    = ishft(2_int_kind**expwidth_exp - 1_int_kind, expbit_exp)
  !< IEEE +Inf bit pattern
integer(kind=int_kind), parameter :: neg_inf_bits &
    = ior(pos_inf_bits, ishft(-1_int_kind, signbit_exp))
  !< IEEE -Inf bit pattern

! Fast integer rounding offset, for exp_repro()
real, parameter :: round_bias = 1.5 * 2_int_kind**(digits(real_mold) - 1)
  !< Binary offset used to trigger rounding of fractional values
integer(int_kind), parameter :: round_bias_bits &
    = transfer(round_bias, int_mold)
  !< Bit representation of rounding bias

!> Remez coefficients for (exp(x) - 1) / x on [-ln2/256, ln2/256] [nondim]
!! Hoisted to module scope (rather than function-local, as in the original
!! submodule): a function-local parameter array is implemented as a static
!! variable, which nvfortran's -Minline refuses to cross-file inline
!! ("subprogram not inlined -- static variable during crossing files"),
!! independent of submodule structure -- found while testing experiment E4.
real, parameter :: exp_remez_c(0:4) = [ &
    1.0, &
    0.4999999999999766853164828717126511037349700927734375, &
    0.166666666666670015839457619222230277955532073974609375, &
    4.1666679392304360740606483659576042555272579193115234375e-2, &
    8.3333340579158539374038383584775147028267383575439453125e-3 &
]

contains

!> Evaluate the inverse cosh, either using a math library or an
!! equivalent expression
function invcosh(x)
  real, intent(in) :: x !< The argument of the inverse of cosh [nondim].  NaNs will
                        !! occur if x<1, but there is no error checking
  real :: invcosh  ! The inverse of cosh of x [nondim]

#ifdef __INTEL_COMPILER
  invcosh = acosh(x)
#else
  invcosh = log(x+sqrt(x*x-1))
#endif

end function invcosh


!> Returns the cube root of a real argument at roundoff accuracy, in a form that works properly with
!! rescaling of the argument by integer powers of 8.  If the argument is a NaN, a NaN is returned.
elemental function cuberoot(x) result(root)
  !$omp declare target
  real, intent(in) :: x !< The argument of cuberoot in arbitrary units cubed [A3]
  real :: root !< The real cube root of x in arbitrary units [A]

  real :: asx ! The absolute value of x rescaled by an integer power of 8 to put it into
              ! the range from 0.125 < asx <= 1.0, in ambiguous units cubed [B3]
  real :: root_asx ! The cube root of asx [B]
  real :: ra_3 ! root_asx cubed [B3]
  real :: num ! The numerator of an expression for the evolving estimate of the cube root of asx
              ! in arbitrary units that can grow or shrink with each iteration [B C]
  real :: den ! The denominator of an expression for the evolving estimate of the cube root of asx
              ! in arbitrary units that can grow or shrink with each iteration [C]
  real :: num_prev ! The numerator of an expression for the previous iteration of the evolving estimate
              ! of the cube root of asx in arbitrary units that can grow or shrink with each iteration [B D]
  real :: np_3 ! num_prev cubed  [B3 D3]
  real :: den_prev ! The denominator of an expression for the previous iteration of the evolving estimate of
              ! the cube root of asx in arbitrary units that can grow or shrink with each iteration [D]
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
end function cuberoot


!> Bit-stable n-th root of x for x in (0, +inf) and integer n >= 1, suitable
!! for evaluation inside `!$omp target` / `do concurrent` offloaded regions.
!!
!! Lowering `x**(1.0/n)` via the compiler produces `exp((1.0/n)*log(x))` — two
!! transcendentals whose last-bit rounding differs between host libm and CUDA
!! libdevice. This routine avoids that path entirely: it uses fixed-iteration
!! Newton on y^n - x = 0, with y^(n-1) evaluated as repeated multiplication,
!! and one bit-precision-polishing iteration at the end.
!!
!! For x in [0, 1] (the case in `MOM_barotropic.F90`'s `bt_rem = av_rem**Instep`)
!! convergence is rapid because the linear initial guess y0 = 1 - (1-x)/n is
!! already within a few percent of the true root.
elemental function nth_root(x, n) result(root)
  !$omp declare target
  real,    intent(in) :: x  !< Argument, x >= 0 [arbitrary]
  integer, intent(in) :: n  !< Root index, n >= 1
  real :: root              !< x**(1/n) in the same units as x

  integer, parameter :: maxitt = 20  ! Fixed (deterministic) iteration count
  real    :: y, ypow_nm1
  real    :: x_n_r, x_nm1_r
  integer :: itt, k

  ! Trivial cases — return early to keep loop tight and avoid 0/0 below.
  if (n <= 1) then
    root = x
    return
  endif
  if (x == 0.0) then
    root = 0.0
    return
  endif

  x_n_r   = real(n)
  x_nm1_r = real(n - 1)

  ! Linear initial guess valid for x in [0, 1] and decent for moderate x > 1.
  ! For our caller (av_rem in [0, 1], typically near 1), this is within ~1%.
  y = 1.0 - (1.0 - x) / x_n_r

  ! Newton iteration:  y_{k+1} = ((n-1)*y_k + x / y_k^{n-1}) / n
  ! All ops are *, +, /. Integer power y^{n-1} is repeated multiplication,
  ! so there is no `pow`/`exp/log` lowering anywhere in the iteration.
  do itt = 1, maxitt
    ypow_nm1 = 1.0
    do k = 1, n - 1
      ypow_nm1 = ypow_nm1 * y
    enddo
    y = (x_nm1_r * y + x / ypow_nm1) / x_n_r
  enddo

  root = y
end function nth_root


!> Rescale `a` to the range [0.125, 1) and compute its cube-root exponent.
!!
!! This function decomposes `a` into the form `s * x * 2**e` so that `x` is
!! in the desired range.  This is accomplished by computing the integral cube
!! root of `e` (as a division) and applying the residual to `x`.
pure subroutine rescale_cbrt(a, x, e_r, s_a)
  !$omp declare target
  real, intent(in) :: a
    !< The number to be rescaled for cube-root computation [A3]
  real, intent(out) :: x
    !< The rescaled value of `a` in the range [0.125, 1) [B3]
  integer(kind=int64), intent(out) :: e_r
    !< The integral component of the cube-root exponent of `a`.
  integer(kind=int64), intent(out) :: s_a
    !< Sign bit of `a`.  A nonzero value indicates negative sign.

  integer(kind=int64) :: xb
    ! Floating point integer representation of `a`
  integer(kind=int64) :: e_a
    ! Exponent of `a`
  integer(kind=int64) :: e_x
    ! Exponent of `x`

  ! Pack bits of a into xb and extract its exponent and sign.
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
end subroutine rescale_cbrt


!> Undo the rescaling of a real number back to its original base.
pure function descale(x, e_a, s_a) result(a)
  !$omp declare target
  real, intent(in) :: x
    !< The rescaled value which is to be restored in ambiguous units [B]
  integer(kind=int64), intent(in) :: e_a
    !< Exponent of the unscaled value
  integer(kind=int64), intent(in) :: s_a
    !< Sign bit of the unscaled value
  real :: a
    !< Restored value with the corrected exponent and sign in arbitrary units [A]

  integer(kind=int64) :: xb
    ! Bit-packed real number into integer form
  integer(kind=int64) :: e_x
    ! Biased exponent of x

  ! Apply the corrected exponent and sign to x.
  xb = transfer(x, 1_int64)
  e_x = ibits(xb, expbit, explen)
  call mvbits(e_a + e_x, 0, explen, xb, expbit)
  call mvbits(s_a, 0, 1, xb, signbit)
  a = transfer(xb, 1.)
end function descale


!> Reproducible exponential function
!!
!! Compute exp(x) with bitwise reproducibility across platforms.
!! (E3: moved here from the former submodule MOM_exp, to test whether
!! -Mextract/-Minline cross-file inlining requires an ordinary module
!! procedure rather than a submodule procedure -- see MOM_exp_data_n128 for
!! the lookup tables this depends on.)
elemental function exp_repro(x) result(a)
  !$omp declare target
  real, intent(in) :: x
    !< Input value
  real :: a
    !< exp(x)

  ! ln2 estimates
  real, parameter :: ln2 = 0.693147180559945309417232121458176568
    !< ln2: 0.693147180559945309417232... [nondim]
  real, parameter :: I_ln2 = 1.44269504088896340735992468100189214
    !< 1 / ln2: 1.4426950408889634073599... [nondim]

  ! The max and min x values between which exp(x) remains finite.
  real, parameter :: xmax &
      = real(maxexponent(real_mold) + 1) * ln2
    !< Largest x value before exp(x) overflow [nondim]
  real, parameter :: xmin &
      = real(minexponent(real_mold) - digits(real_mold) - 1) * ln2
    !< Smallest x value before exp(x) underflow [nondim]

  ! Double-real precision of ln2 used in Cody-Waite range reduction
  ! NOTE: This split assumes real64 precision.
  real, parameter :: ln2_hi = 0.69314718036912381649017333984375
    !< Upper 32 bits of ln2: 6.93147180369123816490e-01 [nondim]
  real, parameter :: ln2_lo = 1.90821492927058770002e-10
    !< Lower precision bits of ln2: 1.90821492927058770002e-10 [nondim]

  ! Subdivide [-ln2/2, ln2/2] into ndiv subintervals to reduce approximation range.
  ! This allows for a smaller polynomial at the cost of a lookup table.
  ! ndiv and the lookup tables are imported from MOM_exp_data_n128.
  real, parameter :: I_ndiv = 1. / real(ndiv)
    !< 1 / ndiv [nondim]
  real, parameter :: n_ln2 = ndiv * I_ln2
    !< ndiv / ln2 [nondim]
  real, parameter :: ln2_ndiv_hi = I_ndiv * ln2_hi
    !< Upper 32 bits of ln2 / ndiv [nondim]
  real, parameter :: ln2_ndiv_lo = I_ndiv * ln2_lo
    !< Lower precision bits of ln2 / ndiv [nondim]
  integer(int_kind), parameter :: idiv_mask = int(ndiv - 1, int_kind)
    !< Used for fast modulo of ndiv

  ! Range of K = nint(x / ln2) for which direct exponent scaling is safe.
  ! Beyond this range, a bias is applied to handle subnormals and overflow.
  ! NOTE: Fortran exponent is defined as one less than IEEE exponent.
  integer, parameter :: Kmin = minexponent(real_mold) + 1
    !< Minimum K before subnormal scaling is needed
    !! Kmin = (minexponent() - 1) + 1 (min exp(r)) (+1 safety buffer)
  integer, parameter :: Kmax = maxexponent(real_mold) - 2
    !< Maximum K before overflow scaling is needed
    !! Kmax = (maxexponent() - 1) - 0 (max exp(r)) (-1 safety buffer)
  integer(kind=int_kind), parameter :: Kbias = maxexponent(real_mold) - 2
    !< Exponent adjustment used for overflow and subnormal scaling
    !! Any bias which rescales 2**K exp(r) to O(1) works here.

  ! Nonfinite testing
  logical :: nonfinite
    ! True if input is a nonfinite float (+/-Inf, NaN)
  integer(kind=int_kind) :: xb
    ! Bit representation of x

  ! Range-reduction variables

  real :: xc
    ! x clamped between xmin and xmax [nondim]
  integer(kind=int_kind) :: K
    ! Nearest IEEE-rounded integer to x/ln2 [nondim]
  real :: Z
    ! Nearest integer to ndiv * K + idiv for ndiv subdivisions
    ! NOTE: Z is stored as real to avoid int-real type conversions.
  integer(kind=int_kind) :: Zi
    ! Integer representation of Z
  integer :: idiv
    ! Subdivision index
  real :: r
    ! Range-reduced input, r = x - Z ln2 / ndiv [nondim]

  ! Polynomial estimation variables

  real :: e
    ! Exponent of range-reduced input, e = 2**(idiv/ndiv) exp(r) [nondim]
  real :: expm1_r
    ! Approximation to exp(r) - 1 [nondim]
  real :: idiv_scale
    ! Estimate of 2**(idiv/ndiv) [nondim]
  real :: idiv_residual
    ! Relative residual, 2**(j/N) = idiv_scale * (1 + idiv_residual) [nondim]

  ! Descaling and subnormal handling

  integer(kind=int_kind) :: eb
    ! Bit representations of e
  integer(kind=int_kind) :: j
    ! Bias added to K to compensate for exponent K beyond {-1022,..,+1023}.
  integer(kind=int_kind) :: fb
    ! Bit representation 2**j, the K exponent rescale

  ! 1. Nonfinite handling
  ! ---------------------
  ! Nonfinites must be handled first to prevent their appearance in
  ! calculations, which may raise unwanted floating point signals.

  xb = transfer(x, int_mold)
  nonfinite = iand(xb, pos_inf_bits) == pos_inf_bits

  if (nonfinite) then
    ! exp(-Inf) = 0, otherwise pass-through +Inf and +/-NaN values
    ! Compute x + x to trigger `Invalid` for signaled NaNs.
    a = merge(0., x + x, xb == neg_inf_bits)
    return
  endif

  ! 2. Range Reduction
  ! ------------------
  ! Apply a range reduction of r = x - K ln2 - (idiv / ndiv) ln2, so that
  !     exp(x) = 2**K 2**(idiv / ndiv) exp(r).
  ! If K = nint(x / ln2) then r is in [-ln2/(2*ndiv), ln2/(2*ndiv)] and exp(r)
  ! can be estimated by a sufficiently accurate polynomial.

  ! Clamp x to [xmin,xmax] to avoid extreme exponents in subnormal handler.
  xc = min(max(x, xmin), xmax)

  ! Compute Z = ndiv K + idiv, where r = x - Z (ln 2 / ndiv)
  Z = NEAREST_INT(xc * n_ln2)
  Zi = transfer(Z + round_bias, int_mold) - round_bias_bits

  ! Compute the subdivision index and integer offset K
  idiv = iand(Zi, idiv_mask)
  K = (Zi - int(idiv, int_kind)) / ndiv

  ! Since Z ~ x N / ln2, the terms in r will nearly cancel and there is some
  ! expected loss of precision.  To compensate, we use a Cody-Waite correction.
  r = (xc - Z * ln2_ndiv_hi) - Z * ln2_ndiv_lo

  ! 3. Polynomial approximation
  ! ---------------------------
  ! a = 2**K e where e = 2**(idiv/ndiv) exp(r), and exp(r) = 1 + r * P(r) where
  ! P(r) is an order-5 Remez minimax polynomial of expm1(r) / r.

  idiv_scale = idiv_scale_lookup(idiv)
  idiv_residual = idiv_residual_lookup(idiv)
  expm1_r = exp_remez_expm1_estrin_4(r)

  ! Evaluate the small correction before the final addition to idiv_scale
  e = idiv_scale + idiv_scale * (idiv_residual + r * expm1_r)

  ! 4. Unscaling
  ! ------------
  ! Compute a = 2**K e, an exact power-of-2 calculation.
  ! Adjust scaling to compensate for subnormal output.

  ! exp(r) has range [0.707/ndiv, 1.414/ndiv], so K shifts by either 0 or -1.
  ! Resolved exponent are in the range {-1022..1023}, so for a to be resolved,
  ! K must be in {-1021,1023}.  (For safety, we actually do {-1020,1022}.)

  ! Determine if K is outside the supported exponent range.
  ! If so, then apply a bias j to normalize the exponent.
  ! Kbias is chosen so that the exponent is "something near 1".
  j = merge(Kbias, 0_int_kind, K < Kmin) + merge(-Kbias, 0_int_kind, K > Kmax)

  ! Get the bit representation of e
  eb = transfer(e, int_mold)

  ! Rescale to e to exp(x), possibly including the j bias.
  eb = eb + ishft(k + j, expbit_exp)
  a = transfer(eb, real_mold)

  ! Undo the 2**j bias as floating point multiplication.
  ! - For "normals", this has no effect.
  ! - For subnormals, this will force subnormal estimation (if enabled).
  ! - For resolvable K beyond this range, it triggers an over/underflow.
  ! - Extreme values of K have already been filtered out by the min/max step.
  fb = ishft(int(expbias_exp, int_kind) - j, expbit_exp)
  a = a * transfer(fb, real_mold)
end function exp_repro


!> Remez polynomial estimate of (exp(x) - 1) / x over [-ln2/256, ln2/256].
!! Coefficients are generated by Sollya 8, and evaluation is in Estrin form.
!! Coefficients (exp_remez_c) are module-scope, not function-local: a
!! function-local parameter array is what blocked -Mextract/-Minline
!! cross-file inlining in experiment E4 (see exp_remez_c's declaration above).
pure function exp_remez_expm1_estrin_4(x) result(e)
  !$omp declare target
  real, intent(in) :: x
    !< Input value; expected range is [-ln2/256, ln2/256] [nondim]
  real :: e
    !< Approximation of (exp(x) - 1) / x [nondim]

  real :: x2, x4
    !< Powers of x [nondim]
  real :: p01, p23
    !< Polynomial partial sums [nondim]

  x2 = x * x
  x4 = x2 * x2

  p01 = exp_remez_c(0) + exp_remez_c(1) * x
  p23 = exp_remez_c(2) + exp_remez_c(3) * x

  ! Final assembly: (exp(x) - 1) / x
  e = (p01 + x2 * p23) + x4 * exp_remez_c(4)
end function exp_remez_expm1_estrin_4


!> An optimized nearest-integer function for floating point reals.
!!
!! The value x is shifted from 2**K (1 + a) to 2**(digits-1+K) (1.5 + [a]),
!! causing the fractional part to be rounded according to the current IEEE
!! settings.  In almost all cases, this is nearest ties-to-even.
!!
!! The +0.5 ensures that the biased exponent of negative numbers does not drop
!! by one, which can cause half-value rounding.
!!
!! The behavior of this function does not match nint() or anint().  The nint()
!! function always ties away from zero, e.g. nint(2.5) = 3.
!!
!! It is essential that compilers not reduce (x+b)-b to x.  This can typically
!! be ensured as long as parentheses are respected.  This is managed by the
!! ENABLE_FAST_RINT macro in MOM_exp.h and assigned to NEAREST_INT().  If
!! unset, then ieee_rint() is used.
pure function fast_rint(x) result(n)
  !$omp declare target
  real, intent(in) :: x
    !< Real value to be rounded to the nearest integer
  real :: n
    !< Nearest integer to x, stored as a real

  n = (x + round_bias) - round_bias
end function fast_rint


!> Returns true if any unit test of intrinsic_functions fails, or false if they all pass.
function intrinsic_functions_unit_tests(verbose) result(fail)
  logical, intent(in) :: verbose !< If true, write results to stdout
  logical :: fail !< True if any of the unit tests fail

  ! Local variables
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
end function intrinsic_functions_unit_tests

!> True if the cube of cuberoot(val) does not closely match val. False otherwise.
logical function Test_cuberoot(verbose, val)
  logical, intent(in) :: verbose !< If true, write results to stdout
  real, intent(in) :: val  !< The real value to test, in arbitrary units [A]
  ! Local variables
  real :: diff ! The difference between val and the cube root of its cube [A].

  diff = val - cuberoot(val)**3
  Test_cuberoot = (abs(diff) > 2.0e-15*abs(val))

  if (Test_cuberoot) then
    write(stdout, '("For val = ",ES22.15,", (val - cuberoot(val**3))) = ",ES9.2," <-- FAIL")') val, diff
  elseif (verbose) then
    write(stdout, '("For val = ",ES22.15,", (val - cuberoot(val**3))) = ",ES9.2)') val, diff

  endif
end function Test_cuberoot

end module MOM_intrinsic_functions
