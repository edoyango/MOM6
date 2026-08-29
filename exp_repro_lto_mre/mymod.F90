module mymod
  implicit none
contains

  !> Top-level function, meant to be cross-inlined into a foreign TU.
  !! Calls a second, private helper function internally -- mirrors
  !! exp_repro() calling exp_remez_expm1_estrin_4().
  elemental function foo(x) result(y)
    !$omp declare target
    real, intent(in) :: x
    real :: y
    y = 2.0 * helper(x) + 1.0
  end function foo

  !> Private helper, no explicit declare target -- mirrors
  !! exp_remez_expm1_estrin_4() in the real code.
  pure function helper(x) result(h)
    real, intent(in) :: x
    real :: h
    h = x * x + x
  end function helper

end module mymod
