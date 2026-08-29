! Minimal reproducer for exp_repro_gpu_perf_investigation.md, section 8.4.
!
! Two-file OpenMP-target program with a genuine cross-TU device call
! (mymod.F90's `foo`, called here in driver2.F90), structurally identical
! to exp_repro() calling exp_remez_expm1_estrin_4() in the real
! MOM_set_viscosity.F90 investigation: an elemental function carrying
! `!$omp declare target` that itself calls a second, private helper with
! no explicit `declare target`.
!
! Confirms -gpu=lto,rdc compiles, links, and runs correctly at this scale --
! the LTO link failure documented in section 8.4 is scale-dependent, not a
! blanket incompatibility between LTO, RDC, and OpenMP-target-GPU.
!
! Build and run (nvhpc 26.3 or 26.5, both confirmed):
!   nvfortran -O4 -mp=gpu -gpu=lto,rdc -c mymod.F90 -o mymod_lto.o
!   nvfortran -O4 -mp=gpu -gpu=lto,rdc -c driver2.F90 -o driver2_lto.o
!   nvfortran -mp=gpu -gpu=lto,rdc -o lto_test mymod_lto.o driver2_lto.o
!   ./lto_test
! Expected output: " sum:    686900.0"

program driver2
  use mymod, only : foo
  implicit none
  integer, parameter :: npts = 100
  real :: x(npts), val(npts)
  integer :: i
  do i = 1, npts
    x(i) = real(i)
  enddo
  !$omp target teams distribute parallel do map(to:x) map(from:val)
  do i = 1, npts
    val(i) = foo(x(i))
  enddo
  !$omp end target teams distribute parallel do
  print *, "sum:", sum(val)
end program driver2
