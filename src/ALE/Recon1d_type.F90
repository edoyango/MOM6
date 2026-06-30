! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> A generic type for vertical 1D reconstructions
module Recon1d_type

use numerical_testing_type, only : testing

implicit none ; private

public Recon1d
public testing

!> The base class for implementations of 1D reconstructions
type, abstract :: Recon1d

  integer :: n = 0 !< Number of cells in column
  real, allocatable, dimension(:) :: u_mean !< Cell mean [A]
  real :: h_neglect = 0. !< A negligibly small width used in cell reconstructions in the same units as h [H]
  real :: x_tolerance = 1. * epsilon(1.) !< Solver tolerance for x in element (0,1) [nondim]
  logical :: check = .false. !< If true, enable some consistency checking

  logical :: debug = .false. !< If true, dump info as calculations are made (do not enable)
contains

  ! The following functions/subroutines are deferred and must be provided specifically by each scheme

  !> Deferred implementation of initialization
  procedure(i_init), deferred :: init
  !> Deferred implementation of reconstruction function
  procedure(i_reconstruct), deferred :: reconstruct
  !> Deferred implementation of the average over an interval
  procedure(i_average), deferred :: average
  !> Deferred implementation of evaluating the reconstruction at a point
  procedure(i_f), deferred :: f
  !> Deferred implementation of the derivative of the reconstruction at a point
  procedure(i_dfdx), deferred :: dfdx
  !> Deferred implementation of check_reconstruction
  !!
  !! Returns True if a check fails. Returns False if all checks pass.
  !! Checks are about internal, or inferred, state for arbitrary inputs.
  !! Checks should cover all the expected properties of a reconstruction.
  procedure(i_check_reconstruction), deferred :: check_reconstruction
  !> Deferred implementation of unit tests for the reconstruction
  !!
  !! Returns True if a test fails. Returns False if all tests pass.
  !! Tests in unit_tests() are usually checks against known (e.g. analytic) solutions.
  procedure(i_unit_tests), deferred :: unit_tests
  !> Deferred implementation of deallocation
  procedure(i_destroy), deferred :: destroy

  ! The following functions/subroutines are shared across all reconstructions and provided by this module
  ! unless replaced for the purpose of optimization

  !> Solves for x such that f(x)=t
  procedure :: x => x
  !> Remaps the column to subgrid h_sub
  procedure :: remap_to_sub_grid => remap_to_sub_grid
  !> Set debugging
  procedure :: set_debug => a_set_debug

  ! The following functions usually point to the same implementation as above but
  ! for derived secondary children these allow invocation of the parent class function.

  !> Second interface to init(), used to reach the primary class if derived from a primary implementation
  procedure(i_init_parent), deferred :: init_parent
  !> Second interface to reconstruct(), used to reach the primary class if derived from a primary implementation
  procedure(i_reconstruct_parent), deferred :: reconstruct_parent

end type Recon1d

interface

  !> Initialize a 1D reconstruction for n cells
  subroutine i_init(this, n, h_neglect, check)
    import :: Recon1d
    class(Recon1d),    intent(out) :: this !< This reconstruction
    integer,           intent(in)  :: n    !< Number of cells in this column
    real, optional,    intent(in)  :: h_neglect !< A negligibly small width used in cell reconstructions [H]
    logical, optional, intent(in)  :: check !< If true, enable some consistency checking
  end subroutine i_init

  !> Calculate a 1D reconstructions based on h(:) and u(:)
  subroutine i_reconstruct(this, h, u)
    import :: Recon1d
    class(Recon1d), intent(inout) :: this !< This reconstruction
    real,           intent(in)    :: h(*) !< Grid spacing (thickness), typically in [H]
    real,           intent(in)    :: u(*) !< Cell mean values [A]
  end subroutine i_reconstruct

  !> Average between xa and xb for cell k of a 1D reconstruction [A]
  !!
  !! It is assumed that 0<=xa<=1, 0<=xb<=1, and xa<=xb
  real function i_average(this, k, xa, xb)
    import :: Recon1d
    class(Recon1d), intent(in) :: this !< This reconstruction
    integer,        intent(in) :: k    !< Cell number
    real,           intent(in) :: xa   !< Start of averaging interval on element (0 to 1)
    real,           intent(in) :: xb   !< End of averaging interval on element (0 to 1)
  end function i_average

  !> Point-wise value of reconstruction [A]
  !!
  !! The function is only valid for 0 <= x <= 1. x is effectively clipped to this range.
  real function i_f(this, k, x)
    import :: Recon1d
    class(Recon1d), intent(in) :: this !< This reconstruction
    integer,        intent(in) :: k    !< Cell number
    real,           intent(in) :: x    !< Non-dimensional position within element [nondim]
  end function i_f

  !> Point-wise value of derivative reconstruction [A]
  !!
  !! The function is only valid for 0 <= x <= 1. x is effectively clipped to this range.
  real function i_dfdx(this, k, x)
    import :: Recon1d
    class(Recon1d), intent(in) :: this !< This reconstruction
    integer,        intent(in) :: k    !< Cell number
    real,           intent(in) :: x    !< Non-dimensional position within element [nondim]
  end function i_dfdx

  !> Point-wise solver for x: f(x)=t [nondim]
  !!
  !! The function solves for the non-dimensional position x within the cell where
  !! the reconstruction f(x)=t. The solver returns x=0 or x=1 if the target, t,
  !! is outside of the cell.
  real function i_x(this, k, t)
    import :: Recon1d
    class(Recon1d), intent(in) :: this !< This reconstruction
    integer,        intent(in) :: k    !< Cell number
    real,           intent(in) :: t    !< Value to solve for [A]
  end function i_x

  !> Returns true if some inconsistency is detected, false otherwise
  !!
  !! The nature of "consistency" is defined by the implementations
  !! and might be no-ops.
  logical function i_check_reconstruction(this, h, u)
    import :: Recon1d
    class(Recon1d), intent(in) :: this !< This reconstruction
    real,           intent(in) :: h(*) !< Grid spacing (thickness), typically in [H]
    real,           intent(in) :: u(*) !< Cell mean values [A]
  end function i_check_reconstruction

  !> Deallocate a 1D reconstruction
  subroutine i_destroy(this)
    import :: Recon1d
    class(Recon1d), intent(inout) :: this !< This reconstruction
  end subroutine i_destroy

  !> Second interface to init(), or to parent init()
  subroutine i_init_parent(this, n, h_neglect, check)
    import :: Recon1d
    class(Recon1d), intent(out) :: this !< This reconstruction
    integer,        intent(in)  :: n    !< Number of cells in this column
    real, optional, intent(in)  :: h_neglect !< A negligibly small width used in cell reconstructions [H]
    logical, optional, intent(in)  :: check !< If true, enable some consistency checking
  end subroutine i_init_parent

  !> Second interface to reconstruct(), or to parent reconstruct()
  subroutine i_reconstruct_parent(this, h, u)
    import :: Recon1d
    class(Recon1d), intent(inout) :: this !< This reconstruction
    real,           intent(in)    :: h(*) !< Grid spacing (thickness), typically in [H]
    real,           intent(in)    :: u(*) !< Cell mean values [A]
  end subroutine i_reconstruct_parent

  !> Runs reconstruction unit tests and returns True for any fails, False otherwise
  !!
  !! Assumes single process/thread context
  logical function i_unit_tests(this, verbose, stdout, stderr)
    import :: Recon1d
    class(Recon1d), intent(inout) :: this    !< This reconstruction
    logical,        intent(in)    :: verbose !< True, if verbose
    integer,        intent(in)    :: stdout  !< I/O channel for stdout
    integer,        intent(in)    :: stderr  !< I/O channel for stderr
  end function i_unit_tests

end interface


  interface
real module function x(this, k, t)
  class(Recon1d), intent(in) :: this !< This reconstruction
  integer,        intent(in) :: k    !< Cell number
  real,           intent(in) :: t    !< Value to solve for [A]

end function x
module subroutine remap_to_sub_grid(this, h0, u0, n1, h_sub, &
                                   isrc_start, isrc_end, isrc_max, isub_src, &
                                   u_sub, uh_sub, u02_err)
  class(Recon1d), intent(in) :: this !< 1-D reconstruction type
  real,    intent(in)  :: h0(*)  !< Source grid widths (size n0) [H]
  real,    intent(in)  :: u0(*)  !< Source grid widths (size n0) [H]
  integer, intent(in)  :: n1      !< Number of cells in target grid
  real,    intent(in)  :: h_sub(*) !< Overlapping sub-cell thicknesses, h_sub [H]
  integer, intent(in)  :: isrc_start(*) !< Index of first sub-cell within each source cell
  integer, intent(in)  :: isrc_end(*) !< Index of last sub-cell within each source cell
  integer, intent(in)  :: isrc_max(*) !< Index of thickest sub-cell within each source cell
  integer, intent(in)  :: isub_src(*) !< Index of source cell for each sub-cell
  real,    intent(out) :: u_sub(*) !< Sub-cell cell averages (size n1) [A]
  real,    intent(out) :: uh_sub(*) !< Sub-cell cell integrals (size n1) [A H]
  real,    intent(out) :: u02_err !< Integrated reconstruction error estimates [A H]
  ! Local variables
! real :: u0_min(this%n), u0_max(this%n) ! Min/max of u0 for each source cell [A]
! real :: ul,ur ! Left/right edge values [A]

end subroutine remap_to_sub_grid
module subroutine a_set_debug(this)
  class(Recon1d), intent(inout) :: this !< 1-D reconstruction type

end subroutine a_set_debug
  end interface

end module Recon1d_type
