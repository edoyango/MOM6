! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> This module contains the hybgen remapping routines from HYCOM, with minor
!! modifications to follow the MOM6 coding conventions
module MOM_hybgen_remap

implicit none ; private

public hybgen_plm_coefs, hybgen_ppm_coefs, hybgen_weno_coefs


  interface
module subroutine hybgen_plm_coefs(si, dpi, slope, nk, ns, thin, PCM_lay)
  integer, intent(in)  :: nk        !< The number of input layers
  integer, intent(in)  :: ns        !< The number of scalar fields to work on
  real,    intent(in)  :: si(nk,ns) !< The cell-averaged input scalar fields [A]
  real,    intent(in)  :: dpi(nk)   !< The input grid layer thicknesses [H ~> m or kg m-2]
  real,    intent(out) :: slope(nk,ns) !< The PLM slope times cell width [A]
  real,    intent(in)  :: thin      !< A negligible layer thickness that can be ignored [H ~> m or kg m-2]
  logical, optional, intent(in)  :: PCM_lay(nk) !< If true for a layer, use PCM remapping for that layer

!-----------------------------------------------------------------------
!  1) coefficients for remapping from one set of vertical cells to another.
!     method: piecewise linear across each input cell with
!             monotonized central-difference limiter.
!
!     van Leer, B., 1977, J. Comp. Phys., 23 276-299.
!
!  2) input arguments:
!       si    - initial scalar fields in pi-layer space
!       dpi   - initial layer thicknesses (dpi(k) = pi(k+1)-pi(k))
!       nk    - number of layers
!       ns    - number of fields
!       thin  - layer thickness (>0) that can be ignored
!       PCM_lay - use PCM for selected layers (optional)
!
!  3) output arguments:
!       slope - coefficients for hybgen_plm_remap
!                profile(y) = si+slope*(y-1),  -0.5 <= y <= 0.5
!
!  4) Tim Campbell, Mississippi State University, October 2002.
!     Alan J. Wallcraft,  Naval Research Laboratory,  Aug. 2007.
!-----------------------------------------------------------------------
!
                 ! of the adjacent cells, usually ~0.5, but always <= 1 [nondim]

end subroutine hybgen_plm_coefs
module subroutine hybgen_ppm_coefs(s, h_src, edges, nk, ns, thin, PCM_lay)
  integer, intent(in)  :: nk        !< The number of input layers
  integer, intent(in)  :: ns        !< The scalar fields to work on
  real,    intent(in)  :: s(nk,ns)  !< The input scalar fields [A]
  real,    intent(in)  :: h_src(nk) !< The input grid layer thicknesses [H ~> m or kg m-2]
  real,    intent(out) :: edges(nk,2,ns) !< The PPM interpolation edge values of the scalar fields [A]
  real,    intent(in)  :: thin      !< A negligible layer thickness that can be ignored [H ~> m or kg m-2]
  logical, optional, intent(in)  :: PCM_lay(nk) !< If true for a layer, use PCM remapping for that layer

!-----------------------------------------------------------------------
!  1) coefficients for remapping from one set of vertical cells to another.
!     method: monotonic piecewise parabolic across each input cell
!
!     Colella, P. & P.R. Woodward, 1984, J. Comp. Phys., 54, 174-201.
!
!  2) input arguments:
!       s     - initial scalar fields in pi-layer space
!       h_src - initial layer thicknesses (>=0)
!       nk    - number of layers
!       ns    - number of fields
!       thin  - layer thickness (>0) that can be ignored
!       PCM_lay - use PCM for selected layers (optional)
!
!  3) output arguments:
!       edges - cell edge scalar values for the PPM reconstruction
!                edges.1 is value at interface above
!                edges.2 is value at interface below
!
!  4) Tim Campbell, Mississippi State University, October 2002.
!     Alan J. Wallcraft,  Naval Research Laboratory,  Aug. 2007.
!-----------------------------------------------------------------------
!
                           ! very thin, or because this is specified by PCM_lay.

  ! This PPM remapper is not currently written to work with massless layers, so set
  ! the thicknesses for very thin layers to some minimum value.
end subroutine hybgen_ppm_coefs
module subroutine hybgen_weno_coefs(s, h_src, edges, nk, ns, thin, PCM_lay)
  integer, intent(in)  :: nk        !< The number of input layers
  integer, intent(in)  :: ns        !< The number of scalar fields to work on
  real,    intent(in)  :: s(nk,ns)  !< The input scalar fields [A]
  real,    intent(in)  :: h_src(nk) !< The input grid layer thicknesses [H ~> m or kg m-2]
  real,    intent(out) :: edges(nk,2,ns) !< The WENO interpolation edge values of the scalar fields [A]
  real,    intent(in)  :: thin      !< A negligible layer thickness that can be ignored [H ~> m or kg m-2]
  logical, optional, intent(in)  :: PCM_lay(nk) !< If true for a layer, use PCM remapping for that layer

!-----------------------------------------------------------------------
!  1) coefficients for remapping from one set of vertical cells to another.
!     method: monotonic WENO-like alternative to PPM across each input cell
!             a second order polynomial approximation of the profiles
!             using a WENO reconciliation of the slopes to compute the
!             interfacial values
!
!     This scheme might have ben developed by Shchepetkin. A.F., personal communication.
!     See also Engwirda, D., and M. Kelley, A WENO-type slope-limiter for a family of piecewise
!       polynomial methods, arXive:1606.08188v1, 27 June 2016.
!
!  2) input arguments:
!       s     - initial scalar fields in pi-layer space
!       h_src - initial layer thicknesses (>=0)
!       nk    - number of layers
!       ns    - number of fields
!       thin  - layer thickness (>0) that can be ignored
!       PCM_lay - use PCM for selected layers (optional)
!
!  3) output arguments:
!       edges - cell edge scalar values for the WENO reconstruction
!                edges.1 is value at interface above
!                edges.2 is value at interface below
!
!  4) Laurent Debreu, Grenoble.
!     Alan J. Wallcraft,  Naval Research Laboratory,  July 2008.
!-----------------------------------------------------------------------
!
!  real, parameter :: dsmll=1.0e-8  ! This has units of [A2], and hence can not be a parameter.
!
                      ! spacing [A H-1 ~> A m-1 or A m2 kg-1]
                      ! very thin, or because this is specified by PCM_lay.
                      ! concentrations and the left and right edges [A2]

end subroutine hybgen_weno_coefs
  end interface

end module MOM_hybgen_remap
