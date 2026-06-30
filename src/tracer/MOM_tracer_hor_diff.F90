! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> Main routine for lateral (along surface or neutral) diffusion of tracers
module MOM_tracer_hor_diff

use MOM_cpu_clock,                only : cpu_clock_id, cpu_clock_begin, cpu_clock_end
use MOM_cpu_clock,                only : CLOCK_MODULE, CLOCK_ROUTINE
use MOM_diag_mediator,            only : post_data, diag_ctrl
use MOM_diag_mediator,            only : register_diag_field, safe_alloc_ptr, time_type
use MOM_domains,                  only : sum_across_PEs, max_across_PEs
use MOM_domains,                  only : create_group_pass, do_group_pass, group_pass_type
use MOM_domains,                  only : pass_vector
use MOM_debugging,                only : hchksum, uvchksum
use MOM_diabatic_driver,          only : diabatic_CS
use MOM_EOS,                      only : calculate_density, EOS_type, EOS_domain
use MOM_error_handler,            only : MOM_error, FATAL, WARNING, MOM_mesg, is_root_pe
use MOM_error_handler,            only : MOM_set_verbosity, callTree_showQuery
use MOM_error_handler,            only : callTree_enter, callTree_leave, callTree_waypoint
use MOM_file_parser,              only : get_param, log_version, param_file_type
use MOM_grid,                     only : ocean_grid_type
use MOM_lateral_mixing_coeffs,    only : VarMix_CS
use MOM_MEKE_types,               only : MEKE_type
use MOM_neutral_diffusion,        only : neutral_diffusion_init, neutral_diffusion_end
use MOM_neutral_diffusion,        only : neutral_diffusion_CS
use MOM_neutral_diffusion,        only : neutral_diffusion_calc_coeffs, neutral_diffusion
use MOM_hor_bnd_diffusion,        only : hbd_CS, hor_bnd_diffusion_init
use MOM_hor_bnd_diffusion,        only : hor_bnd_diffusion, hor_bnd_diffusion_end
use MOM_tracer_registry,          only : tracer_registry_type, tracer_type, MOM_tracer_chksum
use MOM_unit_scaling,             only : unit_scale_type
use MOM_variables,                only : thermo_var_ptrs, vertvisc_type
use MOM_verticalGrid,             only : verticalGrid_type

implicit none ; private

#include <MOM_memory.h>

public tracer_hordiff, tracer_hor_diff_init, tracer_hor_diff_end

!> The control structure for along-layer and epineutral tracer diffusion
type, public :: tracer_hor_diff_CS ; private
  real    :: KhTr           !< The along-isopycnal tracer diffusivity [L2 T-1 ~> m2 s-1].
  real    :: KhTr_Slope_Cff !< The non-dimensional coefficient in KhTr formula [nondim]
  real    :: KhTr_min       !< Minimum along-isopycnal tracer diffusivity [L2 T-1 ~> m2 s-1].
  real    :: KhTr_max       !< Maximum along-isopycnal tracer diffusivity [L2 T-1 ~> m2 s-1].
  real    :: KhTr_passivity_coeff !< Passivity coefficient that scales Rd/dx (default = 0)
                                  !! where passivity is the ratio between along-isopycnal
                                  !! tracer mixing and thickness mixing [nondim]
  real    :: KhTr_passivity_min   !< Passivity minimum (default = 1/2) [nondim]
  real    :: ML_KhTR_scale        !< With Diffuse_ML_interior, the ratio of the
                                  !! truly horizontal diffusivity in the mixed
                                  !! layer to the epipycnal diffusivity [nondim].
  real    :: max_diff_CFL         !< If positive, locally limit the along-isopycnal
                                  !! tracer diffusivity to keep the diffusive CFL
                                  !! locally at or below this value [nondim].
  logical :: KhTr_use_vert_struct  !< If true, uses the equivalent barotropic structure
                                  !! as the vertical structure of tracer diffusivity.
  logical :: full_depth_khtr_min  !< If true, KHTR_MIN is enforced throughout the whole water column.
                                  !! Otherwise, KHTR_MIN is only enforced at the surface. This parameter
                                  !! is only available when KHTR_USE_EBT_STRUCT=True and KHTR_MIN>0.
  logical :: Diffuse_ML_interior  !< If true, diffuse along isopycnals between
                                  !! the mixed layer and the interior.
  logical :: check_diffusive_CFL  !< If true, automatically iterate the diffusion
                                  !! to ensure that the diffusive equivalent of
                                  !! the CFL limit is not violated.
  logical :: use_neutral_diffusion !< If true, use the neutral_diffusion module from within
                                   !! tracer_hor_diff.
  logical :: use_hor_bnd_diffusion !< If true, use the hor_bnd_diffusion module from within
                                   !! tracer_hor_diff.
  logical :: recalc_neutral_surf   !< If true, recalculate the neutral surfaces if CFL has been
                                   !! exceeded
  logical :: limit_bug             !< If true and the answer date is 20240330 or below, use a
                                   !! rotational symmetry breaking bug when limiting the tracer
                                   !! properties in tracer_epipycnal_ML_diff.
  integer :: answer_date           !< The vintage of the order of arithmetic to use for the tracer
                                   !! diffusion.  Values of 20240330 or below recover the answers
                                   !! from the original form of this code, while higher values use
                                   !! mathematically equivalent expressions that recover rotational symmetry
                                   !! when DIFFUSE_ML_TO_INTERIOR is true.
  type(neutral_diffusion_CS), pointer :: neutral_diffusion_CSp => NULL() !< Control structure for neutral diffusion.
  type(hbd_CS), pointer    :: hor_bnd_diffusion_CSp => NULL() !< Control structure for
                                                              !! horizontal boundary diffusion.
  type(diag_ctrl), pointer :: diag => NULL() !< A structure that is used to
                                   !! regulate the timing of diagnostic output.
  logical :: debug                 !< If true, write verbose checksums for debugging purposes.
  logical :: show_call_tree        !< Display the call tree while running. Set by VERBOSITY level.
  logical :: first_call = .true.   !< This is true until after the first call
  !>@{ Diagnostic IDs
  integer :: id_KhTr_u  = -1
  integer :: id_KhTr_v  = -1
  integer :: id_KhTr_h  = -1
  integer :: id_CFL     = -1
  integer :: id_khdt_x  = -1
  integer :: id_khdt_y  = -1
  !>@}

  type(group_pass_type) :: pass_t !< For group halo pass, used in both
                                  !! tracer_hordiff and tracer_epipycnal_ML_diff
end type tracer_hor_diff_CS

!> A type that can be used to create arrays of pointers to 2D arrays
type p2d
  real, dimension(:,:), pointer :: p => NULL() !< A pointer to a 2D array of reals [various]
end type p2d
!> A type that can be used to create arrays of pointers to 2D integer arrays
type p2di
  integer, dimension(:,:), pointer :: p => NULL() !< A pointer to a 2D array of integers
end type p2di

!>@{ CPU time clocks
integer :: id_clock_diffuse, id_clock_epimix, id_clock_pass, id_clock_sync
!>@}


  interface
module subroutine tracer_hordiff(h, dt, MEKE, VarMix, visc, G, GV, US, CS, Reg, tv, do_online_flag, read_khdt_x, read_khdt_y)
  type(ocean_grid_type),      intent(inout) :: G       !< Grid type
  type(verticalGrid_type),    intent(in)    :: GV      !< ocean vertical grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                              intent(in)    :: h       !< Layer thickness [H ~> m or kg m-2]
  real,                       intent(in)    :: dt      !< time step [T ~> s]
  type(MEKE_type),            intent(in)    :: MEKE    !< MEKE fields
  type(VarMix_CS),            intent(in)    :: VarMix  !< Variable mixing type
  type(vertvisc_type),        intent(in)    :: visc    !< Structure with vertical viscosities,
                                                       !! boundary layer properties and related fields
  type(unit_scale_type),      intent(in)    :: US      !< A dimensional unit scaling type
  type(tracer_hor_diff_CS),   pointer       :: CS      !< module control structure
  type(tracer_registry_type), pointer       :: Reg     !< registered tracers
  type(thermo_var_ptrs),      intent(in)    :: tv      !< A structure containing pointers to any available
                                                       !! thermodynamic fields, including potential temp and
                                                       !! salinity or mixed layer density. Absent fields have
                                                       !! NULL ptrs, and these may (probably will) point to
                                                       !! some of the same arrays as Tr does.  tv is required
                                                       !! for epipycnal mixing between mixed layer and the interior.
  ! Optional inputs for offline tracer transport
  logical,          optional, intent(in)    :: do_online_flag !< If present and true, do online
                                                       !! tracer transport with stored velocities.
  ! The next two arguments do not appear to be used anywhere.
  real, dimension(SZIB_(G),SZJ_(G)), &
                    optional, intent(in)    :: read_khdt_x !< If present, these are the zonal diffusivities
                                                       !! times a timestep from a previous run [L2 ~> m2]
  real, dimension(SZI_(G),SZJB_(G)), &
                    optional, intent(in)    :: read_khdt_y !< If present, these are the meridional diffusivities
                                                       !! times a timestep from a previous run [L2 ~> m2]


end subroutine tracer_hordiff
module subroutine tracer_epipycnal_ML_diff(h, dt, Tr, ntr, khdt_epi_x, khdt_epi_y, G, &
                                    GV, US, CS, tv, num_itts)
  type(ocean_grid_type),                    intent(inout) :: G          !< ocean grid structure
  type(verticalGrid_type),                  intent(in)    :: GV         !< ocean vertical grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in)   :: h          !< layer thickness [H ~> m or kg m-2]
  real,                                     intent(in)    :: dt         !< time step [T ~> s]
  type(tracer_type),                        intent(inout) :: Tr(:)      !< tracer array
  integer,                                  intent(in)    :: ntr        !< number of tracers
  real, dimension(SZIB_(G),SZJ_(G)),        intent(in)    :: khdt_epi_x !< Zonal epipycnal diffusivity times
                                                           !! a time step and the ratio of the open face width over
                                                           !! the distance between adjacent tracer points [L2 ~> m2]
  real, dimension(SZI_(G),SZJB_(G)),        intent(in)    :: khdt_epi_y !< Meridional epipycnal diffusivity times
                                                           !! a time step and the ratio of the open face width over
                                                           !! the distance between adjacent tracer points [L2 ~> m2]
  type(unit_scale_type),                    intent(in)    :: US         !< A dimensional unit scaling type
  type(tracer_hor_diff_CS),                 intent(inout) :: CS         !< module control structure
  type(thermo_var_ptrs),                    intent(in)    :: tv         !< thermodynamic structure
  integer,                                  intent(in)    :: num_itts   !< number of iterations (usually=1)



  ! The naming mnemonic is a=above,b=below,L=Left,R=Right,u=u-point,v=v-point.
  ! These are 1-D arrays of pointers to 2-d arrays to minimize memory usage.




  ! The following 3-d arrays were created in 2014 in MOM6 PR#12 to facilitate openMP threading
  ! on an i-loop, which might have been ill advised.
end subroutine tracer_epipycnal_ML_diff
module subroutine tracer_hor_diff_init(Time, G, GV, US, param_file, diag, EOS, diabatic_CSp, CS)
  type(time_type), target,    intent(in)    :: Time       !< current model time
  type(ocean_grid_type),      intent(in)    :: G          !< ocean grid structure
  type(verticalGrid_type),    intent(in)    :: GV         !< ocean vertical grid structure
  type(unit_scale_type),      intent(in)    :: US         !< A dimensional unit scaling type
  type(diag_ctrl), target,    intent(inout) :: diag       !< diagnostic control
  type(EOS_type),  target,    intent(in)    :: EOS        !< Equation of state CS
  type(diabatic_CS), pointer, intent(in)    :: diabatic_CSp !< Equation of state CS
  type(param_file_type),      intent(in)    :: param_file !< parameter file
  type(tracer_hor_diff_CS),   pointer       :: CS         !< horz diffusion control structure

  ! This include declares and sets the variable "version".

end subroutine tracer_hor_diff_init
module subroutine tracer_hor_diff_end(CS)
  type(tracer_hor_diff_CS), pointer :: CS !< module control structure

end subroutine tracer_hor_diff_end
  end interface

end module MOM_tracer_hor_diff
