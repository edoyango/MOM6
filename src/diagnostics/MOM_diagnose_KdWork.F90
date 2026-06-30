! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> Provides diagnostics of work due to a given diffusivity
module MOM_diagnose_kdwork

use MOM_diag_mediator, only : diag_ctrl, time_type, post_data, register_diag_field
use MOM_diag_mediator, only : register_scalar_field
use MOM_error_handler, only : MOM_error, FATAL, WARNING
use MOM_grid,          only : ocean_grid_type
use MOM_unit_scaling,  only : unit_scale_type
use MOM_variables,     only : thermo_var_ptrs
use MOM_verticalGrid,  only : verticalGrid_type
use MOM_spatial_means, only : global_area_integral

implicit none ; private

#include <MOM_memory.h>

public vbf_CS
public kdwork_diagnostics
public Allocate_VBF_CS
public Deallocate_VBF_CS
public KdWork_init
public KdWork_end

!> This structure has memory for used in calculating diagnostics of diffusivity
!! many of the diffusivity diagnostics are copies of other 3d arrays.  It could
!! be written more efficiently, but it is less intrusive to copy into this structure
!! and do all calculations in this module.  These diagnostics may be expensive for
!! routine use.
type vbf_CS
  ! 3d varying Kd contributions
  real, pointer, dimension(:,:,:) :: &
    Bflx_salt => NULL(), & !< Salinity contribution to buoyancy flux at interfaces
                           !! [H Z T-3 ~> m2 s-3 or W m-3]
    Bflx_temp => NULL(), & !< Temperature contribution to buoyancy flux at interfaces
                           !! [H Z T-3 ~> m2 s-3 or W m-3]
    Bflx_salt_dz => NULL(), & !< Salinity contribution to integral of buoyancy flux over layer
                              !! [H Z2 T-3 ~> m3 s-3 or W m-2]
    Bflx_temp_dz => NULL(), & !< Temperature contribution to integral of buoyancy flux over layer
                              !! [H Z2 T-3 ~> m3 s-3 or W m-2]
    ! The following are all allocatable arrays that store copies of process driven Kd, so that
    ! the process driven buoyancy flux and work can be derived at the end of the time step.
    Kd_salt => NULL(), &   !< total diapycnal diffusivity of salt at interfaces [H Z T-1 ~> m2 s-1 or kg m-1 s-1]
    Kd_temp => NULL(), &   !< total diapycnal diffusivity of heat at interfaces [H Z T-1 ~> m2 s-1 or kg m-1 s-1]
    Kd_BBL => NULL(), &    !< diapycnal diffusivity due to BBL at interfaces [H Z T-1 ~> m2 s-1 or kg m-1 s-1]
    Kd_ePBL => NULL(), &   !< diapycnal diffusivity due to ePBL at interfaces [H Z T-1 ~> m2 s-1 or kg m-1 s-1]
    Kd_KS => NULL(), &     !< diapycnal diffusivity due to Kappa Shear at interfaces [H Z T-1 ~> m2 s-1 or kg m-1 s-1]
    Kd_bkgnd => NULL(), &  !< diapycnal diffusivity due to Kd_bkgnd at interfaces [H Z T-1 ~> m2 s-1 or kg m-1 s-1]
    Kd_ddiff_S => NULL(), &!< diapycnal diffusivity due to double diffusion of salt at interfaces
                           !! [H Z T-1 ~> m2 s-1 or kg m-1 s-1]
    Kd_ddiff_T => NULL(), &!< diapycnal diffusivity due to double diffusion of heat at interfaces
                           !![H Z T-1 ~> m2 s-1 or kg m-1 s-1]
    Kd_leak => NULL(), &   !< diapycnal diffusivity due to Kd_leak at interfaces [H Z T-1 ~> m2 s-1 or kg m-1 s-1]
    Kd_quad => NULL(), &   !< diapycnal diffusivity due to Kd_quad at interfaces [H Z T-1 ~> m2 s-1 or kg m-1 s-1]
    Kd_itidal => NULL(), & !< diapycnal diffusivity due to Kd_itidal at interfaces [H Z T-1 ~> m2 s-1 or kg m-1 s-1]
    Kd_Froude => NULL(), & !< diapycnal diffusivity due to Kd_Froude at interfaces [H Z T-1 ~> m2 s-1 or kg m-1 s-1]
    Kd_slope => NULL(), &  !< diapycnal diffusivity due to Kd_slope at interfaces [H Z T-1 ~> m2 s-1 or kg m-1 s-1]
    Kd_lowmode => NULL(), &!< diapycnal diffusivity due to Kd_lowmode at interfaces [H Z T-1 ~> m2 s-1 or kg m-1 s-1]
    Kd_Niku => NULL(), &   !< diapycnal diffusivity due to Kd_Niku at interfaces [H Z T-1 ~> m2 s-1 or kg m-1 s-1]
    Kd_itides => NULL()    !< diapycnal diffusivity due to Kd_itides at interfaces [H Z T-1 ~> m2 s-1 or kg m-1 s-1]

  ! Constant Kd contributions
  real :: Kd_add !< spatially uniform additional diapycnal diffusivity at interfaces [H Z T-1 ~> m2 s-1 or kg m-1 s-1]
                 !! a diagnostic for this diffusivity is not yet included, but this makes it straightforward to add

  !>@{ Diagnostic IDs
  integer :: id_Bdif     = -1, id_Bdif_salt  = -1, id_Bdif_temp  = -1
  integer :: id_Bdif_dz  = -1, id_Bdif_salt_dz  = -1, id_Bdif_temp_dz  = -1
  integer :: id_Bdif_idz  = -1, id_Bdif_salt_idz  = -1, id_Bdif_temp_idz  = -1
  integer :: id_Bdif_idV  = -1, id_Bdif_salt_idV  = -1, id_Bdif_temp_idV  = -1
  integer :: id_Bdif_ePBL  = -1, id_Bdif_dz_ePBL  = -1, id_Bdif_idz_ePBL  = -1, id_Bdif_idV_ePBL  = -1
  integer :: id_Bdif_BBL  = -1, id_Bdif_dz_BBL  = -1, id_Bdif_idz_BBL  = -1, id_Bdif_idV_BBL  = -1
  integer :: id_Bdif_KS  = -1, id_Bdif_dz_KS  = -1, id_Bdif_idz_KS  = -1, id_Bdif_idV_KS  = -1
  integer :: id_Bdif_bkgnd  = -1, id_Bdif_dz_bkgnd  = -1, id_Bdif_idz_bkgnd  = -1, id_Bdif_idV_bkgnd  = -1
  integer :: id_Bdif_ddiff_temp  = -1, id_Bdif_ddiff_salt  = -1
  integer :: id_Bdif_dz_ddiff_temp  = -1, id_Bdif_dz_ddiff_salt  = -1
  integer :: id_Bdif_idz_ddiff_temp  = -1, id_Bdif_idz_ddiff_salt  = -1
  integer :: id_Bdif_idV_ddiff_temp  = -1, id_Bdif_idV_ddiff_salt  = -1
  integer :: id_Bdif_leak  = -1, id_Bdif_dz_leak  = -1, id_Bdif_idz_leak  = -1, id_Bdif_idV_leak  = -1
  integer :: id_Bdif_quad  = -1, id_Bdif_dz_quad  = -1, id_Bdif_idz_quad  = -1, id_Bdif_idV_quad  = -1
  integer :: id_Bdif_itidal  = -1, id_Bdif_dz_itidal  = -1, id_Bdif_idz_itidal  = -1, id_Bdif_idV_itidal  = -1
  integer :: id_Bdif_Froude  = -1, id_Bdif_dz_Froude  = -1, id_Bdif_idz_Froude  = -1, id_Bdif_idV_Froude  = -1
  integer :: id_Bdif_slope  = -1, id_Bdif_dz_slope  = -1, id_Bdif_idz_slope  = -1, id_Bdif_idV_slope  = -1
  integer :: id_Bdif_lowmode  = -1, id_Bdif_dz_lowmode  = -1, id_Bdif_idz_lowmode  = -1, id_Bdif_idV_lowmode  = -1
  integer :: id_Bdif_Niku  = -1, id_Bdif_dz_Niku  = -1, id_Bdif_idz_Niku  = -1, id_Bdif_idV_Niku  = -1
  integer :: id_Bdif_itides  = -1, id_Bdif_dz_itides  = -1, id_Bdif_idz_itides  = -1, id_Bdif_idV_itides  = -1
  !>@}

  logical :: do_bflx_salt = .false.  !< Logical flag to indicate if N2_salt should be computed
  logical :: do_bflx_temp = .false.  !< Logical flag to indicate if N2_temp should be computed
  logical :: do_bflx_salt_dz = .false.  !< Logical flag to indicate if N2_salt should be computed
  logical :: do_bflx_temp_dz = .false.  !< Logical flag to indicate if N2_temp should be computed

end type vbf_CS

! A note on unit descriptions in comments: MOM6 uses units that can be rescaled for dimensional
! consistency testing. These are noted in comments with units like Z, H, L, and T, along with
! their mks counterparts with notation like "a velocity [Z T-1 ~> m s-1]".  If the units
! vary with the Boussinesq approximation, the Boussinesq variant is given first.


  interface
module subroutine KdWork_Diagnostics(G,GV,US,diag,VBF,N2_Salt,N2_Temp,dz)
  type(ocean_grid_type),      intent(in)    :: G       !< Grid type
  type(verticalGrid_type),    intent(in)    :: GV      !< ocean vertical grid structure
  type(unit_scale_type),      intent(in)    :: US      !< A dimensional unit scaling type
  type(diag_ctrl), target,    intent(inout) :: diag    !< regulates diagnostic output
  type (vbf_CS),              intent(inout) :: VBF     !< Vertical buoyancy flux structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)+1), &
                              intent(in)    :: N2_Salt !< Buoyancy frequency [T-2 ~> s-2]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)+1), &
                              intent(in)    :: N2_Temp !< Buoyancy frequency [T-2 ~> s-2]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                              intent(in)    :: dz      !< Grid spacing [Z ~> m]

  ! Work arrays for computing buoyancy flux integrals


end subroutine KdWork_Diagnostics
module subroutine diagnoseKdWork(G, GV, N2, Kd, Bdif_flx, dz, Bdif_flx_dz)
  type(ocean_grid_type),   intent(in)  :: G    !< Grid type
  type(verticalGrid_type), intent(in)  :: GV   !< ocean vertical grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)+1), &
                           intent(in)  :: N2   !< Buoyancy frequency [T-2 ~> s-2]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)+1), &
                           intent(in)  :: Kd   !< Diffusivity [H Z T-1 ~> m2 s-1 or kg m-1 s-1]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)+1), &
                           intent(out) :: Bdif_flx !< Buoyancy flux [H Z T-3 ~> m2 s-3 or W m-3]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                 intent(in), optional :: dz    !< Grid spacing [Z ~> m]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                intent(out), optional :: Bdif_flx_dz !< Buoyancy flux over layer [H Z2 T-3 ~> m3 s-3 or W m-2]


end subroutine diagnoseKdWork
module subroutine Allocate_VBF_CS(G, GV, VBF)
  type(ocean_grid_type),   intent(in) :: G   !< ocean grid structure
  type(verticalGrid_type), intent(in) :: GV  !< ocean vertical grid structure
  type (vbf_CS),        intent(inout) :: VBF !< Vertical buoyancy flux structure


end subroutine Allocate_VBF_CS
module subroutine Deallocate_VBF_CS(VBF)
  type (vbf_CS), intent(inout) :: VBF !< Vertical buoyancy flux structure

end subroutine Deallocate_VBF_CS
module subroutine KdWork_init(Time, G,GV,US,diag,VBF,Use_KdWork_diag)
  type(time_type), target                :: Time             !< model time
  type(ocean_grid_type),   intent(in)    :: G        !< ocean grid structure
  type(verticalGrid_type), intent(in)    :: GV       !< ocean vertical grid structure
  type(unit_scale_type),   intent(in)    :: US       !< A dimensional unit scaling type
  type(diag_ctrl), target, intent(inout) :: diag     !< regulates diagnostic output
  type (vbf_CS), pointer,  intent(inout) :: VBF      !< Vertical buoyancy flux structure
  logical,                 intent(out)   :: Use_KdWork_diag !< Flag if any output was turned on

end subroutine KdWork_init
module subroutine KdWork_end(VBF)
  type (vbf_CS), pointer,  intent(inout) :: VBF      !< Vertical buoyancy flux structure

end subroutine KdWork_end
  end interface

end module MOM_diagnose_kdwork
