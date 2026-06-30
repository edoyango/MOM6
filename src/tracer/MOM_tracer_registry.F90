! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> This module contains subroutines that handle registration of tracers
!! and related subroutines. The primary subroutine, register_tracer, is
!! called to indicate the tracers advected and diffused.
!! It also makes public the types defined in MOM_tracer_types.
module MOM_tracer_registry

! use MOM_diag_mediator, only : diag_ctrl
use MOM_coms,          only : reproducing_sum
use MOM_debugging,     only : hchksum
use MOM_diag_mediator, only : diag_ctrl, register_diag_field, post_data, safe_alloc_ptr
use MOM_diag_mediator, only : diag_grid_storage
use MOM_diag_mediator, only : diag_copy_storage_to_diag, diag_save_grids, diag_restore_grids
use MOM_error_handler, only : MOM_error, FATAL, WARNING, MOM_mesg, is_root_pe
use MOM_file_parser,   only : get_param, log_version, param_file_type
use MOM_hor_index,     only : hor_index_type
use MOM_grid,          only : ocean_grid_type
use MOM_io,            only : vardesc, query_vardesc, cmor_long_std
use MOM_restart,       only : register_restart_field, MOM_restart_CS
use MOM_string_functions, only : lowercase
use MOM_time_manager,  only : time_type
use MOM_unit_scaling,  only : unit_scale_type
use MOM_verticalGrid,  only : verticalGrid_type
use MOM_tracer_types,  only : tracer_type, tracer_registry_type

implicit none ; private

#include <MOM_memory.h>

public register_tracer
public MOM_tracer_chksum, MOM_tracer_chkinv
public register_tracer_diagnostics
public post_tracer_diagnostics_at_sync, post_tracer_transport_diagnostics
public preALE_tracer_diagnostics, postALE_tracer_diagnostics
public tracer_registry_init, lock_tracer_registry, tracer_registry_end
public tracer_name_lookup
public tracer_type, tracer_registry_type

!> Write out checksums for registered tracers
interface MOM_tracer_chksum
  module procedure tracer_array_chksum, tracer_Reg_chksum
end interface MOM_tracer_chksum

!> Calculate and print the global inventories of registered tracers
interface MOM_tracer_chkinv
  module procedure tracer_array_chkinv, tracer_Reg_chkinv
end interface MOM_tracer_chkinv


  interface
module subroutine register_tracer(tr_ptr, Reg, param_file, HI, GV, name, longname, units, &
                           cmor_name, cmor_units, cmor_longname, net_surfflux_name, &
                           NLT_budget_name, net_surfflux_longname, tr_desc, OBC_inflow, &
                           OBC_in_u, OBC_in_v, ad_x, ad_y, df_x, df_y, ad_2d_x, ad_2d_y, &
                           df_2d_x, df_2d_y, advection_xy, registry_diags, &
                           conc_scale, flux_nameroot, flux_longname, flux_units, flux_scale, &
                           convergence_units, convergence_scale, cmor_tendprefix, diag_form, &
                           restart_CS, mandatory, underflow_conc, Tr_out, advect_scheme)
  type(hor_index_type),           intent(in)    :: HI           !< horizontal index type
  type(verticalGrid_type),        intent(in)    :: GV           !< ocean vertical grid structure
  type(tracer_registry_type),     pointer       :: Reg          !< pointer to the tracer registry
  real, dimension(SZI_(HI),SZJ_(HI),SZK_(GV)), &
                                  target        :: tr_ptr       !< target or pointer to the tracer array [CU ~> conc]
  type(param_file_type), intent(in)             :: param_file   !< file to parse for model parameter values
  character(len=*),     optional, intent(in)    :: name         !< Short tracer name
  character(len=*),     optional, intent(in)    :: longname     !< The long tracer name
  character(len=*),     optional, intent(in)    :: units        !< The units of this tracer
  character(len=*),     optional, intent(in)    :: cmor_name    !< CMOR name
  character(len=*),     optional, intent(in)    :: cmor_units   !< CMOR physical dimensions of variable
  character(len=*),     optional, intent(in)    :: cmor_longname !< CMOR long name
  character(len=*),     optional, intent(in)    :: net_surfflux_name     !< Name for net_surfflux diag
  character(len=*),     optional, intent(in)    :: NLT_budget_name       !< Name for NLT_budget diag
  character(len=*),     optional, intent(in)    :: net_surfflux_longname !< Long name for net_surfflux diag
  type(vardesc),        optional, intent(in)    :: tr_desc      !< A structure with metadata about the tracer

  real,                 optional, intent(in)    :: OBC_inflow   !< the tracer for all inflows via OBC for which OBC_in_u
                                                                !! or OBC_in_v are not specified [CU ~> conc]
  real, dimension(:,:,:), optional, pointer     :: OBC_in_u     !< tracer at inflows through u-faces of
                                                                !! tracer cells [CU ~> conc]
  real, dimension(:,:,:), optional, pointer     :: OBC_in_v     !< tracer at inflows through v-faces of
                                                                !! tracer cells [CU ~> conc]

  ! The following are probably not necessary if registry_diags is present and true.
  real, dimension(:,:,:), optional, pointer     :: ad_x         !< diagnostic x-advective flux
                                                                !! [CU H L2 T-1 ~> conc m3 s-1 or conc kg s-1]
  real, dimension(:,:,:), optional, pointer     :: ad_y         !< diagnostic y-advective flux
                                                                !! [CU H L2 T-1 ~> conc m3 s-1 or conc kg s-1]
  real, dimension(:,:,:), optional, pointer     :: df_x         !< diagnostic x-diffusive flux
                                                                !! [CU H L2 T-1 ~> conc m3 s-1 or conc kg s-1]
  real, dimension(:,:,:), optional, pointer     :: df_y         !< diagnostic y-diffusive flux
                                                                !! [CU H L2 T-1 ~> conc m3 s-1 or conc kg s-1]
  real, dimension(:,:),   optional, pointer     :: ad_2d_x      !< vert sum of diagnostic x-advect flux
                                                                !! [CU H L2 T-1 ~> conc m3 s-1 or conc kg s-1]
  real, dimension(:,:),   optional, pointer     :: ad_2d_y      !< vert sum of diagnostic y-advect flux
                                                                !! [CU H L2 T-1 ~> conc m3 s-1 or conc kg s-1]
  real, dimension(:,:),   optional, pointer     :: df_2d_x      !< vert sum of diagnostic x-diffuse flux
                                                                !! [CU H L2 T-1 ~> conc m3 s-1 or conc kg s-1]
  real, dimension(:,:),   optional, pointer     :: df_2d_y      !< vert sum of diagnostic y-diffuse flux
                                                                !! [CU H L2 T-1 ~> conc m3 s-1 or conc kg s-1]

  real, dimension(:,:,:), optional, pointer     :: advection_xy !< convergence of lateral advective tracer fluxes
                                                                !! [CU H T-1 ~> conc m s-1 or conc kg m-2 s-1]
  logical,              optional, intent(in)    :: registry_diags !< If present and true, use the registry for
                                                                !! the diagnostics of this tracer.
  real,                 optional, intent(in)    :: conc_scale   !< A scaling factor used to convert the concentration
                                                                !! of this tracer to its desired units [conc CU-1 ~> 1]
  character(len=*),     optional, intent(in)    :: flux_nameroot !< Short tracer name snippet used construct the
                                                                !! names of flux diagnostics.
  character(len=*),     optional, intent(in)    :: flux_longname !< A word or phrase used construct the long
                                                                !! names of flux diagnostics.
  character(len=*),     optional, intent(in)    :: flux_units   !< The units for the fluxes of this tracer.
  real,                 optional, intent(in)    :: flux_scale   !< A scaling factor used to convert the fluxes
                                                                !! of this tracer to its desired units
                                                                !! [conc m CU-1 H-1 ~> 1] or [conc kg m-2 CU-1 H-1 ~> 1]
  character(len=*),     optional, intent(in)    :: convergence_units !< The units for the flux convergence of
                                                                !! this tracer.
  real,                 optional, intent(in)    :: convergence_scale !< A scaling factor used to convert the flux
                                                                !! convergence of this tracer to its desired units.
                                                                !! [conc m CU-1 H-1 ~> 1] or [conc kg m-2 CU-1 H-1 ~> 1]
  character(len=*),     optional, intent(in)    :: cmor_tendprefix !< The CMOR name for the layer-integrated
                                                                !! tendencies of this tracer.
  integer,              optional, intent(in)    :: diag_form    !< An integer (1 or 2, 1 by default) indicating the
                                                                !! character string template to use in
                                                                !! labeling diagnostics
  type(MOM_restart_CS), optional, intent(inout) :: restart_CS   !< MOM restart control struct
  logical,              optional, intent(in)    :: mandatory    !< If true, this tracer must be read
                                                                !! from a restart file.
  real,                 optional, intent(in)    :: underflow_conc !< A tiny concentration, below which the tracer
                                                                !! concentration underflows to 0 [CU ~> conc].
  type(tracer_type),    optional, pointer       :: Tr_out       !< If present, returns pointer into registry

  integer,                 optional, intent(in) :: advect_scheme !< Advection scheme for this tracer, the default is -1
                                                                !! indicating to use the scheme from MOM_tracer_advect


end subroutine register_tracer
module subroutine lock_tracer_registry(Reg)
  type(tracer_registry_type), pointer    :: Reg    !< pointer to the tracer registry

end subroutine lock_tracer_registry
module subroutine register_tracer_diagnostics(Reg, h, Time, diag, G, GV, US, use_ALE, use_KPP)
  type(ocean_grid_type),      intent(in) :: G    !< The ocean's grid structure
  type(verticalGrid_type),    intent(in) :: GV   !< The ocean's vertical grid structure
  type(unit_scale_type),      intent(in) :: US   !< A dimensional unit scaling type
  type(tracer_registry_type), pointer    :: Reg  !< pointer to the tracer registry
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                              intent(in) :: h    !< Layer thicknesses [H ~> m or kg m-2]
  type(time_type),            intent(in) :: Time !< current model time
  type(diag_ctrl),            intent(in) :: diag !< structure to regulate diagnostic output
  logical,                    intent(in) :: use_ALE !< If true active diagnostics that only
                                                 !! apply to ALE configurations
  logical,                    intent(in) :: use_KPP !< If true active diagnostics that only
                                                 !! apply to CVMix KPP mixings

                                 ! creating additional diagnostics.
                                 ! [units] m3 s-1 or [units] kg s-1.
                                 ! [units] m s-1 or [units] kg m-2 s-1.
end subroutine register_tracer_diagnostics
module subroutine preALE_tracer_diagnostics(Reg, G, GV)
  type(tracer_registry_type), pointer    :: Reg  !< pointer to the tracer registry
  type(ocean_grid_type),      intent(in) :: G    !< The ocean's grid structure
  type(verticalGrid_type),    intent(in) :: GV   !< ocean vertical grid structure

end subroutine preALE_tracer_diagnostics
module subroutine postALE_tracer_diagnostics(Reg, G, GV, diag, dt)
  type(tracer_registry_type), pointer    :: Reg  !< pointer to the tracer registry
  type(ocean_grid_type),      intent(in) :: G    !< The ocean's grid structure
  type(verticalGrid_type),    intent(in) :: GV   !< ocean vertical grid structure
  type(diag_ctrl),            intent(in) :: diag !< regulates diagnostic output
  real,                       intent(in) :: dt   !< total time interval for these diagnostics [T ~> s]

end subroutine postALE_tracer_diagnostics
module subroutine post_tracer_diagnostics_at_sync(Reg, h, diag_prev, diag, G, GV, dt)
  type(ocean_grid_type),      intent(in) :: G    !< The ocean's grid structure
  type(verticalGrid_type),    intent(in) :: GV   !< The ocean's vertical grid structure
  type(tracer_registry_type), pointer    :: Reg  !< pointer to the tracer registry
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                              intent(in) :: h    !< Layer thicknesses [H ~> m or kg m-2]
  type(diag_grid_storage),    intent(in) :: diag_prev !< Contains diagnostic grids from previous timestep
  type(diag_ctrl),            intent(inout) :: diag !< structure to regulate diagnostic output
  real,                       intent(in) :: dt   !< total time step for tracer updates [T ~> s]

                                     ! in [CU H T-1 ~> conc m s-1 or conc kg m-2 s-1]
end subroutine post_tracer_diagnostics_at_sync
module subroutine post_tracer_transport_diagnostics(G, GV, Reg, h_diag, diag)
  type(ocean_grid_type),      intent(in) :: G    !< The ocean's grid structure
  type(verticalGrid_type),    intent(in) :: GV   !< The ocean's vertical grid structure
  type(tracer_registry_type), pointer    :: Reg  !< pointer to the tracer registry
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                              intent(in) :: h_diag !< Layer thicknesses on which to post fields [H ~> m or kg m-2]
  type(diag_ctrl),            intent(in) :: diag !< structure to regulate diagnostic output

                                          ! tracer fluxes [CU H T-1 ~> conc m s-1 or conc kg m-2 s-1]

end subroutine post_tracer_transport_diagnostics
module subroutine tracer_array_chksum(mesg, Tr, ntr, G)
  character(len=*),         intent(in) :: mesg   !< message that appears on the chksum lines
  type(tracer_type),        intent(in) :: Tr(:)  !< array of all of registered tracers
  integer,                  intent(in) :: ntr    !< number of registered tracers
  type(ocean_grid_type),    intent(in) :: G      !< ocean grid structure


end subroutine tracer_array_chksum
module subroutine tracer_Reg_chksum(mesg, Reg, G)
  character(len=*),           intent(in) :: mesg !< message that appears on the chksum lines
  type(tracer_registry_type), pointer    :: Reg  !< pointer to the tracer registry
  type(ocean_grid_type),      intent(in) :: G    !< ocean grid structure


end subroutine tracer_Reg_chksum
module subroutine tracer_array_chkinv(mesg, G, GV, h, Tr, ntr)
  character(len=*),                          intent(in) :: mesg !< message that appears on the chksum lines
  type(ocean_grid_type),                     intent(in) :: G    !< ocean grid structure
  type(verticalGrid_type),                   intent(in) :: GV   !< The ocean's vertical grid structure
  type(tracer_type), dimension(:),           intent(in) :: Tr   !< array of all of registered tracers
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in) :: h    !< Layer thicknesses [H ~> m or kg m-2]
  integer,                                   intent(in) :: ntr  !< number of registered tracers

  ! Local variables
                    ! masses to kg [kg H-1 L-2 ~> 1], depending on whether the Boussinesq approximation is used
                    ! each cell [conc m3] or [conc kg]

end subroutine tracer_array_chkinv
module subroutine tracer_Reg_chkinv(mesg, G, GV, h, Reg)
  character(len=*),                          intent(in) :: mesg !< message that appears on the chksum lines
  type(ocean_grid_type),                     intent(in) :: G    !< ocean grid structure
  type(verticalGrid_type),                   intent(in) :: GV   !< The ocean's vertical grid structure
  type(tracer_registry_type),                pointer    :: Reg  !< pointer to the tracer registry
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in) :: h    !< Layer thicknesses [H ~> m or kg m-2]

  ! Local variables
                    ! masses to kg [kg H-1 L-2 ~> 1], depending on whether the Boussinesq approximation is used
                    ! each cell [conc m3] or [conc kg]

end subroutine tracer_Reg_chkinv
module subroutine tracer_name_lookup(Reg, n, tr_ptr, name)
  type(tracer_registry_type), pointer    :: Reg     !< pointer to tracer registry
  type(tracer_type), pointer             :: tr_ptr  !< target or pointer to the tracer array
  character(len=32), intent(in)          :: name    !< tracer name
  integer, intent(out)                   :: n       !< index to tracer registery

end subroutine tracer_name_lookup
module subroutine tracer_registry_init(param_file, Reg)
  type(param_file_type),      intent(in) :: param_file !< open file to parse for model parameters
  type(tracer_registry_type), pointer    :: Reg        !< pointer to tracer registry


! This include declares and sets the variable "version".

end subroutine tracer_registry_init
module subroutine tracer_registry_end(Reg)
  type(tracer_registry_type), pointer :: Reg !< The tracer registry that will be deallocated
end subroutine tracer_registry_end
  end interface

end module MOM_tracer_registry
