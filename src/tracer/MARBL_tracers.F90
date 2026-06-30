! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> A tracer package for tracers computed in the MARBL library
!!
!! Currently configured for use with marbl0.36.0
!! https://github.com/marbl-ecosys/MARBL/releases/tag/marbl0.36.0
!! (clone entire repo into pkg/MARBL)
module MARBL_tracers

use MOM_coms,            only : EFP_type, root_PE, broadcast
use MOM_debugging,       only : hchksum
use MOM_diag_mediator,   only : diag_ctrl
use MOM_error_handler,   only : is_root_PE, MOM_error, FATAL, WARNING, NOTE
use MOM_file_parser,     only : get_param, log_param, log_version, param_file_type
use MOM_forcing_type,    only : forcing
use MOM_grid,            only : ocean_grid_type
use MOM_interpolate,     only : external_field, init_external_field, time_interp_external
use MOM_CVMix_KPP,       only : KPP_NonLocalTransport, KPP_CS
use MOM_hor_index,       only : hor_index_type
use MOM_interpolate,     only : forcing_timeseries_dataset
use MOM_interpolate,     only : forcing_timeseries_set_time_type_vars
use MOM_interpolate,     only : map_model_time_to_forcing_time
use MOM_io,              only : file_exists, MOM_read_data, slasher, vardesc, var_desc, query_vardesc
use MOM_open_boundary,   only : ocean_OBC_type
use MOM_remapping,       only : reintegrate_column
use MOM_remapping,       only : remapping_CS, initialize_remapping, remapping_core_h
use MOM_restart,         only : query_initialized, MOM_restart_CS, register_restart_field
use MOM_spatial_means,   only : global_mass_int_EFP
use MOM_sponge,          only : set_up_sponge_field, sponge_CS
use MOM_time_manager,    only : time_type
use MOM_tracer_registry, only : register_tracer
use MOM_tracer_types,    only : tracer_type, tracer_registry_type
use MOM_tracer_diabatic, only : tracer_vertdiff, applyTracerBoundaryFluxesInOut
use MOM_tracer_initialization_from_Z, only : MOM_initialize_tracer_from_Z
use MOM_tracer_Z_init,   only : read_Z_edges
use MOM_unit_scaling,    only : unit_scale_type
use MOM_variables,       only : surface, thermo_var_ptrs
use MOM_verticalGrid,    only : verticalGrid_type
use MOM_diag_mediator,   only : register_diag_field, post_data!, safe_alloc_ptr

use MARBL_interface,              only : MARBL_interface_class
use MARBL_interface_public_types, only : marbl_diagnostics_type, marbl_saved_state_type

use atmos_ocean_fluxes_mod, only : aof_set_coupler_flux

implicit none ; private

#include <MOM_memory.h>

public register_MARBL_tracers, initialize_MARBL_tracers
public MARBL_tracers_column_physics, MARBL_tracers_surface_state
public MARBL_tracers_set_forcing
public MARBL_tracers_stock, MARBL_tracers_get, MARBL_tracers_end

! A note on unit descriptions in comments: MOM6 uses units that can be rescaled for dimensional
! consistency testing. These are noted in comments with units like Z, H, L, and T, along with
! their mks counterparts with notation like "a velocity [Z T-1 ~> m s-1]".  If the units
! vary with the Boussinesq approximation, the Boussinesq variant is given first.

!> Temporary type for diagnostic variables coming from MARBL
!! Allocate exactly one of field_[23]d
type :: temp_MARBL_diag
  integer :: id !< index into MOM diagnostic structure
  real, allocatable :: field_2d(:,:)   !< memory for 2D field
  real, allocatable :: field_3d(:,:,:) !< memory for 3D field
end type temp_MARBL_diag

!> MOM6 needs to know the index of some MARBL tracers to properly apply river fluxes
type :: tracer_ind_type
  integer :: no3_ind  !< NO3 index
  integer :: po4_ind  !< PO4 index
  integer :: don_ind  !< DON index
  integer :: donr_ind  !< DONr index
  integer :: dop_ind  !< DOP index
  integer :: dopr_ind  !< DOPr index
  integer :: sio3_ind  !< SiO3 index
  integer :: fe_ind  !< Fe index
  integer :: doc_ind  !< DOC index
  integer :: docr_ind  !< DOCr index
  integer :: alk_ind  !< ALK index
  integer :: alk_alt_co2_ind  !< ALK_ALT_CO2 index
  integer :: dic_ind  !< DIC index
  integer :: dic_alt_co2_ind  !< DIC_ALT_CO2 index
  integer :: abio_dic_ind  !< ABIO_DIC index
  integer :: abio_di14c_ind  !< ABIO_DI14C index
end type tracer_ind_type

!> MOM needs to store some information about saved_state; besides providing these
!! fields to MARBL, they are also written to restart files
type :: saved_state_for_MARBL_type
  character(len=200) :: short_name !< name of variable being saved
  character(len=200) :: file_varname !< name of variable in restart file
  character(len=200) :: units !< variable units
  real, pointer :: field_2d(:,:) => NULL()   !< memory for 2D field
  real, pointer :: field_3d(:,:,:) => NULL() !< memory for 3D field
end type saved_state_for_MARBL_type

!> All calls to MARBL are done via the interface class
type(MARBL_interface_class) :: MARBL_instances

!> Pointer to tracer concentration and to tracer_type in tracer registry
type, private :: MARBL_tracer_data
  real, pointer              :: tr(:,:,:) => NULL() !< Array of tracers used in this subroutine [CU ~> conc]
                                                    !! (ALK tracers use meq m-3 instead of mmol m-3)
  type(tracer_type), pointer :: tr_ptr    => NULL() !< pointer to tracer inside Tr_reg
end type MARBL_tracer_data

!> The control structure for the MARBL tracer package
type, public :: MARBL_tracers_CS ; private
  integer :: ntr          !< The number of tracers that are actually used.
  logical :: debug        !< If true, write verbose checksums for debugging purposes.
  logical :: base_bio_on  !< Will MARBL use base biotic tracers?
  logical :: abio_dic_on  !< Will MARBL use abiotic DIC / DI14C tracers?
  logical :: ciso_on      !< Will MARBL use isotopic tracers?

  integer :: restore_count              !< The number of tracers MARBL is configured to restore
  logical :: coupled_tracers = .false.  !< These tracers are not offered to the coupler.
  logical :: use_ice_category_fields    !< Forcing will include multiple ice categories for ice_frac and shortwave
  logical :: request_Chl_from_MARBL     !< MARBL can provide Chl to use in set_pen_shortwave()
  integer :: ice_ncat                   !< Number of ice categories when use_ice_category_fields = True
  real    :: IC_min                     !< Minimum value for tracer initial conditions [CU ~> conc]
  character(len=200) :: IC_file !< The file in which the age-tracer initial values cam be found.
  logical :: ongrid                     !< True if IC_file is already interpolated to MOM grid
  type(tracer_registry_type), pointer :: tr_Reg => NULL() !< A pointer to the tracer registry
  type(MARBL_tracer_data), dimension(:), allocatable :: tracer_data  !< type containing tracer data and pointer
                                                                     !! into tracer registry

  integer, allocatable, dimension(:) :: ind_tr !< Indices returned by aof_set_coupler_flux if it is used and the
                                               !! surface tracer concentrations are to be provided to the coupler.

  type(diag_ctrl), pointer :: diag => NULL() !< A structure that is used to
                                   !! regulate the timing of diagnostic output.
  type(MOM_restart_CS), pointer :: restart_CSp => NULL() !< A pointer to the restart control structure

  type(vardesc), allocatable :: tr_desc(:) !< Descriptions and metadata for the tracers
  logical :: tracers_may_reinit            !< If true the tracers may be initialized if not found in a restart file

  character(len=200) :: fesedflux_file   !< name of [netCDF] file containing iron sediment flux
  character(len=200) :: feventflux_file  !< name of [netCDF] file containing iron vent flux
  type(forcing_timeseries_dataset) :: d14c_dataset(3) !< File and time axis information for d14c forcing
  real, dimension(3) :: d14c_bands       !< forcing is organized into bands: [30 N, 90 N]; [30 S, 30 N]; [90 S, 30 S]
                                         !! This variable contains D14C for each band [CU ~> conc]
  integer :: d14c_id                     !< id for diagnostic field with d14c forcing
  logical :: read_riv_fluxes             !< If true, use river fluxes supplied from an input file.
                                         !! This is temporary, we will always read river fluxes
  type(forcing_timeseries_dataset) :: riv_flux_dataset !< File and time axis information for river fluxes
  character(len=4) :: restoring_source !< location of tracer restoring data
                                       !! valid values: file, none
  integer :: restoring_nz  !< number of levels in tracer restoring file
  real, allocatable, dimension(:) :: &
      restoring_z_edges  !< The depths of the cell interfaces in the tracer restoring file [Z ~> m]
  real, allocatable, dimension(:) :: &
      restoring_dz  !< The thickness of the cell layers in the tracer restoring file [H ~> m]
  integer :: restoring_timescale_nz  !< number of levels in tracer restoring timescale file
  real, allocatable, dimension(:) :: &
      restoring_timescale_z_edges  !< The depths of the cell interfaces in the tracer restoring timescale file [Z ~> m]
  real, allocatable, dimension(:) :: &
      restoring_timescale_dz  !< The thickness of the cell layers in the tracer restoring timescale file [H ~> m]
  character(len=14) :: restoring_I_tau_source !< location of inverse restoring timescale data
                                              !! valid values: file, grid_dependent
  character(len=200) :: restoring_file !< name of [netCDF] file containing tracer restoring data
  type(remapping_CS) :: restoring_remapCS !< Remapping parameters and work arrays for tracer restoring / timescale
  character(len=200) :: restoring_I_tau_file !< name of [netCDF] file containing inverse restoring timescale
  character(len=200) :: restoring_I_tau_var_name !< name of field containing inverse restoring timescale
  character(len=35) :: marbl_settings_file  !< name of [text] file containing MARBL settings

  real :: bot_flux_mix_thickness !< for bottom flux -> tendency conversion, assume uniform mixing over
                                 !! bottom layer of prescribed thickness [Z ~> m]
  real :: Ibfmt                  !< Reciprocal of bot_flux_mix_thickness [Z-1 ~> m-1]

  type(temp_MARBL_diag), allocatable :: surface_flux_diags(:)  !< collect surface flux diagnostics from all columns
                                                               !! before posting
  type(temp_MARBL_diag), allocatable :: interior_tendency_diags(:)  !< collect tendency diagnostics from all columns
                                                                    !! before posting
  type(saved_state_for_MARBL_type), allocatable :: surface_flux_saved_state(:)  !< surface_flux saved state
  type(saved_state_for_MARBL_type), allocatable :: interior_tendency_saved_state(:)  !< interior_tendency saved state

  ! TODO: If we can post data column by column, all we need are integer arrays for ids
  ! integer, allocatable :: id_surface_flux_diags(:)  !< array of indices for surface_flux diagnostics
  ! integer, allocatable :: id_interior_tendency_diags(:)  !< array of indices for interior_tendency diagnostics

  type(tracer_ind_type) :: tracer_inds  !< Indices to tracers that will have river fluxes added to STF

  !> Need to store global output from both marbl_instance%surface_flux_compute() and
  !! marbl_instance%interior_tendency_compute(). For the former, just need id to register
  !! because we already copy data into CS%STF; latter requires copying data and indices
  !! so currently using temp_MARBL_diag for that.
  integer, allocatable :: id_surface_flux_out(:)  !< register_diag indices for surface_flux output
  integer, allocatable :: id_surface_flux_from_salt_flux(:)  !< register_diag indices for surface_flux from salt_flux
  type(temp_MARBL_diag), allocatable :: interior_tendency_out(:)  !< collect interior tendencies for diagnostic output
  type(temp_MARBL_diag), allocatable :: interior_tendency_out_zint(:)  !< vertical integral of interior tendencies
                                                                       !! (full column)
  type(temp_MARBL_diag), allocatable :: interior_tendency_out_zint_100m(:)  !< vertical integral of interior tendencies
                                                                            !! (top 100m)
  integer :: bot_flux_to_tend_id  !< register_diag index for BOT_FLUX_TO_TEND
  integer, allocatable :: fracr_cat_id(:) !< register_diag index for per-category ice fraction
  integer, allocatable :: qsw_cat_id(:)   !< register_diag index for per-category shortwave

  real :: DIC_salt_ratio !< ratio to convert salt surface flux to DIC surface flux [conc ppt-1]
  real :: ALK_salt_ratio !< ratio to convert salt surface flux to ALK surface flux [conc ppt-1]

  real, allocatable :: STF(:,:,:)          !< surface fluxes returned from MARBL to use in tracer_vertdiff
                                           !! (dims: i, j, tracer) [conc Z T-1 ~> conc m s-1]
  real, pointer :: SFO(:,:,:) => NULL()    !< surface flux output returned from MARBL for use in GCM
                                           !! e.g. CO2 flux to pass to atmosphere (dims: i, j, num_sfo)
                                           !! Units vary based on index of num_sfo dimension
  real, pointer :: ITO(:,:,:,:) => NULL()  !< interior tendency output returned from MARBL for use in GCM
                                           !! e.g. total chlorophyll to use in shortwave penetration
                                           !! (dims: i, j, k, num_ito)
                                           !! Units vary based on index of num_ito dimension

  integer :: u10_sqr_ind   !< index of MARBL forcing field array to copy 10-m wind (squared) into
  integer :: sss_ind       !< index of MARBL forcing field array to copy sea surface salinity into
  integer :: sst_ind       !< index of MARBL forcing field array to copy sea surface temperature into
  integer :: ifrac_ind     !< index of MARBL forcing field array to copy ice fraction into
  integer :: dust_dep_ind  !< index of MARBL forcing field array to copy dust flux into
  integer :: fe_dep_ind    !< index of MARBL forcing field array to copy iron flux into
  integer :: nox_flux_ind  !< index of MARBL forcing field array to copy NOx flux into
  integer :: nhy_flux_ind  !< index of MARBL forcing field array to copy NHy flux into
  integer :: atmpress_ind  !< index of MARBL forcing field array to copy atmospheric pressure into
  integer :: xco2_ind      !< index of MARBL forcing field array to copy CO2 flux into
  integer :: xco2_alt_ind  !< index of MARBL forcing field array to copy CO2 flux (alternate CO2) into
  integer :: d14c_ind      !< index of MARBL forcing field array to copy d14C into

  !> external_field types for river fluxes (added to surface fluxes)
  type(external_field) :: id_din_riv     !< id for time_interp_external.
  type(external_field) :: id_don_riv     !< id for time_interp_external.
  type(external_field) :: id_dip_riv     !< id for time_interp_external.
  type(external_field) :: id_dop_riv     !< id for time_interp_external.
  type(external_field) :: id_dsi_riv     !< id for time_interp_external.
  type(external_field) :: id_dfe_riv     !< id for time_interp_external.
  type(external_field) :: id_dic_riv     !< id for time_interp_external.
  type(external_field) :: id_alk_riv     !< id for time_interp_external.
  type(external_field) :: id_doc_riv     !< id for time_interp_external.

  !> external_field type for d14c (needed if abio_dic_on is True)
  type(external_field) :: id_d14c(3)        !< id for time_interp_external.

  !> Indices for river fluxes (diagnostics)
  integer :: no3_riv_flux          !< NO3 riverine flux
  integer :: po4_riv_flux          !< PO4 riverine flux
  integer :: don_riv_flux          !< DON riverine flux
  integer :: donr_riv_flux         !< DONr riverine flux
  integer :: dop_riv_flux          !< DOP riverine flux
  integer :: dopr_riv_flux         !< DOPr riverine flux
  integer :: sio3_riv_flux         !< SiO3 riverine flux
  integer :: fe_riv_flux           !< Fe riverine flux
  integer :: doc_riv_flux          !< DOC riverine flux
  integer :: docr_riv_flux         !< DOCr riverine flux
  integer :: alk_riv_flux          !< ALK riverine flux
  integer :: alk_alt_co2_riv_flux  !< ALK (alternate CO2) riverine flux
  integer :: dic_riv_flux          !< DIC riverine flux
  integer :: dic_alt_co2_riv_flux  !< DIC (alternate CO2) riverine flux

  !> Indices for forcing fields required to compute interior tendencies
  integer :: dustflux_ind  !< index of MARBL forcing field array to copy dust flux into
  integer :: PAR_col_frac_ind  !< index of MARBL forcing field array to copy PAR column fraction into
  integer :: surf_shortwave_ind  !< index of MARBL forcing field array to copy surface shortwave into
  integer :: potemp_ind  !< index of MARBL forcing field array to copy potential temperature into
  integer :: salinity_ind  !< index of MARBL forcing field array to copy salinity into
  integer :: pressure_ind  !< index of MARBL forcing field array to copy pressure into
  integer :: fesedflux_ind  !< index of MARBL forcing field array to copy iron sediment flux into
  integer :: o2_scalef_ind  !< index of MARBL forcing field array to copy O2 scale length into
  integer :: remin_scalef_ind  !< index of MARBL forcing field array to copy remin scale length into
  type(external_field), allocatable :: id_tracer_restoring(:) !< id number for time_interp_external
  integer, allocatable :: tracer_restoring_ind(:) !< index of MARBL forcing field to copy
                                                  !! per-tracer restoring field into
  integer, allocatable :: tracer_I_tau_ind(:) !< index of MARBL forcing field to copy per-tracer
                                              !! inverse restoring timescale into

  !> Memory for storing river fluxes, tracer restoring fields, and abiotic forcing
  real, allocatable :: d14c(:,:)         !< d14c forcing for abiotic DIC and carbon isotope tracer modules
                                         !! [mmol m-3 s-1]
  real, allocatable :: RIV_FLUXES(:,:,:) !< river flux forcing for applyTracerBoundaryFluxesInOut
                                         !! (needs to be time-integrated when passed to function!)
                                         !! (dims: i, j, tracer) [conc m s-1]
  character(len=15), allocatable :: tracer_restoring_varname(:) !< name of variable being restored
  real, allocatable :: I_tau(:,:,:)  !< inverse restoring timescale for marbl tracers (dims: i, j, k) [s-1]
  real, allocatable, dimension(:,:,:,:) :: restoring_in  !< Restoring fields read from file
                                                         !! (dims: i, j, restoring_nz, restoring_cnt) [tracer units]

  !> Number of surface flux outputs as well as specific indices for each one
  integer :: sfo_cnt       !< number of surface flux outputs from MARBL
  integer :: ito_cnt       !< number of interior tendency outputs from MARBL
  integer :: flux_co2_ind  !< index to co2 flux surface flux output
  integer :: total_Chl_ind !< index to total chlorophyll interior tendency output

  ! TODO: create generic 3D forcing input type to read z coordinate + values
  real    :: fesedflux_scale_factor !< scale factor for iron sediment flux [mmol umol-1 d s-1]
  integer :: fesedflux_nz  !< number of levels in iron sediment flux file
  real, allocatable, dimension(:,:,:) :: fesedflux_in  !< Field to read iron sediment flux into [conc m s-1]
  real, allocatable, dimension(:,:,:) :: feventflux_in  !< Field to read iron vent flux into [conc m s-1]
  real, allocatable, dimension(:) :: &
    fesedflux_z_edges  !< The depths of the cell interfaces in the input data [Z ~> m]
  ! TODO: this thickness does not need to be 3D, but it is easier to make thickness 0
  !       below the surface on a per-column basis (could save memory by storing 1D
  !       thickness from file and then computing a second 1D thickness array in (i,j) loop)
  real, allocatable, dimension(:,:,:) :: &
    fesedflux_dz  !< The thickness of the cell layers in the input data [H ~> m]
end type MARBL_tracers_CS

! Module parameters
real, parameter :: atm_per_Pa = 1./101325.  !< convert from Pa -> atm [atm Pa-1]


  interface
module subroutine configure_MARBL_tracers(GV, US, param_file, CS)
  type(verticalGrid_type),    intent(in) :: GV   !< The ocean's vertical grid structure
  type(unit_scale_type),      intent(in) :: US   !< A dimensional unit scaling type
  type(param_file_type),      intent(in) :: param_file !< A structure to parse for run-time parameters
  type(MARBL_tracers_CS),     pointer    :: CS   !< A pointer that is set to point to the control
                                                 !! structure for this module

end subroutine configure_MARBL_tracers
module function register_MARBL_tracers(HI, GV, US, param_file, CS, tr_Reg, restart_CS, MARBL_computes_chl)
  type(hor_index_type),       intent(in) :: HI   !< A horizontal index type structure.
  type(verticalGrid_type),    intent(in) :: GV   !< The ocean's vertical grid structure
  type(unit_scale_type),      intent(in) :: US   !< A dimensional unit scaling type
  type(param_file_type),      intent(in) :: param_file !< A structure to parse for run-time parameters
  type(MARBL_tracers_CS),     pointer    :: CS   !< A pointer that is set to point to the control
                                                 !! structure for this module
  type(tracer_registry_type), pointer    :: tr_Reg !< A pointer that is set to point to the control
                                                 !! structure for the tracer advection and diffusion module.
  type(MOM_restart_CS), target, intent(inout) :: restart_CS !< MOM restart control struct
  logical,                      intent(out)   :: MARBL_computes_chl  !< If MARBL is computing chlorophyll, MOM
                                                                     !! may use it to compute SW penetration

! Local variables
! This include declares and sets the variable "version".
                                           ! (ALK tracers use meq m-3 instead of mmol m-3)
  logical :: register_MARBL_tracers
  ! read_Z_edges() has several mandatory arguments that we do not use given our expectation
  ! of how the file being read in was created
end function register_MARBL_tracers
module subroutine initialize_MARBL_tracers(restart, day, G, GV, US, h, param_file, diag, OBC, CS, sponge_CSp)
  logical,                               intent(in)    :: restart      !< .true. if the fields have already been
                                                                       !! read from a restart file.
  type(time_type), target,               intent(in)    :: day          !< Time of the start of the run.
  type(ocean_grid_type),                 intent(inout) :: G            !< The ocean's grid structure
  type(verticalGrid_type),               intent(in)    :: GV           !< The ocean's vertical grid structure
  type(unit_scale_type),                 intent(in)    :: US           !< A dimensional unit scaling type
  real, dimension(NIMEM_,NJMEM_,NKMEM_), intent(in)    :: h            !< Layer thicknesses [H ~> m or kg m-2]
  type(param_file_type),                 intent(in)    :: param_file   !< A structure to parse for run-time parameters
  type(diag_ctrl), target,               intent(in)    :: diag         !< Structure used to regulate diagnostic output.
  type(ocean_OBC_type),                  pointer       :: OBC          !< This open boundary condition type specifies
                                                                       !! whether, where, and what open boundary
                                                                       !! conditions are used.
  type(MARBL_tracers_CS),                pointer       :: CS           !< The control structure returned by a previous
                                                                       !! call to register_MARBL_tracers.
  type(sponge_CS),                       pointer       :: sponge_CSp   !< A pointer to the control structure
                                                                       !! for the sponges, if they are in use.

  ! Local variables
                                  ! years m3 s-1 or years kg s-1.

end subroutine initialize_MARBL_tracers
module subroutine register_MARBL_diags(MARBL_diags, diag, day, G, id_diags)

  type(marbl_diagnostics_type), intent(in)    :: MARBL_diags !< MARBL diagnostics from MARBL_instances
  type(time_type), target,      intent(in)    :: day  !< Time of the start of the run.
  type(diag_ctrl), target,      intent(in)    :: diag !< Structure used to regulate diagnostic output.
  !integer, allocatable,         intent(inout) :: id_diags(:) !< allocatable array storing diagnostic index number
  type(ocean_grid_type),              intent(in) :: G    !< The ocean's grid structure
  type(temp_marbl_diag), allocatable, intent(inout) :: id_diags(:) !< allocatable array storing diagnostic index
                                                                   !! number and buffer space for collecting diags
                                                                   !! from all columns


end subroutine register_MARBL_diags
module subroutine setup_saved_state(MARBL_saved_state, HI, GV, restart_CS, tracers_may_reinit, &
    local_saved_state)

  type(marbl_saved_state_type),                  intent(in)    :: MARBL_saved_state !< MARBL saved state from
                                                                                    !! MARBL_instances
  type(hor_index_type),                          intent(in)    :: HI   !< A horizontal index type structure.
  type(verticalGrid_type),                       intent(in)    :: GV   !< The ocean's vertical grid structure
  type(MOM_restart_CS), pointer,                 intent(in)    :: restart_CS !< control structure to add saved state
                                                                             !! to restarts
  logical,                                       intent(in)    :: tracers_may_reinit  !< used to determine mandatory
                                                                                      !! flag in restart
  type(saved_state_for_MARBL_type), allocatable, intent(inout) :: local_saved_state(:) !< allocatable array for local
                                                                                       !! saved state


end subroutine setup_saved_state
module subroutine MARBL_tracers_column_physics(h_old, h_new, ea, eb, fluxes, dt, G, GV, US, CS, tv, &
    KPP_CSp, nonLocalTrans, evap_CFL_limit, minimum_forcing_depth)

  type(ocean_grid_type),   intent(in) :: G    !< The ocean's grid structure
  type(verticalGrid_type), intent(in) :: GV   !< The ocean's vertical grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), &
                           intent(in) :: h_old !< Layer thickness before entrainment [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), &
                           intent(in) :: h_new !< Layer thickness after entrainment [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), &
                           intent(in) :: ea   !< an array to which the amount of fluid entrained
                                              !! from the layer above during this call will be
                                              !! added [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), &
                           intent(in) :: eb   !< an array to which the amount of fluid entrained
                                              !! from the layer below during this call will be
                                              !! added [H ~> m or kg m-2].
  type(forcing),           intent(in) :: fluxes !< A structure containing pointers to thermodynamic
                                              !! and tracer forcing fields.  Unused fields have NULL ptrs.
  real,                    intent(in) :: dt   !< The amount of time covered by this call [T ~> s]
  type(unit_scale_type),   intent(in) :: US   !< A dimensional unit scaling type
  type(MARBL_tracers_CS),     pointer :: CS   !< The control structure returned by a previous
                                              !! call to register_MARBL_tracers.
  type(thermo_var_ptrs),   intent(in) :: tv   !< A structure pointing to various thermodynamic variables
  type(KPP_CS),  optional, pointer    :: KPP_CSp  !< KPP control structure
  real,          optional, intent(in) :: nonLocalTrans(:,:,:) !< Non-local transport [1]
  real,          optional, intent(in) :: evap_CFL_limit !< Limit on the fraction of the water that can
                                              !! be fluxed out of the top layer in a timestep [1]
  real,          optional, intent(in) :: minimum_forcing_depth !< The smallest depth over which
                                              !! fluxes can be applied [m]

  ! Local variables
                                                     ! [S H T-1 ~> ppt m s-1 or ppt kg m-2 s-1].
                                                          ! [conc Z T-1 ~> conc m s-1].
                                                                ! [Z-1 ~> m-1]

end subroutine MARBL_tracers_column_physics
module subroutine MARBL_tracers_set_forcing(day_start, G, CS)

  type(time_type),         intent(in)    :: day_start !< Start time of the fluxes.
  type(ocean_grid_type),   intent(in)    :: G         !< The ocean's grid structure.
  type(MARBL_tracers_CS),  pointer       :: CS        !< The control structure returned by a


                                                   !! [mmol m-2 s-1]

end subroutine MARBL_tracers_set_forcing
module function MARBL_tracers_stock(h, stocks, G, GV, CS, names, units, stock_index)
  real, dimension(NIMEM_,NJMEM_,NKMEM_), intent(in)    :: h      !< Layer thicknesses [H ~> m or kg m-2]
  type(EFP_type), dimension(:),          intent(out)   :: stocks !< the mass-weighted integrated amount of
                                                                 !! each tracer, in kg times concentration units
                                                                 !! [kg conc].
  type(ocean_grid_type),                 intent(in)    :: G      !< The ocean's grid structure
  type(verticalGrid_type),               intent(in)    :: GV     !< The ocean's vertical grid structure
  type(MARBL_tracers_CS),                pointer       :: CS     !< The control structure returned by a
                                                                 !! previous call to register_MARBL_tracers.
  character(len=*), dimension(:),        intent(out)   :: names  !< the names of the stocks calculated.
  character(len=*), dimension(:),        intent(out)   :: units  !< the units of the stocks calculated.
  integer, optional,                     intent(in)    :: stock_index !< the coded index of a specific stock
                                                                      !! being sought.
  integer                                              :: MARBL_tracers_stock   !< Return value: the number of stocks
                                                                                !! calculated here.

  ! Local variables
end function MARBL_tracers_stock
module subroutine MARBL_tracers_surface_state(sfc_state, G, US, CS)
  type(ocean_grid_type),   intent(in)    :: G   !< The ocean's grid structure.
  type(surface),           intent(inout) :: sfc_state !< A structure containing fields that
                                                      !! describe the surface state of the ocean.
  type(unit_scale_type),   intent(in)    :: US  !< A dimensional unit scaling type
  type(MARBL_tracers_CS),  pointer       :: CS  !< The control structure returned by a previous
                                                !! call to register_MARBL_tracers.

  ! Local variables

end subroutine MARBL_tracers_surface_state
module subroutine MARBL_tracers_get(name, G, GV, array, CS)

  character(len=*),         intent(in)    :: name   !< Name of requested tracer.
  type(ocean_grid_type),    intent(in)    :: G      !< The ocean's grid structure.
  type(verticalGrid_type),  intent(in)    :: GV     !< The ocean's vertical grid structure.
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                            intent(inout) :: array  !< Array filled by this routine.
  type(MARBL_tracers_CS),   pointer       :: CS     !< Pointer to the control structure for this module.


end subroutine MARBL_tracers_get
module subroutine MARBL_tracers_end(CS)
  type(MARBL_tracers_CS), pointer, intent(inout) :: CS !< The control structure returned by a previous
                                                       !! call to register_MARBL_tracers.


end subroutine MARBL_tracers_end
module subroutine set_riv_flux_tracer_inds(CS)

  type(MARBL_tracers_CS), pointer, intent(inout) :: CS   !< The MARBL tracers control structure


  ! Initialize tracers from file (unless they were initialized by restart file)
  ! Also save indices of tracers that have river fluxes
end subroutine set_riv_flux_tracer_inds
module subroutine print_marbl_log(log_to_print, G, i, j)

  use marbl_logging, only : marbl_status_log_entry_type
  use marbl_logging, only : marbl_log_type
  use MOM_coms,      only : PE_here

  class(marbl_log_type),           intent(in) :: log_to_print  !< MARBL log to include in MOM6 logfile
  type(ocean_grid_type), optional, intent(in) :: G             !< The ocean's grid structure
  integer,               optional, intent(in) :: i             !< i of (i,j) index of column providing the log
  integer,               optional, intent(in) :: j             !< j of (i,j) index of column providing the log


  ! elem_old is used to keep track of whether all messages are coming from the same point
end subroutine print_marbl_log
  end interface

end module MARBL_tracers
