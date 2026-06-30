submodule (MOM_coms_infra) MOM_coms_infra_s
  implicit none
contains
module procedure PE_here
  pe = mpp_pe()
end procedure PE_here
module procedure root_PE
  pe = mpp_root_pe()
end procedure root_PE
module procedure num_PEs
  npes = mpp_npes()
end procedure num_PEs
module procedure set_rootPE
  call mpp_set_root_pe(pe)
end procedure set_rootPE
module procedure Set_PEList
  call mpp_set_current_pelist(pelist, no_sync)
end procedure Set_PEList
module procedure Get_PEList
  call mpp_get_current_pelist(pelist, name, commiD)
end procedure Get_PEList
module procedure sync_PEs
  call mpp_sync(pelist)
end procedure sync_PEs
module procedure broadcast_char
  integer :: src_PE   ! The processor that is sending the data
  logical :: do_block ! If true add synchronizing barriers
  do_block = .false. ; if (present(blocking)) do_block = blocking
  if (present(from_PE)) then ; src_PE = from_PE ; else ; src_PE = root_PE() ; endif

  if (do_block) call mpp_sync(PElist)
  call mpp_broadcast(dat, length, src_PE, PElist)
  if (do_block) call mpp_sync_self(PElist)

end procedure broadcast_char
module procedure broadcast_int64_0D
  integer :: src_PE   ! The processor that is sending the data
  logical :: do_block ! If true add synchronizing barriers
  do_block = .false. ; if (present(blocking)) do_block = blocking
  if (present(from_PE)) then ; src_PE = from_PE ; else ; src_PE = root_PE() ; endif

  if (do_block) call mpp_sync(PElist)
  call mpp_broadcast(dat, src_PE, PElist)
  if (do_block) call mpp_sync_self(PElist)

end procedure broadcast_int64_0D
module procedure broadcast_int32_0D
  integer :: src_PE   ! The processor that is sending the data
  logical :: do_block ! If true add synchronizing barriers
  do_block = .false. ; if (present(blocking)) do_block = blocking
  if (present(from_PE)) then ; src_PE = from_PE ; else ; src_PE = root_PE() ; endif

  if (do_block) call mpp_sync(PElist)
  call mpp_broadcast(dat, src_PE, PElist)
  if (do_block) call mpp_sync_self(PElist)

end procedure broadcast_int32_0D
module procedure broadcast_int1D
  integer :: src_PE   ! The processor that is sending the data
  logical :: do_block ! If true add synchronizing barriers
  do_block = .false. ; if (present(blocking)) do_block = blocking
  if (present(from_PE)) then ; src_PE = from_PE ; else ; src_PE = root_PE() ; endif

  if (do_block) call mpp_sync(PElist)
  call mpp_broadcast(dat, length, src_PE, PElist)
  if (do_block) call mpp_sync_self(PElist)

end procedure broadcast_int1D
module procedure broadcast_real0D
  integer :: src_PE   ! The processor that is sending the data
  logical :: do_block ! If true add synchronizing barriers
  do_block = .false. ; if (present(blocking)) do_block = blocking
  if (present(from_PE)) then ; src_PE = from_PE ; else ; src_PE = root_PE() ; endif

  if (do_block) call mpp_sync(PElist)
  call mpp_broadcast(dat, src_PE, PElist)
  if (do_block) call mpp_sync_self(PElist)

end procedure broadcast_real0D
module procedure broadcast_real1D
  integer :: src_PE   ! The processor that is sending the data
  logical :: do_block ! If true add synchronizing barriers
  do_block = .false. ; if (present(blocking)) do_block = blocking
  if (present(from_PE)) then ; src_PE = from_PE ; else ; src_PE = root_PE() ; endif

  if (do_block) call mpp_sync(PElist)
  call mpp_broadcast(dat, length, src_PE, PElist)
  if (do_block) call mpp_sync_self(PElist)

end procedure broadcast_real1D
module procedure broadcast_real2D
  integer :: src_PE   ! The processor that is sending the data
  logical :: do_block ! If true add synchronizing barriers
  do_block = .false. ; if (present(blocking)) do_block = blocking
  if (present(from_PE)) then ; src_PE = from_PE ; else ; src_PE = root_PE() ; endif

  if (do_block) call mpp_sync(PElist)
  call mpp_broadcast(dat, length, src_PE, PElist)
  if (do_block) call mpp_sync_self(PElist)

end procedure broadcast_real2D
module procedure broadcast_real3D
  integer :: src_PE   ! The processor that is sending the data
  logical :: do_block ! If true add synchronizing barriers
  do_block = .false. ; if (present(blocking)) do_block = blocking
  if (present(from_PE)) then ; src_PE = from_PE ; else ; src_PE = root_PE() ; endif

  if (do_block) call mpp_sync(PElist)
  call mpp_broadcast(dat, length, src_PE, PElist)
  if (do_block) call mpp_sync_self(PElist)

end procedure broadcast_real3D
module procedure field_chksum_real_0d
  chksum = mpp_chksum(field, pelist, mask_val)
end procedure field_chksum_real_0d
module procedure field_chksum_real_1d
  chksum = mpp_chksum(field, pelist, mask_val)
end procedure field_chksum_real_1d
module procedure field_chksum_real_2d
  chksum = mpp_chksum(field, pelist, mask_val)
end procedure field_chksum_real_2d
module procedure field_chksum_real_3d
  chksum = mpp_chksum(field, pelist, mask_val)
end procedure field_chksum_real_3d
module procedure field_chksum_real_4d
  chksum = mpp_chksum(field, pelist, mask_val)
end procedure field_chksum_real_4d
module procedure sum_across_PEs_int4_0d
  call mpp_sum(field, pelist)
end procedure sum_across_PEs_int4_0d
module procedure sum_across_PEs_int4_1d
  call mpp_sum(field, length, pelist)
end procedure sum_across_PEs_int4_1d
module procedure sum_across_PEs_int8_0d
  call mpp_sum(field, pelist)
end procedure sum_across_PEs_int8_0d
module procedure sum_across_PEs_int8_1d
  call mpp_sum(field, length, pelist)
end procedure sum_across_PEs_int8_1d
module procedure sum_across_PEs_int8_2d
  call mpp_sum(field, length, pelist)
end procedure sum_across_PEs_int8_2d
module procedure sum_across_PEs_real_0d
  call mpp_sum(field, pelist)
end procedure sum_across_PEs_real_0d
module procedure sum_across_PEs_real_1d
  call mpp_sum(field, length, pelist)
end procedure sum_across_PEs_real_1d
module procedure sum_across_PEs_real_2d
  call mpp_sum(field, length, pelist)
end procedure sum_across_PEs_real_2d
module procedure max_across_PEs_int_0d
  call mpp_max(field, pelist)
end procedure max_across_PEs_int_0d
module procedure max_across_PEs_real_0d
  call mpp_max(field, pelist)
end procedure max_across_PEs_real_0d
module procedure max_across_PEs_real_1d
  call mpp_max(field, length, pelist)
end procedure max_across_PEs_real_1d
module procedure min_across_PEs_int_0d
  call mpp_min(field, pelist)
end procedure min_across_PEs_int_0d
module procedure min_across_PEs_real_0d
  call mpp_min(field, pelist)
end procedure min_across_PEs_real_0d
module procedure min_across_PEs_real_1d
  call mpp_min(field, length, pelist)
end procedure min_across_PEs_real_1d
module procedure any_across_PEs
  integer :: field_flag
  field_flag = 0
  if (field) field_flag = 1
  call max_across_PEs(field_flag, pelist)
  any_across_PEs = (field_flag > 0)
end procedure any_across_PEs
module procedure all_across_PEs
  integer :: field_flag
  field_flag = 0
  if (field) field_flag = 1
  call min_across_PEs(field_flag, pelist)
  all_across_PEs = (field_flag > 0)
end procedure all_across_PEs
module procedure MOM_infra_init
  call fms_init(localcomm)
end procedure MOM_infra_init
module procedure MOM_infra_end
  call print_memuse_stats( 'Memory HiWaterMark', always=.TRUE. )
  call fms_end()
end procedure MOM_infra_end
end submodule MOM_coms_infra_s
