submodule (MOM_spatial_means) MOM_spatial_means_s
#include <MOM_memory.h>
  implicit none
contains
module procedure global_area_mean
  real :: tmpForSumming(SZI_(G),SZJ_(G))  ! An unscaled cell integral in [a L2 ~> a m2] or a
  real :: scalefac   ! A scaling factor for the variable that is not reversed [a A-1 ~> 1] or [B A-1 ~> b a-1] or [1]
  real :: temp_scale ! A temporary scaling factor [a A-1 ~> 1] or [b B-1 ~> 1] or [1]
  integer :: i, j, is, ie, js, je
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec

  temp_scale = 1.0 ; if (present(tmp_scale)) temp_scale = tmp_scale

  scalefac = 1.0
  if (present(unscale)) then ; scalefac = unscale
  elseif (present(scale)) then ; scalefac = scale ; endif

  tmpForSumming(:,:) = 0.
  do j=js,je ; do i=is,ie
    tmpForSumming(i,j) = var(i,j) * (scalefac * G%areaT(i,j) * G%mask2dT(i,j))
  enddo ; enddo

  global_area_mean = reproducing_sum(tmpForSumming, unscale=temp_scale*G%US%L_to_m**2) * G%IareaT_global

end procedure global_area_mean
module procedure global_area_mean_v
  real, dimension(SZI_(G),SZJ_(G)) :: tmpForSumming ! An unscaled cell integral [A L2 ~> a m2]
  real :: temp_scale ! A temporary scaling factor [a A-1 ~> 1] or [1]
  integer :: i, j, is, ie, js, je, isB, ieB, jsB, jeB
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec
  isB = G%iscB ; ieB = G%iecB ; jsB = G%jscB ; jeB = G%jecB

  temp_scale = 1.0 ; if (present(tmp_scale)) temp_scale = tmp_scale

  tmpForSumming(:,:) = 0.
  do j=js,je ; do i=is,ie
    tmpForSumming(i,j) = G%areaT(i,j) * &
             (var(i,J) * G%mask2dCv(i,J) + var(i,J-1) * G%mask2dCv(i,J-1)) / &
             max(1.e-20, G%mask2dCv(i,J)+G%mask2dCv(i,J-1))
  enddo ; enddo
  global_area_mean_v = reproducing_sum(tmpForSumming, unscale=G%US%L_to_m**2*temp_scale) * G%IareaT_global

end procedure global_area_mean_v
module procedure global_area_mean_u
  real, dimension(SZI_(G),SZJ_(G)) :: tmpForSumming ! An unscaled cell integral [A L2 ~> a m2]
  real :: temp_scale ! A temporary scaling factor [a A-1 ~> 1] or [1]
  integer :: i, j, is, ie, js, je, isB, ieB, jsB, jeB
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec
  isB = G%iscB ; ieB = G%iecB ; jsB = G%jscB ; jeB = G%jecB

  temp_scale = 1.0 ; if (present(tmp_scale)) temp_scale = tmp_scale

  tmpForSumming(:,:) = 0.
  do j=js,je ; do i=is,ie
    tmpForSumming(i,j) = G%areaT(i,j) * &
             (var(I,j) * G%mask2dCu(I,j) + var(I-1,j) * G%mask2dCu(I-1,j)) / &
             max(1.e-20, G%mask2dCu(I,j)+G%mask2dCu(I-1,j))
  enddo ; enddo
  global_area_mean_u = reproducing_sum(tmpForSumming, unscale=G%US%L_to_m**2*temp_scale) * G%IareaT_global

end procedure global_area_mean_u
module procedure global_area_integral
  real, dimension(SZI_(G),SZJ_(G)) :: tmpForSumming ! An unscaled cell integral in [a m2] or
  real :: scalefac  ! An overall scaling factor for the areas and variable, in units of [a m2 A-1 L-2 ~> 1]
  real :: temp_scale ! A temporary scaling factor [a m2 L-2 A-1 ~> 1] or [b m2 L-2 B-1 ~> 1] or [1]
  integer :: i, j, is, ie, js, je
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec

  if (present(tmp_scale)) then
    temp_scale = G%US%L_to_m**2 * tmp_scale  ! Units of [a m2 A-1 L-2 ~> 1] or [b m2 B-1 L-2 ~> 1]
    scalefac = 1.0
  else
    temp_scale = 1.0
    scalefac = G%US%L_to_m**2
  endif
  if (present(unscale)) then ; scalefac = scalefac * unscale
  elseif (present(scale)) then ; scalefac = scalefac * scale ; endif

  tmpForSumming(:,:) = 0.
  if (present(area)) then
    do j=js,je ; do i=is,ie
      tmpForSumming(i,j) = var(i,j) * (scalefac * area(i,j))
    enddo ; enddo
  else
    do j=js,je ; do i=is,ie
      tmpForSumming(i,j) = var(i,j) * (scalefac * G%areaT(i,j) * G%mask2dT(i,j))
    enddo ; enddo
  endif

  global_area_integral = reproducing_sum(tmpForSumming, unscale=temp_scale)

end procedure global_area_integral
module procedure global_layer_mean
  real, dimension(SZK_(GV)) :: global_layer_mean  !< The mean of the variable in the arbitrary scaled [A ~> a]
  real :: tmpForSumming(G%isc:G%iec,G%jsc:G%jec,SZK_(GV)) ! An unscaled cell integral in [L2 a m ~> a m3] or
  real :: weight(G%isc:G%iec,G%jsc:G%jec,SZK_(GV))  ! The volume or mass of each cell, depending on whether
  real :: scalefac   ! A scaling factor for the variable [a A-1 ~> 1] or [B A-1 ~> b a-1] or [1]
  type(EFP_type) :: laysums(2*SZK_(GV)) ! A vector of sums with heterogeneous meanings, with the first
  real :: global_temp_scalar   ! The global integral of the tracer over all
  real :: global_weight_scalar ! The global integral of the volume or mass over all
  real :: temp_scale ! A temporary scaling factor [a A-1 ~> 1] or [b B-1 ~> 1] or [1]
  integer :: i, j, k, is, ie, js, je, nz
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke

  temp_scale = 1.0 ; if (present(tmp_scale)) temp_scale = tmp_scale

  scalefac = 1.0
  if (present(unscale)) then ; scalefac = unscale
  elseif (present(scale)) then ; scalefac = scale ; endif
  tmpForSumming(:,:,:) = 0. ; weight(:,:,:) = 0.

  do k=1,nz ; do j=js,je ; do i=is,ie
    weight(i,j,k)  =  (GV%H_to_MKS * h(i,j,k)) * (G%areaT(i,j) * G%mask2dT(i,j))
    tmpForSumming(i,j,k) =  scalefac * var(i,j,k) * weight(i,j,k)
  enddo ; enddo ; enddo

  global_temp_scalar = reproducing_sum(tmpForSumming, EFP_lay_sums=laysums(1:nz), only_on_PE=.true., &
                                       unscale=temp_scale*G%US%L_to_m**2)
  global_weight_scalar = reproducing_sum(weight, EFP_lay_sums=laysums(nz+1:2*nz), only_on_PE=.true., &
                                         unscale=G%US%L_to_m**2)
  call EFP_sum_across_PEs(laysums, 2*nz)

  ! Note that temp_scale appears in the denominator here because the variables returned via the
  ! EFP_lay_sums arguments to reproducing sums stay in unscaled mks units.
  do k=1,nz
    global_layer_mean(k) = EFP_to_real(laysums(k)) / (temp_scale*EFP_to_real(laysums(nz+k)))
  enddo

end procedure global_layer_mean
module procedure global_volume_mean
  real :: temp_scale ! A temporary scaling factor [a A-1 ~> 1] or [b B-1 ~> 1] or [1]
  real :: scalefac   ! A scaling factor for the variable [a A-1 ~> 1] or [B A-1 ~> b a-1] or [1]
  real :: weight_here ! The volume or mass of a grid cell [L2 m ~> m3] or [L2 kg m-2 ~> kg]
  real, dimension(SZI_(G),SZJ_(G)) :: tmpForSumming ! The volume or mass integral of the variable in a column
  real, dimension(SZI_(G),SZJ_(G)) :: sum_weight  ! The volume or mass of each column of water
  integer :: i, j, k, is, ie, js, je, nz
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke

  temp_scale = 1.0 ; if (present(tmp_scale)) temp_scale = tmp_scale

  scalefac = 1.0
  if (present(unscale)) then ; scalefac = unscale
  elseif (present(scale)) then ; scalefac = scale ; endif
  tmpForSumming(:,:) = 0. ; sum_weight(:,:) = 0.

  do k=1,nz ; do j=js,je ; do i=is,ie
    weight_here  =  (GV%H_to_MKS * h(i,j,k)) * (G%areaT(i,j) * G%mask2dT(i,j))
    tmpForSumming(i,j) = tmpForSumming(i,j) + scalefac * var(i,j,k) * weight_here
    sum_weight(i,j) = sum_weight(i,j) + weight_here
  enddo ; enddo ; enddo
  global_volume_mean = (reproducing_sum(tmpForSumming, unscale=temp_scale*G%US%L_to_m**2)) / &
                       (reproducing_sum(sum_weight, unscale=G%US%L_to_m**2))

end procedure global_volume_mean
module procedure global_mass_integral
  real :: tmpForSumming(SZI_(G),SZJ_(G)) ! The mass-weighted integral of the variable in a column in
  real :: scalefac   ! An overall scaling factor for the cell mass and variable in [a kg A-1 R-1 Z-1 L-2 ~> 1]
  real :: temp_scale ! A temporary scaling factor [1] or if tmp_scale is present this could be in
  logical :: global_sum ! If true do the sum globally, but if false only do the sum on the current PE.
  integer :: i, j, k, is, ie, js, je, nz
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke

  if (present(tmp_scale)) then
    temp_scale = G%US%RZL2_to_kg * tmp_scale
    if (.not.present(var)) temp_scale = G%US%RZL2_to_kg
    scalefac = 1.0
  else
    temp_scale = 1.0
    scalefac = G%US%RZL2_to_kg
  endif
  if (present(var)) then
    if (present(unscale)) then ; scalefac = scalefac * unscale
    elseif (present(scale)) then ; scalefac = scalefac * scale ; endif
  endif

  tmpForSumming(:,:) = 0.0
  if (present(var)) then
    do k=1,nz ; do j=js,je ; do i=is,ie
      tmpForSumming(i,j) = tmpForSumming(i,j) + var(i,j,k) * &
                ((GV%H_to_RZ * h(i,j,k)) * (scalefac*G%areaT(i,j) * G%mask2dT(i,j)))
    enddo ; enddo ; enddo
  else
    do k=1,nz ; do j=js,je ; do i=is,ie
      tmpForSumming(i,j) = tmpForSumming(i,j) + &
                ((GV%H_to_RZ * h(i,j,k)) * (scalefac*G%areaT(i,j) * G%mask2dT(i,j)))
    enddo ; enddo ; enddo
  endif
  global_sum = .true. ; if (present(on_PE_only)) global_sum = .not.on_PE_only
  if (global_sum) then
    global_mass_integral = reproducing_sum(tmpForSumming, unscale=temp_scale)
  else
    global_mass_integral = 0.0
    do j=js,je ; do i=is,ie
      global_mass_integral = global_mass_integral + tmpForSumming(i,j)
    enddo ; enddo
  endif

end procedure global_mass_integral
module procedure global_mass_int_EFP
  real :: tmpForSum(SZI_(G),SZJ_(G)) ! The mass-weighted integral of the variable in a column [kg a] or [kg]
  real :: scalefac  ! An overall scaling factor for the cell mass and variable [a kg A-1 H-1 L-2 ~> kg m-3 or 1]
  integer :: i, j, k, is, ie, js, je, nz, isr, ier, jsr, jer
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke
  isr = is - (G%isd-1) ; ier = ie - (G%isd-1) ; jsr = js - (G%jsd-1) ; jer = je - (G%jsd-1)

  scalefac = GV%H_to_kg_m2 * G%US%L_to_m**2
  if (present(unscale)) then ; scalefac = unscale * scalefac
  elseif (present(scale)) then ; scalefac = scale * scalefac ; endif

  tmpForSum(:,:) = 0.0
  if (present(var)) then
    do k=1,nz ; do j=js,je ; do i=is,ie
      tmpForSum(i,j) = tmpForSum(i,j) + var(i,j,k) * &
                ((scalefac * h(i,j,k)) * (G%areaT(i,j) * G%mask2dT(i,j)))
    enddo ; enddo ; enddo
  else
    do k=1,nz ; do j=js,je ; do i=is,ie
      tmpForSum(i,j) = tmpForSum(i,j) + &
                ((scalefac * h(i,j,k)) * (G%areaT(i,j) * G%mask2dT(i,j)))
    enddo ; enddo ; enddo
  endif

  global_mass_int_EFP = reproducing_sum_EFP(tmpForSum, isr, ier, jsr, jer, only_on_PE=on_PE_only)

end procedure global_mass_int_EFP
module procedure global_i_mean
  type(EFP_type), allocatable, dimension(:) :: asum      ! The masked sum of the variable in each row [a]
  type(EFP_type), allocatable, dimension(:) :: mask_sum  ! The sum of the mask values in each row [nondim]
  real :: scalefac   ! A scaling factor for the variable [a A-1 ~> 1]
  real :: rescale    ! A factor for redoing any internal rescaling before output [A a-1 ~> 1]
  real :: mask_sum_r ! The sum of the mask values in a row [nondim]
  integer :: is, ie, js, je, idg_off, jdg_off
  integer :: i, j
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec
  idg_off = G%idg_offset ; jdg_off = G%jdg_offset

  scalefac = 1.0
  if (present(unscale)) then ; scalefac = unscale
  elseif (present(scale)) then ; scalefac = scale ; endif

  rescale = 1.0
  if (present(tmp_scale)) then ; if (tmp_scale /= 0.0) then
    scalefac = scalefac * tmp_scale ; rescale = 1.0 / tmp_scale
  endif ; endif
  call reset_EFP_overflow_error()

  allocate(asum(G%jsg:G%jeg))
  if (present(mask)) then
    allocate(mask_sum(G%jsg:G%jeg))

    do j=G%jsg,G%jeg
      asum(j) = real_to_EFP(0.0) ; mask_sum(j) = real_to_EFP(0.0)
    enddo

    do j=js,je ; do i=is,ie
      asum(j+jdg_off) = asum(j+jdg_off) + real_to_EFP(scalefac*array(i,j)*mask(i,j))
      mask_sum(j+jdg_off) = mask_sum(j+jdg_off) + real_to_EFP(mask(i,j))
    enddo ; enddo

    if (query_EFP_overflow_error()) call MOM_error(FATAL, &
      "global_i_mean overflow error occurred before sums across PEs.")

    call EFP_sum_across_PEs(asum(G%jsg:G%jeg), G%jeg-G%jsg+1)
    call EFP_sum_across_PEs(mask_sum(G%jsg:G%jeg), G%jeg-G%jsg+1)

    if (query_EFP_overflow_error()) call MOM_error(FATAL, &
      "global_i_mean overflow error occurred during sums across PEs.")

    do j=js,je
      mask_sum_r = EFP_to_real(mask_sum(j+jdg_off))
      if (mask_sum_r == 0.0 ) then ; i_mean(j) = 0.0 ; else
        i_mean(j) = EFP_to_real(asum(j+jdg_off)) / mask_sum_r
      endif
    enddo

    deallocate(mask_sum)
  else
    do j=G%jsg,G%jeg ; asum(j) = real_to_EFP(0.0) ; enddo

    do j=js,je ; do i=is,ie
      asum(j+jdg_off) = asum(j+jdg_off) + real_to_EFP(scalefac*array(i,j))
    enddo ; enddo

    if (query_EFP_overflow_error()) call MOM_error(FATAL, &
      "global_i_mean overflow error occurred before sum across PEs.")

    call EFP_sum_across_PEs(asum(G%jsg:G%jeg), G%jeg-G%jsg+1)

    if (query_EFP_overflow_error()) call MOM_error(FATAL, &
      "global_i_mean overflow error occurred during sum across PEs.")

    do j=js,je
      i_mean(j) = EFP_to_real(asum(j+jdg_off)) / real(G%ieg-G%isg+1)
    enddo
  endif

  if (rescale /= 1.0) then ; do j=js,je ; i_mean(j) = rescale*i_mean(j) ; enddo ; endif

  deallocate(asum)

end procedure global_i_mean
module procedure global_j_mean
  type(EFP_type), allocatable, dimension(:) :: asum      ! The masked sum of the variable in each row [a]
  type(EFP_type), allocatable, dimension(:) :: mask_sum  ! The sum of the mask values in each row [nondim]
  real :: mask_sum_r ! The sum of the mask values in a row [nondim]
  real :: scalefac   ! A scaling factor for the variable [a A-1 ~> 1]
  real :: rescale    ! A factor for redoing any internal rescaling before output [A a-1 ~> 1]
  integer :: is, ie, js, je, idg_off, jdg_off
  integer :: i, j
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec
  idg_off = G%idg_offset ; jdg_off = G%jdg_offset

  scalefac = 1.0
  if (present(unscale)) then ; scalefac = unscale
  elseif (present(scale)) then ; scalefac = scale ; endif

  rescale = 1.0
  if (present(tmp_scale)) then ; if (tmp_scale /= 0.0) then
    scalefac = scalefac * tmp_scale ; rescale = 1.0 / tmp_scale
  endif ; endif
  call reset_EFP_overflow_error()

  allocate(asum(G%isg:G%ieg))
  if (present(mask)) then
    allocate (mask_sum(G%isg:G%ieg))

    do i=G%isg,G%ieg
      asum(i) = real_to_EFP(0.0) ; mask_sum(i) = real_to_EFP(0.0)
    enddo

    do i=is,ie ; do j=js,je
      asum(i+idg_off) = asum(i+idg_off) + real_to_EFP(scalefac*array(i,j)*mask(i,j))
      mask_sum(i+idg_off) = mask_sum(i+idg_off) + real_to_EFP(mask(i,j))
    enddo ; enddo

    if (query_EFP_overflow_error()) call MOM_error(FATAL, &
      "global_j_mean overflow error occurred before sums across PEs.")

    call EFP_sum_across_PEs(asum(G%isg:G%ieg), G%ieg-G%isg+1)
    call EFP_sum_across_PEs(mask_sum(G%isg:G%ieg), G%ieg-G%isg+1)

    if (query_EFP_overflow_error()) call MOM_error(FATAL, &
      "global_j_mean overflow error occurred during sums across PEs.")

    do i=is,ie
      mask_sum_r = EFP_to_real(mask_sum(i+idg_off))
      if (mask_sum_r == 0.0 ) then ; j_mean(i) = 0.0 ; else
        j_mean(i) = EFP_to_real(asum(i+idg_off)) / mask_sum_r
      endif
    enddo

    deallocate(mask_sum)
  else
    do i=G%isg,G%ieg ; asum(i) = real_to_EFP(0.0) ; enddo

    do i=is,ie ; do j=js,je
      asum(i+idg_off) = asum(i+idg_off) + real_to_EFP(scalefac*array(i,j))
    enddo ; enddo

    if (query_EFP_overflow_error()) call MOM_error(FATAL, &
      "global_j_mean overflow error occurred before sum across PEs.")

    call EFP_sum_across_PEs(asum(G%isg:G%ieg), G%ieg-G%isg+1)

    if (query_EFP_overflow_error()) call MOM_error(FATAL, &
      "global_j_mean overflow error occurred during sum across PEs.")

    do i=is,ie
      j_mean(i) = EFP_to_real(asum(i+idg_off)) / real(G%jeg-G%jsg+1)
    enddo
  endif

  if (rescale /= 1.0) then ; do i=is,ie ; j_mean(i) = rescale*j_mean(i) ; enddo ; endif

  deallocate(asum)

end procedure global_j_mean
module procedure adjust_area_mean_to_zero
  real :: posVals(G%isc:G%iec, G%jsc:G%jec) ! The positive values in a cell or 0 [A ~> a]
  real :: negVals(G%isc:G%iec, G%jsc:G%jec) ! The negative values in a cell or 0 [A ~> a]
  real :: areaXposVals(G%isc:G%iec, G%jsc:G%jec) ! The cell area integral of the positive values [L2 A ~> m2 a]
  real :: areaXnegVals(G%isc:G%iec, G%jsc:G%jec) ! The cell area integral of the negative values [L2 A ~> m2 a]
  type(EFP_type), dimension(2) :: areaInt_EFP ! An EFP version integral of the values on the current PE [m2 a]
  real :: scalefac  ! A scaling factor for the variable [a A-1 ~> 1]
  real :: areaIntPosVals, areaIntNegVals ! The global area integral of the positive and negative values [m2 a]
  real :: posScale, negScale ! The scaling factor to apply to positive or negative values [nondim]
  integer :: i,j
  scalefac = 1.0
  if (present(unscale)) then ; scalefac = unscale
  elseif (present(unit_scale)) then ; scalefac = unit_scale ; endif

  ! areaXposVals(:,:) = 0.  ! This zeros out halo points.
  ! areaXnegVals(:,:) = 0.  ! This zeros out halo points.

  do j=G%jsc,G%jec ; do i=G%isc,G%iec
    posVals(i,j) = max(0., array(i,j))
    areaXposVals(i,j) = G%areaT(i,j) * posVals(i,j)
    negVals(i,j) = min(0., array(i,j))
    areaXnegVals(i,j) = G%areaT(i,j) * negVals(i,j)
  enddo ; enddo

  ! Combining the sums like this avoids separate blocking global sums.
  areaInt_EFP(1) = reproducing_sum_EFP( areaXposVals, only_on_PE=.true., unscale=scalefac*G%US%L_to_m**2 )
  areaInt_EFP(2) = reproducing_sum_EFP( areaXnegVals, only_on_PE=.true., unscale=scalefac*G%US%L_to_m**2 )
  call EFP_sum_across_PEs(areaInt_EFP, 2)
  areaIntPosVals = EFP_to_real( areaInt_EFP(1) )
  areaIntNegVals = EFP_to_real( areaInt_EFP(2) )

  posScale = 0.0 ; negScale = 0.0
  if ((areaIntPosVals>0.).and.(areaIntNegVals<0.)) then ! Only adjust if possible
    if (areaIntPosVals>-areaIntNegVals) then ! Scale down positive values
      posScale = - areaIntNegVals / areaIntPosVals
      do j=G%jsc,G%jec ; do i=G%isc,G%iec
        array(i,j) = (posScale * posVals(i,j)) + negVals(i,j)
      enddo ; enddo
    elseif (areaIntPosVals<-areaIntNegVals) then ! Scale down negative values
      negScale = - areaIntPosVals / areaIntNegVals
      do j=G%jsc,G%jec ; do i=G%isc,G%iec
        array(i,j) = posVals(i,j) + (negScale * negVals(i,j))
      enddo ; enddo
    endif
  endif
  if (present(scaling)) scaling = posScale - negScale

end procedure adjust_area_mean_to_zero
module procedure array_global_min_max
  real    :: tmax, tmin      ! Maximum and minimum tracer values, in the same units as tr_array [CU ~> conc]
  integer :: ijk_min_max(2)  ! Integers encoding the global grid positions of the global minimum and maximum values
  real    :: xyz_min_max(6)  ! A single array with the x-, y- and z-positions of the minimum and
  logical :: valid_PE        ! True if there are any valid points on the local PE.
  logical :: find_location   ! If true, report the locations of the extrema
  integer :: ijk_loc_max     ! An integer encoding the global grid position of the maximum tracer value on this PE
  integer :: ijk_loc_min     ! An integer encoding the global grid position of the minimum tracer value on this PE
  integer :: ijk_loc_here    ! An integer encoding the global grid position of the current grid point
  integer :: itmax, jtmax, ktmax, itmin, jtmin, ktmin
  integer :: i, j, k, isc, iec, jsc, jec
  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec

  find_location = (present(xgmin) .or. present(ygmin) .or. present(zgmin) .or. &
                   present(xgmax) .or. present(ygmax) .or. present(zgmax))

  ! The initial values set here are never used if there are any valid points.
  tmax = -huge(tmax) ; tmin = huge(tmin)

  if (find_location) then
    ! Find the maximum and minimum tracer values on this PE and their locations.
    valid_PE = .false.
    itmax = 0 ; jtmax = 0 ; ktmax = 0 ; ijk_loc_max = 0
    itmin = 0 ; jtmin = 0 ; ktmin = 0 ; ijk_loc_min = 0
    do k=1,nk ; do j=jsc,jec ; do i=isc,iec ; if (G%mask2dT(i,j) > 0.0) then
      valid_PE = .true.
      if (tr_array(i,j,k) > tmax) then
        tmax = tr_array(i,j,k)
        itmax = i ; jtmax = j ; ktmax = k
        ijk_loc_max = ijk_loc(i, j, k, nk, G%HI)
      elseif ((tr_array(i,j,k) == tmax) .and. (k <= ktmax)) then
        ijk_loc_here = ijk_loc(i, j, k, nk, G%HI)
        if (ijk_loc_here > ijk_loc_max) then
          itmax = i ; jtmax = j ; ktmax = k
          ijk_loc_max = ijk_loc_here
        endif
      endif
      if (tr_array(i,j,k) < tmin) then
        tmin = tr_array(i,j,k)
        itmin = i ; jtmin = j ; ktmin = k
        ijk_loc_min = ijk_loc(i, j, k, nk, G%HI)
      elseif ((tr_array(i,j,k) == tmin) .and. (k <= ktmin)) then
        ijk_loc_here = ijk_loc(i, j, k, nk, G%HI)
        if (ijk_loc_here > ijk_loc_min) then
          itmin = i ; jtmin = j ; ktmin = k
          ijk_loc_min = ijk_loc_here
        endif
      endif
    endif ; enddo ; enddo ; enddo
  else
    ! Only the maximum and minimum values are needed, and not their positions.
    do k=1,nk ; do j=jsc,jec ; do i=isc,iec ; if (G%mask2dT(i,j) > 0.0) then
      if (tr_array(i,j,k) > tmax) tmax = tr_array(i,j,k)
      if (tr_array(i,j,k) < tmin) tmin = tr_array(i,j,k)
    endif ; enddo ; enddo ; enddo
  endif

  ! Find the global maximum and minimum tracer values.
  g_max = tmax ; g_min = tmin
  call max_across_PEs(g_max)
  call min_across_PEs(g_min)

  if (find_location) then
    if (g_max < g_min) then
      ! This only occurs if there are no unmasked points anywhere in the domain.
      xyz_min_max(:) = 0.0
    else
      ! Find the global indices of the maximum and minimum locations.  This can
      ! occur on multiple PEs.
      ijk_min_max(1:2) = 0
      if (valid_PE) then
        if (g_min == tmin) ijk_min_max(1) = ijk_loc_min
        if (g_max == tmax) ijk_min_max(2) = ijk_loc_max
      endif
      ! If MOM6 supported taking maxima on arrays of integers, these could be combined as:
      ! call max_across_PEs(ijk_min_max, 2)
      call max_across_PEs(ijk_min_max(1))
      call max_across_PEs(ijk_min_max(2))

      ! Set the positions of the extrema if they occur on this PE.  This will only
      ! occur on a single PE.
      xyz_min_max(1:6) = -huge(xyz_min_max)  ! These huge negative values are never selected by max_across_PEs.
      if (valid_PE) then
        if (ijk_min_max(1) == ijk_loc_min) then
          xyz_min_max(1) = G%geoLonT(itmin,jtmin)
          xyz_min_max(2) = G%geoLatT(itmin,jtmin)
          xyz_min_max(3) = real(ktmin)
        endif
        if (ijk_min_max(2) == ijk_loc_max) then
          xyz_min_max(4) = G%geoLonT(itmax,jtmax)
          xyz_min_max(5) = G%geoLatT(itmax,jtmax)
          xyz_min_max(6) = real(ktmax)
        endif
      endif

      call max_across_PEs(xyz_min_max, 6)
    endif

    if (present(xgmin)) xgmin = xyz_min_max(1)
    if (present(ygmin)) ygmin = xyz_min_max(2)
    if (present(zgmin)) zgmin = xyz_min_max(3)
    if (present(xgmax)) xgmax = xyz_min_max(4)
    if (present(ygmax)) ygmax = xyz_min_max(5)
    if (present(zgmax)) zgmax = xyz_min_max(6)
  endif

  if (g_max < g_min) then
    ! There are no unmasked points anywhere in the domain.
    g_max = 0.0 ; g_min = 0.0
  endif

  if (present(unscale)) then
    ! Rescale g_min and g_max, perhaps changing their units from [CU ~> conc] to [conc]
    g_max = unscale * g_max
    g_min = unscale * g_min
  endif

end procedure array_global_min_max
module procedure ijk_loc
  integer :: ig, jg  ! Global index values with a global computational domain start value of 1.
  integer :: ij_loc  ! The encoding of the horizontal position
  integer :: qturns  ! The number of counter-clockwise quarter turns of the grid that have to be undone
  ig = i + HI%idg_offset + (1 - HI%isg)
  jg = j + HI%jdg_offset + (1 - HI%jsg)

  ! Compensate for the rotation of the model grid to give a rotationally invariant encoding.
  qturns = modulo(HI%turns, 4)
  if (qturns == 0) then
    ij_loc = ig + HI%niglobal * jg
  elseif (qturns == 1) then
    ij_loc = jg + HI%njglobal * ((HI%niglobal+1)-ig)
  elseif (qturns == 2) then
    ij_loc = ((HI%niglobal+1)-ig) + HI%niglobal * ((HI%njglobal+1)-jg)
  elseif (qturns == 3) then
    ij_loc = ((HI%njglobal+1)-jg) + HI%njglobal * ig
  endif

  ijk_loc = ij_loc + (HI%niglobal*HI%njglobal) * (nk-k)

end procedure ijk_loc
end submodule MOM_spatial_means_s
