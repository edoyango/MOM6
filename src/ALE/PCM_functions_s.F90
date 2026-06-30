submodule (PCM_functions) PCM_functions_s
  implicit none
contains
module procedure PCM_reconstruction
  integer :: k
  ppoly_coef(:,1) = u(:)

  ! The edge values are equal to the cell average
  do k = 1,N
    edge_values(k,:) = u(k)
  enddo

end procedure PCM_reconstruction
end submodule PCM_functions_s
