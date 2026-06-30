! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> Implements the general purpose Artificial Neural Network (ANN).
module MOM_ANN

! This file is part of MOM6. See LICENSE.md for the license

use MOM_io, only : MOM_read_data, field_exists
use MOM_error_handler, only : MOM_error, FATAL, MOM_mesg
use numerical_testing_type, only : testing

implicit none ; private

!#include <MOM_memory.h>

public ANN_init, ANN_allocate, ANN_apply, ANN_end, ANN_unit_tests
public ANN_apply_vector_orig, ANN_apply_vector_oi, ANN_apply_array_sio
public set_layer, set_input_normalization, set_output_normalization
public ANN_random, randomize_layer

!> Applies ANN to x, returning results in y
interface ANN_apply
  module procedure ANN_apply_vector_oi
  module procedure ANN_apply_array_sio
end interface ANN_apply

!> Type for a single Linear layer of ANN,
!! i.e. stores the matrix A and bias b
!! for matrix-vector multiplication
!! y = A*x + b.
type, private :: layer_type ; private
  integer :: output_width        !< Number of rows in matrix A
  integer :: input_width         !< Number of columns in matrix A
  logical :: activation = .True. !< If true, apply the default activation function

  real, allocatable :: A(:,:) !< Matrix in column-major order
                              !! of size A(output_width, input_width) [nondim]
  real, allocatable :: b(:)   !< bias vector of size output_width [nondim]
end type layer_type

!> Control structure/type for ANN
type, public :: ANN_CS ; private
  ! Parameters
  integer :: num_layers          !< Number of layers in the ANN, including the input and output.
                                 !! For example, for ANN with one hidden layer, num_layers = 3.
  integer, allocatable &
          :: layer_sizes(:)      !< Array of length num_layers, storing the number of neurons in
                                 !! each layer.

  type(layer_type), allocatable &
          :: layers(:)           !< Array of length num_layers-1, where each element is the Linear
                                 !! transformation between layers defined by Matrix A and vias b.

  real, allocatable :: &
    input_means(:), &  !< Array of length layer_sizes(1) containing the mean of each input feature
                       !! prior to normalization by input_norms [arbitrary].
    input_norms(:), &  !< Array of length layer_sizes(1) containing the *inverse* of the standard
                       !! deviation for each input feature used to normalize (multiply) before
                       !! feeding into the ANN [arbitrary]
    output_means(:), & !< Array of length layer_sizes(num_layers) containing the mean of each
                       !! output prior to normalization by output_norms [arbitrary].
    output_norms(:)    !< Array of length layer_sizes(num_layers) containing the standard deviation
                       !! each output of the ANN will be multiplied [arbitrary]

  integer, public :: parameters = 0 !< Count of number of parameters
end type ANN_CS


  interface
module subroutine ANN_init(CS, NNfile)
  type(ANN_CS), intent(inout)  :: CS     !< ANN control structure.
  character(*), intent(in)     :: NNfile !< The name of NetCDF file having neural network parameters
  ! Local variables

end subroutine ANN_init
module subroutine ANN_allocate(CS, num_layers, layer_sizes)
  type(ANN_CS), intent(inout) :: CS !< ANN control structure
  integer,      intent(in)    :: num_layers !< The number of layers, including the input and output layer
  integer,      intent(in)    :: layer_sizes(num_layers) !< The number of neurons in each layer
  ! Local variables

  ! Assert that there is always an input and output layer
end subroutine ANN_allocate
module subroutine ANN_test(CS, NNfile)
  type(ANN_CS), intent(inout) :: CS     !< ANN control structure.
  character(*), intent(in)    :: NNfile !< The name of NetCDF file having neural network parameters
  ! Local variables

  ! Allocate data
end subroutine ANN_test
module subroutine ANN_end(CS)
  type(ANN_CS), intent(inout) :: CS !< ANN control structure.
  ! Local variables

end subroutine ANN_end
pure elemental module function activation_fn(x) result (y)
  real, intent(in) :: x !< Scalar input value [nondim]
  real             :: y !< Scalar output value [nondim]

end function activation_fn
module subroutine ANN_apply_vector_orig(x, y, CS)
  type(ANN_CS), intent(in)    :: CS                               !< ANN instance
  real,         intent(in)    :: x(CS%layer_sizes(1))             !< Inputs [arbitrary]
  real,         intent(inout) :: y(CS%layer_sizes(CS%num_layers)) !< Outputs [arbitrary]
  ! Local variables

  ! Normalize input
end subroutine ANN_apply_vector_orig
module subroutine ANN_apply_vector_oi(x, y, CS)
  type(ANN_CS), intent(in)    :: CS                               !< ANN instance
  real,         intent(in)    :: x(CS%layer_sizes(1))             !< Inputs [arbitrary]
  real,         intent(inout) :: y(CS%layer_sizes(CS%num_layers)) !< Outputs [arbitrary]
  ! Local variables

end subroutine ANN_apply_vector_oi
module subroutine ANN_apply_array_sio(nij, x, y, CS)
  type(ANN_CS), intent(in)    :: CS !< ANN control structure
  integer,      intent(in)    :: nij !< Size of spatial dimension
  real,         intent(in)    :: x(nij, CS%layer_sizes(1)) !< input [arbitrary]
  real,         intent(inout) :: y(nij, CS%layer_sizes(CS%num_layers)) !< output [arbitrary]
  ! Local variables

end subroutine ANN_apply_array_sio
module subroutine set_layer(ANN, layer, weights, biases, activation)
  type(ANN_CS), intent(inout) :: ANN !< ANN control structure
  integer,      intent(in)    :: layer !< The number of the layer being adjusted
  real,         intent(in)    :: weights(:,:) !< The weights to assign
  real,         intent(in)    :: biases(:) !< The biases to assign
  logical,      intent(in)    :: activation !< Turn on the activation function

end subroutine set_layer
module subroutine set_input_normalization(ANN, means, norms)
  type(ANN_CS),   intent(inout) :: ANN !< ANN control structure
  real, optional, intent(in)    :: means(:) !< The mean of each input
  real, optional, intent(in)    :: norms(:) !< The standard deviation of each input

end subroutine set_input_normalization
module subroutine set_output_normalization(ANN, means, norms)
  type(ANN_CS),   intent(inout) :: ANN !< ANN control structure
  real, optional, intent(in)    :: means(:) !< The mean of each output
  real, optional, intent(in)    :: norms(:) !< The standard deviation of each output

end subroutine set_output_normalization
module subroutine ANN_random(ANN, nlayers, widths)
  type(ANN_CS), intent(inout) :: ANN !< ANN control structure
  integer,      intent(in)    :: nlayers !< Number of layers
  integer,      intent(in)    :: widths(nlayers) !< Width of each layer
  ! Local variables

end subroutine ANN_random
module subroutine randomize_layer(ANN, nlayers, layer, widths)
  type(ANN_CS), intent(inout) :: ANN !< ANN control structure
  integer,      intent(in)    :: nlayers !< Number of layers
  integer,      intent(in)    :: layer !< Layer number to randomize
  integer,      intent(in)    :: widths(nlayers) !< Width of each layer
  ! Local variables

end subroutine randomize_layer
logical module function ANN_unit_tests(verbose)
  logical, intent(in) :: verbose !< If true, write results to stdout
  ! Local variables

end function ANN_unit_tests
  end interface

end module MOM_ANN
