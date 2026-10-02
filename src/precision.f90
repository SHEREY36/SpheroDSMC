!> Kinds and constants shared by every module.
module precision
  use, intrinsic :: iso_fortran_env, only: real64, int64
  implicit none
  private
  public :: dp, i8, pi

  integer, parameter :: dp = real64
  integer, parameter :: i8 = int64
  real(dp), parameter :: pi = 3.141592653589793238462643383279502884_dp
end module precision
