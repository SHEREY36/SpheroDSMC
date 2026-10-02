!> Random variates on top of the intrinsic generator (xoshiro256** in gfortran).
!> Uniform numbers are drawn in blocks for speed.
module random
  use precision, only: dp, i8
  implicit none
  private
  public :: rng_seed, uniform, normal, gamma_variate, beta_variate, unit_vector, &
            unit_vector_perpendicular

  integer, parameter :: nblock = 8192
  real(dp), save :: block(nblock)
  integer, save :: next = nblock + 1
  logical, save :: have_spare = .false.
  real(dp), save :: spare = 0.0_dp

contains

  !> Seed the intrinsic generator reproducibly from one integer.
  subroutine rng_seed(seed)
    integer, intent(in) :: seed
    integer :: n, i
    integer, allocatable :: s(:)
    integer(i8) :: x
    call random_seed(size=n)
    allocate(s(n))
    x = modulo(int(seed, i8) * 7919_i8 + 104729_i8, 2147483647_i8)
    do i = 1, n
      x = modulo(48271_i8 * x + 11_i8, 2147483647_i8)       ! Park-Miller step
      s(i) = int(x)
    end do
    call random_seed(put=s)
    next = nblock + 1
    have_spare = .false.
  end subroutine rng_seed

  !> Uniform on [0, 1).
  real(dp) function uniform()
    if (next > nblock) then
      call random_number(block)
      next = 1
    end if
    uniform = block(next)
    next = next + 1
  end function uniform

  !> Standard normal (Marsaglia polar method).
  real(dp) function normal()
    real(dp) :: x, y, s
    if (have_spare) then
      have_spare = .false.
      normal = spare
      return
    end if
    do
      x = 2.0_dp * uniform() - 1.0_dp
      y = 2.0_dp * uniform() - 1.0_dp
      s = x * x + y * y
      if (s > 0.0_dp .and. s < 1.0_dp) exit
    end do
    s = sqrt(-2.0_dp * log(s) / s)
    spare = y * s
    have_spare = .true.
    normal = x * s
  end function normal

  !> Gamma(shape a, scale 1) for a >= 1 (Marsaglia and Tsang, 2000).
  real(dp) function gamma_variate(a)
    real(dp), intent(in) :: a
    real(dp) :: d, c, x, v, u
    d = a - 1.0_dp / 3.0_dp
    c = 1.0_dp / sqrt(9.0_dp * d)
    do
      x = normal()
      v = 1.0_dp + c * x
      if (v <= 0.0_dp) cycle
      v = v * v * v
      u = 1.0_dp - uniform()                                  ! (0, 1]
      if (u < 1.0_dp - 0.0331_dp * x**4) exit
      if (log(u) < 0.5_dp * x * x + d * (1.0_dp - v + log(v))) exit
    end do
    gamma_variate = d * v
  end function gamma_variate

  !> Beta(a, b) from two Gamma variates.
  real(dp) function beta_variate(a, b)
    real(dp), intent(in) :: a, b
    real(dp) :: x, y
    x = gamma_variate(a)
    y = gamma_variate(b)
    beta_variate = x / (x + y)
  end function beta_variate

  !> Unit vector uniform on the sphere.
  function unit_vector() result(e)
    real(dp) :: e(3)
    e = [normal(), normal(), normal()]
    e = e / norm2(e)
  end function unit_vector

  !> Unit vector uniform on the circle perpendicular to the unit vector u.
  function unit_vector_perpendicular(u) result(e)
    real(dp), intent(in) :: u(3)
    real(dp) :: e(3)
    do
      e = [normal(), normal(), normal()]
      e = e - dot_product(e, u) * u
      if (norm2(e) > 1.0e-12_dp) exit
    end do
    e = e / norm2(e)
  end function unit_vector_perpendicular

end module random
