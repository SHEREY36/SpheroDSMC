!> Simulator particles of a homogeneous gas of smooth spherocylinders.
!>
!> Each particle carries its velocity c, unit axis u, angular velocity w
!> (always perpendicular to u: the spin about the axis of a smooth rod never
!> changes and is zero) and its rotational energy E_rot = I |w|^2 / 2.
module particles
  use precision, only: dp, pi
  use random, only: normal, unit_vector
  implicit none
  private
  public :: gas_t, setup_gas, initialize, temperatures, rotate_axes, non_gaussianity, &
            shear_step, kinetic_stress, thermal_energy

  type :: gas_t
    integer :: n = 0                      ! simulator particles (one per real particle)
    real(dp) :: diameter = 1, length = 0  ! d and cylinder length L = (AR - 1) d
    real(dp) :: aspect_ratio = 1
    real(dp) :: mass = 1, inertia = 1     ! m and transverse moment of inertia I
    real(dp) :: volume = 1                ! cell volume V = N / n_density
    real(dp) :: number_density = 1
    real(dp) :: a_bar = 0, a_sup = 0      ! isotropic mean and maximum of A_perp
    real(dp), allocatable :: c(:, :), u(:, :), w(:, :), erot(:)
  end type gas_t

contains

  !> Geometry, mass and inertia of a uniform spherocylinder of material density rho.
  subroutine setup_gas(g, n, aspect_ratio, diameter, rho, volume_fraction)
    type(gas_t), intent(inout) :: g
    integer, intent(in) :: n
    real(dp), intent(in) :: aspect_ratio, diameter, rho, volume_fraction
    real(dp) :: d, l, vp
    g%n = n
    g%aspect_ratio = aspect_ratio
    g%diameter = diameter
    d = diameter
    l = (aspect_ratio - 1) * d
    g%length = l
    vp = pi * d**3 / 6 + pi * d**2 * l / 4
    g%mass = rho * vp
    g%inertia = pi * rho * d**2 * l**3 / 48 + 3 * pi * rho * d**4 * l / 64 &
              + pi * rho * d**5 / 60 + pi * rho * d**3 * l**2 / 24
    g%number_density = volume_fraction / vp
    g%volume = n / g%number_density
    g%a_bar = pi * (d**2 + d * l + l**2 / 8)
    g%a_sup = pi * d**2 + 4 * d * l + l**2
    allocate(g%c(3, n), g%u(3, n), g%w(3, n), g%erot(n))
  end subroutine setup_gas

  !> Maxwellian velocities at T_tr, isotropic axes, two-dimensional Maxwellian
  !> spins at T_rot.  With exact = .true. the realised temperatures are
  !> rescaled to the requested values (zero total momentum is kept).
  subroutine initialize(g, ttr, trot, exact)
    type(gas_t), intent(inout) :: g
    real(dp), intent(in) :: ttr, trot
    logical, intent(in) :: exact
    integer :: i
    real(dp) :: mean(3), t1, t2
    do i = 1, g%n
      g%c(:, i) = [normal(), normal(), normal()] * sqrt(ttr / g%mass)
      g%u(:, i) = unit_vector()
      g%w(:, i) = [normal(), normal(), normal()] * sqrt(trot / g%inertia)
      g%w(:, i) = g%w(:, i) - dot_product(g%w(:, i), g%u(:, i)) * g%u(:, i)
      g%erot(i) = 0.5_dp * g%inertia * dot_product(g%w(:, i), g%w(:, i))
    end do
    mean = sum(g%c, dim=2) / g%n
    do i = 1, 3
      g%c(i, :) = g%c(i, :) - mean(i)
    end do
    if (exact) then
      call temperatures(g, t1, t2)
      g%c = g%c * sqrt(ttr / t1)
      g%w = g%w * sqrt(trot / t2)
      g%erot = g%erot * (trot / t2)
    end if
  end subroutine initialize

  !> Translational (peculiar) and rotational temperatures.
  subroutine temperatures(g, ttr, trot)
    type(gas_t), intent(in) :: g
    real(dp), intent(out) :: ttr, trot
    real(dp) :: mean(3), s2(3)
    integer :: i
    mean = sum(g%c, dim=2) / g%n
    do i = 1, 3
      s2(i) = sum((g%c(i, :) - mean(i))**2)
    end do
    ttr = g%mass * sum(s2) / (3 * g%n)
    trot = sum(g%erot) / g%n
  end subroutine temperatures

  !> Torque-free motion over dt: each axis precesses about its fixed spin
  !> (exact Rodrigues rotation).
  subroutine rotate_axes(g, dt)
    type(gas_t), intent(inout) :: g
    real(dp), intent(in) :: dt
    integer :: i
    real(dp) :: s, k(3), u(3), ca, sa
    do i = 1, g%n
      s = norm2(g%w(:, i))
      if (s * dt <= 1.0e-14_dp) cycle
      k = g%w(:, i) / s
      u = g%u(:, i)
      ca = cos(s * dt)
      sa = sin(s * dt)
      u = u * ca + [k(2) * u(3) - k(3) * u(2), k(3) * u(1) - k(1) * u(3), &
                    k(1) * u(2) - k(2) * u(1)] * sa + k * dot_product(k, u) * (1 - ca)
      g%u(:, i) = u / norm2(u)
    end do
  end subroutine rotate_axes

  !> Kurtosis coefficients a2 of the translational and rotational distributions:
  !> a2_tr = <(m c^2)^2>/(15 T_tr^2) - 1 (three DOF), a2_rot = <E_rot^2>/(2 T_rot^2) - 1 (two DOF).
  subroutine non_gaussianity(g, ttr, trot, a2_tr, a2_rot)
    type(gas_t), intent(in) :: g
    real(dp), intent(in) :: ttr, trot
    real(dp), intent(out) :: a2_tr, a2_rot
    real(dp) :: mean(3), e
    integer :: i
    mean = sum(g%c, dim=2) / g%n
    a2_tr = 0
    do i = 1, g%n
      e = g%mass * sum((g%c(:, i) - mean)**2)
      a2_tr = a2_tr + e * e
    end do
    a2_tr = a2_tr / (g%n * 15 * ttr**2) - 1
    a2_rot = sum(g%erot**2) / (g%n * 2 * trot**2) - 1
  end subroutine non_gaussianity

  !> Uniform shear flow in the homogeneous (Lagrangian) frame: over a strain a = gdot dt
  !> the peculiar velocities follow dc_x/dt = -gdot c_y, i.e. c_x <- c_x - a c_y exactly
  !> (axes and spins are not affected).  work accumulates the exact change of the kinetic
  !> energy, m sum(-a c_x c_y + a^2 c_y^2 / 2).
  subroutine shear_step(g, a, work)
    type(gas_t), intent(inout) :: g
    real(dp), intent(in) :: a
    real(dp), intent(inout) :: work
    integer :: i
    real(dp) :: s
    s = 0
    do i = 1, g%n
      s = s - a * g%c(1, i) * g%c(2, i) + 0.5_dp * a * a * g%c(2, i)**2
      g%c(1, i) = g%c(1, i) - a * g%c(2, i)
    end do
    work = work + g%mass * s
  end subroutine shear_step

  !> Kinetic stress reduced by n T_tr: P^k_ij / (n T_tr) = m <c_i c_j> / T_tr, in the
  !> order xx, yy, zz, xy, xz, yz (peculiar velocities).
  subroutine kinetic_stress(g, ttr, p)
    type(gas_t), intent(in) :: g
    real(dp), intent(in) :: ttr
    real(dp), intent(out) :: p(6)
    real(dp) :: mean(3), d(3)
    integer :: i
    mean = sum(g%c, dim=2) / g%n
    p = 0
    do i = 1, g%n
      d = g%c(:, i) - mean
      p = p + [d(1) * d(1), d(2) * d(2), d(3) * d(3), d(1) * d(2), d(1) * d(3), d(2) * d(3)]
    end do
    p = g%mass * p / (g%n * ttr)
  end subroutine kinetic_stress

  !> Total kinetic plus rotational energy of the simulator particles.
  real(dp) function thermal_energy(g) result(e)
    type(gas_t), intent(in) :: g
    e = 0.5_dp * g%mass * sum(g%c**2) + sum(g%erot)
  end function thermal_energy

end module particles
