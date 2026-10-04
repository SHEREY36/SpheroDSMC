!> One collision step of the DSMC in one cell (Algorithm 3 of the paper).
!>
!> Candidates:  N_cand = floor(N(N-1) F sigma_c(theta) v_max dt / V + U),
!>              F = A_sup G_max / A_bar.
!> Stage 1 (speed screen): accept if |g.e| >= U v_max, e uniform on the sphere;
!>              the mean of |g.e|/v_max over e is g/(2 v_max).
!> Stage 2 (orientation and spin screen): accept with A_perp G(Xi)/(A_sup G_max),
!>              Xi = (|w_i| + |w_j|)(L + d)/g.
!> Encounter:   loss eps = E[eps|z] B/<B>, E[eps|z] = c_t z + c_r (1 - z), exchange z' ~ Q(U; a),
!>              uniform split of E'_rot, deflection cos(chi) ~ exp(a c + b P2(c)),
!>              spins reset uniformly in the plane perpendicular to each axis.
module collisions
  use precision, only: dp, i8, pi
  use random, only: uniform, unit_vector, unit_vector_perpendicular, beta_variate
  use particles, only: gas_t
  use closure_model, only: model_t, stencil_t, dynamic_factor, spin_factor, &
                           sample_partition, sample_cosine
  implicit none
  private
  public :: counters_t, collision_step

  real(dp), parameter :: loss_cap = 1 - 1.0e-9_dp    ! a draw cannot remove all of E

  type :: counters_t
    integer(i8) :: candidates = 0, speed_accepted = 0, encounters = 0
    integer(i8) :: majorant_violations = 0, energy_clamps = 0, theta_outside_steps = 0
    integer(i8) :: loss_caps = 0
    real(dp) :: contacts = 0                 ! sum of <k>(theta) over encounters
    real(dp) :: dissipated = 0               ! sum of eps E over encounters
  end type counters_t

contains

  subroutine collision_step(g, m, s, theta, k_contacts, dt, vmax, vseen, cnt)
    type(gas_t), intent(inout) :: g
    type(model_t), intent(in) :: m
    type(stencil_t), intent(in) :: s
    real(dp), intent(in) :: theta, k_contacts, dt, vmax
    real(dp), intent(inout) :: vseen
    type(counters_t), intent(inout) :: cnt
    real(dp) :: sigma, inflation, mean_candidates, reach
    real(dp) :: gv(3), e(3), ge, speed, ghat(3), xi, area
    integer(i8) :: ncand, icand
    integer :: i, j

    sigma = g%a_bar * dynamic_factor(m, theta)
    inflation = g%a_sup * m%g_max / g%a_bar
    mean_candidates = real(g%n, dp) * real(g%n - 1, dp) * inflation * sigma * vmax * dt / g%volume
    ncand = int(mean_candidates + uniform(), i8)
    reach = g%length + g%diameter
    cnt%candidates = cnt%candidates + ncand

    do icand = 1, ncand
      i = 1 + int(uniform() * g%n)
      j = 1 + int(uniform() * (g%n - 1))
      if (j >= i) j = j + 1

      ! stage 1: flux of relative motion through a random plane (rate ~ g)
      gv = g%c(:, i) - g%c(:, j)
      e = unit_vector()
      ge = abs(dot_product(gv, e))
      vseen = max(vseen, ge)
      if (ge > vmax) cnt%majorant_violations = cnt%majorant_violations + 1
      if (ge < uniform() * vmax) cycle
      cnt%speed_accepted = cnt%speed_accepted + 1

      ! stage 2: projected excluded area and its enhancement by rotation
      speed = norm2(gv)
      if (speed <= 0) cycle
      ghat = gv / speed
      xi = (norm2(g%w(:, i)) + norm2(g%w(:, j))) * reach / speed
      area = projected_area(g, g%u(:, i), g%u(:, j), ghat)
      if (uniform() * g%a_sup * m%g_max > area * spin_factor(m, s, xi)) cycle

      call encounter(g, m, s, i, j, gv, cnt)
      cnt%contacts = cnt%contacts + k_contacts
    end do
  end subroutine collision_step

  !> Projected excluded area of a non-rotating pair seen along ghat (Steiner):
  !> pi d^2 + 2 d L (|u_i x g| + |u_j x g|) + L^2 |g.(u_i x u_j)|.
  pure real(dp) function projected_area(g, ui, uj, gh) result(a)
    type(gas_t), intent(in) :: g
    real(dp), intent(in) :: ui(3), uj(3), gh(3)
    real(dp) :: d, l
    d = g%diameter
    l = g%length
    a = pi * d * d + 2 * d * l * (norm2(cross(ui, gh)) + norm2(cross(uj, gh))) &
      + l * l * abs(dot_product(gh, cross(ui, uj)))
  end function projected_area

  subroutine encounter(g, m, s, i, j, gv, cnt)
    type(gas_t), intent(inout) :: g
    type(model_t), intent(in) :: m
    type(stencil_t), intent(in) :: s
    integer, intent(in) :: i, j
    real(dp), intent(in) :: gv(3)
    type(counters_t), intent(inout) :: cnt
    real(dp) :: etr, erot, e, z, eps, zp, ep, etr_out, erot_out, split
    real(dp) :: ghat(3), e1(3), e2(3), trial(3), cosx, sinx, phi, gout(3), vcm(3)

    ! incoming pair energy and partition (centre-of-mass motion is untouched)
    etr = 0.25_dp * g%mass * dot_product(gv, gv)
    erot = g%erot(i) + g%erot(j)
    e = etr + erot
    if (e <= 0) return
    z = etr / e

    ! loss: universal Beta shape on the measured conditional mean loss
    ! E[eps|z] = c_t z + c_r (1 - z) (mixture of the node means); older model
    ! files without the loss table use the unconditional node mean
    if (m%loss_memory) then
      eps = max(s%loss_ct * z + s%loss_cr * (1 - z), 0.0_dp) &
            * beta_variate(m%beta_a, m%beta_b) / m%mean_beta
    else
      eps = s%mean_loss * beta_variate(m%beta_a, m%beta_b) / m%mean_beta
    end if
    if (eps >= loss_cap) then
      eps = loss_cap
      cnt%loss_caps = cnt%loss_caps + 1
    end if

    ! outgoing partition of the surviving energy
    zp = sample_partition(m, s, z, eps, cnt%energy_clamps)
    ep = (1 - eps) * e
    etr_out = zp * ep
    erot_out = (1 - zp) * ep
    if (etr_out <= 0 .or. erot_out <= 0) error stop 'collisions: non-positive modal energy'

    ! rotational energy split between the two rods (2 + 2 DOF: uniform)
    split = uniform()
    g%erot(i) = split * erot_out
    g%erot(j) = (1 - split) * erot_out

    ! deflection of the relative velocity, uniform azimuth about ghat
    cosx = sample_cosine(m, s, z, zp)
    sinx = sqrt(max(0.0_dp, 1 - cosx * cosx))
    ghat = gv / norm2(gv)
    trial = [1.0_dp, 0.0_dp, 0.0_dp]
    if (abs(ghat(1)) >= 0.9_dp) trial = [0.0_dp, 1.0_dp, 0.0_dp]
    e1 = trial - dot_product(trial, ghat) * ghat
    e1 = e1 / norm2(e1)
    e2 = cross(ghat, e1)
    phi = 2 * pi * uniform()
    gout = 2 * sqrt(etr_out / g%mass) * (cosx * ghat + sinx * (cos(phi) * e1 + sin(phi) * e2))
    vcm = 0.5_dp * (g%c(:, i) + g%c(:, j))
    g%c(:, i) = vcm + 0.5_dp * gout
    g%c(:, j) = vcm - 0.5_dp * gout

    ! spins: equilibrium direction in the plane perpendicular to each axis
    g%w(:, i) = sqrt(2 * g%erot(i) / g%inertia) * unit_vector_perpendicular(g%u(:, i))
    g%w(:, j) = sqrt(2 * g%erot(j) / g%inertia) * unit_vector_perpendicular(g%u(:, j))

    cnt%encounters = cnt%encounters + 1
    cnt%dissipated = cnt%dissipated + eps * e
  end subroutine encounter

  pure function cross(a, b) result(c)
    real(dp), intent(in) :: a(3), b(3)
    real(dp) :: c(3)
    c = [a(2) * b(3) - a(3) * b(2), a(3) * b(1) - a(1) * b(3), a(1) * b(2) - a(2) * b(1)]
  end function cross

end module collisions
