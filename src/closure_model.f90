!> Collision closures of one (alpha, AR) plane, as functions of theta = T_tr/T_rot.
!>
!> The model file (tools/export_model.py) holds the closures at a few theta
!> nodes.  Between nodes they are combined with linear weights in theta:
!>   * scalars (mean loss, spin factor G) are weighted sums,
!>   * the exchange law is the weighted sum of node quantiles at a shared U
!>     (displacement interpolation of the node laws),
!>   * the angular law is the mixture of the node laws.
!> The encounter cross-section is sigma_c = A_bar D(theta), with D linear in
!> log(theta).  Outside the node range every closure is held at its end value
!> and the excursion is counted.
module closure_model
  use precision, only: dp, i8, pi
  use npz_io, only: npz_t, npz_read, npz_real1, npz_real2, npz_int1, npz_scalar, npz_has
  use random, only: uniform
  implicit none
  private
  public :: model_t, stencil_t, load_model, make_stencil, dynamic_factor, &
            contacts_per_encounter, spin_factor, sample_partition, sample_cosine

  integer, parameter :: max_xi = 64

  type :: model_t
    real(dp) :: alpha = 0, aspect_ratio = 0
    real(dp) :: beta_a = 1.21_dp, beta_b = 3.67_dp      ! loss shape B ~ Beta(a, b)
    real(dp) :: mean_beta = 0                           ! <B> = a/(a+b)
    integer :: nodes = 0
    real(dp), allocatable :: theta(:)                   ! (K) node temperature ratios
    real(dp), allocatable :: mean_loss(:)               ! (K) eps_bar per encounter
    logical :: loss_memory = .false.                   ! loss conditioned on z
    real(dp), allocatable :: loss_ct(:), loss_cr(:)    ! (K) E[eps|z] = c_t z + c_r (1 - z)
    logical :: loss_bounded = .false.                  ! eps ~ Beta(kappa mu, kappa (1 - mu))
    real(dp), allocatable :: loss_kappa(:)             ! (K) Beta concentration (0 at alpha = 1)
    character(len=16) :: contact_model = 'R1-legacy'   ! C1 = contact-point damping, Hertz-Tsuji CTC
    real(dp), allocatable :: xi_grid(:)                 ! (NX)
    real(dp), allocatable :: xi_enhancement(:, :)       ! (NX, K) spin factor G
    real(dp) :: g_max = 0                               ! max of G over the plane
    real(dp), allocatable :: lambda1(:), lambda4(:)     ! (K) dissipative tilt
    integer, allocatable :: memory_form(:)              ! (K) 0 linear, 1 logit cubic
    real(dp), allocatable :: memory_coef(:, :)          ! (3, K)
    real(dp), allocatable :: memory_center(:), memory_scale(:)
    integer :: nu = 0                                   ! quantile grid size
    integer, allocatable :: row_first(:), row_last(:)   ! (K) rows of each node
    real(dp), allocatable :: a_grid(:)                  ! (rows)
    real(dp), allocatable :: quantile(:, :)             ! (NU, rows): Q(u_j; a_r)
    real(dp), allocatable :: angular(:, :)              ! (10, K)
    real(dp), allocatable :: cs_log_theta(:), cs_factor(:)
    real(dp), allocatable :: k_theta(:), k_mean(:)
  end type model_t

  !> The closure state at one theta: at most two nodes with their weights.
  type :: stencil_t
    integer :: n = 1
    integer :: node(2) = 1
    real(dp) :: w(2) = [1.0_dp, 0.0_dp]
    real(dp) :: mean_loss = 0                           ! eps_bar(theta)
    real(dp) :: loss_ct = 0, loss_cr = 0                ! loss rates at theta
    real(dp) :: loss_kappa = 0                          ! Beta concentration at theta
    integer :: nx = 0
    real(dp) :: g(max_xi) = 1.0_dp                      ! G(Xi; theta) on xi_grid
    logical :: outside = .false.                        ! theta beyond the node range
  end type stencil_t

contains

  subroutine load_model(path, m)
    character(len=*), intent(in) :: path
    type(model_t), intent(out) :: m
    type(npz_t) :: z
    real(dp), allocatable :: tmp(:), cs_theta(:)
    integer, allocatable :: offsets(:)
    integer :: k

    call npz_read(path, z)
    m%alpha = npz_scalar(z, 'alpha')
    m%aspect_ratio = npz_scalar(z, 'aspect_ratio')
    m%beta_a = npz_scalar(z, 'loss_beta_a')
    m%beta_b = npz_scalar(z, 'loss_beta_b')
    m%mean_beta = m%beta_a / (m%beta_a + m%beta_b)

    call npz_real1(z, 'theta_nodes', m%theta)
    m%nodes = size(m%theta)
    call npz_real1(z, 'mean_loss', m%mean_loss)
    m%loss_memory = npz_has(z, 'loss_ct')
    if (m%loss_memory) then
      call npz_real1(z, 'loss_ct', m%loss_ct)
      call npz_real1(z, 'loss_cr', m%loss_cr)
    end if
    m%loss_bounded = npz_has(z, 'loss_kappa')
    if (m%loss_bounded) then
      if (.not. m%loss_memory) error stop 'closure_model: loss_kappa without loss_ct/loss_cr'
      call npz_real1(z, 'loss_kappa', m%loss_kappa)
    end if
    if (npz_has(z, 'contact_model_c1')) then
      if (npz_scalar(z, 'contact_model_c1') > 0.5_dp) m%contact_model = 'C1'
    end if
    call npz_real1(z, 'xi_grid', m%xi_grid)
    call npz_real2(z, 'xi_enhancement', m%xi_enhancement)
    if (size(m%xi_grid) > max_xi) error stop 'closure_model: xi grid larger than max_xi'
    m%g_max = maxval(m%xi_enhancement)
    call npz_real1(z, 'lambda1', m%lambda1)
    call npz_real1(z, 'lambda4', m%lambda4)
    call npz_int1(z, 'memory_form', m%memory_form)
    call npz_real2(z, 'memory_coef', m%memory_coef)
    call npz_real1(z, 'memory_center', m%memory_center)
    call npz_real1(z, 'memory_scale', m%memory_scale)

    call npz_real1(z, 'quantile_probability', tmp)
    m%nu = size(tmp)
    call npz_int1(z, 'energy_a_offsets', offsets)
    allocate(m%row_first(m%nodes), m%row_last(m%nodes))
    do k = 1, m%nodes
      m%row_first(k) = offsets(k) + 1
      m%row_last(k) = offsets(k + 1)
    end do
    call npz_real1(z, 'energy_a_grid', m%a_grid)
    call npz_real2(z, 'energy_quantiles', m%quantile)
    call npz_real2(z, 'angular', m%angular)

    call npz_real1(z, 'cs_theta', cs_theta)
    m%cs_log_theta = log(cs_theta)
    call npz_real1(z, 'cs_dynamic_factor', m%cs_factor)
    call npz_real1(z, 'k_theta', m%k_theta)
    call npz_real1(z, 'k_mean', m%k_mean)
  end subroutine load_model

  !> Node weights and scalar closures at temperature ratio theta.
  subroutine make_stencil(m, theta, s)
    type(model_t), intent(in) :: m
    real(dp), intent(in) :: theta
    type(stencil_t), intent(out) :: s
    integer :: k
    real(dp) :: f

    s%outside = theta < m%theta(1) .or. theta > m%theta(m%nodes)
    if (theta <= m%theta(1)) then
      s%n = 1; s%node(1) = 1; s%w = [1.0_dp, 0.0_dp]
    else if (theta >= m%theta(m%nodes)) then
      s%n = 1; s%node(1) = m%nodes; s%w = [1.0_dp, 0.0_dp]
    else
      k = 1
      do while (m%theta(k + 1) < theta)
        k = k + 1
      end do
      f = (theta - m%theta(k)) / (m%theta(k + 1) - m%theta(k))
      s%n = 2; s%node = [k, k + 1]; s%w = [1.0_dp - f, f]
    end if
    s%nx = size(m%xi_grid)
    s%mean_loss = 0
    s%loss_ct = 0
    s%loss_cr = 0
    s%loss_kappa = 0
    s%g(1:s%nx) = 0
    do k = 1, s%n
      s%mean_loss = s%mean_loss + s%w(k) * m%mean_loss(s%node(k))
      if (m%loss_memory) then
        s%loss_ct = s%loss_ct + s%w(k) * m%loss_ct(s%node(k))
        s%loss_cr = s%loss_cr + s%w(k) * m%loss_cr(s%node(k))
      end if
      if (m%loss_bounded) s%loss_kappa = s%loss_kappa + s%w(k) * m%loss_kappa(s%node(k))
      s%g(1:s%nx) = s%g(1:s%nx) + s%w(k) * m%xi_enhancement(:, s%node(k))
    end do
  end subroutine make_stencil

  !> Dynamic factor D(theta) of the encounter cross-section, linear in log(theta).
  real(dp) function dynamic_factor(m, theta)
    type(model_t), intent(in) :: m
    real(dp), intent(in) :: theta
    dynamic_factor = interp(log(theta), m%cs_log_theta, m%cs_factor)
  end function dynamic_factor

  !> Mean number of contacts in one encounter (diagnostic: DEM counts contacts).
  real(dp) function contacts_per_encounter(m, theta)
    type(model_t), intent(in) :: m
    real(dp), intent(in) :: theta
    contacts_per_encounter = interp(theta, m%k_theta, m%k_mean)
  end function contacts_per_encounter

  !> Spin factor G(Xi) of the encounter area at the current stencil.
  real(dp) function spin_factor(m, s, xi)
    type(model_t), intent(in) :: m
    type(stencil_t), intent(in) :: s
    real(dp), intent(in) :: xi
    spin_factor = interp(xi, m%xi_grid, s%g(1:s%nx))
  end function spin_factor

  !> Outgoing translational fraction z' given the incoming fraction z and the
  !> loss eps:  z' = sum_k w_k Q_k(U; a_k),  a_k = lambda1 + memory(z) + lambda4 eps_k,
  !> with eps_k = eps * eps_bar_k / eps_bar(theta) on the loss scale of node k.
  real(dp) function sample_partition(m, s, z, eps, clamps) result(zp)
    type(model_t), intent(in) :: m
    type(stencil_t), intent(in) :: s
    real(dp), intent(in) :: z, eps
    integer(i8), intent(inout) :: clamps
    real(dp) :: u, x, fu, a, memory, cov, blend, q_lo, q_hi
    integer :: j, k, nd, lo, hi, r

    u = uniform()
    x = u * real(m%nu - 1, dp)                       ! uniform probability grid
    j = min(int(x) + 1, m%nu - 1)
    fu = x - real(j - 1, dp)
    zp = 0
    do k = 1, s%n
      nd = s%node(k)
      memory = memory_term(m, nd, z)
      cov = eps
      if (s%mean_loss > 0 .and. m%mean_loss(nd) > 0) cov = eps * m%mean_loss(nd) / s%mean_loss
      a = m%lambda1(nd) + memory + m%lambda4(nd) * cov
      lo = m%row_first(nd)
      hi = m%row_last(nd)
      if (a < m%a_grid(lo) .or. a > m%a_grid(hi)) then
        clamps = clamps + 1
        a = min(max(a, m%a_grid(lo)), m%a_grid(hi))
      end if
      r = bracket(m%a_grid(lo:hi), a) + lo - 1
      blend = (a - m%a_grid(r)) / (m%a_grid(r + 1) - m%a_grid(r))
      q_lo = (1 - fu) * m%quantile(j, r) + fu * m%quantile(j + 1, r)
      q_hi = (1 - fu) * m%quantile(j, r + 1) + fu * m%quantile(j + 1, r + 1)
      zp = zp + s%w(k) * ((1 - blend) * q_lo + blend * q_hi)
    end do
  end function sample_partition

  !> Memory of the incoming partition in the natural parameter of node nd.
  real(dp) function memory_term(m, nd, z)
    type(model_t), intent(in) :: m
    integer, intent(in) :: nd
    real(dp), intent(in) :: z
    real(dp) :: zc, x
    if (m%memory_form(nd) == 1) then              ! bounded logit cubic (low theta, near spheres)
      zc = min(max(z, 1.0e-12_dp), 1.0_dp - 1.0e-12_dp)
      x = tanh(((log(zc / (1 - zc)) - m%memory_center(nd)) / max(m%memory_scale(nd), 1.0e-8_dp)) / 2)
      memory_term = m%memory_coef(1, nd) * x + m%memory_coef(2, nd) * x**2 + m%memory_coef(3, nd) * x**3
    else                                          ! Sinkhorn bridge: lambda3 z
      memory_term = m%memory_coef(1, nd) * z
    end if
  end function memory_term

  !> Scattering cosine c = ghat.ghat' from the node mixture of
  !> p(c | z, z') ~ exp[a c + b P2(c)],  a = eta10 + eta11 z' + eta12 z,
  !> b = eta20 + eta21 z' + eta22 z, each clamped to the range of its node.
  real(dp) function sample_cosine(m, s, z, zp) result(c)
    type(model_t), intent(in) :: m
    type(stencil_t), intent(in) :: s
    real(dp), intent(in) :: z, zp
    real(dp) :: a, b, fmax, vertex, r
    integer :: nd
    nd = s%node(1)
    if (s%n == 2) then
      if (uniform() >= s%w(1)) nd = s%node(2)
    end if
    associate (p => m%angular(:, nd))
      a = min(max(p(1) + p(2) * zp + p(3) * z, p(7)), p(8))
      b = min(max(p(4) + p(5) * zp + p(6) * z, p(9)), p(10))
    end associate
    fmax = max(a + b, -a + b)                      ! f(c) = a c + b P2(c) at c = +-1
    if (b < 0) then
      vertex = -a / (3 * b)
      if (abs(vertex) < 1) fmax = max(fmax, a * vertex + b * 0.5_dp * (3 * vertex**2 - 1))
    end if
    do
      c = 2 * uniform() - 1
      r = a * c + b * 0.5_dp * (3 * c * c - 1) - fmax
      if (uniform() <= exp(r)) exit
    end do
  end function sample_cosine

  !> Linear interpolation, held constant outside [x(1), x(n)].
  pure real(dp) function interp(x0, x, y)
    real(dp), intent(in) :: x0, x(:), y(:)
    integer :: i, n
    n = size(x)
    if (x0 <= x(1)) then
      interp = y(1)
    else if (x0 >= x(n)) then
      interp = y(n)
    else
      i = bracket(x, x0)
      interp = y(i) + (y(i + 1) - y(i)) * (x0 - x(i)) / (x(i + 1) - x(i))
    end if
  end function interp

  !> Index i with x(i) <= x0 < x(i+1), restricted to 1..n-1 (binary search).
  pure integer function bracket(x, x0) result(lo)
    real(dp), intent(in) :: x(:), x0
    integer :: hi, mid
    lo = 1
    hi = size(x)
    do while (hi - lo > 1)
      mid = (lo + hi) / 2
      if (x(mid) <= x0) then
        lo = mid
      else
        hi = mid
      end if
    end do
  end function bracket

end module closure_model
