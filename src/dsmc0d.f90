!> Zero-dimensional DSMC of a dilute gas of inelastic smooth spherocylinders.
!>
!>   dsmc0d input.nml
!>
!> The collision model (one restitution coefficient and aspect ratio) is read
!> from a single .npz file made by tools/export_model.py.  The program follows
!> the homogeneous cooling state (flow = 'hcs') or the uniform shear flow
!> (flow = 'usf', shear rate gdot) from T_tr = T_tr0, T_rot = T_rot0 and writes
!>   <prefix>.ts     time series (one row every output_dtau contacts/particle)
!>   <prefix>.stress USF only: kinetic stress P^k_ij/(n T_tr), T* = T_tr/(m gdot^2 d^2)
!>                   and the energy ledger (thermal energy, shear work, collisional change)
!>   <prefix>.speed  reduced speed distribution     (if sampling is enabled)
!>   <prefix>.rot    reduced rotational-energy distribution
!>   <prefix>.log    parameters, counters and wall-clock times
!>
!> Clocks: tau_enc counts encounters per particle (one DSMC event is one
!> binary encounter); tau_c counts contacts per particle, tau_c = sum <k>(theta)
!> over encounters, which is what a DEM simulation counts.
program dsmc0d
  use precision, only: dp, i8, pi
  use random, only: rng_seed
  use closure_model, only: model_t, stencil_t, load_model, make_stencil, &
                           dynamic_factor, contacts_per_encounter, ablation_mode, ablation_code
  use particles, only: gas_t, setup_gas, initialize, temperatures, rotate_axes, &
                       non_gaussianity, shear_step, kinetic_stress, thermal_energy
  use collisions, only: counters_t, collision_step
  use distributions, only: histogram_t, init_histograms, sample_histograms, write_histograms
  implicit none

  ! ---- input (namelist /run/)
  character(len=512) :: model_file = 'model.npz', output_prefix = 'run'
  integer :: n_particles = 100000, seed = 1
  real(dp) :: volume_fraction = 0.01_dp, density = 0.763944_dp, diameter = 1.0_dp
  real(dp) :: T_tr0 = 1.0_dp, T_rot0 = 1.0_dp
  logical :: exact_start = .true.
  real(dp) :: tau_end = 25.0_dp              ! contacts per particle
  real(dp) :: t_end = huge(1.0_dp)
  real(dp) :: collision_fraction = 0.02_dp   ! encounters per particle per time step
  real(dp) :: output_dtau = 0.05_dp          ! contacts per particle between rows
  real(dp) :: sample_tau_start = huge(1.0_dp), sample_dtau = 0.5_dp
  real(dp) :: rescale_below = 0.0_dp         ! similarity rescaling when T < rescale_below T(0); 0 = off
  character(len=8) :: flow = 'hcs'           ! 'hcs' (cooling) or 'usf' (uniform shear flow)
  real(dp) :: shear_rate = 0.0_dp            ! USF: gdot
  real(dp) :: stress_tau_start = huge(1.0_dp) ! USF: steady averages over tau_c >= this
  real(dp) :: max_shear_step = 0.01_dp       ! USF: largest strain gdot dt per time step
  character(len=32) :: ablation = 'none'     ! none, no_spin_screen, loss_unconditioned,
                                             ! lambda4_off, isotropic_scattering, memoryless_exchange
  namelist /run/ model_file, output_prefix, n_particles, seed, volume_fraction, density, &
                 diameter, T_tr0, T_rot0, exact_start, tau_end, t_end, collision_fraction, &
                 output_dtau, sample_tau_start, sample_dtau, rescale_below, &
                 flow, shear_rate, stress_tau_start, max_shear_step, ablation

  type(model_t) :: model
  type(stencil_t) :: st
  type(gas_t) :: gas
  type(counters_t) :: cnt
  type(histogram_t) :: hspeed, hrot
  character(len=512) :: input_file
  integer :: u, ts, ios
  integer(i8) :: clock0, clock1, rate, steps, samples
  real(dp) :: t, ttr, trot, theta, tau_c, tau_enc, dt, nu_enc, vmax, vseen, vprev
  real(dp) :: kc, a2_tr, a2_rot, next_out, next_sample
  real(dp) :: wall_tau_c1, wall_tau_enc1, wall_total, theta_min, theta_max
  real(dp) :: t_ref, factor, vmean(3)
  integer :: k
  integer(i8) :: rescalings
  ! USF
  logical :: usf, in_window
  integer :: sf, nrows
  real(dp) :: pk(6), tstar, eth0, eth, e1, w_shear, w_coll
  real(dp) :: win_t0, win_tau0, win_wsh0, win_wcoll0, win_eth0, smean(10), sse(10)
  real(dp), allocatable :: steady(:, :)

  ! ---------------------------------------------------------------- input
  if (command_argument_count() < 1) then
    write(*, '(a)') 'usage: dsmc0d input.nml'
    stop 1
  end if
  call get_command_argument(1, input_file)
  open(newunit=u, file=trim(input_file), status='old', action='read', iostat=ios)
  if (ios /= 0) error stop 'cannot open the input file'
  read(u, nml=run, iostat=ios)
  if (ios /= 0) error stop 'cannot read namelist /run/'
  close(u)
  usf = trim(flow) == 'usf'
  if (.not. usf .and. trim(flow) /= 'hcs') error stop 'flow must be hcs or usf'
  if (usf .and. shear_rate <= 0) error stop 'flow = usf needs shear_rate > 0'
  if (usf .and. rescale_below > 0) error stop 'similarity rescaling is for the cooling state only'
  ablation_mode = ablation_code(ablation)
  if (ablation_mode < 0) error stop 'unknown ablation'

  call load_model(trim(model_file), model)
  call rng_seed(seed)
  call setup_gas(gas, n_particles, model%aspect_ratio, diameter, density, volume_fraction)
  call initialize(gas, T_tr0, T_rot0, exact_start)
  call init_histograms(hspeed, hrot)

  open(newunit=ts, file=trim(output_prefix)//'.ts', status='replace')
  write(ts, '(a)') '# t  tau_enc  tau_c  T_tr  T_rot  theta  T  a2_tr  a2_rot  dt'
  if (usf) then
    open(newunit=sf, file=trim(output_prefix)//'.stress', status='replace')
    write(sf, '(a)') '# t  tau_c  T_tr  T_rot  theta  Tstar  Pk_xx  Pk_yy  Pk_zz  Pk_xy  Pk_xz  Pk_yz' // &
                     '  E_th  W_shear  W_coll  resid   (P^k_ij reduced by n T_tr; W cumulative)'
  end if

  ! ---------------------------------------------------------------- time march
  t = 0; dt = 0; steps = 0; samples = 0
  next_out = 0; next_sample = sample_tau_start
  vprev = 0
  wall_tau_c1 = -1; wall_tau_enc1 = -1
  rescalings = 0
  call temperatures(gas, ttr, trot)
  t_ref = (3 * ttr + 2 * trot) / 5
  theta_min = huge(1.0_dp); theta_max = 0
  eth0 = thermal_energy(gas); w_shear = 0; w_coll = 0
  in_window = .false.; nrows = 0
  win_t0 = 0; win_tau0 = 0; win_wsh0 = 0; win_wcoll0 = 0; win_eth0 = 0
  allocate(steady(10, 1024))
  call system_clock(clock0, rate)
  do
    call temperatures(gas, ttr, trot)
    ! Similarity rescaling: in the hard-contact limit an encounter is unchanged by a
    ! uniform rescaling of all velocities and spins, so the cooling state can be
    ! followed to long times at a fixed temperature scale.  The time column is then
    ! in rescaled units; the collision count and the reduced quantities are not affected.
    if (rescale_below > 0 .and. (3 * ttr + 2 * trot) / 5 < rescale_below * t_ref) then
      ! The mean velocity is zero up to round-off and collisions conserve it, so it
      ! does not cool; it is removed here, or repeated rescalings would amplify it
      ! until it swamped the thermal velocities.
      factor = sqrt(t_ref / ((3 * ttr + 2 * trot) / 5))
      vmean = sum(gas%c, dim=2) / gas%n
      do k = 1, gas%n
        gas%c(:, k) = factor * (gas%c(:, k) - vmean)
      end do
      gas%w = factor * gas%w
      gas%erot = factor**2 * gas%erot
      vprev = factor * vprev
      rescalings = rescalings + 1
      call temperatures(gas, ttr, trot)
    end if
    theta = ttr / trot
    tau_enc = 2 * real(cnt%encounters, dp) / gas%n
    tau_c = 2 * cnt%contacts / gas%n

    if (tau_c >= next_out - 1.0e-12_dp) then
      call non_gaussianity(gas, ttr, trot, a2_tr, a2_rot)
      write(ts, '(10es17.9)') t, tau_enc, tau_c, ttr, trot, theta, (3 * ttr + 2 * trot) / 5, &
                               a2_tr, a2_rot, dt
      if (usf) then
        call kinetic_stress(gas, ttr, pk)
        tstar = ttr / (gas%mass * (shear_rate * gas%diameter)**2)
        eth = thermal_energy(gas)
        write(sf, '(16es17.9)') t, tau_c, ttr, trot, theta, tstar, pk, eth, w_shear, w_coll, &
                                eth - eth0 - w_shear - w_coll
        if (tau_c >= stress_tau_start .and. tau_c <= tau_end) then
          if (.not. in_window) then
            in_window = .true.
            win_t0 = t; win_tau0 = tau_c; win_wsh0 = w_shear; win_wcoll0 = w_coll; win_eth0 = eth
          end if
          call append_row(steady, nrows, [tstar, theta, pk(1), pk(2), pk(3), pk(4), &
                                          pk(1) - pk(2), pk(2) - pk(3), a2_tr, a2_rot])
        end if
      end if
      next_out = next_out + output_dtau
    end if
    if (tau_c >= next_sample .and. tau_c <= tau_end) then
      call sample_histograms(gas, ttr, trot, hspeed, hrot)
      samples = samples + 1
      next_sample = next_sample + sample_dtau
    end if
    if (tau_c >= tau_end .or. t >= t_end) exit

    ! closures and time step at the current temperature ratio
    theta_min = min(theta_min, theta)
    theta_max = max(theta_max, theta)
    call make_stencil(model, theta, st)
    if (st%outside) cnt%theta_outside_steps = cnt%theta_outside_steps + 1
    kc = contacts_per_encounter(model, theta)
    nu_enc = gas%number_density * gas%a_bar * dynamic_factor(model, theta) &
           * 4 * sqrt(ttr / (pi * gas%mass))
    dt = collision_fraction / nu_enc
    if (usf) dt = min(dt, max_shear_step / shear_rate)
    vmax = max(5 * sqrt(2.0_dp) * sqrt((ttr + 2 * trot / 3) / gas%mass), vprev)

    ! Strang splitting: half free rotation, collisions, half free rotation
    ! (USF: the exact shear map is split the same way around the collisions)
    call rotate_axes(gas, 0.5_dp * dt)
    if (usf) call shear_step(gas, 0.5_dp * shear_rate * dt, w_shear)
    vseen = 0
    if (usf) e1 = thermal_energy(gas)
    call collision_step(gas, model, st, theta, kc, dt, vmax, vseen, cnt)
    if (usf) w_coll = w_coll + thermal_energy(gas) - e1
    vprev = vseen
    if (usf) call shear_step(gas, 0.5_dp * shear_rate * dt, w_shear)
    call rotate_axes(gas, 0.5_dp * dt)
    t = t + dt
    steps = steps + 1

    call system_clock(clock1)
    if (wall_tau_c1 < 0 .and. 2 * cnt%contacts / gas%n >= 1) &
      wall_tau_c1 = real(clock1 - clock0, dp) / rate
    if (wall_tau_enc1 < 0 .and. 2 * real(cnt%encounters, dp) / gas%n >= 1) &
      wall_tau_enc1 = real(clock1 - clock0, dp) / rate
  end do
  call system_clock(clock1)
  wall_total = real(clock1 - clock0, dp) / rate
  close(ts)
  if (usf) close(sf)
  if (samples > 0) call write_histograms(trim(output_prefix), hspeed, hrot)

  ! ---------------------------------------------------------------- log
  open(newunit=u, file=trim(output_prefix)//'.log', status='replace')
  if (usf) then
    write(u, '(a)') 'dsmc0d: uniform shear flow of smooth spherocylinders'
  else
    write(u, '(a)') 'dsmc0d: homogeneous cooling of smooth spherocylinders'
  end if
  write(u, nml=run)
  write(u, '(a,f8.4,a,f6.3)') 'model: alpha = ', model%alpha, '  AR = ', model%aspect_ratio
  write(u, '(a,*(es12.4))') 'theta nodes: ', model%theta
  write(u, '(a,4es15.7)') 'mass, inertia, number density, volume: ', gas%mass, gas%inertia, &
                          gas%number_density, gas%volume
  write(u, '(a,2es15.7)') 'A_bar, A_sup: ', gas%a_bar, gas%a_sup
  write(u, '(a,es15.7)') 'G_max: ', model%g_max
  write(u, '(a,i0)') 'steps: ', steps
  write(u, '(a,i0)') 'candidates: ', cnt%candidates
  write(u, '(a,i0)') 'speed-screen accepted: ', cnt%speed_accepted
  write(u, '(a,i0)') 'encounters: ', cnt%encounters
  write(u, '(a,f12.5)') 'tau_enc (encounters/particle): ', 2 * real(cnt%encounters, dp) / gas%n
  write(u, '(a,f12.5)') 'tau_c (contacts/particle): ', 2 * cnt%contacts / gas%n
  write(u, '(a,es15.7)') 'final time: ', t
  write(u, '(a,2f10.5)') 'theta range: ', theta_min, theta_max
  write(u, '(a,i0)') 'majorant violations: ', cnt%majorant_violations
  write(u, '(a,i0)') 'energy-axis clamps: ', cnt%energy_clamps
  write(u, '(a,i0)') 'loss caps: ', cnt%loss_caps
  write(u, '(a,l1)') 'loss conditioned on z: ', model%loss_memory
  write(u, '(a,l1)') 'bounded loss law (Beta with mean E[eps|z]): ', model%loss_bounded
  write(u, '(2a)') 'contact model: ', trim(model%contact_model)
  write(u, '(2a)') 'ablation: ', trim(ablation)
  write(u, '(a,i0)') 'similarity rescalings: ', rescalings
  write(u, '(a,i0)') 'steps with theta outside the nodes: ', cnt%theta_outside_steps
  write(u, '(a,i0)') 'distribution samples: ', samples
  write(u, '(a,f14.4)') 'wall seconds to tau_c = 1: ', wall_tau_c1
  write(u, '(a,f14.4)') 'wall seconds to tau_enc = 1: ', wall_tau_enc1
  if (usf) then
    write(u, '(a,es15.7)') 'shear rate: ', shear_rate
    write(u, '(a,i0,a,2f12.4)') 'steady window rows: ', nrows, '  tau_c from, to: ', win_tau0, tau_c
    if (nrows >= 20) then
      call block_stats(steady(:, 1:nrows), smean, sse)
      write(u, '(a)') 'steady quantities: Tstar theta Pk_xx Pk_yy Pk_zz Pk_xy N1 N2 a2_tr a2_rot' // &
                      ' (P^k reduced by n T_tr; N1 = Pk_xx - Pk_yy, N2 = Pk_yy - Pk_zz)'
      write(u, '(a,10es15.7)') 'steady mean: ', smean
      write(u, '(a,10es15.7)') 'steady se (10 blocks): ', sse
      write(u, '(a,es15.7)') 'contacts per particle per unit strain: ', &
                             (tau_c - win_tau0) / (shear_rate * (t - win_t0))
      write(u, '(a,4es15.7)') 'steady window energy: shear work, collisional change, thermal change, ' // &
                              '(shear + collisional)/shear: ', w_shear - win_wsh0, w_coll - win_wcoll0, &
                              thermal_energy(gas) - win_eth0, &
                              (w_shear - win_wsh0 + w_coll - win_wcoll0) / max(abs(w_shear - win_wsh0), tiny(1.0_dp))
    end if
    write(u, '(a,es15.7)') 'energy bookkeeping residual / |shear work|: ', &
                           (thermal_energy(gas) - eth0 - w_shear - w_coll) / max(abs(w_shear), tiny(1.0_dp))
  end if
  write(u, '(a,f14.4)') 'wall seconds total (time march): ', wall_total
  close(u)

contains

  !> Append one row to a growing column-major table.
  subroutine append_row(a, n, row)
    real(dp), allocatable, intent(inout) :: a(:, :)
    integer, intent(inout) :: n
    real(dp), intent(in) :: row(:)
    real(dp), allocatable :: tmp(:, :)
    if (n == size(a, 2)) then
      allocate(tmp(size(a, 1), 2 * size(a, 2)))
      tmp(:, 1:n) = a(:, 1:n)
      call move_alloc(tmp, a)
    end if
    n = n + 1
    a(:, n) = row
  end subroutine append_row

  !> Means of the rows and standard errors from 10 consecutive blocks.
  subroutine block_stats(x, mean, se)
    real(dp), intent(in) :: x(:, :)
    real(dp), intent(out) :: mean(:), se(:)
    integer, parameter :: nb = 10
    integer :: b, bl
    real(dp) :: bm(size(x, 1), nb), bbar(size(x, 1))
    bl = size(x, 2) / nb
    do b = 1, nb
      bm(:, b) = sum(x(:, (b - 1) * bl + 1:b * bl), dim=2) / bl
    end do
    bbar = sum(bm, dim=2) / nb
    mean = sum(x, dim=2) / size(x, 2)
    se = sqrt(sum((bm - spread(bbar, 2, nb))**2, dim=2) / (nb * (nb - 1)))
  end subroutine block_stats

end program dsmc0d
