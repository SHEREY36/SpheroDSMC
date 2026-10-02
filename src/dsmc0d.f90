!> Zero-dimensional DSMC of a dilute gas of inelastic smooth spherocylinders.
!>
!>   dsmc0d input.nml
!>
!> The collision model (one restitution coefficient and aspect ratio) is read
!> from a single .npz file made by tools/export_model.py.  The program follows
!> the homogeneous cooling state from T_tr = T_tr0, T_rot = T_rot0 and writes
!>   <prefix>.ts     time series (one row every output_dtau contacts/particle)
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
                           dynamic_factor, contacts_per_encounter
  use particles, only: gas_t, setup_gas, initialize, temperatures, rotate_axes, &
                       non_gaussianity
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
  namelist /run/ model_file, output_prefix, n_particles, seed, volume_fraction, density, &
                 diameter, T_tr0, T_rot0, exact_start, tau_end, t_end, collision_fraction, &
                 output_dtau, sample_tau_start, sample_dtau

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

  call load_model(trim(model_file), model)
  call rng_seed(seed)
  call setup_gas(gas, n_particles, model%aspect_ratio, diameter, density, volume_fraction)
  call initialize(gas, T_tr0, T_rot0, exact_start)
  call init_histograms(hspeed, hrot)

  open(newunit=ts, file=trim(output_prefix)//'.ts', status='replace')
  write(ts, '(a)') '# t  tau_enc  tau_c  T_tr  T_rot  theta  T  a2_tr  a2_rot  dt'

  ! ---------------------------------------------------------------- time march
  t = 0; dt = 0; steps = 0; samples = 0
  next_out = 0; next_sample = sample_tau_start
  vprev = 0
  wall_tau_c1 = -1; wall_tau_enc1 = -1
  theta_min = huge(1.0_dp); theta_max = 0
  call system_clock(clock0, rate)
  do
    call temperatures(gas, ttr, trot)
    theta = ttr / trot
    tau_enc = 2 * real(cnt%encounters, dp) / gas%n
    tau_c = 2 * cnt%contacts / gas%n

    if (tau_c >= next_out - 1.0e-12_dp) then
      call non_gaussianity(gas, ttr, trot, a2_tr, a2_rot)
      write(ts, '(10es17.9)') t, tau_enc, tau_c, ttr, trot, theta, (3 * ttr + 2 * trot) / 5, &
                               a2_tr, a2_rot, dt
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
    vmax = max(5 * sqrt(2.0_dp) * sqrt((ttr + 2 * trot / 3) / gas%mass), vprev)

    ! Strang splitting: half free rotation, collisions, half free rotation
    call rotate_axes(gas, 0.5_dp * dt)
    vseen = 0
    call collision_step(gas, model, st, theta, kc, dt, vmax, vseen, cnt)
    vprev = vseen
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
  if (samples > 0) call write_histograms(trim(output_prefix), hspeed, hrot)

  ! ---------------------------------------------------------------- log
  open(newunit=u, file=trim(output_prefix)//'.log', status='replace')
  write(u, '(a)') 'dsmc0d: homogeneous cooling of smooth spherocylinders'
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
  write(u, '(a,i0)') 'steps with theta outside the nodes: ', cnt%theta_outside_steps
  write(u, '(a,i0)') 'distribution samples: ', samples
  write(u, '(a,f14.4)') 'wall seconds to tau_c = 1: ', wall_tau_c1
  write(u, '(a,f14.4)') 'wall seconds to tau_enc = 1: ', wall_tau_enc1
  write(u, '(a,f14.4)') 'wall seconds total (time march): ', wall_total
  close(u)
end program dsmc0d
