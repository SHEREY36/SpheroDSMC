# SpheroDSMC

Zero-dimensional direct simulation Monte Carlo (DSMC) of a dilute gas of
inelastic, smooth spherocylinders, with the collision-data-informed relaxation
closures of

> M. Gbolasere and A. Morris, *Collision-data-informed relaxation closures for
> DSMC of granular gases of non-spherical particles* (submitted to Physics of Fluids).

The code is the reference implementation of Algorithm 3 of the paper. It is
small on purpose (about 1100 lines of Fortran 2008, no dependencies) so that every
step of the method can be read next to its equation.

## Build and run

```bash
make                                   # gfortran -O3; `make debug` for a checked build
./bin/dsmc0d examples/hcs_AR2.0_a0.80.nml
```

The example follows the homogeneous cooling state (HCS) of rods with aspect
ratio AR = 2 and restitution α = 0.8 for 25 collisions per particle with
2.5×10⁵ simulator particles (about one minute on one core).

## What the program does

| paper | code |
|---|---|
| encounter cross-section σ_c(θ) = Ā⊥ D(θ), Eq. (8) | `closure_model::dynamic_factor`, `particles::setup_gas` |
| speed screen \|g·e\| ≥ U v_max (Sec. 4.2) | `collisions::collision_step` |
| orientation and spin screen A⊥ G(Ξ)/(A_sup G_max), Eqs. (7), (9), (10) | `collisions::projected_area`, `closure_model::spin_factor` |
| candidate count (no-time-counter), Eq. (B.8) | `collisions::collision_step` |
| loss ε = [c_t z + c_r (1 − z)] B/⟨B⟩, B ~ Beta(ζ_n/2, ζ_s/2), Eqs. (14)–(15) | `collisions::encounter` |
| exchange z′ = Σ_k w_k Q_k(U; a_k), Eqs. (16), (17), (19) | `closure_model::sample_partition` |
| rotational split and spin reset (Sec. 5.6) | `collisions::encounter` |
| angular law exp[a cos χ + b P₂(cos χ)], Eq. (18) | `closure_model::sample_cosine` |
| velocity update, Eq. (F.2) and Algorithm 3 | `collisions::encounter` |
| free rotation (exact Rodrigues, Strang split) | `particles::rotate_axes` |

At fixed α and AR every closure depends on the gas only through the temperature
ratio θ = T_tr/T_rot. The closures are tabulated at a few θ nodes and combined
with linear weights in θ (the exchange law by quantile averaging, the angular
law as a mixture). The time step is chosen each step so that a particle has a
probability `collision_fraction` (default 0.02) of an encounter per step.

### Clocks

One DSMC event is one binary *encounter*: everything from first contact to final
separation, possibly several contacts. The output gives both

* `tau_enc`: encounters per particle (the DSMC clock), and
* `tau_c`: contacts per particle, the sum of the mean contacts per encounter
  ⟨k⟩(θ) measured in the trajectory calculations. This is what a DEM simulation
  counts, and `tau_end` is expressed in it.

## Input (`namelist /run/`)

| name | default | meaning |
|---|---|---|
| `model_file` | | closure model `.npz` (one α, one AR) |
| `output_prefix` | `run` | prefix of the output files |
| `n_particles` | 100000 | simulator particles |
| `seed` | 1 | random seed |
| `volume_fraction` | 0.01 | sets the number density n = φ/V_p |
| `density` | 0.763944 | material density (m = ρ V_p; m = 1 at AR = 2) |
| `diameter` | 1 | rod diameter d |
| `T_tr0`, `T_rot0` | 1, 1 | initial temperatures |
| `exact_start` | `.true.` | rescale the initial draw to exactly T_tr0, T_rot0 |
| `tau_end`, `t_end` | 25, ∞ | stop at this many contacts per particle or this time |
| `collision_fraction` | 0.02 | encounters per particle per time step |
| `output_dtau` | 0.05 | contacts per particle between time-series rows |
| `sample_tau_start`, `sample_dtau` | ∞, 0.5 | distribution sampling |
| `rescale_below` | 0 (off) | rescale all velocities and spins back to the initial temperature whenever T falls below this fraction of it (long cooling runs; the time column is then in rescaled units) |
| `flow` | `hcs` | `hcs`: homogeneous cooling; `usf`: uniform shear flow |
| `shear_rate` | 0 | USF shear rate γ̇ |
| `stress_tau_start` | ∞ | USF: steady averages over τ_c ≥ this |
| `max_shear_step` | 0.01 | USF: largest strain γ̇Δt per time step |

### Uniform shear flow

With `flow = 'usf'` the peculiar velocities are sheared in the homogeneous (Lagrangian)
frame, c_x ← c_x − γ̇ c_y Δt, exactly; axes and spins are unaffected. The shear map is
split around the collisions in half steps, like the free rotation. The time step is the
smaller of the collision step and `max_shear_step`/γ̇. The exact shear work and the
energy change of the collisions are accumulated, so the energy ledger closes to round-off,
and in the steady state the shear work balances the collisional loss.

## Output

* `<prefix>.ts`: t, tau_enc, tau_c, T_tr, T_rot, θ, T = (3T_tr + 2T_rot)/5,
  a₂ of the translational and rotational distributions, dt.
* `<prefix>.speed`, `<prefix>.rot`: distributions of c* = |c − U|/√(2T_tr/m)
  and x = E_rot/T_rot with the Maxwellian references.
* `<prefix>.stress` (USF): t, τ_c, T_tr, T_rot, θ, T* = T_tr/(m γ̇² d²), the kinetic
  stress P^k_ij/(n T_tr) (xx, yy, zz, xy, xz, yz), and the energy ledger (thermal energy,
  cumulative shear work and collisional change, bookkeeping residual).
* `<prefix>.log`: parameters, counters (candidates, encounters, majorant
  violations, table clamps), the contact model of the closure and wall-clock times,
  including the time to the first collision per particle. USF runs add the steady
  means of T*, θ, P^k_ij, N₁, N₂ and a₂ with 10-block standard errors, the contacts per
  particle per unit strain, and the energy balance over the averaging window.

## Model files

`models/` holds the final closures of the paper for AR = 1.1, 1.5, 2, 2.5, 3 and
α = 0.65, 0.8, 0.95 (about 3 MB each). They are exported from the released model
folder (`closure_v2.npz` with its three tables) with

```bash
python3 tools/export_model.py --artifact closure_v2.npz \
    --angular angular_memory.json --cross-section encounter_cross_section.json \
    --loss-memory loss_memory.json --alpha 0.8 --AR 2.0 --out models/model_AR2.0_a0.80.npz
```

A model file is an uncompressed NumPy archive (`np.savez`); `src/npz_io.f90`
reads it directly. Its contents are listed at the top of `tools/export_model.py`.
`tools/export_planes.sh model_dir [AR ...]` exports every fitted plane of a model
folder (restricted to some restitutions with `ALPHAS="0.65 0.8 0.95"`). With the
loss table the loss of an encounter follows the incoming partition,
E[ε|z] = c_t z + c_r (1 − z); model files without it use the unconditional node mean.

## Reproducing the paper

```bash
runs/make_paper_campaign.sh            # HCS against DEM, long-time distributions, benchmark
```

runs four realizations of 2.5×10⁵ particles at every (AR, α) of `models/` to 25
collisions per particle (100 at AR = 1.1), ten realizations of 10⁵ particles at
AR = 2 to 1500 collisions per particle with the distributions sampled from 500 on,
and times the first collision per particle on one core (`runs/bench/summary.txt`).
The results land in `runs/hcs_paper/` and `runs/hcs_ng/`. On five cores the
campaign takes about two hours.

### Contact model C1

`runs/make_c1_campaign.sh` runs the campaigns of the contact-model C1 closure (damping of
the normal velocity of the contact points) from planes exported with
`CONTACT_MODEL=C1 OUT_DIR=models_C1 ELASTIC=1 tools/export_planes.sh <model folder>`:
HCS against DEM, the elastic-equilibrium test on the learned α = 1 planes, the uniform
shear flow at the 40 DEM cases of `runs/usf_protocol.csv` (from a cold and a hot start)
and the long-time distributions. `tools/usf_summary.py runs/C1/usf --dem <DEM table>`
reduces the USF runs (cold and hot runs of one seed are averaged first).

## Verification

`tests/compare_reference.py` compares an ensemble of runs with the original
Python implementation of the method (DSMC_V2). With the final model, the
asymptotic temperature ratio θ* of the Fortran code agrees with the Python
validation campaign to within 0.8 % at AR = 1.1 (α = 0.8), 1.5 (0.65), 2 (0.8)
and 3 (0.95).
