#!/bin/bash
# Near-sphere HCS (AR = 1.1): exchange takes ~20 encounters, so theta settles
# only after ~50-100 contacts per particle.  4 seeds of N = 250000 to tau_c = 100.
JOBS=${1:-2}
mkdir -p runs/hcs_long
: > runs/hcs_long/list.txt
for model in models/model_AR1.1_a*.npz; do
  tag=$(basename $model .npz | sed 's/model_//')
  for s in 1 2 3 4; do
    cat > runs/hcs_long/${tag}_s$s.nml <<NML
&run
  model_file = '$model'
  output_prefix = 'runs/hcs_long/${tag}_s$s'
  n_particles = 250000
  seed = $((s + 10))
  tau_end = 100.0
  output_dtau = 0.1
  sample_tau_start = 60.0
  sample_dtau = 1.0
/
NML
    echo runs/hcs_long/${tag}_s$s.nml >> runs/hcs_long/list.txt
  done
done
xargs -a runs/hcs_long/list.txt -P $JOBS -I{} ./bin/dsmc0d {}
