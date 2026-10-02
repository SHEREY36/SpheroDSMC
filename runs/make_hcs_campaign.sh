#!/bin/bash
# HCS campaign for the paper: every exported (alpha, AR) model, 4 seeds of
# N = 250000 (10^6 particles per case in all), to 25 contacts per particle.
# Usage: runs/make_hcs_campaign.sh [JOBS]   (run from the repository root)
JOBS=${1:-2}
mkdir -p runs/hcs
: > runs/hcs/list.txt
for model in models/model_AR*_a*.npz; do
  tag=$(basename $model .npz | sed 's/model_//')
  for s in 1 2 3 4; do
    cat > runs/hcs/${tag}_s$s.nml <<NML
&run
  model_file = '$model'
  output_prefix = 'runs/hcs/${tag}_s$s'
  n_particles = 250000
  seed = $s
  volume_fraction = 0.01
  density = 0.763944
  T_tr0 = 1.0
  T_rot0 = 1.0
  tau_end = 25.0
  output_dtau = 0.05
  sample_tau_start = 5.0
  sample_dtau = 0.5
/
NML
    echo runs/hcs/${tag}_s$s.nml >> runs/hcs/list.txt
  done
done
xargs -a runs/hcs/list.txt -P $JOBS -I{} ./bin/dsmc0d {}
