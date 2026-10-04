#!/bin/bash
# Campaigns of the paper, run from the repository root with the final closure
# planes in models/ (tools/export_planes.sh on the final model folder):
#
#   runs/make_paper_campaign.sh [JOBS]
#
#  1. HCS against DEM: AR 1.1, 1.5, 2, 2.5, 3 at alpha 0.65, 0.8, 0.95, four
#     realizations of 2.5x10^5 particles, to 25 collisions per particle (100 at
#     AR 1.1, where translation-rotation exchange is slow)        -> runs/hcs_paper
#  2. Long-time distributions at AR 2: ten realizations of 10^5 particles to
#     1500 collisions per particle, sampled every 5 collisions per particle from
#     500 on, with similarity rescaling of the temperature scale  -> runs/hcs_ng
#  3. The serial timing benchmark on the same planes (runs/bench), run alone.
# JOBS concurrent runs (default 5: one core stays free).
set -e
JOBS=${1:-5}
make -s
mkdir -p runs/hcs_paper runs/hcs_ng
LIST=runs/paper_list.txt
: > $LIST
for AR in 1.1 1.5 2.0 2.5 3.0; do
  for a in 0.65 0.80 0.95; do
    model=models/model_AR${AR}_a${a}.npz
    [ -f $model ] || { echo "missing $model"; exit 1; }
    tau_end=25.0
    [ "$AR" = "1.1" ] && tau_end=100.0
    for s in 1 2 3 4; do
      p=runs/hcs_paper/AR${AR}_a${a}_s$s
      [ -f $p.log ] && grep -q "wall seconds total" $p.log && continue
      cat > $p.nml <<NML
&run
  model_file = '$model'
  output_prefix = '$p'
  n_particles = 250000
  seed = $s
  volume_fraction = 0.01
  density = 0.763944
  tau_end = $tau_end
  output_dtau = 0.05
/
NML
      echo $p.nml >> $LIST
    done
  done
done
for a in 0.65 0.80 0.95; do
  for s in $(seq 1 10); do
    p=runs/hcs_ng/AR2.0_a${a}_s$s
    [ -f $p.log ] && grep -q "wall seconds total" $p.log && continue
    cat > $p.nml <<NML
&run
  model_file = 'models/model_AR2.0_a${a}.npz'
  output_prefix = '$p'
  n_particles = 100000
  seed = $((100 + s))
  volume_fraction = 0.01
  density = 0.763944
  tau_end = 1500.0
  output_dtau = 1.0
  sample_tau_start = 500.0
  sample_dtau = 5.0
  rescale_below = 1.0e-6
/
NML
    echo $p.nml >> $LIST
  done
done
# long runs first so that the short ones fill the gaps
{ grep hcs_ng $LIST; grep -v hcs_ng $LIST; } > $LIST.ordered
echo "$(wc -l < $LIST.ordered) runs, $JOBS at a time"
xargs -a $LIST.ordered -P $JOBS -I{} sh -c './bin/dsmc0d {} > /dev/null'

ARS="1.1 1.5 2.0 2.5 3.0" ALPHAS="0.65 0.80 0.95" runs/bench_serial.sh
