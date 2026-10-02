#!/bin/bash
# HCS campaign of the paper from one closure model directory (the layout of
# DSMC_V2/models/<dir>: closure_v2.npz, angular_memory.json,
# encounter_cross_section.json). Run from the repository root:
#
#   runs/make_paper_campaign.sh /path/to/model_dir [JOBS]
#
#  1. exports AR {1.25, 1.5, 2, 2.5, 3} x alpha {0.6, 0.8, 0.95} to models/, and
#     AR 1.25 at the Rubio-Largo restitutions 0.88, 0.90, 0.96;
#  2. runs 4 x 250000 particles of every case to 25 contacts per particle
#     (100 for AR 1.25, where translation-rotation exchange is slow);
#  3. reruns the serial timing benchmark on the same models (runs/bench).
# Results land in runs/hcs_paper/, which the paper's analysis script reads.
# JOBS concurrent runs (default 5: one core stays free); the benchmark runs alone.
set -e
DIR=${1:?usage: runs/make_paper_campaign.sh model_dir [JOBS]}
JOBS=${2:-5}
make -s

ALPHAS="0.6 0.8 0.95" tools/export_planes.sh "$DIR" 1.25 1.5 2.0 2.5 3.0
ALPHAS="0.88 0.9 0.96" tools/export_planes.sh "$DIR" 1.25

OUT=runs/hcs_paper
mkdir -p $OUT
: > $OUT/list.txt
cases="1.25:0.60 1.25:0.80 1.25:0.95 1.25:0.88 1.25:0.90 1.25:0.96"
for AR in 1.5 2.0 2.5 3.0; do
  for a in 0.60 0.80 0.95; do cases="$cases $AR:$a"; done
done
for c in $cases; do
  AR=${c%%:*}; a=${c##*:}
  model=models/model_AR${AR}_a${a}.npz
  [ -f $model ] || { echo "missing $model (plane not in the artifact?)"; continue; }
  tau_end=25.0
  [ "$AR" = "1.25" ] && tau_end=100.0
  for s in 1 2 3 4; do
    p=$OUT/AR${AR}_a${a}_s$s
    [ -f $p.log ] && grep -q "wall seconds" $p.log && continue      # finished
    cat > $p.nml <<NML
&run
  model_file = '$model'
  output_prefix = '$p'
  n_particles = 250000
  seed = $s
  volume_fraction = 0.01
  density = 0.763944
  T_tr0 = 1.0
  T_rot0 = 1.0
  tau_end = $tau_end
  output_dtau = 0.05
  sample_tau_start = 5.0
  sample_dtau = 0.5
/
NML
    echo $p.nml >> $OUT/list.txt
  done
done
echo "$(wc -l < $OUT/list.txt) runs, $JOBS at a time"
xargs -a $OUT/list.txt -P $JOBS -I{} sh -c './bin/dsmc0d {} > /dev/null'

ARS="1.25 1.5 2.0 2.5 3.0" ALPHAS="0.60 0.80 0.95" runs/bench_serial.sh
