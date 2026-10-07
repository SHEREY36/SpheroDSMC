#!/bin/bash
# DSMC campaigns of the contact-model C1 closure, run from the repository root:
#
#   CONTACT_MODEL=C1 OUT_DIR=models_C1 ELASTIC=1 tools/export_planes.sh <C1 model folder>
#   runs/make_c1_campaign.sh [JOBS]            (PARTS="hcs elastic ng" by default)
#
#  hcs      HCS against the C1 DEM: AR 1.1, 1.5, 2, 2.5, 3 at alpha 0.5, 0.65, 0.8, 0.95,
#           eight realizations of 2.5x10^5 particles to 25 collisions per particle
#           (100 at AR 1.1, where translation-rotation exchange is slow)   -> $OUT/hcs_paper
#  elastic  elastic equilibrium on the learned alpha = 1 planes at AR 1.1, 2, 3: four
#           realizations of 10^5 particles to 200 collisions per particle from
#           theta = 1, 0.2 and 5 at T = (3 T_tr + 2 T_rot)/5 = 1     -> $OUT/elastic
#  usf      uniform shear flow at the 40 DEM cases of runs/usf_protocol.csv (AR 1.5, 2, 2.5, 3 x
#           alpha 0.5..0.95; shear rate, cold and hot starts, run length and averaging window
#           of the Python validation protocol usf-encounter-v1): four realizations of 10^5
#           particles per start                                       -> $OUT/usf
#  ablation the five ablations of collisions.f90/closure_model.f90 (no refit): HCS at AR 1.1, 2, 3
#           x alpha 0.65, 0.8, four realizations of 2.5x10^5; USF at AR 2 x alpha 0.65, 0.8,
#           0.95, two realizations of 10^5 per start                  -> $OUT/ablation
#  ng       long-time distributions at AR 2: ten realizations of 10^5 particles to 1500
#           collisions per particle, sampled every 5 from 500 on, similarity rescaling
#                                                                    -> $OUT/hcs_ng
#  bench    serial timing (runs/bench_serial.sh on $MODELS); run it alone, not with others
#
# MODELS (default models_C1) and OUT (default runs/C1) select the planes and the output.
# RUNNER=slurm submits the runs as one Slurm array (one core per run; ACCOUNT, QOS, PARTITION,
# RUN_TIME from the environment) instead of running them here; check the logs once the array
# has finished with CHECK_ONLY=1 runs/make_c1_campaign.sh.
# Finished runs (log with "wall seconds total") are skipped, so the script can be rerun.
# Every run must report "contact model: C1"; the script checks the logs at the end.
set -e
JOBS=${1:-10}
MODELS=${MODELS:-models_C1}
OUT=${OUT:-runs/C1}
PARTS=${PARTS:-"hcs elastic ng"}
RUNNER=${RUNNER:-local}
[ -n "${CHECK_ONLY:-}" ] && PARTS=""
make -s
LIST=$OUT/list.txt
mkdir -p "$OUT"
: > "$LIST"

nml() {   # nml PREFIX MODEL N SEED TAU_END [extra namelist lines]
  local p=$1 model=$2 n=$3 seed=$4 tau=$5
  shift 5
  [ -f "$model" ] || { echo "missing $model" >&2; exit 1; }
  [ -f "$p.log" ] && grep -q "wall seconds total" "$p.log" && return 0
  {
    echo "&run"
    echo "  model_file = '$model'"
    echo "  output_prefix = '$p'"
    echo "  n_particles = $n"
    echo "  seed = $seed"
    echo "  volume_fraction = 0.01"
    echo "  density = 0.763944"
    echo "  tau_end = $tau"
    for line in "$@"; do echo "  $line"; done
    echo "/"
  } > "$p.nml"
  echo "$p.nml" >> "$LIST"
}

for part in $PARTS; do
  case $part in
    ng)
      mkdir -p "$OUT/hcs_ng"
      for a in 0.65 0.80 0.95; do
        for s in $(seq 1 10); do
          nml "$OUT/hcs_ng/AR2.0_a${a}_s$s" "$MODELS/model_AR2.0_a${a}.npz" 100000 $((100 + s)) 1500.0 \
              "output_dtau = 1.0" "sample_tau_start = 500.0" "sample_dtau = 5.0" "rescale_below = 1.0e-6"
        done
      done ;;
    elastic)
      mkdir -p "$OUT/elastic"
      for AR in 1.1 2.0 3.0; do
        # theta0 T_tr0 T_rot0 with (3 T_tr0 + 2 T_rot0)/5 = 1
        for start in "1.0 1.0 1.0" "0.2 0.38461538 1.92307692" "5.0 1.47058824 0.29411765"; do
          set -- $start
          for s in 1 2 3 4; do
            nml "$OUT/elastic/AR${AR}_theta$1_s$s" "$MODELS/model_AR${AR}_a1.00.npz" 100000 $((200 + s)) 200.0 \
                "T_tr0 = $2" "T_rot0 = $3" "output_dtau = 0.5"
          done
        done
      done ;;
    usf)
      mkdir -p "$OUT/usf"
      tail -n +4 runs/usf_protocol.csv | while IFS=, read -r ar a gd tc rc th rh tend tstart; do
        for br in cold hot; do
          if [ $br = cold ]; then t0=$tc; r0=$rc; else t0=$th; r0=$rh; fi
          for s in 1 2 3 4; do
            nml "$OUT/usf/AR${ar}_a${a}_${br}_s$s" "$MODELS/model_AR${ar}_a${a}.npz" 100000 $((300 + s)) "$tend" \
                "flow = 'usf'" "shear_rate = $gd" "T_tr0 = $t0" "T_rot0 = $r0" \
                "stress_tau_start = $tstart" "output_dtau = 0.05"
          done
        done
      done ;;
    ablation)
      mkdir -p "$OUT/ablation"
      for ab in no_spin_screen loss_unconditioned lambda4_off isotropic_scattering memoryless_exchange; do
        for AR in 1.1 2.0 3.0; do
          for a in 0.65 0.80; do
            tau=25.0
            [ "$AR" = 1.1 ] && tau=100.0
            for s in 1 2 3 4; do
              nml "$OUT/ablation/hcs_${ab}_AR${AR}_a${a}_s$s" "$MODELS/model_AR${AR}_a${a}.npz" 250000 $s $tau \
                  "output_dtau = 0.05" "ablation = '$ab'"
            done
          done
        done
        grep "^2.0," runs/usf_protocol.csv | while IFS=, read -r ar a gd tc rc th rh tend tstart; do
          case $a in 0.65|0.80|0.95) ;; *) continue ;; esac
          for br in cold hot; do
            if [ $br = cold ]; then t0=$tc; r0=$rc; else t0=$th; r0=$rh; fi
            for s in 1 2; do
              nml "$OUT/ablation/usf_${ab}_AR${ar}_a${a}_${br}_s$s" "$MODELS/model_AR${ar}_a${a}.npz" 100000 \
                  $((300 + s)) "$tend" "flow = 'usf'" "shear_rate = $gd" "T_tr0 = $t0" "T_rot0 = $r0" \
                  "stress_tau_start = $tstart" "output_dtau = 0.05" "ablation = '$ab'"
            done
          done
        done
      done ;;
    hcs)
      mkdir -p "$OUT/hcs_paper"
      for AR in 1.1 1.5 2.0 2.5 3.0; do
        for a in 0.50 0.65 0.80 0.95; do
          tau=25.0
          [ "$AR" = 1.1 ] && tau=100.0
          for s in $(seq 1 8); do
            nml "$OUT/hcs_paper/AR${AR}_a${a}_s$s" "$MODELS/model_AR${AR}_a${a}.npz" 250000 $s $tau \
                "output_dtau = 0.05"
          done
        done
      done ;;
    bench)
      sed "s#models/model_#$MODELS/model_#g; s#runs/bench#$OUT/bench#g" runs/bench_serial.sh \
        | ARS="1.1 1.5 2.0 2.5 3.0" ALPHAS="0.65 0.80 0.95" bash ;;
    *) echo "unknown part $part" >&2; exit 2 ;;
  esac
done

if [ -s "$LIST" ]; then
  # long runs first so that the short ones fill the gaps
  { grep hcs_ng "$LIST" || true; grep elastic "$LIST" || true; grep usf "$LIST" || true;
    grep -v "hcs_ng\|elastic\|usf" "$LIST" || true; } \
    > "$LIST.ordered"
  n=$(wc -l < "$LIST.ordered")
  if [ "$RUNNER" = slurm ]; then
    stamp=$(date +%Y%m%dT%H%M%S)
    cp "$LIST.ordered" "$OUT/list_$stamp.txt"
    mkdir -p "$OUT/slurm"
    sbatch -A "${ACCOUNT:-morri353}" -p "${PARTITION:-cpu}" -q "${QOS:-normal}" -J sdsmc_c1 \
      --time="${RUN_TIME:-01:00:00}" --array=0-$((n - 1)) -o "$OUT/slurm/%A_%a.out" \
      runs/run_list.sbatch "$OUT/list_$stamp.txt"
    echo "$n runs submitted as one array; when it has finished: CHECK_ONLY=1 OUT=$OUT runs/make_c1_campaign.sh"
    exit 0
  fi
  echo "$n runs, $JOBS at a time"
  xargs -a "$LIST.ordered" -P "$JOBS" -I{} sh -c './bin/dsmc0d {} > /dev/null'
fi

# every finished run must be contact model C1 with no loss cap
bad=0
for log in $(ls "$OUT"/*/*.log 2>/dev/null); do
  grep -q "wall seconds total" "$log" || { echo "unfinished: $log"; bad=1; continue; }
  grep -q "contact model: C1" "$log" || { echo "not contact model C1: $log"; bad=1; }
  grep -q "loss caps: 0" "$log" || { echo "loss caps: $log"; bad=1; }
  case $log in
    */ablation/*) ;;
    *) if grep -q "^ablation: " "$log" && ! grep -q "^ablation: none" "$log"; then
         echo "ablation set: $log"; bad=1
       fi ;;
  esac
done
[ $bad = 0 ] && echo "all runs finished, contact model C1, no loss caps"
exit $bad
