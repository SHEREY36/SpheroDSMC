#!/bin/bash
# Export every fitted (alpha < 1, AR) plane of a closure model directory.
#   tools/export_planes.sh /path/to/model_dir [AR ...]
# model_dir must hold closure_v2.npz, angular_memory.json and
# encounter_cross_section.json (the layout of DSMC_V2/models/*).
# Without AR arguments all aspect ratios with a full plane are exported; the
# environment variable ALPHAS (e.g. ALPHAS="0.6 0.8 0.95") limits the restitutions.
set -e
DIR=${1:?usage: tools/export_planes.sh model_dir [AR ...]}
shift
PY=${PYTHON:-python3}
ART=$DIR/closure_v2.npz
ANG=$DIR/angular_memory.json
CS=$DIR/encounter_cross_section.json
LOSS=""
[ -f $DIR/loss_memory.json ] && LOSS="--loss-memory $DIR/loss_memory.json"
PLANES=$(ALPHAS="${ALPHAS:-}" $PY - "$ART" "$@" <<'EOF'
import os, sys, numpy as np
c = np.load(sys.argv[1])["surface_coordinates"]
want = [float(a) for a in sys.argv[2:]]
alphas = [float(a) for a in os.environ["ALPHAS"].split()]
for a, ar in sorted({(round(x[0], 4), round(x[2], 4)) for x in c}):
    if a < 1.0 and (not want or any(abs(ar - w) < 1e-9 for w in want)) \
            and (not alphas or any(abs(a - w) < 1e-9 for w in alphas)):
        tag = f"{ar:.2f}"
        tag = tag[:-1] if tag.endswith("0") else tag          # 2.0, 1.25, 1.1
        print(a, ar, f"model_AR{tag}_a{a:.2f}.npz")
EOF
)
mkdir -p models
while read a ar name; do
  $PY tools/export_model.py --artifact $ART --angular $ANG --cross-section $CS $LOSS \
      --alpha $a --AR $ar --out models/$name \
      || echo "skipped alpha=$a AR=$ar"
done <<< "$PLANES"
