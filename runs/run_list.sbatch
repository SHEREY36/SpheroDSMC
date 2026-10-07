#!/bin/bash
# One SpheroDSMC run per array task: line SLURM_ARRAY_TASK_ID + 1 of a list of .nml files.
#   submitted by runs/make_c1_campaign.sh with RUNNER=slurm (from the repository root)
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1900M
LIST=${1:?list of namelist files}
cd "${SLURM_SUBMIT_DIR:?}"
nml=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "$LIST")
[ -n "$nml" ] || exit 0
p=${nml%.nml}
if [ -f "$p.log" ] && grep -q "wall seconds total" "$p.log"; then exit 0; fi
./bin/dsmc0d "$nml" > /dev/null
grep -q "wall seconds total" "$p.log"
