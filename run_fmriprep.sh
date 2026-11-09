#!/bin/bash
#SBATCH --job-name=fmriprep
#SBATCH --array=1-3                 # 1-N, N = number of lines in subjects.txt
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=logs/fmriprep_%A_%a.out
#SBATCH --error=logs/fmriprep_%A_%a.err

# ---------------- EDIT THESE ----------------
BIDS=/scratch5/$USER/<project>/bids
OUT=/scratch5/$USER/<project>/derivatives/fmriprep
WORK=/scratch5/$USER/fmriprep_work
SUBJECTS=subjects.txt
# --------------------------------------------

set -euo pipefail

module use /software/modules/CINN
module load fmriprep/25.2.7

SUB=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$SUBJECTS")
if [ -z "$SUB" ]; then
    echo "No subject on line ${SLURM_ARRAY_TASK_ID} of $SUBJECTS" >&2
    exit 1
fi

mkdir -p "$OUT" "$WORK/sub-$SUB"
echo "Running fMRIPrep for sub-$SUB on $(hostname) at $(date)"

fmriprep "$BIDS" "$OUT" participant \
    --participant-label "$SUB" \
    -w "$WORK/sub-$SUB" \
    --nprocs "$SLURM_CPUS_PER_TASK" \
    --mem-mb 30000 \
    --output-spaces MNI152NLin2009cAsym:res-2 T1w \
    --skip-bids-validation

echo "Finished sub-$SUB at $(date)"
