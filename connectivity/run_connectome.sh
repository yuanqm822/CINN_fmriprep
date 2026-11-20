#!/bin/bash
#SBATCH --job-name=connectome
#SBATCH --array=1-3                 # 1-N, N = number of lines in subjects.txt
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=1:00:00
#SBATCH --output=logs/connectome_%A_%a.out
#SBATCH --error=logs/connectome_%A_%a.err

# ---------------- EDIT THESE ----------------
FMRIPREP=/scratch5/$USER/<project>/derivatives/fmriprep
OUT=/scratch5/$USER/<project>/derivatives/connectome
SUBJECTS=/scratch5/$USER/<project>/subjects.txt
# --------------------------------------------

# ---------------- OPTIONS -------------------
# Denoising options passed to connectome.py. Examples:
#   (--strategy baseline)                            high-pass only (reference)
#   (--strategy simple)                              motion + WM/CSF + high-pass
#   (--strategy simple --global-signal basic)        ... + global signal regression
#   (--strategy scrubbing --fd-threshold 0.5)        volume censoring
#   (--strategy compcor)                             anatomical CompCor
#   (--atlas 200 --task rest)                        Schaefer 200, task label
OPTIONS=(--task rest --strategy simple --global-signal basic --atlas 400)
# --------------------------------------------

set -euo pipefail

module use /software/modules/CINN
module load nilearn/1.0

SUB=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$SUBJECTS")
if [ -z "$SUB" ]; then
    echo "No subject on line ${SLURM_ARRAY_TASK_ID} of $SUBJECTS" >&2
    exit 1
fi

echo "Connectome for sub-$SUB on $(hostname) at $(date)"
echo "Options: ${OPTIONS[*]} $*"

python "${SLURM_SUBMIT_DIR:-.}/connectome.py" "$FMRIPREP" "$SUB" "$OUT" "${OPTIONS[@]}" "$@"
