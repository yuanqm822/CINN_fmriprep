#!/bin/bash
#SBATCH --job-name=mriqc
#SBATCH --array=1-3                 # 1-N, N = number of lines in subjects.txt
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=logs/mriqc_%A_%a.out
#SBATCH --error=logs/mriqc_%A_%a.err

# ---------------- EDIT THESE ----------------
BIDS=/scratch5/$USER/<project>/bids
OUT=//scratch5/$USER/<project>/derivatives/mriqc
WORK=/scratch5/$USER/mriqc_work
SUBJECTS=subjects.txt

# ---------------- OPTIONS -------------------
# Any extra MRIQC flags, one per item. Examples:
#   (--modalities T1w bold)       only these modalities
#   (--fd_thres 0.5)              framewise-displacement threshold (mm)
#   (--no-datalad-get)            don't try to fetch data with DataLad
EXTRA_ARGS=()
# --------------------------------------------

set -euo pipefail

module use /software/modules/CINN
module load mriqc/26.1.0

SUB=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$SUBJECTS")
if [ -z "$SUB" ]; then
    echo "No subject on line ${SLURM_ARRAY_TASK_ID} of $SUBJECTS" >&2
    exit 1
fi

mkdir -p "$OUT" "$WORK/sub-$SUB"
echo "Running MRIQC for sub-$SUB on $(hostname) at $(date)"
echo "Extra args: ${EXTRA_ARGS[*]} $*"

mriqc "$BIDS" "$OUT" participant \
    --participant-label "$SUB" \
    -w "$WORK/sub-$SUB" \
    --nprocs "$SLURM_CPUS_PER_TASK" \
    --mem-gb 14 \
    --no-sub \
    "${EXTRA_ARGS[@]}" \
    "$@"

echo "Finished sub-$SUB at $(date)"
