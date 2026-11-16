#!/bin/bash
#SBATCH --job-name=group_metrics
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=logs/group_metrics_%j.out
#SBATCH --error=logs/group_metrics_%j.err

# Group-level denoising benchmark (Wang et al. 2024): loss of tDOF, QC-FC,
# DM-FC and network modularity for every strategy in the connectome folder.
# Run once, after motion_qc.py and run_connectome.sh have finished for ALL subjects.

# ---------------- EDIT THESE ----------------
CONNECTOMES=/scratch5/$USER/<project>/derivatives/connectome
MOTION_QC=/scratch5/$USER/<project>/derivatives/fmriprep/motion_qc_task-rest.tsv
PARTICIPANTS=/scratch5/$USER/<project>/bids/participants.tsv
AGE_COL=""          # leave empty to use the first column starting with "age"
SEX_COL=sex
# --------------------------------------------

set -euo pipefail

module use /software/modules/CINN
module load nilearn/1.0

ARGS=(--sex-col "$SEX_COL")
[ -n "$AGE_COL" ] && ARGS+=(--age-col "$AGE_COL")

python "${SLURM_SUBMIT_DIR:-.}/group_metrics.py" "$CONNECTOMES" "$MOTION_QC" "$PARTICIPANTS" \
    "${ARGS[@]}" "$@"
