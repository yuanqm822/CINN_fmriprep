# Functional connectivity with Nilearn (in development)

> This chapter extends the [main guide](../README.md): it starts from fMRIPrep outputs
> and is still being developed. Scripts and defaults may change — check the
> [roadmap](#roadmap) and contact Simon before relying on it for a study.

All commands below are run from this folder:

```bash
cd CINN_fmriprep/connectivity
mkdir -p logs
```

The job scripts read the subject list from `../subjects.txt` (the same file used for
fMRIPrep); change `SUBJECTS=` in a script to use another list.

## Overview

After fMRIPrep, `run_connectome.sh` denoises the preprocessed resting-state data and
computes a functional connectome for each run, following the workflow of
Wang et al. (2024). It uses `nilearn/1.0`, a CINN Python environment (built with
pixi) that contains Nilearn and TemplateFlow, so there is no conda environment to
set up. The Schaefer 2018 atlases are pre-downloaded to the shared TemplateFlow
folder.

For each run it:

1. selects confounds from fMRIPrep's `desc-confounds_timeseries.tsv` with
   `load_confounds_strategy`
2. extracts denoised Schaefer-atlas time series (`NiftiLabelsMasker`)
3. computes a correlation matrix (`ConnectivityMeasure`)

It needs fMRIPrep outputs in `MNI152NLin2009cAsym:res-2`, which the default
`OUTPUT_SPACES` in `run_fmriprep.sh` already produces.

## Step 1: motion quality control

Before computing connectomes, check which runs pass the motion criteria used by
Wang et al. (2024). A run is flagged if mean framewise displacement (FD) is above
0.25 mm, if scrubbing at FD 0.2 mm removes more than 80% of volumes, or if less than
1 minute of data is left. It only reads fMRIPrep's confounds files, so it runs in
seconds on the login node:

```bash
module load nilearn/1.0
python motion_qc.py /scratch5/$USER/<project>/derivatives/fmriprep --task rest
```

It prints a summary and writes `motion_qc_task-rest.tsv` (mean/max FD, % volumes
removed, minutes left, `exclude` and the reason) into the fMRIPrep folder. Leave
flagged runs out of later analyses. The paper's data had a TR of 2–2.5 s; with
multiband data (e.g. TR 0.8 s), respiration inflates FD, so check the numbers rather
than relying only on the flag.

## Step 2: connectomes

Edit the paths under `EDIT THESE` in `run_connectome.sh`, set `--array`, then:

```bash
mkdir -p logs
sbatch run_connectome.sh
```

Choose the denoising strategy in the `OPTIONS` line:

| Strategy | `OPTIONS` | What it removes |
|----------|-----------|-----------------|
| baseline | `--strategy baseline` | low-frequency drift only (the paper's reference) |
| simple + GSR (default) | `--strategy simple --global-signal basic` | motion, WM/CSF, global signal, low-frequency drift |
| simple | `--strategy simple` | as above, without global signal |
| scrubbing | `--strategy scrubbing --fd-threshold 0.5` | full motion + WM/CSF, plus censoring of volumes with FD > 0.5 mm |
| scrubbing (strict) | `--strategy scrubbing --fd-threshold 0.2` | as above with FD > 0.2 mm |
| CompCor | `--strategy compcor` | motion and anatomical CompCor components |

Other options: `--atlas 100|200|300|400|...` (Schaefer parcels, 7 networks) and
`--task` (default `rest`). One-off changes can also be passed at submission, e.g.
`sbatch run_connectome.sh --atlas 200`.

Outputs in your connectome folder, per run:

- `*_timeseries.tsv` — denoised time series, one column per region
- `*_connectome.tsv` — region × region correlation matrix, with region names
- `*_desc.json` — strategy, number of confound regressors and volumes kept

To compare strategies (needed for Step 3), submit the job once per strategy. Flags
given at submission override the `OPTIONS` line:

```bash
for s in "--strategy baseline" \
         "--strategy simple" "--strategy simple --global-signal basic" \
         "--strategy scrubbing --fd-threshold 0.5" "--strategy scrubbing --fd-threshold 0.2" \
         "--strategy compcor"; do
    sbatch run_connectome.sh $s
done
```

## Step 3: group benchmark (whole sample)

Once motion QC and connectomes are done for **all** subjects, `run_group_metrics.sh`
compares the strategies with the four metrics from Wang et al. (2024):

| Metric | What it measures | Better |
|--------|------------------|--------|
| Loss of tDOF | confound regressors + censored volumes | lower |
| QC-FC | partial correlation between subjects' mean FD and each edge, controlling for age and sex (% significant edges, uncorrected and FDR; median \|r\|) | lower |
| DM-FC | Spearman correlation between region distance and QC-FC | closer to 0 |
| Modularity | mean Louvain modularity of each subject's connectome | higher |

Runs flagged by `motion_qc.py` are left out, and a subject's runs are averaged.
Age and sex come from the BIDS `participants.tsv`; the script uses the first column
whose name starts with `age` (e.g. `age@ses-01` from BIDScoin) unless you set
`AGE_COL`. Edit the paths at the top of `run_group_metrics.sh`, then:

```bash
sbatch run_group_metrics.sh
```

Outputs in `<connectome folder>/group/`:

- `group_metrics.tsv` — one row per atlas × strategy, with ranks per metric (as in the paper's Fig 12)
- `group_metrics.png` — bar charts of each metric
- `qcfc_<atlas>_<strategy>.tsv` — per-edge QC-FC, distance and p-value

These are group statistics: with fewer than ~20 subjects the script warns that the
numbers are unreliable. As in the paper, no single "best" strategy is computed —
weigh motion removal (QC-FC, DM-FC) against data lost (tDOF), and choose before
looking at your effects of interest.

Notes:

- **Simple + GSR** keeps every volume, which you need for any analysis that relies on
  a continuous time series. Wang et al. found it to be a good choice in that case.
- **Scrubbing with multiband data:** at short TRs, respiration inflates framewise
  displacement, so a strict threshold (e.g. 0.2 mm) can remove a large share of
  volumes. The script defaults to 0.5 mm and warns if more than half of a run is
  removed. Check `volumes_kept` in `*_desc.json`.
- **Scrubbing uses FD only**, as in the paper (no DVARS threshold). Add
  `--dvars-threshold 1.5` to also censor on DVARS. Runs with fewer than 10 volumes
  left after scrubbing are skipped.
- **ICA-AROMA** (in the paper) is not available: fMRIPrep stopped producing it from
  version 23.1.
- To adapt the analysis, copy `connectome.py` and edit it. Run it with
  `module load nilearn/1.0` then `python connectome.py --help`.

<details>
<summary>For the CINN admin: how <code>nilearn/1.0</code> is installed</summary>

The environment is defined in [`env/pixi.toml`](env/pixi.toml) and
built with CINN's own pixi under `/software/CINN/pixi`:

```bash
export PIXI_CACHE_DIR=/scratch5/$USER/pixi_cache   # keep the cache out of /home
mkdir -p /software/CINN/envs/nilearn-1.0
cp env/pixi.toml /software/CINN/envs/nilearn-1.0/
cd /software/CINN/envs/nilearn-1.0 && /software/CINN/pixi/bin/pixi install
chmod -R a+rX /software/CINN/envs/nilearn-1.0
```

The module `/software/modules/CINN/nilearn/1.0` puts
`/software/CINN/envs/nilearn-1.0/.pixi/envs/default/bin` on `PATH` and sets
`TEMPLATEFLOW_HOME=/software/CINN/templateflow`. New versions go in a new folder
(e.g. `nilearn-1.1`) so existing analyses keep the same package versions.
</details>

## Roadmap

Planned or possible next steps (edit as the chapter develops):

- [ ] Test on the first CINN group dataset and settle default thresholds for multiband data (TR 0.8 s)
- [ ] More atlases from the paper: Gordon 333, MIST, DiFuMo (probabilistic, `NiftiMapsMasker`)
- [ ] Optional physiological / respiration-related confounds for multiband data
- [ ] Task-based connectivity (e.g. task-regressed residuals)
- [ ] Packaged option: giga-connectome BIDS app as a CINN module

## References

- Wang H-T et al. (2024) Continuous evaluation of denoising strategies in resting-state fMRI connectivity using fMRIPrep and Nilearn. *PLOS Comput Biol* 20(3): e1011942. <https://doi.org/10.1371/journal.pcbi.1011942>
- Nilearn: <https://nilearn.github.io>
- TemplateFlow: <https://www.templateflow.org>
