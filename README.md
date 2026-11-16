# Running fMRIPrep and MRIQC on RACC (CINN)

fMRIPrep and MRIQC are installed on RACC (University of Reading) as CINN (Centre for Integrative Neuroscience and Neurodynamics) modules (Neurodesk containers via
Transparent Singularity). You do **not** need to clone a package, build a conda
environment, or download containers yourself — load the module and submit a Slurm job.

> This replaces the older conda/Python-package approach. If you followed that guide,
> remove any `CINN_fmriprep` / fmriprep path lines it added to your `~/.bashrc`.

| Tool     | Module            |
|----------|-------------------|
| fMRIPrep | `fmriprep/25.2.7` |
| MRIQC    | `mriqc/26.1.0`    |
| Nilearn (connectivity) | `nilearn/1.0` |

---

## 1. One-time setup

Add this line to the end of your `~/.bashrc`:

```bash
module use /software/modules/CINN
```

Then reload and check:

```bash
source ~/.bashrc
module avail fmriprep mriqc
```

The modules already set:

- the shared FreeSurfer licence (`FS_LICENSE`) — no need for `--fs-license-file`
- container bind paths for `/software`, `/home`, `/storage`, `/scratch2`, `/scratch3`, `/scratch4`, `/scratch5`, `/scratch6`.

## 2. Prepare your data

- A valid BIDS dataset (e.g. made with bidscoin) <https://bidscoin.readthedocs.io/en/stable/preparation.html>
- A `subjects.txt` file with one subject ID per line, **without** the `sub-` prefix:

```
101
102
103
```

## 3. Edit and submit the job script

Clone this repository into your project folder on RACC and move into it:

```bash
cd /scratch5/$USER/<project>
git clone https://github.com/yuanqm822/CINN_fmriprep.git
cd CINN_fmriprep
```

Put your `subjects.txt` in this folder (or copy `examples/subjects.txt` and edit it).

Open `run_fmriprep.sh` and change the four lines under `EDIT THESE`:
your BIDS folder, output folder, work folder and subject list. Then set
`#SBATCH --array=1-N` where `N` is the number of lines in `subjects.txt`:

```bash
wc -l < subjects.txt
```

Submit:

```bash
mkdir -p logs
sbatch run_fmriprep.sh
```

Each subject runs as a separate array task, in parallel, as resources become available.

### Adjusting fMRIPrep options

The `OPTIONS` block in `run_fmriprep.sh` controls output spaces and any extra
fMRIPrep flags. For example, to get FreeSurfer surface outputs and force BBR
coregistration:

```bash
OUTPUT_SPACES=(fsnative fsaverage)
EXTRA_ARGS=(--force bbr)
```

Other common settings:

| Goal | Setting |
|------|---------|
| Volume + surface outputs | `OUTPUT_SPACES=(MNI152NLin2009cAsym:res-2 fsaverage5)` |
| Skip FreeSurfer (faster, no surfaces) | `EXTRA_ARGS=(--fs-no-reconall)` |
| CIFTI grayordinates | `EXTRA_ARGS=(--cifti-output 91k)` |
| Drop dummy scans | `EXTRA_ARGS=(--dummy-scans 4)` |
| Skip distortion correction | `EXTRA_ARGS=(--ignore fieldmaps)` |
| Several at once | `EXTRA_ARGS=(--force bbr --cifti-output 91k)` |

For a one-off change you can also pass flags at submission without editing the
script. They are added to the end of the fMRIPrep command:

```bash
sbatch run_fmriprep.sh --force bbr --dummy-scans 4
```

Surface spaces (`fsnative`, `fsaverage*`, `fsLR`) need FreeSurfer reconstruction,
so don't combine them with `--fs-no-reconall`. The full list of flags is at
<https://fmriprep.org/en/stable/usage.html>. `run_mriqc.sh` has the same
`EXTRA_ARGS` option.

MRIQC works the same way with `run_mriqc.sh`. When all participant jobs have
finished, run the group report once:

```bash
module load mriqc/26.1.0
mriqc /storage/<project>/bids /storage/<project>/derivatives/mriqc group --no-sub
```

## 4. Monitor jobs

```bash
squeue -u $USER                     # running / pending jobs
sacct -j <jobid> --format=JobID,State,Elapsed,MaxRSS
tail -f logs/fmriprep_<jobid>_1.out # follow one subject's log
scancel <jobid>                     # cancel all tasks in an array
scancel <jobid>_1                   # cancel one task
```

Outputs appear in your output folder; the HTML report for each subject
(`sub-XX.html`) is the first thing to check.

## 5. Resources and timing

The example requests 8 CPUs, 32 GB RAM and 12 hours per subject. As a rough guide,
one T1w plus one ~200-volume BOLD run takes around 45 minutes; more runs, higher
resolution or `--cifti-output` / FreeSurfer reconstruction take longer. Request
more time rather than less — a job killed at the time limit has to restart.

## 6. Tips and troubleshooting

- **Work directory:** keep `-w` on a scratch volume (e.g. `/scratch5`), not `/home` —
  it gets large. Delete it once a subject has finished successfully.
- **Log files:** jobs will fail straight away, without writing any output, if the
  `logs/` folder doesn't exist, so run `mkdir -p logs` before `sbatch`. Each subject
  gets `logs/fmriprep_<jobid>_<task>.out` (progress) and `.err` (errors). Check the
  `.err` file first when a job fails.
- **"File not found" for a path that exists:** the container can only see bound
  folders. If your data live outside `/software`, `/home`, `/storage` or
  `/scratch2`–`/scratch6`, add the path before running:
  `export APPTAINER_BINDPATH="$APPTAINER_BINDPATH,/other/path"`
- **MRIQC hangs or fails at the end:** always use `--no-sub`; compute nodes cannot
  upload metrics to the MRIQC web API.
- **Re-running a subject:** reusing the same work directory lets fMRIPrep pick up
  where it stopped.
- **Help text:** `module help fmriprep/25.2.7`

## 7. Next step: functional connectivity with Nilearn

After fMRIPrep, `run_connectome.sh` denoises the preprocessed resting-state data and
computes a functional connectome for each run, following the workflow of
Wang et al. (2024). It uses `nilearn/1.0`, a CINN Python environment (built with
pixi) that contains Nilearn and TemplateFlow, so there is no conda environment to
set up. The Schaefer 2018 atlases are pre-downloaded to the shared TemplateFlow
folder, so jobs work on compute nodes without internet.

For each run it:

1. selects confounds from fMRIPrep's `desc-confounds_timeseries.tsv` with
   `load_confounds_strategy`
2. extracts denoised Schaefer-atlas time series (`NiftiLabelsMasker`)
3. computes a correlation matrix (`ConnectivityMeasure`)

It needs fMRIPrep outputs in `MNI152NLin2009cAsym:res-2`, which the default
`OUTPUT_SPACES` in `run_fmriprep.sh` already produces.

### Step 1: motion quality control

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

### Step 2: connectomes

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

The environment is defined in [`envs/nilearn/pixi.toml`](envs/nilearn/pixi.toml) and
built with CINN's own pixi under `/software/CINN/pixi`:

```bash
export PIXI_CACHE_DIR=/scratch5/$USER/pixi_cache   # keep the cache out of /home
mkdir -p /software/CINN/envs/nilearn-1.0
cp envs/nilearn/pixi.toml /software/CINN/envs/nilearn-1.0/
cd /software/CINN/envs/nilearn-1.0 && /software/CINN/pixi/bin/pixi install
chmod -R a+rX /software/CINN/envs/nilearn-1.0
```

The module `/software/modules/CINN/nilearn/1.0` puts
`/software/CINN/envs/nilearn-1.0/.pixi/envs/default/bin` on `PATH` and sets
`TEMPLATEFLOW_HOME=/software/CINN/templateflow`. New versions go in a new folder
(e.g. `nilearn-1.1`) so existing analyses keep the same package versions.
</details>

## Documentation

- fMRIPrep: <https://fmriprep.org/en/stable/>
- MRIQC: <https://mriqc.readthedocs.io>
- Nilearn: <https://nilearn.github.io>
- Wang H-T et al. (2024) Continuous evaluation of denoising strategies in resting-state fMRI connectivity using fMRIPrep and Nilearn. *PLOS Comput Biol* 20(3): e1011942. <https://doi.org/10.1371/journal.pcbi.1011942>
- Transparent Singularity: <https://neurodesk.org/developers/architecture/transparent-singularity/>

## Acknowledgements

Based on Nick Hedger's original guide and `CINN_fmriprep` package for running
fMRIPrep on RACC2.

Questions or problems: contact Dr Simon Yuan (simon.yuan@reading.ac.uk).
