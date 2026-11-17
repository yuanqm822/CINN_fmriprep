#!/usr/bin/env python
"""Motion quality control on fMRIPrep outputs, following Wang et al. (2024).

A run is flagged for exclusion if ANY of these is true (paper, Materials and methods,
"Participant exclusion based on motion"):
  1. mean framewise displacement (FD) > 0.25 mm
  2. more than 80% of volumes are removed by scrubbing at FD 0.2 mm
  3. less than 1 minute of data remains after that scrubbing

Scrubbing is done exactly as in the paper: Nilearn's 'scrubbing' strategy with
fd_threshold=0.2 and no DVARS threshold (it also removes segments shorter than
5 volumes between censored volumes).

Note: the paper's data had TR 2-2.5 s. With multiband data (e.g. TR 0.8 s),
respiration inflates FD, so more runs may be flagged; the paper itself lists this
as an open question. Look at the numbers, not only the flag.

Usage (inside `module load nilearn/1.0`):
    python motion_qc.py <fmriprep_dir> [--task rest] [--out motion_qc.tsv]

Writes one row per run and prints a summary; runs for all subjects found.

Optionally (--save-confounds DIR) it also saves, for every run, the confound
regressors chosen with Nilearn's load_confounds, e.g. for your own GLM:
    python motion_qc.py <fmriprep_dir> --save-confounds <dir> \
        --motion full --wm-csf full --global-signal basic
Levels for --motion / --wm-csf / --global-signal:
    basic        original parameters                    (motion: 6, wm_csf: 2)
    power2       + quadratic terms                      (12, 4)
    derivatives  + 1st temporal derivatives             (12, 4)
    full         + derivatives + quadratic terms of     (24, 8)
                   both parameters and derivatives
High-pass (cosine) regressors are included unless --no-high-pass.
"""
import argparse
import glob
import json
import os
import re
import warnings

import numpy as np
import pandas as pd
from nilearn.interfaces.fmriprep import load_confounds, load_confounds_strategy

SPACE = "MNI152NLin2009cAsym"
warnings.filterwarnings("ignore", category=DeprecationWarning)
warnings.filterwarnings("ignore", category=FutureWarning)

p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
p.add_argument("fmriprep_dir")
p.add_argument("--task", default="rest", help="BIDS task label (default: rest)")
p.add_argument("--subjects", nargs="*", help="only these subject IDs (default: all)")
p.add_argument("--out", default=None, help="output TSV (default: <fmriprep_dir>/motion_qc_task-<task>.tsv)")
p.add_argument("--mean-fd", type=float, default=0.25, help="max mean FD in mm (default: 0.25)")
p.add_argument("--scrub-fd", type=float, default=0.2, help="FD threshold for scrubbing check (default: 0.2)")
p.add_argument("--max-removed", type=float, default=0.8, help="max fraction removed (default: 0.8)")
p.add_argument("--min-minutes", type=float, default=1.0, help="min data left in minutes (default: 1)")
LEVELS = ["basic", "power2", "derivatives", "full"]
g = p.add_argument_group("saving confound regressors (optional)")
g.add_argument("--save-confounds", metavar="DIR", default=None,
               help="also save the selected confounds of every run to DIR")
g.add_argument("--motion", choices=LEVELS, default="full", help="motion parameters (default: full = 24)")
g.add_argument("--wm-csf", choices=LEVELS, default="basic", help="WM/CSF signals (default: basic = 2)")
g.add_argument("--global-signal", choices=LEVELS, default=None, help="global signal (default: not included)")
g.add_argument("--no-high-pass", action="store_true", help="leave out the cosine high-pass regressors")
a = p.parse_args()

if a.save_confounds:
    os.makedirs(a.save_confounds, exist_ok=True)
    conf_strategy = (["high_pass"] if not a.no_high_pass else []) + ["motion", "wm_csf"]
    conf_kw = dict(motion=a.motion, wm_csf=a.wm_csf)
    if a.global_signal:
        conf_strategy.append("global_signal"); conf_kw["global_signal"] = a.global_signal
    conf_desc = f"motion{a.motion}wmcsf{a.wm_csf}" + (f"gs{a.global_signal}" if a.global_signal else "")

bolds = sorted(glob.glob(f"{a.fmriprep_dir}/sub-*/**/func/*_task-{a.task}_*"
                         f"space-{SPACE}*desc-preproc_bold.nii.gz", recursive=True))
# one preprocessed image per run is enough (several resolutions may exist)
runs = {}
for b in bolds:
    key = os.path.basename(b).split("_space-")[0]
    runs.setdefault(key, b)
if a.subjects:
    keep = {s.removeprefix("sub-") for s in a.subjects}
    runs = {k: v for k, v in runs.items() if re.match(r"sub-([^_]+)", k)[1] in keep}
if not runs:
    raise SystemExit(f"No task-{a.task} preprocessed BOLD found in {a.fmriprep_dir}")

rows = []
for key, bold in runs.items():
    conf_tsv = re.sub(r"_space-.*$", "_desc-confounds_timeseries.tsv", bold)
    fd = pd.read_csv(conf_tsv, sep="\t")["framewise_displacement"]
    meta = bold.replace(".nii.gz", ".json")
    tr = json.load(open(meta)).get("RepetitionTime") if os.path.exists(meta) else None

    _, sample_mask = load_confounds_strategy(bold, denoise_strategy="scrubbing",
                                             fd_threshold=a.scrub_fd, std_dvars_threshold=None)
    n = len(fd)
    kept = n if sample_mask is None else len(sample_mask)
    removed = 1 - kept / n
    minutes_left = kept * tr / 60 if tr else np.nan

    reasons = []
    if fd.mean() > a.mean_fd:
        reasons.append(f"mean FD {fd.mean():.3f} > {a.mean_fd}")
    if removed > a.max_removed:
        reasons.append(f"{removed:.0%} removed at FD {a.scrub_fd}")
    if tr and minutes_left < a.min_minutes:
        reasons.append(f"{minutes_left:.2f} min left")

    if a.save_confounds:
        conf, _ = load_confounds(bold, strategy=conf_strategy, **conf_kw)
        stem = os.path.join(a.save_confounds, f"{key}_desc-{conf_desc}_confounds")
        conf.to_csv(f"{stem}.tsv", sep="\t", index=False)
        with open(f"{stem}.json", "w") as fh:
            json.dump(dict(source=os.path.basename(conf_tsv), strategy=conf_strategy, **conf_kw,
                           n_regressors=conf.shape[1], columns=list(conf.columns)), fh, indent=2)

    rows.append(dict(run=key, TR=tr, volumes=n, mean_FD=round(fd.mean(), 4),
                     max_FD=round(fd.max(), 4),
                     pct_FD_over_thresh=round(100 * (fd > a.scrub_fd).mean(), 1),
                     volumes_kept_scrubbing=kept, pct_removed_scrubbing=round(100 * removed, 1),
                     minutes_left=round(minutes_left, 2), exclude=bool(reasons),
                     reason="; ".join(reasons)))

df = pd.DataFrame(rows)
out = a.out or os.path.join(a.fmriprep_dir, f"motion_qc_task-{a.task}.tsv")
df.to_csv(out, sep="\t", index=False)

with pd.option_context("display.width", 200, "display.max_columns", 20):
    print(df[["run", "mean_FD", "pct_removed_scrubbing", "minutes_left", "exclude", "reason"]]
          .to_string(index=False))
print(f"\n{df.exclude.sum()} of {len(df)} runs flagged for exclusion. Full table: {out}")
if a.save_confounds:
    print(f"Saved {conf_desc} confounds ({', '.join(conf_strategy)}) for {len(df)} runs in {a.save_confounds}")
