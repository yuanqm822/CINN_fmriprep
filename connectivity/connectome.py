#!/usr/bin/env python
"""Denoise fMRIPrep outputs and compute functional connectomes with Nilearn.

Implements the post-fMRIPrep workflow of Wang et al. (2024, PLOS Comput Biol,
doi:10.1371/journal.pcbi.1011942): confounds are selected with
nilearn.interfaces.fmriprep.load_confounds_strategy, applied while extracting
atlas time series (NiftiLabelsMasker), and correlated with ConnectivityMeasure.

Run motion_qc.py first to see which runs meet the paper's motion criteria.

Usage (inside `module load nilearn/1.0`):
    python connectome.py <fmriprep_dir> <subject> <output_dir> [options]

Example:
    python connectome.py /scratch5/$USER/proj/derivatives/fmriprep 101 \
        /scratch5/$USER/proj/derivatives/connectome --task rest \
        --strategy simple --global-signal basic --atlas 400
"""
import argparse
import glob
import json
import os
import warnings

import numpy as np
import pandas as pd
import templateflow.api as tf
from nilearn.connectome import ConnectivityMeasure
from nilearn.interfaces.fmriprep import load_confounds, load_confounds_strategy
from nilearn.maskers import NiftiLabelsMasker

SPACE = "MNI152NLin2009cAsym"
warnings.filterwarnings("ignore", category=DeprecationWarning)
warnings.filterwarnings("ignore", category=FutureWarning)

p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
p.add_argument("fmriprep_dir")
p.add_argument("subject", help="subject ID, with or without 'sub-'")
p.add_argument("output_dir")
p.add_argument("--task", default="rest", help="BIDS task label (default: rest)")
p.add_argument("--strategy", default="simple",
               choices=["baseline", "simple", "scrubbing", "compcor"],
               help="denoising strategy (Table 2 of Wang et al. 2024; default: simple). "
                    "'baseline' = high-pass filter only, the paper's reference strategy. "
                    "ICA-AROMA is not offered: fMRIPrep >= 23.1 no longer produces it")
p.add_argument("--global-signal", default=None, choices=["basic", "power2", "derivatives", "full"],
               help="add global signal regression (not used with compcor; default: none)")
p.add_argument("--fd-threshold", type=float, default=0.5,
               help="framewise displacement threshold in mm for scrubbing (default: 0.5)")
p.add_argument("--dvars-threshold", type=float, default=None,
               help="standardised DVARS threshold for scrubbing (default: none, i.e. FD only, "
                    "as in Wang et al. 2024)")
p.add_argument("--atlas", type=int, default=400, choices=[100, 200, 300, 400, 500, 600, 800, 1000],
               help="Schaefer 2018 parcellation size, 7 networks (default: 400)")
p.add_argument("--resolution", type=int, default=2, help="fMRIPrep res- of the MNI output (default: 2)")
a = p.parse_args()

sub = a.subject.removeprefix("sub-")
os.makedirs(a.output_dir, exist_ok=True)

# ---- atlas from the shared TemplateFlow folder (downloaded there if missing) ----
atlas = tf.get(SPACE, resolution=a.resolution, atlas="Schaefer2018",
               desc=f"{a.atlas}Parcels7Networks", suffix="dseg", extension=".nii.gz")
if not atlas:
    raise SystemExit(f"Schaefer {a.atlas} atlas not found in $TEMPLATEFLOW_HOME and could not be downloaded")
labels_tsv = str(atlas).replace("_dseg.nii.gz", "_dseg.tsv")
labels = (pd.read_csv(labels_tsv, sep="\t")["name"].tolist()
          if os.path.exists(labels_tsv) else [f"roi{i + 1}" for i in range(a.atlas)])

# ---- denoising options --------------------------------------------------------
if a.strategy == "baseline":
    if a.global_signal:
        print("Note: --global-signal is ignored for the baseline strategy")
    kw = dict(strategy=["high_pass"])

    def get_confounds(bold):
        return load_confounds(bold, **kw)
else:
    kw = dict(denoise_strategy=a.strategy)
    if a.global_signal and a.strategy != "compcor":
        kw["global_signal"] = a.global_signal
    if a.strategy == "scrubbing":
        kw.update(fd_threshold=a.fd_threshold, std_dvars_threshold=a.dvars_threshold)

    def get_confounds(bold):
        return load_confounds_strategy(bold, **kw)

pattern = (f"{a.fmriprep_dir}/sub-{sub}/**/func/sub-{sub}*_task-{a.task}*"
           f"_space-{SPACE}_res-{a.resolution}_desc-preproc_bold.nii.gz")
bolds = sorted(glob.glob(pattern, recursive=True))
if not bolds:
    raise SystemExit(f"No preprocessed BOLD found matching:\n  {pattern}\n"
                     f"Was fMRIPrep run with --output-spaces {SPACE}:res-{a.resolution}?")

tag = f"atlas-Schaefer{a.atlas}_denoise-{a.strategy}" + ("GSR" if "global_signal" in kw else "")
if a.strategy == "scrubbing":
    tag += f"FD{a.fd_threshold:g}".replace(".", "p")
summary = []
for bold in bolds:
    mask = bold.replace("desc-preproc_bold", "desc-brain_mask")
    confounds, sample_mask = get_confounds(bold)
    n_total = len(confounds)
    n_kept = n_total if sample_mask is None else len(sample_mask)
    if n_kept < 10:
        print(f"SKIPPED {os.path.basename(bold)}: only {n_kept}/{n_total} volumes left after "
              f"scrubbing - exclude this run (see motion_qc.py)")
        summary.append(dict(volumes_kept=n_kept, volumes_total=n_total))
        continue

    masker = NiftiLabelsMasker(atlas, mask_img=mask, standardize="zscore_sample")
    ts = masker.fit_transform(bold, confounds=confounds, sample_mask=sample_mask)
    fc = ConnectivityMeasure(kind="correlation", standardize="zscore_sample").fit_transform([ts])[0]

    names = labels[:ts.shape[1]]
    stem = os.path.join(a.output_dir, os.path.basename(bold).split("_space-")[0] + "_" + tag)
    pd.DataFrame(ts, columns=names).to_csv(f"{stem}_timeseries.tsv", sep="\t", index=False)
    pd.DataFrame(fc, index=names, columns=names).to_csv(f"{stem}_connectome.tsv", sep="\t")

    info = dict(bold=os.path.basename(bold), atlas=os.path.basename(str(atlas)),
                denoise=a.strategy,
                **{k: v for k, v in kw.items() if k not in ("denoise_strategy", "strategy")},
                n_confounds=confounds.shape[1], volumes_total=n_total, volumes_kept=n_kept)
    with open(f"{stem}_desc.json", "w") as fh:
        json.dump(info, fh, indent=2)
    summary.append(info)
    print(f"{os.path.basename(stem)}: {ts.shape[1]} ROIs, {n_kept}/{n_total} volumes kept, "
          f"{confounds.shape[1]} confound regressors")

if any(s["volumes_kept"] < 0.5 * s["volumes_total"] for s in summary):
    print("WARNING: scrubbing removed more than half the volumes in at least one run")
