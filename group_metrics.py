#!/usr/bin/env python
"""Group-level denoising benchmark metrics, following Wang et al. (2024).

For every denoising strategy found in the connectome folder, computes:

  loss of tDOF   number of confound regressors + number of censored volumes
                 (lower = more data preserved)
  QC-FC          partial correlation, across subjects, between mean framewise
                 displacement and each connectivity edge, controlling for age and
                 sex; reported as % of edges with p < 0.05 (uncorrected and FDR)
                 and median |r| (lower = less residual motion)
  DM-FC          Spearman correlation between the Euclidean distance of each pair
                 of regions and their QC-FC (closer to 0 = less distance-dependent
                 motion artefact)
  modularity     average Louvain modularity Q of each subject's connectome
                 (positive edges only; higher = clearer network structure)

Runs flagged by motion_qc.py are left out. If a subject has several runs, their
connectomes (Fisher z) and mean FD are averaged. These metrics need a group:
results with fewer than ~20 subjects are not meaningful.

Usage (inside `module load nilearn/1.0`):
    python group_metrics.py <connectome_dir> <motion_qc.tsv> <participants.tsv> \
        [--age-col age] [--sex-col sex] [--out <dir>]
"""
import argparse
import glob
import json
import os
import re
import warnings
from collections import defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
import templateflow.api as tf
import nibabel as nib
from scipy import stats

warnings.filterwarnings("ignore", category=DeprecationWarning)
warnings.filterwarnings("ignore", category=FutureWarning)
SPACE = "MNI152NLin2009cAsym"

p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
p.add_argument("connectome_dir", help="output folder of run_connectome.sh")
p.add_argument("motion_qc", help="TSV written by motion_qc.py")
p.add_argument("participants", help="BIDS participants.tsv (needs age and sex)")
p.add_argument("--age-col", default=None, help="age column (default: first column starting with 'age')")
p.add_argument("--sex-col", default="sex", help="sex column (default: sex)")
p.add_argument("--resolution", type=int, default=2, help="atlas resolution used (default: 2)")
p.add_argument("--out", default=None, help="output folder (default: <connectome_dir>/group)")
p.add_argument("--min-subjects", type=int, default=20, help="warn below this many subjects (default: 20)")
a = p.parse_args()
out = a.out or os.path.join(a.connectome_dir, "group")
os.makedirs(out, exist_ok=True)

# ---- inputs -------------------------------------------------------------------
qc = pd.read_csv(a.motion_qc, sep="\t").set_index("run")
part = pd.read_csv(a.participants, sep="\t", dtype={"participant_id": str})
part["participant_id"] = part["participant_id"].str.removeprefix("sub-")
age_col = a.age_col or next((c for c in part.columns if c.lower().startswith("age")), None)
covars = []
if age_col in part:
    part["_age"] = pd.to_numeric(part[age_col], errors="coerce"); covars.append("_age")
if a.sex_col in part:
    part["_sex"] = pd.factorize(part[a.sex_col].astype(str).str.upper().str[0])[0].astype(float)
    covars.append("_sex")
if not covars:
    print("WARNING: no age/sex columns found - QC-FC will be a plain correlation")
part = part.set_index("participant_id")

files = sorted(glob.glob(os.path.join(a.connectome_dir, "*_connectome.tsv")))
pat = re.compile(r"^(?P<run>sub-(?P<sub>[^_]+)_.*?)_atlas-(?P<atlas>Schaefer\d+)_denoise-(?P<strat>.+)_connectome\.tsv$")
groups = defaultdict(list)
for f in files:
    m = pat.match(os.path.basename(f))
    if m:
        groups[(m["atlas"], m["strat"])].append((m["run"], m["sub"], f))
if not groups:
    raise SystemExit(f"No *_connectome.tsv files found in {a.connectome_dir}")

_coords = {}


def region_coords(atlas, names):
    """MNI centroids of the atlas regions, in the order of the connectome columns."""
    if atlas not in _coords:
        n = int(atlas.removeprefix("Schaefer"))
        img = tf.get(SPACE, resolution=a.resolution, atlas="Schaefer2018",
                     desc=f"{n}Parcels7Networks", suffix="dseg", extension=".nii.gz")
        lut = pd.read_csv(str(img).replace("_dseg.nii.gz", "_dseg.tsv"), sep="\t")
        nii = nib.load(str(img))
        lab = np.asarray(nii.dataobj).astype(int)
        by_name = {}
        for idx, name in zip(lut["index"], lut["name"]):
            vox = np.argwhere(lab == idx)
            if len(vox):
                by_name[name] = nib.affines.apply_affine(nii.affine, vox.mean(0))
        _coords[atlas] = by_name
    return np.array([_coords[atlas][nm] for nm in names])


def fdr(pvals, q=0.05):
    """Benjamini-Hochberg: boolean mask of discoveries."""
    pv = np.asarray(pvals); n = pv.size
    order = np.argsort(pv)
    passed = pv[order] <= q * np.arange(1, n + 1) / n
    mask = np.zeros(n, bool)
    if passed.any():
        mask[order[: np.max(np.nonzero(passed)) + 1]] = True
    return mask


def residualise(y, X):
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    return y - X @ beta


# ---- metrics per strategy ---------------------------------------------------------
rows = []
for (atlas, strat), items in sorted(groups.items()):
    subj_z, subj_fd, subj_dof, names = defaultdict(list), defaultdict(list), defaultdict(list), None
    n_excluded = 0
    for run, sub, f in items:
        if run not in qc.index:
            print(f"  {run}: not in motion QC table - skipped"); continue
        if bool(qc.loc[run, "exclude"]):
            n_excluded += 1; continue
        fc = pd.read_csv(f, sep="\t", index_col=0)
        names = list(fc.columns)
        z = np.arctanh(np.clip(fc.values, -0.999999, 0.999999))
        np.fill_diagonal(z, 0)
        d = json.load(open(f.replace("_connectome.tsv", "_desc.json")))
        subj_z[sub].append(z)
        subj_fd[sub].append(qc.loc[run, "mean_FD"])
        subj_dof[sub].append(d["n_confounds"] + d["volumes_total"] - d["volumes_kept"])

    subs = [s for s in subj_z if s in part.index and not part.loc[s, covars].isna().any()] \
        if covars else list(subj_z)
    dropped = set(subj_z) - set(subs)
    if dropped:
        print(f"  {strat}: no age/sex for {sorted(dropped)} - left out")
    n = len(subs)
    if n < 4:
        print(f"{atlas} {strat}: only {n} usable subjects - skipped"); continue

    Z = np.stack([np.mean(subj_z[s], axis=0) for s in subs])           # subjects x R x R
    fd = np.array([np.mean(subj_fd[s]) for s in subs])
    dof = np.array([np.mean(subj_dof[s]) for s in subs])
    iu = np.triu_indices(Z.shape[1], k=1)
    edges = Z[:, iu[0], iu[1]]                                          # subjects x edges

    # QC-FC: partial correlation of FD with each edge, controlling covariates
    X = np.column_stack([np.ones(n)] + [part.loc[subs, c].values for c in covars])
    fd_r = residualise(fd, X)
    ed_r = residualise(edges, X)
    ed_r -= ed_r.mean(0); fd_c = fd_r - fd_r.mean()
    r = (fd_c @ ed_r) / (np.linalg.norm(fd_c) * np.linalg.norm(ed_r, axis=0) + 1e-12)
    df = n - 2 - len(covars)
    t = r * np.sqrt(df / np.clip(1 - r ** 2, 1e-12, None))
    pv = 2 * stats.t.sf(np.abs(t), df)

    # DM-FC: distance dependence of QC-FC
    xyz = region_coords(atlas, names)
    dist = np.linalg.norm(xyz[iu[0]] - xyz[iu[1]], axis=1)
    dmfc = stats.spearmanr(dist, r).statistic

    # Modularity (Louvain on positive weights), per subject
    q = []
    for zmat in Z:
        w = np.where(zmat > 0, zmat, 0)
        g = nx.from_numpy_array(w)
        comms = nx.community.louvain_communities(g, weight="weight", seed=0)
        q.append(nx.community.modularity(g, comms, weight="weight"))

    rows.append(dict(atlas=atlas, strategy=strat, n_subjects=n, n_runs_excluded=n_excluded,
                     mean_FD=round(fd.mean(), 4), loss_tDOF_mean=round(dof.mean(), 1),
                     qcfc_pct_sig_uncorrected=round(100 * (pv < 0.05).mean(), 2),
                     qcfc_pct_sig_fdr=round(100 * fdr(pv).mean(), 2),
                     qcfc_median_abs_r=round(np.median(np.abs(r)), 4),
                     dmfc_spearman_r=round(dmfc, 4),
                     modularity_mean=round(np.mean(q), 4)))
    np.savetxt(os.path.join(out, f"qcfc_{atlas}_{strat}.tsv"),
               np.column_stack([iu[0], iu[1], dist, r, pv]), delimiter="\t",
               header="i\tj\tdistance_mm\tqcfc_r\tp", comments="", fmt=["%d", "%d", "%.2f", "%.5f", "%.3g"])
    print(f"{atlas} {strat}: {n} subjects done")

res = pd.DataFrame(rows)
if res.empty:
    raise SystemExit("No strategy had enough subjects.")
if res.n_subjects.min() < a.min_subjects:
    print(f"\nWARNING: some strategies have fewer than {a.min_subjects} subjects - "
          "treat the metrics as unreliable")

# Ranks within each atlas (1 = best), as in the paper's Fig 12
better_low = ["loss_tDOF_mean", "qcfc_pct_sig_uncorrected", "qcfc_median_abs_r"]
for c in better_low:
    res[f"rank_{c}"] = res.groupby("atlas")[c].rank(method="min")
res["rank_dmfc"] = res.groupby("atlas")["dmfc_spearman_r"].transform(lambda s: s.abs().rank(method="min"))
res["rank_modularity"] = res.groupby("atlas")["modularity_mean"].rank(ascending=False, method="min")
res.to_csv(os.path.join(out, "group_metrics.tsv"), sep="\t", index=False)

# Figure: one row of bar charts per atlas
metrics = [("loss_tDOF_mean", "Loss of tDOF"), ("qcfc_pct_sig_uncorrected", "QC-FC: % sig. edges"),
           ("qcfc_median_abs_r", "QC-FC: median |r|"), ("dmfc_spearman_r", "DM-FC (r)"),
           ("modularity_mean", "Modularity Q")]
atlases = res.atlas.unique()
fig, axes = plt.subplots(len(atlases), len(metrics), figsize=(4 * len(metrics), 3.2 * len(atlases)),
                         squeeze=False)
for i, at in enumerate(atlases):
    sub = res[res.atlas == at]
    for j, (col, title) in enumerate(metrics):
        ax = axes[i, j]
        ax.barh(sub.strategy, sub[col], color="#4C72B0")
        ax.set_title(f"{title}\n{at}" if i == 0 else at, fontsize=9)
        ax.tick_params(labelsize=8)
        if j: ax.set_yticklabels([])
fig.tight_layout()
fig.savefig(os.path.join(out, "group_metrics.png"), dpi=150)

with pd.option_context("display.width", 200, "display.max_columns", 20):
    print("\n" + res[["atlas", "strategy", "n_subjects", "loss_tDOF_mean", "qcfc_pct_sig_uncorrected",
                      "qcfc_median_abs_r", "dmfc_spearman_r", "modularity_mean"]].to_string(index=False))
print(f"\nSaved {out}/group_metrics.tsv, group_metrics.png and per-edge QC-FC tables")
