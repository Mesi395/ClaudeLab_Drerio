# ruff: noqa: E501
"""Supervisor check of C007 claim C7-3: cell-cycle signal only in Hr9 and Hr8.

Independent, score-free measure: per heart, share of ventricle-called cells (agent_chamber from
results/c005_agent_labels.csv) with >= 2 distinct genes of the Li SD2 proliferating set detected,
and the same after thinning every cell to 500 UMIs (seed 7).
Usage: uv run python -I research/cm-state-axes/checks/c007_verify.py data/clean/hu2022_cm_allgenes_sanitized.h5ad proposals/loop/li2025_sd2_sets.tsv results/c005_agent_labels.csv results/c001_library_qc.csv
"""

import sys

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

a = ad.read_h5ad(sys.argv[1])
sets = pd.read_csv(sys.argv[2], sep="\t")
axcol = next(c for c in sets.columns if c.lower() in ("axis", "state"))
gcol = next(c for c in sets.columns if c.lower() in ("gene", "symbol"))
prolif = [g for g in sets.loc[sets[axcol].str.contains("prolif"), gcol] if g in a.var_names]
print("proliferating genes used:", len(prolif), prolif)
lab = pd.read_csv(sys.argv[3]).set_index("cell").reindex(a.obs_names)
vent = (lab["agent_chamber"] == "ventricle").to_numpy()
X = sp.csr_matrix(a.X)[:, [list(a.var_names).index(g) for g in prolif]]
tot = np.asarray(sp.csr_matrix(a.X).sum(1)).ravel()
heart = lab["heart"].to_numpy()
time = a.obs["time"].astype(str).to_numpy()

rng = np.random.default_rng(7)
Xd = X.toarray()
p = np.minimum(1.0, 500 / tot)[:, None]
Xt = rng.binomial(Xd.astype(int), p)
res = pd.DataFrame(
    {
        "heart": heart,
        "time": time,
        "vent": vent,
        "cc2": (Xd > 0).sum(1) >= 2,
        "cc2_thin": (Xt > 0).sum(1) >= 2,
    }
)[vent]
out = (
    res.groupby(["heart", "time"])
    .agg(n=("cc2", "size"), frac_cc2=("cc2", "mean"), frac_cc2_thin=("cc2_thin", "mean"))
    .reset_index()
)
qc = pd.read_csv(sys.argv[4])
qc["heart"] = qc["lib"].str.replace(r"^(Hr2|H8|Hr6|Hr7)[abav]$", r"\1", regex=True)
q = qc.groupby("heart").agg(
    median_depth=("median_total_counts", "median"), mt=("median_mt_frac", "median")
)
out = out.join(q, on="heart").sort_values("frac_cc2", ascending=False)
print(out.round(4).to_string(index=False))
