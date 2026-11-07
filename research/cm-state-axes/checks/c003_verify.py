# ruff: noqa: E501
"""Supervisor check of C003: do axis scores track chamber (atrial share) rather than injury?

Uses only the sanitised file plus the researcher's per-heart table for the score values.
Usage: uv run python -I research/cm-state-axes/checks/c003_verify.py data/clean/hu2022_cm_allgenes_sanitized.h5ad results/c003_heart_scores.csv
"""

import sys

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
from scipy.stats import spearmanr

a = ad.read_h5ad(sys.argv[1])
X = sp.csr_matrix(a.X).astype(np.float32)
lib = a.obs["orig.ident"].astype(str).to_numpy()
av = a.obs["AV"].astype(str).to_numpy()
tot = np.asarray(X.sum(1)).ravel()
g = {s: i for i, s in enumerate(a.var_names)}
heart = pd.Series(lib).str.replace(r"^(Hr2|H8|Hr6|Hr7)[abav]$", r"\1", regex=True).to_numpy()


def share(sym):
    return np.asarray(X[:, g[sym]].todense()).ravel() / tot


# atrial index per cell: share of myh6 (atrial myosin) among myh6 + myh7l (likely ventricular paralog)
m6, m7l = share("myh6"), share("myh7l")
atrial = (m6 > m7l) & (m6 > 0)

b = a.copy()
b.X = X.copy()
sc.pp.normalize_total(b, target_sum=1e4)
sc.pp.log1p(b)
sets = {
    "border_zone_like": [
        "nppa",
        "mustn1b",
        "myh7",
        "hk1",
        "pkma",
        "pdk2a",
        "ldha",
        "nppb",
        "mmp14b",
        "ankrd1a",
    ],
    "injury_responsive": ["nppa", "nppb", "tnfrsf11b", "hspb6", "hspb11", "ankrd1a"],
}
for k, gl in sets.items():
    sc.tl.score_genes(b, gl, score_name=k, random_state=0)

print("Chamber-only libraries (AV from metadata, constant per library): median score_genes")
df = pd.DataFrame(
    {"lib": lib, "AV": av, **{k: b.obs[k].to_numpy() for k in sets}, "atrial_call": atrial}
)
print(
    df[df.AV != "Wholeheart"]
    .groupby(["lib", "AV"])[list(sets) + ["atrial_call"]]
    .median(numeric_only=True)
    .round(3)
    .to_string()
)
print(
    "\nfraction of cells called atrial by myh6>myh7l, by AV:",
    df.groupby("AV").atrial_call.mean().round(3).to_dict(),
)

hs = pd.read_csv(sys.argv[2])
af = pd.Series(atrial).groupby(heart).mean()
print(
    "\nAcross 18 hearts: Spearman of per-heart median thinned score (researcher table) with atrial-call fraction"
)
for k in sets:
    t = hs[hs.state == k].set_index("heart")
    r = spearmanr(t["median_thinned_score"], af.reindex(t.index)).statistic
    print(f"  {k}: rho {r:.2f}")
print("\natrial-call fraction per heart:", af.round(3).sort_values().to_dict())

# within ventricle-called cells only: injured vs control hearts, border axis per-heart median
vent = ~atrial
per = pd.DataFrame(
    {
        "heart": heart,
        "time": a.obs["time"].astype(str).to_numpy(),
        "s": b.obs["border_zone_like"].to_numpy(),
    }
)[vent]
ph = per.groupby(["heart", "time"]).s.median().reset_index()
print("\nborder axis, ventricle-called cells only, per-heart median (unthinned score_genes):")
print(ph.sort_values("s").to_string(index=False))
