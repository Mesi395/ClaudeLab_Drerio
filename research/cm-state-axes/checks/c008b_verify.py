# ruff: noqa: E501
"""Supervisor check of C008: C8-4 (one nppa program) and C8-5 (Hr11 ambient pattern) in ventricle cells.

Own disjoint sets without nppa/nppb (border: rest of its set minus shared genes; injury: tnfrsf11b, hspb6, hspb11;
dediff: desma, nkx2.5, tbx20, gata4, hand2, mef2cb, cacybp), own depth quintiles of Hr11.
Usage: uv run python -I research/cm-state-axes/checks/c008b_verify.py data/clean/hu2022_cm_allgenes_sanitized.h5ad results/c005_agent_labels.csv
"""

import sys

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
from scipy.stats import spearmanr

a = ad.read_h5ad(sys.argv[1])
lab = pd.read_csv(sys.argv[2]).set_index("cell").reindex(a.obs_names)
a = a[(lab["agent_chamber"] == "ventricle").to_numpy()].copy()
X = sp.csr_matrix(a.X).astype(np.float32)
lib = a.obs["orig.ident"].astype(str).to_numpy()
tot = np.asarray(X.sum(1)).ravel()
names = list(a.var_names)
b = a.copy()
b.X = X.copy()
sc.pp.normalize_total(b, target_sum=1e4)
sc.pp.log1p(b)
sets = {
    "border": ["myh7", "hk1", "pkma", "pdk2a", "ldha", "mmp14b"],
    "injury": ["tnfrsf11b", "hspb6", "hspb11"],
    "dediff": ["desma", "nkx2.5", "tbx20", "gata4", "hand2", "mef2cb", "cacybp"],
}
for k, gl in sets.items():
    sc.tl.score_genes(b, gl, score_name=k, random_state=0)
pairs = [("border", "injury"), ("border", "dediff"), ("injury", "dediff")]
rows = []
for lb in np.unique(lib):
    m = lib == lb
    if m.sum() < 50:
        continue
    rows.append({p: spearmanr(b.obs[p[0]][m], b.obs[p[1]][m]).statistic for p in pairs})
df = pd.DataFrame(rows)
print("C8-4 disjoint sets without nppa/nppb, ventricle cells: median within-library rho per pair")
print({f"{x}-{y}": round(float(df[(x, y)].median()), 3) for x, y in pairs}, "over", len(df), "libraries")

def share(g):
    return np.asarray(X[:, names.index(g)].todense()).ravel() / tot

ctrl = np.isin(lib, ["H5", "H6", "H7", "H8v"])
h11 = lib == "Hr11"
edges = np.quantile(tot[h11], [0, 0.2, 0.4, 0.6, 0.8, 1])
print("\nC8-5 Hr11 vs control ventricle cells, Hr11 depth quintiles: ratio of mean share; excess molecules per cell")
for g in ["nppa", "myh6", "hbba1.1"]:
    sh = share(g)
    out = []
    for q in range(5):
        lo, hi = edges[q], edges[q + 1] + (1 if q == 4 else 0)
        mh = h11 & (tot >= lo) & (tot < hi)
        mc = ctrl & (tot >= lo) & (tot < hi)
        r = sh[mh].mean() / sh[mc].mean() if sh[mc].mean() > 0 else float("inf")
        exc = (sh[mh].mean() - sh[mc].mean()) * np.median(tot[mh])
        out.append(f"q{q}: {r:.2f}x, +{exc:.2f} mol (n_ctrl {mc.sum()})")
    print(f"  {g}: " + " | ".join(out))
