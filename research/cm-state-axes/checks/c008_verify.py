# ruff: noqa: E501
"""Supervisor check after the user's question: what do the axes track WITHIN ventricle-called cells?

Per library (ventricle-called cells only, agent_chamber from results/c005_agent_labels.csv), Spearman of each
axis score (scanpy score_genes; UCell for border) with log10 depth, mito share incl. NC-002333 features, and
the nppa share alone. Also per-heart median border score vs heart (variance share).
Usage: uv run python -I research/cm-state-axes/checks/c008_verify.py data/clean/hu2022_cm_allgenes_sanitized.h5ad results/c005_agent_labels.csv
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
v = (lab["agent_chamber"] == "ventricle").to_numpy()
a = a[v].copy()
X = sp.csr_matrix(a.X).astype(np.float32)
lib = a.obs["orig.ident"].astype(str).to_numpy()
tot = np.asarray(X.sum(1)).ravel()
names = np.array(a.var_names)
mito = np.array([n.startswith("mt-") or n.startswith("NC-002333") for n in names])
mito_share = np.asarray(X[:, mito].sum(1)).ravel() / tot
nppa = np.asarray(X[:, list(names).index("nppa")].todense()).ravel() / tot
print("ventricle-called cells:", a.n_obs)

b = a.copy()
b.X = X.copy()
sc.pp.normalize_total(b, target_sum=1e4)
sc.pp.log1p(b)
SETS = {
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
    "dedifferentiated": [
        "mustn1b",
        "nppa",
        "nppb",
        "desma",
        "nkx2.5",
        "tbx20",
        "gata4",
        "hand2",
        "mef2cb",
        "cacybp",
    ],
    "proliferating": ["pcna", "cdk1", "ccng1", "ccng2", "foxm1", "cenpf"],
}
for k, gl in SETS.items():
    sc.tl.score_genes(b, gl, score_name=k, random_state=0)


def ucell(cols, max_rank=1500):
    out = np.empty(X.shape[0])
    n = len(cols)
    for c in range(X.shape[0]):
        s, e = X.indptr[c], X.indptr[c + 1]
        r = pd.Series(-X.data[s:e]).rank(method="average").to_numpy()
        ro = dict(zip(X.indices[s:e], r, strict=True))
        rk = np.array([min(ro.get(g, max_rank + 1), max_rank + 1) for g in cols])
        out[c] = 1 - (rk.sum() - n * (n + 1) / 2) / (n * max_rank)
    return out


uc = ucell([list(names).index(g) for g in SETS["border_zone_like"]])
rows = []
for lb in np.unique(lib):
    m = lib == lb
    if m.sum() < 50:
        continue
    r = {"lib": lb, "n": int(m.sum())}
    for k in SETS:
        s = b.obs[k].to_numpy()[m]
        r[f"{k}~depth"] = spearmanr(s, np.log10(tot[m])).statistic
        r[f"{k}~mito"] = spearmanr(s, mito_share[m]).statistic
        r[f"{k}~nppa"] = spearmanr(s, nppa[m]).statistic
    r["border_ucell~depth"] = spearmanr(uc[m], np.log10(tot[m])).statistic
    rows.append(r)
df = pd.DataFrame(rows)
print(f"libraries with >=50 ventricle cells: {len(df)}")
print("\nMedian over libraries of within-library Spearman rho (ventricle cells only):")
print(df.drop(columns=["lib", "n"]).median().round(3).to_string())
print(
    "\nF1 within ventricle: libraries where |rho(UCell, depth)| > |rho(score_genes, depth)| for border:",
    int((df["border_ucell~depth"].abs() > df["border_zone_like~depth"].abs()).sum()),
    "of",
    len(df),
)
heart = pd.Series(lib).str.replace(r"^(Hr2|H8|Hr6|Hr7)[abav]$", r"\1", regex=True).to_numpy()
s = pd.Series(b.obs["border_zone_like"].to_numpy())
eta2 = s.groupby(heart).transform("mean").var() / s.var()
print(
    "\nborder score: share of cell-level variance explained by heart (eta^2):",
    round(float(eta2), 3),
)
