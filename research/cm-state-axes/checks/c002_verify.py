# ruff: noqa: E501
"""Supervisor checks of C002 claims P1 (border axis vs matched null), P4 (rank vs depth) and P7 (Hr11 nppa).

Independent of the researcher's code: scanpy score_genes, own UCell and AUCell-style rank scores,
own random gene sets, own depth bins.
Usage: uv run python -I research/cm-state-axes/checks/c002_verify.py data/clean/hu2022_cm_allgenes_sanitized.h5ad
"""

import sys

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
from scipy.stats import spearmanr

rng = np.random.default_rng(1)
a = ad.read_h5ad(sys.argv[1])
X = sp.csr_matrix(a.X).astype(np.float32)
lib = a.obs["orig.ident"].astype(str).to_numpy()
tot = np.asarray(X.sum(1)).ravel()
logtot = np.log10(tot)
genes = np.array(a.var_names)
gi = {g: i for i, g in enumerate(genes)}
BORDER = ["nppa", "mustn1b", "myh7", "hk1", "pkma", "pdk2a", "ldha", "nppb", "mmp14b", "ankrd1a"]
idx = [gi[g] for g in BORDER]

b = a.copy()
b.X = X.copy()
sc.pp.normalize_total(b, target_sum=1e4)
sc.pp.log1p(b)


def score_genes(gl, seed=0):
    sc.tl.score_genes(b, gl, score_name="_s", random_state=seed)
    return b.obs["_s"].to_numpy()


def rank_scores(cols, max_rank=1500, k_frac=0.05):
    """Per cell: UCell (zeros at max_rank+1) and AUCell-style AUC with zeros tie-broken at random."""
    n_all = X.shape[1]
    k = int(np.ceil(k_frac * n_all))
    u = np.empty(X.shape[0])
    auc = np.empty(X.shape[0])
    n = len(cols)
    colset = np.array(cols)
    for c in range(X.shape[0]):
        s, e = X.indptr[c], X.indptr[c + 1]
        ind, val = X.indices[s:e], X.data[s:e]
        r = pd.Series(-val).rank(method="average").to_numpy()  # 1 = highest
        rank_of = dict(zip(ind, r, strict=True))
        nd = len(ind)
        rk_u = np.array([min(rank_of.get(g, max_rank + 1), max_rank + 1) for g in colset])
        u[c] = 1 - (rk_u.sum() - n * (n + 1) / 2) / (n * max_rank)
        # AUCell-style: zero genes take random ranks in nd+1..n_all
        rk_a = []
        for g in colset:
            if g in rank_of:
                rk_a.append(rank_of[g])
            else:
                rk_a.append(nd + 1 + rng.integers(0, n_all - nd))
        rk_a = np.array(rk_a)
        auc[c] = np.sum(np.clip(k - rk_a + 1, 0, None)) / (n * k)
    return u, auc


def within_lib_rho(score):
    out = {}
    for lb in np.unique(lib):
        m = lib == lb
        out[lb] = spearmanr(score[m], logtot[m]).statistic
    return pd.Series(out)


sg = score_genes(BORDER)
uc, au = rank_scores(idx)
res = pd.DataFrame(
    {
        "score_genes": within_lib_rho(sg).abs(),
        "ucell": within_lib_rho(uc).abs(),
        "aucell_style": within_lib_rho(au).abs(),
    }
)
print("P4 border axis: median within-library |rho| with log10 total counts")
print(res.median().round(3).to_string())
print(
    "libraries where UCell |rho| > score_genes |rho|:",
    int((res.ucell > res.score_genes).sum()),
    "of 22",
)
print(
    "libraries where AUCell-style |rho| > score_genes |rho|:",
    int((res.aucell_style > res.score_genes).sum()),
    "of 22",
)

# P1: border P90 per library vs 50 random sets matched on mean-expression bin
mean = np.asarray(b.X.mean(0)).ravel()
bins = pd.qcut(mean, 25, labels=False, duplicates="drop")
real_p90 = pd.Series(sg).groupby(lib).quantile(0.9)
null = []
for s in range(50):
    gl = [
        genes[rng.choice(np.where((bins == bins[i]) & ~np.isin(np.arange(len(genes)), idx))[0])]
        for i in idx
    ]
    null.append(pd.Series(score_genes(gl, seed=s + 1)).groupby(lib).quantile(0.9))
null = pd.concat(null, axis=1)
above = real_p90 > null.quantile(0.95, axis=1)
print(
    "\nP1 border score_genes P90 above 95th pct of 50 matched random sets:",
    int(above.sum()),
    "of 22 libraries",
)

# P7: nppa share in Hr11 vs uninjured cells, matched on total-count quintiles of Hr11
share = np.asarray(X[:, gi["nppa"]].todense()).ravel() / tot
ctrl = np.isin(lib, ["H5", "H6", "H7", "H8a", "H8v"])
h11 = lib == "Hr11"
edges = np.quantile(tot[h11], [0, 0.2, 0.4, 0.6, 0.8, 1])
print(
    "\nP7 nppa share of UMIs, Hr11 vs uninjured, matched on Hr11 depth quintiles (mean, mean, n_ctrl):"
)
for q in range(5):
    lo, hi = edges[q], edges[q + 1] + (1 if q == 4 else 0)
    mh = h11 & (tot >= lo) & (tot < hi)
    mc = ctrl & (tot >= lo) & (tot < hi)
    print(
        f"  bin {q} [{lo:.0f},{hi:.0f}): Hr11 {share[mh].mean():.4f}  ctrl {share[mc].mean():.4f}  ratio {share[mh].mean() / share[mc].mean():.2f}  n_ctrl {mc.sum()}"
    )
