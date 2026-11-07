"""Labeling protocol v1: a fixed per-cell, per-axis call with artefact vetoes.

Pre-registered in JOURNAL.md and stated as rules in docs/labeling-protocol.md. For one gene set:

1. Thin every cell to T UMIs (binomial). Cells below T keep their counts and get `low_depth`;
   they receive only the rank (UCell) call.
2. Score twice: UCell on raw ranks, mean log-normalised expression on thinned counts.
3. Null: B decoy sets, each gene replaced by a random gene from the same detection-rate bin of
   the same library. p per method, scDRS-style: the cell's z against its own decoys, ranked
   against the pooled decoy z's of all cells in the library (a plain (1+#)/(1+B) p is floored
   at 1/501, which BH within a heart of thousands of cells can never call).
4. Call per heart: Benjamini-Hochberg on each method's p; positive when q <= Q_POS in both
   methods, negative when p >= P_NEG in both, otherwise uncertain.
5. Vetoes turn any call into uncertain: doublet flag, >= 2 non-CM (ambient) genes detected,
   hemoglobin fraction above the library 95th percentile, pct_mito > MITO_MAX.
6. Gate: per heart the positive fraction must beat the decoy false-positive rate (one-sided
   binomial test, p < 0.05). A timepoint is supported when >= floor(n/2)+1 of its n hearts pass.
   A positive is high-confidence only if its timepoint passes and still passes with any one
   heart left out.
7. Eligibility: if the within-heart partial Spearman rho of the thinned score with log UMI or
   with hemoglobin fraction has |rho| >= RHO_MAX, the axis is `depth_confounded` and emits no
   positives.

`label_axis(adata, genes)` is the one public entry point (team lane agreement, issue #10).
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy import stats

# Pre-registered constants (JOURNAL.md, 2026-11-07). Never tuned on results.
T = 600
B = 500
N_BINS = 20
MAX_RANK = 1500
Q_POS = 0.05
P_NEG = 0.10
MITO_MAX = 30.0
HB_PCT = 95
RHO_MAX = 0.10
GATE_P = 0.05
N_DECOY_AXES = 20  # decoy sets re-called as pseudo-axes to estimate the false-positive rate
POOL = 50  # decoy columns whose z's are pooled across the library's cells for the p-value
RATE_FLOOR = 1e-3  # decoy rate used when no decoy cell is called positive
SEED = 0


def thin(counts: sp.csr_matrix, target: int, rng: np.random.Generator):
    """Binomial thinning of each cell to `target` UMIs; cells below it are left as they are."""
    x = sp.csr_matrix(counts, dtype=np.float64, copy=True)
    total = np.asarray(x.sum(1)).ravel()
    keep = np.where(total > target, target / np.maximum(total, 1), 1.0)
    p = np.repeat(keep, np.diff(x.indptr))
    x.data = rng.binomial(x.data.astype(np.int64), p).astype(np.float64)
    x.eliminate_zeros()
    return x, total < target


def lognorm(x: sp.csr_matrix) -> sp.csr_matrix:
    total = np.asarray(x.sum(1)).ravel()
    y = sp.csr_matrix(x.multiply(1e4 / np.maximum(total, 1)[:, None]))
    return y.log1p()


def ucell_d(x: sp.csr_matrix, max_rank: int = MAX_RANK) -> sp.csr_matrix:
    """Per cell, max_rank - min(rank, max_rank) on expressed genes (rank 1 = highest).

    UCell(S) = (D @ 1_S + n(n+1)/2) / (n * max_rank) with n = |S|, because an undetected gene
    sits at max_rank and contributes 0 to D.
    """
    x = sp.csr_matrix(x)
    d = x.copy().astype(np.float64)
    for i in range(x.shape[0]):
        lo, hi = x.indptr[i], x.indptr[i + 1]
        if hi > lo:
            r = stats.rankdata(-x.data[lo:hi], method="average")
            d.data[lo:hi] = max_rank - np.minimum(r, max_rank)
    d.eliminate_zeros()
    return d


def _indicator(sets: np.ndarray, n_genes: int) -> sp.csc_matrix:
    """Gene sets (sets x slots, gene indices) -> genes x sets 0/1 matrix."""
    k, n = sets.shape
    rows = sets.ravel()
    cols = np.repeat(np.arange(k), n)
    return sp.csc_matrix((np.ones(k * n), (rows, cols)), shape=(n_genes, k))


def decoy_sets(detect: np.ndarray, idx: np.ndarray, rng, b: int = B, bins: int = N_BINS):
    """b sets of len(idx) genes, each slot drawn from the same detection-rate bin as idx."""
    expressed = np.flatnonzero(detect > 0)
    edges = np.quantile(detect[expressed], np.linspace(0, 1, bins + 1)[1:-1])
    gene_bin = np.full(len(detect), -1)
    gene_bin[expressed] = np.searchsorted(edges, detect[expressed], side="right")
    out = np.empty((b, len(idx)), dtype=np.int64)
    for j, g in enumerate(idx):
        # An undetected set gene is matched by undetected genes. Drawing detected genes for it made the
        # decoys outscore the set wherever set genes are absent (found 2026-11-07, plan section H check).
        pool = np.flatnonzero(gene_bin == gene_bin[g]) if gene_bin[g] >= 0 else np.flatnonzero(detect == 0)
        pool = pool[pool != g]
        out[:, j] = rng.choice(pool, size=b)
    return out


def _bh(p: np.ndarray) -> np.ndarray:
    return stats.false_discovery_control(p, method="bh") if len(p) else p


def _pooled_p(obs: np.ndarray, null: np.ndarray, rng, pool: int = POOL) -> np.ndarray:
    """scDRS-style p: each cell's z against its own decoys, ranked against the decoy z's of
    every cell in the library (pooled), so p is not floored at 1 / (1 + B)."""
    mu = null.mean(1, keepdims=True)
    sd = null.std(1, keepdims=True) + 1e-12
    z = (obs[:, None] - mu)[:, 0] / sd[:, 0]
    cols = rng.choice(null.shape[1], size=min(pool, null.shape[1]), replace=False)
    zn = np.sort(((null[:, cols] - mu) / sd).ravel())
    return (1 + len(zn) - np.searchsorted(zn, z, side="left")) / (1 + len(zn))


def _partial_spearman(y, x, z, group) -> float:
    """Spearman rho of y and x controlling for z, on ranks taken within each group."""
    df = pd.DataFrame({"y": y, "x": x, "z": z, "g": np.asarray(group)})
    r = df.groupby("g", observed=True)[["y", "x", "z"]].rank(pct=True)
    r = (r - r.groupby(df["g"]).transform("mean")).to_numpy()
    sd = r.std(0)
    if sd[0] == 0 or sd[1] == 0:
        return 0.0
    c = np.eye(3)
    for i, j in [(0, 1), (0, 2), (1, 2)]:
        if sd[i] > 0 and sd[j] > 0:
            c[i, j] = c[j, i] = np.corrcoef(r[:, i], r[:, j])[0, 1]
    num = c[0, 1] - c[0, 2] * c[1, 2]
    den = np.sqrt(max((1 - c[0, 2] ** 2) * (1 - c[1, 2] ** 2), 1e-12))
    return float(num / den) if np.isfinite(num) else 0.0


def _call(p_rank, p_score, low_depth, heart):
    """BH within heart on each method; both must agree (rank only for low-depth cells)."""
    q_rank = np.empty_like(p_rank)
    q_score = np.empty_like(p_score)
    for h in np.unique(heart):
        m = heart == h
        q_rank[m] = _bh(p_rank[m])
        q_score[m] = _bh(p_score[m])
    pos = (q_rank <= Q_POS) & ((q_score <= Q_POS) | low_depth)
    neg = (p_rank >= P_NEG) & ((p_score >= P_NEG) | low_depth)
    return np.where(pos, "positive", np.where(neg, "negative", "uncertain")), q_rank, q_score


def gate(passed: pd.Series) -> bool:
    """>= floor(n/2)+1 of n hearts pass."""
    n = len(passed)
    return n > 0 and int(passed.sum()) >= n // 2 + 1


def _veto(obs: pd.DataFrame, counts: sp.csr_matrix, var_names, ambient_genes) -> pd.DataFrame:
    v = pd.DataFrame(index=obs.index)
    v["veto_doublet"] = obs.get("doublet_flag", pd.Series(False, index=obs.index)).astype(bool)
    amb = [var_names.get_loc(g) for g in ambient_genes if g in var_names]
    n_amb = np.diff(sp.csr_matrix(counts[:, amb]).indptr) if amb else np.zeros(len(obs))
    v["n_nonCM_genes"] = n_amb
    v["veto_nonCM"] = n_amb >= 2
    hb = np.flatnonzero(var_names.str.match(r"^hb[ab]"))
    total = np.asarray(counts.sum(1)).ravel()
    hb_frac = np.asarray(counts[:, hb].sum(1)).ravel() / np.maximum(total, 1)
    v["hb_frac"] = hb_frac
    lib = obs["orig_ident"].astype(str).to_numpy()
    cut = pd.Series(hb_frac).groupby(lib).transform(lambda s: np.percentile(s, HB_PCT))
    v["veto_hb"] = (hb_frac > cut.to_numpy()) & (hb_frac > 0)
    mito = obs["pct_mito"] if "pct_mito" in obs else pd.Series(0.0, index=obs.index)
    v["veto_mito"] = mito.to_numpy() > MITO_MAX
    v["vetoed"] = v[["veto_doublet", "veto_nonCM", "veto_hb", "veto_mito"]].any(axis=1)
    return v


def label_axis(
    adata,
    genes,
    *,
    seed: int = SEED,
    ambient_genes=(),
    target: int = T,
    n_decoys: int = B,
    heart_key: str = "heart",
    time_key: str = "condition",
    library_key: str = "orig_ident",
    precomputed: dict | None = None,
) -> pd.DataFrame:
    """Per-cell call, both scores, empirical p, flags; axis verdict in `.attrs`.

    adata: cells x genes, raw counts in layers['counts'] (or X); obs needs library, heart and
    time columns (`doublet_flag` and `pct_mito` are used when present).
    `.attrs['verdict']` is one of 'supported', 'not supported', 'depth_confounded';
    `.attrs['gate']` holds the per-heart table as records; `.attrs['rho']` the eligibility rhos.
    """
    rng = np.random.default_rng(seed)
    obs = adata.obs
    var_names = pd.Index(adata.var_names)
    idx = np.array([var_names.get_loc(g) for g in genes if g in var_names])
    if len(idx) < 2:
        raise ValueError(f"only {len(idx)} of {len(genes)} genes are in the data")
    pre = precomputed or prepare(adata, seed=seed, target=target, library_key=library_key)
    counts, thinned_log, d_rank, low = pre["counts"], pre["thinned_log"], pre["d_rank"], pre["low"]
    lib = obs[library_key].astype(str).to_numpy()
    heart = obs[heart_key].astype(str).to_numpy()
    n_cells, n_genes = counts.shape
    n = len(idx)

    s_rank = np.zeros(n_cells)
    s_score = np.zeros(n_cells)
    p_rank = np.zeros(n_cells)
    p_score = np.zeros(n_cells)
    null_rank = np.zeros((n_cells, N_DECOY_AXES))
    null_score = np.zeros((n_cells, N_DECOY_AXES))
    pn_rank = np.zeros((n_cells, N_DECOY_AXES))
    pn_score = np.zeros((n_cells, N_DECOY_AXES))
    const = n * (n + 1) / 2
    for L in np.unique(lib):
        m = np.flatnonzero(lib == L)
        detect = np.asarray((counts[m] > 0).mean(0)).ravel()
        sets = np.vstack([idx[None, :], decoy_sets(detect, idx, rng, b=n_decoys)])
        ind = _indicator(sets, n_genes)
        r = (np.asarray((d_rank[m] @ ind).todense()) + const) / (n * MAX_RANK)
        s = np.asarray((thinned_log[m] @ ind).todense()) / n
        s_rank[m], s_score[m] = r[:, 0], s[:, 0]
        p_rank[m] = _pooled_p(r[:, 0], r[:, 1:], rng)
        p_score[m] = _pooled_p(s[:, 0], s[:, 1:], rng)
        # Decoy pseudo-axes: decoy k against the other decoys, for the false-positive rate.
        for k in range(N_DECOY_AXES):
            others = np.r_[1 : 1 + k, 2 + k : r.shape[1]]
            null_rank[m, k], null_score[m, k] = r[:, 1 + k], s[:, 1 + k]
            pn_rank[m, k] = _pooled_p(r[:, 1 + k], r[:, others], rng)
            pn_score[m, k] = _pooled_p(s[:, 1 + k], s[:, others], rng)

    call, q_rank, q_score = _call(p_rank, p_score, low, heart)
    veto = _veto(obs, counts, var_names, ambient_genes)
    call = np.where(veto["vetoed"].to_numpy(), "uncertain", call)

    decoy_pos = np.zeros((n_cells, N_DECOY_AXES), dtype=bool)
    for k in range(N_DECOY_AXES):
        c, _, _ = _call(pn_rank[:, k], pn_score[:, k], low, heart)
        decoy_pos[:, k] = (c == "positive") & ~veto["vetoed"].to_numpy()

    log_umi = np.log10(np.maximum(np.asarray(counts.sum(1)).ravel(), 1))
    hb = veto["hb_frac"].to_numpy()
    rho_umi = _partial_spearman(s_score, log_umi, hb, heart)
    rho_hb = _partial_spearman(s_score, hb, log_umi, heart)
    confounded = max(abs(rho_umi), abs(rho_hb)) >= RHO_MAX

    out = pd.DataFrame(
        {
            "call": call,
            "score_ucell_raw": s_rank,
            "score_mean_thinned": s_score,
            "p_ucell": p_rank,
            "p_score": p_score,
            "q_ucell": q_rank,
            "q_score": q_score,
            "low_depth": low,
        },
        index=obs.index,
    ).join(veto)
    g = heart_gate(out, decoy_pos, obs[heart_key], obs[time_key])
    if confounded:
        out["call"] = np.where(out["call"] == "positive", "uncertain", out["call"])
        verdict = "depth_confounded"
    else:
        verdict = "supported" if g.groupby("time")["supported"].first().any() else "not supported"
    t_ok = obs[time_key].astype(str).map(g.groupby("time")["robust"].first()).fillna(False)
    out["high_confidence"] = (out["call"] == "positive") & t_ok.to_numpy().astype(bool)
    out.attrs.update(
        verdict=verdict, gate=g.to_dict("records"), rho={"log_umi": rho_umi, "hb_frac": rho_hb}
    )
    return out


def heart_gate(out: pd.DataFrame, decoy_pos: np.ndarray, heart, time) -> pd.DataFrame:
    """Per heart: positive fraction vs decoy rate (binomial); per time: gate and LOHO."""
    df = pd.DataFrame(
        {
            "heart": np.asarray(heart).astype(str),
            "time": np.asarray(time).astype(str),
            "pos": out["call"].eq("positive").to_numpy(),
        }
    )
    df["decoy_rate"] = decoy_pos.mean(1)
    rows = []
    for (h, t), d in df.groupby(["heart", "time"]):
        k, n = int(d["pos"].sum()), len(d)
        rate = max(float(d["decoy_rate"].mean()), RATE_FLOOR)
        p = stats.binomtest(k, n, rate, alternative="greater").pvalue
        rows.append(
            {
                "heart": h,
                "time": t,
                "n_cells": n,
                "pos_frac": k / n,
                "decoy_rate": float(d["decoy_rate"].mean()),
                "binom_p": p,
                "pass": p < GATE_P,
            }
        )
    g = pd.DataFrame(rows)
    sup, rob = {}, {}
    for t, d in g.groupby("time"):
        sup[t] = gate(d["pass"])
        rob[t] = sup[t] and all(gate(d["pass"].drop(i)) for i in d.index) if len(d) > 1 else False
    g["supported"] = g["time"].map(sup)
    g["robust"] = g["time"].map(rob)
    return g


def prepare(adata, *, seed: int = SEED, target: int | None = T, library_key: str = "orig_ident"):
    """Shared per-cell matrices: counts, thinned log-normalised counts, UCell rank offsets."""
    counts = adata.layers.get("counts", adata.X)
    counts = sp.csr_matrix(counts)
    if target is None:  # no thinning: the reference for the thinning-agreement metric
        thinned, low = counts.astype(np.float64), np.zeros(counts.shape[0], dtype=bool)
    else:
        thinned, low = thin(counts, target, np.random.default_rng(seed + 1))
    return {
        "counts": counts,
        "thinned": thinned,
        "thinned_log": lognorm(thinned),
        "d_rank": ucell_d(counts),
        "low": low,
    }
