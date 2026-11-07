"""Axis scoring: UCell on sparse counts and expression-matched random gene sets.

UCell follows Andreatta & Carmona 2021 (doi:10.1016/j.csbj.2021.06.043): genes are ranked
per cell in decreasing order, ties averaged, ranks above `max_rank` set to `max_rank + 1`, and
score = 1 - (sum of signature ranks - n(n+1)/2) / (n * max_rank).
Departure, stated: undetected genes are all given rank `max_rank + 1`, which equals UCell's
result whenever a cell has fewer than `max_rank` detected genes (true for nearly all cells here).
Re-implemented rather than installed (pyucell) to precompute ranks once for many null sets.
"""

from __future__ import annotations

import numpy as np
import scipy.sparse as sp
from scipy.stats import rankdata

from zf.panel import Panel


def rank_matrix(X) -> sp.csr_matrix:
    """Per-cell descending ranks of the non-zero entries (ties averaged), same sparsity as X."""
    X = sp.csr_matrix(X, dtype=float)
    R = X.copy()
    for i in range(X.shape[0]):
        s, e = X.indptr[i], X.indptr[i + 1]
        if e > s:
            R.data[s:e] = rankdata(-X.data[s:e], method="average")
    return R


def ucell_from_ranks(R: sp.csr_matrix, idx, max_rank: int = 1500) -> np.ndarray:
    """UCell score per cell for the gene columns `idx`, from `rank_matrix` output."""
    idx = np.asarray(idx)
    r = R[:, idx].toarray()
    r[r == 0] = max_rank + 1
    r = np.minimum(r, max_rank + 1)
    n = len(idx)
    return 1 - (r.sum(1) - n * (n + 1) / 2) / (n * max_rank)


def matched_random_sets(means, gene_idx, n_sets=100, n_bins=25, seed=0, exclude=()):
    """Random gene sets, each gene drawn from the same mean-expression quantile bin as the axis
    gene it replaces, without replacement within a set, excluding `exclude`."""
    means = np.asarray(means)
    edges = np.quantile(means, np.linspace(0, 1, n_bins + 1)[1:-1])
    bins = np.digitize(means, edges)
    excl = np.zeros(len(means), bool)
    excl[np.asarray(exclude, dtype=int)] = True
    pools = {b: np.flatnonzero((bins == b) & ~excl) for b in np.unique(bins[gene_idx])}
    rng = np.random.default_rng(seed)
    out = np.empty((n_sets, len(gene_idx)), dtype=int)
    for k in range(n_sets):
        used: set[int] = set()
        for j, g in enumerate(gene_idx):
            pool = np.setdiff1d(pools[bins[g]], list(used))
            out[k, j] = rng.choice(pool)
            used.add(int(out[k, j]))
    return out


# Panel-based scores used by zf.checks: a mean score with matched controls, and a rank score.

N_BINS = 25


def expression_bins(panel: Panel) -> np.ndarray:
    if panel.bins is None:
        means = panel.logn.mean(axis=0)
        edges = np.quantile(means, np.linspace(0, 1, N_BINS + 1)[1:-1])
        panel.bins = np.digitize(means, edges)
    return panel.bins


def score_mean(panel: Panel, genes) -> np.ndarray:
    """Set mean minus the mean of the other genes in the same expression bins (Tirosh 2016)."""
    idx = panel.idx(genes)
    bins = expression_bins(panel)
    ctrl = np.setdiff1d(np.flatnonzero(np.isin(bins, bins[idx])), idx)
    return panel.logn[:, idx].mean(axis=1) - panel.logn[:, ctrl].mean(axis=1)


def score_rank(panel: Panel, genes) -> np.ndarray:
    """Mann-Whitney U of the set's within-cell expression ranks, scaled to [0, 1] (UCell-like)."""
    idx = panel.idx(genes)
    order = panel.logn.argsort(axis=1, kind="stable")
    ranks = np.empty_like(order)
    rows = np.arange(order.shape[0])[:, None]
    ranks[rows, order] = np.arange(order.shape[1])[None, :]  # 0 = lowest
    n, m = len(idx), panel.logn.shape[1] - len(idx)
    u = ranks[:, idx].sum(axis=1) - n * (n - 1) / 2
    return u / (n * m)
