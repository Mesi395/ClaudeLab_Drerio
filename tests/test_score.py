"""Tests for zf.score: UCell on sparse counts and matched random gene sets (C002)."""

import numpy as np
import scipy.sparse as sp
from scipy.stats import rankdata

from zf.score import matched_random_sets, rank_matrix, ucell_from_ranks


def _ucell_bruteforce(x, idx, max_rank):
    # UCell (Andreatta 2021): rank genes descending, ties averaged, ranks > max_rank -> max_rank+1;
    # score = 1 - (sum ranks - n(n+1)/2) / (n * max_rank)
    r = rankdata(-x, method="average")
    r[x == 0] = max_rank + 1  # all undetected genes tie beyond max_rank
    r = np.minimum(r, max_rank + 1)
    n = len(idx)
    return 1 - (r[idx].sum() - n * (n + 1) / 2) / (n * max_rank)


def test_ucell_matches_bruteforce():
    rng = np.random.default_rng(1)
    dense = rng.poisson(0.5, size=(30, 200)).astype(float)
    X = sp.csr_matrix(dense)
    idx = np.array([0, 5, 17, 40, 199])
    R = rank_matrix(X)
    got = ucell_from_ranks(R, idx, max_rank=50)
    want = np.array([_ucell_bruteforce(dense[i], idx, 50) for i in range(30)])
    np.testing.assert_allclose(got, want, atol=1e-12)


def test_ucell_invariant_to_cell_scaling():
    rng = np.random.default_rng(2)
    dense = rng.poisson(1.0, size=(10, 100)).astype(float)
    idx = np.arange(5)
    a = ucell_from_ranks(rank_matrix(sp.csr_matrix(dense)), idx, 40)
    b = ucell_from_ranks(rank_matrix(sp.csr_matrix(dense * 7.0)), idx, 40)
    np.testing.assert_allclose(a, b)


def test_matched_random_sets_bins_and_exclusion():
    means = np.linspace(0, 1, 1000)
    gene_idx = np.array([10, 500, 990])
    sets = matched_random_sets(means, gene_idx, n_sets=20, n_bins=10, seed=0, exclude=gene_idx)
    assert sets.shape == (20, 3)
    bins = np.digitize(means, np.quantile(means, np.linspace(0, 1, 11)[1:-1]))
    for s in sets:
        assert list(bins[s]) == list(bins[gene_idx])
        assert not set(s) & set(gene_idx)
        assert len(set(s)) == 3
