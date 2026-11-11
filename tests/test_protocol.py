"""Labeling protocol v1 on synthetic data: null calibration, planted defects, a true state."""

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy import stats

from zf import protocol

N_GENES = 400


def _synthetic(depth_genes=(), state_genes=(), state_frac=0.0, seed=0):
    """6 hearts x 300 cells, Poisson counts; depth varies 4x within each heart.

    depth_genes: their share rises with depth (a fake marker that tracks library size).
    state_genes: up 6x in a random state_frac of cells in hearts h0-h3 (time 'late').
    """
    rng = np.random.default_rng(seed)
    base = rng.gamma(0.6, 1.0, N_GENES)
    base /= base.sum()
    rows, obs = [], []
    for h in range(6):
        depth = rng.uniform(600, 2400, 300)
        for d in depth:
            w = base.copy()
            if len(depth_genes):
                w[list(depth_genes)] *= (d / 600) ** 1.5
            if len(state_genes) and h < 4 and rng.random() < state_frac:
                w[list(state_genes)] *= 6
            rows.append(rng.poisson(d * w / w.sum()))
            obs.append(
                {"orig_ident": f"L{h}", "heart": f"h{h}", "condition": "late" if h < 4 else "early"}
            )
    x = sp.csr_matrix(np.array(rows, dtype=np.float32))
    a = ad.AnnData(X=x, obs=pd.DataFrame(obs, index=[f"c{i}" for i in range(x.shape[0])]))
    a.var_names = [f"g{i}" for i in range(N_GENES)]
    a.layers["counts"] = x
    return a


def test_null_gene_set_has_uniform_p_and_is_not_supported():
    a = _synthetic()
    genes = [f"g{i}" for i in range(20, 40)]
    out = protocol.label_axis(a, genes, n_decoys=200)
    # Matched decoys: a random gene set's empirical p is uniform (KS p > 0.01).
    assert stats.kstest(out["p_score"], "uniform").pvalue > 0.01
    assert out.attrs["verdict"] == "not supported"


def test_planted_depth_tracking_markers_are_depth_confounded():
    depth = list(range(50, 70))  # 20 genes picked to track depth
    a = _synthetic(depth_genes=depth)
    out = protocol.label_axis(a, [f"g{i}" for i in depth], n_decoys=200)
    assert out.attrs["verdict"] == "depth_confounded", out.attrs["rho"]
    assert not out["call"].eq("positive").any()


def test_true_state_is_supported_and_called_in_its_cells():
    state = list(range(100, 115))
    a = _synthetic(state_genes=state, state_frac=0.3)
    out = protocol.label_axis(a, [f"g{i}" for i in state])
    assert out.attrs["verdict"] == "supported", out.attrs["gate"]
    late = a.obs["condition"].eq("late").to_numpy()
    # Conservative by design: BH over 300 cells with p >= 1/501 finds part of the 30%.
    assert out["call"][late].eq("positive").mean() > 0.08
    assert out["call"][~late].eq("positive").mean() < 0.02


def test_thinning_hits_target_and_flags_shallow_cells():
    x = sp.csr_matrix(np.array([[300, 300, 400], [100, 100, 0]], dtype=float))
    t, low = protocol.thin(x, 600, np.random.default_rng(0))
    assert t[0].sum() <= 1000 and abs(t[0].sum() - 600) < 80
    assert list(low) == [False, True]
    assert t[1].sum() == 200


def test_gate_rule_is_majority_plus_one():
    assert protocol.gate(pd.Series([True, True, True, False]))
    assert not protocol.gate(pd.Series([True, True, False, False]))
    assert protocol.gate(pd.Series([True, True, False]))
    assert protocol.gate(pd.Series([True] * 4 + [False] * 3))


def test_undetected_set_gene_is_matched_by_undetected_decoys():
    # A set gene absent from a heart adds 0 to every cell score; a detected decoy in its slot would
    # outscore the set and make every presence test conservative (found 2026-11-07).
    detect = np.r_[np.zeros(50), np.linspace(0.01, 0.5, 200)]
    idx = np.array([0, 60, 120])
    dec = protocol.decoy_sets(detect, idx, np.random.default_rng(0), b=100)
    assert (detect[dec[:, 0]] == 0).all()
    assert (detect[dec[:, 1:]] > 0).all()
