"""Do the four states co-occur in the same cells? Card B's last step (TASKS.md score-3).

Per heart and axis pair (A, B), on the label-free cardiomyocytes:
- rho: within-heart partial Spearman of the two thinned scores (protocol v1) given log UMI.
- null: the same rho with B replaced by K decoy sets of B, matched on detection rate within each
  library exactly as the protocol's decoys. Two disjoint gene sets correlate a little through
  composition alone, so the test is rho against its own decoys:
  p_hi = (1 + #null >= rho) / (1 + K), p_lo = (1 + #null <= rho) / (1 + K).
- joint positives: cells positive on both axes vs the product of the two positive fractions
  (calls from data/derived/labels_v1.tsv.gz, written by scripts/run_protocol.py).
Per timepoint, the protocol's gate rule: the pair co-occurs when >= floor(n/2)+1 of its n hearts
have p_hi < 0.05, and comes apart when that many have p_lo < 0.05; otherwise "not shown".
Secondary (added after the first run, JOURNAL 2026-11-07 16:40): sign test of the per-heart excess
over hearts, per timepoint and over all hearts (time "all"). Positive control: the dissociation
control set split into two halves (pair "dissociation_half_a~dissociation_half_b"), which share
one program and must co-occur if the test can see co-occurrence at all.

Usage: uv run python scripts/cooccurrence.py   (after scripts/run_protocol.py)
"""

import itertools
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp
from run_protocol import axis_genes, dissoc_genes
from scipy import stats

from zf import io, protocol, provenance

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "hk"
LABELS = io.DATA / "derived" / "labels_v1.tsv.gz"
K = 100  # decoys of B per pair; p floor 1/101 is below the 0.05 per-heart cut
P_HEART = 0.05


def main() -> None:
    a = io.load()
    counts = sp.csr_matrix(a.layers["counts"])
    labels = pd.read_csv(LABELS, sep="\t", index_col=0).loc[a.obs_names]
    # Same thinning and seed as protocol.prepare, so the real scores match labels_v1 exactly.
    thinned, _ = protocol.thin(counts, protocol.T, np.random.default_rng(protocol.SEED + 1))
    tlog = protocol.lognorm(thinned)
    var = pd.Index(a.var_names)
    lib = a.obs["orig_ident"].astype(str).to_numpy()
    heart = a.obs["heart"].astype(str).to_numpy()
    time = a.obs["condition"].astype(str).to_numpy()
    log_umi = np.log10(np.maximum(np.asarray(counts.sum(1)).ravel(), 1))
    axes = axis_genes()
    scored = {g for gs in axes.values() for g in gs}
    dis = dissoc_genes(a.var_names, scored)
    control = {"dissociation_half_a": dis[0::2], "dissociation_half_b": dis[1::2]}
    sets = {**axes, **control}
    genes = {ax: np.array([var.get_loc(g) for g in gs if g in var]) for ax, gs in sets.items()}
    score = {ax: np.asarray(tlog[:, idx].mean(1)).ravel() for ax, idx in genes.items()}
    for ax in axes:
        if not np.allclose(score[ax], labels[f"{ax}.score_mean_thinned"].to_numpy(), atol=1e-9):
            raise ValueError(f"{ax}: recomputed score differs from {LABELS}; rerun run_protocol")
    none = np.zeros(a.n_obs, dtype=bool)  # control halves have no protocol calls
    pos = {ax: labels[f"{ax}.call"].eq("positive").to_numpy() for ax in axes} | dict.fromkeys(
        control, none
    )
    pairs = [*itertools.combinations(axes, 2), tuple(control)]

    rng = np.random.default_rng(protocol.SEED)
    decoy = {}  # axis -> cells x K decoy scores, decoys drawn per library
    for ax, idx in genes.items():
        d = np.zeros((a.n_obs, K))
        for L in np.unique(lib):
            m = np.flatnonzero(lib == L)
            detect = np.asarray((counts[m] > 0).mean(0)).ravel()
            ind = protocol._indicator(protocol.decoy_sets(detect, idx, rng, b=K), a.n_vars)
            d[m] = np.asarray((tlog[m] @ ind).todense()) / len(idx)
        decoy[ax] = d

    rows = []
    for ax_a, ax_b in pairs:
        for h in np.unique(heart):
            m = heart == h
            g = heart[m]
            rho = protocol._partial_spearman(score[ax_a][m], score[ax_b][m], log_umi[m], g)
            null = np.array(
                [
                    protocol._partial_spearman(score[ax_a][m], decoy[ax_b][m, k], log_umi[m], g)
                    for k in range(K)
                ]
            )
            pa, pb = pos[ax_a][m].mean(), pos[ax_b][m].mean()
            rows.append(
                {
                    "pair": f"{ax_a}~{ax_b}",
                    "heart": h,
                    "time": time[m][0],
                    "n_cells": int(m.sum()),
                    "rho": rho,
                    "null_median": float(np.median(null)),
                    "excess": rho - float(np.median(null)),
                    "p_hi": (1 + (null >= rho).sum()) / (1 + K),
                    "p_lo": (1 + (null <= rho).sum()) / (1 + K),
                    "pos_a": int(pos[ax_a][m].sum()),
                    "pos_b": int(pos[ax_b][m].sum()),
                    "pos_both": int((pos[ax_a][m] & pos[ax_b][m]).sum()),
                    "pos_both_expected": float(pa * pb * m.sum()),
                }
            )
    ph = pd.DataFrame(rows)
    ph.to_csv(OUT / "cooccurrence_per_heart.csv", index=False)

    summ = []
    groups = [*ph.groupby(["pair", "time"]), *(((k, "all"), d) for k, d in ph.groupby("pair"))]
    for (pair, t), d in groups:
        hi, lo = d["p_hi"] < P_HEART, d["p_lo"] < P_HEART
        verdict = (
            "co-occur" if protocol.gate(hi) else "come apart" if protocol.gate(lo) else "not shown"
        )
        up = int((d["excess"] > 0).sum())
        summ.append(
            {
                "pair": pair,
                "time": t,
                "n_hearts": len(d),
                "rho_median": d["rho"].median(),
                "excess_median": d["excess"].median(),
                "hearts_p_hi": int(hi.sum()),
                "hearts_p_lo": int(lo.sum()),
                "hearts_excess_pos": up,
                "sign_p_two_sided": stats.binomtest(up, len(d), 0.5).pvalue,
                "pos_both": int(d["pos_both"].sum()),
                "pos_both_expected": d["pos_both_expected"].sum(),
                "verdict": verdict,
            }
        )
    s = pd.DataFrame(summ)
    s.to_csv(OUT / "cooccurrence.csv", index=False)
    provenance.write_manifest(
        ROOT / "results" / "runs" / "hk-cooccurrence",
        command=" ".join(["uv run python scripts/cooccurrence.py", *sys.argv[1:]]),
        inputs=[io.LABEL_FREE, io.FLAGS, LABELS],
        params={"K": K, "p_heart": P_HEART, "T": protocol.T, "seed": protocol.SEED},
        outcomes={},
    )
    print(s.to_string())


if __name__ == "__main__":
    main()
