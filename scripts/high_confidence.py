"""Statistical plan v2, section H: high-confidence presence without asking every heart to match.

Pre-registered: committed with plan section H before its first run.

The task owner (2026-11-07): cutoffs need not demand that every library matches; define cutoffs of
high confidence, and find where the states form. Three levels, each with its own null:
- cell: "high" if its set score (mean log-normalised count over the set) is above the 95th
  percentile of the scores of 1,000 detection-matched decoy gene sets in the same heart and chamber
  (cells thinned to 600 UMIs, veto applied; 15% of ventricle cells are below 600 and drop out);
- heart: the share of high cells, compared with each decoy set's own share of high cells;
  p = (1 + decoys at or above the set) / 1001, floor 0.001. The null is competitive: a co-regulated
  set can show high cells in a healthy heart, so a positive heart means "more high cells than random
  genes", not "injury state";
- condition: partial conjunction (Benjamini & Heller 2008), Simes combination,
  P(r/n) = min over r <= i <= n of (n - r + 1) / (i - r + 1) * p(i), and the 95% lower bound
  r_hat = max{r : P(r/n) <= 0.05} on the number of hearts with high cells in excess (Wang & Owen
  2019). Hearts are independent animals.
Reported per set, condition and chamber: r_hat, and P at the majority r = floor(n/2) + 1 with BH
over the six primary sets x four conditions x two chambers. "Injury-specific presence": r_hat >= 4
of 7 at
7 dpi, r_hat = 0 in the uninjured hearts, and the section A test positive at 7 dpi.
Controls: dissociation set (positive: handling stress is present in every heart); 20 random 10-gene
sets (negative: should reach r_hat = 0 almost always).

Usage: PYTHONPATH=src:scripts ZF_DATA=<data> .venv/bin/python scripts/high_confidence.py
"""

import numpy as np
import pandas as pd
import scipy.sparse as sp
import yaml
from pseudobulk_de import OUT
from pseudobulk_set_score import PRIMARY
from pseudobulk_set_test import marker_sets
from run_protocol import ROOT, strata
from scipy.stats import false_discovery_control

from zf import io, protocol

N_DECOY, N_RANDOM, SEED, ALPHA, MIN_CELLS = 1000, 20, 20261109, 0.05, 30


def heart_p(cnt, lg, idx, rng):
    """Excess share of high cells and its decoy p-value for one set in one heart x chamber."""
    detect = np.asarray((cnt > 0).mean(0)).ravel()
    sets = np.vstack([idx[None, :], protocol.decoy_sets(detect, idx, rng, b=N_DECOY)])
    s = np.asarray((lg @ protocol._indicator(sets, cnt.shape[1])).todense()) / len(idx)
    high = (s > np.quantile(s[:, 1:], 0.95)).mean(0)
    return high[0] - high[1:].mean(), (1 + np.sum(high[1:] >= high[0])) / (1 + N_DECOY)


def partial_conjunction(p, r):
    """Simes partial-conjunction p for 'at least r of n' (Benjamini & Heller 2008)."""
    p, n = np.sort(p), len(p)
    i = np.arange(r, n + 1)
    return float(min(1.0, np.min((n - r + 1) / (i - r + 1) * p[i - 1])))


def main():
    a = io.load()
    counts = sp.csr_matrix(a.layers["counts"])
    st = strata(a)
    controls = yaml.safe_load((ROOT / "config" / "axes.yaml").read_text())["controls"]
    fib = [
        g["symbol"] for g in next(c for c in controls if c["name"] == "ambient_fibroblast")["genes"]
    ]
    veto = protocol._veto(a.obs, counts, pd.Index(a.var_names), fib)["vetoed"].to_numpy()
    thinned, low = protocol.thin(counts, protocol.T, np.random.default_rng(protocol.SEED + 1))
    lgall = protocol.lognorm(thinned)
    var = pd.Index(a.var_names)
    ms = marker_sets()
    rng = np.random.default_rng(SEED)
    heart = a.obs["heart"].astype(str).to_numpy()
    cond = a.obs["condition"].astype(str).to_numpy()
    sets = {n: ms[n] for n in PRIMARY + ["dissociation_control"]}
    detected = var[np.asarray((thinned[~veto & ~low] > 0).mean(0)).ravel() > 0.01]
    sets |= {
        f"random_{i:02d}": list(rng.choice(detected, 10, replace=False)) for i in range(N_RANDOM)
    }

    rows = []
    for chamber in ["ventricle", "atrium"]:
        use = st[chamber] & ~veto & ~low
        for h in sorted(set(heart[use])):
            m = np.flatnonzero(use & (heart == h))
            if len(m) < MIN_CELLS:
                continue
            cnt, lg = sp.csr_matrix(thinned[m]), sp.csr_matrix(lgall[m])
            for name, genes in sets.items():
                idx = var.get_indexer([g for g in genes if g in var])
                ex, p = heart_p(cnt, lg, idx, rng)
                rows.append(
                    {
                        "chamber": chamber,
                        "heart": h,
                        "condition": cond[m][0],
                        "n_cells": len(m),
                        "set": name,
                        "excess_high": ex,
                        "p_heart": p,
                    }
                )
    per = pd.DataFrame(rows)

    res = []
    for (chamber, c, name), g in per.groupby(["chamber", "condition", "set"]):
        p, n = g["p_heart"].to_numpy(), len(g)
        pcs = {r: partial_conjunction(p, r) for r in range(1, n + 1)}
        res.append(
            {
                "chamber": chamber,
                "condition": c,
                "set": name,
                "primary": name in PRIMARY,
                "n_hearts": n,
                "r_majority": n // 2 + 1,
                "p_majority": pcs[n // 2 + 1],
                "r_hat": max([r for r in pcs if pcs[r] <= ALPHA], default=0),
                "hearts_p_le_0.05": int((p <= ALPHA).sum()),
                "median_excess": g["excess_high"].median(),
                **{f"p_r{r}": v for r, v in pcs.items()},
            }
        )
    res = pd.DataFrame(res)
    prim = res["primary"]
    res.loc[prim, "q_majority"] = false_discovery_control(res.loc[prim, "p_majority"].to_numpy())
    res["present_majority"] = res["q_majority"] <= ALPHA
    per.to_csv(OUT / "high_confidence_per_heart.csv", index=False)
    res.to_csv(OUT / "high_confidence.csv", index=False)
    rnd = res[res["set"].str.startswith("random_")]
    print(
        f"negative control: random sets with r_hat >= 1 in {int((rnd['r_hat'] >= 1).sum())}"
        f" of {len(rnd)} condition x chamber cells"
    )
    cols = [
        "chamber",
        "condition",
        "set",
        "n_hearts",
        "hearts_p_le_0.05",
        "r_hat",
        "p_majority",
        "q_majority",
        "median_excess",
    ]
    with pd.option_context("display.width", 200, "display.max_columns", 20):
        print(res[~res["set"].str.startswith("random_")][cols].round(4).to_string(index=False))


if __name__ == "__main__":
    main()
