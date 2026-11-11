"""Statistical plan v2, section 2F: does injury shift every cell a little, or switch on a subgroup?

Pre-registered: committed with the plan's section 2F before its first run.

Per ventricle heart, cells at >= 600 UMIs thinned to exactly 600 (equal detection), veto applied.
Per gene set, against 200 detection-matched decoy sets drawn within that heart:
- tail excess = share of cells whose set score (mean log-normalised count over the set) is above
  the 95th percentile of the heart's decoy scores, minus the decoys' own share above it (null 0);
- dispersion excess = log(variance / mean) of the set's summed counts across cells, minus the
  median of the same over decoys. A subgroup switching on adds variance between cells; a small
  shift of every cell barely changes it.
Both measures are shape measures: decoys are matched on detection within the heart, so a uniform
level shift of every cell is absorbed (synthetic check: null tail -0.005 / disp +0.02; 10% of cells
x5 tail +0.07 / disp +1.58; every cell x1.5 tail +0.006 / disp +0.06). The level is item 18's test.
Hearts are compared with the plan's test (Welch t, one-sided up, exact relabelling): 7 dpi vs
control is the primary family (BH over the six primary sets), 3 dpi vs control the second.
Decision: "a subgroup switches on" if tail excess q <= 0.05 and dispersion excess p <= 0.05 at the
same contrast; "every cell shifts" if the heart-level set score shifted (item 18) and the
dispersion excess did not; otherwise "cannot tell at 600 UMIs".

Usage: PYTHONPATH=src:scripts ZF_DATA=<data> .venv/bin/python scripts/cell_tail_test.py
"""

import numpy as np
import pandas as pd
import scipy.sparse as sp
import yaml
from pseudobulk_de import OUT
from pseudobulk_set_plan import CONTRASTS, STEP, test
from pseudobulk_set_score import PRIMARY
from pseudobulk_set_test import marker_sets
from run_protocol import ROOT, strata
from scipy.stats import false_discovery_control

from zf import io, protocol

N_DECOY, SEED, ALPHA = 200, 20261107, 0.05


def shape(cnt, lg, idx, rng):
    """Tail excess and dispersion excess of one set in one heart's cells."""
    detect = np.asarray((cnt > 0).mean(0)).ravel()
    sets = np.vstack([idx[None, :], protocol.decoy_sets(detect, idx, rng, b=N_DECOY)])
    ind = protocol._indicator(sets, cnt.shape[1])
    s = np.asarray((lg @ ind).todense()) / len(idx)
    c = np.asarray((cnt @ ind).todense())
    t = np.quantile(s[:, 1:], 0.95)
    tail = (s[:, 0] > t).mean() - (s[:, 1:] > t).mean()
    m, v = c.mean(0), c.var(0, ddof=1)
    ok = (m > 0) & (v > 0)
    d = np.full(len(m), np.nan)
    d[ok] = np.log(v[ok] / m[ok])
    return tail, d[0] - np.nanmedian(d[1:])


def main():
    a = io.load()
    counts = sp.csr_matrix(a.layers["counts"])
    st = strata(a)
    controls = yaml.safe_load((ROOT / "config" / "axes.yaml").read_text())["controls"]
    fib = [g["symbol"] for g in next(c for c in controls if c["name"] == "ambient_fibroblast")["genes"]]
    veto = protocol._veto(a.obs, counts, pd.Index(a.var_names), fib)["vetoed"].to_numpy()
    thinned, low = protocol.thin(counts, protocol.T, np.random.default_rng(protocol.SEED + 1))
    lg = protocol.lognorm(thinned)
    var = pd.Index(a.var_names)
    sets = marker_sets()
    rng = np.random.default_rng(SEED)
    use = st["ventricle"] & ~veto & ~low
    heart = a.obs["heart"].astype(str).to_numpy()
    cond = a.obs["condition"].astype(str).to_numpy()

    rows = []
    for h in sorted(set(heart[use])):
        m = np.flatnonzero(use & (heart == h))
        if len(m) < 30:
            continue
        cnt, l = sp.csr_matrix(thinned[m]), sp.csr_matrix(lg[m])
        for name in PRIMARY + ["dissociation_control"]:
            idx = np.array([var.get_loc(g) for g in sets[name] if g in var])
            tail, disp = shape(cnt, l, idx, rng)
            rows.append({"heart": h, "condition": cond[m][0], "n_cells": len(m), "set": name, "tail_excess": tail, "disp_excess": disp})
    per = pd.DataFrame(rows)
    per.to_csv(OUT / "cell_tail_per_heart.csv", index=False)

    res = []
    for name, g in per.groupby("set"):
        meta = g.drop_duplicates("heart").set_index("heart")[["condition"]].assign(heart=lambda d: d.index)
        for measure in ["tail_excess", "disp_excess"]:
            score = g.set_index("heart")[measure]
            if score.isna().any():  # set genes undetected in a heart: the measure is not computable
                res.append({"set": name, "primary": name in PRIMARY, "measure": measure, "contrast": "not computable",
                            "n_case": int(score.isna().sum())})
                continue
            for cname, (case, ctrl, only) in CONTRASTS.items():
                res.append({"set": name, "primary": name in PRIMARY, "measure": measure, "contrast": cname,
                            **{k: v for k, v in test(score, meta, case, ctrl, only).items() if k != "loo_flag"}})
    res = pd.DataFrame(res)
    for (measure, cname), g in res[res["primary"] & res["contrast"].isin(STEP)].groupby(["measure", "contrast"]):
        res.loc[g.index, "q_bh"] = false_discovery_control(g["p"].to_numpy())
    res.to_csv(OUT / "cell_tail.csv", index=False)

    item18 = pd.read_csv(OUT / "pseudobulk_set_plan.csv")
    item18 = item18[item18["run"] == "raw"].set_index(["set", "contrast"])
    calls = []
    for name in PRIMARY:
        out = {"set": name}
        for cname in STEP:
            t = res[(res["set"] == name) & (res["contrast"] == cname)].set_index("measure")
            if len(t) < 2:
                out[cname] = "cannot tell at 600 UMIs (not computable)"
                continue
            sub = t.loc["tail_excess", "q_bh"] <= ALPHA and t.loc["disp_excess", "p"] <= ALPHA
            shifted = item18.loc[(name, cname), "q_bh"] <= ALPHA
            out[cname] = "a subgroup switches on" if sub else ("every cell shifts" if shifted and t.loc["disp_excess", "p"] > ALPHA else "cannot tell at 600 UMIs")
        calls.append(out)
    calls = pd.DataFrame(calls)
    calls.to_csv(OUT / "cell_tail_calls.csv", index=False)
    with pd.option_context("display.width", 220, "display.max_columns", 20):
        print(res[res["contrast"].isin(["7v0", "3v0", "3v30"])][["set", "measure", "contrast", "n_case", "n_ctrl", "est", "lo", "hi", "p", "q_bh"]].round(3).to_string(index=False))
        print(calls.to_string(index=False))


if __name__ == "__main__":
    main()
