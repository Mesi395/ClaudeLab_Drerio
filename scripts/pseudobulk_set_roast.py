"""Section A re-run with the strategy section G chose (docs/statistical-plan.md, G1 decision).

G1 (scripts/strategy_benchmark.py, 9fe5edd) chose M6, the rotation gene-set test on pseudobulk log
expression, over M4 (exact relabelling of the per-heart z-score): both calibrated, M6 had 7 points
more mean shift power at 3 vs 4 hearts. The plan says: if the choice is not M4, section A is re-run
with it. This script swaps only the p-value: same pseudobulks, sets, contrasts, BH families and
qualifiers as scripts/pseudobulk_set_plan.py (estimates, leave-one-heart-out and the handling /
ambient correlations are read from its committed output). The rotation test is run as benchmarked
(strategy_benchmark.roast, 999 rotations, the contrast's hearts only).

Usage: PYTHONPATH=src:scripts ZF_DATA=<data> .venv/bin/python scripts/pseudobulk_set_roast.py
"""

import numpy as np
import pandas as pd
import scipy.sparse as sp
import yaml
from pseudobulk_de import MIN_CELLS, OUT, pseudobulk
from pseudobulk_set_plan import ALPHA, CONTRASTS, R_MAX, SMALL, STEP, groups, plan_sets
from pseudobulk_set_score import PRIMARY
from run_protocol import ROOT, strata
from scipy.stats import false_discovery_control
from strategy_benchmark import roast

from zf import io, protocol

SEED, N_RANDOM, RANDOM_SIZE = 20261110, 200, 10


def main():
    a = io.load()
    counts = sp.csr_matrix(a.layers["counts"])
    st = strata(a)
    cham = np.where(st["ventricle"], "ventricle", np.where(st["atrium"], "atrium", "uncertain"))
    controls = yaml.safe_load((ROOT / "config" / "axes.yaml").read_text())["controls"]
    fib = [
        g["symbol"] for g in next(c for c in controls if c["name"] == "ambient_fibroblast")["genes"]
    ]
    var = pd.Index(a.var_names)
    keep = pd.Series(
        ~protocol._veto(a.obs, counts, var, fib)["vetoed"].to_numpy(), index=a.obs_names
    )
    thinned, _ = protocol.thin(counts, protocol.T, np.random.default_rng(protocol.SEED + 1))
    sets = plan_sets(fib)
    rng = np.random.default_rng(SEED)

    rows, rand_rows = [], []
    for run, cm in {"raw": counts, "thinned600": sp.csr_matrix(thinned.astype(np.int64))}.items():
        xa, ma = pseudobulk(a, cm, cham, keep)
        v = (ma["n_cells"] >= MIN_CELLS) & ma["chamber"].eq("ventricle")
        x, meta = xa[v], ma[v]
        pool = x.columns[(x > 0).all(axis=0)]
        rand = [list(rng.choice(pool, RANDOM_SIZE, replace=False)) for _ in range(N_RANDOM)]
        for cname, (case, ctrl, only) in CONTRASTS.items():
            ic, ik = groups(meta, case, ctrl, only)
            xc = x.loc[ic.append(ik)]
            small = xc.index[meta.loc[xc.index, "n_cells"] >= SMALL]
            for name, genes in sets.items():
                g = [k for k in genes if k in var]
                row = {
                    "run": run,
                    "set": name,
                    "primary": name in PRIMARY,
                    "contrast": cname,
                    "p_roast": roast(xc, ic, g, rng),
                }
                if cname in STEP and run == "raw":
                    row["p_roast_small_out"] = roast(xc.loc[small], ic.intersection(small), g, rng)
                rows.append(row)
            ps = np.array([roast(xc, ic, g, rng) for g in rand])
            rand_rows.append(
                {
                    "run": run,
                    "contrast": cname,
                    "n_sets": len(ps),
                    "share_p05": float(np.mean(ps <= ALPHA)),
                }
            )

    res = pd.DataFrame(rows)
    for _, g in res[res["primary"] & res["contrast"].isin(STEP)].groupby(["run", "contrast"]):
        res.loc[g.index, "q_roast"] = false_discovery_control(g["p_roast"].to_numpy())
    old = pd.read_csv(OUT / "pseudobulk_set_plan.csv")
    res = res.merge(
        old[
            [
                "run",
                "set",
                "contrast",
                "est",
                "lo",
                "hi",
                "p",
                "q_bh",
                "loo_flag",
                "r_handling",
                "r_ambient",
            ]
        ].rename(columns={"p": "p_exact", "q_bh": "q_exact"}),
        on=["run", "set", "contrast"],
        how="left",
    )
    rand = pd.DataFrame(rand_rows)

    calls = []
    for name in PRIMARY:
        g = res[(res["run"] == "raw") & (res["set"] == name)].set_index("contrast")
        th = res[(res["run"] == "thinned600") & (res["set"] == name)].set_index("contrast")
        out = {"set": name}
        for cname in STEP:
            r = g.loc[cname]
            if r["q_roast"] > ALPHA:
                out[cname] = "not shown"
                continue
            q = []
            if r["r_handling"] >= R_MAX:
                q.append("cannot be separated from handling")
            if r["r_ambient"] >= R_MAX:
                q.append("cannot be separated from ambient RNA")
            if r["loo_flag"]:
                q.append("rests on one heart")
            if th.loc[cname, "p_roast"] > ALPHA:
                q.append("depth-sensitive")
            if r["p_roast_small_out"] > ALPHA:
                q.append("small-heart-sensitive")
            kind = (
                "cardiomyocyte programme" if not any(k.startswith("cannot") for k in q) else "shift"
            )
            out[cname] = "; ".join([kind] + q) + (" [3 dpi confounded]" if cname == "3v0" else "")
        calls.append(out)
    calls = pd.DataFrame(calls)

    res.to_csv(OUT / "pseudobulk_set_roast.csv", index=False)
    rand.to_csv(OUT / "pseudobulk_set_roast_random.csv", index=False)
    calls.to_csv(OUT / "pseudobulk_set_roast_calls.csv", index=False)
    cols = [
        "set",
        "contrast",
        "est",
        "p_exact",
        "q_exact",
        "p_roast",
        "q_roast",
        "p_roast_small_out",
    ]
    with pd.option_context(
        "display.width", 250, "display.max_columns", 30, "display.max_colwidth", 90
    ):
        print(res[(res["run"] == "raw")][cols].round(4).to_string(index=False))
        print(
            res[(res["run"] == "thinned600") & res["primary"]][cols[:7]]
            .round(4)
            .to_string(index=False)
        )
        print(rand.round(3).to_string(index=False))
        print(calls.to_string(index=False))


if __name__ == "__main__":
    main()
