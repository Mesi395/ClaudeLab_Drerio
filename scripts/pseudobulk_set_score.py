"""Per-heart gene-set scores tested by exact relabelling of hearts (a self-contained set test).

Pre-registered 2026-11-07: this file, its sets, statistic and decision rule were committed before
it was first run.

Why: the competitive set test (pseudobulk_set_test.py) compares a set's mean Wald statistic with
decoy gene sets, so it treats genes as independent, and correlated set genes inflate it (Wu & Smyth
2012, camera). Here the unit is the heart. Each ventricle pseudobulk gets one score per set and the
hearts are relabelled, so gene-gene correlation stays intact (self-contained null, Goeman &
Buehlmann 2007). It also gives the 30 dpi control a formal test, which until now was judged by eye.

Score: DESeq2 median-of-ratios size factors over all kept ventricle pseudobulks, log2(count/sf + 1),
z-scored per gene across those samples, averaged over the set's genes present. No step uses the
labels, so the relabelling null is exact.

Tests: one-sided, Welch t as the statistic, p = share of all relabellings with t >= observed
(observed included). dpi_3 > uninjured (4 vs 4, 70 relabellings, min p 0.014); dpi_30 > uninjured
(3 vs 4, 35, min p 0.029); dpi_3 > dpi_30 (4 vs 3, 35, min p 0.029); dpi_7 > uninjured as an
extra. Effect size: Hedges' g. BH across the six PRIMARY sets within each run and contrast.

Decision rule (fixed before running; on raw permutation p, BH q reported beside it):
- fades (injury-like): 3>U p<=0.05 and 3>30 p<=0.05
- persists:            3>U p<=0.05 and 30>U p<=0.05 and 3>30 p>0.05
- up at 3 dpi only:    3>U p<=0.05 and neither of the above (unresolved)
- not shown:           3>U p>0.05
A call whose deciding tests do not all reach BH q<=0.05 is labelled nominal. With 35 relabellings
and six sets, a 30 dpi test reaches q<=0.05 only if several sets sit at the floor.
Raw counts are primary; counts thinned to 600 UMIs per cell are the sensitivity run.
memory_bz_core_no_ap1 and dissociation_control are read beside the primary sets: a pattern that
the dissociation set shares, or that vanishes without AP-1 genes, is not an injury programme.

Usage: PYTHONPATH=src:scripts ZF_DATA=<data> .venv/bin/python scripts/pseudobulk_set_score.py
"""

from itertools import combinations

import numpy as np
import pandas as pd
import scipy.sparse as sp
import yaml
from pseudobulk_de import MIN_CELLS, OUT, pseudobulk
from pseudobulk_set_test import marker_sets
from run_protocol import ROOT, strata
from scipy.stats import false_discovery_control

from zf import io, protocol

PRIMARY = [
    "memory_bz_core",
    "border_zone",
    "dedifferentiated",
    "injury_responsive",
    "proliferating",
    "glycolysis_honkoop2019",
]
CONTRASTS = [("dpi_3", "uninjured"), ("dpi_30", "uninjured"), ("dpi_3", "dpi_30"), ("dpi_7", "uninjured")]
ALPHA = 0.05


def set_scores(x, sets):
    """Per-sample score of each set: size-factor log expression, z per gene, mean over set genes."""
    pos = x.loc[:, (x > 0).all(axis=0)]
    lg = np.log(pos)
    sf = np.exp(lg.sub(lg.mean(axis=0), axis=1).median(axis=1))
    y = np.log2(x.div(sf, axis=0) + 1)
    z = (y - y.mean()) / y.std(ddof=1)  # genes with no variance become NaN and drop out of the mean
    return pd.DataFrame({n: z[[g for g in gs if g in z.columns]].mean(axis=1) for n, gs in sets.items()})


def welch(a, b):
    return (a.mean() - b.mean()) / np.sqrt(a.var(ddof=1) / len(a) + b.var(ddof=1) / len(b))


def hedges(a, b):
    na, nb = len(a), len(b)
    pooled = np.sqrt(((na - 1) * a.var(ddof=1) + (nb - 1) * b.var(ddof=1)) / (na + nb - 2))
    return (a.mean() - b.mean()) / pooled * (1 - 3 / (4 * (na + nb) - 9))


def exact_perm(v, k):
    """One-sided exact p of Welch t for the first k of v against the rest, over all relabellings."""
    obs = welch(v[:k], v[k:])
    idx = np.arange(len(v))
    null = np.array([welch(v[list(c)], v[np.setdiff1d(idx, c)]) for c in combinations(idx, k)])
    return obs, float(np.mean(null >= obs - 1e-12)), len(null)


def call(p):
    up3, up30, fade = (p[c] <= ALPHA for c in ["dpi_3>uninjured", "dpi_30>uninjured", "dpi_3>dpi_30"])
    if not up3:
        return "not shown", ["dpi_3>uninjured"]
    if fade:
        return "fades", ["dpi_3>uninjured", "dpi_3>dpi_30"]
    if up30:
        return "persists", ["dpi_3>uninjured", "dpi_30>uninjured"]
    return "up at 3 dpi only", ["dpi_3>uninjured"]


def main():
    a = io.load()
    counts = sp.csr_matrix(a.layers["counts"])
    st = strata(a)
    cham = np.where(st["ventricle"], "ventricle", np.where(st["atrium"], "atrium", "uncertain"))
    controls = yaml.safe_load((ROOT / "config" / "axes.yaml").read_text())["controls"]
    fib = [g["symbol"] for g in next(c for c in controls if c["name"] == "ambient_fibroblast")["genes"]]
    keep = pd.Series(~protocol._veto(a.obs, counts, pd.Index(a.var_names), fib)["vetoed"].to_numpy(), index=a.obs_names)
    thinned, _ = protocol.thin(counts, protocol.T, np.random.default_rng(protocol.SEED + 1))
    sets = marker_sets()

    rows, per_sample = [], []
    for run, cm in {"raw": counts, "thinned600": sp.csr_matrix(thinned.astype(np.int64))}.items():
        x, meta = pseudobulk(a, cm, cham, keep)
        v = meta["chamber"].eq("ventricle") & (meta["n_cells"] >= MIN_CELLS)
        x, meta = x[v], meta[v]
        s = set_scores(x, sets)
        per_sample.append(s.join(meta[["heart", "condition", "n_cells"]]).assign(run=run))
        for case, ctrl in CONTRASTS:
            ic, ik = meta.index[meta["condition"].eq(case)], meta.index[meta["condition"].eq(ctrl)]
            for name in sets:
                vals = s.loc[ic.append(ik), name].to_numpy()
                t, p, n = exact_perm(vals, len(ic))
                rows.append(
                    {
                        "run": run,
                        "set": name,
                        "primary": name in PRIMARY,
                        "contrast": f"{case}>{ctrl}",
                        "n_case": len(ic),
                        "n_ctrl": len(ik),
                        "n_genes": int(s[name].notna().any() and sum(g in x.columns for g in sets[name])),
                        "mean_case": vals[: len(ic)].mean(),
                        "mean_ctrl": vals[len(ic) :].mean(),
                        "welch_t": t,
                        "hedges_g": hedges(vals[: len(ic)], vals[len(ic) :]),
                        "p_perm": p,
                        "n_perm": n,
                    }
                )
    res = pd.DataFrame(rows)
    res["q_bh"] = np.nan
    prim = res["primary"]
    for _, g in res[prim].groupby(["run", "contrast"]):
        res.loc[g.index, "q_bh"] = false_discovery_control(g["p_perm"].to_numpy())

    calls = []
    for (run, name), g in res.groupby(["run", "set"]):
        p, q = g.set_index("contrast")["p_perm"], g.set_index("contrast")["q_bh"]
        c, used = call(p)
        nominal = name in PRIMARY and c != "not shown" and not (q[used] <= ALPHA).all()
        calls.append({"run": run, "set": name, "primary": name in PRIMARY, "call": c, "nominal": nominal, **{f"p[{k}]": p[k] for k in p.index}})
    calls = pd.DataFrame(calls)

    res.to_csv(OUT / "pseudobulk_set_score.csv", index=False)
    calls.to_csv(OUT / "pseudobulk_set_score_calls.csv", index=False)
    pd.concat(per_sample).to_csv(OUT / "pseudobulk_set_score_samples.csv")
    with pd.option_context("display.width", 200, "display.max_columns", 20):
        print(res.round(3).to_string(index=False))
        print(calls.round(3).to_string(index=False))


if __name__ == "__main__":
    main()
