"""Statistical plan v2, section G: choose the strategy by testing the candidates on this data.

Pre-registered: committed with plan section G before its first full run. A two-split smoke run
(BENCH_SPLITS=2) checked only that the code runs; it writes and prints no rates.

Null: the seven 7 dpi ventricle hearts split 3 vs 4 in all 35 ways. Real heart-to-heart variation,
no condition difference. Planted signal by binomial thinning (Gerard 2020): the set genes' counts
are thinned to 1/f in every cell except the "case" cells that carry the signal, so those carry an
f-fold higher set expression. Shift: every case cell (f = 1.25, 1.5, 2). Subgroup: a share pi of the
case cells (pi = 0.1, 0.2; f = 3).

G1, does a state programme shift between groups of hearts? Candidates, one-sided, alpha 0.05:
  M1 cell Wilcoxon: per-cell set score (mean log-normalised count), case cells vs control cells.
  M2 mixed model: per-cell score ~ group + (1 | heart), Wald z (statsmodels MixedLM, REML).
  M3 heart t: per-heart score (plan 2A, set_scores), Welch t with its t-distribution p.
  M4 heart exact: the same statistic, p by exact relabelling (35 relabellings, floor 0.029).
  M5 competitive: mean per-gene Welch t of the set vs 500 expression-matched decoy gene sets.
  M6 rotation: gene-level ROAST-style test on pseudobulk log expression (Wu et al. 2010), mean gene
     t, 999 rotations; unweighted (limma fry with voomByGroup weights needs R, absent on both
     hosts).
  M7 PLAGE: per-heart score = first singular vector of the set's per-gene z matrix (Tomfohr et al.
     2005), exact relabelling.
  M8 singscore: per-heart mean rank of the set genes within that heart, scaled (Foroutan et al.
     2018), exact relabelling.
Candidates come from the literature survey (research/s4_m1, s4_m2): M7 and M8 are the per-sample
scores it names as alternatives to the z-score; M2 is the mixed-model class it says is poorly
supported at 4-7 hearts, kept to test that on this data. False-positive rate: share of null tests
with p <= 0.05 over 35 splits x (6 primary + 50 random 10-gene sets), Clopper-Pearson 95% interval
(splits share hearts, so the interval is too narrow). Power: share of 35 splits x 6 primary sets
detected, per planted scenario. Decision: a candidate is calibrated if its false-positive interval
reaches 0.05 or lower. The strategy is the calibrated candidate with the highest mean power over the
three shift scenarios; within two percentage points the one with fewer distributional assumptions
wins (order M4, M8, M7, M3, M6, M5, M2, M1: exact relabelling of a simple score first, model-based
tests last).

G2, every cell or a subgroup? Plan 2F's shape test (tail excess and dispersion excess against
within-heart decoys, cells thinned to 600 UMIs) on the same splits and plants, raw exact p per set.
Reported per set and scenario: share of splits called "a subgroup switches on" (both p <= 0.05).
Decision (fixed, after the methods survey s4): plan 2F's "every cell shifts" is an absence-of-
evidence call, so it stands for a set only if that set's subgroup-call rate is >= 0.8 under both
subgroup plants and <= 0.1 under the null and the shift plants; otherwise the set's 2F call becomes
"cannot tell at 600 UMIs".

Usage: PYTHONPATH=src:scripts ZF_DATA=<data> .venv/bin/python scripts/strategy_benchmark.py
"""

import os
import warnings
from itertools import combinations
from multiprocessing import get_context

import numpy as np
import pandas as pd
import scipy.sparse as sp
import statsmodels.formula.api as smf
import yaml
from cell_tail_test import shape
from pseudobulk_de import OUT
from pseudobulk_set_score import PRIMARY, exact_perm, set_scores
from pseudobulk_set_test import marker_sets
from run_protocol import ROOT, strata
from scipy import stats

from zf import io, protocol

SEED, ALPHA, N_RANDOM, N_DECOY, N_ROT = 20261108, 0.05, 50, 500, 999
SCENARIOS = {
    "null": None,
    "shift x1.25": (1.0, 1.25),
    "shift x1.5": (1.0, 1.5),
    "shift x2": (1.0, 2.0),
    "subgroup 10% x3": (0.1, 3.0),
    "subgroup 20% x3": (0.2, 3.0),
}
SHAPE_SCENARIOS = ["null", "shift x1.5", "shift x2", "subgroup 10% x3", "subgroup 20% x3"]
ORDER = [
    "M4 heart exact",
    "M8 singscore exact",
    "M7 PLAGE exact",
    "M3 heart t",
    "M6 rotation (ROAST)",
    "M5 competitive",
    "M2 mixed model",
    "M1 cell Wilcoxon",
]
D = {}  # filled once in main, shared with forked workers


def plant(sub, carrier, scen, rng):
    """Every cell except the signal carriers keeps 1/f of its set counts (Gerard 2020 thinning)."""
    if scen is None:
        return sub
    pi, f = scen
    full = carrier & (rng.random(len(sub)) < pi)
    out = sub.copy()
    out[~full] = rng.binomial(sub[~full], 1 / f)
    return out


def mixed(score, case, heart):
    cells = pd.DataFrame({"s": score, "g": case.astype(float), "h": heart})
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        r = smf.mixedlm("s ~ g", cells, groups=cells["h"]).fit(reml=True)
    ok = bool(r.converged) and not any("converge" in str(x.message).lower() for x in w)
    return float(stats.norm.sf(r.params["g"] / r.bse["g"])), ok


def logexpr(x):
    """log2 expression after median-of-ratios size factors, as in set_scores."""
    lg = np.log(x.loc[:, (x > 0).all(axis=0)])
    sf = np.exp(lg.sub(lg.mean(axis=0), axis=1).median(axis=1))
    return np.log2(x.div(sf, axis=0) + 1)


def roast(x, case, genes, rng):
    """Rotation gene-set test (Wu et al. 2010): mean gene t, null by rotating the residual space."""
    y = logexpr(x)[genes]
    c = y.index.isin(case).astype(float)
    q = np.linalg.qr(np.column_stack([np.ones(len(c)), c]), mode="complete")[0]
    w = q.T[1:] @ y.to_numpy()  # row 0: group effect, rows 1..d: residuals (d = hearts - 2)
    w[0] *= np.sign(q[:, 1] @ (c - c.mean()))
    ss, d = (w**2).sum(0), len(c) - 2
    w, ss = w[:, ss > 0], ss[ss > 0]
    if not ss.size:  # no set gene varies across these hearts: no test (the floor p would be false)
        return np.nan

    def stat(b):
        return (b / np.sqrt((ss - b**2) / d)).mean(-1)

    r = rng.standard_normal((N_ROT, d + 1))
    null = stat((r / np.linalg.norm(r, axis=1, keepdims=True)) @ w)
    return (1 + np.sum(null >= stat(w[0]))) / (1 + N_ROT)


def plage(y):
    """First singular vector of the set's per-gene z matrix, signed to follow the mean z."""
    z = ((y - y.mean()) / y.std(ddof=1)).dropna(axis=1)
    u, s, _ = np.linalg.svd(z.to_numpy(), full_matrices=False)
    score = u[:, 0] * s[0]
    return score * np.sign(np.corrcoef(score, z.mean(axis=1))[0, 1])


def singscore(y, genes):
    """Mean within-heart rank of the set genes, scaled to [-0.5, 0.5] (unidirectional singscore)."""
    r = y.rank(axis=1).loc[:, genes].mean(axis=1).to_numpy()
    k, g = len(genes), y.shape[1]
    lo, hi = (k + 1) / 2, (2 * g - k + 1) / 2
    return (r - lo) / (hi - lo) - 0.5


def competitive(x, case, genes, rng):
    y = logexpr(x)
    c = y.index.isin(case)
    a, b = y[c], y[~c]
    t = (
        (a.mean() - b.mean()) / np.sqrt(a.var(ddof=1) / len(a) + b.var(ddof=1) / len(b))
    ).to_numpy()
    ok = np.isfinite(t)
    idx = np.flatnonzero(y.columns.isin(genes) & ok)
    dec = protocol.decoy_sets(np.where(ok, y.mean().to_numpy(), 0), idx, rng, b=N_DECOY)
    t = np.nan_to_num(t)
    return (1 + np.sum(t[dec].mean(1) >= t[idx].mean())) / (1 + N_DECOY)


def heart_p(v):
    """Welch t-distribution p and exact relabelling p, case hearts first (3 of 7)."""
    if np.isnan(v).any():
        return np.nan, np.nan
    return stats.ttest_ind(v[:3], v[3:], equal_var=False, alternative="greater").pvalue, exact_perm(
        v, 3
    )[1]


def g1(si, case_h, rng):
    case = np.isin(D["heart"], case_h)
    hearts = list(case_h) + [h for h in D["x"].index if h not in case_h]
    rows = []
    for scen_name, scen in SCENARIOS.items():
        for name in D["sets"] if scen is None else PRIMARY:
            genes = D["sets"][name]
            sub = plant(D["sub"][name], case, scen, rng)
            total = D["total"] - (D["sub"][name] - sub).sum(1)
            cs = np.log1p(sub * (1e4 / np.maximum(total, 1))[:, None]).mean(1)
            x = D["x"].copy()
            x[genes] = pd.DataFrame(sub, columns=genes).groupby(D["heart"]).sum().loc[x.index]
            p3, p4 = heart_p(set_scores(x, {name: genes})[name].loc[hearts].to_numpy())
            y = logexpr(x).loc[hearts]
            p7 = exact_perm(plage(y[genes]), 3)[1]
            p8 = exact_perm(singscore(y, genes), 3)[1]
            p2, ok2 = mixed(cs, case, D["heart"])
            ps = {
                "M1 cell Wilcoxon": (
                    stats.mannwhitneyu(cs[case], cs[~case], alternative="greater").pvalue,
                    True,
                ),
                "M2 mixed model": (p2, ok2),
                "M3 heart t": (p3, True),
                "M4 heart exact": (p4, True),
                "M5 competitive": (competitive(x, case_h, genes, rng), True),
                "M6 rotation (ROAST)": (roast(x, case_h, genes, rng), True),
                "M7 PLAGE exact": (p7, True),
                "M8 singscore exact": (p8, True),
            }
            rows += [
                {
                    "split": si,
                    "scenario": scen_name,
                    "set": name,
                    "random": name not in PRIMARY,
                    "method": m,
                    "p": p,
                    "converged": ok,
                }
                for m, (p, ok) in ps.items()
            ]
    return rows


def g2(si, case_h, rng):
    tcase = np.isin(D["theart"], case_h)
    hearts = list(case_h) + [h for h in D["x"].index if h not in case_h]
    rows = []
    for scen_name in SHAPE_SCENARIOS:
        for name in PRIMARY:
            idx = D["idx"][name]
            new = plant(D["tsub"][name], tcase, SCENARIOS[scen_name], rng)
            delta = (new - D["tsub"][name]).astype(float)
            r, c = np.nonzero(delta)
            cnt = (
                D["thin"] + sp.csr_matrix((delta[r, c], (r, idx[c])), shape=D["thin"].shape)
            ).tocsr()
            cnt.eliminate_zeros()
            per = {}
            for h in hearts:
                m = np.flatnonzero(D["theart"] == h)
                ch = cnt[m]
                per[h] = shape(ch, protocol.lognorm(ch), idx, rng)
            pt = heart_p(np.array([per[h][0] for h in hearts]))[1]
            pd_ = heart_p(np.array([per[h][1] for h in hearts]))[1]
            rows.append(
                {
                    "split": si,
                    "scenario": scen_name,
                    "set": name,
                    "p_tail": pt,
                    "p_disp": pd_,
                    "subgroup_call": bool(pt <= ALPHA and pd_ <= ALPHA),
                }
            )
    return rows


def run(job):
    si, case_h = job
    rng = np.random.default_rng([SEED, si])
    return g1(si, case_h, rng), g2(si, case_h, rng)


def rate(hit, n):
    k = int(hit)
    ci = stats.binomtest(k, n).proportion_ci(0.95, method="exact") if n else None
    return {
        "n": n,
        "k": k,
        "rate": k / n if n else np.nan,
        "lo": ci.low if n else np.nan,
        "hi": ci.high if n else np.nan,
    }


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
    use = st["ventricle"] & ~veto & (a.obs["condition"].astype(str).to_numpy() == "dpi_7")
    var = pd.Index(a.var_names)
    heart = a.obs["heart"].astype(str).to_numpy()
    cnt = counts[np.flatnonzero(use)]
    D["heart"], D["total"] = heart[use], np.asarray(cnt.sum(1)).ravel()
    hs = sorted(set(D["heart"]))
    assert len(hs) == 7, hs
    ind = sp.csr_matrix(
        (
            np.ones(len(D["heart"])),
            (pd.Index(hs).get_indexer(D["heart"]), np.arange(len(D["heart"]))),
        )
    )
    D["x"] = pd.DataFrame((ind @ cnt).toarray().astype(int), index=hs, columns=var)
    ms = marker_sets()
    sets = {n: [g for g in ms[n] if g in var] for n in PRIMARY}
    pool = var[(D["x"] > 0).all(axis=0).to_numpy()]
    rng = np.random.default_rng(SEED)
    sets |= {f"random_{i:02d}": list(rng.choice(pool, 10, replace=False)) for i in range(N_RANDOM)}
    D["sets"] = sets
    D["sub"] = {n: cnt[:, var.get_indexer(g)].toarray().astype(np.int64) for n, g in sets.items()}
    tm = np.flatnonzero(use & ~low)
    D["thin"] = sp.csr_matrix(thinned[tm])
    D["theart"] = heart[tm]
    D["idx"] = {n: var.get_indexer(sets[n]) for n in PRIMARY}
    D["tsub"] = {n: D["thin"][:, D["idx"][n]].toarray().astype(np.int64) for n in PRIMARY}

    jobs = list(enumerate(combinations(hs, 3)))[: int(os.environ.get("BENCH_SPLITS", "35"))]
    with get_context("fork").Pool(min(16, len(jobs))) as p:  # fork: workers share D without copying
        out = p.map(run, jobs)
    tests = pd.DataFrame([r for o in out for r in o[0]])
    shapes = pd.DataFrame([r for o in out for r in o[1]])
    print(
        f"tests {len(tests)}, shape tests {len(shapes)}, "
        f"mixed model not converged {int((~tests['converged']).sum())}, "
        f"NaN p {int(tests['p'].isna().sum())}, "
        f"NaN shape p {int(shapes[['p_tail', 'p_disp']].isna().any(axis=1).sum())}"
    )
    if len(jobs) < 35:
        return  # smoke run: no rates

    summ = []
    for (m, s), g in tests.groupby(["method", "scenario"]):
        g = g if s == "null" else g[~g["random"]]
        summ.append(
            {
                "method": m,
                "scenario": s,
                **rate(((g["p"] <= ALPHA) & g["converged"]).sum(), len(g)),
                "not_converged": int((~g["converged"]).sum()),
            }
        )  # a fit that did not converge counts as not detected
    summ = pd.DataFrame(summ)
    fpr = summ[summ["scenario"] == "null"].set_index("method")
    power = summ[summ["scenario"].str.startswith("shift")].groupby("method")["rate"].mean()
    cal = [m for m in ORDER if fpr.loc[m, "lo"] <= ALPHA]
    best = max(power[cal]) if cal else np.nan
    choice = next((m for m in ORDER if m in cal and power[m] >= best - 0.02), "none calibrated")
    summ["calibrated"] = summ["method"].isin(cal)
    summ["mean_shift_power"] = summ["method"].map(power)
    summ["chosen"] = summ["method"].eq(choice)
    sh = pd.DataFrame(
        [
            {
                "set": n,
                "scenario": s,
                **rate(g["subgroup_call"].sum(), len(g)),
                "tail_rate": (g["p_tail"] <= ALPHA).mean(),
                "disp_rate": (g["p_disp"] <= ALPHA).mean(),
            }
            for (n, s), g in shapes.groupby(["set", "scenario"])
        ]
    )
    r = sh.pivot_table(index="set", columns="scenario", values="rate")
    fit = (r[[s for s in r if s.startswith("subgroup")]].min(axis=1) >= 0.8) & (
        r[["null"] + [s for s in r if s.startswith("shift")]].max(axis=1) <= 0.1
    )
    sh["fit_for_2F"] = sh["set"].map(fit)
    tests.to_csv(OUT / "strategy_benchmark_tests.csv", index=False)
    summ.to_csv(OUT / "strategy_benchmark.csv", index=False)
    shapes.to_csv(OUT / "strategy_benchmark_shape_tests.csv", index=False)
    sh.to_csv(OUT / "strategy_benchmark_shape.csv", index=False)
    with pd.option_context("display.width", 200, "display.max_columns", 20):
        print(summ.round(3).to_string(index=False))
        print(sh.round(3).to_string(index=False))
    print(f"chosen strategy: {choice}")
    print(f"sets where 2F can tell every-cell from subgroup: {sorted(fit[fit].index)}")


if __name__ == "__main__":
    main()
