"""Statistical plan v2, question A (docs/statistical-plan.md sections 2A, 3, 5).

Pre-registered: committed together with the plan before its first run.

Heart-level test of each gene set: one score per ventricle heart (set_scores of
pseudobulk_set_score.py), Welch t, one-sided, p by exact relabelling of hearts. 7 dpi vs control is
the primary family (BH over six sets), 3 dpi vs control the second (confounded by chemistry and
series; the Hr11+Hr12 estimate is shown beside it), 30 dpi and 3 vs 30 dpi are estimates. Gene
splits, handling / ambient correlations, leave-one-heart-out, thinning, small hearts out and a
ventricle-minus-atrium line qualify the calls.

Usage: PYTHONPATH=src:scripts ZF_DATA=<data> .venv/bin/python scripts/pseudobulk_set_plan.py
"""

import numpy as np
import pandas as pd
import scipy.sparse as sp
import yaml
from pseudobulk_de import MIN_CELLS, OUT, pseudobulk
from pseudobulk_set_score import PRIMARY, exact_perm, set_scores
from pseudobulk_set_test import marker_sets
from run_protocol import ROOT, strata
from scipy import stats
from scipy.stats import false_discovery_control

from zf import io, protocol

GLYCOLYTIC = {"hk1", "pkma", "pdk2a", "ldha", "pgk1", "eno1a", "tpi1b", "pgam2", "gapdh"}
LIBSIZE = {"ttn.1", "ttn.2", "bves", "synpo2lb"}
V2_3DPI = ["Hr11", "Hr12"]
N_RANDOM, RANDOM_SIZE, SMALL, R_MAX, ALPHA, SEED = 200, 10, 100, 0.6, 0.05, 20261107
CONTRASTS = {  # name: (case condition, control condition, case hearts restricted to or None)
    "7v0": ("dpi_7", "uninjured", None),
    "3v0": ("dpi_3", "uninjured", None),
    "3v0_v2only": ("dpi_3", "uninjured", V2_3DPI),
    "30v0": ("dpi_30", "uninjured", None),
    "3v30": ("dpi_3", "dpi_30", None),
}
STEP = {"7v0": 1, "3v0": 2}


def plan_sets(fib):
    s = marker_sets()
    s["dedifferentiated_no_glycolysis"] = [g for g in s["dedifferentiated"] if g not in GLYCOLYTIC]
    s["dedifferentiated_no_glycolysis_no_libsize"] = [
        g for g in s["dedifferentiated_no_glycolysis"] if g not in LIBSIZE
    ]
    s["injury_responsive_no_dissociation"] = [
        g for g in s["injury_responsive"] if g not in s["dissociation_control"]
    ]
    s["ambient_fibroblast"] = fib
    s["pdk2a"] = ["pdk2a"]
    return s


def welch_ci(a, b):
    va, vb = a.var(ddof=1) / len(a), b.var(ddof=1) / len(b)
    df = (va + vb) ** 2 / (va**2 / (len(a) - 1) + vb**2 / (len(b) - 1))
    d, h = a.mean() - b.mean(), stats.t.ppf(0.975, df) * np.sqrt(va + vb)
    return d, d - h, d + h


def groups(meta, case, ctrl, only):
    ic = meta.index[meta["condition"].eq(case) & (meta["heart"].isin(only) if only else True)]
    return ic, meta.index[meta["condition"].eq(ctrl)]


def test(score, meta, case, ctrl, only):
    ic, ik = groups(meta, case, ctrl, only)
    v = score.loc[ic.append(ik)].to_numpy()
    t, p, n = exact_perm(v, len(ic))
    d, lo, hi = welch_ci(v[: len(ic)], v[len(ic) :])
    k, loo = len(ic), []
    for i in range(len(v)):  # leave heart i out; the case group shrinks when i < k
        w, kk = np.delete(v, i), k - (i < k)
        loo.append(w[:kk].mean() - w[kk:].mean())
    flag = bool(np.any(np.sign(loo) != np.sign(d)) or np.any(np.abs(np.array(loo) - d) > abs(d) / 2))
    return {"n_case": len(ic), "n_ctrl": len(ik), "est": d, "lo": lo, "hi": hi, "welch_t": t,
            "p": p, "n_perm": n, "loo_min": min(loo), "loo_max": max(loo), "loo_flag": flag}


def within_r(s, meta, a, b):
    """Correlation of two set scores after centring each on its condition mean (all hearts)."""
    c = lambda x: x - x.groupby(meta["condition"]).transform("mean")  # noqa: E731
    return float(np.corrcoef(c(s[a]), c(s[b]))[0, 1])


def main():
    a = io.load()
    counts = sp.csr_matrix(a.layers["counts"])
    st = strata(a)
    cham = np.where(st["ventricle"], "ventricle", np.where(st["atrium"], "atrium", "uncertain"))
    controls = yaml.safe_load((ROOT / "config" / "axes.yaml").read_text())["controls"]
    fib = [g["symbol"] for g in next(c for c in controls if c["name"] == "ambient_fibroblast")["genes"]]
    keep = pd.Series(~protocol._veto(a.obs, counts, pd.Index(a.var_names), fib)["vetoed"].to_numpy(), index=a.obs_names)
    thinned, _ = protocol.thin(counts, protocol.T, np.random.default_rng(protocol.SEED + 1))
    sets = plan_sets(fib)
    named = list(sets)
    rng = np.random.default_rng(SEED)

    rows, rand_rows, cham_rows = [], [], []
    for run, cm in {"raw": counts, "thinned600": sp.csr_matrix(thinned.astype(np.int64))}.items():
        xa, ma = pseudobulk(a, cm, cham, keep)
        big = ma["n_cells"] >= MIN_CELLS
        xa, ma = xa[big], ma[big]
        v = ma["chamber"].eq("ventricle")
        x, meta = xa[v], ma[v]
        pool = x.columns[(x > 0).all(axis=0)]
        rand = {f"random_{i}": list(rng.choice(pool, RANDOM_SIZE, replace=False)) for i in range(N_RANDOM)}
        s = set_scores(x, {**sets, **rand})
        s.join(meta[["heart", "condition", "n_cells"]]).assign(run=run).to_csv(OUT / f"pseudobulk_set_plan_scores_{run}.csv")
        small_out = meta["n_cells"] >= SMALL
        for name in named:
            r_h = within_r(s, meta, name, "dissociation_control")
            r_a = within_r(s, meta, name, "ambient_fibroblast")
            for cname, (case, ctrl, only) in CONTRASTS.items():
                row = {"run": run, "set": name, "primary": name in PRIMARY, "contrast": cname,
                       **test(s[name], meta, case, ctrl, only), "r_handling": r_h, "r_ambient": r_a}
                if cname in STEP:
                    row["p_small_out"] = test(s.loc[small_out, name], meta[small_out], case, ctrl, only)["p"]
                rows.append(row)
        for cname, (case, ctrl, only) in CONTRASTS.items():
            ps = np.array([test(s[k], meta, case, ctrl, only)["p"] for k in rand])
            rand_rows.append({"run": run, "contrast": cname, "n_sets": len(ps), "share_p05": float(np.mean(ps <= ALPHA))})
        # Chamber line: score over all chamber pseudobulks, ventricle minus atrium per heart.
        sa = set_scores(xa, sets)
        both = ma.groupby("heart")["chamber"].nunique().pipe(lambda c: c[c == 2].index)
        hv = ma[v & ma["heart"].isin(both)].reset_index().set_index("heart")
        ha = ma[~v & ma["heart"].isin(both)].reset_index().set_index("heart")
        diff = sa.loc[hv["sample"]].set_axis(hv.index) - sa.loc[ha.loc[hv.index, "sample"]].set_axis(hv.index)
        mh = hv[["condition"]].assign(heart=hv.index)
        for name in PRIMARY + ["dedifferentiated_no_glycolysis", "ambient_fibroblast", "dissociation_control"]:
            for cname in STEP:
                case, ctrl, only = CONTRASTS[cname]
                cham_rows.append({"run": run, "set": name, "contrast": cname, **test(diff[name], mh, case, ctrl, only)})

    res = pd.DataFrame(rows)
    for (run, cname), g in res[res["primary"] & res["contrast"].isin(STEP)].groupby(["run", "contrast"]):
        res.loc[g.index, "q_bh"] = false_discovery_control(g["p"].to_numpy())
    rand = pd.DataFrame(rand_rows)
    chamber = pd.DataFrame(cham_rows)

    calls = []
    for name in PRIMARY:
        g = res[(res["run"] == "raw") & (res["set"] == name)].set_index("contrast")
        th = res[(res["run"] == "thinned600") & (res["set"] == name)].set_index("contrast")
        out = {"set": name}
        for cname in STEP:
            r = g.loc[cname]
            if r["q_bh"] > ALPHA:
                out[cname] = "not shown"
                continue
            q = []
            if r["r_handling"] >= R_MAX:
                q.append("cannot be separated from handling")
            if r["r_ambient"] >= R_MAX:
                q.append("cannot be separated from ambient RNA")
            if r["loo_flag"]:
                q.append("rests on one heart")
            if th.loc[cname, "p"] > ALPHA:
                q.append("depth-sensitive")
            if r["p_small_out"] > ALPHA:
                q.append("small-heart-sensitive")
            kind = "cardiomyocyte programme" if not any(k.startswith("cannot") for k in q) else "shift"
            out[cname] = "; ".join([kind] + q) + (" [3 dpi confounded]" if cname == "3v0" else "")
        f = g.loc["3v30"]
        out["3v30"] = f"{f['est']:+.2f} [{f['lo']:+.2f}, {f['hi']:+.2f}] p {f['p']:.3f}" + (" fades" if f["p"] <= 1 / 35 + 1e-9 else "")
        calls.append(out)
    for split, parent in [("dedifferentiated_no_glycolysis", "dedifferentiated"),
                          ("memory_bz_core_no_ap1", "memory_bz_core"),
                          ("injury_responsive_no_dissociation", "injury_responsive")]:
        g = res[(res["run"] == "raw") & (res["set"] == split)].set_index("contrast")
        calls.append({"set": f"{split} (split of {parent})", **{c: f"p {g.loc[c, 'p']:.3f}" for c in STEP}})
    calls = pd.DataFrame(calls)

    res.to_csv(OUT / "pseudobulk_set_plan.csv", index=False)
    rand.to_csv(OUT / "pseudobulk_set_plan_random.csv", index=False)
    chamber.to_csv(OUT / "pseudobulk_set_plan_chamber.csv", index=False)
    calls.to_csv(OUT / "pseudobulk_set_plan_calls.csv", index=False)
    cols = ["run", "set", "contrast", "n_case", "n_ctrl", "est", "lo", "hi", "p", "q_bh", "loo_flag", "r_handling", "r_ambient", "p_small_out"]
    with pd.option_context("display.width", 250, "display.max_columns", 30, "display.max_colwidth", 90):
        print(res[res["run"] == "raw"][cols].round(3).to_string(index=False))
        print(rand.round(3).to_string(index=False))
        print(chamber[chamber["run"] == "raw"][["set", "contrast", "n_case", "n_ctrl", "est", "lo", "hi", "p"]].round(3).to_string(index=False))
        print(calls.to_string(index=False))


if __name__ == "__main__":
    main()
