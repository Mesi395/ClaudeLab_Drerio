"""C004: chamber call + validation, agent label file, ventricle-only positive control, G4 gate.

Rules: research/cm-state-axes/cycles/C004/reading-rules.md (fixed 15:43 CET, before computing).
uv run python scripts/c004_probe.py data/clean/hu2022_cm_allgenes_sanitized.h5ad markers results
"""

import json
import re
import subprocess
import sys
from datetime import datetime
from itertools import combinations, product
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import yaml
from scipy.stats import spearmanr

from zf.io import ALLOWED_OBS
from zf.score import matched_random_sets

TARGET, SEED, N_NULL, CIRCULAR, MIN_CELLS = 500, 0, 100, {"ttn.2"}, 50
HEARTS = {
    "Ctrl": ["H5", "H6", "H7", "H8"],
    "3dpi": ["Hr11", "Hr12", "Hr24", "Hr27"],
    "7dpi": ["Hr1", "Hr2", "Hr6", "Hr7", "Hr8", "Hr9", "Hr13"],
    "30dpi": ["Hr3", "Hr4", "Hr19"],
}
PRED = {
    "border_zone_like": ["7dpi"],
    "proliferating": ["7dpi"],
    "injury_responsive": ["3dpi"],
    "dedifferentiated": ["3dpi", "7dpi"],
}
CHAMBER_LIBS = {
    "H8a": "atrium",
    "Hr6a": "atrium",
    "Hr7a": "atrium",
    "H8v": "ventricle",
    "Hr6v": "ventricle",
    "Hr7v": "ventricle",
}


def heart(lib):
    return re.sub(r"^(Hr2|H8|Hr6|Hr7)[abv]$", r"\1", lib)


def auc(inj, ctl):
    return float(np.mean([1.0 if i > c else 0.5 if i == c else 0.0 for i, c in product(inj, ctl)]))


def score_all(b, axes, vn, means, all_idx):
    """score_genes for each axis and its 100 matched null sets on AnnData b (log1p CP10k)."""
    vidx = {g: i for i, g in enumerate(vn)}
    real, null = {}, {}
    for k, (st, genes) in enumerate(axes.items()):
        gi = np.array([vidx[g] for g in genes])
        sc.tl.score_genes(b, genes, random_state=SEED, score_name="_s", use_raw=False)
        real[st] = b.obs.pop("_s").to_numpy()
        sets = matched_random_sets(means, gi, N_NULL, 25, SEED + k, exclude=all_idx)
        cols = []
        for j, s in enumerate(sets):
            sc.tl.score_genes(
                b, [vn[i] for i in s], random_state=SEED + j, score_name="_s", use_raw=False
            )
            cols.append(b.obs.pop("_s").to_numpy())
        null[st] = np.stack(cols, 1)
    return real, null


def main(src, marker_dir, out_dir):
    out = Path(out_dir)
    a = ad.read_h5ad(src)
    assert set(a.obs.columns) == set(ALLOWED_OBS)
    counts = a.layers["counts"].tocsr().astype(np.float32)
    tot = np.asarray(counts.sum(1)).ravel()
    lib = a.obs["orig.ident"].astype(str).to_numpy()
    hrt = np.array([heart(x) for x in lib])
    time_of = {h: t for t, hs in HEARTS.items() for h in hs}
    vn = list(a.var_names)
    vidx = {g: i for i, g in enumerate(vn)}

    # ---- 1. chamber call ----
    def col(g):
        return counts[:, vidx[g]].toarray().ravel()

    atr = col("myh6")
    calls = {"A": col("myh7"), "B": col("myh7") + col("myh7l")}
    val, chamber = [], {}
    for rule, ven in calls.items():
        c = np.where(atr > ven, "atrium", np.where(ven > atr, "ventricle", "unassigned"))
        chamber[rule] = c
        for lb, truth in CHAMBER_LIBS.items():
            m = lib == lb
            assigned = c[m] != "unassigned"
            val.append(
                {
                    "rule": rule,
                    "lib": lb,
                    "truth": truth,
                    "n_cells": int(m.sum()),
                    "coverage": float(assigned.mean()),
                    "accuracy_assigned": float((c[m][assigned] == truth).mean())
                    if assigned.any()
                    else np.nan,
                }
            )
    vdf = pd.DataFrame(val)
    bal = vdf.groupby(["rule", "truth"]).accuracy_assigned.mean().groupby("rule").mean()
    cov = vdf.groupby("rule").coverage.mean()
    vdf["rule_balanced_accuracy"] = vdf.rule.map(bal)
    vdf["rule_mean_coverage"] = vdf.rule.map(cov)
    vdf.to_csv(out / "c004_chamber_validation.csv", index=False)
    best = bal.idxmax()
    usable = bool(bal[best] >= 0.9 and cov[best] >= 0.5)
    ch = chamber[best]
    print(
        "chamber",
        bal.round(3).to_dict(),
        cov.round(3).to_dict(),
        "chosen",
        best,
        usable,
        flush=True,
    )

    # ---- thinning (identical to C003) ----
    rng = np.random.default_rng(SEED)
    p = np.minimum(1.0, TARGET / tot)
    th = counts.copy()
    th.data = rng.binomial(th.data.astype(np.int64), np.repeat(p, np.diff(th.indptr))).astype(
        np.float32
    )
    th.eliminate_zeros()
    b = ad.AnnData(X=th, var=pd.DataFrame(index=vn), obs=pd.DataFrame(index=a.obs_names))
    sc.pp.normalize_total(b, target_sum=1e4)
    sc.pp.log1p(b)
    full = ad.AnnData(X=counts.copy(), var=pd.DataFrame(index=vn))
    sc.pp.normalize_total(full, target_sum=1e4)
    sc.pp.log1p(full)
    means = np.asarray(full.X.mean(0)).ravel()
    axes = {}
    for f in sorted(Path(marker_dir).glob("*.yaml")):
        spec = yaml.safe_load(f.read_text())
        axes[spec["state"]] = [
            m["symbol"]
            for m in spec["genes"]
            if m["symbol"] in vidx and m["symbol"] not in CIRCULAR
        ]
    all_idx = np.array(sorted({vidx[g] for gs in axes.values() for g in gs}))

    # ---- 4. agent labels (all cells) ----
    real, null = score_all(b, axes, vn, means, all_idx)
    states = list(axes)
    z = np.full((len(lib), len(states)), -np.inf)
    for j, st in enumerate(states):
        for h in set(hrt):
            m = hrt == h
            thr = np.percentile(null[st][m].ravel(), 99)
            sd = null[st][m].ravel().std()
            exc = (real[st][m] - thr) / sd
            z[m, j] = np.where(real[st][m] > thr, exc, -np.inf)
    high = np.isfinite(z)
    primary = np.where(high.any(1), np.array(states)[np.argmax(z, 1)], "unassigned")
    labels = pd.DataFrame(
        {
            "cell": a.obs_names,
            "library": lib,
            "heart": hrt,
            "agent_state": primary,
            "agent_states_all": [";".join(np.array(states)[r]) for r in high],
            "agent_chamber": ch,
        }
    )
    for st in states:
        labels[f"score_{st}"] = np.round(real[st], 5)
    labels.to_csv(out / "c004_agent_labels.csv", index=False)
    lab_summary = labels.groupby("agent_state").size().to_dict()
    print("labels", lab_summary, flush=True)

    # ---- 2. positive control, ventricle-called cells ----
    vm = ch == "ventricle"
    bv = b[vm].copy()
    rv, nv = score_all(bv, axes, vn, means, all_idx)
    hv = hrt[vm]
    hrows, crow = [], []
    for st in states:
        for h in sorted(time_of):
            m = hv == h
            n = int(m.sum())
            thr = np.percentile(nv[st][m].ravel(), 99) if n else np.nan
            hrows.append(
                {
                    "state": st,
                    "heart": h,
                    "time": time_of[h],
                    "n_ventricle_cells": n,
                    "included": n >= MIN_CELLS,
                    "median_thinned_score": float(np.median(rv[st][m])) if n else np.nan,
                    "tail_fraction": float((rv[st][m] > thr).mean()) if n else np.nan,
                }
            )
        hd = pd.DataFrame([r for r in hrows if r["state"] == st and r["included"]])
        ctl = hd[hd.time == "Ctrl"]
        for t in PRED[st]:
            inj = hd[hd.time == t]
            r = {
                "state": st,
                "contrast": f"{t}_vs_Ctrl",
                "n_injured_hearts": len(inj),
                "n_ctrl_hearts": len(ctl),
            }
            for stat in ("tail_fraction", "median_thinned_score"):
                r[f"{stat}_diff_of_medians"] = float(inj[stat].median() - ctl[stat].median())
                r[f"{stat}_auc"] = auc(inj[stat].tolist(), ctl[stat].tolist())
                r[f"{stat}_pass"] = bool(
                    r[f"{stat}_auc"] >= 0.8 and r[f"{stat}_diff_of_medians"] > 0
                )
            crow.append(r)
    pd.DataFrame(hrows).to_csv(out / "c004_ventricle_heart_scores.csv", index=False)
    pc = pd.DataFrame(crow)
    pc.to_csv(out / "c004_ventricle_positive_control.csv", index=False)
    passing = sorted(set(pc[pc.tail_fraction_pass].state))

    # ---- 3. G4 gate ----
    from collections import Counter

    cnt = Counter(g for gs in axes.values() for g in gs)
    shared = sorted(g for g, c in cnt.items() if c > 1)
    disjoint = {st: [g for g in gs if cnt[g] == 1] for st, gs in axes.items()}
    g4 = {"axes_passing_tail": passing, "shared_genes_removed": shared, "disjoint_sets": disjoint}
    if passing:
        rd, nd = score_all(
            bv, {s: g for s, g in disjoint.items() if len(g) >= 2}, vn, means, all_idx
        )
        rows = []
        for s1, s2 in combinations(rd, 2):
            for h in sorted(set(hv)):
                m = hv == h
                if m.sum() < MIN_CELLS:
                    continue
                rr = spearmanr(rd[s1][m], rd[s2][m])[0]
                nr = [spearmanr(nd[s1][m, j], nd[s2][m, j])[0] for j in range(N_NULL)]
                rows.append(
                    {
                        "pair": f"{s1}|{s2}",
                        "heart": h,
                        "rho": rr,
                        "null_p05": np.percentile(nr, 5),
                        "null_p95": np.percentile(nr, 95),
                    }
                )
        gd = pd.DataFrame(rows)
        gd.to_csv(out / "c004_g4_cooccurrence.csv", index=False)
        summ = {}
        for pair, x in gd.groupby("pair"):
            up, dn = (x.rho > x.null_p95).mean(), (x.rho < x.null_p05).mean()
            summ[pair] = {
                "frac_hearts_above_null95": up,
                "frac_hearts_below_null05": dn,
                "verdict": "co-occur" if up >= 0.8 else "come apart" if dn >= 0.8 else "mixed",
            }
        g4["pairs"] = summ
    else:
        g4["verdict"] = "cannot be answered with these sets: no axis passes the positive control"
    (out / "c004_g4.json").write_text(json.dumps(g4, indent=2, default=str))

    meta = {
        "commit": subprocess.check_output(["git", "rev-parse", "HEAD"], text=True).strip(),
        "dirty_tree": True,
        "command": "uv run python " + " ".join(sys.argv),
        "input": src,
        "chamber_rule_chosen": best,
        "chamber_rule_usable": usable,
        "chamber_balanced_accuracy": bal.round(4).to_dict(),
        "chamber_mean_coverage": cov.round(4).to_dict(),
        "agent_state_counts": {k: int(v) for k, v in lab_summary.items()},
        "agent_chamber_counts": {
            k: int(v) for k, v in pd.Series(ch).value_counts().to_dict().items()
        },
        "label_file_format": {
            "file": "c004_agent_labels.csv",
            "key": "cell = obs_names of the sanitised h5ad (one row per cell, 32,044 rows)",
            "agent_state": "one of " + ", ".join(states) + ", unassigned",
            "agent_states_all": "semicolon-separated high axes (empty if none)",
            "agent_chamber": "atrium | ventricle | unassigned",
            "score_<axis>": "thinned score_genes per axis",
            "rule": "research/cm-state-axes/cycles/C004/reading-rules.md section 4",
        },
        "params": {"target_umis": TARGET, "seed": SEED, "n_null": N_NULL, "min_cells": MIN_CELLS},
        "time": datetime.now().astimezone().isoformat(timespec="seconds"),
    }
    (out / "c004_probe.json").write_text(json.dumps(meta, indent=2))
    print(json.dumps({k: meta[k] for k in list(meta)[5:10]}, indent=1))
    print(pc.round(3).to_string(index=False))
    print(json.dumps(g4, indent=1, default=str)[:1500])


if __name__ == "__main__":
    main(*sys.argv[1:4])
