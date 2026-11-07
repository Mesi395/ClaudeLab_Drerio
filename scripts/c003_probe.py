"""C003 Probe: positive control by time per heart (thinned scores, tail vs null), and the Hr11
ambient-RNA test. Rules: research/cm-state-axes/cycles/C003/reading-rules.md (fixed first).

uv run python scripts/c003_probe.py data/clean/hu2022_cm_allgenes_sanitized.h5ad \
    proposals/loop results
"""

import json
import re
import subprocess
import sys
from datetime import datetime
from itertools import product
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import yaml

from zf.io import ALLOWED_OBS
from zf.score import matched_random_sets

TARGET, SEED, N_NULL, CIRCULAR = 500, 0, 100, {"ttn.2"}
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
NON_CM = [
    "hbba1.1",
    "hbaa1.1",
    "hbba1",
    "col1a1a",
    "col1a2",
    "postnb",
    "fn1b",
    "kdrl",
    "cdh5",
    "fli1a",
    "mpeg1.1",
    "lyz",
    "tcf21",
]
BINS = [0, 600, 900, 1350, np.inf]


def heart(lib):
    return re.sub(r"^(Hr2|H8|Hr6|Hr7)[abv]$", r"\1", lib)


def sg(adata, genes, seed):
    sc.tl.score_genes(adata, list(genes), random_state=seed, score_name="_s", use_raw=False)
    return adata.obs.pop("_s").to_numpy()


def auc(inj, ctl):
    pairs = list(product(inj, ctl))
    return float(np.mean([1.0 if i > c else 0.5 if i == c else 0.0 for i, c in pairs]))


def main(src, marker_dir, out_dir):
    out = Path(out_dir)
    a = ad.read_h5ad(src)
    assert set(a.obs.columns) == set(ALLOWED_OBS)
    counts = a.layers["counts"].tocsr().astype(np.float32)
    tot = np.asarray(counts.sum(1)).ravel()
    lib = a.obs["orig.ident"].astype(str).to_numpy()
    hrt = np.array([heart(x) for x in lib])
    time_of = {h: t for t, hs in HEARTS.items() for h in hs}
    assert set(hrt) == set(time_of), "heart map mismatch"
    vn = list(a.var_names)
    vidx = {g: i for i, g in enumerate(vn)}

    # ---------- Q1: thinned score_genes, tail vs matched null, per heart ----------
    rng = np.random.default_rng(SEED)
    p = np.minimum(1.0, TARGET / tot)
    th = counts.copy()
    th.data = rng.binomial(th.data.astype(np.int64), np.repeat(p, np.diff(th.indptr))).astype(
        np.float32
    )
    th.eliminate_zeros()
    b = ad.AnnData(X=th, var=pd.DataFrame(index=vn))
    sc.pp.normalize_total(b, target_sum=1e4)
    sc.pp.log1p(b)
    # null bins from unthinned means, as in C002
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
    hrows, crow = [], []
    for k, (st, genes) in enumerate(axes.items()):
        gi = np.array([vidx[g] for g in genes])
        real = sg(b, genes, SEED)
        nulls = matched_random_sets(means, gi, N_NULL, 25, SEED + k, exclude=all_idx)
        nul = np.stack([sg(b, [vn[i] for i in s], SEED + j) for j, s in enumerate(nulls)], 1)
        for h in sorted(time_of):
            m = hrt == h
            thr = np.percentile(nul[m].ravel(), 99)
            hrows.append(
                {
                    "state": st,
                    "heart": h,
                    "time": time_of[h],
                    "n_cells": int(m.sum()),
                    "median_thinned_score": float(np.median(real[m])),
                    "null99": float(thr),
                    "tail_fraction": float((real[m] > thr).mean()),
                }
            )
        hd = pd.DataFrame([r for r in hrows if r["state"] == st])
        for excl in (False, True):
            d = hd[hd.heart != "Hr11"] if excl else hd
            ctl = d[d.time == "Ctrl"]
            for t in PRED[st]:
                inj = d[d.time == t]
                r = {
                    "state": st,
                    "contrast": f"{t}_vs_Ctrl",
                    "hr11_excluded": excl,
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
        print(st, "done", flush=True)
    pd.DataFrame(hrows).to_csv(out / "c003_heart_scores.csv", index=False)
    pd.DataFrame(crow).to_csv(out / "c003_positive_control.csv", index=False)

    # ---------- Q2: ambient test ----------
    t_of_lib = a.obs["time"].astype(str).to_numpy()
    group = np.where(
        np.isin(lib, ["Hr11", "Hr12", "Hr24", "Hr27"]),
        lib,
        np.where(t_of_lib == "Ctrl", "Ctrl", ""),
    )
    ctrl_m = group == "Ctrl"
    excl_pat = re.compile(r"^(mt-|NC-|rp[ls]|RP[LS]|hb)")
    share_ctrl = np.asarray(counts[ctrl_m].sum(0)).ravel() / tot[ctrl_m].sum()
    order = np.argsort(-share_ctrl)
    abundant = [
        vn[i] for i in order if not excl_pat.match(vn[i]) and vn[i] not in ("nppa", "ttn.2", "myh6")
    ][:10]
    panels = {
        "target": ["nppa", "ttn.2"],
        "atrial": ["myh6"],
        "abundant_cm": abundant,
        "non_cm": [g for g in NON_CM if g in vidx],
    }
    mito = np.array([bool(re.match(r"^(mt-|NC-002333)", g)) for g in vn])
    mito_share = np.asarray(counts[:, mito].sum(1)).ravel() / tot
    bin_of = np.digitize(tot, BINS[1:-1])
    rows = []
    genes = [g for gs in panels.values() for g in gs]
    gidx = [vidx[g] for g in genes]
    S = counts[:, gidx].toarray() / tot[:, None]
    for grp in ["Hr11", "Hr12", "Hr24", "Hr27"]:
        for bb in range(4):
            gm, cm = (group == grp) & (bin_of == bb), ctrl_m & (bin_of == bb)
            if gm.sum() < 20 or cm.sum() < 20:
                continue
            med_depth = float(np.median(tot[gm]))
            for pan, gs in panels.items():
                for g in gs:
                    j = genes.index(g)
                    sg_, sc_ = S[gm, j].mean(), S[cm, j].mean()
                    rows.append(
                        {
                            "group": grp,
                            "depth_bin": bb,
                            "n_group": int(gm.sum()),
                            "n_ctrl": int(cm.sum()),
                            "panel": pan,
                            "gene": g,
                            "mean_share_group": sg_,
                            "mean_share_ctrl": sc_,
                            "ratio": sg_ / sc_ if sc_ > 0 else np.inf,
                            "excess_molecules_per_cell": (sg_ - sc_) * med_depth,
                            "bin_median_depth": med_depth,
                        }
                    )
            rows.append(
                {
                    "group": grp,
                    "depth_bin": bb,
                    "n_group": int(gm.sum()),
                    "n_ctrl": int(cm.sum()),
                    "panel": "technical",
                    "gene": "mito_incl_NC",
                    "mean_share_group": mito_share[gm].mean(),
                    "mean_share_ctrl": mito_share[cm].mean(),
                    "ratio": mito_share[gm].mean() / mito_share[cm].mean(),
                    "excess_molecules_per_cell": np.nan,
                    "bin_median_depth": med_depth,
                }
            )
    amb = pd.DataFrame(rows)
    amb.to_csv(out / "c003_ambient_by_bin.csv", index=False)

    # verdict components (rules fixed in advance)
    sh = amb[amb.depth_bin.isin([0, 1])]
    v = {}
    for grp in ["Hr11", "Hr12", "Hr24", "Hr27"]:
        x = sh[sh.group == grp]
        v[grp] = {
            "nppa_ratio_bins01_min": float(x[x.gene == "nppa"].ratio.min()) if len(x) else None,
            "nppa_ge2x_bins01": bool(len(x) and (x[x.gene == "nppa"].ratio >= 2).all()),
            "non_cm_frac_ge2x_bins01": float(
                x[x.panel == "non_cm"].groupby("gene").ratio.min().ge(2).mean()
            )
            if len(x)
            else None,
            "abundant_cm_frac_ge2x_bins01": float(
                x[x.panel == "abundant_cm"].groupby("gene").ratio.min().ge(2).mean()
            )
            if len(x)
            else None,
        }
    hn = amb[(amb.group == "Hr11") & (amb.gene == "nppa")].sort_values("depth_bin")
    em = hn.excess_molecules_per_cell.to_numpy()
    v["Hr11_nppa_ratio_by_bin"] = hn.ratio.round(3).tolist()
    v["Hr11_nppa_excess_molecules_by_bin"] = np.round(em, 3).tolist()
    v["Hr11_nppa_ratio_falls_with_depth"] = bool(np.all(np.diff(hn.ratio.to_numpy()) < 0))
    v["Hr11_nppa_excess_molecules_within_2x"] = bool(em.min() > 0 and em.max() / em.min() <= 2)
    others = sum(v[g]["nppa_ge2x_bins01"] for g in ["Hr12", "Hr24", "Hr27"])
    v["n_other_3dpi_with_nppa_ge2x"] = int(others)
    hr11_specific = others <= 1
    noncm = (v["Hr11"]["non_cm_frac_ge2x_bins01"] or 0) >= 0.5
    if (
        hr11_specific
        and noncm
        and v["Hr11_nppa_ratio_falls_with_depth"]
        and v["Hr11_nppa_excess_molecules_within_2x"]
    ):
        verdict = "ambient favoured"
    elif others >= 2 and not noncm:
        verdict = "3 dpi biology favoured"
    else:
        verdict = "inconclusive"
    v["verdict"] = verdict
    v["parts"] = {"hr11_specific": hr11_specific, "non_cm_raised": noncm}
    v["panels"] = panels
    (out / "c003_ambient_verdict.json").write_text(json.dumps(v, indent=2))

    meta = {
        "commit": subprocess.check_output(["git", "rev-parse", "HEAD"], text=True).strip(),
        "dirty_tree": True,
        "command": "uv run python " + " ".join(sys.argv),
        "input": src,
        "params": {"target_umis": TARGET, "seed": SEED, "n_null": N_NULL, "depth_bins": BINS[:-1]},
        "time": datetime.now().astimezone().isoformat(timespec="seconds"),
        "outputs": [
            "c003_heart_scores.csv",
            "c003_positive_control.csv",
            "c003_ambient_by_bin.csv",
            "c003_ambient_verdict.json",
        ],
    }
    (out / "c003_probe.json").write_text(json.dumps(meta, indent=2))
    print(json.dumps(v, indent=1, default=str))


if __name__ == "__main__":
    main(*sys.argv[1:4])
