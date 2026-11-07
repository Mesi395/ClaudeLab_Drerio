"""C007: ventricle-only positive control with control-calibrated tails and a null-AUC test, for the
C004 sourced sets and the Li 2025 SD2 sets, side by side; leave-Hr9-out.

Rules: research/cm-state-axes/cycles/C007/reading-rules.md (fixed 16:06 CET).
uv run python scripts/c007_probe.py data/clean/hu2022_cm_allgenes_sanitized.h5ad markers \
    proposals/loop/li2025_sd2_sets.tsv results [tag]
"""

import json
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

sys.path.insert(0, str(Path(__file__).parent))
from c004_probe import HEARTS, heart  # noqa: E402

TARGET, SEED, N_NULL, CIRCULAR = 500, 0, 100, {"ttn.2"}
PRED = {
    "injury_responsive": ["3dpi"],
    "dedifferentiated": ["3dpi", "7dpi"],
    "proliferating": ["7dpi"],
    "border_zone_like": ["7dpi", "3dpi"],
}
CTRL = HEARTS["Ctrl"]


def auc(inj, ctl):
    return float(np.mean([1.0 if i > c else 0.5 if i == c else 0.0 for i, c in product(inj, ctl)]))


def sg(b, genes, seed):
    sc.tl.score_genes(b, list(genes), random_state=seed, score_name="_s", use_raw=False)
    return b.obs.pop("_s").to_numpy()


def tails(score, hv):
    """Control-calibrated tail per heart; control hearts use leave-one-out thresholds."""
    out = {}
    ctrl_m = np.isin(hv, CTRL)
    t_all = np.percentile(score[ctrl_m], 99)
    for h in np.unique(hv):
        t = np.percentile(score[ctrl_m & (hv != h)], 99) if h in CTRL else t_all
        out[h] = float((score[hv == h] > t).mean())
    return out


def contrast(tl, t, drop=()):
    hs = [h for h in tl if h not in drop]
    inj = [tl[h] for h in hs if h in HEARTS[t]]
    ctl = [tl[h] for h in hs if h in CTRL]
    return auc(inj, ctl), float(np.median(inj) - np.median(ctl))


def main(src, marker_dir, li_tsv, out_dir, tag="c007"):
    out = Path(out_dir)
    a = ad.read_h5ad(src)
    assert set(a.obs.columns) == set(ALLOWED_OBS)
    lab = pd.read_csv(out / "c004_agent_labels.csv")
    assert (lab.cell.to_numpy() == a.obs_names.to_numpy()).all()
    vm = lab.agent_chamber.to_numpy() == "ventricle"
    counts = a.layers["counts"].tocsr().astype(np.float32)
    tot = np.asarray(counts.sum(1)).ravel()
    lib = a.obs["orig.ident"].astype(str).to_numpy()
    hrt = np.array([heart(x) for x in lib])
    vn = list(a.var_names)
    vidx = {g: i for i, g in enumerate(vn)}
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
    bv = b[vm].copy()
    hv = hrt[vm]

    sets = {}
    for f in sorted(Path(marker_dir).glob("*.yaml")):
        spec = yaml.safe_load(f.read_text())
        sets[("C004", spec["state"])] = [
            m["symbol"]
            for m in spec["genes"]
            if m["symbol"] in vidx and m["symbol"] not in CIRCULAR
        ]
    li = pd.read_csv(li_tsv, sep="\t")
    for ax, x in li.groupby("axis"):
        sets[("Li_SD2", ax)] = [g for g in x.sort_values("rank").gene if g in vidx]
    excl = np.array(sorted({vidx[g] for gs in sets.values() for g in gs}))

    rows, heart_rows = [], []
    for k, ((family, ax), genes) in enumerate(sets.items()):
        gi = np.array([vidx[g] for g in genes])
        real = sg(bv, genes, SEED)
        tl = tails(real, hv)
        nulls = matched_random_sets(means, gi, N_NULL, 25, SEED + k, exclude=excl)
        ntl = [tails(sg(bv, [vn[i] for i in s], SEED + j), hv) for j, s in enumerate(nulls)]
        for h, v in tl.items():
            heart_rows.append(
                {
                    "family": family,
                    "axis": ax,
                    "heart": h,
                    "time": next(t for t, hs in HEARTS.items() if h in hs),
                    "n_ventricle_cells": int((hv == h).sum()),
                    "tail_fraction": v,
                    "null_median_tail": float(np.median([n[h] for n in ntl])),
                }
            )
        for t in PRED[ax]:
            for drop in ((), ("Hr9",)):
                a_, d_ = contrast(tl, t, drop)
                na = [contrast(n, t, drop)[0] for n in ntl]
                r = {
                    "family": family,
                    "axis": ax,
                    "n_genes": len(genes),
                    "contrast": f"{t}_vs_Ctrl",
                    "hr9_excluded": bool(drop),
                    "auc": a_,
                    "diff_of_medians": d_,
                    "null_auc_p95": float(np.percentile(na, 95)),
                    "null_auc_median": float(np.median(na)),
                }
                r["pass"] = bool(a_ >= 0.8 and d_ > 0 and a_ > r["null_auc_p95"])
                rows.append(r)
        print(family, ax, "done", flush=True)
    res = pd.DataFrame(rows)
    res.to_csv(out / f"{tag}_positive_control.csv", index=False)
    pd.DataFrame(heart_rows).to_csv(out / f"{tag}_heart_tails.csv", index=False)
    meta = {
        "commit": subprocess.check_output(["git", "rev-parse", "HEAD"], text=True).strip(),
        "dirty_tree": True,
        "command": "uv run python " + " ".join(sys.argv),
        "input": src,
        "sets": {f"{f}:{a_}": g for (f, a_), g in sets.items()},
        "params": {
            "target_umis": TARGET,
            "seed": SEED,
            "n_null": N_NULL,
            "threshold": "99th pct of control-heart ventricle cells, LOO for controls",
        },
        "time": datetime.now().astimezone().isoformat(timespec="seconds"),
    }
    (out / f"{tag}_probe.json").write_text(json.dumps(meta, indent=2))
    pd.set_option("display.width", 250)
    print(res.round(3).to_string(index=False))


if __name__ == "__main__":
    main(*sys.argv[1:6])
