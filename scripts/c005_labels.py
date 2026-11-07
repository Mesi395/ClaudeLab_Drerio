"""C005 fix F10: agent state labels for ventricle-called cells only; other cells not_scored.

Rules: research/cm-state-axes/cycles/C005/reading-rules.md section 1 (fixed before computing).
uv run python scripts/c005_labels.py data/clean/hu2022_cm_allgenes_sanitized.h5ad markers results
"""

import json
import subprocess
import sys
from datetime import datetime
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import yaml

sys.path.insert(0, str(Path(__file__).parent))
from c004_probe import CIRCULAR, SEED, TARGET, heart, score_all  # noqa: E402

from zf.io import ALLOWED_OBS  # noqa: E402


def main(src, marker_dir, out_dir):
    out = Path(out_dir)
    a = ad.read_h5ad(src)
    assert set(a.obs.columns) == set(ALLOWED_OBS)
    old = pd.read_csv(out / "c004_agent_labels.csv")
    assert (old.cell.to_numpy() == a.obs_names.to_numpy()).all(), "cell order mismatch"
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
    axes = {}
    for f in sorted(Path(marker_dir).glob("*.yaml")):
        spec = yaml.safe_load(f.read_text())
        axes[spec["state"]] = [
            m["symbol"]
            for m in spec["genes"]
            if m["symbol"] in vidx and m["symbol"] not in CIRCULAR
        ]
    all_idx = np.array(sorted({vidx[g] for gs in axes.values() for g in gs}))

    vm = old.agent_chamber.to_numpy() == "ventricle"
    bv = b[vm].copy()
    real, null = score_all(bv, axes, vn, means, all_idx)
    states = list(axes)
    hv = hrt[vm]
    z = np.full((vm.sum(), len(states)), -np.inf)
    for j, st in enumerate(states):
        for h in set(hv):
            m = hv == h
            nl = null[st][m].ravel()
            thr, sd = np.percentile(nl, 99), nl.std()
            z[m, j] = np.where(real[st][m] > thr, (real[st][m] - thr) / sd, -np.inf)
    high = np.isfinite(z)
    prim = np.where(high.any(1), np.array(states)[np.argmax(z, 1)], "unassigned")

    new = old[["cell", "library", "heart", "agent_chamber"]].copy()
    new["agent_state"] = "not_scored"
    new["agent_states_all"] = ""
    new.loc[vm, "agent_state"] = prim
    new.loc[vm, "agent_states_all"] = [";".join(np.array(states)[r]) for r in high]
    for st in states:
        new[f"score_{st}"] = np.nan
        new.loc[vm, f"score_{st}"] = np.round(real[st], 5)
    new = new[
        ["cell", "library", "heart", "agent_state", "agent_states_all", "agent_chamber"]
        + [f"score_{s}" for s in states]
    ]
    new.to_csv(out / "c005_agent_labels.csv", index=False)

    n_scored = int(vm.sum())
    cmp_ = (
        pd.DataFrame(
            {
                "c004_count": old.agent_state.value_counts(),
                "c005_count": new.agent_state.value_counts(),
            }
        )
        .fillna(0)
        .astype(int)
    )
    cmp_["c005_chance_expectation"] = np.where(cmp_.index.isin(states), round(0.01 * n_scored), 0)
    cmp_["c005_ratio_to_chance"] = np.where(
        cmp_.index.isin(states), cmp_.c005_count / (0.01 * n_scored), np.nan
    )
    cmp_["c005_at_chance_0.5_2x"] = cmp_.c005_ratio_to_chance.between(0.5, 2)
    cmp_.to_csv(out / "c005_label_counts_vs_c004.csv")
    by_heart = pd.crosstab(new.heart[vm], new.agent_state[vm])
    by_heart.to_csv(out / "c005_labels_by_heart.csv")
    old_v = old.agent_state[vm].to_numpy()
    meta = {
        "commit": subprocess.check_output(["git", "rev-parse", "HEAD"], text=True).strip(),
        "dirty_tree": True,
        "command": "uv run python " + " ".join(sys.argv),
        "input": src,
        "supersedes": "c004_agent_labels.csv (kept unchanged); correction for F10 (chamber)",
        "n_cells": int(len(new)),
        "n_scored_ventricle": n_scored,
        "agreement_with_c004_on_ventricle_cells": float((old_v == prim).mean()),
        "label_file_format": {
            "file": "c005_agent_labels.csv",
            "cell": "obs_names of the sanitised h5ad; one row per cell (32,044)",
            "agent_state": "border_zone_like | dedifferentiated | injury_responsive | "
            "proliferating | unassigned (ventricle-called cells) | not_scored "
            "(atrium-called or chamber-unassigned cells)",
            "agent_states_all": "semicolon-separated high axes; empty if none or not scored",
            "agent_chamber": "atrium | ventricle | unassigned; sourced rule B, see C004",
            "score_<axis>": "thinned score_genes within ventricle-called cells; NaN if not scored",
            "rule": "research/cm-state-axes/cycles/C005/reading-rules.md section 1",
        },
        "time": datetime.now().astimezone().isoformat(timespec="seconds"),
    }
    (out / "c005_agent_labels.json").write_text(json.dumps(meta, indent=2))
    print(cmp_.round(3).to_string())
    print(
        json.dumps(
            {k: meta[k] for k in ("n_scored_ventricle", "agreement_with_c004_on_ventricle_cells")}
        )
    )
    print(by_heart.to_string())


if __name__ == "__main__":
    main(*sys.argv[1:4])
