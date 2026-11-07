"""C002: binomial depth thinning to equal depth, then rescore (R3.6, stronger form).

uv run python scripts/c002_thinning.py data/clean/hu2022_cm_allgenes_sanitized.h5ad markers results
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
import scipy.sparse as sp
import yaml
from scipy.stats import spearmanr

from zf.io import ALLOWED_OBS
from zf.score import rank_matrix, ucell_from_ranks

TARGET, SEED, MAX_RANK, CIRCULAR = 500, 0, 1500, {"ttn.2"}


def scores(counts, var_names, axes, seed):
    a = ad.AnnData(X=counts.copy(), var=pd.DataFrame(index=var_names))
    sc.pp.normalize_total(a, target_sum=1e4)
    sc.pp.log1p(a)
    R = rank_matrix(counts)
    vidx = {g: i for i, g in enumerate(var_names)}
    out = {}
    for st, genes in axes.items():
        out[(st, "ucell")] = ucell_from_ranks(R, [vidx[g] for g in genes], MAX_RANK)
        sc.tl.score_genes(a, genes, random_state=seed, score_name="_s", use_raw=False)
        out[(st, "score_genes")] = a.obs.pop("_s").to_numpy()
    return out


def main(src, marker_dir, out_dir):
    a = ad.read_h5ad(src)
    assert set(a.obs.columns) == set(ALLOWED_OBS)
    counts = a.layers["counts"].tocsr().astype(np.float32)
    tot = np.asarray(counts.sum(1)).ravel()
    lib = a.obs["orig.ident"].astype(str).to_numpy()
    vn = list(a.var_names)
    axes = {}
    for f in sorted(Path(marker_dir).glob("*.yaml")):
        spec = yaml.safe_load(f.read_text())
        axes[spec["state"]] = [
            m["symbol"] for m in spec["genes"] if m["symbol"] in vn and m["symbol"] not in CIRCULAR
        ]
    rng = np.random.default_rng(SEED)
    p = np.minimum(1.0, TARGET / tot)
    th = counts.copy()
    th.data = rng.binomial(th.data.astype(np.int64), np.repeat(p, np.diff(th.indptr))).astype(
        np.float32
    )
    th.eliminate_zeros()
    th = sp.csr_matrix(th)
    tot_th = np.asarray(th.sum(1)).ravel()
    s0, s1 = scores(counts, vn, axes, SEED), scores(th, vn, axes, SEED)
    logtot = np.log10(tot)
    rows = []
    for lb in sorted(set(lib)):
        m = lib == lb
        for (st, meth), v0 in s0.items():
            v1 = s1[(st, meth)]
            rows.append(
                {
                    "state": st,
                    "method": meth,
                    "lib": lb,
                    "n_cells": int(m.sum()),
                    "n_below_target": int((tot[m] < TARGET).sum()),
                    "rho_orig_logtot": float(spearmanr(v0[m], logtot[m])[0]),
                    "rho_thinned_logtot": float(spearmanr(v1[m], logtot[m])[0]),
                    "rho_orig_vs_thinned": float(spearmanr(v0[m], v1[m])[0]),
                }
            )
    d = pd.DataFrame(rows)
    d.to_csv(Path(out_dir) / "c002_thinning_by_library.csv", index=False)
    s = (
        d.assign(a0=d.rho_orig_logtot.abs(), a1=d.rho_thinned_logtot.abs())
        .groupby(["state", "method"])
        .agg(
            median_abs_rho_orig=("a0", "median"),
            median_abs_rho_thinned=("a1", "median"),
            median_rho_orig_vs_thinned=("rho_orig_vs_thinned", "median"),
        )
    )
    s["fraction_of_tie_left"] = s.median_abs_rho_thinned / s.median_abs_rho_orig
    s.to_csv(Path(out_dir) / "c002_thinning_summary.csv")
    meta = {
        "commit": subprocess.check_output(["git", "rev-parse", "HEAD"], text=True).strip(),
        "dirty_tree": True,
        "command": "uv run python " + " ".join(sys.argv),
        "input": src,
        "target_umis": TARGET,
        "seed": SEED,
        "cells_below_target": int((tot < TARGET).sum()),
        "median_total_after_thinning": float(np.median(tot_th)),
        "time": datetime.now().astimezone().isoformat(timespec="seconds"),
    }
    (Path(out_dir) / "c002_thinning.json").write_text(json.dumps(meta, indent=2))
    print(json.dumps(meta, indent=1))
    print(s.round(3).to_string())


if __name__ == "__main__":
    main(*sys.argv[1:4])
