"""Per-library depth, genes detected, mito fraction and chemistry; marker presence and detection.

Reads ONLY the sanitised copy (see zf.io.sanitize_h5ad). Study cm-state-axes, cycle C001.

uv run python scripts/data_facts.py data/clean/hu2022_cm_allgenes_sanitized.h5ad markers results
"""

import json
import subprocess
import sys
from datetime import datetime
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import yaml

from zf.io import ALLOWED_OBS

# 10x chemistry per GEO (docs/data.md): Hr24 and Hr27 are v3, the others v2.
V3_LIBS = {"Hr24", "Hr27"}


def main(src: str, marker_dir: str, out_dir: str) -> None:
    a = ad.read_h5ad(src)
    assert set(a.obs.columns) == set(ALLOWED_OBS), "input is not the sanitised file"
    out = Path(out_dir)
    out.mkdir(exist_ok=True)
    X = a.layers["counts"].tocsr()
    total = np.asarray(X.sum(1)).ravel()
    ngenes = np.asarray((X > 0).sum(1)).ravel()
    mt = a.var_names.str.startswith("mt-")
    mtfrac = np.asarray(X[:, mt].sum(1)).ravel() / np.maximum(total, 1)
    df = pd.DataFrame(
        {
            "lib": a.obs["orig.ident"].astype(str).to_numpy(),
            "time": a.obs["time"].astype(str).to_numpy(),
            "AV": a.obs["AV"].astype(str).to_numpy(),
            "total_counts": total,
            "n_genes": ngenes,
            "mt_frac": mtfrac,
            "nCount_RNA": a.obs["nCount_RNA"].to_numpy(),
            "nFeature_RNA": a.obs["nFeature_RNA"].to_numpy(),
            "percent_mito": a.obs["percent.mito"].to_numpy(),
        }
    )
    g = df.groupby("lib")
    lib = pd.DataFrame(
        {
            "time": g["time"].first(),
            "AV": g["AV"].first(),
            "chemistry": [("v3" if lb in V3_LIBS else "v2") for lb in g.groups],
            "n_cells": g.size(),
            "median_total_counts": g["total_counts"].median(),
            "median_n_genes": g["n_genes"].median(),
            "median_mt_frac": g["mt_frac"].median(),
            "median_percent_mito_supplied": g["percent_mito"].median(),
            "max_abs_count_diff_vs_nCount_RNA": g.apply(
                lambda d: float(np.max(np.abs(d.total_counts - d.nCount_RNA)))
            ),
        }
    ).sort_values(["time", "median_total_counts"])
    lib.to_csv(out / "c001_library_qc.csv")

    # marker presence and per-library detection fraction
    rows = []
    vn = set(a.var_names)
    libs = df["lib"].to_numpy()
    for f in sorted(Path(marker_dir).glob("*.yaml")):
        spec = yaml.safe_load(f.read_text())
        for m in spec["genes"]:
            gene = m["symbol"]
            row = {"state": spec["state"], "gene": gene, "present": gene in vn}
            if gene in vn:
                det = np.asarray((X[:, a.var_names.get_loc(gene)] > 0).todense()).ravel()
                for lb in sorted(set(libs)):
                    row[lb] = round(float(det[libs == lb].mean()), 4)
                row["all"] = round(float(det.mean()), 4)
            rows.append(row)
    pd.DataFrame(rows).to_csv(out / "c001_marker_detection.csv", index=False)

    meta = {
        "commit": subprocess.check_output(["git", "rev-parse", "HEAD"], text=True).strip(),
        "dirty_tree": bool(subprocess.check_output(["git", "status", "--porcelain"], text=True)),
        "command": "uv run python " + " ".join(sys.argv),
        "input": src,
        "time": datetime.now().astimezone().isoformat(timespec="seconds"),
        "n_cells": int(a.n_obs),
        "n_genes": int(a.n_vars),
        "n_libraries": int(df["lib"].nunique()),
        "n_mt_genes": int(mt.sum()),
        "missing_values_in_covariates": int(df.isna().sum().sum()),
        "counts_integer": bool(np.all(X.data == np.round(X.data))),
        "outputs": ["c001_library_qc.csv", "c001_marker_detection.csv"],
    }
    (out / "c001_data_facts.json").write_text(json.dumps(meta, indent=2))
    print(json.dumps(meta, indent=1))
    print(lib.round(3).to_string())


if __name__ == "__main__":
    main(*sys.argv[1:4])
