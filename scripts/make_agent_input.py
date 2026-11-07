"""Make the only data file agents may read, from the processed cardiomyocyte AnnData.

Keeps raw counts and an allowlist of technical and design columns; drops everything else
(other cell columns, embeddings, gene annotations, unstructured metadata). Run once by the
data owner, on their machine:

    uv run python scripts/make_agent_input.py PATH/TO/processed.h5ad

Writes:
    data/processed/cardiomyocytes.h5ad   counts + allowlisted columns (git-ignored)
    config/cells.csv.gz                  cell barcode and library (committed)
"""

import hashlib
import json
import sys
from pathlib import Path

import anndata as ad
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
KEEP_OBS = ["orig.ident", "nCount_RNA", "nFeature_RNA", "percent.mito", "time", "AV"]
OUT_H5AD = ROOT / "data/processed/cardiomyocytes.h5ad"
OUT_CELLS = ROOT / "config/cells.csv.gz"


def strip(adata: ad.AnnData) -> ad.AnnData:
    counts = adata.layers["counts"]
    obs = adata.obs[[c for c in KEEP_OBS if c in adata.obs.columns]].copy()
    return ad.AnnData(X=counts, obs=obs, var=pd.DataFrame(index=adata.var_names.copy()))


def cell_list(adata: ad.AnnData) -> pd.DataFrame:
    return pd.DataFrame({"cell": adata.obs_names, "library": adata.obs["orig.ident"].to_numpy()})


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main(source: str):
    out = strip(ad.read_h5ad(source))
    OUT_H5AD.parent.mkdir(parents=True, exist_ok=True)
    out.write_h5ad(OUT_H5AD, compression="gzip")
    cell_list(out).to_csv(OUT_CELLS, index=False, compression={"method": "gzip", "mtime": 0})
    summary = {
        "cells": out.n_obs,
        "genes": out.n_vars,
        "libraries": int(out.obs["orig.ident"].nunique()),
        "obs_columns": list(out.obs.columns),
        "cells_csv_sha256": sha256(OUT_CELLS),
        "h5ad_sha256": sha256(OUT_H5AD),
    }
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main(sys.argv[1])
