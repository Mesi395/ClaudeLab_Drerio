"""Data loading: counts, metadata and per-cell technical covariates.

`sanitize_h5ad` is the leakage guard (study cm-state-axes, requirement R5.3). The supplied
Hu et al. 2022 files embed the authors' cell labels. This function opens the source with h5py
and reads **only** allowlisted elements, so the values of dropped columns are never loaded.

Usage:
    uv run python -m zf.io data/hu2022/<in>.h5ad data/clean/<out>.h5ad
"""

from __future__ import annotations

import json
import subprocess
import sys
from datetime import datetime
from pathlib import Path

import anndata as ad
import h5py
import pandas as pd
from anndata.io import read_elem

ALLOWED_OBS = (
    "orig.ident",
    "nCount_RNA",
    "nFeature_RNA",
    "percent.mito",
    "time",
    "morphine",
    "AV",
)


def _commit() -> str:
    try:
        return subprocess.check_output(["git", "rev-parse", "HEAD"], text=True).strip()
    except Exception:  # noqa: BLE001
        return "unknown"


def sanitize_h5ad(src: str | Path, dst: str | Path) -> dict:
    """Copy `src` to `dst`, keeping only the allowlisted obs columns, counts, var names and stats.

    X in the output is the integer `counts` layer (the authors' normalisation is dropped).
    Raises KeyError if an allowlisted column or the counts layer is missing.
    Returns a provenance record (names only, never values of dropped columns).
    """
    src, dst = Path(src), Path(dst)
    with h5py.File(src, "r") as f:
        obs_grp = f["obs"]
        idx_key = obs_grp.attrs["_index"]
        obs_names = read_elem(obs_grp[idx_key])
        missing = [c for c in ALLOWED_OBS if c not in obs_grp]
        if missing:
            raise KeyError(f"allowlisted obs columns missing: {missing}")
        obs = pd.DataFrame({c: read_elem(obs_grp[c]) for c in ALLOWED_OBS})
        obs.index = pd.Index(obs_names, name=None).astype(str)
        dropped_obs = [c for c in obs_grp.attrs.get("column-order", []) if c not in ALLOWED_OBS]

        var = read_elem(f["var"])  # gene names and per-gene stats only
        if "layers" not in f or "counts" not in f["layers"]:
            raise KeyError("layers['counts'] missing")
        counts = read_elem(f["layers/counts"])
        dropped_other = {
            "uns": list(f["uns"].keys()) if "uns" in f else [],
            "obsm": list(f["obsm"].keys()) if "obsm" in f else [],
            "obsp": list(f["obsp"].keys()) if "obsp" in f else [],
            "varm": list(f["varm"].keys()) if "varm" in f else [],
            "raw": "raw" in f,
        }

    out = ad.AnnData(X=counts.copy(), obs=obs, var=var)
    out.layers["counts"] = counts
    dst.parent.mkdir(parents=True, exist_ok=True)
    out.write_h5ad(dst, compression="gzip")
    record = {
        "source": str(src),
        "output": str(dst),
        "commit": _commit(),
        "command": " ".join(sys.argv),
        "time": datetime.now().astimezone().isoformat(timespec="seconds"),
        "n_obs": int(out.n_obs),
        "n_vars": int(out.n_vars),
        "kept_obs": list(ALLOWED_OBS),
        "dropped_obs_names": dropped_obs,
        "dropped_other_names": dropped_other,
        "X": "layers['counts'] copied (authors' X dropped)",
    }
    Path(str(dst) + ".json").write_text(json.dumps(record, indent=2))
    return record


if __name__ == "__main__":
    rec = sanitize_h5ad(sys.argv[1], sys.argv[2])
    print(json.dumps({k: rec[k] for k in ("output", "n_obs", "n_vars", "kept_obs")}, indent=1))
