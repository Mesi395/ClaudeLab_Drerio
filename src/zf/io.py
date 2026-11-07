"""Data loading: counts, metadata and per-cell technical covariates.

`build()` reads the team's stripped file from `scripts/make_agent_input.py` (raw counts plus an
allowlist of technical columns; the labelled original stays with its owner), adds per-cell
covariates and the heart mapping, and writes the analysis copy to data/. Analysis code calls
`load()`, which opens that copy and joins the per-cell flags written by `zf.confounds`.
"""

import os
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

from zf import config

ROOT = Path(__file__).resolve().parents[2]
DATA = Path(os.environ.get("ZF_DATA", ROOT / "data"))
SOURCE = Path(os.environ.get("ZF_SOURCE_H5AD", DATA / "processed" / "cardiomyocytes.h5ad"))
LABEL_FREE = DATA / "cm_labelfree.h5ad"
FLAGS = DATA / "derived" / "cm_flags.tsv"

# Metadata kept from the stripped file (the same allowlist as make_agent_input.KEEP_OBS).
KEEP = {
    "orig.ident": "orig_ident",
    "nCount_RNA": "ncount_source",
    "nFeature_RNA": "nfeature_source",
    "percent.mito": "pct_mito_source",
    "time": "time_source",
    "AV": "av_source",
}
TIME = {"Ctrl": "uninjured", "3dpi": "dpi_3", "7dpi": "dpi_7", "30dpi": "dpi_30"}
AV = {"Wholeheart": "whole", "Atrium": "atrium", "Ventricle": "ventricle"}
# Two libraries from one heart (GEO sample records; teammate review on PR #8).
SAME_HEART = {"Hr2a": "Hr2", "Hr2b": "Hr2", "H8a": "H8", "H8v": "H8"}
SAME_HEART |= {"Hr6a": "Hr6", "Hr6v": "Hr6", "Hr7a": "Hr7", "Hr7v": "Hr7"}
# GRCz11 / Ensembl 92 names: vmhc is myh7 (and vmhcl is myh7l), amhc is myh6.
VENTRICULAR, ATRIAL = ["myh7", "myh7l"], ["myh6"]


def heart_of(library: str) -> str:
    return SAME_HEART.get(library, library)


def library_table() -> pd.DataFrame:
    ds = config.load_dataset(ROOT / "config" / "dataset.yaml")
    t = pd.DataFrame([lib.model_dump() for lib in ds.libraries]).set_index("library")
    t["heart"] = [heart_of(lib) for lib in t.index]
    # Chamber is known only for libraries GEO titles as atrium or ventricle.
    t["chamber"] = t["chamber"].replace({"whole": "unknown"})
    t["morphine_study"] = t["study"].eq("GSE158919")
    return t


def _va_call(a: ad.AnnData) -> pd.DataFrame:
    """Per-cell ventricular minus atrial myosin expression (log1p CP10k of raw counts)."""
    x = sp.csr_matrix(a.layers["counts"])
    scale = 1e4 / np.asarray(x.sum(1)).ravel()

    def mean_log(genes):
        cols = [a.var_names.get_loc(g) for g in genes]
        sub = sp.csr_matrix(x[:, cols].multiply(scale[:, None]))
        return np.asarray(sub.log1p().mean(1)).ravel()

    score = mean_log(VENTRICULAR) - mean_log(ATRIAL)
    call = np.where(score > 0, "ventricle", np.where(score < 0, "atrium", "unknown"))
    return pd.DataFrame({"va_score": score, "va_call": call}, index=a.obs_names)


def va_accuracy(obs: pd.DataFrame) -> pd.DataFrame:
    """Accuracy of the per-cell chamber call on libraries whose chamber GEO states."""
    known = obs[obs["chamber"].isin(["atrium", "ventricle"])]
    hit = known["va_call"].astype(str).eq(known["chamber"].astype(str))
    return hit.groupby(known["orig_ident"].astype(str), observed=True).agg(["mean", "size"])


def build(source: Path = SOURCE) -> ad.AnnData:
    a = ad.read_h5ad(source)
    libs = library_table()

    obs = a.obs[list(KEEP)].rename(columns=KEEP)
    obs["orig_ident"] = obs["orig_ident"].astype(str)
    counts = sp.csr_matrix(a.X)  # the stripped file stores raw counts in X
    total = np.asarray(counts.sum(1)).ravel()
    mito = np.asarray(a.var_names.str.startswith("mt-"))
    obs["total_counts"] = total
    obs["log_umi"] = np.log10(total)
    obs["n_genes"] = np.diff(counts.indptr)
    obs["pct_mito"] = 100 * np.asarray(counts[:, mito].sum(1)).ravel() / total
    for col in ["heart", "condition", "chamber", "chemistry", "study", "morphine_study"]:
        obs[col] = obs["orig_ident"].map(libs[col]).to_numpy()

    # Two sources for the same fact: GEO records (config/dataset.yaml) and the owner's file.
    assert obs["time_source"].astype(str).map(TIME).eq(obs["condition"]).all(), "time differs"
    geo_av = obs["chamber"].replace({"unknown": "whole"})
    assert obs["av_source"].astype(str).map(AV).eq(geo_av).all(), "chamber differs"
    assert np.allclose(obs["ncount_source"], total), "nCount_RNA is not the raw count sum"

    out = ad.AnnData(X=counts, obs=obs, var=pd.DataFrame(index=a.var_names))
    out.layers["counts"] = counts
    out.obs = out.obs.join(_va_call(out))
    for col in ["orig_ident", "heart", "condition", "chamber", "chemistry", "study", "va_call"]:
        out.obs[col] = out.obs[col].astype("category")
    DATA.mkdir(parents=True, exist_ok=True)
    out.write_h5ad(LABEL_FREE, compression="gzip")
    return out


def load(flags: bool = True) -> ad.AnnData:
    """The label-free cardiomyocytes. X and layers['counts'] hold raw counts."""
    a = ad.read_h5ad(LABEL_FREE)
    if flags:
        f = pd.read_csv(FLAGS, sep="\t", index_col=0)
        missing = a.obs_names.difference(f.index)
        if len(missing):
            raise ValueError(
                f"{len(missing)} cells have no doublet/ambient flags; run zf.confounds"
            )
        a.obs = a.obs.join(f)
    return a


if __name__ == "__main__":
    from zf import provenance

    a = build()
    print(a)
    acc = va_accuracy(a.obs)
    print(acc)
    out = ROOT / "results" / "hk"
    out.mkdir(parents=True, exist_ok=True)
    acc.to_csv(out / "va_accuracy.tsv", sep="\t")
    provenance.write_manifest(
        ROOT / "results" / "runs" / "hk-io-build",
        command="uv run python -m zf.io",
        inputs=[SOURCE],
        params={"keep": list(KEEP)},
        outcomes={},
    )
