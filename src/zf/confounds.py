"""Per-cell technical flags from the full 10x libraries: doublets and ambient RNA.

The cardiomyocyte file holds cardiomyocytes only, but doublets and ambient RNA are properties
of the whole library. So both are estimated on every cell of the 22 filtered GEO matrices
(data/raw/, scripts/get_data.sh) and then mapped to the cardiomyocyte barcodes.

- Doublets: scanpy pp.scrublet (the doublet covariate of config/checks.yaml) on each library
  separately, all cells, expected rate 0.8% per 1,000 cells recovered (10x user guide).
- Ambient RNA: the profile of each library is the pooled expression of its non-cardiomyocyte
  cells. Genes specific to non-cardiomyocytes (pooled share at least 10 times their share in
  cardiomyocytes; stricter ratios drop the abundant soup genes, because ambient RNA itself
  puts them into cardiomyocytes at a few percent) should be absent from a clean cardiomyocyte,
  so their counts in one are contamination. The contamination fraction of a cell is the
  non-negative least-squares fit of those counts to the library profile scaled by its depth.
  Noisy below ~1,000 UMIs: one soup-gene count can give 0.1 to 0.3, so use ambient_frac as a
  continuous covariate and ambient_flag only as a screen.

Flags only: no cell is removed. Output: data/derived/cm_flags.tsv (one row per cardiomyocyte)
and data/derived/library_flags.tsv.
"""

import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp

from zf import io, provenance

RAW = io.DATA / "raw"
DERIVED = io.DATA / "derived"
RESULTS = io.ROOT / "results" / "hk"
SPECIFIC_RATIO = 10  # non-CM share / CM share for a gene to count as non-CM specific
MIN_SHARE = 1e-4  # pooled non-CM share below which a gene is too rare to estimate from
MAX_GENES = 300
FLAG_MADS = 3  # ambient outlier: above the library median by this many scaled MADs
FLAG_FLOOR = 0.1  # ...and above 10% contamination; most libraries have median and MAD near 0
MAX_CALL_RATIO = 4  # auto doublet threshold is rejected if it calls > 4x the expected rate


def ambient_fraction(x_g, total, profile_g):
    """1-D NNLS: x_ig ~ rho_i * total_i * p_g over the non-CM genes g, rho_i >= 0.

    x_g: cells x genes counts on the non-CM genes; total: counts per cell; profile_g: the
    library's ambient share of each of those genes (fraction of all its non-CM counts).
    """
    x_g = sp.csr_matrix(x_g)
    num = np.asarray(x_g @ profile_g).ravel()
    den = total * float(profile_g @ profile_g)
    return np.clip(num / den, 0, 1)


def _mad_flag(v):
    med = np.median(v)
    mad = 1.4826 * np.median(np.abs(v - med))
    return v > max(med + FLAG_MADS * mad, FLAG_FLOOR)


def _read(lib):
    a = sc.read_10x_h5(RAW / f"{lib}_filtered_matrix.h5")
    a.var_names_make_unique()
    a.obs_names = [f"{lib}_{b.split('-')[0]}" for b in a.obs_names]
    return a


def run(seed: int = 0) -> pd.DataFrame:
    cm = io.load(flags=False)
    libs = cm.obs["orig_ident"].astype(str)
    full = {lib: _read(lib) for lib in libs.unique()}

    # Non-CM-specific genes, chosen once over all libraries.
    genes = cm.var_names
    for a in full.values():
        genes = genes.intersection(a.var_names)
    cm_sum = np.zeros(len(genes))
    non_sum = np.zeros(len(genes))
    for lib, a in full.items():
        is_cm = a.obs_names.isin(cm.obs_names[libs == lib])
        x = sp.csr_matrix(a[:, genes].X)
        cm_sum += np.asarray(x[is_cm].sum(0)).ravel()
        non_sum += np.asarray(x[~is_cm].sum(0)).ravel()
    cm_share, non_share = cm_sum / cm_sum.sum(), non_sum / non_sum.sum()
    ok = (non_share > MIN_SHARE) & (non_share >= SPECIFIC_RATIO * cm_share)
    order = np.argsort(-non_share * ok)[: min(MAX_GENES, int(ok.sum()))]
    specific = genes[order]

    rows, lib_rows, scores = [], [], {}
    for lib, a in full.items():
        mine = cm.obs_names[libs == lib]
        missing = mine.difference(a.obs_names)
        if len(missing):
            raise ValueError(f"{lib}: {len(missing)} cardiomyocyte barcodes not in the matrix")
        is_cm = a.obs_names.isin(mine)
        x = sp.csr_matrix(a.X)
        n = a.n_obs
        rate = 0.008 * n / 1000
        sc.pp.scrublet(a, expected_doublet_rate=rate, random_state=seed)
        score = a.obs["doublet_score"].to_numpy()
        threshold = a.uns["scrublet"].get("threshold", np.nan)
        scores[lib] = (score, is_cm)

        xs = sp.csr_matrix(a[:, specific].X)
        profile = np.asarray(xs[~is_cm].sum(0)).ravel() / x[~is_cm].sum()
        total = np.asarray(x[is_cm].sum(1)).ravel()
        rho = ambient_fraction(xs[is_cm], total, profile)
        rows.append(
            pd.DataFrame(
                {
                    "doublet_score": score[is_cm],
                    "ambient_frac": rho,
                    "ambient_umi": np.asarray(xs[is_cm].sum(1)).ravel(),
                    "ambient_flag": _mad_flag(rho),
                },
                index=a.obs_names[is_cm],
            )
        )
        lib_rows.append(
            {
                "library": lib,
                "cells_full": n,
                "cells_cm": int(is_cm.sum()),
                "doublet_expected_rate": rate,
                "doublet_threshold_auto": threshold,
                "ambient_share_specific": float(profile.sum()),
                "ambient_frac_median": float(np.median(rho)),
            }
        )
    # scrublet's automatic threshold fails on a unimodal score histogram (seen on the shallow
    # Hr11): it then calls most cells doublets. Such a library gets the median threshold of the
    # libraries whose automatic call rate is plausible.
    lt = pd.DataFrame(lib_rows).set_index("library")
    lt["auto_rate"] = [
        (scores[lib][0] > lt.loc[lib, "doublet_threshold_auto"]).mean() for lib in lt.index
    ]
    # A missing automatic threshold (NaN) would silently call no doublets: treat it as failed.
    sane = lt["doublet_threshold_auto"].notna() & (
        lt["auto_rate"] <= MAX_CALL_RATIO * lt["doublet_expected_rate"]
    )
    lt["doublet_threshold"] = lt["doublet_threshold_auto"].where(
        sane, lt.loc[sane, "doublet_threshold_auto"].median()
    )
    lt["doublet_threshold_source"] = np.where(sane, "auto", "median_of_sane_libraries")
    calls = []
    for lib, (score, is_cm) in scores.items():
        call = score > lt.loc[lib, "doublet_threshold"]
        lt.loc[lib, "doublet_rate_full"] = call.mean()
        lt.loc[lib, "doublet_rate_cm"] = call[is_cm].mean()
        calls.append(call[is_cm])
    flags = pd.concat(rows)
    flags.insert(1, "doublet_flag", np.concatenate(calls))
    flags = flags.loc[cm.obs_names]
    DERIVED.mkdir(parents=True, exist_ok=True)
    flags.to_csv(DERIVED / "cm_flags.tsv", sep="\t")
    lt.to_csv(DERIVED / "library_flags.tsv", sep="\t")
    pd.Series(specific, name="gene").to_csv(DERIVED / "ambient_genes.tsv", index=False)
    RESULTS.mkdir(parents=True, exist_ok=True)
    lt.to_csv(RESULTS / "library_flags.tsv", sep="\t")
    provenance.write_manifest(
        io.ROOT / "results" / "runs" / "hk-confounds",
        command="uv run python -m zf.confounds",
        inputs=[io.LABEL_FREE, *(RAW / f"{lib}_filtered_matrix.h5" for lib in full)],
        params={"seed": seed, "specific_ratio": SPECIFIC_RATIO, "max_call_ratio": MAX_CALL_RATIO},
        outcomes={},
    )
    return flags


if __name__ == "__main__":
    f = run()
    print(f.describe().T)
    print(pd.read_csv(DERIVED / "library_flags.tsv", sep="\t").to_string())
