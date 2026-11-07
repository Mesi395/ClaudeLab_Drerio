"""Calibration of the ventricle 3 dpi vs Ctrl pseudobulk model by label permutation.

Octo 2026-11-07 16:43. 4 vs 4 hearts -> all C(8,4) = 70 label assignments (observed included).
Each one refits DESeq2 ~condition; we record the DEG count (same rule and exclusions as
pseudobulk_de.py) and the mean Wald stat of every marker set (mitochondrial genes dropped).
Negative control: 30 dpi vs Ctrl (many DEGs there would point at batch, not injury).

Usage: uv run python scripts/pseudobulk_permutation.py
"""

import os
from itertools import combinations
from multiprocessing import Pool

import numpy as np
import pandas as pd
import scipy.sparse as sp
import yaml
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
from pseudobulk_de import MIN_CELLS, OUT, degs, pseudobulk
from pseudobulk_set_test import MITO, marker_sets
from run_protocol import ROOT, axis_genes, dissoc_genes, strata

from zf import io, protocol

X = META = EXCLUDE = None


def run(case_samples):
    meta = META.assign(condition=np.where(META.index.isin(case_samples), "case", "ctrl"))
    x = X.loc[:, X.sum(0) >= 10]
    dds = DeseqDataSet(counts=x, metadata=meta, design="~condition", quiet=True, n_cpus=1)
    dds.deseq2()
    ds = DeseqStats(dds, contrast=["condition", "case", "ctrl"], quiet=True, n_cpus=1)
    ds.summary()
    res = ds.results_df
    d = degs(res, EXCLUDE)
    r = res[~res.index.str.match(MITO)].dropna(subset=["stat"])
    row = {"case": ";".join(sorted(case_samples)), "degs": len(d),
           "up": int((d["log2FoldChange"] > 0).sum())}
    for name, g in marker_sets().items():
        row[f"mean_stat_{name}"] = r.loc[r.index.intersection(g), "stat"].mean()
    return row


def main():
    global X, META, EXCLUDE
    a = io.load()
    counts = sp.csr_matrix(a.layers["counts"])
    st = strata(a)
    cham = np.where(st["ventricle"], "ventricle", np.where(st["atrium"], "atrium", "uncertain"))
    controls = yaml.safe_load((ROOT / "config" / "axes.yaml").read_text())["controls"]
    fib = [g["symbol"] for g in next(c for c in controls if c["name"] == "ambient_fibroblast")["genes"]]
    veto = protocol._veto(a.obs, counts, pd.Index(a.var_names), fib)["vetoed"].to_numpy()
    keep = pd.Series(~veto, index=a.obs_names)
    scored = {g for gs in axis_genes().values() for g in gs}
    ambient = set(pd.read_csv(io.DATA / "derived" / "ambient_genes.tsv")["gene"])
    EXCLUDE = set(dissoc_genes(a.var_names, scored)) | ambient
    x, meta = pseudobulk(a, counts, cham, keep)
    del a, counts  # ponytail: forked workers need only the pseudobulk matrix
    ok = meta["n_cells"].ge(MIN_CELLS) & meta["chamber"].eq("ventricle")
    rows = []
    for t in ["dpi_3", "dpi_30"]:
        v = ok & meta["condition"].isin(["uninjured", t])
        X, META = x[v], meta[v][["heart", "condition"]]
        observed = tuple(META.index[META["condition"].eq(t)])
        assigns = [observed] if t == "dpi_30" else list(combinations(META.index, len(observed)))
        with Pool(min(len(assigns), 12)) as p:  # ponytail: 12 keeps RAM on the 31 GB ROG
            out = p.map(run, assigns)
        for o, c in zip(out, assigns):
            o.update(contrast=f"ventricle_{t}_vs_ctrl", n_ctrl=int((~META.index.isin(c)).sum()),
                     n_case=len(c), observed=set(c) == set(observed))
        rows += out
        print(t, len(out), flush=True)
    df = pd.DataFrame(rows)
    perm = df[df["contrast"].eq("ventricle_dpi_3_vs_ctrl")]
    obs = perm[perm["observed"]].iloc[0]
    summary = []
    for col in ["degs"] + [c for c in df if c.startswith("mean_stat_")]:
        summary.append({"contrast": "ventricle_dpi_3_vs_ctrl", "case": "PERMUTATION_SUMMARY",
                        "statistic": col, "observed_value": obs[col],
                        "perm_median": perm[col].median(), "perm_max": perm[col].max(),
                        "perm_p_ge": (perm[col] >= obs[col]).mean(), "n_perm": len(perm)})
    pd.concat([df, pd.DataFrame(summary)]).to_csv(OUT / "pseudobulk_permutation.csv", index=False)
    print(pd.DataFrame(summary).to_string())
    print(df[df["contrast"].eq("ventricle_dpi_30_vs_ctrl")].T.to_string())


if __name__ == "__main__":
    main()
