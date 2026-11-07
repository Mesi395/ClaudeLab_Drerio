"""Pseudobulk DE per heart x chamber with pydeseq2, and a data-derived injury gene set.

Owner directive 2026-11-07 15:11: "Downsampling robustness checks and sample-level pseudobulking
are very important. The injury-responsiveness score could be derived from a DESeq2 DEG list
from sample pseudobulks of ventricle subsets." Hearts are the replicates (Squair 2021).

1. Pseudobulk = sum of raw counts per heart x chamber over non-vetoed cardiomyocytes; chamber
   from GEO for a/v libraries, the per-cell call for whole-heart ones; < MIN_CELLS dropped.
2. Ventricle: ~condition (3 dpi vs uninjured, 7 dpi vs uninjured). Both chambers:
   ~heart + chamber + chamber:condition; the interaction is the ventricle-specific effect.
3. Clean: drop dissociation/heat-shock and ambient (non-CM) genes; padj < 0.05, |log2FC| > 1;
   top TOP_N up genes = the data-derived injury set.
4. Leave-one-heart-out DEG sets and their Jaccard overlap.
5. Thinned (600 UMIs per cell before summing) re-run: logFC correlation and DEG overlap.
6. Overlap with the literature injury set.

Usage: uv run python scripts/pseudobulk_de.py
"""

import numpy as np
import pandas as pd
import scipy.sparse as sp
import yaml
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
from run_protocol import ROOT, axis_genes, dissoc_genes, strata

from zf import io, protocol

OUT = ROOT / "results" / "hk"
MIN_CELLS = 30
PADJ, LFC, TOP_N = 0.05, 1.0, 30


def pseudobulk(a, counts, cham, keep):
    obs = a.obs
    key = obs["heart"].astype(str) + "|" + pd.Series(cham, index=obs.index)
    key = key.where(keep & (cham != "uncertain"))
    groups = key.dropna().unique()
    rows, meta = [], []
    for g in sorted(groups):
        m = (key == g).to_numpy()
        heart, ch = g.split("|")
        meta.append(
            {
                "sample": g,
                "heart": heart,
                "chamber": ch,
                "condition": str(obs.loc[m, "condition"].iloc[0]),
                "n_cells": int(m.sum()),
            }
        )
        rows.append(np.asarray(counts[m].sum(0)).ravel())
    meta = pd.DataFrame(meta).set_index("sample")
    x = pd.DataFrame(np.vstack(rows).astype(int), index=meta.index, columns=a.var_names)
    return x, meta


def fit(x, meta, design):
    dds = DeseqDataSet(counts=x.loc[:, x.sum(0) >= 10], metadata=meta, design=design, quiet=True)
    dds.deseq2()
    return dds


def stats(dds, contrast):
    """contrast: [factor, case, ref], or one design-matrix column name (interaction term)."""
    if isinstance(contrast, str):
        cols = list(dds.obsm["design_matrix"].columns)
        contrast = (np.array(cols) == contrast).astype(float)
    ds = DeseqStats(dds, contrast=contrast, quiet=True)
    ds.summary()
    return ds.results_df


def deseq(x, meta, design, contrast):
    return stats(fit(x, meta, design), contrast)


def degs(res, exclude):
    r = res.dropna(subset=["padj"])
    r = r[(r["padj"] < PADJ) & (r["log2FoldChange"].abs() > LFC) & ~r.index.isin(exclude)]
    return r


def jaccard(a, b):
    a, b = set(a), set(b)
    return len(a & b) / len(a | b) if a | b else np.nan


def main():
    a = io.load()
    counts = sp.csr_matrix(a.layers["counts"])
    st = strata(a)
    cham = np.where(st["ventricle"], "ventricle", np.where(st["atrium"], "atrium", "uncertain"))
    controls = yaml.safe_load((ROOT / "config" / "axes.yaml").read_text())["controls"]
    fib = [
        g["symbol"] for g in next(c for c in controls if c["name"] == "ambient_fibroblast")["genes"]
    ]
    veto = protocol._veto(a.obs, counts, pd.Index(a.var_names), fib)["vetoed"].to_numpy()
    keep = pd.Series(~veto, index=a.obs_names)
    genes = axis_genes()
    scored = {g for gs in genes.values() for g in gs}
    dissoc = set(dissoc_genes(a.var_names, scored))
    ambient = set(pd.read_csv(io.DATA / "derived" / "ambient_genes.tsv")["gene"])
    exclude = dissoc | ambient

    thinned, _ = protocol.thin(counts, protocol.T, np.random.default_rng(protocol.SEED + 1))
    runs = {"raw": counts, "thinned600": sp.csr_matrix(thinned.astype(np.int64))}
    summary, sets = [], {}
    for run, cm in runs.items():
        x, meta = pseudobulk(a, cm, cham, keep)
        dropped = meta[meta["n_cells"] < MIN_CELLS]
        x, meta = x.loc[meta["n_cells"] >= MIN_CELLS], meta[meta["n_cells"] >= MIN_CELLS]
        if run == "raw":
            meta.assign(dropped=False).reset_index().to_csv(
                OUT / "pseudobulk_samples.csv", index=False
            )
            dropped.reset_index().to_csv(OUT / "pseudobulk_dropped.csv", index=False)
        for t in ["dpi_3", "dpi_7"]:
            arms = meta["condition"].isin(["uninjured", t])
            v = arms & meta["chamber"].eq("ventricle")
            res = deseq(x[v], meta[v], "~condition", ["condition", t, "uninjured"])
            res.to_csv(OUT / f"pseudobulk_de_{run}_ventricle_{t}_vs_ctrl.csv.gz")
            d = degs(res, exclude)
            up = d[d["log2FoldChange"] > 0].sort_values("padj").index[:TOP_N].tolist()
            sets[(run, "ventricle", t)] = (res, d.index.tolist(), up)
            n_arm = meta[v].groupby("condition")["heart"].nunique().to_dict()
            row = {
                "run": run,
                "model": "ventricle ~condition",
                "contrast": f"{t}_vs_uninjured",
                "hearts_ctrl": n_arm.get("uninjured", 0),
                "hearts_case": n_arm.get(t, 0),
                "degs": len(d),
                "up": int((d["log2FoldChange"] > 0).sum()),
                "top_up": ";".join(up),
            }
            summary.append(row)
            print(row, flush=True)
            # Interaction: hearts with both chambers.
            both = meta[arms].groupby("heart")["chamber"].nunique()
            hb = both[both == 2].index
            m2 = arms & meta["heart"].isin(hb)
            mm = meta[m2].copy()
            # One explicit chamber:condition column: formulaic's own interaction coding adds a
            # column per chamber level, which is collinear with heart (condition is nested in it).
            mm["vent_injured"] = (mm["chamber"].eq("ventricle") & mm["condition"].eq(t)).astype(
                float
            )
            try:
                dds = fit(x[m2], mm, "~heart + chamber + vent_injured")
                term = "vent_injured"
                res2 = stats(dds, term)
            except Exception as e:  # noqa: BLE001 - reported in the summary row, not silent
                summary.append(
                    {
                        "run": run,
                        "model": "interaction",
                        "contrast": t,
                        "top_up": f"not run: {e!r}"[:300],
                    }
                )
            else:
                res2.to_csv(OUT / f"pseudobulk_de_{run}_interaction_{t}.csv.gz")
                d2 = degs(res2, exclude)
                up2 = d2[d2["log2FoldChange"] > 0].sort_values("padj").index[:TOP_N].tolist()
                sets[(run, "interaction", t)] = (res2, d2.index.tolist(), up2)
                n2 = mm.groupby("condition", observed=True)["heart"].nunique().to_dict()
                summary.append(
                    {
                        "run": run,
                        "model": f"~heart+chamber+chamber:condition [{term}]",
                        "contrast": f"{t}_vs_uninjured",
                        "hearts_ctrl": n2.get("uninjured", 0),
                        "hearts_case": n2.get(t, 0),
                        "degs": len(d2),
                        "up": int((d2["log2FoldChange"] > 0).sum()),
                        "top_up": ";".join(up2),
                    }
                )
            print(summary[-1], flush=True)

    s = pd.DataFrame(summary)
    s.to_csv(OUT / "pseudobulk_de_summary.csv", index=False)

    # 5. Thinning robustness of the DE itself.
    rob = []
    for (run, model, t), (res, d, up) in sets.items():
        if run != "raw" or ("thinned600", model, t) not in sets:
            continue
        rt, dt, upt = sets[("thinned600", model, t)]
        j = res.index.intersection(rt.index)
        rho = pd.Series(res.loc[j, "log2FoldChange"]).corr(
            rt.loc[j, "log2FoldChange"], method="spearman"
        )
        rob.append(
            {
                "model": model,
                "contrast": t,
                "lfc_spearman_raw_vs_thinned": rho,
                "deg_jaccard": jaccard(d, dt),
                "top_up_jaccard": jaccard(up, upt),
                "degs_raw": len(d),
                "degs_thinned": len(dt),
            }
        )
    pd.DataFrame(rob).to_csv(OUT / "pseudobulk_thinning_robustness.csv", index=False)

    # 6. Corroboration with the literature injury set.
    lit = set(genes["injury_responsive"])
    cor = [
        {
            "model": m,
            "contrast": t,
            "run": r,
            "data_up": len(up),
            "literature": len(lit),
            "overlap": ";".join(sorted(lit & set(up))),
            "jaccard": jaccard(lit, up),
        }
        for (r, m, t), (_, _, up) in sets.items()
    ]
    pd.DataFrame(cor).to_csv(OUT / "pseudobulk_vs_literature.csv", index=False)

    # 4. Leave-one-heart-out on the raw ventricle 3 dpi model: DEG-set stability.
    x, meta = pseudobulk(a, counts, cham, keep)
    x, meta = x.loc[meta["n_cells"] >= MIN_CELLS], meta[meta["n_cells"] >= MIN_CELLS]
    v = meta["condition"].isin(["uninjured", "dpi_3"]) & meta["chamber"].eq("ventricle")
    full_up = sets[("raw", "ventricle", "dpi_3")][2]
    lo = []
    for h in meta[v]["heart"].unique():
        vv = v & meta["heart"].ne(h)
        if meta[vv].groupby("condition")["heart"].nunique().min() < 2:
            continue
        res = deseq(x[vv], meta[vv], "~condition", ["condition", "dpi_3", "uninjured"])
        d = degs(res, exclude)
        up = d[d["log2FoldChange"] > 0].sort_values("padj").index[:TOP_N].tolist()
        lo.append(
            {
                "left_out": h,
                "degs": len(d),
                "top_up_jaccard_vs_full": jaccard(up, full_up),
                "top_up": ";".join(up),
            }
        )
        print(lo[-1], flush=True)
    pd.DataFrame(lo).to_csv(OUT / "pseudobulk_loho_dpi3.csv", index=False)
    print(s.to_string())


if __name__ == "__main__":
    main()
