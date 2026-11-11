"""Leave-one-heart-out stability of the ventricle 3 dpi up genes, with size-free metrics.

Why: the Jaccard of top-30 lists in pseudobulk_de.py section 4 is capped when the lists differ in
size (15 vs 30 genes caps it at 0.50). Here each leave-one-heart-out (LOHO) fit is compared with
the full fit by rank and by containment, and the same numbers are computed for label shuffles of
the 8 hearts (4 vs 4, true split and its mirror excluded), each shuffle against its own full fit.

Same model as section 4: ventricle, uninjured vs dpi_3, raw counts, ~condition. Lists as there:
degs() (padj < 0.05, |log2FC| > 1, dissociation/ambient genes out), up only, top 30 by padj.
`*_rank` variants use the top 30 / top 100 by Wald stat with no threshold, so they are defined
for every fit, including shuffles whose full fit has few or no DEGs.

Usage: PYTHONPATH=src ZF_DATA=<data> .venv/bin/python scripts/pseudobulk_loho_rank.py
"""

from itertools import combinations

import numpy as np
import pandas as pd
import scipy.sparse as sp
import yaml
from pseudobulk_de import MIN_CELLS, OUT, TOP_N, degs, deseq, pseudobulk
from run_protocol import ROOT, axis_genes, dissoc_genes, strata

from zf import io, protocol

N_SHUFFLE, SEED, WIDE = 10, 20261107, 100


def tops(res, exclude):
    """Full DEG table, DEG up list (top 30 by padj), and rank-based top-30 up (by stat)."""
    d = degs(res, exclude)
    up = d[d["log2FoldChange"] > 0].sort_values("padj").index[:TOP_N].tolist()
    r = res.dropna(subset=["stat"])
    r = r[~r.index.isin(exclude)].sort_values("stat", ascending=False)
    return d, up, r.index[:TOP_N].tolist(), r.index[:WIDE].tolist()


def spearman(a, b, genes):
    g = [x for x in genes if x in a.index and x in b.index]
    g = [x for x in g if pd.notna(a.loc[x, "stat"]) and pd.notna(b.loc[x, "stat"])]
    return a.loc[g, "stat"].corr(b.loc[g, "stat"], method="spearman") if len(g) > 2 else np.nan


def compare(full, sub, exclude):
    """Metrics of one left-out fit `sub` against its full fit `full` (results tables)."""
    d, up, up_rank, _ = tops(full, exclude)
    ds, ups, ups_rank, wide = tops(sub, exclude)
    full_up = d[d["log2FoldChange"] > 0].index
    lfc = sub["log2FoldChange"].reindex(full_up)
    overlap = lambda x, y: len(set(x) & set(y)) / min(len(x), len(y)) if x and y else np.nan  # noqa: E731
    return {
        "degs_full": len(d),
        "up_full": len(full_up),
        "degs_sub": len(ds),
        "up_sub": len(ups),
        "spearman_all": spearman(full, sub, full.index),
        "spearman_degs": spearman(full, sub, d.index),
        "sign_agree": float((lfc > 0).mean()) if len(full_up) else np.nan,
        "contain_top100": len(set(up) & set(wide)) / len(up) if up else np.nan,
        "overlap_coef_top30": overlap(up, ups),
        "jaccard_top30": len(set(up) & set(ups)) / len(set(up) | set(ups)) if up or ups else np.nan,
        "contain_top100_rank": len(set(up_rank) & set(wide)) / len(up_rank),
        "overlap_coef_top30_rank": overlap(up_rank, ups_rank),
    }


def loho(x, meta, exclude, case, tag, full):
    """One row per left-out heart; `case` = the hearts labelled injured in this labelling."""
    rows = []
    for h in sorted(meta["heart"].unique()):
        vv = meta["heart"].ne(h)
        m = meta[vv].assign(condition=np.where(meta[vv]["heart"].isin(case), "case", "ctrl"))
        if m.groupby("condition")["heart"].nunique().min() < 2:
            continue
        res = deseq(x[vv], m, "~condition", ["condition", "case", "ctrl"])
        rows.append(
            {"kind": tag[0], "shuffle": tag[1], "left_out": h,
             "case_hearts": ";".join(sorted(case)), **compare(full, res, exclude)}
        )
        print(rows[-1], flush=True)
    return rows


def main():
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
    exclude = set(dissoc_genes(a.var_names, scored)) | ambient

    x, meta = pseudobulk(a, counts, cham, keep)
    v = meta["n_cells"].ge(MIN_CELLS) & meta["condition"].isin(["uninjured", "dpi_3"])
    v &= meta["chamber"].eq("ventricle")
    x, meta = x[v], meta[v]
    hearts = sorted(meta["heart"])
    assert len(hearts) == 8 and meta["heart"].is_unique, hearts
    true_case = set(meta.loc[meta["condition"].eq("dpi_3"), "heart"])
    assert len(true_case) == 4, true_case

    def full_fit(case):
        m = meta.assign(condition=np.where(meta["heart"].isin(case), "case", "ctrl"))
        return deseq(x, m, "~condition", ["condition", "case", "ctrl"])

    full = full_fit(true_case)
    rows = loho(x, meta, exclude, true_case, ("observed", 0), full)

    # 35 balanced splits; drop the true one (its mirror is the same split) -> 34 candidates.
    mirror = set(hearts) - true_case
    cands = [set(c) for c in combinations(hearts, 4) if set(c) not in (true_case, mirror)]
    cands = [c for c in cands if min(c) < min(set(hearts) - c)]  # one labelling per split
    assert len(cands) == 34, len(cands)
    rng = np.random.default_rng(SEED)
    for i, k in enumerate(rng.choice(len(cands), N_SHUFFLE, replace=False), 1):
        case = cands[k] if rng.random() < 0.5 else set(hearts) - cands[k]  # random orientation
        rows += loho(x, meta, exclude, case, ("shuffle", i), full_fit(case))
    rows_df = pd.DataFrame(rows)
    rows_df.to_csv(OUT / "pseudobulk_loho_rank.csv", index=False)

    metrics = [c for c in rows_df if c not in ("kind", "shuffle", "left_out", "case_hearts")]
    summ = rows_df.groupby("kind")[metrics].median().T
    print(summ.to_string())
    print("shuffles with >=1 full up DEG:",
          rows_df[rows_df["kind"].eq("shuffle")].groupby("shuffle")["up_full"].first().gt(0).sum())


if __name__ == "__main__":
    main()
