"""C002 Probe: four candidate axes vs matched random sets and vs depth, within library.

Reading rules: research/cm-state-axes/cycles/C002/reading-rules.md (fixed before running).
uv run python scripts/c002_probe.py data/clean/hu2022_cm_allgenes_sanitized.h5ad \
    proposals/loop results
"""

import json
import os
import re
import subprocess
import sys
import time
from datetime import datetime
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import yaml
from scipy.stats import spearmanr

from zf.io import ALLOWED_OBS
from zf.score import matched_random_sets, rank_matrix, ucell_from_ranks

N_NULL = int(sys.argv[4]) if len(sys.argv) > 4 else 100
N_SG_NULL = int(sys.argv[5]) if len(sys.argv) > 5 else 50  # score_genes null sets (slower)
MAX_RANK = 1500
SEED = 0
CIRCULAR = {"ttn.2"}
# Exploratory options (C002, added after the primary run): finer null bins and removal of
# mitochondrial, NC_002333 and ribosomal genes from the null pool.
NULL_BINS = int(os.environ.get("C002_NULL_BINS", "25"))
NULL_EXCLUDE_TECH = os.environ.get("C002_NULL_EXCLUDE_TECH", "0") == "1"


def sp_rho(a, b):
    return float(spearmanr(a, b)[0])


def coherence(M):
    """Mean pairwise Spearman among columns of dense M (cells x genes); NaN-safe."""
    keep = M.std(0) > 0
    if keep.sum() < 2:
        return np.nan
    c = spearmanr(M[:, keep])[0]
    c = np.atleast_2d(c)
    iu = np.triu_indices_from(c, 1)
    return float(np.nanmean(c[iu]))


def score_genes_idx(adata, names, seed):
    sc.tl.score_genes(adata, list(names), random_state=seed, score_name="_s", use_raw=False)
    return adata.obs.pop("_s").to_numpy()


def main(src, marker_dir, out_dir):
    t0 = time.time()
    out = Path(out_dir)
    a = ad.read_h5ad(src)
    assert set(a.obs.columns) == set(ALLOWED_OBS), "input is not the sanitised file"
    counts = a.layers["counts"].tocsr().astype(np.float32)
    tot = np.asarray(counts.sum(1)).ravel()
    logtot = np.log10(tot)
    ngenes = np.asarray((counts > 0).sum(1)).ravel()
    lib = a.obs["orig.ident"].astype(str).to_numpy()
    libs = sorted(set(lib))
    a.X = counts.copy()
    sc.pp.normalize_total(a, target_sum=1e4)
    sc.pp.log1p(a)
    Xn = a.X.tocsr()
    means = np.asarray(Xn.mean(0)).ravel()
    vidx = {g: i for i, g in enumerate(a.var_names)}
    R = rank_matrix(counts)
    print(f"ranks done {time.time() - t0:.0f}s", flush=True)

    axes = {}
    for f in sorted(Path(marker_dir).glob("*.yaml")):
        spec = yaml.safe_load(f.read_text())
        genes = [m["symbol"] for m in spec["genes"]]
        axes[spec["state"]] = (f, [g for g in genes if g in vidx and g not in CIRCULAR])
    all_axis_idx = np.array(sorted({vidx[g] for _, gs in axes.values() for g in gs}))
    null_excl = all_axis_idx
    if NULL_EXCLUDE_TECH:
        tech = [i for i, g in enumerate(a.var_names) if re.match(r"^(mt-|NC-|rp[ls]|RP[LS])", g)]
        null_excl = np.union1d(all_axis_idx, tech)

    cells = pd.DataFrame({"lib": lib, "log10_total": logtot, "n_genes": ngenes}, index=a.obs_names)
    rows, summ = [], []
    for k, (state, (_, genes)) in enumerate(axes.items()):
        gi = np.array([vidx[g] for g in genes])
        uc = ucell_from_ranks(R, gi, MAX_RANK)
        sg = score_genes_idx(a, genes, SEED)
        cells[f"{state}_ucell"], cells[f"{state}_score_genes"] = uc, sg
        nulls = matched_random_sets(means, gi, N_NULL, NULL_BINS, SEED + k, exclude=null_excl)
        null_uc = np.stack([ucell_from_ranks(R, s, MAX_RANK) for s in nulls], 1)
        null_sg = np.stack(
            [score_genes_idx(a, a.var_names[s], SEED + j) for j, s in enumerate(nulls[:N_SG_NULL])],
            1,
        )
        print(f"{state}: nulls done {time.time() - t0:.0f}s", flush=True)
        for lb in libs:
            m = lib == lb
            G = Xn[m][:, gi].toarray()
            coh = coherence(G)
            null_coh = np.array([coherence(Xn[m][:, s].toarray()) for s in nulls])
            r = {"state": state, "lib": lb, "n_cells": int(m.sum()), "n_genes_axis": len(gi)}
            r["coherence"] = coh
            r["coherence_null95"] = float(np.nanpercentile(null_coh, 95))
            r["coherence_above_null"] = bool(coh > r["coherence_null95"])
            for meth, s, ns in (("ucell", uc, null_uc), ("score_genes", sg, null_sg)):
                p90 = np.percentile(s[m], 90)
                np90 = np.percentile(ns[m], 90, axis=0)
                r[f"{meth}_p90"] = float(p90)
                r[f"{meth}_p90_null95"] = float(np.percentile(np90, 95))
                r[f"{meth}_p90_above_null"] = bool(p90 > r[f"{meth}_p90_null95"])
                r[f"{meth}_median"] = float(np.median(s[m]))
                r[f"{meth}_rho_logtot"] = sp_rho(s[m], logtot[m])
                r[f"{meth}_rho_ngenes"] = sp_rho(s[m], ngenes[m])
                nr = [abs(sp_rho(ns[m, j], logtot[m])) for j in range(ns.shape[1])]
                r[f"{meth}_null_median_abs_rho_logtot"] = float(np.nanmedian(nr))
            r["rho_ucell_vs_score_genes"] = sp_rho(uc[m], sg[m])
            rows.append(r)
        df = pd.DataFrame([x for x in rows if x["state"] == state])
        s = {"state": state, "genes": ",".join(genes), "n_genes": len(genes)}
        s["libs_coherence_above_null"] = int(df["coherence_above_null"].sum())
        for meth in ("ucell", "score_genes"):
            s[f"libs_{meth}_p90_above_null"] = int(df[f"{meth}_p90_above_null"].sum())
            s[f"{meth}_median_abs_rho_logtot"] = float(df[f"{meth}_rho_logtot"].abs().median())
            s[f"{meth}_median_abs_rho_ngenes"] = float(df[f"{meth}_rho_ngenes"].abs().median())
            s[f"{meth}_null_median_abs_rho_logtot"] = float(
                df[f"{meth}_null_median_abs_rho_logtot"].median()
            )
        s["libs_methods_agree_ge0.5"] = int((df["rho_ucell_vs_score_genes"] >= 0.5).sum())
        s["median_rho_methods"] = float(df["rho_ucell_vs_score_genes"].median())
        lm = df.set_index("lib")
        for lb in ("Hr11", "Hr7a", "Hr2a", "Hr2b", "H7"):
            s[f"{lb}_ucell_median"] = float(lm.loc[lb, "ucell_median"])
            s[f"{lb}_score_genes_median"] = float(lm.loc[lb, "score_genes_median"])
        libdepth = cells.groupby("lib")["log10_total"].median()
        for meth in ("ucell", "score_genes"):
            s[f"between_lib_rho_{meth}_median_vs_depth"] = sp_rho(
                lm[f"{meth}_median"], libdepth.loc[lm.index]
            )
        summ.append(s)

    pd.DataFrame(rows).to_csv(out / "c002_axis_by_library.csv", index=False)
    pd.DataFrame(summ).to_csv(out / "c002_axis_summary.csv", index=False)
    cells.to_csv(out / "c002_scores.csv.gz")
    print(f"axes done {time.time() - t0:.0f}s", flush=True)

    # per-gene depth association vs same-bin background genes
    edges = np.quantile(means, np.linspace(0, 1, 26)[1:-1])
    bins = np.digitize(means, edges)
    rng = np.random.default_rng(SEED)
    bg = rng.choice(np.setdiff1d(np.arange(len(means)), all_axis_idx), 2000, replace=False)
    test_idx = np.concatenate([all_axis_idx, [vidx[g] for g in CIRCULAR if g in vidx], bg])
    rho = np.full((len(test_idx), len(libs)), np.nan)
    for j, lb in enumerate(libs):
        m = lib == lb
        G = Xn[m][:, test_idx].toarray()
        ok = G.std(0) > 0
        if ok.sum():
            c = spearmanr(np.column_stack([logtot[m], G[:, ok]]))[0]
            rho[ok, j] = c[0, 1:]
    med = np.nanmedian(rho, 1)
    gdf = pd.DataFrame(
        {"gene": a.var_names[test_idx], "bin": bins[test_idx], "median_rho_logtot": med}
    )
    gdf["role"] = (
        ["axis"] * len(all_axis_idx)
        + ["circular"] * (len(test_idx) - len(all_axis_idx) - len(bg))
        + ["background"] * len(bg)
    )
    bgd = gdf[gdf.role == "background"]
    lo = bgd.groupby("bin")["median_rho_logtot"].quantile(0.05)
    hi = bgd.groupby("bin")["median_rho_logtot"].quantile(0.95)
    nbg = bgd.groupby("bin").size()
    gdf["bg_p05"], gdf["bg_p95"] = gdf["bin"].map(lo), gdf["bin"].map(hi)
    gdf["bg_n"] = gdf["bin"].map(nbg)
    gdf["tracks_depth"] = (gdf.median_rho_logtot > gdf.bg_p95) | (
        gdf.median_rho_logtot < gdf.bg_p05
    )
    gdf["direction"] = np.where(
        gdf.median_rho_logtot > gdf.bg_p95,
        "up_with_depth",
        np.where(gdf.median_rho_logtot < gdf.bg_p05, "down_with_depth", ""),
    )
    gdf[gdf.role != "background"].to_csv(out / "c002_gene_depth.csv", index=False)
    gdf[gdf.role == "background"].to_csv(out / "c002_gene_depth_background.csv", index=False)

    # C9: depth-matched share of UMIs, Hr11 vs control libraries
    ctrl = cells.index[a.obs["time"].astype(str).to_numpy() == "Ctrl"]
    sel = (lib == "Hr11") | a.obs_names.isin(ctrl)
    qb = np.quantile(tot[sel], np.linspace(0, 1, 6))
    c9 = []
    for g in ("nppa", "ttn.2", "nppb"):
        share = counts[:, vidx[g]].toarray().ravel() / tot
        for b in range(5):
            inb = sel & (tot >= qb[b]) & (tot <= qb[b + 1])
            h, c = inb & (lib == "Hr11"), inb & a.obs_names.isin(ctrl)
            c9.append(
                {
                    "gene": g,
                    "depth_bin": b,
                    "tot_lo": qb[b],
                    "tot_hi": qb[b + 1],
                    "n_hr11": int(h.sum()),
                    "n_ctrl": int(c.sum()),
                    "hr11_median_share": float(np.median(share[h])),
                    "ctrl_median_share": float(np.median(share[c])),
                    "hr11_mean_share": float(share[h].mean()),
                    "ctrl_mean_share": float(share[c].mean()),
                }
            )
    pd.DataFrame(c9).to_csv(out / "c002_c9_depth_matched.csv", index=False)

    meta = {
        "commit": subprocess.check_output(["git", "rev-parse", "HEAD"], text=True).strip(),
        "dirty_tree": bool(subprocess.check_output(["git", "status", "--porcelain"], text=True)),
        "command": "uv run python " + " ".join(sys.argv),
        "input": src,
        "markers": {s: str(f) for s, (f, _) in axes.items()},
        "params": {
            "n_null_ucell": N_NULL,
            "n_null_score_genes": N_SG_NULL,
            "max_rank": MAX_RANK,
            "seed": SEED,
            "null_bins": NULL_BINS,
            "null_exclude_tech": NULL_EXCLUDE_TECH,
            "gene_depth_bins": 25,
            "excluded_circular": sorted(CIRCULAR),
            "normalisation": "normalize_total 1e4 + log1p for score_genes and coherence",
        },
        "time": datetime.now().astimezone().isoformat(timespec="seconds"),
        "runtime_s": round(time.time() - t0),
        "outputs": [
            "c002_axis_by_library.csv",
            "c002_axis_summary.csv",
            "c002_scores.csv.gz",
            "c002_gene_depth.csv",
            "c002_gene_depth_background.csv",
            "c002_c9_depth_matched.csv",
        ],
    }
    (out / "c002_probe.json").write_text(json.dumps(meta, indent=2))
    print(json.dumps(meta, indent=1))


if __name__ == "__main__":
    main(*sys.argv[1:4])
