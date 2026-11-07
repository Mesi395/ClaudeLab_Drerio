"""Set-level pseudobulk test on the saved per-gene DESeq2 tables (no refit).

Octo 2026-11-07 16:43 (two reviews: Fable + paper-memory reviewer). Per marker set and model:
mean DESeq2 Wald stat vs 500 decoy sets matched gene-for-gene on baseMean bin -> empirical p.
Mitochondrial genes (mt-*, NC-002333*) are dropped first. Plus the hypergeometric overlap of the
3 dpi up set (raw, thinned600) with every set, on the 18,794-gene background and on the tested genes.

Usage: uv run python scripts/pseudobulk_set_test.py
"""

import re

import numpy as np
import pandas as pd
import yaml
from scipy.stats import hypergeom

from run_protocol import DISSOC_EXTRA, ROOT

OUT = ROOT / "results" / "hk"
MK = OUT / "markers_v1"
N_DECOY, N_BINS, SEED = 500, 20, 20261107
BACKGROUND = 18794  # genes in the count matrix (directive); tested genes reported alongside
MITO = re.compile(r"^(mt-|NC-002333)")
PADJ, LFC = 0.05, 1.0  # same DEG rule as pseudobulk_de.py


def yaml_set(name, roles=None):
    g = yaml.safe_load((MK / f"{name}.yaml").read_text())["genes"]
    return [x["symbol"] for x in g if roles is None or x["role"] in roles]


def marker_sets():
    core = ["nppb", "mustn1b", "myh7", "junba", "junbb", "fosab", "fosl1a"]  # nppa out (Emi)
    return {
        "memory_bz_core": core,
        "memory_bz_core_no_ap1": [g for g in core if not g.startswith(("jun", "fos"))],
        "border_zone": yaml_set("border_zone"),
        "border_zone_core": yaml_set("border_zone", {"core"}),
        "dedifferentiated": yaml_set("dedifferentiated"),
        "injury_responsive": yaml_set("injury_responsive"),
        "proliferating": yaml_set("proliferating"),
        "glycolysis_honkoop2019": ["hk1", "pkma", "pdk2a"],
        "dissociation_control": yaml_set("controls"),
    }


def load(path):
    r = pd.read_csv(path, index_col=0)
    return r[~r.index.str.match(MITO)].dropna(subset=["stat"])


def decoy_test(res, genes, rng):
    """Mean Wald stat of the set vs N_DECOY sets drawn bin-for-bin on log baseMean."""
    s = [g for g in genes if g in res.index]
    if not s:
        return {"n_in_table": 0}
    bins = pd.qcut(np.log1p(res["baseMean"]), N_BINS, labels=False, duplicates="drop")
    pool = {b: np.setdiff1d(np.flatnonzero(bins.to_numpy() == b), res.index.get_indexer(s))
            for b in bins.unique()}
    stat = res["stat"].to_numpy()
    obs = stat[res.index.get_indexer(s)].mean()
    want = bins.loc[s].value_counts()
    null = np.array([
        np.concatenate([stat[rng.choice(pool[b], k, replace=False)] for b, k in want.items()]).mean()
        for _ in range(N_DECOY)
    ])
    return {
        "n_in_table": len(s),
        "mean_stat": obs,
        "decoy_mean": null.mean(),
        "decoy_sd": null.std(ddof=1),
        "z_vs_decoy": (obs - null.mean()) / null.std(ddof=1),
        "p_up": (1 + (null >= obs).sum()) / (N_DECOY + 1),
        "p_two_sided": (1 + (np.abs(null - null.mean()) >= abs(obs - null.mean())).sum())
        / (N_DECOY + 1),
        "genes_in_table": ";".join(s),
    }


def up_set(res):
    ambient_path = ROOT / "data" / "derived" / "ambient_genes.tsv"
    ambient = set(pd.read_csv(ambient_path)["gene"]) if ambient_path.exists() else set()
    dissoc = set(yaml_set("controls")) | set(res.index[res.index.str.match(DISSOC_EXTRA)])
    r = res.dropna(subset=["padj"])
    r = r[(r["padj"] < PADJ) & (r["log2FoldChange"] > LFC) & ~r.index.isin(dissoc | ambient)]
    return set(r.index), ambient_path.exists()


def main():
    rng = np.random.default_rng(SEED)
    sets = marker_sets()
    rows = []
    for run in ["raw", "thinned600"]:
        for model, f in [
            ("ventricle_dpi3_vs_ctrl", f"pseudobulk_de_{run}_ventricle_dpi_3_vs_ctrl.csv.gz"),
            ("interaction_dpi3", f"pseudobulk_de_{run}_interaction_dpi_3.csv.gz"),
        ]:
            res = load(OUT / f)
            for name, genes in sets.items():
                rows.append({"test": "decoy_mean_wald", "run": run, "model": model, "set": name,
                             "set_size": len(genes), "tested_genes": len(res),
                             **decoy_test(res, genes, rng)})
            if model != "ventricle_dpi3_vs_ctrl":
                continue
            up, amb = up_set(res)
            for name, genes in sets.items():
                k = up & set(genes)
                for bg, n_bg, n_set in [("matrix_18794", BACKGROUND, len(genes)),
                                        ("tested", len(res), len(set(genes) & set(res.index)))]:
                    rows.append({"test": f"hypergeom_up_{bg}", "run": run, "model": model,
                                 "set": name, "set_size": n_set, "tested_genes": len(res),
                                 "up_set_size": len(up), "overlap": len(k),
                                 "overlap_genes": ";".join(sorted(k)),
                                 "expected": len(up) * n_set / n_bg,
                                 "p_up": hypergeom.sf(len(k) - 1, n_bg, n_set, len(up)),
                                 "ambient_list_used": amb})
    out = pd.DataFrame(rows)
    out.to_csv(OUT / "pseudobulk_set_test.csv", index=False)
    cols = ["test", "run", "model", "set", "n_in_table", "mean_stat", "z_vs_decoy", "p_up",
            "p_two_sided", "overlap", "up_set_size"]
    print(out[[c for c in cols if c in out]].to_string())


if __name__ == "__main__":
    main()
