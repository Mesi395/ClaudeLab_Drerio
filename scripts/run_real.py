"""Score the cardiomyocyte states on the real data and run every check on each.

    uv run python scripts/run_real.py

Inputs: GEO counts in data/raw/ (scripts/get_data.sh), the cell list config/cells.csv.gz,
marker sets markers/*.yaml (core genes; roles per config/axes.yaml), conditions from
config/dataset.yaml. No labels are read.

Writes results/real/checks.tsv (axis x check outcome and statistic),
results/real/summary.tsv (score mean and 90th percentile per axis and condition)
and results/runs/real-v0/manifest.yaml.
"""

import json
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp

from zf import checks, config, provenance
from zf.panel import Panel
from zf.score import score_mean

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data/raw"
CELLS = ROOT / "config/cells.csv.gz"
OUT = ROOT / "results/real"
N_TOP_GENES = 3000
SEED = 0


def read_library(lib: str, barcodes: list[str]):
    with h5py.File(RAW / f"{lib}_filtered_matrix.h5") as f:
        m = f["matrix"]
        x = sp.csc_matrix((m["data"][:], m["indices"][:], m["indptr"][:]), shape=m["shape"][:])
        names = [n.decode() for n in m["features/name"][:]]
        bc = [b.decode().split("-")[0] for b in m["barcodes"][:]]
    pos = {b: i for i, b in enumerate(bc)}
    cols = [pos[b.split("_", 1)[1]] for b in barcodes]
    return x[:, cols].T.tocsr(), names


def load():
    cells = pd.read_csv(CELLS)
    blocks, names = [], None
    for lib, group in cells.groupby("library", sort=False):
        x, n = read_library(lib, group.cell.tolist())
        if names is None:
            names = n
        assert n == names, "gene order differs between libraries"
        blocks.append(x)
    counts = sp.vstack(blocks).tocsr()
    order = np.concatenate([g.index.to_numpy() for _, g in cells.groupby("library", sort=False)])
    cells = cells.loc[order].reset_index(drop=True)
    first = ~pd.Index(names).duplicated()
    return counts[:, np.flatnonzero(first)], list(pd.Index(names)[first]), cells


def gene_sets(axes):
    sets, by_role = {}, {}
    for axis in axes.axes:
        markers = config.load_markers(ROOT / axis.markers)
        sets[axis.name] = [g.symbol for g in markers.genes if g.role in axis.score_roles]
        for g in markers.genes:
            by_role.setdefault(g.role, []).append(g.symbol)
    controls = {
        c.name: [g.symbol for g in c.genes]
        or sorted({s for r in c.from_roles for s in by_role.get(r, [])})
        for c in axes.controls
    }
    return sets, controls


def main():
    axes = config.load_axes(ROOT / "config/axes.yaml")
    libs = {
        lib.library: lib.condition
        for lib in config.load_dataset(ROOT / "config/dataset.yaml").libraries
    }
    counts, genes, cells = load()
    total = np.asarray(counts.sum(axis=1)).ravel()

    sets, controls = gene_sets(axes)
    tested = dict(sets, library_size_suspect=controls.pop("library_size_suspect"))
    needed = {g for s in [*tested.values(), *controls.values()] for g in s}
    missing = sorted(needed - set(genes))
    detection = np.asarray((counts > 0).mean(axis=0)).ravel()
    keep = set(np.argsort(detection)[::-1][:N_TOP_GENES]) | {
        genes.index(g) for g in needed if g in genes
    }
    keep = sorted(keep)
    sub_genes = [genes[i] for i in keep]

    condition = cells.library.map(libs).to_numpy()
    panel = Panel(
        counts=counts[:, keep],
        genes=sub_genes,
        axis_genes=[],
        library=cells.library.to_numpy(),
        condition=condition,
        covariates={},
        total=total,
    )
    panel.covariates = {
        "depth": np.log(total),
        "ambient": score_mean(
            panel, [g for g in controls["ambient_fibroblast"] if g in panel.index]
        ),
        "dissociation": score_mean(
            panel, [g for g in controls["dissociation"] if g in panel.index]
        ),
    }

    rows, summary, outcomes = [], [], {}
    expect = {a.name: [(e.higher, e.than) for e in a.expect] for a in axes.axes}
    for axis, genes_ in tested.items():
        panel.axis_genes = [g for g in genes_ if g in panel.index]
        panel.expect = expect.get(axis, [])
        result = checks.run_all(panel, seed=SEED)
        outcomes[axis] = result
        for check_id, outcome in result.items():
            stat = {} if outcome == "not_run" and check_id not in ("doublet",) else None
            if outcome != "not_run":
                stat = checks.evaluate(check_id, panel, SEED)
                stat.pop("outcome", None)
            rows.append(
                {"axis": axis, "check": check_id, "outcome": outcome, "statistic": json.dumps(stat)}
            )
        score = score_mean(panel, panel.axis_genes)
        for cond in ["uninjured", "dpi_3", "dpi_7", "dpi_30"]:
            s = score[condition == cond]
            summary.append(
                {
                    "axis": axis,
                    "condition": cond,
                    "n_cells": len(s),
                    "mean": s.mean(),
                    "p90": np.quantile(s, 0.9),
                }
            )
        print(axis, result)

    OUT.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(OUT / "checks.tsv", sep="\t", index=False)
    pd.DataFrame(summary).to_csv(OUT / "summary.tsv", sep="\t", index=False, float_format="%.4f")
    provenance.write_manifest(
        ROOT / "results/runs/real-v0",
        command="uv run python scripts/run_real.py",
        inputs=[
            CELLS,
            ROOT / "config/axes.yaml",
            ROOT / "config/checks.yaml",
            *sorted((ROOT / "markers").glob("*.yaml")),
        ],
        params={
            "seed": SEED,
            "n_cells": int(counts.shape[0]),
            "n_genes_scored": len(sub_genes),
            "genes_missing_from_data": missing,
            "doublet": "not_run: no doublet covariate in this run",
            "expectations_status": "proposed, not yet confirmed by the domain expert",
        },
        outcomes=outcomes,
    )


if __name__ == "__main__":
    main()
