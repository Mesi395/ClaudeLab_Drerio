# ruff: noqa: E501
"""Supervisor checks of C001 claims C1, C8, C9 and the myh7 anomaly.

Usage: uv run python -I research/cm-state-axes/checks/c001_verify.py data/clean/hu2022_cm_allgenes_sanitized.h5ad
Reads only the sanitised file. Prints summaries, never label values (there are none to print).
"""

import sys

import anndata as ad
import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp

path = sys.argv[1]

# C1: structure of the sanitised file
with h5py.File(path, "r") as f:
    print("top-level keys:", list(f.keys()))
    print("obs keys:", list(f["obs"].keys()))
    print("var keys:", list(f["var"].keys()))
    for k in ("uns", "obsm", "obsp", "varm", "layers"):
        if k in f:
            print(k, "keys:", list(f[k].keys()))
    hits = []
    f.visit(
        lambda n: (
            hits.append(n)
            if any(
                s in n.lower()
                for s in (
                    "ident",
                    "label",
                    "annot",
                    "cluster",
                    "subtype",
                    "celltype",
                    "cell_type",
                )  # token list generalised for the public repo
            )
            else None
        )
    )
    print("label-like HDF5 paths:", hits)

a = ad.read_h5ad(path)
X = a.X.tocsr() if sp.issparse(a.X) else sp.csr_matrix(a.X)
lib = a.obs["orig.ident"].astype(str).to_numpy()
tot = np.asarray(X.sum(1)).ravel()
print("integer counts:", bool(np.all(X.data == np.round(X.data))), "| shape", a.shape)

# C8: per-library median depth, recomputed
med = pd.Series(tot).groupby(lib).median().sort_values()
print("\nC8 median total counts per CM, lowest 6:\n", med.head(6).to_string())
print("highest 3:\n", med.tail(3).to_string())

# C9: detection of nppa and ttn.2 in Hr11 vs H7, plus share of UMIs
genes = list(a.var_names)


def det(g, lb):
    j = genes.index(g)
    col = X[:, j].toarray().ravel()
    m = lib == lb
    return (col[m] > 0).mean(), np.median(col[m] / tot[m])


for g in ("nppa", "ttn.2"):
    for lb in ("Hr11", "H7"):
        d, s = det(g, lb)
        print(f"C9 {g} {lb}: detected in {d:.4f} of cells; median share of cell UMIs {s:.5f}")

# myh7 anomaly: detection of myosin heavy chain symbols and top genes by total counts
gs = np.asarray(X.sum(0)).ravel()
detfrac = np.asarray((X > 0).sum(0)).ravel() / X.shape[0]
myh = [i for i, g in enumerate(genes) if g.lower().startswith(("myh", "vmhc", "amhc"))]
print("\nmyosin heavy chain symbols (gene, detected fraction, share of all UMIs):")
for i in sorted(myh, key=lambda i: -gs[i]):
    print(f"  {genes[i]:12s} {detfrac[i]:.4f} {gs[i] / gs.sum():.5f}")
top = np.argsort(-gs)[:15]
print("\ntop 15 genes by total counts:", [(genes[i], round(gs[i] / gs.sum(), 4)) for i in top])
