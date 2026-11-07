"""C003 exploratory (not in the fixed rules): within each 3 dpi / Ctrl library, Spearman of a cell's
nppa, ttn.2 and hsp70l share with its non-CM panel share, and with its mito share (incl. NC).

uv run python scripts/c003_exploratory_cellcorr.py \
    data/clean/hu2022_cm_allgenes_sanitized.h5ad results
"""

import json
import re
import subprocess
import sys
from datetime import datetime
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

NON_CM = [
    "hbba1.1",
    "hbba1",
    "col1a1a",
    "col1a2",
    "postnb",
    "fn1b",
    "kdrl",
    "cdh5",
    "fli1a",
    "mpeg1.1",
    "lyz",
    "tcf21",
]


def main(src, out_dir):
    a = ad.read_h5ad(src)
    X = a.layers["counts"].tocsr()
    tot = np.asarray(X.sum(1)).ravel()
    vn = list(a.var_names)
    vidx = {g: i for i, g in enumerate(vn)}
    noncm = np.asarray(X[:, [vidx[g] for g in NON_CM]].sum(1)).ravel() / tot
    mito = (
        np.asarray(
            X[:, [i for i, g in enumerate(vn) if re.match(r"^(mt-|NC-002333)", g)]].sum(1)
        ).ravel()
        / tot
    )
    lib = a.obs["orig.ident"].astype(str).to_numpy()
    rows = []
    for lb in ["Hr11", "Hr12", "Hr24", "Hr27", "H5", "H6", "H7"]:
        m = lib == lb
        for g in ["nppa", "ttn.2", "hsp70l"]:
            s = X[m][:, vidx[g]].toarray().ravel() / tot[m]
            rows.append(
                {
                    "lib": lb,
                    "gene": g,
                    "n_cells": int(m.sum()),
                    "rho_with_noncm_share": float(spearmanr(s, noncm[m])[0]),
                    "rho_with_mito_share": float(spearmanr(s, mito[m])[0]),
                }
            )
    d = pd.DataFrame(rows)
    d.to_csv(Path(out_dir) / "c003_exploratory_cellcorr.csv", index=False)
    meta = {
        "commit": subprocess.check_output(["git", "rev-parse", "HEAD"], text=True).strip(),
        "dirty_tree": True,
        "command": "uv run python " + " ".join(sys.argv),
        "input": src,
        "label": "exploratory, added after seeing C003 results",
        "time": datetime.now().astimezone().isoformat(timespec="seconds"),
    }
    (Path(out_dir) / "c003_exploratory_cellcorr.json").write_text(json.dumps(meta, indent=2))
    print(d.round(3).to_string(index=False))


if __name__ == "__main__":
    main(*sys.argv[1:3])
