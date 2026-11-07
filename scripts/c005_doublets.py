"""C005 R3.2: scrublet doublet scores on the full GEO per-library matrices, mapped to CM cells.

Rules: research/cm-state-axes/cycles/C005/reading-rules.md section 2 (fixed before computing).
Inputs: data/raw/<lib>_filtered_matrix.h5 (scripts/get_data.sh, GEO GSE159032 / GSE158919),
results/c004_agent_labels.csv (axis scores), results/c005_agent_labels.csv (state calls).

uv run python -I scripts/c005_doublets.py data/raw results
"""

import json
import subprocess
import sys
import time
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import spearmanr

STATES = ["border_zone_like", "dedifferentiated", "injury_responsive", "proliferating"]


def main(raw_dir, out_dir):
    out = Path(out_dir)
    lab4 = pd.read_csv(out / "c004_agent_labels.csv")
    lab5 = pd.read_csv(out / "c005_agent_labels.csv")
    assert (lab4.cell == lab5.cell).all()
    lab = lab4[["cell", "library"] + [f"score_{s}" for s in STATES]].copy()
    lab["agent_state_c005"] = lab5.agent_state
    lab["barcode"] = lab.cell.str.split("_", n=1).str[1] + "-1"
    per_cell, per_lib = [], []
    t0 = time.time()
    for lb in sorted(lab.library.unique()):
        f = Path(raw_dir) / f"{lb}_filtered_matrix.h5"
        row = {"library": lb, "status": "not run"}
        try:
            r = sc.read_10x_h5(f)
            r.var_names_make_unique()
            sc.pp.scrublet(r, random_state=0)
            sub = lab[lab.library == lb]
            mapped = sub.barcode.isin(r.obs_names)
            row.update(
                {
                    "n_all_cells": int(r.n_obs),
                    "n_cm": int(len(sub)),
                    "mapping_rate": float(mapped.mean()),
                    "all_cells_doublet_rate": float(r.obs.predicted_doublet.mean()),
                    "threshold": float(r.uns["scrublet"]["threshold"])
                    if "threshold" in r.uns.get("scrublet", {})
                    else np.nan,
                }
            )
            d = r.obs.loc[sub.barcode[mapped], ["doublet_score", "predicted_doublet"]]
            s = sub[mapped].copy()
            s["doublet_score"] = d.doublet_score.to_numpy()
            s["predicted_doublet"] = d.predicted_doublet.to_numpy().astype(bool)
            per_cell.append(s)
            row["status"] = "ok" if row["mapping_rate"] > 0.9 else "excluded (mapping <= 0.9)"
            row["cm_doublet_rate"] = float(s.predicted_doublet.mean())
            for st in STATES:
                row[f"rho_doublet_{st}"] = float(spearmanr(s.doublet_score, s[f"score_{st}"])[0])
        except Exception as e:  # noqa: BLE001
            row["status"] = f"not run: {type(e).__name__}: {e}"[:200]
        per_lib.append(row)
        print(lb, row["status"], f"{time.time() - t0:.0f}s", flush=True)
        pd.DataFrame(per_lib).to_csv(out / "c005_doublets_by_library.csv", index=False)
    pl = pd.DataFrame(per_lib)
    pc = pd.concat(per_cell)
    ok = pl.status == "ok"
    pc = pc[pc.library.isin(pl.library[ok])]
    pc[["cell", "library", "doublet_score", "predicted_doublet", "agent_state_c005"]].to_csv(
        out / "c005_doublets_cells.csv", index=False
    )
    summ = {
        "libraries_ok": int(ok.sum()),
        "libraries_total": len(pl),
        "not_ok": pl.loc[~ok, ["library", "status"]].to_dict("records"),
        "median_abs_rho_doublet": {
            st: float(pl.loc[ok, f"rho_doublet_{st}"].abs().median()) for st in STATES
        },
        "tracks_doublets_rule_ge_0.1": {},
        "state_doublet_share": {},
    }
    for st in STATES:
        summ["tracks_doublets_rule_ge_0.1"][st] = summ["median_abs_rho_doublet"][st] >= 0.1
    ref = pc[pc.agent_state_c005 == "unassigned"].predicted_doublet.mean()
    summ["unassigned_ventricle_doublet_share"] = float(ref)
    for st in STATES + ["not_scored"]:
        x = pc[pc.agent_state_c005 == st]
        share = float(x.predicted_doublet.mean()) if len(x) else np.nan
        summ["state_doublet_share"][st] = {
            "n": int(len(x)),
            "share": share,
            "ratio_to_unassigned": share / ref if ref else None,
            "enriched_ge_2x": bool(ref and share / ref >= 2),
        }
    summ["meta"] = {
        "commit": subprocess.check_output(["git", "rev-parse", "HEAD"], text=True).strip(),
        "dirty_tree": True,
        "command": "uv run python -I " + " ".join(sys.argv),
        "inputs": [raw_dir, "c004_agent_labels.csv", "c005_agent_labels.csv"],
        "scrublet": "scanpy.pp.scrublet defaults, random_state=0, per library, all cell types",
        "time": datetime.now().astimezone().isoformat(timespec="seconds"),
    }
    (out / "c005_doublets_summary.json").write_text(json.dumps(summ, indent=2, default=str))
    print(json.dumps({k: v for k, v in summ.items() if k != "meta"}, indent=1, default=str))


if __name__ == "__main__":
    main(*sys.argv[1:3])
