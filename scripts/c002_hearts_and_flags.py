"""C002: heart-level summary of per-library tests, and per-gene depth flags
into proposals/loop/*.yaml (moved there from markers/ after the runs).

Hr2a/Hr2b, H8a/H8v, Hr6a/Hr6v and Hr7a/Hr7v are one heart each (supervisor note, C002), so the
22 libraries are 18 hearts. A heart passes a test only if all of its libraries pass (conservative).

uv run python scripts/c002_hearts_and_flags.py results markers
"""

import json
import re
import subprocess
import sys
from datetime import datetime
from pathlib import Path

import pandas as pd
import yaml


def heart(lib: str) -> str:
    return re.sub(r"^(Hr2|H8|Hr6|Hr7)[abv]$", r"\1", lib)


def main(res_dir: str, marker_dir: str) -> None:
    res = Path(res_dir)
    d = pd.read_csv(res / "c002_axis_by_library.csv")
    d["heart"] = d["lib"].map(heart)
    d["methods_agree"] = d["rho_ucell_vs_score_genes"] >= 0.5
    tests = [
        "coherence_above_null",
        "ucell_p90_above_null",
        "score_genes_p90_above_null",
        "methods_agree",
    ]
    h = d.groupby(["state", "heart"])[tests].all().groupby("state").sum()
    h.columns = [f"hearts_{c}" for c in h.columns]
    h["n_hearts"] = d.groupby("state")["heart"].nunique()
    h.to_csv(res / "c002_heart_summary.csv")
    print(h.to_string())

    g = pd.read_csv(res / "c002_gene_depth.csv").set_index("gene")
    for f in sorted(Path(marker_dir).glob("*.yaml")):
        spec = yaml.safe_load(f.read_text())
        for m in spec["genes"]:
            if m["symbol"] in g.index:
                r = g.loc[m["symbol"]]
                m["c002_depth"] = {
                    "median_rho_logtot": round(float(r.median_rho_logtot), 3),
                    "background_p05_p95": [round(float(r.bg_p05), 3), round(float(r.bg_p95), 3)],
                    "tracks_depth": bool(r.tracks_depth),
                    "direction": "" if pd.isna(r.direction) else str(r.direction),
                    "evidence": "results/c002_gene_depth.csv (R-E11)",
                }
            if m["symbol"] == "myh7":
                m["c002_symbol_check"] = (
                    "ZFIN: vmhc is a previous name of myh7 (ZDB-GENE-991123-5); myh7l is a "
                    "different gene, previously vmhcl (ZDB-GENE-061027-393). Kept myh7 as "
                    "sourced; it is detected in few CMs here (R-E9). Not substituted."
                )
        f.write_text(yaml.safe_dump(spec, sort_keys=False, width=100, allow_unicode=True))
    meta = {
        "commit": subprocess.check_output(["git", "rev-parse", "HEAD"], text=True).strip(),
        "dirty_tree": bool(subprocess.check_output(["git", "status", "--porcelain"], text=True)),
        "command": "uv run python " + " ".join(sys.argv),
        "inputs": ["c002_axis_by_library.csv", "c002_gene_depth.csv"],
        "time": datetime.now().astimezone().isoformat(timespec="seconds"),
        "rule": "heart passes only if all its libraries pass; 22 libraries = 18 hearts",
    }
    (res / "c002_heart_summary.json").write_text(json.dumps(meta, indent=2))


if __name__ == "__main__":
    main(*sys.argv[1:3])
