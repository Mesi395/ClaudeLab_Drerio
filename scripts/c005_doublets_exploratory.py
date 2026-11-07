"""C005 exploratory (not in the fixed rules; added after scrublet's automatic threshold called no
or few doublets): per library, the top 10% of CM cells by doublet score are "high-score"; the share
of each c005 state among high-score cells vs among unassigned ventricle cells.

uv run python scripts/c005_doublets_exploratory.py results
"""

import json
import subprocess
import sys
from datetime import datetime
from pathlib import Path

import pandas as pd

STATES = ["border_zone_like", "dedifferentiated", "injury_responsive", "proliferating"]


def main(out_dir):
    out = Path(out_dir)
    c = pd.read_csv(out / "c005_doublets_cells.csv")
    c["high_score"] = c.groupby("library").doublet_score.transform(lambda s: s >= s.quantile(0.9))
    ref = c[c.agent_state_c005 == "unassigned"].high_score.mean()
    rows = []
    for st in STATES + ["unassigned", "not_scored"]:
        x = c[c.agent_state_c005 == st]
        rows.append(
            {
                "agent_state": st,
                "n": len(x),
                "share_high_score": x.high_score.mean(),
                "ratio_to_unassigned": x.high_score.mean() / ref,
            }
        )
    d = pd.DataFrame(rows)
    d.to_csv(out / "c005_doublets_exploratory_top10.csv", index=False)
    meta = {
        "commit": subprocess.check_output(["git", "rev-parse", "HEAD"], text=True).strip(),
        "dirty_tree": True,
        "command": "uv run python " + " ".join(sys.argv),
        "input": "c005_doublets_cells.csv",
        "label": "exploratory: top-10% within-library CM doublet score, added after results",
        "time": datetime.now().astimezone().isoformat(timespec="seconds"),
    }
    (out / "c005_doublets_exploratory_top10.json").write_text(json.dumps(meta, indent=2))
    print(d.round(3).to_string(index=False))


if __name__ == "__main__":
    main(sys.argv[1])
