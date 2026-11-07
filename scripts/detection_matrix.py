"""Evaluate the checks on planted defects: how often each check flags each defect.

    uv run python scripts/detection_matrix.py [n_seeds]

Writes results/detection/matrix.tsv (rate of `flag` per defect and check, signal checks first)
and results/runs/detection-matrix/manifest.yaml.
"""

import sys
from pathlib import Path

import pandas as pd

from zf import checks, config, provenance, synthetic

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results/detection/matrix.tsv"


def main(n_seeds: int = 20):
    gating = [
        c for c in config.load_checks(ROOT / "config/checks.yaml").checks if c.family != "report"
    ]
    defects = [None] + [c.detects for c in gating]
    rows = []
    for defect in defects:
        for seed in range(n_seeds):
            outcomes = checks.run_all(synthetic.make(defect, seed=seed), seed=seed)
            rows += [
                {"defect": defect or "clean", "seed": seed, "check": k, "outcome": v}
                for k, v in outcomes.items()
            ]
    long = pd.DataFrame(rows)
    rate = long.assign(flag=long.outcome.eq("flag")).pivot_table(
        index="defect", columns="check", values="flag", aggfunc="mean", sort=False
    )[[c.id for c in gating]]
    OUT.parent.mkdir(parents=True, exist_ok=True)
    rate.to_csv(OUT, sep="\t", float_format="%.2f")

    own = {c.detects: c.id for c in gating}
    detection = {d: float(rate.loc[d, own[d]]) for d in own}
    off = rate.drop(index="clean").copy()
    for d, c in own.items():
        off.loc[d, c] = float("nan")
    outcomes = {
        "summary": {
            "min_detection": "pass" if min(detection.values()) >= 0.9 else "flag",
            "false_alarms_clean": "pass" if rate.loc["clean"].max() == 0 else "flag",
            "max_off_diagonal": "pass" if float(off.max().max()) <= 0.1 else "flag",
        }
    }
    provenance.write_manifest(
        ROOT / "results/runs/detection-matrix",
        command=f"uv run python scripts/detection_matrix.py {n_seeds}",
        inputs=[ROOT / "config/checks.yaml", OUT],
        params={
            "n_seeds": n_seeds,
            "detection": detection,
            "max_off_diagonal": float(off.max().max()),
            "false_alarm_clean_max": float(rate.loc["clean"].max()),
        },
        outcomes=outcomes,
    )
    print(rate.to_string(float_format="%.2f"))


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 20)
