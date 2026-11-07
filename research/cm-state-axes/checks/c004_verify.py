# ruff: noqa: E501
"""Supervisor checks of C004: label file integrity, chamber accuracy from the file, injury calls by chamber.

Usage: uv run python -I research/cm-state-axes/checks/c004_verify.py data/clean/hu2022_cm_allgenes_sanitized.h5ad results/c004_agent_labels.csv
"""

import sys

import anndata as ad
import pandas as pd

a = ad.read_h5ad(sys.argv[1], backed="r")
obs = a.obs[["orig.ident", "AV", "time"]].astype(str)
lab = pd.read_csv(sys.argv[2])
print("columns:", list(lab.columns))
print(
    "rows:",
    len(lab),
    "| unique cells:",
    lab["cell"].nunique(),
    "| all cells in data:",
    set(lab["cell"]) == set(obs.index),
)
lab = lab.set_index("cell").join(obs)
print("agent_state counts:", lab["agent_state"].value_counts().to_dict())

six = lab[lab["AV"] != "Wholeheart"]
ct = pd.crosstab([six["orig.ident"], six["AV"]], six["agent_chamber"])
print("\nchamber call vs metadata in the six chamber-only libraries:\n", ct.to_string())
called = six[six["agent_chamber"].isin(["atrium", "ventricle"])]
acc = (
    (called["agent_chamber"].str[0] == called["AV"].str.lower().str[0]).groupby(called["AV"]).mean()
)
print(
    "accuracy among called cells by true chamber:",
    acc.round(3).to_dict(),
    "| balanced:",
    round(acc.mean(), 3),
    "| coverage:",
    round(len(called) / len(six), 3),
)

inj = lab[lab["agent_state"] == "injury_responsive"]
print("\ninjury_responsive calls by agent_chamber:", inj["agent_chamber"].value_counts().to_dict())
print("injury_responsive calls by time:", inj["time"].value_counts().to_dict())
print("all cells by time:", lab["time"].value_counts().to_dict())
for s in ["injury_responsive", "proliferating", "dedifferentiated", "border_zone_like"]:
    sub = lab[lab["agent_state"] == s]
    print(
        f"{s}: share by time {(sub['time'].value_counts() / lab['time'].value_counts()).round(4).to_dict()}"
    )
