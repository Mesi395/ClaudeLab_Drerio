# ruff: noqa: E501
"""Supervisor checks of C005: corrected label file, and the doublet score distribution.

Usage: uv run python -I research/cm-state-axes/checks/c005_verify.py results/c005_agent_labels.csv results/c005_doublets_by_library.csv
"""

import sys

import pandas as pd

lab = pd.read_csv(sys.argv[1])
print("rows", len(lab), "unique", lab["cell"].nunique())
print(
    "agent_state by agent_chamber:\n",
    pd.crosstab(lab["agent_state"], lab["agent_chamber"]).to_string(),
)
st = lab[~lab["agent_state"].isin(["unassigned", "not_scored"])]
print("\nstate calls per heart (top 6):", st["heart"].value_counts().head(6).to_dict())
d = pd.read_csv(sys.argv[2])
print("\ndoublet table columns:", list(d.columns))
print(d.head(25).to_string())
