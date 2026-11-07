# Note for the human who runs `src/zf/evaluate.py`

Written by the researcher, 2026-11-07 15:55 CET. The agent never read `holdout/`.

**Use `results/c005_agent_labels.csv`.** It supersedes `c004_agent_labels.csv`, which is kept
unchanged. The c004 file scored all cells regardless of chamber; most of its injury-responsive
calls were atrial cells (supervisor check, F10).

One row per cell (32,044 rows); `cell` = obs_names of the sanitised file. The format is in
`results/c005_agent_labels.json` (`label_file_format`).

| Column | What it claims | How far to trust it |
|---|---|---|
| `agent_chamber` | atrium, ventricle or unassigned, from a rule on the sourced chamber genes *myh6* (atrial) and *myh7* plus `myh7l` (ventricular). | **The most defensible column.** On the six chamber-only libraries it has balanced accuracy 0.991 at coverage 0.735 (`results/c004_chamber_validation.csv`). The figure is optimistic: the rule was chosen on those same libraries, and `myh7l`'s ventricular role is not sourced. |
| `agent_state` | For ventricle-called cells: the strongest axis above its per-heart null 99th percentile, or `unassigned`. Atrial and chamber-unassigned cells are `not_scored`. | **At or below chance.** Chance alone gives about 1% of the 14,024 scored cells per axis (about 140). The counts are border-zone 67 (0.48×), dedifferentiated 53 (0.38×), injury-responsive 145 (1.03×), proliferating 183 (1.30×) (`results/c005_label_counts_vs_c004.csv`). |
| `agent_states_all` | Every axis above its threshold. | Same caveat. |
| `score_<axis>` | Depth-thinned score_genes within ventricle-called cells. | These are continuous scores, not states. They track depth (C002), and no axis passed a convincing positive control (C003, C004). |

How to read an agreement result:

- Agreement on `agent_chamber` tests the chamber rule.
- Agreement on `agent_state` mostly tests whether the expert's states coincide with the few cells
  above a near-chance threshold. Those cells concentrate in hearts Hr9 and Hr11
  (`results/c005_labels_by_heart.csv`).
- Low agreement on states is the expected result, and it is consistent with the study's finding:
  the sourced marker sets do not mark separable states in these data.
