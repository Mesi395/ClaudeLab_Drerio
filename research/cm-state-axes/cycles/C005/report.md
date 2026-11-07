# C005 report: fix F10, then Probe R3.2 (researcher)

Written 2026-11-07 15:55 CET. Rules fixed at 15:48, before computing: [reading-rules.md](reading-rules.md).
Evidence R-E17 and R-E18.

## 1. F10 fixed: corrected label file

- New file: [results/c005_agent_labels.csv](../../../../results/c005_agent_labels.csv), 32,044
  rows. Format and provenance are in [c005_agent_labels.json](../../../../results/c005_agent_labels.json).
  `c004_agent_labels.csv` is kept unchanged.
- Rule: the C004 state rule, re-applied within the 14,024 ventricle-called cells (scoring and null
  thresholds computed on that subset). The 18,020 atrial or chamber-unassigned cells are
  `not_scored`.
- Effect ([c005_label_counts_vs_c004.csv](../../../../results/c005_label_counts_vs_c004.csv)):

| agent_state | c004 | c005 | c005 / chance (about 140) |
|---|---|---|---|
| border_zone_like | 186 | 67 | 0.48 (below chance) |
| dedifferentiated | 199 | 53 | 0.38 (below chance) |
| injury_responsive | 1,188 | 145 | 1.03 (at chance) |
| proliferating | 407 | 183 | 1.30 (at chance) |
| unassigned | 30,064 | 13,576 | |
| not_scored | 0 | 18,020 | |

- On ventricle cells, c005 and c004 agree for 0.994 of cells (mostly `unassigned`). The F10
  excess of injury-responsive calls was atrial: from 1,188 calls to 145.
- The calls concentrate in Hr9 (35, 20, 42, 52) and Hr11 (20, 21, 20, 8); control hearts have
  almost none ([c005_labels_by_heart.csv](../../../../results/c005_labels_by_heart.csv)).
- Note for the human who runs evaluate.py: [note-for-evaluate.md](note-for-evaluate.md).

**Claim C5-1.** After chamber correction, no state label class rises above the chance
expectation; two are below it. The label file honestly encodes "no separable states found".
`agent_chamber` is its defensible content.

## 2. Doublets (R3.2): run, with a threshold failure

- **Run:** scanpy `pp.scrublet` (defaults, random_state 0) on all 22 GEO filtered matrices, all cell
  types, 3,745–16,587 cells per library. CM barcodes mapped at rate 1.0 in every library
  ([c005_doublets_by_library.csv](../../../../results/c005_doublets_by_library.csv), summary
  [c005_doublets_summary.json](../../../../results/c005_doublets_summary.json), per-cell
  [c005_doublets_cells.csv](../../../../results/c005_doublets_cells.csv)).
- The first attempt failed in every library because scikit-image was missing. It is kept as
  `results/c005_doublets_attempt1_failed.*`. Fixed with `uv add scikit-image`.

**C5-2, threshold failure (observed).** Scrublet's automatic threshold (0.45–0.67) called a doublet
rate of 0.000 among all cells in all 22 libraries, to three decimals. Libraries of this size are
expected to hold a few percent doublets or more (an expectation, not measured here). So the
automatic call failed to separate doublets, and the pre-declared "share of states among likely
doublets" is **uninformative**: zero doublets in every class. That part of R3.2 is **not
answered**, which is not the same as "no doublets".

**C5-3, continuous score (pre-declared rule).** Median within-library |rho| between the doublet
score and the axis scores:

| Axis | Median \|rho\| | Tracks doublets (≥ 0.1)? |
|---|---|---|
| border-zone | 0.067 | no |
| dedifferentiated | 0.079 | no |
| injury-responsive | 0.091 | no |
| proliferating | 0.050 | no |

By the rule, no axis tracks doublets. Two patterns qualify that reading:

- **Atrium-only libraries are positive:** H8a 0.32–0.34, Hr7a 0.13–0.14, Hr6a 0.12–0.16 for
  border-zone, dedifferentiated and injury-responsive. In atrial libraries, high nppa-axis cells
  have higher doublet scores. The *nppa*-led axes are mixed with whatever drives scrublet there.
- **Hr11 and Hr27 are negative** (−0.20 and −0.20 for border-zone). This is the opposite of a
  doublet explanation for Hr11's high scores. It also fits C003: contamination dilutes the CM
  genes rather than adding *nppa*.

**C5-4, exploratory (top 10% of CM doublet scores per library; added after results).** The share of
each state's cells that are high-score, relative to unassigned ventricle cells: border-zone 0.32,
dedifferentiated 0.81, injury-responsive 0.74, proliferating 0.77, and not_scored (atrial or
unassigned) 1.28 ([c005_doublets_exploratory_top10.csv](../../../../results/c005_doublets_exploratory_top10.csv)).
No state is enriched among high doublet scores; all four are slightly depleted.

**Reading.** Doublets do not explain the agent's state calls or the axis scores within
libraries. The check is weaker than intended for two reasons:

- the automatic doublet call failed;
- scrublet detects mainly heterotypic doublets, and its scores were not validated here (no
  expected-rate threshold, no second method such as scDblFinder).

## Requirement coverage

- R5.2: corrected file delivered, with a note for evaluate.py.
- R3.2: run in 22/22 libraries.
  - Continuous association: answered (no axis tracks doublets by the rule).
  - Share of states among called doublets: not answered (automatic threshold failure).
  - Exploratory top-10%: no enrichment.

## Weaknesses

- The scrublet threshold was not recalibrated. The exploratory threshold is relative within CMs,
  not against all cells.
- One doublet method.
- Doublet scores are compared with the c004 all-cell axis scores, as pre-declared, not with the
  ventricle-only c005 scores.

## Status

**Review requested.** Item 1 is done. Item 2 is run, with its limits stated; the doublet-share
question is not answered.
