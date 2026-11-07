# C005 reading rules (researcher), fixed 2026-11-07 15:48 CET before computing

## 1. Corrected label file (F10)

- Cells: only those called ventricle by the C004 chamber rule B (`agent_chamber == ventricle` in
  `results/c004_agent_labels.csv`). All other cells get agent_state = `not_scored`.
- The state rule is C004 section 4, re-applied **within ventricle-called cells**, as in the C004
  ventricle positive control:
  - cells thinned to 500 UMIs (seed 0), log1p CP10k, score_genes computed on the ventricle subset;
  - per heart and axis, the threshold is the 99th percentile of 100 matched null-set scores over
    that heart's ventricle cells;
  - the primary state is the high axis with the largest (score − threshold) / null SD;
    `unassigned` if no axis is high.
- Report the label counts from c004 and c005 side by side. The c004 file is not changed.
- Chance expectation per axis: about 1% of scored cells. A label class is "at chance" if its count
  is within 0.5–2× of that expectation.

## 2. Doublets (R3.2), only if the GEO download completes

- Per library: scanpy `pp.scrublet` on the full GEO filtered matrix (all cell types), with default
  parameters and random_state 0. CM barcodes are mapped by barcode sequence within the library;
  the mapping rate is reported and must exceed 90% per library, or that library is excluded and
  named.
- "Likely doublet" = scrublet's predicted_doublet call (its automatic threshold). The score is also
  used continuously.
- Association: within-library Spearman between the doublet score and each axis score (the C004
  all-cell thinned score_genes, `score_<axis>` in c004). A score "tracks doublets" if the median
  |rho| over libraries is ≥ 0.1.
- State calls among doublets: per state in the c005 file, the share of its cells that are likely
  doublets, against the share among `unassigned` ventricle cells. "Enriched" if ≥ 2×.
- If scrublet fails in a library, it is reported as not run, never as passed.

## Addendum, 2026-11-07 15:52 CET, exploratory, after partial results

Scrublet's automatic threshold called 0 doublets in the first libraries (H5, H6), so the
pre-declared share among likely doublets may be empty. Exploratory fallback: the top 10% of CM cells
by doublet score within each library are high-score (scripts/c005_doublets_exploratory.py).
It is reported as exploratory and is not part of the R3.2 verdict.
