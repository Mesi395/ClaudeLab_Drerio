# C004 reading rules (researcher), fixed 2026-11-07 15:43 CET before computing

## 1. Chamber call (per cell, raw counts)

- Sourced genes:
  - atrial: *myh6* (amhc; Berdougo et al. 2003, doi:10.1242/dev.00838, via ZFIN ZDB-PUB-031103-20,
    a search-level check)
  - ventricular: *myh7* (vmhc; Yelon et al. 1999, doi:10.1006/dbio.1999.9406; ZFIN rename R-E7)
- **Rule A** (sourced genes only): atrium if count(myh6) > count(myh7); ventricle if
  count(myh7) > count(myh6); unassigned if they are equal (including both zero).
- **Rule B** adds `myh7l` (formerly *vmhcl*, "ventricular myosin heavy chain-like"). The name
  suggests ventricle, but no source read states it, so this rule is flagged. The ventricular count
  is myh7 + myh7l.
- Validation is on the six chamber-only libraries, which carry metadata AV labels, not expert
  labels. Per library: coverage (the share of cells assigned) and accuracy among assigned cells.
  Across the six: balanced accuracy over the atrium and ventricle libraries.
- The rule with the higher balanced accuracy is used downstream. Because the rule is selected on
  the same six libraries, its accuracy is optimistic; this is stated in the report.
- A rule is usable if its balanced accuracy is ≥ 0.9 and its coverage is ≥ 0.5.

## 2. Positive control in ventricle-called cells

- Same method as C003:
  - cells thinned to 500 UMIs, seed 0
  - score_genes, with 100 matched null sets (25 bins)
  - per heart, the tail fraction above the null's 99th percentile
  - pass = AUC ≥ 0.8 and a positive difference of heart medians, with the same predicted contrasts
- It is computed **only on ventricle-called cells**, and scoring is done on that subset.
- Hearts with fewer than 50 ventricle-called cells are reported but left out of the contrasts.

## 3. G4 (only if an axis passes on the tail)

- Disjoint sets: any gene that is on more than one axis is removed from all axes.
- Per heart, the Spearman correlation between each pair of axis scores is compared with the
  distribution of Spearman correlations between pairs of matched null sets (100 pairs).
- "Co-occur" if the real rho exceeds the null 95th percentile in ≥ 80% of hearts. "Come apart" if
  it lies below the null 5th percentile in ≥ 80% of hearts. Otherwise "mixed".
- If no axis passes, G4 is declared "cannot be answered with these sets", with the evidence: the
  shared genes and the chamber and depth ties.

## 4. Agent label file (R5.2)

- All 32,044 cells.
- Scores: thinned score_genes on all cells, with the sourced axis sets minus ttn.2, as in C003.
  For each heart and axis, the threshold is the 99th percentile of the null-set scores over that
  heart's cells.
- An axis is "high" in a cell if its score is above the threshold.
- agent_state is "unassigned" if no axis is high. Otherwise it is the high axis with the largest
  (score − threshold) / SD of the heart's null scores. agent_states_all lists every high axis.
- agent_chamber comes from the chamber rule chosen in step 1.
- This rule is fixed before any look at any annotation; holdout/ is never read.
- The labels carry the caveat that no axis passed a positive control (C003), unless C004 changes
  that.
