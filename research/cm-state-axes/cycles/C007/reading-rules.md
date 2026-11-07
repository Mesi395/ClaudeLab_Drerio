# C007 reading rules (researcher), fixed 2026-11-07 16:06 CET before building sets or computing

## 1. Set-building rule (Li 2025 SD2 = `sources/li2025_supp/MOESM5.xlsx`)

The source cluster for each axis:

| Axis | Sheet | Cluster(s) |
|---|---|---|
| injury_responsive | Stereo-seq regeneration | `1:BZm (activ.)` |
| border_zone_like | Stereo-seq regeneration | genes that are markers of **at least 3 of the 5** BZm domains (`1:` activ., `2:` dediff., `3:` dediff.-prolif., `4:` prolif., `5:` re-diff.), ranked by mean avg_log2FC over the domains where they are markers |
| dedifferentiated | scRNA-seq regeneration | `Cardiomyocytes (dediff.)` |
| proliferating | scRNA-seq regeneration | `Proliferating cells` |

A gene counts as a marker if it passes all of these:

- p_val_adj < 0.05;
- pct.1 ≥ 0.2;
- avg_log2FC ≥ 0.5 for scRNA-seq rows, or ≥ 0.25 for Stereo-seq rows (spot data have smaller fold
  changes);
- present in the sanitised data;
- not a technical gene: `mt-`, `NC-`, `rp[ls]`, `RP[LS]` or `hb` prefixes.

CM-internal filter: drop a gene if, in the **scRNA-seq chamber** sheet,

- it is a marker (p_val_adj < 0.05, avg_log2FC ≥ 0.5) of **any non-CM cluster** (every cluster
  except `Cardiomyocytes (A)` and `Cardiomyocytes (V)`); or
- it is a marker of `Cardiomyocytes (A)` but not of `Cardiomyocytes (V)` (atrial-specific).

Disjoint sets: rank each axis's surviving genes by avg_log2FC (border-zone by mean avg_log2FC). A
gene on several axes goes to the axis where its rank is best. Ties go in this order: proliferating,
dedifferentiated, injury_responsive, border_zone_like.

Gene count: the **top 25** per axis after the disjoint step. If fewer than 25 survive, use all
and report the number.

Output: `markers/li2025_sd2_sets.tsv` (axis, gene, sheet, cluster, rank, avg_log2FC, pct.1,
pct.2, p_val_adj).

## 2. Predicted time (Li used resection; its stages run from 6 hpa to 28 dpa; ours are 3, 7 and 30 dpi)

| Axis | Li stage | Predicted contrast here |
|---|---|---|
| injury_responsive | 6 hpa (activ.) | 3 dpi > Ctrl. Weak expectation: the program may have faded by 3 dpi |
| dedifferentiated | 12 hpa–3 dpa | 3 dpi > Ctrl and 7 dpi > Ctrl |
| proliferating | 7 dpa (peak 7–14 dpi) | 7 dpi > Ctrl |
| border_zone_like | 6 hpa–14 dpa | 7 dpi > Ctrl (3 dpi reported as well) |

## 3. Pipeline (identical for the C004 sets and the Li sets)

- Ventricle-called cells (C004 rule B, now sourced by Li SD2, SUP-E10).
- Binomial thinning to 500 UMIs (seed 0).
- score_genes on log1p CP10k, scored within the ventricle subset.
- 100 expression-matched null sets per axis (25 bins).

## 4. Calibrated tail (declared because the C004 tail rule was miscalibrated)

- Per axis, threshold T = the 99th percentile of the axis score over the ventricle cells of the
  **control hearts**. For each control heart, T is computed from the **other three** control
  hearts only (leave-one-out).
- Tail fraction per heart = share of its ventricle cells above T.
- Contrast: AUC (injured hearts > control hearts) and the difference of heart medians.
- **Null calibration:** the same tail procedure is applied to each of the 100 matched null sets,
  which gives a null distribution of AUCs. This tests whether random genes of the same expression
  show the same time contrast, a contrast that global shifts (depth, contamination, stress)
  between injured and control hearts could produce.
- **Pass:** AUC ≥ 0.8, a positive difference of medians, **and** the real AUC above the 95th
  percentile of the null AUCs.
- **Leave-Hr9-out:** the same with Hr9 removed. A pass that survives without Hr9 is "not driven
  by one heart".
- **Clearer than C004:** more axes pass under this rule, or larger effects (difference of
  medians), on the same pipeline.

Anything changed after results is labelled exploratory.

## Addendum, 2026-11-07 16:08 CET: exploratory variant, after seeing the gene lists, before any scores

The fixed CM-internal filter keeps genes that mark CMs against all cell types, such as sarcomere
genes (dedifferentiated: myh7l, myl7, actc1a; border-zone: actc1, tnnt2a, OXPHOS genes). Exploratory
variant: also drop every marker of Cardiomyocytes (A) or (V) in the chamber sheet (p_adj < 0.05,
log2FC >= 0.5). Same pipeline and pass rule. Labelled exploratory; it is not the C007 verdict.
