# C007 report: Probe of direction A (researcher)

Written 2026-11-07 16:12 CET. Rules fixed at 16:06, before building sets or computing:
[reading-rules.md](reading-rules.md). Evidence R-E20.

## The sets

[markers/li2025_sd2_sets.tsv](../../../../markers/li2025_sd2_sets.tsv) gives the sheet, cluster,
rank, log2FC, pct and adjusted p for every gene. There are 25 genes per axis, disjoint across axes,
after dropping 3,304 genes that mark non-CM or atrial-only clusters.

- **proliferating** (scRNA-seq "Proliferating cells"): a clean cell-cycle set (*tk1*, *cenpf*,
  *aspm*, *tpx2*, *birc5a*, *plk1*, *cdk1*, *ccnb1*, *nusap1*, *ndc80*…).
- **injury_responsive** (Stereo-seq BZm activ.): *tnfrsf11b*, *nppb*, *hspb11*, *flnca*, *rln3a*,
  *hspb9*, *hspb6*… It overlaps the C001 injury set (*tnfrsf11b*, *hspb6*, *hspb11*, *nppb*).
- **dedifferentiated** (scRNA-seq "Cardiomyocytes (dediff.)"): mostly **CM-identity and sarcomere
  genes** (*myh7*, *csrp3*, *desma*, *mustn1b*, *actn2b*, *myh7l*, *nppa*, *myl7*, *actc1a*…).
- **border_zone_like** (markers of at least 3 of 5 BZm domains): mostly **sarcomere and OXPHOS
  genes** (*ACTC1*, *acta1a*, *xirp1*, *cox6a2*, *tnnt2a*, *atp5f1b*…).

Observation: Li's markers are each cluster against **all** cell types. For the dedifferentiated
and border-zone clusters, they are dominated by genes that separate CMs from non-CMs, not one CM
state from another. In the exploratory variant that also removes every CM-identity marker,
dedifferentiated keeps 3 genes and border-zone 1
([sets](../../../../markers/li2025_sd2_sets_exploratory_noCMid.tsv)). **Li's SD2 does not contain
CM-internal marker lists for these two states.**

## Side-by-side positive control ([c007_positive_control.csv](../../../../results/c007_positive_control.csv))

Same pipeline for both families: ventricle-called cells, cells thinned to 500 UMIs, score_genes,
18 hearts. The tail threshold is calibrated on control hearts (leave-one-out for controls). A pass
also requires the AUC to beat the 95th percentile of AUCs from 100 matched random sets.

| Axis | Contrast | C004 AUC (diff) | Li SD2 AUC (diff) | Null AUC p95 (C004 / Li) | Without Hr9 (C004 / Li) | Pass? |
|---|---|---|---|---|---|---|
| border-zone | 7 dpi > Ctrl | 0.714 (+0.016) | 0.536 (+0.004) | 0.789 / 0.821 | 0.667 / 0.500 | No / No |
| border-zone | 3 dpi > Ctrl | 0.750 (+0.016) | 0.312 (−0.012) | 0.750 / 0.628 | 0.750 / 0.312 | No / No |
| dedifferentiated | 3 dpi > Ctrl | 0.438 (−0.002) | 0.188 (−0.006) | 0.750 / 0.625 | same | No / No |
| dedifferentiated | 7 dpi > Ctrl | 0.536 (+0.009) | 0.464 (0.000) | 0.789 / 0.821 | 0.458 / 0.417 | No / No |
| injury-responsive | 3 dpi > Ctrl | 0.750 (+0.010) | 0.438 (−0.005) | 0.812 / 0.720 | same | No / No |
| proliferating | 7 dpi > Ctrl | 0.571 (+0.006) | 0.643 (+0.002) | 0.821 / 0.804 | 0.542 / 0.583 | No / No |

**Claim C7-1. No axis passes, in either family.** The Li SD2 sets do not beat the C004 sets.
They are lower on 5 of 6 contrasts; proliferating is the exception (0.643 against 0.571), and
it is still well below its null.

**Claim C7-2. The two C004 "passes" do not survive calibration.** C004 had injury-responsive at
3 dpi (0.812) and dedifferentiated at 7 dpi (0.857) passing on an uncalibrated tail.

- With a control-calibrated threshold, they become 0.750 and 0.536.
- Random gene sets of matched expression reach AUCs of 0.63–0.82 at their 95th percentile. So
  with 4 control hearts, an AUC near 0.8 is within what random genes give. The C004 rule
  (AUC ≥ 0.8 alone) was too lenient.

**Claim C7-3. Per-heart values show heart-level outliers, not a time effect**
([c007_heart_tails.csv](../../../../results/c007_heart_tails.csv)):

- Hr9 and Hr8 (both 7 dpi) are high on nearly every axis in both families. Hr9 has C004
  border-zone 0.241 and Li proliferating 0.041; Hr8 has Li border-zone 0.170. Hr11 is high only on
  the C004 border-zone set (0.124).
- The other injured hearts sit at control level.
- The clean Li cell-cycle set finds a tail in Hr9 (0.041) and Hr8 (0.016), while the other 7 dpi
  hearts sit at 0.000–0.007 and the controls at 0.002–0.007. Without Hr9 its AUC falls to 0.583.
- For comparison, the 2026 resection study reports 1–2% cycling CMs at 7 dpi (abstract, R-E19).
  Most of our 7 dpi hearts do not show that.

## Exploratory, after the gene lists (no-CM-identity sets)

[results/c007_exploratory_noCMid_positive_control.csv](../../../../results/c007_exploratory_noCMid_positive_control.csv).
It finished (exit 0). **No pass either.** Results for the Li sets:

- border-zone (1 gene, *tmem199*): 7 dpi AUC 0.821, under its null p95 of 0.858;
- dedifferentiated (3 genes): 0.688 and 0.518;
- injury-responsive (8 genes): 0.125;
- proliferating (unchanged): 0.643.

Removing CM-identity genes leaves too few genes to form an axis.

## What this means for the direction

- **Direction A, as built from Li SD2, fails clearly.** For dedifferentiated and border-zone, the
  published lists are cluster-versus-all-cells markers. They hold almost no CM-internal genes,
  and what they mark in our CMs is CM identity, depth and heart-level differences.
- For proliferating and injury-responsive, the sets are clean, but our hearts do not show the
  predicted time contrast beyond random genes, except Hr9 (and partly Hr8).
- **The signal sits in a few hearts, not in timepoints.** Hr9 and Hr8 are high across unrelated
  axes. A heart-level factor is the likely cause (injury size, dissection, Hr8's morphine-study
  origin, or library quality). That points to **direction C's question** ("what structure do these
  ventricular CMs hold?"), with count splitting, and to per-heart covariates.
- A different dataset with more border-zone CMs (sorted *nppa*+ cells, GSE139218) is the other
  candidate.

## Weaknesses

- Reproduced byte-identically by a rerun after a style-only lint fix (R-E20).
- With 4 control hearts, AUC is coarse, and the null AUC p95 varies from 0.63 to 0.82 between sets.
- The CM-internal filter in the fixed rule was too weak, a known risk I named in C006 but did not
  fix before fixing the rules.
- The border-zone rule (markers of at least 3 of 5 domains) is one of several defensible choices.
  The others were not tried.
- The Hr9/Hr8 factor is unexplained.

## Status

**Review requested.** The verdict under the fixed rules is a clear failure of direction A as
built. The exploratory variant is reported separately.
