# Statistical plan v2

Written 2026-11-07, after the hackathon. v1 grew one test at a time (results/FINDINGS.md items
1-17). This plan starts from the design, sets one primary test per question, and states what each
earlier analysis becomes. A senior design review (Fable, 2026-11-07) changed the first draft: the
chemistry term, the handling covariate and the omnibus gate were removed, 7 dpi became the primary
contrast, and ambient and gene-overlap checks were added. Sections 2A, 3 and 5 are fixed before the
v2 run (scripts/pseudobulk_set_plan.py is committed with this file and has not been run).

## 1. The design, and what it can and cannot identify

Unit of replication: the heart. Cells from one heart are not independent replicates, so every
injury claim is tested with hearts as the unit (pseudobulk; Squair et al. 2021, Crowell et al.
2020). Cells only inform estimates within a heart.

Ventricle pseudobulks with at least 30 cells (18 hearts):

| condition | hearts | chemistry | source series |
|---|---|---|---|
| uninjured | H5, H6, H7, H8 | 10x v2 | main |
| 3 dpi | Hr11, Hr12, Hr24, Hr27 (47 cells) | Hr24, Hr27: 10x v3; others v2 | Hr11: second series |
| 7 dpi | Hr1, Hr2, Hr6, Hr7, Hr8, Hr9, Hr13 | 10x v2 | Hr8: second series |
| 30 dpi | Hr3, Hr4, Hr19 (55 cells) | 10x v2 | main |

Sources: chemistry from the GEO sample records of GSE159032 ("10x chemistry version: v3" for Hr24
and Hr27, v2 for every other heart used here). The second series is GSE158919, a morphine pain study
from the same lab; Hr11 and Hr8 are its untreated hearts ("morphine treatment: no"), so they carry
no drug but come from a separate experiment. Sequencing depth does not separate v2 from v3 hearts,
so chemistry is known only from metadata.

What follows:
- **7 dpi vs control carries the primary injury claim**: 7 hearts vs 4, one chemistry, one heart
  from the second series. Its exact relabelling floor (1/330) leaves room for multiplicity.
- **3 dpi vs control is confounded.** Chemistry is nested inside 3 dpi, so a chemistry term does
  not adjust anything: with it, the 3 dpi estimate is exactly Hr11 and Hr12 against the controls.
  Every 3 dpi number is labelled confounded, and the two-v2-heart estimate is shown beside it.
- **30 dpi** has 3 hearts (one with 55 cells): estimates with intervals, no calls of equality.
- **The control arm is the least variable** (per-heart score SD 0.1-0.5 vs 0.4-1.0 after injury),
  so group variances are not pooled (Welch statistic).
- **Two heart-level nuisance axes exist**: handling stress (the dissociation / heat-shock set; high in
  Hr11, Hr24, Hr3, Hr9) and ambient fibroblast RNA (col1a1a and postnb rise 2^4-2^5 in the 3 dpi
  cardiomyocyte pseudobulk). Neither can be regressed out: the handling set shares genes with the
  sets it would correct, and handling also rises with injury itself (a mediator). They are handled by
  gene splits and a correlation criterion (2A).

## 2. Questions and their tests

### A. Does injury shift a state's programme? (heart level, primary)
- **Measure:** one score per heart and gene set (`set_scores` in scripts/pseudobulk_set_score.py:
  median-of-ratios size factors, log2, z per gene across the 18 hearts, mean over set genes). No step
  uses the labels, so the relabelling null is exact. The z-scaling compresses very large shifts; it
  does not bias the test.
- **Primary sets (six):** memory border-zone core, border zone, dedifferentiated, injury-responsive,
  proliferating, glycolysis (hk1, pkma, pdk2a).
- **Gene splits** (they decide what a positive set means, not whether it is positive):
  dedifferentiated without its nine glycolytic genes (hk1, pkma, pdk2a, ldha, pgk1, eno1a, tpi1b,
  pgam2, gapdh), and also without its four library-size genes (ttn.1, ttn.2, bves, synpo2lb);
  memory core without AP-1 genes; injury-responsive without its dissociation gene (socs3a).
  Overlap is declared, not counted as replication: glycolysis is inside dedifferentiated; junba,
  junbb, fosab are in both the memory core and the handling set; nppb and mustn1b are in the memory
  core and border zone.
- **Controls:** handling (dissociation set), ambient fibroblast (col1a1a, col1a2, postnb, dcn),
  pdk2a alone, and 200 random 10-gene sets from genes detected in every heart. The random sets
  are a genome-wide shift diagnostic, not a false-positive rate: the share of random sets with
  p ≤ 0.05 at a contrast says how much the hearts differ overall.
- **Statistic:** Welch t, one-sided (injured above control); p by exact relabelling of hearts,
  observed labelling included. Effect: difference in mean score with a 95% Welch interval.
- **Contrasts:** 7 dpi vs control (330 relabellings, floor 0.003); 3 dpi vs control (70, floor
  0.014); 3 dpi Hr11+Hr12 vs control (15, floor 0.067, estimate only); 30 dpi vs control (35, floor
  0.029); 3 vs 30 dpi (35, floor 0.029).
- **Criteria attached to any positive set** (reported as qualifiers, not gates):
  1. handling: within-condition correlation of the set score with the handling score
     (both centred on their condition means, 18 hearts); r ≥ 0.6 means the set cannot be separated
     from handling with this design;
  2. ambient: the same with the ambient fibroblast score; r ≥ 0.6 means it cannot be separated from
     fibroblast RNA;
  3. one heart: leaving any single heart out changes the estimate's sign or by more than half;
  4. depth: the call is repeated on counts thinned to 600 UMIs per cell;
  5. small hearts: the call is repeated without the two hearts under 100 cells (Hr27, Hr19).
- **Chamber line (robustness, 3 and 7 dpi):** per heart, ventricle score minus atrium score (sets
  scored over all chamber pseudobulks with at least 30 cells). Additive heart-level nuisances
  (handling, chemistry, series, depth) cancel. It measures the ventricle-local part only, so low
  power is expected; it supports, it does not decide.

### B. Which genes change? (gene level, supporting)
DESeq2 on the ventricle pseudobulks, `~ condition`, 7 dpi vs control as the primary contrast, 3 dpi
reported with the confounding label and a v2-only refit. Conclusions about states rest on A.

### C. Can a state be called cell by cell? (cell level)
A negative per-cell result means something only with its power. Plant each signature into a fraction
of real cells thinned to 600 UMIs by scaling the real set genes (gene-gene correlation kept), at the
per-heart decoy-calibrated threshold actually used; report power by heart depth, set size and
fraction. State-level test without per-cell labels: Milo neighbourhood differential abundance with
hearts as replicates, `~ condition`, plus a v2-only sensitivity (Dann et al. 2022).

### D. Do two states co-occur?
Only where C shows power. Within-heart association of two scores at neighbourhood or metacell level,
decoy-pair null, the heart as the unit, and the smallest detectable correlation reported.

### E. Which annotation to trust?
Without ground truth, an annotation earns trust on four criteria, in this order:
1. **Calibration:** positive and negative controls behave (the handling set is found where handling
   differs; decoy sets are not).
2. **Replication:** the call holds across hearts and survives leaving one heart out.
3. **Specificity:** it follows injury time and is separable from handling, ambient RNA and chemistry.
4. **Independent evidence:** it agrees with data that share no genes or cells with it, such as
   spatial (tomo-seq) border-zone positions, in situ data or lineage data.

Also required: a table of set overlaps (above); a statement of how cardiomyocytes were selected and
that the filter uses no gene of a tested set; the published annotation's own markers scored per
heart as a sanity check of that annotation; one table of every set, control and floor, including
the non-significant rows. Labellers count as agreeing only if they share no marker genes; the
latent-class model assumes independent errors, so kappa 0 by construction is not evidence.

### F. Does injury shift every cell a little, or switch on a subgroup? (added 2026-11-07, before its run)
The question from the task owner's explanation: a population mean shift, or a high tail that grows.
Per ventricle heart, cells at 600 UMIs or more thinned to exactly 600, two shape measures per set
against 200 detection-matched decoy sets drawn within that heart: tail excess (share of cells above
the decoys' 95th percentile, minus the decoys' share) and dispersion excess (log variance-to-mean of
the set's summed counts across cells, minus the decoys' median). Both ignore a uniform level shift,
which section A measures. Same test, contrasts and BH families as A
(scripts/cell_tail_test.py). Decision: "a subgroup switches on" if tail excess q <= 0.05 and
dispersion excess p <= 0.05 at the same contrast; "every cell shifts" if A found a shift and the
dispersion excess did not move; otherwise "cannot tell at 600 UMIs".

## 3. Multiplicity
Fixed sequence, floors stated beside every test:
1. **7 dpi vs control, six primary sets, BH q ≤ 0.05** (the primary family).
2. **3 dpi vs control, six primary sets, BH q ≤ 0.05**, all six pre-declared; floor 0.014, so the
   best reachable q is about 0.04; every result carries the confounding label.
3. **30 dpi vs control and 3 vs 30 dpi:** estimate, interval and one-sided exact p; no multiplicity
   claim and no "persistent" category (3 hearts cannot show equality).
Gene splits, controls and the chamber line are reported with raw p; they qualify calls, they do
not add calls.

## 4. What the v1 analyses become

| v1 analysis (FINDINGS item) | status in v2 |
|---|---|
| per-cell protocol: UCell, decoys, per-heart gate (1-9) | supporting: specificity calibration; power added in C |
| 4 vs 4 DESeq2 at 3 dpi (10-11) | superseded by B (7 dpi primary, 3 dpi labelled) |
| competitive decoy set test (13) and its ordering | retired: it treats genes as independent (handling z 9.1 there, p 0.14 when hearts are relabelled) |
| relabelling calibration of the DEG count (14) | supporting |
| 30 dpi read by eye (15) | replaced by the 30 dpi estimates of A |
| leave-one-heart-out Jaccard and rank stability (12, 16) | retired: at 8 hearts it cannot separate signal from noise |
| per-heart score, pairwise exact relabelling (17) | the test of A, extended to 7 dpi |
| within-heart partial Spearman co-occurrence | supporting, until D adds power |
| latent-class agreement | retired as a measure of trust (E) |

## 5. Decision rules for A (fixed before the v2 run)
- **Injury shift at 7 dpi:** step 1 q ≤ 0.05.
- **Shift at 3 dpi (confounded):** step 2 q ≤ 0.05, always with the label and the Hr11+Hr12 estimate.
- **Cardiomyocyte programme:** a shift whose handling and ambient correlations are both below 0.6.
  Otherwise "shift, cannot be separated from handling / ambient RNA".
- **Rests on one heart:** the leave-one-heart-out flag fires. **Depth- or size-sensitive:** the
  thinned or small-hearts-out repeat misses p ≤ 0.05.
- **Dedifferentiation beyond glycolysis:** claimed only if dedifferentiated without its glycolytic
  genes reaches p ≤ 0.05 at the same contrast; otherwise the dedifferentiated result is reported as
  the glycolytic shift.
- **Returns toward control by 30 dpi:** reported as the 3 vs 30 dpi estimate and interval, with
  "fades" only when the one-sided exact p reaches its floor of 0.029.
- Anything else: **not shown with this design**, with its interval.

## References (checked against Crossref by the design review, 2026-11-07)
- Squair JW et al. (2021) Nat Commun 12:5692. doi:10.1038/s41467-021-25960-2
- Crowell HL et al. (2020) Nat Commun 11:6077. doi:10.1038/s41467-020-19894-4
- Goeman JJ, Bühlmann P (2007) Bioinformatics 23:980-987. doi:10.1093/bioinformatics/btm051
- Wu D, Smyth GK (2012) camera. Nucleic Acids Res 40:e133. doi:10.1093/nar/gks461
- Wu D et al. (2010) ROAST. Bioinformatics 26:2176-2182. doi:10.1093/bioinformatics/btq401
- Janssen A (1997) Stat Probab Lett 36:9-21. doi:10.1016/S0167-7152(97)00043-6
- Dann E et al. (2022) Milo. Nat Biotechnol 40:245-253. doi:10.1038/s41587-021-01033-z
- Risso D et al. (2014) RUV. Nat Biotechnol 32:896-902. doi:10.1038/nbt.2931
- van den Brink SC et al. (2017) Nat Methods 14:935-936. doi:10.1038/nmeth.4437
- Honkoop H et al. (2019) eLife 8:e50163. doi:10.7554/eLife.50163 (hk1, pkma, pdk2a confirmed by
  in situ hybridisation in border-zone cardiomyocytes at 7 dpi, border zone vs remote; our result
  is "consistent with", not a replication)
- Benjamini Y, Bogomolov M (2014) J R Stat Soc B 76:297-318. doi:10.1111/rssb.12028
