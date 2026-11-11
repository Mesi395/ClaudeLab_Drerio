# Statistical plan v2

Written 2026-11-07, after the hackathon. v1 grew one test at a time (results/FINDINGS.md items
1-17). This plan starts from the design, sets one primary test per question, and states what each
earlier analysis becomes. A senior design review (Fable, 2026-11-07) changed the first draft: the
chemistry term, the handling covariate and the omnibus gate were removed, 7 dpi became the primary
contrast, and ambient and gene-overlap checks were added. Sections 2A, 3 and 5 are fixed before the
v2 run (scripts/pseudobulk_set_plan.py is committed with this file and has not been run).
Sections G and H were added after a literature survey of statistical strategies, each committed
with its script (scripts/strategy_benchmark.py, scripts/high_confidence.py) before its first full run.

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
- **A label discrepancy, noted.** Hu et al. Supplementary Data 1 lists Hr1, Hr2a and Hr2b at 0 dpi;
  their GEO records (GSM4817933, GSM4817934, GSM4817938: "Injured heart ... dpi: 7"), the task data
  and our config say 7 dpi. We keep 7 dpi and flag it for the task owner.

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

### G. Which strategy? Candidates chosen from the literature, then tested on this data (added 2026-11-07, before its run)
Sections A and F were chosen by reasoning from the design and from benchmarks on other datasets
(Squair et al. 2021, Crowell et al. 2020). G makes the choice on this data. The candidate list comes
from a literature survey (research/s4_* in the run record): four reports, two on methods papers and
two on what similar injury single-cell studies did (10 zebrafish heart studies, mammalian MI studies).
What the survey found:
- Similar studies: 9 of 10 zebrafish heart studies (Hu 2022 included) test cells, not hearts; 5 of 10
  pool hearts per library, so they cannot test hearts at all. In mammalian MI only Kuppe 2022 and the
  Squair-lab studies treat the animal or patient as the unit. None uses an "in at least u of n animals"
  rule, none asks whether every cell shifts or a subgroup switches on, none corrects for chemistry.
- Methods papers: pseudobulk beats per-cell tests for false positives (Squair 2021, Murphy & Skene
  2022, Junttila 2022); mixed models inflate false positives at 5 animals per group (Zimmerman 2021,
  He 2021), and we have 3-4; per-sample scores (z-score, PLAGE) can be unstable at small n, singscore
  scores each sample alone (Foroutan 2018); rotation tests avoid the coarse p-values of relabelling at
  small n (Wu 2010); limma fry with group-specific variance weights (voomByGroup, You 2023) is the one
  option the survey rates above ours at this n, and it needs R, which neither host has.
- Not candidates, with the reason: Milo and CNA (abundance, not programme; CNA cannot permute across
  the 3 dpi chemistry batch), propeller (proportions of clusters we do not call), NEBULA / MAST mixed
  (the M2 class), IDEAS and waddR (distribution tests; G2 checks our 2F test instead), batch correction
  (chemistry is nested in time point, so it cannot be separated, Leek 2010, Nygaard 2016).

**Null and planted signal.** The seven 7 dpi ventricle hearts split 3 vs 4 in all 35 ways: real
heart-to-heart variation (including the second-series heart Hr8), no condition difference. A known
signal is planted by binomial thinning (Gerard 2020): the set genes' counts are thinned to 1/f in
every cell except the cells that carry the signal. Shift: every cell of the 3 "case" hearts carries it
(f = 1.25, 1.5, 2). Subgroup: a share pi of the case cells carries it (pi = 0.1, 0.2; f = 3).
3 vs 4 matches the 3 and 30 dpi arms; the primary 7 vs 4 contrast has more hearts, so its power is
higher than measured here.

**G1, heart-level shift (question A).** Candidates, one-sided, alpha 0.05:
| | strategy | used by | assumption that can fail here |
|---|---|---|---|
| M1 | per-cell Wilcoxon on a cell set score | 9 of 10 zebrafish heart studies surveyed (incl. Hu 2022); mouse MI (Yamada 2022, Calcagno 2022, Cui 2020) | cells of one heart are independent |
| M2 | per-cell mixed model, heart random intercept | Tabulae Paralytica (lmerTest on module scores), Kathe 2022 (NB mixed model), Zimmerman et al. 2021 | a normal random effect estimated from 7 hearts |
| M3 | per-heart z-score, Welch t, t-distribution p | sample-level tests in Kuppe 2022, Nicin 2022 | normality at 3 vs 4 hearts |
| M4 | per-heart z-score, Welch t, exact relabelling | plan 2A (current choice); Squair et al. 2021 (pseudobulk) | exchangeable hearts; floor 0.029 at 3 vs 4 |
| M5 | competitive: set genes vs expression-matched decoy genes | v1 item 13; camera-type tests | genes independent (Wu & Smyth 2012) |
| M6 | rotation gene-set test (ROAST-style) on pseudobulk log expression | Wu et al. 2010 (limma `roast`/`fry`) | normal residuals; 5 residual df; unweighted |
| M7 | per-heart PLAGE score (first singular vector), exact relabelling | Tomfohr et al. 2005; GSVA package | one dominant direction across 7 hearts |
| M8 | per-heart singscore (mean within-heart rank), exact relabelling | Foroutan et al. 2018 | rank ties among low genes |
False-positive rate: share of null tests with p <= 0.05 over 35 splits x (6 primary + 50 random
10-gene sets), with a Clopper-Pearson 95% interval (splits share hearts, so the interval is too
narrow). Power: share of 35 splits x 6 primary sets detected per planted scenario.
**Decision (fixed):** a candidate is calibrated if its false-positive interval reaches 0.05 or lower.
The chosen strategy is the calibrated candidate with the highest mean power over the three shift
scenarios; within two percentage points the one with fewer distributional assumptions wins
(M4, M8, M7, M3, M6, M5, M2, M1). Mixed-model fits that do not converge count as not detected and are
reported. If the choice is not M4, section A is re-run with it and FINDINGS item 18 is restated.

**G2, every cell or a subgroup (question F).** Plan 2F's shape test on the same splits and plants
(raw exact p per set). Reported: the share of tests called "a subgroup switches on" per set and
scenario. 2F's "every cell shifts" is an absence-of-evidence call, so it needs power (survey s4_m2):
it stands for a set only if that set's subgroup-call rate is >= 0.8 under both subgroup plants and
<= 0.1 under the null and the shift plants; otherwise item 19's call for that set becomes "cannot
tell at 600 UMIs". A shape p that cannot be computed (no set gene detected in some heart after
thinning) counts as no call, so such a set cannot pass.

Synthetic check before the run (scripts/strategy_benchmark.py functions, 300 simulated 7-heart
datasets each): no effect and independent genes, ROAST 0.060 and competitive 0.057 at p <= 0.05;
no effect with a shared heart effect on the set genes, ROAST 0.043 and competitive 0.283.
Two-split smoke runs on the real data checked only that the code runs.

### H. High-confidence presence without asking every heart to match (added 2026-11-07, before its run)
The task owner: cutoffs need not demand that every library matches; define cutoffs of high
confidence, and find where the states form. In the similar zebrafish and mammalian injury studies
(literature survey, research/s4_d1 and s4_d2 in the run record) a state is "present" when a cluster
exists and staining in about three animals agrees; replication across animals is judged by eye and no
study uses a formal "in at least u of n animals" rule. The statistical structure for that statement is
partial conjunction (Benjamini & Heller 2008; Wang, Gui, Su, Sabatti & Owen 2022), with a lower
confidence bound on the number of animals (Wang & Owen 2019). Three levels, each with its own null
(scripts/high_confidence.py):
- **cell:** "high" if its set score (mean log-normalised count over the set) is above the 95th
  percentile of the scores of 1,000 detection-matched decoy gene sets in the same heart and chamber
  (cells thinned to 600 UMIs, veto applied; hearts with fewer than 30 cells left are skipped);
- **heart:** the share of high cells against each decoy set's own share; p = (1 + decoys at or above
  the set) / 1001, floor 0.001. The null is competitive: a positive heart means "more high cells than
  random genes of the same detection", not "injury state";
- **condition:** Simes partial conjunction for "high cells in excess in at least r of the n hearts",
  P(r/n) = min over r <= i <= n of (n - r + 1) / (i - r + 1) * p(i) (hearts are independent animals),
  and r_hat = max{r : P(r/n) <= 0.05}, a 95% lower bound on how many hearts carry the state.
Why 1,000 decoys: with u hearts at the floor and the rest null, P = (n - u + 1) x floor; at 200 decoys
a lone set could pass BH over six sets only if every heart sat at the floor, at 1,000 any u can.
Why Simes over Fisher: Fisher's combination answers "at least one heart" more strongly than "at least r"
and is dominated by one extreme heart; Simes is the form the replicability papers give with r_hat.
Reported per set, condition and chamber: r_hat, P at the majority r = floor(n/2) + 1 (the majority the
v1 gate asked for, now one test instead of n separate cutoffs) with BH over six primary sets x four
conditions x two chambers, and the per-heart excess (where the state shows).
**Injury-specific presence** (fixed): r_hat >= 4 of 7 at 7 dpi, r_hat = 0 in the uninjured hearts, and
the section A test positive at 7 dpi. Presence alone is not injury: a set present in uninjured hearts
is baseline variation.
Controls: dissociation set (positive: handling stress is present in every heart); 20 random 10-gene
sets (negative: should reach r_hat = 0 almost always).
Checks before the run (synthetic, 80 hearts each, after the decoy fix): no signal, 6% of hearts at
p <= 0.05 (median p 0.50); 10% of cells in a state, 61%; 20%, 80%. Simes with 4 of 7 hearts at the
floor: 0.004 at 1,000 decoys, 0.020 at 200. The first synthetic check found the decoy bug fixed in
9f349b3 (an undetected set gene was matched by detected genes, so decoys outscored the set wherever
set genes were absent).

## 3. Multiplicity
Separate BH families, one per time point (renamed 2026-11-07: v2 called this a "fixed sequence";
separate BH per family controls the error averaged over families, Benjamini & Bogomolov 2014).
Floors stated beside every test:
1. **7 dpi vs control, six primary sets, BH q ≤ 0.05** (the primary family).
2. **3 dpi vs control, six primary sets, BH q ≤ 0.05**, all six pre-declared; floor 1/70 = 0.014. A lone set
   at the floor gets q = 6/70 = 0.086 and cannot pass; at least two sets at the floor are needed
   (q <= 0.043; corrected 2026-11-07, v2 said "about 0.04"); every result carries the confounding label.
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

Added with sections G and H (checked against Crossref 2026-11-07):
- Benjamini Y, Heller R (2008) Biometrics 64:1215-1222. doi:10.1111/j.1541-0420.2007.00984.x
- Wang J, Owen AB (2019) J Am Stat Assoc 114:158-168. doi:10.1080/01621459.2017.1385465
- Wang J, Gui L, Su WJ, Sabatti C, Owen AB (2022) Ann Stat 50:1890-1909. doi:10.1214/21-aos2139
- Gerard D (2020) BMC Bioinformatics 21:206. doi:10.1186/s12859-020-3450-9
- Tomfohr J et al. (2005) PLAGE. BMC Bioinformatics 6:225. doi:10.1186/1471-2105-6-225
- Foroutan M et al. (2018) singscore. BMC Bioinformatics 19:404. doi:10.1186/s12859-018-2435-4
- Lee E et al. (2008) PLoS Comput Biol 4:e1000217. doi:10.1371/journal.pcbi.1000217
- Zimmerman KD et al. (2021) Nat Commun 12:738. doi:10.1038/s41467-021-21038-1
- Murphy AE, Skene NG (2022) Nat Commun 13:7851. doi:10.1038/s41467-022-35519-4
- Junttila S et al. (2022) Brief Bioinform 23:bbac286. doi:10.1093/bib/bbac286
- He L et al. (2021) NEBULA. Commun Biol 4:629. doi:10.1038/s42003-021-02146-6
- You Y et al. (2023) voomByGroup. Genome Biol 24:107. doi:10.1186/s13059-023-02949-2
- Leek JT et al. (2010) Nat Rev Genet 11:733-739. doi:10.1038/nrg2825
- Nygaard V et al. (2016) Biostatistics 17:29-39. doi:10.1093/biostatistics/kxv027
