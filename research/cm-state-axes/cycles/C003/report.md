# C003 report: Probe (researcher)

Written 2026-11-07 15:45 CET. Rules fixed at 15:36, before computing: [reading-rules.md](reading-rules.md).
The numbers come from the files named in each row. "Exploratory" marks anything not in the rules.

## C002 corrections carried in

- **F6 (P7 wording):** the Hr11 *nppa* excess is largest in shallow cells and smaller in deep
  cells. It does not disappear. With C003's fixed bins, the ratio in the top bin is 1.04 and the
  supervisor's bins give 1.53, so the size of the top-bin effect depends on the binning.

## Q1. Positive control (F5)

Depth-controlled score: score_genes on cells thinned to 500 UMIs. Tail fraction: share of a
heart's cells above the 99th percentile of 100 matched null sets (null expectation 0.01).
Hearts: 4 Ctrl, 4 at 3 dpi, 7 at 7 dpi, 3 at 30 dpi. Per-heart values are in
[c003_heart_scores.csv](../../../../results/c003_heart_scores.csv); contrasts are in
[c003_positive_control.csv](../../../../results/c003_positive_control.csv).

| Axis | Predicted contrast | Tail fraction (primary): AUC, difference of heart medians | Median thinned score: AUC, difference | Pass? |
|---|---|---|---|---|
| border-zone | 7 dpi > Ctrl | 0.571, +0.003 | 0.500, +0.004 | **No** |
| proliferating | 7 dpi > Ctrl | 0.464, −0.000 | 0.821, +0.022 | **No** on the primary; yes on the secondary |
| injury-responsive | 3 dpi > Ctrl | 0.438, −0.006 | 0.375, +0.001 | **No** |
| dedifferentiated | 3 dpi > Ctrl | 0.562, +0.003 | 0.562, +0.054 | **No** |
| dedifferentiated | 7 dpi > Ctrl | 0.464, 0.000 | 0.464, −0.006 | **No** |

Excluding Hr11 changes no verdict.

**Claims.**

- **Q1-a.** No axis passes its fair positive control on the tail statistic. Tail fractions mostly
  sit near or below the null's 1%; for example, border-zone ranges from 0.000 to 0.024 per heart.
  The C002 failure is therefore not only the P90 statistic's insensitivity. With these sets and
  depth-controlled scores, the predicted time pattern is not visible per heart.
- **Q1-b.** The proliferating "pass" on the secondary median is carried by 3 of 7 hearts. Hr6, Hr8
  and Hr9 have median thinned scores of 0.35–0.37; the other four 7 dpi hearts and all controls
  sit between −0.07 and −0.10. Those three hearts also have the highest *ccng1* detection (0.79–0.84;
  others 0.53–0.72; [c001_marker_detection.csv](../../../../results/c001_marker_detection.csv)).
  *ccng1* was flagged in C001 for its conflicting meaning (a p53 or stress target). This is an
  inference: no leave-one-gene-out was run. Treat the pass as unverified.
- **Q1-c.** Control heart H8 (pooled atrium and ventricle libraries) has the highest tails:
  injury-responsive 0.180, dedifferentiated 0.026, border-zone 0.019. A chamber effect is the
  likely reason, since *nppa* and *nppb* are chamber-biased and H8a is atrium-only. That is a
  reason to control for chamber before any time contrast. It was not tested.
- **Limit (F5 not fully closed):** this control cannot tell "markers carry no state" from "these CMs
  hold too few border-zone or proliferating cells to move a heart-level statistic". A positive
  control with a known strong time effect (a gene set validated in cryoinjury single cells) is
  still missing.

## Q2. The Hr11 excess: ambient RNA or 3 dpi biology?

Shallow bins are 0 (<600 UMIs) and 1 (600–900). Ratios compare each group with Ctrl in the same
bin ([c003_ambient_by_bin.csv](../../../../results/c003_ambient_by_bin.csv)). Rule components are
in [c003_ambient_verdict.json](../../../../results/c003_ambient_verdict.json).

| Test (fixed rule) | Result | Points to |
|---|---|---|
| *nppa* ≥ 2× in shallow bins, in the other 3 dpi hearts | Hr12 1.77 and 2.09; Hr24 0.11; Hr27 1.63. **0 of 3** meet the bar. Hr11: 2.67 and 3.10 | Hr11-specific by the rule; Hr12 is borderline |
| Non-CM genes ≥ 2× in Hr11's shallow bins | **10 of 12** (fibroblast *col1a1a* 20.4 and 16.5, *postnb* 12.4 and 48.4, *fn1b*; endothelial *cdh5*, *kdrl*, *fli1a*; immune *lyz*, *mpeg1.1*; *tcf21*). Blood genes are not raised | Foreign RNA in Hr11 CMs |
| …the same in the other 3 dpi hearts | Hr12 2, Hr24 5, Hr27 7 of 12 | Non-CM RNA is raised in **all** 3 dpi hearts, most in Hr11 |
| *nppa* ratio falls with depth; excess molecules within 2× | Ratio by bin 2.67, 3.10, 2.36, 1.04: it does not fall steadily. Excess molecules 1.53, 2.16, 2.43, 0.28: not within 2× | Not the fixed-amount ambient pattern |
| **Verdict by the fixed rule** | **Inconclusive.** Hr11-specific and non-CM raised, but the depth pattern of *nppa* fails | |

Further evidence, not part of the verdict rule:

- **Atrial soup does not explain *nppa*.** *myh6*, the atrial myosin, is 1.13 and 1.58 in Hr11's
  shallow bins. An atrial CM soup rich enough to triple *nppa* should raise *myh6* too. *nppa* is a
  CM gene, so fibroblast or endothelial soup cannot supply it either.
- **Abundant CM genes are diluted in Hr11** (*myl7* 0.24–0.27, *gapdh* 0.23–0.25, *tnnt2a*
  0.40–0.46). This fits a large foreign-RNA fraction in Hr11 CMs. *hsp70l* is 5.89 and 4.83 in
  Hr11 only (Hr12 0.84 and 0.51), a heat-shock or stress gene, which suggests Hr11-specific stress.
  Mito share including NC is 1.48× in Hr11's shallow bins, and 2.1× in Hr24 and Hr27.
- **Exploratory, cell level:** within each library, a cell's *nppa* share does not rise with its
  non-CM share. The Spearman rho is −0.05 in Hr11, −0.43 in Hr12 and −0.13 to −0.20 elsewhere
  ([c003_exploratory_cellcorr.csv](../../../../results/c003_exploratory_cellcorr.csv)). The *nppa*
  excess does not travel with the contamination.
- ***ttn.2*** (hinted, FU-6): 3.58 and 4.18 in Hr11, but also 2.38 and 3.63 in Hr12 and 6.11 in Hr27.
  Its excess is **shared by the 3 dpi hearts**, not Hr11-specific.

**Which explanation the evidence favours.** The fixed rule gives *inconclusive*. My reading
(inference) has three parts:

- Hr11 has a strong technical signature: foreign non-CM RNA, CM-gene dilution, heat-shock stress
  and higher mito share.
- The *nppa* excess is probably **not** that ambient RNA. The likely soups (non-CM, atrial) do not
  carry it, and it does not co-vary with contamination in single cells.
- The *nppa* excess looks more like a 3 dpi response exaggerated in Hr11, since Hr12 is about 2×
  and Hr27 1.6×. For *ttn.2*, the shared 3 dpi excess points to time, not to one library.

**The trap the slide names** (a marker that tracks a technical artefact) is better matched by the
non-CM and stress genes than by *nppa*. Any axis containing fibroblast or endothelial-expressed
genes (*mmp14b*, which C001 already flagged) or stress genes (*hsp* genes; injury-responsive
includes *hspb6* and *hspb11*) is exposed to 3 dpi contamination.

**What would overturn this.** An empty-droplet soup profile for Hr11 rich in *nppa* would; that
needs the raw matrices, which we don't have. So would a doublet score showing that the *nppa*-high
Hr11 cells are CM–CM doublets. Or a leave-one-library-out showing Hr12's 2× vanishes with
better depth matching.

## Requirement coverage

- R2.3 a tail statistic: done. F5 positive control: run, all axes fail on the primary statistic;
  the control's own sensitivity is unproven.
- R3.3 ambient proxy per library: done for the 3 dpi hearts and controls (non-CM panel and mito
  including NC). Not done: the 7 and 30 dpi libraries, and doublets (R3.2).
- G4 co-occurrence: not started. Given Q1, scoring these sets for co-occurrence would measure
  shared genes and contamination.

## Weaknesses

- One thinning seed. A rerun reproduced all four result files byte-identically (R-E14).
- Shallow-bin ratios for Hr24 and Hr27 rest on bin 1 only.
- The non-CM genes have low base rates, so their ratios are noisy, and `inf` appears where Ctrl is 0.
- The heart contrasts have 4 Ctrl hearts, so AUC is descriptive.
- Chamber (AV) is not controlled.
- The *ccng1* explanation of the proliferating pass is untested.

## Next step proposed

After the freeze, for the demo and README: state the trap result honestly. The agent found
contamination and stress in 3 dpi libraries, strongest in Hr11, on its own. The *nppa* lead was
found unaided, but its cause is unresolved. The *ttn.2* result was hinted. Future work:
1. A positive control with a validated set.
2. Leave-one-gene-out per axis (*ccng1*).
3. A chamber-controlled contrast.
4. Doublet scores.

## Status

**Review requested.** Q1 and Q2 were both answered under the fixed rules: Q1 negative, Q2
inconclusive, with a stated reading.
