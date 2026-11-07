# Research map: cardiomyocyte state axes after cryoinjury, without ground truth

Study folder: `research/cm-state-axes/`. Started 2026-11-07 14:12 CET.
Repository rules (`AGENTS.md`, `CONTRIBUTING.md`) apply on top of this map.

## 1. Request

Word for word, in order:

> hello

> first move the data from the Dcownloads folder to the current workspace
> *(screenshot of ~/Downloads: HR_1.csv, CM_zotero_papers.csv, Integration_zotero_papers.csv,
> singlecellgenomics_zotero_papers.csv, 20261107_ClaudeLabDay.zip,
> Hu2022_CM_Ctrl_CMo…ocessed_allgenes.h5ad, Hu2022_CM_Ctrl_CMo…_processed_HVG.h5ad)*

> here is my  project, and here is the paper related to it https://www.nature.com/articles/s41467-025-59070-0
>
> For cell scoring on the state axis:
> Scanpy's built-in gene scoring function is scanpy.tl.score_genes() (which computes average expression subtraction relative to a reference pool), AUCell is an external ranking-based scoring algorithm originally designed for R.
> Because AUCell is ranking-based, it is completely independent of normalization units and library size. You can run AUCell on a Scanpy AnnData object using two primary Python tools: pyUCell or the OmicVerse wrapper.
>
> Now nanalyse and tell me what we are doing here
> *(attached: the problem slide "Can an agent annotate cell states when there's no ground truth?",
> transcribed in [docs/problem.md](../../docs/problem.md))*

## 2. Goal

**Problem.** After cryoinjury, zebrafish cardiomyocytes are said to take on several states:
injury-responsive, border-zone-like, dedifferentiated and proliferating. These labels come from
clustering and hand-picked marker genes. No single-cell ground truth exists. Two questions are
open. Do the states co-occur in the same cells, or are they separate populations? And what do the
published markers really mark? At least one supposed "dedifferentiation" marker is known to track
library size rather than biology. The dataset makes this risk sharp. Hr11, a 3 dpi library, has
far lower sequencing depth than the others, so any score that follows depth will also look like
an injury response.

**How it is handled today.** The data's authors (Hu et al. 2022) clustered merged libraries
without batch or ambient-RNA correction and labelled cardiomyocyte subtypes by hand. Spatial
studies define the border zone by position: Wu et al. 2016 by tomo-seq along the injury axis,
Honkoop et al. 2019 by sorting *nppa*-reporter cells, and Li et al. 2025 (the paper the user
supplied) with Stereo-seq domains at the wound edge after apex resection. Li et al. report a
stage-specific sequence of border-zone states: activated, dedifferentiated,
dedifferentiated-proliferating, proliferating and re-differentiating. Language-model annotators
(GPTCelltype, CASSIA, LICT) assign labels but rarely test whether a marker reflects a technical
artefact.

**Key idea.** Treat each state as a continuous axis, scored per cell from a marker set sourced
independently of this dataset (literature or spatial data), rather than as a cluster label. Then
each axis has to earn trust through checks the agent runs itself. The two scoring methods use
different assumptions: scanpy `score_genes` (set mean minus an expression-matched control pool)
and a rank-based method (AUCell or UCell), and the two must agree. The other checks are a matched
random gene-set null, partial correlation with library size and genes detected within each
library, doublet scores, an ambient-RNA proxy, and robustness to dropping single genes. The user
proposes that a rank-based score is independent of normalization and library size. The study
treats this as a hypothesis to test. Ranks are unaffected by per-cell scaling, but lower depth
means more zeros, more ties and fewer detected genes, which can still move a rank score. Hr11 is
the natural test case.

**What is to be found and shown.** For each of the four states:

- a sourced marker set, with each marker that tracks a confound flagged rather than silently
  dropped;
- two-method scores for all 32,044 cardiomyocytes;
- null, confound and stability results.

Across states, show whether the axes coincide or come apart within each timepoint, once time,
library and technical covariates are controlled. Then the agent writes its own annotation. A
human compares it with the held-out expert annotation using `src/zf/evaluate.py`, and the study
lists which known traps the agent caught without being told. The two human questions stay in
view: with no ground truth, which annotation should we trust, and which checks did the agent run
on its own versus which would a domain expert add?

**Success, scope and evaluation unit.** Success means every axis has sourced markers, two-method
scores and null, confound and stability results in `results/`, each with its commit, command and
input. Flagged markers carry evidence. The co-occurrence question gets an answer with its scope
and confound controls. A negative or inconclusive answer is acceptable if it is honest. Scope:
cardiomyocytes only, the 22 libraries of Hu et al. 2022 (uninjured, 3, 7 and 30 dpi), and no use
of the authors' labels. The evaluation unit is the single cell on each state axis. The final
agreement check is per cell against the expert annotation and is run by a human only. The
deliverables are the results, the README and a 3-minute demo, due under the hackathon schedule on
2026-11-07: integration 15:30, feature freeze 16:15, submission 17:00. No manuscript was
requested.

### Breakdown

| Part | Requirements | Success means | Evidence needed |
|---|---|---|---|
| **G1 Markers** | R1.1 a marker set per state (4), each gene with a DOI or a result file. R1.2 every gene present in the data under its zebrafish symbol. R1.3 each source's injury model, timepoint and resolution (cell or spot) recorded, so a mismatch with cryoinjury single cells is visible. | Four YAMLs in `proposals/loop/` (was `markers/`), sourced and checked against `var_names`. | Sources opened in full or the access failure stated; a presence check script. |
| **G2 Scores** | R2.1 every cardiomyocyte scored on each axis with `score_genes` and a rank-based method. R2.2 cross-method agreement per axis. R2.3 a matched random gene-set null per axis. | Score table plus agreement and null summaries in `results/`. | Script outputs with provenance. |
| **G3 Confounds** | R3.1 library size and genes detected, within library. R3.2 doublet score. R3.3 ambient-RNA proxy per library. R3.4 library-specific structure (Hr11 and chemistry: Hr24 and Hr27 are v3). R3.5 markers flagged on evidence. R3.6 a direct test of whether the rank-based score is independent of depth. | `results/confounds_*.csv`; flagged genes in the YAMLs. | Within-library statistics, not pooled ones. |
| **G4 Co-occurrence** | R4.1 pairwise correlation and joint-high fractions per axis pair, per timepoint. R4.2 controlled for time, library and confounds. R4.3 a stated answer, coincide or come apart, with its scope. | `results/cooccurrence.csv`, one figure and the answer. | As R4.1 and R4.2, plus a null for joint-high fractions. |
| **G5 Trust without ground truth** | R5.1 which checks the agent ran itself versus what an expert would add. R5.2 an agent label file in the format `evaluate.py` expects (a human runs it). R5.3 no exposure to the held-out labels: no `holdout/`, and the label columns embedded in the supplied `.h5ad` stripped before use. | A list of checks; a label file; a leakage log. | `JOURNAL.md` entries; a sanitised data file with an allowlist of columns. |

Goal history: 2026-11-07 14:12, first version from the request, the slide, `README.md`,
`docs/problem.md`, `docs/data.md` and `TASKS.md`.

## 3. How far we are, per part

After C001 and C002 (researcher; evidence R-E1–R-E13; reports [C001](cycles/C001/report.md), [C002](cycles/C002/report.md)).

- **G1 Markers: drafts exist, evidence uneven.** Four candidate YAMLs in `proposals/loop/` (moved from `markers/` for review), each gene sourced.
  Existing work does not give four separate marker sets: the same genes (*nppa*, *nppb*, *mustn1b*)
  are named for injury response, border zone and dedifferentiation. Border-zone is best supported
  (two or more independent sources, genes detected in a large share of CMs). Injury-responsive genes
  were defined at 6 h after injury, but our earliest sample is 3 dpi. Proliferation genes are
  clear but detected in very few cells. Dedifferentiation is described mostly by proteins and by
  genes going *down* (sarcomeres), so it lacks a positive transcript set. Only one cryoinjury
  single-cell source was read (Honkoop 2019, 352 cells).
- **G2 Scores: done for the four draft axes (C002).** Every CM is scored with score_genes and UCell.
  The two methods agree almost perfectly in every library, so their agreement is no independent
  check. No axis beats random genes of matched expression on its high tail. Only border-zone comes
  close to showing its genes co-vary beyond random sets.
- **G3 Confounds: depth done; doublets and ambient RNA not started.** The rank score is *more* tied
  to depth than score_genes. Thinning cells to equal depth removes most of that tie, so it is
  technical. The user's claim that a rank score is independent of depth does not hold here. Ten
  axis genes are flagged as depth-tracking. Hr11's excess *nppa* survives depth matching but is
  concentrated in low-count cells, which fits ambient RNA (untested). The mito fraction misses the
  `NC-002333.*` mitochondrial features.
- **C003 positive control and ambient test.** No axis shows the change the biology predicts over
  time, counted per heart (18 hearts) on a depth-controlled tail statistic. Proliferating passes
  only on a secondary median, and that rests on three hearts, probably through the stress-linked
  gene *ccng1*. All 3 dpi libraries carry foreign non-CM RNA, Hr11 the most, together with
  dilution of CM genes, heat-shock stress (*hsp70l*) and higher mito share. The Hr11 *nppa* excess
  is **not** well explained by that contamination: the likely soups don't carry *nppa*, and it
  doesn't co-vary with contamination per cell. The fixed rule rates it inconclusive. *ttn.2*
  (hinted) is raised in all 3 dpi hearts.
- **C004 chamber and labels.** A chamber call from *myh6* against *myh7*+*myh7l* recovers the six
  chamber-only libraries with balanced accuracy 0.991. It is optimistic, and *myh7l* is
  unsourced. In ventricular CMs only, injury-responsive (3 dpi) and dedifferentiated (7 dpi) pass
  the tail rule, but weakly: injured hearts sit near the null's 1% while control hearts sit below
  it, and their medians move the wrong way. Border-zone just misses (AUC 0.786). The agent label
  file exists (R5.2). Most cells are unassigned, and three of the four state labels are at chance
  level.
- **C005 label fix (F10).** The corrected label file scores ventricle-called cells only. The
  injury-responsive calls drop from 1,188 to 145, because most c004 calls were atrial. After the
  correction, every state class is at or below chance (`results/c005_label_counts_vs_c004.csv`).
- **C005 doublets (R3.2).** Scrublet on all 22 GEO matrices: no axis score tracks the doublet score
  within libraries, and no state is enriched among high scores. Scrublet's own doublet call failed
  (0 doublets), so the share of states among called doublets is not answered.
- **C007, direction A tested.** Sets built from Li 2025 SD2 (25 genes per axis, disjoint,
  non-CM genes removed) do not beat the C004 sets. No axis passes a control-calibrated,
  null-tested positive control in ventricular CMs, and the two weak C004 passes vanish under
  calibration. Li's dedifferentiated and border-zone markers are CM-versus-all markers, mostly
  sarcomere and OXPHOS genes, with almost no CM-internal genes. The signal sits in a few hearts:
  Hr9 and Hr8 are high on unrelated axes. It does not sit in timepoints.
- **G4 Co-occurrence: cannot be resolved with these sets (C004).** Overlapping sets co-occur by
  construction. Disjoint sets in ventricular CMs are indistinguishable from random gene pairs:
  every pair is "mixed". Scoring the axes as they stand would mostly measure shared
  genes and depth.
- **G5 R5.3: done.** A tested sanitiser keeps only the allowlist; no labels seen.
- Overall: no existing work, nor a combination, already scores these four states in cryoinjured
  single CMs. The gap is real at the "separate the states" step.

## 4. Existing work, by approach

**Cryoinjury single cells.** Honkoop 2019 (eLife, doi:10.7554/eLife.50163): sorted
*nppa*-reporter CMs at 7 dpi. Border zone = *nppa*, *vmhc* (now *myh7*), *mustn1b*, glycolysis genes;
mitochondrial and OXPHOS genes down. Read: first ~45% of full text. Serves G1 border-zone. Reusable.
Constanty 2025, Forman-Rubinsky 2025 (snRNA-seq): abstracts only.

**Spatial, spots or sections.** Li 2025 (Nat Commun, resection, Stereo-seq): five stage-specific
border-zone domains with named genes; only *tnfrsf11b* was checked in cryoinjury. Read: first ~35%
of full text twice; supplements not opened. Wu 2016 (Dev Cell, tomo-seq, cryoinjury 3 and 7 dpi): two
border zones, BMP targets; read via a search summary only. Both serve G1, but the transfer to single
cryoinjury cells is unproven.

**Reporter and functional studies.** *ankrd1a* (Boskovic 2026), *mustn1b* (Meena 2026), *foxm1/cenpf*
(Zuppo 2023), *gata4* (Kikuchi 2010), embryonic myosin protein (Sallin 2015), sarcomere silencing
(Ben-Yair 2019). Abstracts or titles only. They give genes, not single-cell scores.

**Scoring and confound tools.** See [methods-and-tools.md](../../literature/methods-and-tools.md).
No source read claims rank scores are fully independent of depth. pyUCell is UCell, not AUCell.

**Li 2025 supplementary tables (opened in C006).** Full marker tables for chamber, scRNA-seq
regeneration clusters (including "Cardiomyocytes (dediff.)" and "Proliferating cells") and Stereo-seq
border-zone domains, plus markers per day. They source `myh7l` as ventricular. Read: structure and
targeted rows. Reusable for direction A.

**Forward-citation finds (abstracts only).** A 2026 *Cells* paper (resection, scRNA-seq,
an AR-specific CM cluster; 1–2% of CMs cycling at 7 dpi), which gives an expected size for a
proliferating tail. An eLife 2026 paper on *phlda2*+ primordial CMs, a CM subpopulation that could
appear as a "state".

**This dataset's own labels.** *ttn.2* is named in the source paper (Hu et al. 2022), so any *ttn.2* result is not found unaided. It is circular and kept only as a comparison.

## 5. Our experiments and results

**Story so far.** C001 first made the data safe to use: a tested sanitiser removed the authors' labels.
It then read the main marker sources and checked every candidate gene against the data. The main
finding is that the literature does not give four separable marker sets. Border-zone markers are
the best supported. The other three states are thin, overlap border-zone, mismatch our timepoints,
or are barely detected. Hr11 is the lowest-depth library among CMs, though less extreme than the
all-cell numbers suggested. A lead for the next step: two genes are detected *more* often in Hr11
than in a control library, despite its lower depth.

C002 scored all four draft axes with both methods, against matched random gene sets and against
depth, within each library, under rules fixed before computing. Its answer is mostly negative. No
axis stands out from random genes of the same expression on its high tail. The two methods agree so
closely that agreement proves little. The rank method is more depth-bound than score_genes, not
less, though thinning to equal depth shows that this tie is technical and removable. Border-zone is
the only axis whose genes nearly co-vary as a program. The Hr11 *nppa* lead survived depth matching
and now points at ambient RNA as the next test. *ttn.2* behaves the same way; *ttn.2* is named in the source paper (Hu et al. 2022), so any *ttn.2* result is not found unaided.

| Cycle | Tried | Key result (evidence) | Signal | Does not show |
|---|---|---|---|---|
| C001 Understand | Sanitiser; literature; per-library QC and marker detection | Labels removed, test passing (R-E1). Shared genes across 3 states (R-E3, YAMLs). Proliferation and injury genes detected in ≤2% of CMs, except *ccng1* (R-E8, `results/c001_marker_detection.csv`). Hr11 has the lowest median depth per CM, but four libraries are close (`results/c001_library_qc.csv`). | Moderate for a border-zone axis; weak for the other three as separate axes | No scores, nulls or confound tests yet. Detection is not expression. Supplements not read. |
| C002 Probe | score_genes + UCell for 4 axes; 100 matched null sets; within-library depth; per-gene depth flags; binomial thinning; Hr11 depth-matched | P90 above null in 0–3 of 22 libraries, 0–2 of 18 hearts. Coherence above null: border-zone 17, injury-responsive 13, dedifferentiated 9, proliferating 4 of 22 libraries. Methods agree in 22/22. UCell depth tie larger than score_genes in 22/22 libraries; thinning removes most of it. 10 of 26 genes flagged. Hr11 *nppa* excess persists after depth matching, except in the top-depth bin (R-E11–R-E13, `results/c002_*`) | Negative for the axes as defined. Weak-positive for border-zone coherence. Strong against rank-score depth independence. A strong lead for the ambient trap | One rank method (UCell, re-implemented); one thinning seed; P90 may miss rare states; ambient and doublets untested |
| C003 Probe | Positive control by time per heart (thinned score_genes, tail vs 100-set null); Hr11 ambient test vs Hr12, Hr24, Hr27 and Ctrl with CM, atrial and non-CM panels | No axis passes on the tail (AUC 0.44–0.57). Proliferating passes on the median (AUC 0.82) via Hr6, Hr8 and Hr9. Non-CM genes ≥2× in 10/12 for Hr11's shallow cells, and also raised in the other 3 dpi hearts. *nppa* ≥2× in 0/3 other 3 dpi hearts; *myh6* not raised; *nppa* not correlated with contamination per cell (exploratory). Verdict inconclusive (R-E14, R-E15, `results/c003_*`) | Negative for the axes as defined. Strong for 3 dpi contamination as a trap. Weak for ambient as the cause of *nppa* | No validated positive control; chamber not controlled; one seed; no doublet scores |
| C004 Probe | Chamber call validated on the six chamber-only libraries; ventricle-only positive control; G4 with disjoint sets; agent label file | Rule B balanced accuracy 0.991 at coverage 0.735; rule A fails (0.867, 0.358). Ventricle-only tail AUC: injury-responsive 0.812 and dedifferentiated 7 dpi 0.857 pass; border-zone 0.786 and proliferating 0.679 fail. G4: all 6 pairs mixed. Labels: 30,064 unassigned; state counts near chance (R-E16, `results/c004_*`) | Chamber call strong. Axes still weak: the passes rest on control hearts sitting below the null. G4 unresolvable with these sets | `myh7l` unsourced; rule selected on its validation set; the null is uncalibrated per heart; no doublets |
| C005 fix + Probe | Ventricle-only label file; scrublet doublets on the GEO matrices | Labels: border-zone 67, dedifferentiated 53, injury-responsive 145, proliferating 183 against a chance level of about 140 (R-E17). Doublets: scrublet ran in 22/22 libraries (mapping 1.0), but its automatic threshold called 0 doublets everywhere. Median \|rho\| doublet score vs axis is 0.05–0.09 (under 0.1); atrial libraries are positive and Hr11/Hr27 negative. No state is enriched among the top 10% of doublet scores (exploratory) (R-E18) | The label file honestly says "no separable states". Doublets do not explain the scores | Auto-threshold failure; one method; no recalibrated doublet call |
| C007 Probe (direction A) | Li 2025 SD2 sets (25 genes, disjoint, CM-internal filter) against the C004 sets on the same ventricle-only pipeline, with a control-calibrated tail, a null-AUC test and leave-Hr9-out | 0 of 6 contrasts pass for either family. Li lower on 5 of 6. Null AUC p95 0.63–0.82. Hr9 and Hr8 high across axes. Exploratory no-CM-identity sets: no pass (R-E20, `results/c007_*`) | Clear negative for A as built. Points to a heart-level factor and to C or a different dataset | One run; 4 control hearts; one border-zone set rule |

**Supervisor check of C001** (SUP-E5, [cycles/C001/review.md](cycles/C001/review.md)):

- **Checked and confirmed:** the sanitised copy has no labels; the Hr11 depth figures; the Hr11
  lead. The lead is sharper than reported: *nppa*'s share of each cell's RNA is about 4× higher in
  Hr11 than in H7.
- **Correction:** the *vmhc* → `myh7` mapping is probably wrong for this data, because `myh7l` is
  the abundant one. To be fixed in C002.
- **New:** the mitochondrial contig `NC-002333.4` may be missing from the mito fraction.

## 6. Directions compared, and the current decision

**Decision after C001 (2026-11-07 15:30):** no direction is chosen yet. A direction decision needs
the reading gaps in FU-1 closed. Meanwhile, a Probe (C002) runs the core checks (two methods,
null, within-library depth) on all four candidate axes. It tells us which axes carry signal at
all, and it uses the same code as a border-zone-only probe. Doublets, ambient RNA and
co-occurrence follow this signal.

Candidate directions:

- **A. Literature markers, two-method scoring.** This is the user's and the slide's plan.
- **B. Spatially derived markers.** Border-zone and injury-area genes from Wu 2016 tomo-seq,
  Honkoop 2019 or Li 2025 domains. They are independent of this dataset, but the injury model and
  resolution differ.
- **C. Data-driven programs.** For example NMF or clusters, matched to literature afterwards. This
  risks circularity: markers found and tested on the same counts.


### Directions compared (researcher draft, C006, 2026-11-07 16:10; for the supervisor's decision)

What C002–C005 established for every direction:

- A usable axis must pass a positive control in ventricle-called cells with depth thinning, with a
  null calibrated in control hearts.
- It must use disjoint gene sets if co-occurrence (G4) is to mean anything.
- The current sourced sets fail these tests. Three things weigh against them: shared genes, depth
  and chamber ties, and contamination at 3 dpi.

New in C006 (R-E19): Li 2025 Supplementary Data 2 gives full Seurat marker tables for a scRNA-seq
"Cardiomyocytes (dediff.)" cluster (484 rows) and for "Proliferating cells" (477). It also lists
five Stereo-seq border-zone domains: activated 378, dedifferentiated 513, dediff.-prolif. 217,
proliferating 546 and re-differentiating 314 rows. Supplementary Data 5 gives up to 100 markers
per regeneration day. The same table sources `myh7l` as a ventricular CM marker, which closes the
gap in the chamber rule. It also shows `myh7` (*vmhc*) as a border-zone-domain marker
(log2FC 2.65–2.90 in dediff.-prolif. and prolif.), which may explain its rarity in our cells.

| Direction | Expected contribution | Main uncertainty | Strongest fair comparator | Cost |
|---|---|---|---|---|
| **A. Improved literature sets** (Li 2025 SD2/SD5 full lists; a sourced cell-cycle set; dedifferentiation including sarcomere-down genes; full reads of Honkoop, Wu, Ben-Yair, Sallin) | Larger, sourced, independent sets. Top-N genes per Li domain or cluster, made disjoint, can be scored with the C004 machinery. Li's scRNA-seq dediff. and prolif. clusters are single-cell, unlike the C001 spot lists | Model transfer: resection versus our cryoinjury, and Li's timepoints from 6 hpa to 28 dpa. Seurat "markers versus all cell types" include CM-versus-non-CM genes, which must be filtered to CM-internal contrasts. The current sets failed partly for biological reasons (rare states) that bigger lists will not fix | The current sourced sets (C002–C005), with the same positive control and nulls | Low. The tables are downloaded (R-E19); parsing, filtering and one Probe run are about one cycle. Reading the full texts is about one more |
| **B. Markers from raw spatial data** (Wu 2016 tomo-seq GSE74652, cryoinjury 3 and 7 dpi; Li 2025 Stereo-seq, CNGB STOmics STT0000071) | Border-zone genes defined by position, the field's own definition. Wu is **cryoinjury**, matching our model. Re-deriving from raw data lets us choose CM-internal contrasts and control depth | Tomo-seq sections mix cell types, so border-zone genes may be non-CM (contamination again). Stereo-seq spots are multicellular. Download and processing cost; CNGB access untested | Li SD2 domain lists (A), which already summarise B for Li; B adds Wu's cryoinjury | Medium to high. Wu's GEO tables may be small and quick; Li raw Stereo-seq is large. One cycle for Wu alone |
| **C. Data-driven programs in ventricular CMs with count splitting** (Neufeld et al. 2024, doi:10.1093/biostatistics/kxac047; found on the methods sheet, not read in full) | Finds whatever structure the data hold (programs, NMF or clusters) without relying on transferred markers. Count splitting gives valid tests of the programs on independent counts. The literature (A, B) is then used only to *name* programs | Programs may mostly capture depth, chamber, contamination and stress (all shown to be strong here). Naming a program from literature risks the circularity the study avoids. Count splitting assumes Poisson-like counts | A plus B sets scored on the same cells: do the data-driven programs beat sourced sets on the positive control? | Medium: NMF plus count splitting in Python, about one cycle, plus the controls already built |

Researcher's recommendation, for the supervisor to judge: **A first, as a Probe, using Li 2025 SD2.**

- Filter to CM-internal markers (border-zone domains against Vm, and dediff. against CM V).
- Use the top 25 genes per state, disjoint.
- Run the C004 ventricle-only positive control with a calibrated null.
- Read Wu 2016 in full in parallel to add the cryoinjury check (B-light).
- C is the fallback if A fails: it answers "what states do these cells hold at all" but cannot
  by itself say which literature state they are.
- What would change this: if A's sets also fail the positive control, the problem is likely the
  data (rare states, shallow CMs), not the markers. C, or a different dataset with more
  border-zone CMs (for example sorted *nppa*+ cells, GSE139218), then becomes the better path.

## 7. Search trail

Supervisor intake, 2026-11-07:

- The DOI search found Li et al. 2025, *Nat Commun* 16, doi:10.1038/s41467-025-59070-0, PMC12009352.
- nature.com redirected to a login page, so the paper was read from PMC instead, characters
  0–100,000 only. The Data Availability section and the full supplementary list were not read.

Researcher C001, 2026-11-07:

- Zotero list (one file; 384 unique titles), filtered to 82 zebrafish heart-regeneration papers:
  `sources/zotero_zf_heart_regen_states.csv`.
- Opened: eLife Honkoop 2019 (first 100k chars), Europe PMC XML of Li 2025 (first 100k chars),
  and ZFIN for *myh7*/*vmhc*. For Wu 2016, one web search, which returned ZFIN figure pages.
- Not opened: Li 2025 Data Availability and Supplementary Data 2 and 5; the second halves of both
  full texts; GSE74652 and GSE139218.
- No forward-citation search yet (OpenAlex). Not covered: Beisaw 2020 (AP-1, cryoinjury), other
  single-cell regeneration atlases, cell-cycle gene sets mapped to zebrafish.
- Stopped because of the 15:30 deadline, not because searches stopped returning new work.

Supervisor search, 2026-11-07, about 15:33 by the clock (first written as 15:50 in error), preparing the method review that a direction decision will
need:

- **Zebrafish CM states, 2023–2025, cryoinjury single-cell.** Two searches returned only Honkoop
  2019 and Li 2025, already in the map. Nothing new surfaced. Forman-Rubinsky 2025 (snRNA-seq) is
  known from its abstract only.
- **Scoring versus depth.** JASMINE (Noureen 2022, eLife 71994, already cited) attributes scoring
  bias to imbalance in the number of detected genes, with single-cell methods less affected than
  ssGSEA. A Genome Research paper (vol. 36, issue 3, p. 630; 2026) reports bias from control-gene
  binning in `score_genes`-type methods and limited score comparability for Seurat, scanpy, UCell
  and JASMINE. Its full text was behind a login, so it was not read. A Bioconductor support thread
  suggests downsampling (`scuttle::downsampleMatrix`) to test depth artefacts; this is a forum
  suggestion and matches the team's binomial thinning idea.
- **Ambient RNA in heart.** No source found that names *nppa* or *myl7* as ambient markers. The
  general finding is that contamination appears as highly expressed cell-type genes in other
  cells, and heart snRNA-seq contamination is mostly cardiomyocyte cytoplasm. Abstract-level only.
- **Not yet searched:** forward citations of Honkoop 2019 and Li 2025 (part of FU-1).

Researcher C007, 2026-11-07: no new search. Li SD2 (MOESM5) parsed in full for the set rule (R-E20).

Researcher C006, 2026-11-07 (R-E19):

- Forward citations via OpenAlex.
  - Honkoop 2019 (W2995700709, 292 citing works) filtered to "zebrafish cardiomyocyte": 68
    works. New and relevant: a 2026 *Cells* paper (10.3390/cells15161455; resection, scRNA-seq,
    an AR-specific CM cluster; 1–2% of CMs cycling at 7 dpi by flow cytometry; abstract only), and
    an eLife 2026 paper on *phlda2*+ primordial CMs (10.7554/elife.110256; abstract only).
  - Li 2025 (W4409591603, 36 citing works): mostly reviews; no new cryoinjury CM-state dataset.
  - Raw lists are in `sources/openalex/`.
- Li 2025 supplements downloaded from Springer: `sources/li2025_supp/` (MOESM3 description,
  MOESM4–6, 8–10). SD2 = MOESM5 and SD5 = MOESM8, mapped from the order in the MOESM3 legend
  (inferred). Read with `scripts/xlsx_peek.py`: sheet structure, cluster names and targeted rows,
  not every row.
- Li 2025 data: CNGB STOmics STT0000071 (web search; not opened).
- Still not opened: Wu 2016 full text and GSE74652; the second half of Honkoop 2019; Ben-Yair
  2019 and Sallin 2015 full texts; Beisaw 2020; Neufeld 2024 (count splitting) full text; MOESM7
  (SD4, a pdf).

Researcher C004, 2026-11-07: three web searches for chamber-gene sources found Yelon 1999 (vmhc)
and Berdougo 2003 (amhc), both via catalogue and ZFIN records. No source was found for `myh7l`
(vmhcl) ventricular expression.

Researcher C003, 2026-11-07: no new literature search; data-only Probe.

Researcher C002, 2026-11-07:

- ZFIN `myh7l` (ZDB-GENE-061027-393) via a search summary; ZFIN `myh7` opened in C001. NCBI
  E-utilities for NC_002333 (the nuccore web page was behind a browser check).
- Teammate branches not read (the user said optional).

## 8. Current state

- Cycle C007, status: **Continue**. C006 (Understand) was reviewed:
  - R5.1 inventory done ([checks-inventory.md](checks-inventory.md)).
  - FU-1 partly closed.
  - Chamber rule source-backed (SUP-E10).
- Latest review: [cycles/C006/review.md](cycles/C006/review.md).
- Latest assignment: [cycles/C007/assignment.md](cycles/C007/assignment.md) (Probe of direction
  A: sets from Li SD2 against the C004 sets, ventricle-only, leaving out Hr9).
- Answered so far:
  - Rank scores are more tied to depth (SUP-E6).
  - Published border and injury sets mark atrium (SUP-E7).
  - Hr11 carries a technical signature (C003).
  - With confounds removed, the original sets carry at most a weak injury signal (C004, C005).
  - Co-occurrence cannot be resolved with these sets (C004).
  - Agent state labels are at or below chance; chamber labels are validated and source-backed.
- For the human-run agreement check: `results/c005_agent_labels.csv` with
  [cycles/C005/note-for-evaluate.md](cycles/C005/note-for-evaluate.md).
- Not done: the ambient profile from raw droplets; a calibrated doublet call (FU-7); the rest of
  FU-1; a direction decision and its method review; Deepen.
- Scope changes: none. User-facing explainer: [docs/story.html](../../docs/story.html).
