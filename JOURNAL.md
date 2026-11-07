# Journal

Append-only. One entry per result, decision or exchange with the team, with the time and links to outputs.

## 2026-11-07 14:08 CET: input data moved from ~/Downloads into data/ (user request)

- `data/hu2022/`: `Hu2022_CM_Ctrl_CMonly_processed_HVG.h5ad`, `Hu2022_CM_Ctrl_CMonly_processed_allgenes.h5ad`.
- `data/zotero/`: `HR_1.csv`, `CM_zotero_papers.csv`, `Integration_zotero_papers.csv`, `singlecellgenomics_zotero_papers.csv`.
  All four have the same MD5 (`4d38e63f8e7d645f0bde8ba26c8ea926`), so they are one file under four names. Ask the user whether the intended exports are missing.
- `data/claudelabday/20261107_ClaudeLabDay.zip`: moved, not extracted. Its contents are unknown. List it before unpacking and do not extract anything that looks like held-out annotation.
- `data/` is gitignored; nothing committed.

## 2026-11-07 14:51 CET: study started; label columns found in the supplied data

- The Research Supervisor opened the study [research/cm-state-axes/](research/cm-state-axes/research-map.md).
- The supplied `.h5ad` files contain the authors' annotation columns. Their names are not listed in this public repository.
  The supervisor saw only column names and the number of unique values per column, never the label values.
  The agreement check is not void. Evidence: [evidence.md](research/cm-state-axes/evidence.md) SUP-E1 to SUP-E3.
- Decision: no one uses the data until a sanitised copy with an allowlist of columns exists (first step of cycle C001).
- First assignment: [cycles/C001/assignment.md](research/cm-state-axes/cycles/C001/assignment.md) (Understand: marker sources per state, data facts).

## 2026-11-07 15:20 CET: C001 (researcher): leakage guard, candidate markers, data facts

- Leakage guard: `src/zf/io.py::sanitize_h5ad`, test-first ([tests/test_io_sanitize.py](tests/test_io_sanitize.py), 3 passed).
  Sanitised copy `data/clean/hu2022_cm_allgenes_sanitized.h5ad` keeps only the obs allowlist, counts, var names and stats.
  Leakage log: no label value printed or seen; only key *names* of the source file were listed. `holdout/` not opened. The agreement check is not void.
- Data facts: [results/c001_library_qc.csv](results/c001_library_qc.csv), [results/c001_marker_detection.csv](results/c001_marker_detection.csv), provenance [results/c001_data_facts.json](results/c001_data_facts.json).
- Draft marker YAMLs (status `candidate`) in `markers/`.
- Report: [research/cm-state-axes/cycles/C001/report.md](research/cm-state-axes/cycles/C001/report.md), status Checkpoint.

## 2026-11-07 15:45 CET: pulled remote main (force-pushed); trap hint noted

- `origin/main` was rewritten (371e96d → fea5089). Local main was moved with `git reset --keep origin/main`. Uncommitted work is unchanged.
  The new main removes the annotation's column and label names from `docs/data.md` and the annotation download from `scripts/get_data.sh`.
- Trap-test caveat: *ttn.2* is named in the source paper (Hu et al. 2022), so any *ttn.2* result is not found unaided (evidence R-E6).
  Label values were never seen, so the agreement check is not void.
  Recorded as FU-6 in [follow-ups.md](research/cm-state-axes/follow-ups.md).
- Eleven unmerged teammate branches overlap this study (Li 2025 supplementary markers, heart as the unit, the depth-thinning check, a second sanitiser, colliding `markers/*.yaml`). See FU-5.

## 2026-11-07 15:40 CET: C002 (researcher): four axes vs null and depth

- Reading rules fixed before computing: [research/cm-state-axes/cycles/C002/reading-rules.md](research/cm-state-axes/cycles/C002/reading-rules.md).
- Code: `src/zf/score.py` (UCell + matched random sets, test-first: [tests/test_score.py](tests/test_score.py)), `scripts/c002_probe.py`, `scripts/c002_hearts_and_flags.py`, `scripts/c002_thinning.py`.
- Results: [results/c002_axis_summary.csv](results/c002_axis_summary.csv), [results/c002_heart_summary.csv](results/c002_heart_summary.csv), [results/c002_thinning_summary.csv](results/c002_thinning_summary.csv), [results/c002_gene_depth.csv](results/c002_gene_depth.csv), [results/c002_c9_depth_matched.csv](results/c002_c9_depth_matched.csv); exploratory null in `results/c002_exploratory_null/`.
- Depth flags written into `markers/*.yaml` (`c002_depth`); no gene dropped. `myh7` kept after the ZFIN check (F1).
- F2: `percent.mito` excludes the 19 `NC-002333.*` mitochondrial features.
- *ttn.2* is named in the source paper (Hu et al. 2022), so any *ttn.2* result is not found unaided (FU-6). `holdout/` not opened.
- Report: [research/cm-state-axes/cycles/C002/report.md](research/cm-state-axes/cycles/C002/report.md), status Review requested.

## 2026-11-07 15:37 CET: C002 reviewed by the supervisor; C003 assigned

- **Correction of timestamps:** the supervisor's entry above headed "15:45 CET" was written at about 15:32 by the system clock.
  The research map's "15:50" search note was written at about 15:33. The times were estimated, not read from the clock.
- Independent recheck ([checks/c002_verify.out](research/cm-state-axes/checks/c002_verify.out)) confirmed three claims:
  - The border axis is not above its matched null.
  - Rank scores are more tied to depth than `score_genes`. This extends to an AUCell-style score, so the user's rank-independence hypothesis is refuted for these data.
  - Hr11's *nppa* excess survives depth matching. The claim that it vanishes at high depth is overstated: it shrinks.
- Finding F5: C002 had no positive control, so "no axis beats random genes" cannot yet be read as "the markers are fake".
- Review: [research/cm-state-axes/cycles/C002/review.md](research/cm-state-axes/cycles/C002/review.md).
  Next: [cycles/C003/assignment.md](research/cm-state-axes/cycles/C003/assignment.md) (positive control per heart; ambient test of the Hr11 lead).

## 2026-11-07 15:48 CET: C003 reviewed; chamber confound found; C004 assigned

*(Timestamp correction, appended at 15:44 by the clock: this entry and the C003 review and C004 assignment headed 15:46–15:48 were written at about 15:42. The supervisor's estimated times ran a few minutes ahead of the system clock.)*

- Supervisor check ([checks/c003_verify.out](research/cm-state-axes/checks/c003_verify.out)):
  - Within the same heart, atrium-only libraries score 5–10× higher on the border-zone axis than ventricle-only libraries. The *nppa*-led sets track chamber.
  - The atrial fraction varies widely across hearts, which confounds C003's between-heart positive control.
  - In ventricle-called cells (exploratory, crude call), the injured hearts score highest.
- Review: [cycles/C003/review.md](research/cm-state-axes/cycles/C003/review.md).
  Next: [cycles/C004/assignment.md](research/cm-state-axes/cycles/C004/assignment.md) (validated chamber call, a ventricle-only positive control, G4 or "cannot answer", and the agent label file).

## 2026-11-07 15:47 CET: C003 (researcher): positive control and the Hr11 ambient test

- Rules fixed first: [research/cm-state-axes/cycles/C003/reading-rules.md](research/cm-state-axes/cycles/C003/reading-rules.md).
- Positive control per heart (18 hearts): no axis passes on the tail statistic. [results/c003_positive_control.csv](results/c003_positive_control.csv), [results/c003_heart_scores.csv](results/c003_heart_scores.csv).
- Hr11 ambient test: verdict inconclusive by the fixed rule. Non-CM RNA, *hsp70l* stress and CM-gene dilution mark Hr11 and, less strongly, the other 3 dpi hearts. [results/c003_ambient_by_bin.csv](results/c003_ambient_by_bin.csv), [results/c003_ambient_verdict.json](results/c003_ambient_verdict.json); exploratory [results/c003_exploratory_cellcorr.csv](results/c003_exploratory_cellcorr.csv).
- ttn.2 results are hinted (FU-6). `holdout/` was not opened.
- Report: [research/cm-state-axes/cycles/C003/report.md](research/cm-state-axes/cycles/C003/report.md), status Review requested.

## 2026-11-07 15:52 CET: C004 (researcher): chamber call, agent label file, ventricle-only control, G4

- Rules fixed first: [research/cm-state-axes/cycles/C004/reading-rules.md](research/cm-state-axes/cycles/C004/reading-rules.md).
- Agent label file for the human-run agreement check (R5.2): [results/c004_agent_labels.csv](results/c004_agent_labels.csv). Its format is in [results/c004_probe.json](results/c004_probe.json) (`label_file_format`). `holdout/` was never read; the rule was fixed before computing.
- Chamber validation: [results/c004_chamber_validation.csv](results/c004_chamber_validation.csv). Ventricle-only control: [results/c004_ventricle_positive_control.csv](results/c004_ventricle_positive_control.csv). G4: [results/c004_g4.json](results/c004_g4.json).
- Report: [research/cm-state-axes/cycles/C004/report.md](research/cm-state-axes/cycles/C004/report.md), status Review requested.

## 2026-11-07 15:49 CET (clock): C004 reviewed; label file needs a chamber fix; C005 assigned

- Supervisor check ([checks/c004_verify.out](research/cm-state-axes/checks/c004_verify.out)):
  - The label file is complete (32,044 unique cells).
  - Chamber call balanced accuracy 0.986 on the six chamber-only libraries (optimistic, because the rule was chosen there).
  - **F10:** three-quarters of the agent's `injury_responsive` calls are atrial cells, most frequent in uninjured and 30 dpi hearts. The labels inherit the chamber trap.
- G4 accepted as an honest negative: it cannot be resolved with the sourced sets.
- Review: [cycles/C004/review.md](research/cm-state-axes/cycles/C004/review.md).
  Next: [cycles/C005/assignment.md](research/cm-state-axes/cycles/C005/assignment.md) (ventricle-only label correction; doublets on the GEO matrices if time allows).


## 2026-11-07 15:56 CET: C005 (researcher): corrected label file (F10)

- **Use [results/c005_agent_labels.csv](results/c005_agent_labels.csv) for evaluate.py.** It supersedes c004 (kept). States are scored in ventricle-called cells only; other cells are `not_scored`. Format: [results/c005_agent_labels.json](results/c005_agent_labels.json). Note for the human runner: [research/cm-state-axes/cycles/C005/note-for-evaluate.md](research/cm-state-axes/cycles/C005/note-for-evaluate.md).
- Effect: injury-responsive calls fell from 1,188 to 145. All state classes are at or below chance ([results/c005_label_counts_vs_c004.csv](results/c005_label_counts_vs_c004.csv)).
- `uv add scikit-image` (for scanpy's scrublet threshold); pyproject.toml and uv.lock changed.

## 2026-11-07 15:59 CET: C005 (researcher): doublets (R3.2)

- Scrublet on the 22 GEO matrices (downloaded with `scripts/get_data.sh` to `data/raw/`): [results/c005_doublets_by_library.csv](results/c005_doublets_by_library.csv), [results/c005_doublets_summary.json](results/c005_doublets_summary.json).
- No axis tracks the doublet score (median |rho| ≤ 0.091). The automatic doublet call found 0 doublets, so the share of states among doublets is **not answered** (not "passed"). Exploratory: [results/c005_doublets_exploratory_top10.csv](results/c005_doublets_exploratory_top10.csv).
- Report: [research/cm-state-axes/cycles/C005/report.md](research/cm-state-axes/cycles/C005/report.md), status Review requested.

## 2026-11-07 16:01 CET (clock): C005 reviewed; C006 assigned

- Supervisor check ([checks/c005_verify.out](research/cm-state-axes/checks/c005_verify.out)):
  - The corrected labels [results/c005_agent_labels.csv](results/c005_agent_labels.csv) call states only in ventricle cells. F10 is closed.
  - The doublet threshold failed (about 0 called), so the categorical doublet check is **not answered** (FU-7). Continuous doublet scores show no axis tracking doublets.
- Use `results/c005_agent_labels.csv` for the human-run agreement check; see [note-for-evaluate.md](research/cm-state-axes/cycles/C005/note-for-evaluate.md).
- Review: [cycles/C005/review.md](research/cm-state-axes/cycles/C005/review.md).
  Next: [cycles/C006/assignment.md](research/cm-state-axes/cycles/C006/assignment.md) (an inventory of checks for the submission, then closing FU-1 to choose a direction).


## 2026-11-07 16:10 CET: C006 (researcher): checks inventory and direction draft

- R5.1 inventory, for the submission: [research/cm-state-axes/checks-inventory.md](research/cm-state-axes/checks-inventory.md). It lists who initiated each check, its result and the trap caught, plus the checks an expert would add. *ttn.2* is marked as hinted.
- Li 2025 Supplementary Data 2 and 5 downloaded to `research/cm-state-axes/sources/li2025_supp/`. They source `myh7l` as ventricular. Forward citations are in `research/cm-state-axes/sources/openalex/`.
- Direction comparison drafted in research-map section 6. Report: [research/cm-state-axes/cycles/C006/report.md](research/cm-state-axes/cycles/C006/report.md), status Checkpoint.

## 2026-11-07 16:08 CET (clock): C006 reviewed; C007 (Probe of direction A) assigned

- The inventory of checks (R5.1), for the submission: [research/cm-state-axes/checks-inventory.md](research/cm-state-axes/checks-inventory.md).
  It has 35 checks with who initiated each, plus 11 checks a domain expert would add.
- The supervisor confirmed in Li 2025 SD2 that `myh6` marks atrial and `myh7l` ventricular cardiomyocytes (SUP-E10). The chamber rule is now source-backed.
- Review: [cycles/C006/review.md](research/cm-state-axes/cycles/C006/review.md).
  Next: [cycles/C007/assignment.md](research/cm-state-axes/cycles/C007/assignment.md) (sets built from Li SD2 against the C004 sets on the ventricle-only positive control).


## 2026-11-07 16:12 CET: C007 (researcher): direction A (Li 2025 SD2 sets) Probe

- Rules fixed first: [research/cm-state-axes/cycles/C007/reading-rules.md](research/cm-state-axes/cycles/C007/reading-rules.md). Sets: [markers/li2025_sd2_sets.tsv](markers/li2025_sd2_sets.tsv).
- No axis passes for either the Li or the C004 sets: [results/c007_positive_control.csv](results/c007_positive_control.csv), per heart in [results/c007_heart_tails.csv](results/c007_heart_tails.csv). The C004 passes vanish under a calibrated threshold. Hr9 and Hr8 drive most tails.
- Report: [research/cm-state-axes/cycles/C007/report.md](research/cm-state-axes/cycles/C007/report.md), status Review requested.

## 2026-11-07 16:20 CET (clock): C007 reviewed; study stops at submission, not complete

- C007 (Probe of direction A, sets from Li SD2): 0 of 6 contrasts pass for either marker family ([results/c007_positive_control.csv](results/c007_positive_control.csv)).
  The supervisor confirmed that the cell-cycle signal sits in hearts Hr9 and Hr8 ([checks/c007_verify.out](research/cm-state-axes/checks/c007_verify.out)).
- Methodology and findings for submission: [docs/findings-report.md](docs/findings-report.md).
- Status: **not complete**. Restart point: direction C or a validated positive control, after FU-1 and FU-7. See [cycles/C007/review.md](research/cm-state-axes/cycles/C007/review.md).

## 2026-11-07 16:24 CET (clock): methodology and report split into two documents

- `docs/findings-report.md` (named in the 16:20 entry) is replaced by [docs/research-methodology.md](docs/research-methodology.md), covering the two members, their responsibilities and how a cycle runs, and [docs/research-report.md](docs/research-report.md), the progress so far against that methodology.

## 2026-11-07 16:33 CET (clock): review changes for PR #17 (Adish)

- Removed every naming of the source file's label columns from the repository (journal, evidence SUP-E1/SUP-E2, C001 assignment, sanitiser test stand-ins). The repository goes public.
- Replaced the old-doc *ttn.2* reference with: *ttn.2* is named in the source paper (Hu et al. 2022), so any *ttn.2* result is not found unaided (R-E6, FU-6).
- Moved this study's marker files from `markers/` to `proposals/loop/` for Emi's review, per the loop spec. Result sidecars from before the move still record the old `markers/` paths. The contents are unchanged.

## 2026-11-07 16:31 CET (clock): one-page research summary added

- [docs/research-summary.md](docs/research-summary.md): the question, the method in brief, six findings with linked numbers, answers to the problem's questions, and what is not established. The detailed log stays in [docs/research-report.md](docs/research-report.md).

## 2026-11-07 16:39 CET (clock): user question, what do the markers mark inside the ventricle? C008 assigned

- Supervisor check ([checks/c008_verify.out](research/cm-state-axes/checks/c008_verify.out)):
  - In 14,024 ventricle cells, the border, injury and dedifferentiated scores follow the *nppa* share (median rho 0.52–0.58) more than depth (0.22).
  - Proliferating tracks neither. UCell's larger depth tie holds (19/19 libraries).
- Gap found by the user's question (FU-9): the C002, C003 and C005 checks had been run on both chambers mixed.
  Next: [cycles/C008/assignment.md](research/cm-state-axes/cycles/C008/assignment.md).


## 2026-11-07 16:45 CET: C008 (researcher): checks re-run inside the ventricle (user correction FU-9)

- Rules: [research/cm-state-axes/cycles/C008/reading-rules.md](research/cm-state-axes/cycles/C008/reading-rules.md). Summary: [results/c008_summary.json](results/c008_summary.json).
- Inside the ventricle, the border, injury and dedifferentiated axes are one *nppa*/*nppb* signal. Without those genes they do not co-vary ([results/c008_program_by_library.csv](results/c008_program_by_library.csv)).
- The Hr11 *nppa* excess in ventricle cells fits the ambient rule; this changes the C003 conclusion ([results/c008_ambient_ventricle.csv](results/c008_ambient_ventricle.csv)).
- Report: [research/cm-state-axes/cycles/C008/report.md](research/cm-state-axes/cycles/C008/report.md), status Checkpoint.

## 2026-11-07 16:45 CET (clock): C008 reviewed; summary updated; study paused at submission, not complete

- C008 (ventricle re-checks): [report](research/cm-state-axes/cycles/C008/report.md), [review](research/cm-state-axes/cycles/C008/review.md), supervisor check [checks/c008b_verify.out](research/cm-state-axes/checks/c008b_verify.out).
- Inside the ventricle, the three *nppa*-led states are one program; without *nppa* and *nppb* they barely co-vary (confirmed).
- Hr11's *nppa* excess fits ambient RNA, probably atrial (partly confirmed; C003's reading is corrected).
- [docs/research-summary.md](docs/research-summary.md) updated.

