# Tasks

Status: `[ ]` open · `[~]` in progress (name) · `[x]` done (link to PR or result).
Each task ends in a check that a script can run, so the loop can use it as a `/goal`.

## Who does what

| Role | Person | Lanes | Owns the question |
|------|--------|-------|-------------------|
| Scientist, team lead | Eme | `markers`, review of all biology | Which markers, which traps, which expert checks |
| Data and confounds | TBD | `data`, `confound` | Is a score biology or a technical artefact? |
| Scoring and co-occurrence | TBD | `score` | Do the four axes coincide or come apart? |
| Agent loop and evaluation | TBD | `loop`, `eval`, `docs` | Does the agent run the checks itself, and does it agree with the expert? |

## Phase 1, until 14:30: everyone can load the data

- [x] **data-1** (from sanitised copy, `src/zf/io.py`; [results/c001_data_facts.json](results/c001_data_facts.json)) `scripts/get_data.sh` and `src/zf/io.py` load the counts into one AnnData. Check: 32,044 cells, 18,794 genes, 22 libraries.
- [x] **data-2** ([results/c001_library_qc.csv](results/c001_library_qc.csv); 0 missing) Per-cell technical covariates: total counts, genes detected, mitochondrial fraction, library, timepoint. Check: no missing values.
- [~] **markers-1** (researcher C001: candidate YAMLs, presence in [results/c001_marker_detection.csv](results/c001_marker_detection.csv); `slc2a1a` absent) First marker sets for the four states, one YAML per state, each gene with a citation. Check: every gene exists in the data.
- [ ] **loop-1** `PROMPT.md`, `JOURNAL.md` and the `/goal` text; Claude Code can run one step end to end on a 2,000-cell subsample.
- [ ] **eval-1** `holdout/` set up; `evaluate.py` computes agreement (ARI, per-state confusion table) from a file of agent labels. Check: random labels give ARI near 0.

## Phase 2, until 15:30: scores and traps

- [x] **score-1** ([results/c002_axis_summary.csv](results/c002_axis_summary.csv); agree in 22/22 libraries) Score every cardiomyocyte on each axis with two methods (scanpy `score_genes`, AUCell or UCell). Check: the two methods correlate per axis.
- [x] **score-2** ([results/c002_axis_summary.csv](results/c002_axis_summary.csv), 100 matched sets) Null for each axis: random gene sets matched on mean expression. Check: report per axis how far the real score is above the null.
- [~] **confound-1** (within-library rho in [results/c002_axis_by_library.csv](results/c002_axis_by_library.csv); the `confounds_libsize.csv` file name not yet written) Correlation of each axis score with library size and genes detected, within each library. Check: table in `results/confounds_libsize.csv`.
- [~] **confound-2** (scrublet run: [results/c005_doublets_by_library.csv](results/c005_doublets_by_library.csv); auto-threshold called 0 doublets; `confounds_doublets.csv` name not written) Doublet scores (scrublet or scDblFinder); correlation of each axis with doublet score. Check: table in `results/confounds_doublets.csv`.
- [ ] **confound-3** Ambient RNA proxy: non-cardiomyocyte marker expression (for example fibroblast genes) inside cardiomyocytes, per library. Check: table per library.
- [ ] **confound-4** Library-specific clusters: is any cluster made mostly of one library (Hr11, 3 dpi, looks like one)? Check: per-cluster library entropy.
- [ ] **markers-2** Eme flags the known traps, for example the "dedifferentiation" marker that tracks library size. The agent must find these without being told; they go in `holdout/traps.yaml`.

## Phase 3, until 16:15: the loop runs on its own

- [~] **score-3** (C004: [results/c004_g4.json](results/c004_g4.json), all pairs mixed; no figure yet) Co-occurrence: pairwise correlation and joint high-score fractions per axis pair, per timepoint. Check: `results/cooccurrence.csv` and one figure.
- [ ] **loop-2** The loop reads literature (PubMed), proposes markers, scores cells and runs every confound check. It drops or flags markers that fail and logs each decision in `JOURNAL.md`.
- [ ] **eval-2** Agreement with the held-out annotation, and which traps the agent caught on its own. Check: `results/agreement.md`.

## Phase 4, until 17:00: submit

- [ ] **docs-1** `README.md`: problem, how to run, results with links, where it got stuck.
- [ ] **docs-2** Make the repository public; check licenses; fill the submission issue.
- [ ] **docs-3** 3-minute demo: team name, the question, one trap the agent caught, one it missed, what a domain expert adds.
