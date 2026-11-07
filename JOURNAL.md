# Journal

Append-only. One entry per result, decision or exchange with the team, with the time and links to outputs.

## 2026-11-07 15:35 — Card A: label-free data, technical flags, marker sets v1 (hk/data-markers)

- Label-free copy: `src/zf/io.py` `build()` keeps an allowlist of technical metadata, writes the authors' published per-cell annotation and their clustering to the held-out folder, drops their embeddings, and writes `data/cm_labelfree.h5ad` (32,044 cells x 18,794 genes, 22 libraries, 18 hearts: 4 uninjured, 4 at 3 dpi, 7 at 7 dpi, 3 at 30 dpi). Covariates: total_counts, log_umi, n_genes, pct_mito (recomputed from raw counts, matches the source column within 1e-6), library, heart, condition, chamber (GEO titles only; whole-heart libraries = unknown), chemistry, study, morphine_study. Timepoint and chamber are cross-checked between GEO records (`config/dataset.yaml`) and the owner's file; both agree for every cell.
- Per-cell chamber call (myh7 + myh7l vs myh6, sign of the difference) on the six chamber-dissected libraries: accuracy H8a 0.80, H8v 0.87, Hr6a 0.65, Hr6v 0.72, Hr7a 0.57, Hr7v 0.75. Too weak to use as a covariate; recorded as `va_call` only.
- Flags on the full GEO libraries (`src/zf/confounds.py`, all cell types, mapped to cardiomyocyte barcodes; every barcode found): per-library table `results/hk/library_flags.tsv`.
  - Doublets (scanpy pp.scrublet): 99 cardiomyocytes flagged (0.3%). scrublet's automatic threshold collapsed on Hr11 (0.018, calling 96% of the library); Hr11 now uses the median threshold of the other libraries (0.54). Observation: the shallowest library breaks the default doublet caller.
  - Ambient RNA: 106 non-cardiomyocyte-specific genes (col1a1a/b, col1a2, fn1b, sparc, cfd, cd74a/b, mmp2, krt genes). Median contamination per library 0 to 0.18; Hr11 is highest (0.18), Hr24 next (0.10), both at 3 dpi. Observation: ambient load is confounded with 3 dpi. Caveat: at ~700 UMIs one soup-gene count gives a fraction of 0.1 to 0.3, so `ambient_frac` is noisy per cell and `ambient_flag` (above 10% and 3 MADs) marks 11% of cells; use `ambient_frac` as a continuous covariate, not the flag as a filter. No cell is dropped.
- Marker sets v1 (`results/hk/markers_v1/`, README there; claims in `markers/sources/hk_corpus_claims.tsv`; corpus `literature/hk-corpus.md`): 37 core genes, disjoint across axes, each in var_names with a DOI. Injury axis has no immediate-early genes; nppa sits in `ambiguous.yaml`; the published depth-suspect titin genes are `library_size_suspect`, out of every core. Observation (r2 helper): the top-50 genes of the authors' dedifferentiation cluster include dissociation-stress genes (hsp90aa1.1, unc45b, cited4b, her4.1).
- Audio transcript: done by the r5 helper (`.octoally-pro/scratch/research/r5_transcript.md`, project folder only).

## Pre-registered trust card

Written 2026-11-07 15:35 Berlin, before any axis score exists. Not edited after scores exist; changes go in a new dated entry that says why.

Unit of replication: the heart (18). Chemistry (10x v3, Hr24/Hr27) exists only at 3 dpi and is never a separate covariate; library, heart and timepoint are never regressed out.

Per axis, the gating checks are those in `config/checks.yaml` (null_matched_sets, direction, depth, doublet, ambient, dissociation, batch, leave_one_gene_out, method_agreement); outcome `pass` only, `within_noise` and `unstable` do not pass. New parts this card adds to them:

1. **Detection-matched decoys.** Each axis is compared with 1,000 random gene sets matched on mean expression and on detection rate. Per-cell empirical p = (1 + decoys scoring at least as high) / 1,001; Benjamini–Hochberg FDR within each heart.
2. **Marker bootstrap stability.** Resample the core genes with replacement 200 times; a cell is stable if it stays above the FDR cut in at least 80% of resamples.
3. **Clean flags.** A high-confidence cell has `doublet_flag` false and `ambient_frac` below 0.1.
4. **Depth thinning.** Binomially thin every cell to the 3 dpi depth floor (Hr11 median) and rescore; a high-confidence cell must stay above the cut after thinning. The 3 dpi libraries (median UMIs per cardiomyocyte: Hr11 659, Hr27 1,039, Hr24 1,158, Hr12 1,166) are the natural depth-trap test.
5. **Per-heart replication.** Per axis and condition, a one-sided test per heart, combined by partial conjunction (Benjamini–Heller 2008) with Simes combination, threshold u = floor(n/2) + 1 hearts of the n in that condition (3 of 4 uninjured, 3 of 4 at 3 dpi, 4 of 7 at 7 dpi, 2 of 3 at 30 dpi).
6. **Leave one heart out.** The axis-level verdict must not change when any single heart is removed.
7. **Meta-analysis.** Random-effects REML across hearts plus an equal-weight analysis (REML with tau2 = 0 lets Hr11 dominate). tau2 and I2 are reported, never gated on (3 to 4 hearts per condition).

A cell is **high-confidence** for an axis only if 1 to 4 all hold; an axis is **trusted** only if every gating check passes and 5 and 6 hold. Everything else is reported as an observation, not a claim. IDR is not used.

## 2026-11-07 15:55 — Review of card A (reviewer, before commit)

- Chamber-call accuracies above now come from a file: `results/hk/va_accuracy.tsv`, written by `uv run python -m zf.io` (manifest `results/runs/hk-io-build/manifest.yaml`).
- `results/hk/library_flags.tsv` had lost its header row (hand copy). `zf.confounds` now writes it itself, with manifest `results/runs/hk-confounds/manifest.yaml`. The rerun reproduced `data/derived/cm_flags.tsv` byte for byte.
- `zf.confounds` could not run from the locked environment: scrublet's automatic threshold needs scikit-image, which was missing. Added with `uv add scikit-image`.
- A missing automatic scrublet threshold (NaN) used to count as "sane" and call no doublets. It now falls back to the median threshold like Hr11.
- Observation, not fixed: scrublet calls far fewer doublets than the 10x expected rate. In `results/hk/library_flags.tsv`, 8 of 22 libraries call none, and the median ratio `doublet_rate_full / doublet_expected_rate` is far below 1. So `doublet_flag` false is a weak clean-cell test (trust card item 3). Use `doublet_score` as a continuous covariate too, and treat the doublet check as low-power until scDblFinder or a lower threshold is compared.
