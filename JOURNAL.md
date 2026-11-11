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

## 2026-11-07 16:00 — Labeling protocol v1, pre-registered before any axis call is looked at

Owner directive 15:43: the deliverable is a fixed labeling protocol, graded by stability. Design by the senior reviewer (15:39) plus the owner's two statistical upgrades. Code `src/zf/protocol.py` (`label_axis(adata, genes)`), runner `scripts/run_protocol.py`, rules `docs/labeling-protocol.md`. Constants are fixed in the module and are not tuned on results.

Differences from the trust card above, and why:
- Decoys: 500 per axis (not 1,000), matched on detection rate within each library (20 bins), the reviewer's design; per-library matching makes the null follow each library's depth.
- Thinning target T = 600 UMIs; cells below it get `low_depth` and the rank call only.
- Two methods: UCell on raw ranks (max rank 1,500) and mean log-normalised expression on thinned counts; positive only if both agree.
- Call: BH within heart, positive q ≤ 0.05, negative p ≥ 0.10 in both methods, the rest uncertain.
- Vetoes to uncertain: doublet_flag, ≥2 non-CM (ambient) genes, hemoglobin fraction above the library 95th percentile, pct_mito > 30.
- Gate: per heart positive fraction vs decoy false-positive rate (20 decoys re-called as pseudo-axes; one-sided binomial, p < 0.05); timepoint supported when ≥ floor(n/2)+1 hearts pass; high-confidence needs the gate to survive leave-one-heart-out.
- Eligibility: within-heart partial Spearman |rho| ≥ 0.10 of the thinned score with log UMI or hemoglobin fraction ⇒ axis `depth_confounded`, no positives.
- Stability metrics: thinning kappa (thinned vs unthinned calls), marker bootstrap (80% of genes, positive kept in ≥80% of resamples), leave-one-heart-out. Sensitivity: gate without Hr8 and Hr11 (morphine study).
- Not done in v1 (time): meta-analysis, partial conjunction by Simes, the dip test, the taxonomy odds-ratio test, the multi-gene chamber score. Chamber strata use Card A's myosin call (57–87% accurate) for whole-heart libraries.

## 2026-11-07 16:07 — Protocol v1 amendment after run 1 (mechanical defects, not tuning)

Run 1 called 0 positives on every axis in every stratum. Two causes, both mechanical:
- The empirical p (1 + #decoys ≥ score)/(1 + 500) cannot go below 1/501. BH within a heart of N cells then needs at least 4% of the heart's cells at that floor, so any rarer state is uncallable. Fix: scDRS-style pooled p (z against the cell's own decoys, ranked against 50 decoy columns' z's pooled over the library's cells).
- The "≥2 non-CM genes" veto counted any of the 106 soup genes from `zf.confounds` and vetoed 18,947 of 32,044 cells (59%). Fix: the 4 fibroblast control genes of `config/axes.yaml`.
No threshold, T, decoy count, gate or ρ cut changed. Run 1's tables were overwritten by run 2 (same paths).

## 2026-11-07 16:12 — Dissociation control added (domain input 15:09); bootstrap cut to 5

Added before run 3 results: dissociation/heat-shock control axis through the same protocol, per-heart level (`results/hk/dissociation_per_heart.csv`), injury × dissociation partial ρ and injury verdict without the top-10% dissociation cells (`results/hk/injury_vs_dissociation.csv`), ventricle − atrium contrast per heart (`results/hk/chamber_contrast.csv`). Run 2 was stopped (about 1.7 min per axis with 20 bootstraps); run 3 uses 5 bootstrap resamples to fit the 16:25 cut. Decisions in `docs/decisions.md`.
