# C002 reading rules (researcher), fixed 2026-11-07 15:24 CET before any score was computed

All statistics are within library (22 libraries); summaries are counts of libraries.
Scores: scanpy score_genes on log1p(CP10k) of the sanitised counts (default ctrl_size=50,
n_bins=25, fixed seed); UCell (Andreatta 2021 formula, maxRank 1500) on the raw counts.
Null: 100 random gene sets per axis, each gene drawn from the same mean-expression bin
(25 quantile bins of mean log1p CP10k over all cells) as the axis gene it replaces; ttn.2 excluded.

1. **Above null (per library):** the axis's 90th-percentile score in the library exceeds the
   95th percentile of the same statistic over the 100 null sets. **Coherence (per library):** the
   mean pairwise Spearman correlation among the axis genes (log1p CP10k) exceeds the 95th percentile
   of the null sets. An axis is "above null" under a method if it passes in at least 18 of 22
   libraries; "mixed" if 8-17; "not above null" if fewer than 8.
   Caveat fixed in advance: P90 above null can come from genes that are simply high in all CMs;
   only coherence speaks to a shared program.
2. **Tracks depth:** within-library Spearman rho of the score with log10 total counts and with
   genes detected. An axis tracks depth if the median |rho| over libraries is >= 0.3; "weak" if
   0.1-0.3; otherwise "not". Also reported relative to the null sets' median |rho|.
3. **Methods agree:** within-library Spearman between score_genes and UCell >= 0.5 in at least 18
   of 22 libraries.
4. **R3.6, rank less tied to depth:** UCell's median |rho with log counts| is lower than
   score_genes' for at least 3 of 4 axes and also for the null sets. Hr11 vs Hr7a, Hr2a, Hr2b is
   descriptive only (timepoint and chamber differ).
5. **Per-gene depth flag:** within-library Spearman of the gene's log1p CP10k with log10 total
   counts, median over libraries. Flag the gene "tracks_depth" if this median lies outside the
   5th-95th percentile of 2,000 random background genes in the same mean-expression bin.
   Flags are written into the YAMLs; no gene is dropped.
6. **C9 lead:** nppa and ttn.2 share of UMIs in Hr11 vs control cells in the same total-count
   bin (depth-matched). "Not explained by depth" if the Hr11 median share stays at least 2x the
   control median within matched bins.

Any change after seeing results is labelled exploratory.

## Addendum, 2026-11-07 15:30 CET, after the primary results, before the thinning check

- Heart counts (supervisor note): 22 libraries are 18 hearts. A heart passes only if all of its
  libraries pass. The library thresholds above are kept; heart counts are reported alongside.
- Exploratory null (labelled exploratory, run after the primary results): 100 bins, with mt-,
  NC-002333 and ribosomal genes removed from the null pool.
- **Thinning check (R3.6, stronger form):** binomially thin every cell's counts to 500 UMIs (cells
  below 500 are kept unthinned and counted), then rescore both methods. A method's depth tie is
  "technical" if its within-library |rho(score, log10 original depth)| falls by at least half after
  thinning. Rank is "less tied to depth" if, after thinning, its median |rho| is lower than
  score_genes', and if its rho between original and thinned scores is higher.
