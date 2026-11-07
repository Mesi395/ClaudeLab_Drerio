# C002 report: Probe (researcher)

Written 2026-11-07 15:35 CET. Reading rules fixed at 15:24, before any score was computed:
[reading-rules.md](reading-rules.md). The addendum (15:30) covers heart counts, the exploratory
null and the thinning check. All numbers below come from the files in `results/` named in each row.

## Findings resolved first

- **F1 (vmhc).** By ZFIN, *vmhc* is a previous name of `myh7`. `myh7l` is a different paralog,
  formerly *vmhcl* (R-E9). The YAML keeps `myh7` with the check recorded. Its low detection here
  stays unexplained. **Partly closed:** the symbol is settled from a source, but the reason
  for the low detection is not. Ensembl IDs are not in `var`. The teammate TSV was not read
  (optional, per the user).
- **F2 (NC-002333.4).** NC_002333 is the *Danio rerio* mitochondrial genome (NCBI). The 19
  `NC-002333.*` features are excluded from the supplied `percent.mito`, which equals the `mt-`
  share exactly. They carry a median of 0.0628 of each cell's UMIs (R-E10). **Closed:** any mito
  covariate must add them.

## Claims

Library thresholds per the rules: at least 18 of 22 libraries counts as "above null" or "agree",
8–17 as "mixed". Heart counts are out of 18 hearts; a heart passes only if all its libraries pass
([c002_heart_summary.csv](../../../../results/c002_heart_summary.csv)).

| # | Claim | Evidence for | Against / limits |
|---|---|---|---|
| P1 | **No axis is above its matched random null on the 90th-percentile score**, with either method. Libraries passing (UCell / score_genes): border-zone 0/1, dedifferentiated 0/0, injury-responsive 0/3, proliferating 2/1. Hearts: 0–2 of 18. The axes do not mark a high subset of cells any better than random genes of the same expression. | [c002_axis_summary.csv](../../../../results/c002_axis_summary.csv), heart summary | The exploratory null (100 bins, no mito/NC/ribosomal genes in the pool) gives the same reading: 0–8 of 22 libraries ([exploratory](../../../../results/c002_exploratory_null/c002_axis_summary.csv)). P90 is one statistic; a rare state (a few % of cells) may not move P90. |
| P2 | **Coherence (the genes co-vary beyond random sets) is mixed, and weakest for proliferating.** Libraries above null: border-zone 17, injury-responsive 13, dedifferentiated 9, proliferating 4 of 22. Hearts: 13, 10, 6, 2 of 18. Border-zone is closest to a real program, one library short of the rule. | axis summary; heart summary | Coherence values are small in absolute terms (for example border-zone in H7 0.086 against a null 95th percentile of 0.083; [by library](../../../../results/c002_axis_by_library.csv)). The exploratory null raises injury-responsive to 16 and proliferating to 8. |
| P3 | **The two methods agree almost perfectly**: within-library Spearman ≥ 0.5 in 22/22 libraries and 18/18 hearts for every axis, median 0.91–0.96. Agreement therefore cannot serve as an independent check. | axis summary | Agreement for the random sets was not computed. Both methods read the same few detected genes. |
| P4 | **The rank score (UCell) is *more* tied to depth than score_genes, not less.** Median within-library \|rho\| with log10 counts, UCell vs score_genes: border-zone 0.468 vs 0.287; dedifferentiated 0.438 vs 0.272; injury-responsive 0.364 vs 0.185; proliferating 0.276 vs 0.073. UCell's tie is larger in 22/22 libraries for every axis. By rule 2, UCell "tracks depth" on three axes; score_genes is "weak" on three and "not" on proliferating. | axis summary; by-library | Matched random sets are just as tied (UCell null median 0.29–0.51): the tie is generic to scoring sparse CMs, not specific to these markers. |
| P5 | **Thinning to equal depth removes most of UCell's depth tie, so that tie is technical.** After binomial thinning to 500 UMIs, the UCell share of the tie that remains is 0.13–0.30; for score_genes it is 0.27–0.54, and 1.02 for proliferating, which had almost no tie to begin with. After thinning, UCell is less tied than score_genes on 2 of 4 axes. Its scores correlate less with the unthinned scores on 3 of 4 axes (UCell 0.748, 0.728, 0.796, 0.659 against score_genes 0.760, 0.755, 0.804, 0.656, in the order of the table). By rule 4, the rank score is **not** less tied to depth. The user's hypothesis is refuted for these data and this rank method. | [c002_thinning_summary.csv](../../../../results/c002_thinning_summary.csv) | One seed, one target depth. Thinning also removes signal. AUCell was not run, only UCell. |
| P6 | **Per-gene depth flags**: 10 of 26 axis genes (counting genes shared between axes once) lie outside their expression bin's 5th–95th percentile of 2,000 background genes. Above it ("up with depth"): *desma, mef2cb, nkx2.5, pkma, pdk2a, ldha, hspb11, ccng2*. Below it: *myh7, gata4*. Flags are written into the YAMLs; no gene was dropped. | [c002_gene_depth.csv](../../../../results/c002_gene_depth.csv); `markers/*.yaml` (`c002_depth`) | Fragile: top-bin axis genes sit at 0.24–0.28 against a bin 95th percentile of 0.257 (for example *nppb* 0.240 is not flagged, *pkma* 0.276 is). Detection effects in sparse genes inflate these correlations. |
| P7 | **C9 (F3): Hr11's excess *nppa* is not explained by depth, and it is concentrated in low-depth cells.** In depth-matched bins, Hr11's mean *nppa* share is 0.0045, 0.0044, 0.0048 and 0.0041 in bins 0–3, against 0.0017, 0.0014, 0.0015 and 0.0017 in control cells. In the top bin it is 0.0041 against 0.0039. Rule 6 is met in bins 0–3. *nppb* shows no Hr11 excess. | [c002_c9_depth_matched.csv](../../../../results/c002_c9_depth_matched.csv) | Hr11 is 3 dpi and controls are uninjured, so injury biology is not excluded. An excess that is largest in low-count droplets and gone at high depth also fits ambient RNA, which adds a roughly fixed number of molecules per droplet: a **hypothesis**, untested. |
| P8 | ***ttn.2* (named in the source paper, Hu et al. 2022, so not found unaided; FU-6):** the same pattern as *nppa*. The Hr11 mean share is about 3.6–4.4× control in bins 0–3 and 1.1× in the top bin, and the gene is flagged "up with depth" (median rho 0.285). It was kept out of the axis sets. | C9 table; gene-depth table | Hinted. Its role as the slide's trap gene is not confirmed; `holdout/` was not read. |

## What this means for the direction

- None of the four sourced sets behaves as a distinct state axis in these cells by the fixed
  rules. The scores mostly reflect depth or detection, and random genes do as well on P90.
- Border-zone is the only axis near the coherence bar. That fits C001: it is the best-sourced
  axis.
- Proliferating fails coherence. Its genes are barely detected (C001), so a sourced generic
  cell-cycle set (FU-1) is needed before the axis can be judged. Dedifferentiated needs a
  redefinition (for example a sarcomere-down axis), as the supervisor anticipated.
- **score_genes should be the primary score,** with thinning as a standard control. UCell
  adds no independent check here (P3), and its larger depth tie is mostly technical (P5).
- The Hr11 *nppa* and *ttn.2* excess is the sharpest lead for the trap the slide names. The
  ambient-RNA test (R3.3) is the next decisive check: compare each gene's share in empty or
  low-count droplets, or use the non-CM marker share per library.

## Requirement coverage

- R2.1 scores for all 32,044 cells, both methods: done (`c002_scores.csv.gz`).
- R2.2 agreement: done. R2.3 matched null: done, plus an exploratory variant.
- R3.1 depth within library: done. R3.5 flags into YAMLs: done. R3.6: done, with correlation and
  thinning.
- Not done: doublets (R3.2), ambient RNA (R3.3), co-occurrence (G4), the mito covariate with NC
  added. All are out of this Probe's scope.

## Deviations and weaknesses

- UCell was re-implemented rather than run from pyUCell or AUCell. It is checked against a
  brute-force formula, but not against the package. Zeros are set to max_rank + 1, which equals
  UCell when a cell has fewer than 1,500 detected genes.
- `score_genes` uses its own internal control pool; the matched null is a second, outer control.
- Heart counts were added after the primary results, on supervisor instruction and not because of
  the results. The library thresholds were left unchanged.
- The P90 rule may miss rare states. Coherence pools all cells of a library.
- No subagents were used (FU-4).
- The F1 ZFIN check for `myh7l` came from a search summary of the ZFIN page, not an opened page.

## Next step proposed

1. The ambient-RNA test on the Hr11 *nppa* excess (R3.3), and a mito covariate with NC added.
2. For proliferating, a sourced generic cell-cycle set; for dedifferentiated, a sarcomere-down
   definition. Both need the FU-1 reading first.
3. Only then, co-occurrence, scored with score_genes and with thinning as a control.

## Status

**Review requested.**
