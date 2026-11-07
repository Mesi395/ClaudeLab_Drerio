# C008 reading rules (researcher), fixed 2026-11-07 16:38 CET before computing

Cells: ventricle-called (C004 rule B, `agent_chamber == ventricle`). Libraries with at least 50
ventricle cells are analysed (19 expected). Summaries count libraries and **hearts**: a heart passes
if all of its analysed libraries pass. Sets: `proposals/loop/*.yaml` (the C004 sets, minus ttn.2).

1. **Ties (C002 re-run).** Per library, Spearman of each axis score with:
   - log10 total counts (original depth);
   - genes detected;
   - *nppa* share (raw *nppa* / total);
   - mito share including `NC-002333.*`.

   Scores computed:
   - score_genes, unthinned;
   - score_genes on cells thinned to 500 UMIs;
   - UCell on raw counts.

   All scores are computed within the ventricle subset. The C002 thresholds apply: "tracks" if the
   median |rho| is ≥ 0.3, "weak" if 0.1–0.3. "Rank less tied to depth" holds if UCell's median
   |rho| with depth is lower than score_genes' for at least 3 of 4 axes.
2. **Matched null (C002 re-run).** Thinned score_genes. The C002 rules apply (P90 above the null 95th
   percentile; coherence above the null 95th percentile), with 100 sets per axis. "Above null" if at
   least 80% of libraries pass.
3. **Is it one program?** For the three nppa-led axes (border, injury, dedifferentiated), per
   library, take the Spearman between each pair of thinned score_genes axis scores.
   - Variant (a): the sets as they are.
   - Variant (b): *nppa* and *nppb* removed from all sets. The remaining sets are also made
     disjoint, which removes the shared *mustn1b* and *ankrd1a*.
   - Null: the Spearman between paired matched null sets (100 pairs).
   - "Co-vary" if the real rho is above the null 95th percentile in at least 80% of libraries and of
     hearts.
   - Also report each axis's rho with the *nppa* share in variant (b).
   - If (b) still co-varies, the axes share a program beyond *nppa*. If not, their earlier
     agreement was *nppa* itself.
4. **Ambient (C003 re-run, ventricle cells).** The C003 panels and fixed depth bins; Hr11, Hr12,
   Hr24 and Hr27 against control ventricle cells; ratios in bins 0–1. The C003 verdict rule applies.
5. **Doublets (C005 re-run, ventricle cells).** Per library, Spearman of the scrublet score
   (`results/c005_doublets_cells.csv`) with the thinned score_genes axis scores. "Tracks doublets"
   if the median |rho| is ≥ 0.1.

Anything not listed here is labelled exploratory.
