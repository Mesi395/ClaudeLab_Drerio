# C008 report: re-checks within ventricle cells (researcher)

Written 2026-11-07 16:45 CET. Rules fixed at 16:38, before computing: [reading-rules.md](reading-rules.md).
Evidence R-E21. Data: 14,024 ventricle-called cells, 19 libraries with at least 50 cells, 18
hearts. All numbers come from [c008_summary.json](../../../../results/c008_summary.json) and
the per-library tables beside it.

## 1. Ties inside the ventricle ([c008_ties_by_library.csv](../../../../results/c008_ties_by_library.csv))

Median within-library Spearman rho:

| Axis | score_genes vs depth (unthinned) | vs depth (thinned) | vs *nppa* share (unthinned) | vs *nppa* share (thinned) | UCell vs depth | Libraries where UCell's depth tie > score_genes' |
|---|---|---|---|---|---|---|
| border-zone | 0.215 | −0.006 | 0.521 | 0.369 | 0.437 | 19/19 |
| dedifferentiated | 0.232 | −0.008 | 0.524 | 0.361 | 0.443 | 19/19 |
| injury-responsive | 0.227 | 0.058 | 0.577 | 0.361 | 0.388 | 19/19 |
| proliferating | 0.036 | −0.044 | −0.019 | −0.019 | 0.291 | 18/19 |

- **C8-1.** Inside the ventricle, the three nppa-led scores follow the *nppa* share more than depth.
  After thinning, the depth tie is gone (|rho| ≤ 0.06), but the *nppa* tie remains (0.36–0.37).
  This confirms and extends SUP-E12.
- **C8-2.** The C002 rank conclusion **stands** inside the ventricle: UCell is more tied to depth in
  18–19 of 19 libraries.

## 2. Matched null inside the ventricle ([c008_null_by_library.csv](../../../../results/c008_null_by_library.csv))

| Axis | P90 above null (libraries; hearts) | Coherence above null (libraries; hearts) |
|---|---|---|
| border-zone | 0/19; 0/18 | 10/19; 9/18 |
| dedifferentiated | 0/19; 0/18 | 6/19; 6/18 |
| injury-responsive | 0/19; 0/18 | 13/19; 12/18 |
| proliferating | 1/19; 1/18 | 2/19; 2/18 |

- **C8-3.** The C002 reading **stands**: no axis is above null on P90, and coherence is mixed or
  absent. It is slightly weaker for border-zone (10/19 against 17/22 on all CMs).

## 3. One program? ([c008_program_by_library.csv](../../../../results/c008_program_by_library.csv))

| Pair | (a) sets as they are: libraries above null; median rho (null p95) | (b) no *nppa*/*nppb*, disjoint |
|---|---|---|
| border \| dedifferentiated | 19/19 (18/18 hearts); 0.665 (0.154) | 0/19 (0/18); 0.071 (0.141) |
| border \| injury | 19/19 (18/18); 0.428 (0.139) | 0/19 (0/18); 0.054 (0.134) |
| injury \| dedifferentiated | 19/19 (18/18); 0.440 (0.103) | 2/19 (2/18); 0.076 (0.166) |

Variant (b) sets:

- border-zone: *myh7, hk1, pkma, pdk2a, ldha, mmp14b*;
- injury: *tnfrsf11b, hspb6, hspb11*;
- dedifferentiated: *desma, nkx2.5, tbx20, gata4, hand2, mef2cb, cacybp*.

Their median rho with the *nppa* share is 0.006, 0.061 and 0.036.

- **C8-4. Inside the ventricle, the three nppa-led axes are one signal, and that signal is
  *nppa*/*nppb*** (with the shared *mustn1b* and *ankrd1a*). Without those genes they do not
  co-vary beyond random gene sets, and they no longer follow *nppa*.
- In these data, the four "states" reduce to an *nppa*/*nppb* program and a cell-cycle set. The
  cell-cycle set is barely above null itself (§2; C007).
- G4 now has an answer within its scope: with the literature sets, the apparent co-occurrence of
  border, injury and dedifferentiation is shared genes (F7), not shared biology.

## 4. Ambient, ventricle cells ([c008_ambient_ventricle.csv](../../../../results/c008_ambient_ventricle.csv))

For Hr11 against control ventricle cells, by depth bin (0: <600, 1: 600–900, 2: 900–1,350, 3:
≥1,350 UMIs):

| Gene | Ratio by bin 0–3 | Excess molecules per cell, bins 0–3 |
|---|---|---|
| *nppa* | 7.31, 5.42, 4.70, 3.10 | 1.46, 1.75, 2.07, 1.57 |
| *myh6* (atrial) | 13.88, 4.79, 6.20, 3.63 | 0.15, 0.16, 0.28, 0.27 |
| *ttn.2* (hinted) | 2.47, 2.93, 3.64, 2.36 | 0.52, 0.69, 1.17, 1.02 |

Other observations:

- Non-CM genes: at least 2× in 10 of 12 in Hr11's shallow bins.
- Abundant CM genes are diluted (*myl7* 0.22, *gapdh* 0.24, *actc1a* 0.26).
- Other 3 dpi hearts, *nppa* in bins 0–1: Hr12 1.12 and 1.38; Hr24 0.25 (bin 1 only). Hr27 has
  fewer than 20 ventricle cells in bins 0–1.

The C003 verdict rule, applied to ventricle cells as pre-declared:

- Hr11-specific (0 other 3 dpi hearts at ≥2×): **yes**.
- Non-CM genes raised in ≥ half: **yes** (0.83).
- *nppa* ratio falls with depth: **yes** (7.31 → 3.10).
- Excess molecules within 2×: **yes** (1.46–2.07).
- Verdict: **"ambient favoured"**.

**C8-5. The C003 Hr11 conclusion changes.** On all CMs it was "inconclusive", and my reading then
was "probably not ambient; atrial soup does not explain *nppa*". Inside the ventricle:

- the excess has the fixed-amount signature of ambient RNA, about 1.5–2 extra *nppa* molecules
  per cell at every depth;
- the atrial gene *myh6* is raised 3.6–13.9× in Hr11 ventricle cells.

C003's *myh6* check was confounded by the atrial cells in the mixed population. The likely source
is an *nppa*-rich soup (atrial CMs express *nppa* highly). Caveats:

- The chamber call itself uses *myh6*. Strongly contaminated Hr11 cells may be called atrium and
  drop out, so the ratios may be underestimated, not created.
- Without empty-droplet data, this is still inference.
- *ttn.2* (hinted) has a less clean pattern (its ratio does not fall with depth).

## 5. Doublets, ventricle cells

The median |rho| between the scrublet score and the thinned score_genes axes is 0.059–0.069.
**C8-6.** The C005 conclusion stands: no axis tracks the doublet score. The automatic doublet call
still failed (FU-7).

## Where earlier conclusions change

| Earlier claim | Now |
|---|---|
| C002 P4/P5: UCell more depth-tied; thinning removes it | **Stands** in the ventricle |
| C002 P1/P2: not above null; coherence mixed | **Stands**; border-zone coherence is weaker in the ventricle |
| C002 axes "track depth" | **Narrowed**: inside the ventricle, the three nppa-led axes track the *nppa* share. Depth matters only before thinning |
| C003 Q2: Hr11 *nppa* inconclusive, "probably not ambient" | **Changed**: ambient favoured by the pre-declared rule in ventricle cells; the atrial-soup check is reversed |
| C004/C007 G4: cannot be resolved | **Sharpened**: the three nppa-led axes co-vary only through shared *nppa*/*nppb*; without them, no co-variation |
| C005 doublets: no axis tracks | **Stands** |

## Weaknesses

- Reproduced byte-identically by a rerun (R-E21).
- The chamber call depends on *myh6*, which interacts with the ambient test.
- Hr27 is missing from the shallow bins.
- Variant (b) sets are small (3–7 genes).
- The C002 per-gene depth flags were not re-run in the ventricle.

## Status

**Checkpoint** at the submission deadline. Done: items 1–5 of the rules. Missing: the per-gene
depth flags in the ventricle, and a check of the ambient result with a
chamber call that does not use *myh6*.
