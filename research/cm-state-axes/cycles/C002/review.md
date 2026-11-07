# C002 review (supervisor)

Written 2026-11-07 15:36 CET. Report: [report.md](report.md), status Review requested.

## Were the rules fixed first?

[reading-rules.md](reading-rules.md) states that its rules were fixed at 15:24. The primary outputs
are dated 15:27. The addendum dated 15:30 adds heart counts, which I asked for, and the thinning
rule, and both come before the thinning run. The exploratory null is labelled exploratory.
**Accepted.** The file's modification time (15:30) reflects the addendum.

## What I checked myself

Script: [../../checks/c002_verify.py](../../checks/c002_verify.py). Output:
[../../checks/c002_verify.out](../../checks/c002_verify.out). It is independent of the researcher's
code: scanpy `score_genes`, my own UCell, my own AUCell-style score, my own random gene sets and my
own depth bins. Evidence SUP-E6.

| Claim | My result | Verdict |
|---|---|---|
| P1: border axis not above matched null (score_genes P90) | 0 of 22 libraries above the 95th percentile of 50 matched random sets | **Confirmed** |
| P4: rank score more tied to depth than score_genes | Border axis, median within-library \|rho\| with log10 counts: score_genes 0.287, UCell 0.468, AUCell-style 0.427. Both rank scores are more tied than score_genes in 22/22 libraries. | **Confirmed, and extended to an AUCell-style score.** The AUCell-style score breaks zero ties at random, which is my assumption about AUCell. The decoupler or R AUCell package itself was not run. |
| P7: Hr11 *nppa* excess survives depth matching | Mean share against uninjured cells, binned on Hr11 depth quintiles: ratios of 2.53, 2.81, 3.10 and 3.13 in bins 0–3; **1.53 in the top bin** | **Confirmed in bins 0–3.** "Disappears at the top" is **overstated**: with my bins the excess shrinks to 1.53×, while the researcher's top bin gives 1.05×. The bin edges differ, so the top-bin claim depends on binning. The pattern (larger in shallow cells) holds. |

Not checked by me: P2 coherence, P5 thinning, P6 gene flags, F1 and F2 (the evidence entries read
as sound), and the hearts table.

## Findings

- **F5 (G2–G3, interpretation of P1 and P2).** The Probe has **no positive control**. "No axis
  beats random genes" could mean that the sets carry no state, or that the test cannot see a state
  held by a small fraction of cells. A P90 score per library is insensitive to a state present in a
  few percent of cells, and the border zone is a thin rim of tissue. A test is needed that should
  pass if the biology is there. Each axis's depth-controlled score should rise in injured hearts
  over uninjured hearts at its expected time (proliferation and border zone at 7 dpi), counted per
  heart. A tail statistic such as the fraction of cells above the null's 99th percentile also helps.
  This matters because, without the control, P1 cannot tell "fake markers" from "insensitive test".
  **Fix before any claim that the sets carry no state.**
- **F6 (G3, P7).** The top-bin "disappears" claim depends on binning; my bins give 1.53×. Report
  the excess as "largest in shallow cells, smaller in deep cells".
- **F7 (G4, the co-occurrence risk).** Three axes share *nppa* and *nppb*, and all axes are tied to
  depth. Co-occurrence computed on these sets would be inflated by construction. G4 needs disjoint
  sets, depth-controlled scores and a matched-null comparison. This is not urgent until F5 is
  answered.
- **Cleared:** F1 (symbol settled from ZFIN; the low detection of `myh7` stays open as a note) and
  F2 (NC features identified; the mito covariate must include them).

## Step back

**Signal.**

- Strong and checked: in these sparse cardiomyocytes, gene-set scores track depth. Rank methods
  track it *more*, not less. Agreement between methods is near total, so it is not an independent
  check.
- This answers the user's scoring question directly: `score_genes`, with a thinning control, is the
  better primary score here, and rank-based independence from depth is refuted for these data.
- Weak-to-negative, and not yet interpretable without F5: whether the sourced sets mark states.
  Border-zone comes closest on coherence.
- The trap lead sharpened: an *nppa* excess in the shallowest library that is largest in its
  shallowest cells, consistent with ambient RNA, though not tested.

**What we know now that we did not.** The four-axis picture from the literature does not appear by
default in these data, and the obvious checks (method agreement, rank scoring) give false comfort.
The way forward is to make each axis pass a positive control and to explain the Hr11 excess.

**Direction.** Still Probe-level, and no direction decision yet (FU-1 is open). Under the deadline
(freeze 16:15, submission 17:00), the most decision-relevant next work is: (1) the positive
control (F5), which decides whether any axis is usable at all; (2) the ambient test of the Hr11
lead (R3.3), which decides whether the slide's trap was caught. Redefining axes from FU-1 reading
is valuable, but it comes after these two because it cannot finish before the freeze.

## Decision

**Continue.** Next: [../C003/assignment.md](../C003/assignment.md) (Probe).
