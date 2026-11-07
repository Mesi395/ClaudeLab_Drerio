# C008 assignment: Probe (re-check under the chamber confound)

Study: `research/cm-state-axes/`. Assigned 2026-11-07 16:38 CET by the clock, by the Research
Supervisor, after a user correction (FU-9). Read SUP-E12 first.

## Work kind

Probe.

## Question

Within ventricle-called cells only, what do the four axis scores track, and do the C002, C003 and
C005 conclusions still hold there?

## Why

The C002 results (matched null, depth ties, rank against score_genes, per-gene depth flags), the
C003 ambient panel and the C005 doublet correlations were computed on all cardiomyocytes,
although C003 found that the axes track chamber. SUP-E12 shows that inside the ventricle, the
border, injury and dedifferentiated scores follow the *nppa* share (median rho 0.52–0.58) more
than depth (0.22). Heart explains 0.26 of the variance.

## Decision it informs

- Which earlier findings stand as stated, and which must be narrowed.
- Whether, inside the ventricle, the three *nppa*-led axes are one program. If they are, the four
  states reduce to "*nppa* program" and "cell cycle" in these data.

## Constraints

- Reading rules time-stamped first.
- Ventricle-called cells; within library; hearts counted; thinning as the depth control; matched
  nulls.
- Test whether the three axes still co-vary after removing *nppa* and *nppb* from all sets
  (leave-gene-out).
- Report where the earlier conclusions change.

## Evidence expected

- The C002 tables re-run on ventricle cells.
- The leave-*nppa*/*nppb*-out co-variation.
- The ambient panel and doublet correlations on ventricle cells.
- Results with provenance.
- Map sections 3 and 5 updated.

## Timing

Submission is at 17:00. Report what fits by about 16:55 as `Checkpoint`; the rest can follow
after the submission.
