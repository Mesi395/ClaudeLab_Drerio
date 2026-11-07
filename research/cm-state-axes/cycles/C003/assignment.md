# C003 assignment: Probe

Study: `research/cm-state-axes/`. Assigned 2026-11-07 15:37 CET by the Research Supervisor.
Read [../C002/review.md](../C002/review.md) first.

## Work kind

Probe.

## Questions

1. **Positive control (F5).** Does any axis show the change the biology predicts? Its score should
   rise in injured hearts over uninjured hearts at the expected time (border zone and proliferation
   at 7 dpi; injury-responsive at 3 dpi), with depth controlled and **hearts counted** (18 hearts).
   Can a tail statistic see a state that a P90 per library cannot?
2. **The Hr11 lead (R3.3, ambient RNA).** Is Hr11's excess of *nppa* (and *ttn.2*, hinted, FU-6)
   better explained by ambient RNA than by 3 dpi biology?

## Goal part and map gap

- G2–G3 (R2.3, R3.3) and the G5 trap question.
- C002 found no axis above its null, but had no positive control.

## Decision it informs

- Whether any axis is usable as a state axis, which then gates co-occurrence (G4).
- Whether the slide's trap (a marker that tracks a technical artefact) was caught, and by which
  check.

## Constraints

- Reading rules time-stamped before computing. Within-heart or within-library statistics.
  `score_genes` is primary; report thinning as the depth control.
- For the time contrast, compare hearts, not cells. With 3–7 hearts per timepoint, report effect
  sizes and the per-heart values, not just p-values.
- Ambient evidence must be able to fail. For example, if ambient RNA explains the excess, other
  abundant CM genes should show a similar excess in Hr11's shallow cells, and non-CM genes (blood,
  fibroblast) should be raised there too. If only *nppa* and *ttn.2* rise, ambient RNA is the
  weaker explanation. Choose the tests; do not settle for a single plot.
- Hr12, Hr24 and Hr27 are the other 3 dpi hearts. Use them to separate "3 dpi" from "Hr11 only".
- Use only the sanitised data. A mito covariate must include `NC-002333.*` (F2).

## Evidence expected

- A per-heart table of depth-controlled axis scores by timepoint, with the contrast rule fixed in
  advance.
- A tail statistic per axis against the matched null.
- An ambient test table for Hr11 against the other 3 dpi hearts and the controls, covering *nppa*,
  *ttn.2* (hinted), other abundant CM genes and non-CM genes.
- A clear statement of which explanation the evidence favours, and what would overturn it.
- Results in `results/` with provenance.
- Map sections 3, 5 and 7 updated.

## Signal looked for

Either at least one axis passes a fair positive control, or none does. And the Hr11 excess is
either specific to Hr11 or shared with the other 3 dpi hearts.

## Timing

Feature freeze 16:15. Report by about 16:10. Use `Checkpoint` for whatever is unfinished, and state
exactly what is missing.
