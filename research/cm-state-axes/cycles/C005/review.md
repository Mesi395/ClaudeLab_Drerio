# C005 review (supervisor)

Written 2026-11-07 16:00 CET by the clock. Report: [report.md](report.md), status Review requested.
The rules were fixed at 15:48. The exploratory addendum at 15:52 is labelled exploratory. Accepted.

## What I checked myself

Script: [../../checks/c005_verify.py](../../checks/c005_verify.py). Output:
[../../checks/c005_verify.out](../../checks/c005_verify.out). Evidence SUP-E9.

| Claim | Result | Verdict |
|---|---|---|
| F10 fixed: states are scored only in ventricle-called cells | Crosstab of agent_state by agent_chamber. Every state call is in ventricle cells (border 67, dedifferentiated 53, injury 145, proliferating 183, unassigned 13,576). Atrial (8,816) and chamber-unassigned (9,204) cells are all `not_scored`. 32,044 unique cells. | **Confirmed.** F10 is closed. |
| Calls concentrate in Hr9 and Hr11 | State calls per heart: Hr9 149, Hr11 69, Hr8 44, Hr3 28, H6 27, Hr6 26 | **Confirmed.** One heart carries a third of all state calls. |
| Doublets: scrublet's automatic threshold called about 0 doublets in all 22 libraries | Per-library table: the all-cell doublet rate is 0.000–0.0004 at thresholds of 0.45–0.67, in libraries of up to 16,587 cells | **Confirmed as a threshold failure.** At those loadings, 10x doublet rates are expected to be several percent, so "0 doublets" is not credible. The report correctly calls this unanswered, not passed. |

## Findings

- **F11 (R3.2).** The doublet requirement is **partly met**. Continuous doublet scores show no
  axis tracking doublets by the fixed rule (median |rho| 0.050–0.091). The atrium libraries are
  positive (H8a about 0.33), which fits the chamber confound rather than doublets. The
  categorical part needs a threshold calibrated to each library's expected doublet rate, or a
  second method. Logged as FU-7. It does not block the next step.
- **The labels.** At or below chance for every state (border and dedifferentiated below; injury and
  proliferating at about chance). This is an honest representation of the evidence; the note for
  the `evaluate.py` runner says so. `agent_chamber` is the only label backed by validation.

## Step back

**Signal.** C005 closed the label defect and added one more confound test with a null result on
the continuous scale. The study now has a consistent answer: the sourced literature sets, after
depth and chamber control, carry at most a weak injury signal in a few hearts. They do not support
a four-state annotation, and they cannot answer co-occurrence.

**Is this still the best use of the work?** No more Probes of the same sets. Five Probes have
mapped what these sets can and cannot do. The remaining gap toward the goal is G1: better-founded
axes. That needs the direction decision deferred since C001. It requires:

- closing FU-1 (full reads, supplements, a forward-citation search, a sourced cell-cycle set,
  dedifferentiation transcript markers);
- a fair comparison of directions A (improved literature sets), B (spatial-derived from raw
  spatial data) and C (data-driven programs in ventricular CMs with count splitting against
  circularity, Neufeld 2024, matched to literature afterwards).

For the 17:00 submission, the user also needs R5.1: the inventory of checks the agent ran on its
own versus those a domain expert would add.

## Decision

**Continue.** Next: [../C006/assignment.md](../C006/assignment.md) (Understand, with the R5.1
inventory first).
