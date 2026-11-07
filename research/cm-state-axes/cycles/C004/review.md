# C004 review (supervisor)

Written 2026-11-07 15:48 CET by the clock. Report: [report.md](report.md), status Review requested.
The rules were fixed at 15:43, before computing. Accepted.

## What I checked myself

Script: [../../checks/c004_verify.py](../../checks/c004_verify.py). Output:
[../../checks/c004_verify.out](../../checks/c004_verify.out). Evidence SUP-E8.

| Claim | Result | Verdict |
|---|---|---|
| The label file covers every cell, once | 32,044 rows, 32,044 unique, and the set equals the sanitised obs index. Counts match the report: unassigned 30,064, injury 1,188, proliferating 407, dedifferentiated 199, border 186. | **Confirmed** |
| Chamber call: balanced accuracy 0.991, coverage 0.735 | From the label file against AV in the six chamber-only libraries: accuracy 1.000 for atrium and 0.973 for ventricle among called cells, balanced 0.986, coverage 0.701 | **Confirmed in substance.** It differs slightly from the reported figures, probably because the report averages per library and I pool cells. The rule was chosen on the same libraries, so the figure is optimistic, as the report says. |
| The labels can go to the human-run agreement check | `injury_responsive` calls by chamber: atrium 885, ventricle 184, unassigned 119. Share of cells called injury-responsive by time: Ctrl 0.050, 3 dpi 0.027, 7 dpi 0.022, 30 dpi 0.074. | **Problem (F10).** Three-quarters of the injury calls are atrial cells, and they are most frequent in uninjured and 30 dpi hearts. The label rule applies the axes to every chamber, although the positive control was run in ventricle only. The file inherits the chamber trap (F8). |

Not rechecked: the ventricle-only positive-control numbers and the G4 table. The report's own
caveats, that the passes are near null and the median moves the wrong way, already make them
non-decisive.

Descriptive, from the label file: border and dedifferentiated calls are most frequent at 3 dpi
(0.0102 and 0.0112, against 0.0031 and 0.0050 in Ctrl), and proliferating calls at 7 dpi (0.0172
against 0.0116). These are small, in the predicted direction, and not tested.

## Findings

- **F10 (R5.2).** The agent labels must not call atrial cells into the injury states without
  saying so. Fix: apply the state rule only to ventricle-called cells (the injury is ventricular,
  and the positive control was ventricle-only), set atrial and unassigned-chamber cells to
  `not_scored`, and document the rule change as a correction. This matters because otherwise
  the human-run check would mostly score chamber, not states. **Fix before the file goes to the
  human check.**
- **F5 (positive control).** Partly resolved. Removing chamber raised border from 0.571 to 0.786
  and gave weak passes for dedifferentiated (7 dpi) and injury-responsive (3 dpi). All effects are
  near null level, and only Hr9 is clearly above null. Per-heart null calibration is the open
  method issue.
- **G4.** Accepted as "cannot be resolved with these sourced sets". Overlapping sets co-occur by
  construction, and disjoint sets look like random pairs (all 6 pairs mixed). This is an honest
  negative with a stated scope; it is not a finding that the states are independent.
- **R3.2 doublets.** Not run. It is a slide requirement. CM-only matrices are blind to CM–non-CM
  doublets, so it needs the full per-library GEO matrices (`scripts/get_data.sh`).

## Step back

**Signal.** The study's strongest results are now about the markers and methods, not the states:

- Scores track depth, and rank methods more so (C002).
- The published border and injury sets track chamber (C003, F8).
- Hr11 carries a contamination and stress signature (C003).
- With these two confounds removed, the sourced sets carry at most a weak, near-null injury signal
  in ventricular CMs (C004).

So the co-occurrence question cannot be answered with literature sets alone. That is a real
answer for the user's "what do published markers really mark?" and for the slide's question about
trusting an annotation without ground truth.

**How far we are.** G1 sets are drafted with flags. G2 is done. G3 has depth and chamber covered,
ambient partly, and doublets not run. G4 is an honest negative. For G5, the label file exists but
needs F10. FU-1 reading is still open, so no direction decision was made. This study has run as a
series of Probes and has not reached Deepen.

**Best use of the remaining time** (submission 17:00): fix F10, because the human check depends
on it, then run doublets (R3.2) if the GEO download fits.

## Decision

**Continue.** Next: [../C005/assignment.md](../C005/assignment.md).
