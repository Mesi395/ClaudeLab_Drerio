# C003 review (supervisor)

Written 2026-11-07 15:47 CET. Report: [report.md](report.md), status Review requested.
The rules were fixed at 15:36, before computing (reading-rules.md). Accepted.

## What I checked myself

Script: [../../checks/c003_verify.py](../../checks/c003_verify.py). Output:
[../../checks/c003_verify.out](../../checks/c003_verify.out). Evidence SUP-E7.

I chose a check that could overturn the reading of Q1. C003 found that no axis passes the
positive control. That could mean the markers carry no injury state, or that a confound between
hearts hides it. The control heart H8 has the highest tails, and H8 pools an atrium-only and a
ventricle-only library, so I tested chamber.

| Check | Result |
|---|---|
| Axis scores in the three hearts split into atrium-only and ventricle-only libraries (AV is fixed per library from metadata) | Border-zone median score_genes: H8a 0.679 vs H8v 0.068; Hr6a 0.359 vs Hr6v 0.058; Hr7a 0.345 vs Hr7v −0.004. Injury-responsive: 0.876 vs 0.404, 0.313 vs 0.189, 0.371 vs 0.167. **Within the same heart, atrial CMs score far higher.** |
| A crude per-cell chamber call, myh6 share > myh7l share (my assumption, not validated against a source) | Calls 0.661 of cells in atrium-only libraries atrial, against 0.021 in ventricle-only libraries. Across hearts the atrial fraction ranges from 0.043 (Hr19) to 0.827 (Hr27). |
| Per-heart axis score against atrial fraction, 18 hearts | Spearman rho: border 0.26, injury-responsive 0.42 |
| Border axis in ventricle-called cells only, per heart (exploratory, unthinned) | The top four hearts are injured: Hr9 0.637, Hr8 0.408, Hr11 0.397, Hr12 0.244. Controls range 0.072–0.154. By my count the AUC is about 0.57 for 7 dpi against Ctrl and about 0.75 for 3 dpi against Ctrl. |

I did not recheck Q2's numbers (the ambient panel). The reasoning in the report is careful, and its
verdict, "inconclusive, Hr11-specific technical signature, nppa probably not ambient", is not what
the next decision rests on.

## Findings

- **F8 (G1–G4, central).** The border-zone, injury-responsive and dedifferentiated axes all lean on
  *nppa* and *nppb*, and they track **chamber identity**: atrial CMs score 5–10× higher than
  ventricular CMs of the same heart. Whole-heart libraries differ widely in atrial fraction. So
  every comparison between hearts, including C003's positive control, mixes injury with chamber
  composition. Every co-occurrence between these axes would also be inflated by atrial cells.
  This matters because it is the largest signal these axes carry, and because it answers the
  user's question "what do published markers really mark?": in uninjured hearts, the sourced
  "border-zone" genes mostly mark atrium. It is closed by:
  - a per-cell chamber call **validated** against the six chamber-only libraries, which are a
    true labelled set for chamber that comes from metadata, not from the expert annotation;
  - repeating the positive control within ventricular CMs (injury is ventricular), with thinning
    and hearts counted.
- **F5 (partly open).** Read through F8, C003's "no axis passes" is not evidence that the markers
  are fake. The ventricle-only hint (the injured hearts at the top) is exploratory, uses my crude
  chamber call, and is unthinned.
- **F9 (Q2).** Accepted as inconclusive. Hr11 carries a technical signature: non-CM genes raised,
  *hsp70l* about 5×, mito plus NC about 1.5×. The *nppa* excess is better matched by an
  exaggerated 3 dpi response (Hr12 is at 1.8–2.1×) than by fixed-amount ambient RNA. Chamber also
  matters here, because *nppa* is atrial: Hr11's atrial fraction is 0.287, the same as control H7.
  So chamber does not explain the Hr11 excess. Doublets and an empty-droplet soup profile remain
  untested; they need the raw GEO matrices.
- **The proliferating secondary pass** (AUC 0.82 on medians) rests on 3 hearts with high *ccng1*
  and no leave-one-out. Not decision-grade.

## Step back

**Signal.** C003 gave a negative on its own terms. My check turns it into the clearest positive
lead of the study: the published border and injury sets are dominated by an atrial identity
signal, and once restricted to ventricular cells, injured hearts tend to score highest. Two
technical traps are now documented: depth (C002) and chamber composition (C003). Hr11 adds a
third, contamination and stress.

**What changed in how far we are.**

- G3 has its central finding.
- G4 cannot be answered on the raw axes; it needs ventricle-only, disjoint, depth-controlled
  scores.
- G1 needs a chamber column in the marker records.

**Direction.** This is still Probe work under deadline. The next decision-relevant step is to
validate the chamber call and rerun the positive control in ventricular CMs. That one result
decides whether any axis is usable, and therefore whether G4 gets a real answer or a documented
"cannot be answered with these sets". The agent's label file (R5.2) is needed for the
human-run agreement check and must exist by submission.

## Decision

**Continue.** Next: [../C004/assignment.md](../C004/assignment.md) (Probe).
