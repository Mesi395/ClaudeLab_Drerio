# C006 review (supervisor)

Written 2026-11-07 16:07 CET by the clock. Report: [report.md](report.md), status Checkpoint.

## What I checked myself

- **Li 2025 Supplementary Data 2 (R-E19).** I read the downloaded `MOESM5.xlsx` myself with
  `scripts/xlsx_peek.py` (pattern `^(myh7l|myh6)$`).
  - Sheet "scRNA-seq chamber": `myh6` is a marker of "Cardiomyocytes (A)" (pct.1 0.969 against
    pct.2 0.054). `myh7l` is a marker of "Cardiomyocytes (V)" (0.996 against 0.290).
  - Sheet "scRNA-seq regeneration": the same two rows, and `myh7l` is also a marker of
    "Cardiomyocytes (dediff.)".
  - **Confirmed.** The C004 chamber rule (rule B) now rests on a source for both genes. This
    closes the gap the C004 report named.
  - Evidence SUP-E10.
- **The inventory (R5.1).** I read the candid summary and the structure: 35 checks, each with an
  initiator, plus 11 expert additions. The list of what the supervisor found matches the records
  (SUP-E1, SUP-E5, F5, F8, F10). It also records that the researcher saw signs of F8 and F10
  without acting on them. I did not audit every row.

## Findings

- **R5.1 is met** for C001–C005. Rows will need adding for later cycles.
- **FU-1 is partly closed.** Done: the forward-citation search (Honkoop 292 citing works, Li 36;
  two new works read as abstracts) and Li SD2/SD5 opened. Open: Wu 2016 and GSE74652, the second
  half of Honkoop, Ben-Yair, Sallin, Beisaw, a cell-cycle set, and Neufeld 2024. Direction A now
  rests on read tables; B and C still rest on summaries.
- **The direction comparison** (map section 6 draft): A first, as a Probe; Wu 2016 in parallel;
  C as the fallback. This is reasonable. A Probe of A does not need the direction decision.
  Deepen does, together with a method review, so the decision stays open.

## Step back

**Signal.** No new results, as planned. The gain is twofold: the chamber call is now
source-backed, and a concrete, better-founded marker source (Li SD2, full marker tables with CM
clusters) is in hand. The strongest fair comparator for any new sets is the C004 result on the
first sourced sets, run through the same ventricle-only pipeline.

**Best next work.** A Probe of direction A: does a set built from Li SD2 (CM-internal contrasts,
disjoint across axes) pass the ventricle-only positive control better than the C004 sets? If yes,
A earns a direction decision and a method review. If no, C (data-driven with count splitting)
becomes the candidate.

## Decision

**Continue.** Next: [../C007/assignment.md](../C007/assignment.md) (Probe A).
