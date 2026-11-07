# C007 review (supervisor)

Written 2026-11-07 16:20 CET by the clock. Report: [report.md](report.md), status Review requested.
The rules were fixed at 16:06, before the sets or scores existed. Accepted.

## What I checked myself

Script: [../../checks/c007_verify.py](../../checks/c007_verify.py). Output:
[../../checks/c007_verify.out](../../checks/c007_verify.out). Evidence SUP-E11.

The check is score-free and independent: the share of ventricle-called cells with two or more
genes of the Li SD2 proliferating set detected, per heart, raw and after thinning to 500 UMIs.

| Heart | Share, raw | Share, thinned |
|---|---|---|
| Hr9 (7 dpi) | 0.0296 | 0.0068 |
| Hr8 (7 dpi) | 0.0104 | 0.0044 |
| Other 7 dpi hearts | 0.0000–0.0037 | |
| Controls | 0.0000–0.0023 | |

**C7-3 is confirmed:** the cell-cycle signal sits in Hr9 and Hr8. I did not recheck the AUC table
(C7-1, C7-2). It is consistent with this result.

## Decision

Direction A, as built from SD2, fails the Probe. Direction C (data-driven programs with count
splitting, per-heart covariates) is the next candidate.

The study stops here at the hackathon submission, by the user's call, with status **not
complete**. Restart point: the decision between direction C and a validated positive-control
dataset (GSE139218), after FU-1 and FU-7. Summary for submission:
[docs/research-methodology.md](../../../../docs/research-methodology.md) and [docs/research-report.md](../../../../docs/research-report.md).
