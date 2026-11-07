# C008 review (supervisor)

Written 2026-11-07 16:44 CET by the clock. Report: [report.md](report.md), status Checkpoint.
The rules were fixed at 16:38. Accepted.

## What I checked myself

Script: [../../checks/c008b_verify.py](../../checks/c008b_verify.py). Output:
[../../checks/c008b_verify.out](../../checks/c008b_verify.out). Evidence SUP-E13. My own disjoint
sets and my own depth bins.

| Claim | My result | Verdict |
|---|---|---|
| C8-4: without *nppa* and *nppb*, the three axes no longer co-vary | Median within-library rho over 19 libraries: border–injury 0.085, border–dediff 0.083, injury–dediff 0.121 | **Confirmed.** Near the reported null p95 of 0.10–0.15. Inside the ventricle, the three "states" are one *nppa*/*nppb* program. |
| C8-5: Hr11's *nppa* excess in ventricle cells fits ambient RNA | Excess of +1.46, +1.42, +1.67, +2.02 and +1.91 molecules per cell across Hr11 depth quintiles, within 2×. Ratios 7.62, 5.36, 4.71, 6.54 and 4.37 (not monotonic). *myh6* 3.93–21.28×. *hbba1.1* 0.55–1.24× (not raised). | **Partly confirmed.** The fixed-molecule pattern and the atrial *myh6* in ventricle cells support an atrial or CM-derived soup. The "steadily falling" ratio is not reproduced in my bins. Blood RNA is not raised, so this is not a general contamination. |

## Decision

C003's "atrial soup ruled out" is **corrected**: it was confounded by atrial cells in the mixed
population. The study's best current reading is that Hr11's *nppa* excess is ambient RNA,
probably atrial in origin, with moderate confidence. Empty droplets and a chamber call that does
not use *myh6* are needed to settle it.

The submission is at 17:00, so the study pauses here as **not complete**. The open items for C008
(per-gene depth flags in the ventricle, an ambient check with an independent chamber call) carry
over.
