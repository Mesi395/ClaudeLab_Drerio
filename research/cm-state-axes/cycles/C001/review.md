# C001 review (supervisor)

Written 2026-11-07 15:30 CET. Report: [report.md](report.md), status Checkpoint.

## What I checked myself

Script: [../../checks/c001_verify.py](../../checks/c001_verify.py). Output:
[../../checks/c001_verify.out](../../checks/c001_verify.out). Evidence SUP-E5.

| Claim | Check | Result |
|---|---|---|
| C1: the sanitised copy has no labels | HDF5 keys of `data/clean/hu2022_cm_allgenes_sanitized.h5ad`, and a search of all paths for label-like names | **Confirmed.** obs holds only the 7 allowlisted columns. uns, obsm and obsp are empty. The only label-like hit is `orig.ident`, which is the sample name. |
| C8: Hr11 depth | Recomputed per-library median total counts from X | **Confirmed exactly.** Hr11 659, Hr7a 689, Hr2b 751, Hr2a 752, H8v 776. |
| C9: *nppa* and *ttn.2* lead | Recomputed detection, and added the median share of each cell's UMIs | **Confirmed and sharpened.** In Hr11, the median share of *nppa* is 0.00301 against 0.00075 in H7, about 4× higher, even though Hr11 has lower depth. For *ttn.2* it is 0.00140 against 0. So the rise is not only a detection effect. |
| `myh7` anomaly | Detection and UMI share of every myosin heavy chain symbol | **New finding.** `myh7l` is detected in 0.4772 of CMs and `myh6` in 0.3168, but `myh7` in only 0.0340. In this annotation the abundant ventricular myosin is probably `myh7l`, so the YAML mapping *vmhc* → `myh7` is likely wrong here. This is an inference; the symbol history still needs checking. |
| (not claimed) covariates | Top genes by total counts | **New finding.** The top gene is `NC-002333.4` (0.0582 of all UMIs), which looks like the mitochondrial genome contig (NC_002333 is the zebrafish mtDNA RefSeq). The `mt-` fraction may be undercounted if it excludes this contig. Hemoglobin `hbba1.1` is also among the top 15, a possible ambient-RNA sign. |

Not checked by me: the literature claims (C2–C7) beyond reading the evidence entries, and the
marker YAML contents gene by gene.

## Findings

- **F1 (R1.2, G1).** The *vmhc* symbol is likely mapped to the wrong gene (`myh7` instead of
  probably `myh7l`). This shows the "present in data" check can pass while scoring a nearly silent
  gene. It matters because a wrong symbol silently weakens an axis. It is closed by a sourced
  symbol check (ZFIN aliases or Ensembl IDs) for every renamed or ambiguous candidate. **Fix before
  scoring.**
- **F2 (R3.1, G3).** The mitochondrial fraction may miss `NC-002333.4`. It matters because mito
  fraction is a covariate in the confound tests. It is closed by identifying the contig from a
  source and stating whether `percent.mito` and the recomputed fraction include it. **Fix before
  using mito as a covariate.**
- **F3 (G3, the C9 lead).** The 4× higher *nppa* share in a low-depth library is the kind of
  signal the depth and ambient checks must explain. It goes into C002 as a test case, not as a
  conclusion.
- **F4 (G1, Understand gaps).** These remain unread: the Li supplements, the second halves of
  Honkoop and Li, the Wu, Ben-Yair and Sallin full texts, and a forward-citation search. They do
  not block a Probe, but they **block any direction decision or Deepen** (contract: no decision
  rests on an abstract or snippet). Logged as FU-1.

## Step back

**Signal.** Existing work does not give four separable marker sets. Border-zone is moderately
supported. Injury-responsive, dedifferentiated and proliferating are thin, overlapping, mistimed or
barely detected. That is a real and useful finding for the user's question "what do published
markers really mark?". It also means the co-occurrence question (G4) risks being decided by shared
genes rather than biology.

**What we know now that we did not.** The data is safe to use. The depth trap in cardiomyocytes is
narrower than the all-cell numbers suggested. The library with the lowest depth shows *higher*
levels of a top border-zone or stress gene, so depth is not the only technical axis; ambient or
library effects are a live alternative. The account of how far we are moved: G1 has drafts, and
G5 R5.3 is done.

**Best use of the next work.** The biggest open gap toward the goal is G2–G3: no axis has been
scored or tested against confounds. The hackathon deadlines leave limited time (feature freeze
16:15, submission 17:00). The researcher proposed a Probe on the border-zone axis only. I widen it
to all four candidate axes, because the scoring and checks are the same code, and the user's
question is which axes survive, not whether one does. Doublets, ambient RNA and co-occurrence come
after this signal.

## Decision

**Continue.** Next: [../C002/assignment.md](../C002/assignment.md) (Probe).
