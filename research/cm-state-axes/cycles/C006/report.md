# C006 report: Understand (researcher)

Written 2026-11-07 16:10 CET. No new analyses were run.

## Part 1: R5.1 checks inventory, done

[../../checks-inventory.md](../../checks-inventory.md) has 35 checks from C001–C005. For each it
gives who initiated it, the result with its evidence ID, and the trap it caught, followed by 11
checks a domain expert would add.

Candid points:

- The supervisor, not the researcher, found:
  - the embedded answer key;
  - the mito-contig gap;
  - the doubt about `myh7`;
  - the missing positive control;
  - the chamber confound (F8);
  - the atrial label defect (F10).
- The researcher saw signs of F8 and F10 (the *vmhc* chamber flag, H8's high tails, the H8/Hr3
  label concentration) but did not test or fix them.
- *ttn.2* is marked as hinted (FU-6).
- Traps the researcher caught unprompted:
  - method agreement gives false comfort;
  - rank scores are not depth-free (the user's claim, framed as a hypothesis by the supervisor);
  - the injury markers' timing and model mismatch;
  - single-gene meaning conflicts (*ccng1*);
  - the *nppa* excess in Hr11 (its cause is still open);
  - 3 dpi contamination, through a supervisor-designed test.

## Part 2: FU-1 and the direction comparison, partly done

- Forward-citation search done (OpenAlex). Two new relevant works were found, read as abstracts
  only (R-E19).
- Li 2025 Supplementary Data 2 and 5 opened. They give full sourced lists for border-zone
  domains and for scRNA-seq dediff. and prolif. CM clusters. They also source `myh7l` as
  ventricular, closing the chamber rule's gap.
- Direction comparison A, B, C drafted in map section 6, with contribution, uncertainty,
  comparator and cost. The researcher recommends A (Li SD2-based, CM-internal, disjoint sets) as
  the next Probe, Wu 2016 read in parallel, and C as the fallback.
- Not done (FU-1 still partly open): Wu 2016 full text and GSE74652; the second half of Honkoop
  2019; Ben-Yair 2019; Sallin 2015; Beisaw 2020; a sourced cell-cycle set; a full read of Neufeld
  2024.

## Weaknesses

- The citing works were filtered by title keywords, so some relevant works may be missed.
- The supplement sheets were read in part.
- The mapping of SD2 and SD5 to the MOESM file numbers is inferred.
- The recommendation rests on reading the tables, not on testing them.

## Status

**Checkpoint.** Part 1 is complete and ready for the submission. Part 2 is drafted; FU-1 is
partly open, as listed above.
