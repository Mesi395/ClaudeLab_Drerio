# C007 assignment: Probe of direction A

Study: `research/cm-state-axes/`. Assigned 2026-11-07 16:08 CET by the clock, by the Research
Supervisor. Read [../C006/review.md](../C006/review.md) and map section 6 first.

## Work kind

Probe.

## Question

Do axis sets built from Li 2025 Supplementary Data 2 do better than the C004 sourced sets on the
same ventricle-only positive control?

Build the sets from:

- the scRNA-seq "Cardiomyocytes (dediff.)" and "Proliferating cells" clusters;
- the Stereo-seq BZm domains.

Filter them to contrasts internal to cardiomyocytes (drop genes that mark non-CM types), make them
disjoint across axes, and fix a gene count per set in advance.

## Goal part and map gap

G1 and G2: better-founded axes. Direction A in map section 6.

## Decision it informs

Whether direction A earns a direction decision and a method review for Deepen, or whether
direction C (data-driven, with count splitting) becomes the candidate.

## Comparator

C004 on the original sets, run through the same pipeline and rules:

- ventricle-called cells;
- binomial thinning;
- `score_genes`;
- 18 hearts;
- matched null.

Report both side by side.

## Constraints

- Reading rules time-stamped first. That includes the set-building rule (thresholds, filter, gene
  count), so sets are not tuned on this data.
- Li 2025 used resection, and its stages differ from ours (6 hpa to 14 dpa against our 3, 7 and
  30 dpi). State the predicted time for each axis in advance.
- Calibrate the null threshold per heart where the C004 report found the tail rule miscalibrated.
  Declare that in advance too.
- Count hearts. Report per-heart values. Keep Hr9 visible, and show whether any pass survives
  leaving out Hr9.
- Use only the sanitised data. Use the Li files in `sources/li2025_supp/`.

## Evidence expected

- The sets with their provenance (a YAML or TSV with sheet, cluster and rank).
- A side-by-side positive-control table with the C004 sets.
- A leave-Hr9-out result.
- Results in `results/` with provenance.
- Map sections 3, 5 and 7 updated.

## Signal looked for

A clearer pass than C004: larger effects, not driven by one heart. Or a clear failure that sends
the study toward C.

## Timing

Submission is at 17:00. Report by about 16:40. Use `Checkpoint` if the work is unfinished.
