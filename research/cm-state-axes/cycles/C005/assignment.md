# C005 assignment: Probe and fix

Study: `research/cm-state-axes/`. Assigned 2026-11-07 15:49 CET by the clock, by the Research
Supervisor. Read [../C004/review.md](../C004/review.md) first.

## Work kind

A fix (F10), then a Probe (R3.2).

## 1. Fix F10 (first, required)

- Write a corrected label file, `results/c005_agent_labels.csv`. Apply the state rule only to
  ventricle-called cells. Set atrial and chamber-unassigned cells to `not_scored`.
- Keep `c004_agent_labels.csv` unchanged, and record the correction and its effect on the label
  counts.
- Document the format and the rule in the file's provenance.
- Add a short note for the human who runs `evaluate.py`. It should say what each column claims,
  which columns are near chance, and that `agent_chamber` comes from a sourced rule validated on
  six libraries.

## 2. Probe R3.2, doublets (if time allows)

**Question.** Do axis scores or agent state calls track doublet likelihood?

This requires the full per-library matrices (`scripts/get_data.sh`, GEO). CM-only data cannot see
CM–non-CM doublets. If the download or runtime does not fit before about 16:35, stop and report
`Checkpoint` with what is missing. Do not report doublets as checked.

Evidence expected:

- a per-library doublet score mapped to CM barcodes;
- within-library association with each axis score;
- the share of agent state calls among likely doublets.

## Constraints

- Reading rules time-stamped first.
- Use only the sanitised data plus the public GEO matrices. Never `holdout/`.
- Map sections 3, 5 and 7 updated.

## Timing

Submission is at 17:00. Report by about 16:35.
