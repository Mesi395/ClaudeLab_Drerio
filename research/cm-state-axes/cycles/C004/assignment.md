# C004 assignment: Probe

Study: `research/cm-state-axes/`. Assigned 2026-11-07 15:48 CET by the Research Supervisor.
Read [../C003/review.md](../C003/review.md) first (finding F8).

## Work kind

Probe.

## Questions

1. **Chamber call.** Can a per-cell atrium/ventricle call, built from sourced chamber genes, be
   validated on the six chamber-only libraries? H8a, Hr6a and Hr7a are atrium; H8v, Hr6v and
   Hr7v are ventricle. These libraries are a labelled set for chamber from metadata, not from
   the expert annotation.
2. **Positive control, in ventricular CMs.** With cells thinned to equal depth, `score_genes` and
   18 hearts counted: does any axis rise in injured hearts over controls at its predicted time?
3. **If an axis passes,** a first co-occurrence result (G4) within ventricular CMs, using
   **disjoint** gene sets (no gene on two axes), per heart, against a matched null. If none
   passes, state that G4 cannot be answered with these sets, and show why using the shared genes
   and the chamber and depth ties.
4. **Agent label file (R5.2).** Write the agent's own per-cell calls (barcode plus the agent's
   state calls, including "unassigned") in the format `src/zf/evaluate.py` expects. If the stub
   defines no format, write one and document it. Never read `holdout/`.

## Goal part and map gap

G2–G4 and G5 R5.2. F8 showed that the axes track chamber, which confounds every comparison
between hearts.

## Decision it informs

Whether any axis is usable as a state axis. What answer G4 gets. What goes to the human-run
agreement check.

## Constraints

- Reading rules time-stamped before computing.
- Chamber genes must come from a source (DOI or ZFIN), not from this data's clusters.
  Report validation accuracy on the six libraries and the threshold rule, fixed in advance.
- Count hearts. Report effect sizes and per-heart values.
- Labels must follow a rule fixed in advance, from your own scores. Never infer them from, or
  tune them to, any annotation.
- Use only the sanitised data.

## Evidence expected

- A chamber-call validation table.
- A per-heart ventricle-only positive-control table.
- A G4 result or a documented "cannot answer".
- The label file plus its format note.
- Results in `results/` with provenance.
- Map sections 3, 5 and 7 updated.

## Timing

Feature freeze 16:15. Report by about 16:15. If items 3 and 4 do not fit, deliver item 4 first,
because the label file must exist for submission. Use `Checkpoint` for the rest.
