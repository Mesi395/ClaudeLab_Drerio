# C001 assignment: Understand

Study: `research/cm-state-axes/`. Assigned 2026-11-07 14:51 CET by the Research Supervisor.
Read [../../research-map.md](../../research-map.md) first, then `AGENTS.md`, `docs/problem.md`,
`docs/data.md` and `literature/methods-and-tools.md`.

## Work kind

Understand.

## Question

For each of the four states (injury-responsive, border-zone-like, dedifferentiated,
proliferating), what candidate marker genes does existing work offer? On what evidence: injury
model, timepoint, cell or spot, organism? How far does existing work get toward scoring these
states in cryoinjured zebrafish cardiomyocytes, and what does our dataset actually hold to score
them with?

## Goal part and map gap

- G1 (R1.1–R1.3): no marker sets exist yet.
- G5 (R5.3): the supplied data carries the answer key.
- Map sections 3, 4 and 7 are empty.

## Decision it informs

Which direction in map section 6 to take into a Probe: literature markers (A), spatially derived
markers (B), or a combination. It also shows whether any state lacks credible independent markers.

## Hard constraint first: the leakage guard (R5.3)

The `.h5ad` files in `data/hu2022/` embed the authors' labels (evidence SUP-E1).

- Before any other use of the data, write a script in `scripts/` or `src/zf/io.py`, with a test
  written first, that copies the data keeping **only** an allowlist:
  - obs: `orig.ident`, `nCount_RNA`, `nFeature_RNA`, `percent.mito`, `time`, `morphine`, `AV`
  - `layers['counts']`, `var_names` and the var stats
- Drop every other obs column (all annotation and clustering columns; their names are not listed
  in this public repository) and every `uns` entry.
- Never print, plot, tabulate or load the values of a dropped column. Work only from the
  sanitised copy afterwards.
- If any label value is seen, record it in `JOURNAL.md` and in your report, and stop.
- Do not open `holdout/`.

## Findings to resolve

- Li et al. 2025 has been read only as a partial summary (SUP-E4). Read it in full where
  possible, including Supplementary Data 2 and 5, or state the access failure.
- The user's claim that a rank-based score (AUCell) is "completely independent of normalization
  units and library size" is a hypothesis, not a fact. Note what the sources say about it.
  pyUCell implements UCell, not AUCell.
- The four Zotero CSVs in `data/zotero/` are identical copies (MD5). Use the one file as a reading
  list, and say what it covers.

## Evidence expected

- For each state, candidate marker genes, each with a source (DOI) and:
  - injury model, timepoint and resolution
  - whether the gene exists in the sanitised data under that zebrafish symbol
  - how often the gene is detected per library
- Existing work sorted by approach in map section 4, with how deeply each work was read.
- A search trail in map section 7: field terms, citations followed, what could not be opened.
- Basic data facts from the sanitised file, with provenance in `results/`:
  - per-library depth (total counts, genes detected), so the Hr11 depth trap is quantified
  - mitochondrial fraction
  - chemistry notes (the TASKS data-1 and data-2 checks)
- Draft `markers/*.yaml`, one per state, marked as candidates. Every gene needs a source.
- Map sections 3, 4, 5 and 7 updated. Where existing work, or a combination, already reaches a
  part of the goal, say so.

## Signal looked for

Whether each state has markers from at least two independent sources that are detectable in
these cells, and which states rest on thin or conflicting evidence.

This assignment is not after a preferred marker list. Do not choose or drop markers based on how
they behave in this dataset yet; that is later Probe work. The agent is meant to find the known
traps itself, so no candidates are named here.

## Timing

Hackathon schedule today: integration 15:30, feature freeze 16:15, submission 17:00. Report by
about 15:30. If work is unfinished, use status `Checkpoint` and name exactly what is missing.
