# C002 assignment: Probe

Study: `research/cm-state-axes/`. Assigned 2026-11-07 15:30 CET by the Research Supervisor.
Read the map, then [../C001/review.md](../C001/review.md).

## Work kind

Probe.

## Question

Do the four candidate axes, scored as sourced, carry signal beyond random gene sets and beyond
sequencing depth within each library? Is a rank-based score less tied to depth than `score_genes`?

## Goal part and map gap

- G2 R2.1–R2.3 and G3 R3.1, R3.5 and R3.6: nothing scored yet.
- This is the core trap the slide names: a marker that tracks library size.

## Decision it informs

- Which axes are worth carrying into doublet, ambient and co-occurrence tests.
- Whether a border-zone-centred direction (A+B) is viable, or whether the states need to be
  redefined, for example with a sarcomere-down axis for dedifferentiation or a generic cell-cycle
  set for proliferation.
- Whether `score_genes` or a rank method should be the primary score.

## Findings to resolve first

From the C001 review:

- **F1.** Verify the symbol for *vmhc*, and any other renamed or ambiguous candidate, from a
  source (ZFIN or Ensembl). The supervisor saw `myh7l` detected in 48% of CMs against `myh7` in
  3.4%. Correct the YAMLs and record the correction.
- **F2.** Identify `NC-002333.4`, which carries 5.8% of all UMIs and is probably the mitochondrial
  contig. State whether the mito fraction includes it before using mito as a covariate.

## Constraints

- **Fix the reading rules before computing.** Write in your report, time-stamped, what counts as
  "above null", "tracks depth" and "the methods agree". Changes after seeing results are labelled
  exploratory.
- Statistics are within library, never pooled across libraries.
- Do not choose or drop markers by how they behave here without a recorded rule; flag them rather
  than drop them.
- Keep circular genes such as *ttn.2* out of the axis sets. Testing them as single genes is allowed.
- Use only the sanitised data. Use uv only (`uv add`), and write a test first for any code under
  `src/zf/`.

## Evidence expected

- Per-cell scores for each axis with `score_genes` and one rank method (AUCell via decoupler, or
  UCell).
- Per-library agreement between the two methods.
- An expression-matched random gene-set null per axis.
- The within-library association of each axis score with log total counts and genes detected.
- The R3.6 rank-versus-depth test, including Hr11 against libraries of similar depth (Hr7a, Hr2a,
  Hr2b).
- Per-gene depth association for every axis gene, with flags written into the YAMLs.
- The C9 lead (*nppa* and *ttn.2* higher in Hr11) examined against depth.
- Tables in `results/` with provenance.
- Map sections 3, 5 and 7 updated.
- The reading rules stated before the numbers.

## Signal looked for

Which axes stay above their null and stay weakly tied to depth within libraries under both
methods, and which do not. Either answer is useful.

## Timing

Feature freeze is 16:15. Report by about 16:05 with `Review requested` or `Checkpoint`.
