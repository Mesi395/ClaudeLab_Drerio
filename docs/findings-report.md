# Cardiomyocyte state axes after cryoinjury: methodology and findings

Study `cm-state-axes`, run on 2026-11-07 by a two-agent research loop (Research Supervisor and
Researcher, Claude). Status at submission: **not complete**. Seven cycles were run and reviewed,
and no Deepen phase was reached (see "Status and limits"). Every number below is quoted from a file
in `results/` or `research/cm-state-axes/checks/`, linked beside it. The full record is the
[research map](../research/cm-state-axes/research-map.md).

## 1. The question

After cryoinjury, zebrafish cardiomyocytes are said to take on four states: injury-responsive,
border-zone-like, dedifferentiated and proliferating. With no ground truth, can an agent score
every cell on each state, test the scores against technical confounds, and say whether the states
co-occur? Data: Hu et al. 2022, 32,044 cardiomyocytes, 18,794 genes, 22 libraries, which are 18
hearts (uninjured, 3, 7 and 30 dpi).

## 2. Methodology: how the loop works

### Roles and records

- **Research Supervisor.** Holds the goal. Writes one assignment per cycle. Rechecks the claims
  each decision rests on with its own code in `research/cm-state-axes/checks/`. Reviews the
  method. Decides whether to continue.
- **Researcher.** Reads literature, writes code test-first, runs the analyses and reports claims
  with evidence IDs.
- **Shared records.** The [research map](../research/cm-state-axes/research-map.md) (goal, state,
  decisions), per-cycle `assignment.md`, `report.md`, `review.md` and `reading-rules.md`, an
  append-only [evidence log](../research/cm-state-axes/evidence.md),
  [follow-ups](../research/cm-state-axes/follow-ups.md), and [JOURNAL.md](../JOURNAL.md).

### Cycle kinds

- **Understand:** read the existing work and build the map.
- **Probe:** a cheap, faithful test of a direction.
- **Deepen:** serious building or validation. It needs a saved direction decision and a method
  review; it was not reached.

### Rules that guard the evidence

1. **Leakage guard first.** The supplied `.h5ad` carried the authors' labels. Before any analysis,
   a tested allowlist sanitiser (`src/zf/io.py::sanitize_h5ad`, `tests/test_io_sanitize.py`)
   wrote a label-free copy. It reads only the allowlisted columns through h5py, so the labels are
   never loaded. The supervisor verified the output keys (evidence SUP-E5). `holdout/` was never
   opened.
2. **Reading rules fixed before computing**, in each cycle's `reading-rules.md`, time-stamped.
   Changes made after seeing results are labelled exploratory.
3. **Within-library or within-heart statistics**, with hearts as the unit (18), never pooled.
4. **Controls:**
   - expression-matched random gene sets as the null;
   - binomial thinning of every cell to 500 UMIs, as the depth control;
   - a positive control (scores must rise at the biologically predicted time);
   - chamber-only libraries as a labelled set for chamber, taken from metadata, not from the
     expert annotation.
5. **Independent supervisor rechecks** of the decisive claims, with separate code
   (`checks/c00*_verify.py`, outputs `.out`). Results the supervisor confirmed, corrected or
   could not check are marked in each `review.md`.
6. **Provenance.** Every result has a sidecar with the commit, command and input. Primary runs
   were rerun and reproduced byte-identically (researcher).

### Scoring pipeline

- `score_genes` (scanpy) on log1p CP10k as the primary score.
- UCell, re-implemented and tested against a brute-force formula, as the rank comparator.
- AUCell-style scoring in the supervisor check only.
- Marker sets from literature (`markers/*.yaml`, each gene with a DOI) and from Li 2025
  Supplementary Data 2 (`markers/li2025_sd2_sets.tsv`).

### Cycles

| Cycle | Kind | Question | Result (evidence) |
|---|---|---|---|
| C001 | Understand | What markers does existing work offer, and what does the data hold? | Sanitised data; the literature gives overlapping sets; Hr11 has the lowest depth (R-E1–R-E8, SUP-E5) |
| C002 | Probe | Do the four sets beat matched nulls and depth? Is rank scoring less tied to depth? | No set beats null; rank scores are *more* tied to depth (SUP-E6) |
| C003 | Probe | Positive control; is the Hr11 *nppa* excess ambient RNA? | No pass. The supervisor found the chamber confound (SUP-E7) |
| C004 | Probe | Chamber call; ventricle-only positive control; co-occurrence; labels | Chamber call validated; weak passes; co-occurrence unresolved (SUP-E8) |
| C005 | Fix and Probe | Ventricle-only labels; doublets | Labels fixed; doublet threshold failed (SUP-E9) |
| C006 | Understand | Inventory of checks; close the reading gaps | [checks-inventory.md](../research/cm-state-axes/checks-inventory.md); Li SD2 read (SUP-E10) |
| C007 | Probe | Do sets built from Li SD2 beat the C004 sets? | No; signal is concentrated in two hearts |

## 3. Findings

**F1. Rank-based scores are *more* tied to sequencing depth than `score_genes`, not less.**

- Border-zone axis, median within-library |rho| with log10 counts
  ([check output](../research/cm-state-axes/checks/c002_verify.out)): score_genes 0.287, UCell
  0.468, AUCell-style 0.427.
- Both rank scores are more tied in 22 of 22 libraries.
- All four axes: [results/c002_axis_summary.csv](../results/c002_axis_summary.csv).
- Thinning to equal depth removes most of UCell's tie, so the tie is technical
  ([results/c002_thinning_summary.csv](../results/c002_thinning_summary.csv)).
- Confidence: high. The supervisor reproduced it independently. The AUCell package itself was not
  run.

**F2. The published border-zone and injury marker sets mainly mark the atrium.**

- Median border-zone score, atrium-only against ventricle-only library of the same heart
  ([check output](../research/cm-state-axes/checks/c003_verify.out)): H8 0.679 vs 0.068; Hr6
  0.359 vs 0.058; Hr7 0.345 vs −0.004.
- The atrial fraction per heart ranges from 0.043 to 0.827, which confounds every comparison
  between hearts.
- The chamber genes are source-backed: `myh6` atrial and `myh7l` ventricular, from Li 2025 SD2
  (SUP-E10).
- The per-cell chamber call reaches balanced accuracy 0.986 at coverage 0.701 on the six
  chamber-only libraries ([check output](../research/cm-state-axes/checks/c004_verify.out)).
  That is optimistic, because the rule was chosen on the same libraries.
- Confidence: high for the within-heart contrast.

**F3. With depth and chamber controlled, no sourced set shows the predicted time signal above
random genes.**

- Ventricle-only cells, thinned, 18 hearts, a tail threshold calibrated on control hearts, and a
  pass requiring an AUC above the matched-null 95th percentile: 0 of 6 contrasts pass for either
  the C004 sets or the Li SD2 sets ([results/c007_positive_control.csv](../results/c007_positive_control.csv)).
- The two C004 passes did not survive calibration.
- Confidence: moderate. With 4 control hearts the test has little power, and its sensitivity is
  not validated against a known positive (F5 in the C002 review).

**F4. The remaining signal sits in a few hearts, not in timepoints.**

- Share of ventricle cells with two or more cell-cycle genes from the Li SD2 set detected
  ([check output](../research/cm-state-axes/checks/c007_verify.out)): Hr9 0.0296 and Hr8 0.0104.
  Every other 7 dpi heart is at or below 0.0037; three controls are 0.0000 and H6 is 0.0023.
- Hr9 and Hr8 also have the lowest mitochondrial fractions.
- The cause has not been tested. Candidates: injury size, Hr8's origin in a different study, or
  library quality.

**F5. Hr11 (3 dpi, the shallowest library) carries a technical signature.**

- *nppa*'s share of each cell's RNA is 2.53–3.13× the uninjured level in depth-matched bins
  ([check output](../research/cm-state-axes/checks/c002_verify.out)). The excess is largest in
  shallow cells.
- Non-CM genes and *hsp70l* are raised in Hr11 ([results/c003_ambient_by_bin.csv](../results/c003_ambient_by_bin.csv)).
- The *nppa* excess does not behave like fixed-amount ambient RNA, so whether it is biology or
  artefact is inconclusive.
- Teammate PR #14 independently reports the highest ambient load in Hr11 (median 18%).

**F6. Co-occurrence cannot be resolved with literature sets.**

- Overlapping sets co-occur by construction.
- Disjoint sets look like random pairs: all 6 pairs "mixed"
  ([results/c004_g4.json](../results/c004_g4.json)).
- This is not evidence that the states are independent.

**F7. Agent labels for the human-run agreement check:**
[results/c005_agent_labels.csv](../results/c005_agent_labels.csv).

- States are scored only in ventricle-called cells; the others are `not_scored`.
- State calls are at or below chance. `agent_chamber` is the only validated column.
- Read [note-for-evaluate.md](../research/cm-state-axes/cycles/C005/note-for-evaluate.md) before
  running `src/zf/evaluate.py`.

**Traps caught, and by whom** (the full list is in
[checks-inventory.md](../research/cm-state-axes/checks-inventory.md)):

- **The researcher, unprompted:** the false comfort of method agreement, the Hr11 *nppa* excess,
  and the timing mismatch of the injury markers.
- **The supervisor:** the embedded answer key, the uncounted mitochondrial contig, the missing
  positive control, the chamber confound (F2), and atrial cells leaking into the labels.
- ***ttn.2*:** any result about it is **hinted**, not found unaided. An earlier repo document
  named an authors' group "ttn.2 (dedifferentiated)" (FU-6).

## 4. Where it got stuck

- **No validated positive control.** The sensitivity of the time contrast is unknown.
- **Doublets unresolved.** Scrublet's automatic threshold called about 0 doublets in every library
  (FU-7), while PR #14 reports 96% in Hr11 with the same tool. A calibrated threshold or a second
  method is needed.
- **No ambient profile from empty droplets.** It needs the raw (unfiltered) matrices.
- **Reading gaps (FU-1).** Wu 2016 full text, the second half of Honkoop 2019, Ben-Yair 2019,
  Sallin 2015, Beisaw 2020 and a sourced zebrafish cell-cycle set.
- **Errors the loop made and then corrected:**
  - a *vmhc* symbol doubt;
  - atrial cells in the first label file;
  - an overstated "disappears at high depth";
  - supervisor timestamps that ran ahead of the clock.
  All are corrected in the records.

## 5. What it would need next

- **Direction C:** data-driven programs in ventricular CMs with count splitting (Neufeld 2024),
  with per-heart covariates.
- A validated positive-control dataset: Honkoop 2019 *nppa*+ sorted CMs, GSE139218.
- Raw droplet matrices for an ambient-RNA profile, and a second doublet method.
- An expert review of the traps and of the chamber question. The questions for Eme: is the
  chamber confound known, and is per-chamber scoring the right unit?

## 6. How to run

```shell
uv sync
scripts/get_data.sh                       # public GEO matrices (doublet step only)
uv run python -m zf.io <labelled.h5ad> data/clean/hu2022_cm_allgenes_sanitized.h5ad
uv run pytest
uv run python scripts/data_facts.py data/clean/hu2022_cm_allgenes_sanitized.h5ad markers results
uv run python scripts/c002_probe.py data/clean/hu2022_cm_allgenes_sanitized.h5ad markers results 100 100
# further steps: scripts/c00*_*.py; the exact commands are in each results/*.json sidecar
```

Team note: the agreed agent input is `data/processed/cardiomyocytes.h5ad` from
`scripts/make_agent_input.py` (#11, #13). `sanitize_h5ad` predates that and should be replaced by
it (see the PR).

## 7. Data and sources

- Hu et al. 2022, *Nat Genet*, doi:10.1038/s41588-022-01129-5 (GEO GSE159032, GSE158919). Public.
- Li et al. 2025, *Nat Commun*, doi:10.1038/s41467-025-59070-0, Supplementary Data 2 and 5.
- Honkoop et al. 2019, *eLife*, doi:10.7554/eLife.50163.
- Wu et al. 2016 (search summary only).
- ZFIN, NCBI and OpenAlex.

## 8. How Claude was used

- Claude Code with two agents: a Research Supervisor and a Researcher.
- They followed a shared research contract and their own skills (`research-supervision`,
  `evidence-research`).
- The agents read the literature (web, Europe PMC, OpenAlex), wrote code test-first with uv,
  pytest and ruff, ran the analyses, and reviewed each other's claims.

## 9. Status and limits

- **Not complete.** Seven cycles, no direction decision and no Deepen. The supervisor's method
  review is its own, not independent.
- **Manuscript readiness: not assessed.** No manuscript was requested.
- **Domain and statistical review by humans is still needed**, especially for F2 to F4.
