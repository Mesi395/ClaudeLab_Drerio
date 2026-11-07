# Research report: cardiomyocyte state axes after cryoinjury (progress so far)

Study `cm-state-axes`, 2026-11-07. It was run with the method in
[research-methodology.md](research-methodology.md). For a one-page overview, see
[research-summary.md](research-summary.md). Status: **not complete**. Seven cycles were
run and reviewed, and no Deepen phase was reached. Every number is quoted from a file in
`results/` or `research/cm-state-axes/checks/`, linked beside it. The full record is in the
[research map](../research/cm-state-axes/research-map.md).

## 1. Goal

After cryoinjury, zebrafish cardiomyocytes are said to take on four states: injury-responsive,
border-zone-like, dedifferentiated and proliferating. With no ground truth, the goal is to:

- find sourced marker genes for each state (G1);
- score every cell on each state axis (G2);
- test the scores against library size, doublets and ambient RNA (G3);
- test whether the axes co-occur (G4);
- show which checks the agent ran on its own, and hand over labels for a human-run agreement
  check (G5).

Data: Hu et al. 2022, 32,044 cardiomyocytes × 18,794 genes, 22 libraries, which are 18 hearts
(uninjured, 3, 7 and 30 dpi).

## 2. Path of the study

Each row is one cycle: what the supervisor asked, what came back, what the supervisor's own check
showed, and what was decided.

| Cycle | Kind | Question | Result | Supervisor check | Decision |
|---|---|---|---|---|---|
| C001 | Understand | What markers does existing work offer, and what does the data hold? | Labels stripped by a tested sanitiser. The literature gives overlapping sets (*nppa*, *nppb* and *mustn1b* on three states). Hr11 has the lowest median depth. | Sanitiser confirmed. Found that `percent.mito` omits the mitochondrial contig, and a doubt about a gene symbol (SUP-E5). | Probe all four axes |
| C002 | Probe | Do the sets beat matched nulls and depth? Is rank scoring less tied to depth? | No set beats null; rank scores are more tied to depth. | Reproduced with separate code, and extended to an AUCell-style score. Found that there was **no positive control** (SUP-E6). | Probe a positive control and the Hr11 lead |
| C003 | Probe | Do scores rise after injury? Is the Hr11 *nppa* excess ambient RNA? | No positive-control pass. Hr11 has a technical signature; the cause of the *nppa* excess is inconclusive. | Found the **chamber confound**: the sets mark atrium (SUP-E7). | Validate a chamber call; rerun in ventricle only |
| C004 | Probe | Chamber call; ventricle-only control; co-occurrence; labels | Chamber call validated. Weak passes. Co-occurrence unresolved. Labels written. | Found that the **labels inherit the chamber trap** (SUP-E8). | Fix the labels; doublets |
| C005 | Fix and Probe | Ventricle-only labels; doublets | Labels fixed. Doublet threshold failed (about 0 called). | Fix confirmed; threshold failure confirmed (SUP-E9). | No more probes of the same sets; read to choose a direction |
| C006 | Understand | Inventory of checks; close the reading gaps; compare directions | [checks-inventory.md](../research/cm-state-axes/checks-inventory.md); Li 2025 SD2 read; directions A, B and C drafted | Chamber genes now source-backed (SUP-E10) | Probe direction A |
| C007 | Probe | Do sets from Li 2025 SD2 beat the C004 sets? | No: 0 of 6 contrasts pass for either family. The C004 passes vanish under a calibrated null. | Cell-cycle signal confirmed only in hearts Hr9 and Hr8 (SUP-E11) | A fails as built; C is the next candidate |

## 3. Findings

**F1. Rank-based scores are more tied to sequencing depth than `score_genes`.**

- Border-zone axis, median within-library |rho| with log10 counts
  ([checks/c002_verify.out](../research/cm-state-axes/checks/c002_verify.out)):
  score_genes 0.287, UCell 0.468, AUCell-style 0.427.
- Both rank scores are more tied in 22 of 22 libraries.
- All axes: [results/c002_axis_summary.csv](../results/c002_axis_summary.csv). Thinning removes
  most of UCell's tie, so the tie is technical
  ([results/c002_thinning_summary.csv](../results/c002_thinning_summary.csv)).
- The user's hypothesis that rank scores are independent of library size is refuted for these
  data.
- Confidence: high. The AUCell package itself was not run.

**F2. The published border-zone and injury marker sets mainly mark the atrium.**

- Median border-zone score, atrium-only against ventricle-only library of the same heart
  ([checks/c003_verify.out](../research/cm-state-axes/checks/c003_verify.out)): H8 0.679 vs
  0.068; Hr6 0.359 vs 0.058; Hr7 0.345 vs −0.004.
- The atrial fraction per heart ranges from 0.043 to 0.827.
- The chamber call (`myh6` atrial, `myh7l` ventricular, both source-backed by Li 2025 SD2,
  SUP-E10) reaches balanced accuracy 0.986 at coverage 0.701 on the six chamber-only libraries
  ([checks/c004_verify.out](../research/cm-state-axes/checks/c004_verify.out)). That figure is
  optimistic, because the rule was chosen there.
- Confidence: high for the within-heart contrast.

**F3. With depth and chamber controlled, no sourced set shows the predicted time signal above
random genes.**

- Ventricle-only cells, thinned, 18 hearts, a threshold calibrated on control hearts, and a pass
  requiring an AUC above the matched-null 95th percentile: 0 of 6 contrasts pass for either the
  C004 or the Li SD2 sets ([results/c007_positive_control.csv](../results/c007_positive_control.csv)).
- Confidence: moderate. There are 4 control hearts, and the test's sensitivity is not validated.

**F4. The remaining signal sits in a few hearts, not in timepoints.**

- Share of ventricle cells with two or more Li SD2 cell-cycle genes detected
  ([checks/c007_verify.out](../research/cm-state-axes/checks/c007_verify.out)): Hr9 0.0296,
  Hr8 0.0104. The other 7 dpi hearts are at most 0.0037 and the controls at most 0.0023.
- The cause has not been tested: injury size, Hr8's origin in a different study, or library
  quality.

**F5. Hr11 (3 dpi, the shallowest library) carries a technical signature.**

- *nppa*'s share of each cell's RNA is 2.53–3.13× the uninjured level in depth-matched bins
  ([checks/c002_verify.out](../research/cm-state-axes/checks/c002_verify.out)).
- Non-CM genes and *hsp70l* are raised
  ([results/c003_ambient_by_bin.csv](../results/c003_ambient_by_bin.csv)).
- Whether the *nppa* excess is biology or artefact is inconclusive. Teammate PR #14 independently
  reports the highest ambient load in Hr11.

**F6. Co-occurrence (G4) cannot be resolved with literature sets.**

- Overlapping sets co-occur by construction.
- Disjoint sets look like random pairs ([results/c004_g4.json](../results/c004_g4.json)).
- This is not evidence that the states are independent.

**F7. Labels for the human-run agreement check:**
[results/c005_agent_labels.csv](../results/c005_agent_labels.csv).

- States are scored in ventricle-called cells only. State calls are at or below chance.
- `agent_chamber` is the validated column.
- Read [note-for-evaluate.md](../research/cm-state-axes/cycles/C005/note-for-evaluate.md) first.

## 4. Progress per goal part

| Part | Status | What is left |
|---|---|---|
| G1 Markers | Sourced candidate sets with depth flags (`proposals/loop/*.yaml`); Li SD2 sets (`proposals/loop/li2025_sd2_sets.tsv`). Moved from `markers/` for Emi's review | Better-founded axes; reading gaps (FU-1) |
| G2 Scores | Done for all cells with two methods and a matched null | – |
| G3 Confounds | Depth: done. Chamber: done (found). Ambient: partly; Hr11 inconclusive. Doublets: continuous scores show no tracking; threshold failed. | An empty-droplet ambient profile; a calibrated or second doublet method (FU-7) |
| G4 Co-occurrence | An honest "cannot resolve with these sets" | Axes that pass a positive control first |
| G5 Trust | Inventory of checks done; label file delivered; no exposure to the labels | The human-run agreement check |

## 5. Which checks the agent ran on its own

From [checks-inventory.md](../research/cm-state-axes/checks-inventory.md):

- **The researcher, unprompted:** the false comfort of method agreement, the Hr11 *nppa* excess,
  the timing mismatch of the injury markers, the chance-level reading of the labels.
- **The supervisor:** the embedded answer key, the mitochondrial contig, the missing positive
  control, the chamber confound, the atrial labels.
- **An expert would add:** ambient RNA from empty droplets, a second doublet method, the
  origin of the morphine-study libraries, a dissociation-stress axis, a validated positive
  control, spatial validation.
- *ttn.2* is named in the source paper (Hu et al. 2022), so any *ttn.2* result is not found unaided (FU-6).

## 6. Where it got stuck, and what comes next

- **Stuck:**
  - no validated positive control;
  - doublets unresolved (about 0% here against 96% of Hr11 in PR #14 with the same tool);
  - no raw droplets for ambient RNA;
  - reading gaps.
- **Errors made and corrected:**
  - a symbol doubt;
  - atrial cells in the first labels;
  - an overstated depth claim;
  - supervisor timestamps ahead of the clock.
- **Next:**
  - direction C: data-driven programs in ventricular CMs with count splitting (Neufeld 2024) and
    per-heart covariates;
  - a validated positive control (Honkoop 2019 *nppa*+ CMs, GSE139218);
  - an expert decision on chamber as the unit of scoring.
- Manuscript readiness: **not assessed**; no manuscript was requested.

## 7. How to run

```shell
uv sync
scripts/get_data.sh                       # public GEO matrices (doublet step only)
uv run python -m zf.io <labelled.h5ad> data/clean/hu2022_cm_allgenes_sanitized.h5ad
uv run pytest
uv run python scripts/data_facts.py data/clean/hu2022_cm_allgenes_sanitized.h5ad proposals/loop results
uv run python scripts/c002_probe.py data/clean/hu2022_cm_allgenes_sanitized.h5ad proposals/loop results 100 100
# further steps: scripts/c00*_*.py; the exact commands are in each results/*.json sidecar
```

The team's agreed agent input is `data/processed/cardiomyocytes.h5ad` from
`scripts/make_agent_input.py` (#11, #13). `sanitize_h5ad` predates it.

## Sources

- Hu et al. 2022, *Nat Genet*, doi:10.1038/s41588-022-01129-5 (GEO GSE159032, GSE158919)
- Li et al. 2025, *Nat Commun*, doi:10.1038/s41467-025-59070-0
- Honkoop et al. 2019, *eLife*, doi:10.7554/eLife.50163
- ZFIN, NCBI, OpenAlex
