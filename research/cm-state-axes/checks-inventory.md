# Checks inventory (R5.1): what the team checked, who started it, what it caught

Study `cm-state-axes`, cycles C001–C005. Written by the researcher on 2026-11-07 at 16:05 CET.
Evidence IDs point to [evidence.md](evidence.md): R-E* are the researcher's, SUP-E* the
supervisor's. "Who initiated" means who first asked for the check, not who ran it. All checks
were run by the researcher unless marked otherwise.

- **Researcher, unprompted:** the researcher chose the check without being asked.
- **Supervisor-prompted:** an assignment or review finding asked for it.
- **User or team input:** the user's request, the slide, `TASKS.md`, `AGENTS.md`, the repo docs or
  a teammate's branch.

## Candid summary

- **Found by the researcher without prompting:**
  - the method-agreement false comfort;
  - the coherence test;
  - the *nppa* excess in Hr11 (its cause is still open);
  - the timing mismatch of the injury markers;
  - the conflicting meaning of *ccng1* and *ccng2*;
  - the atrial-soup and cell-level contamination tests;
  - the chance-level reading of the label counts;
  - byte-identical reruns.
- **Found by the supervisor, not the researcher:**
  - the embedded answer key (SUP-E1);
  - the mitochondrial contig missing from the mito fraction (F2, SUP-E5);
  - the doubt about `myh7` detection (F1, SUP-E5);
  - the missing positive control (F5);
  - the **chamber confound (F8, SUP-E7)**;
  - the **atrial injury-responsive label defect (F10, SUP-E8)**.
  
  The researcher had flagged *vmhc* as "chamber identity" in C001 and saw control heart H8's high
  tails in C003, but tested neither; the supervisor did.
- **From the user or the team:**
  - the two-method plan and the rank-is-depth-free claim (the user);
  - the Hr11 depth trap and chemistry notes (`docs/data.md`);
  - binomial thinning and counting hearts rather than libraries (the team's research plan, relayed
    by the supervisor).
- **Hinted, not found unaided:** *ttn.2*. An earlier `docs/data.md` named an authors' cluster
  "ttn.2 (dedifferentiated)" (FU-6). Every *ttn.2* result below is hinted.

## Inventory

| # | Cycle | What it tests | Who initiated | Result (evidence) | Trap caught |
|---|---|---|---|---|---|
| 1 | C001 | Strip the authors' labels from the `.h5ad` before any use (allowlist, test-first) | Supervisor-prompted (SUP-E1) | Clean copy with 7 obs columns; no label value seen (R-E1; confirmed SUP-E5) | **Embedded answer key** (found by the supervisor) |
| 2 | all | Never read `holdout/` | User or team (`AGENTS.md`) | Not opened in any cycle (JOURNAL entries) | Leakage |
| 3 | C001 | Every marker gene present under its zebrafish symbol | User or team (`TASKS.md` markers-1; R1.2) | 1 gene absent (`slc2a1a`); *vmhc* renamed `myh7` from ZFIN (R-E7, R-E8) | Renamed symbols |
| 4 | C002 | Is `myh7` the right gene for *vmhc*? | Supervisor-prompted (F1, SUP-E5) | ZFIN: *vmhc* = `myh7`; `myh7l` = former *vmhcl*. `myh7` is barely detected; unresolved (R-E9) | A "present" check passing on a near-silent gene (found by the supervisor) |
| 5 | C001 | Per-library depth, genes detected, mito, chemistry | User or team (`TASKS.md` data-2) and supervisor | Hr11 has the lowest median depth per CM, with four libraries close (R-E8) | **Hr11 depth trap**: quantified, and smaller in CMs than the all-cell figures suggested |
| 6 | C001 | Recomputed total counts against the supplied `nCount_RNA` | Researcher, unprompted | Exact match (R-E8) | – |
| 7 | C002 | What is `NC-002333.*`? Does the mito fraction include it? | Supervisor-prompted (F2, SUP-E5) | The zebrafish mitochondrial genome; excluded from `percent.mito`; median 6.3% of UMIs (R-E10) | **Undercounted mito fraction** (found by the supervisor) |
| 8 | C001 | Marker detection per library | Supervisor-prompted (assignment) | Proliferation and injury genes in ≤2% of CMs (R-E8) | Barely detectable axes |
| 9 | C001 | Source model, time and resolution for each marker | User or team (R1.3) | Li 2025 "activated" genes defined at 6 hpa on resection spots; our earliest timepoint is 3 dpi (R-E3) | **Timing and model mismatch** (researcher) |
| 10 | C001 | Conflicting meaning of individual markers | Researcher, unprompted | *ccng1* and *ccng2* flagged as stress or arrest genes (YAML) | Meaning conflict; later explains the C003 proliferating "pass" |
| 11 | C001 | Circular marker kept out of the axis sets | Researcher, from a hint | *ttn.2* flagged as circular, kept out (R-E6) | Circularity (**hinted**, FU-6) |
| 12 | C001→C003 | *nppa* (and *ttn.2*, hinted) higher in the lowest-depth library | Researcher, unprompted (*nppa*); supervisor sharpened it (SUP-E5) | Survives depth matching (R-E11 P7); not explained by atrial or non-CM soup (R-E14, R-E15); the fixed rule rates it inconclusive | The slide's trap candidate. Cause **unresolved** |
| 13 | C002 | Two scoring methods agree | User or team (the user's plan; `TASKS.md` score-1) | Spearman ≥ 0.5 in 22/22 libraries, median 0.91–0.96 (R-E11) | **Agreement gives false comfort**: it is not an independent check (researcher) |
| 14 | C002 | Axis vs expression-matched random gene sets (P90) | User or team (`TASKS.md` score-2), designed by the researcher | 0–3 of 22 libraries above null (R-E11); confirmed SUP-E6 | The axes look like random genes |
| 15 | C002 | Do the axis genes co-vary beyond random sets (coherence)? | Researcher, unprompted | Border-zone 17, injury-responsive 13, dedifferentiated 9, proliferating 4 of 22 libraries (R-E11) | – |
| 16 | C002 | Null robustness: finer bins, no mito/NC/ribosomal genes | Researcher, unprompted (exploratory) | Same reading (R-E12) | – |
| 17 | C002 | Within-library correlation of each score with depth | User or team (R3.1, `TASKS.md` confound-1) | score_genes weak (0.07–0.29); UCell 0.28–0.47 (R-E11) | Score tracks depth |
| 18 | C002 | Is a rank score independent of depth? (the user's AUCell claim) | User or team (the user's claim), framed as a hypothesis by the supervisor | **Refuted**: UCell is more tied to depth in 22/22 libraries; an AUCell-style score too (R-E11; SUP-E6) | **Rank methods are not depth-free** |
| 19 | C002 | Binomial thinning to equal depth, then rescore | User or team (team research plan), relayed by the supervisor | Thinning removes most of UCell's depth tie: it is technical (R-E13) | Depth tie is technical |
| 20 | C002 | Per-gene depth flags against same-bin background | User or team (R3.5) | 10 of 26 genes flagged (R-E11) | Depth-tracking markers (flags are fragile) |
| 21 | C002 | Count hearts, not libraries | User or team (team research plan), relayed by the supervisor | 22 libraries = 18 hearts (R-E11 heart summary) | Pseudo-replication |
| 22 | all | Reading rules fixed and time-stamped before computing | Supervisor-prompted (C002 on) | Rules files in C002–C005; later changes labelled exploratory | Forking paths |
| 23 | C002–C005 | Byte-identical rerun of each result | Researcher, unprompted | Reproduced in C002, C003 and C004 (R-E11, R-E14, R-E16) | – |
| 24 | C003 | Positive control: axis rises at the predicted time, per heart | Supervisor-prompted (F5) | No axis passes on the tail; proliferating passes only on the median (R-E14) | **No positive control** (found by the supervisor) |
| 25 | C003 | Is the proliferating "pass" one gene? | Researcher, unprompted (an inference) | Carried by 3 hearts with the highest *ccng1* detection; leave-one-out not run (C003 report) | Single-gene artefact (suspected) |
| 26 | C003 | Hr11 excess: ambient RNA or 3 dpi biology, against the other 3 dpi hearts | Supervisor-prompted (assignment); rules and panels by the researcher | Inconclusive by rule. Non-CM genes ≥2× in 10 of 12 in Hr11, and raised in all 3 dpi hearts; *hsp70l* about 5× in Hr11 only (R-E14) | **3 dpi non-CM contamination and Hr11 stress** |
| 27 | C003 | Atrial-soup check (*myh6*) | Researcher, unprompted | Not raised in Hr11 (R-E14) | Rules out atrial soup for *nppa* |
| 28 | C003 | Cell-level: *nppa* vs non-CM share | Researcher, unprompted (exploratory) | No positive association (R-E15) | – |
| 29 | C003→C004 | Chamber confound | **Supervisor (F8, SUP-E7)**. The researcher saw H8's high tails in C003 but did not test them | nppa-led scores are much higher in atrium-only libraries of the same heart (SUP-E7) | **Chamber composition confounds heart contrasts** |
| 30 | C004 | Chamber call from sourced genes, validated on 6 chamber-only libraries | Supervisor-prompted | Balanced accuracy 0.991 at coverage 0.735; optimistic, and `myh7l` is unsourced (R-E16; SUP-E8: 0.986 pooled) | Makes chamber controllable |
| 31 | C004 | Positive control in ventricle-called cells | Supervisor-prompted | Two weak passes that rest on control hearts below the null; border-zone 0.786 (R-E16) | Uncalibrated null per heart |
| 32 | C004 | G4 co-occurrence with disjoint sets against a null | User or team (goal G4), supervisor-prompted | All pairs "mixed": cannot be resolved with these sets (R-E16) | **Shared genes would make co-occurrence true by construction** (F7, supervisor) |
| 33 | C004 | Label counts against chance | Researcher, unprompted | State labels at or below chance (C004 report) | Labels look like signal but are noise |
| 34 | C004→C005 | Injury-responsive labels in atrial cells | **Supervisor (F10, SUP-E8)**. The researcher saw the H8/Hr3 concentration and named chamber as the likely cause, but did not fix the rule | 885 of 1,188 calls were atrial (SUP-E8); corrected file 145 calls (R-E17) | **Label file inherited the chamber confound** |
| 35 | C005 | Doublets (scrublet, full GEO matrices) against the axes and states | User or team (R3.2, `TASKS.md` confound-2), supervisor-prompted | Run in 22/22 libraries. The automatic threshold called 0 doublets (failure, FU-7); no axis tracks the score (median \|rho\| ≤ 0.091); no state enriched (exploratory) (R-E18) | Partly met; the doublet share is not answered |

## Checks a domain expert would add (not run, or could not be run)

| Check | Why it matters here | Why not run |
|---|---|---|
| Ambient RNA from empty droplets (SoupX, CellBender, DecontX) | Decides the Hr11 *nppa* question directly; contamination is visible at 3 dpi | Only filtered matrices on GEO; raw droplets not available |
| Doublets with an expected-rate threshold and a second method (scDblFinder), plus CM–non-CM marker co-expression | The scrublet automatic call failed (FU-7) | Time; scDblFinder is R |
| Morphine-study libraries (Hr8, Hr11, GSE158919) as a treatment or batch confound | Hr11 is the trap library; Hr8 is among the hearts that carry the proliferating "pass" | Not tested. `docs/data.md` notes the issue, and the `morphine` column is constant in the data |
| Chemistry (v3: Hr24, Hr27) as a covariate | Hr24 is the outlier for *nppa* (0.11×) | Noted, not tested |
| Dissociation-stress control axis (immediate-early and heat-shock genes) | *hsp70l* marks Hr11; stress genes sit on the injury axis (*hspb6*, *hspb11*) | Not run. A teammate branch holds a decoy and dissociation axis (FU-5), not read |
| A positive control with a validated set (for example the sorted *nppa*+ border-zone cells, GSE139218) | Separates "markers mark nothing" from "test too weak" | Data not downloaded (FU-1) |
| Spatial validation (in situ, tomo-seq or Stereo-seq) of axis-high cells | Border zone is a position, not a transcript state | Out of scope for computation |
| Known proliferation rates (EdU/PCNA share of CMs at 7 dpi) as the expected tail size | Sets what a proliferating axis should find | Literature value not read (FU-1) |
| Leave-one-gene-out and gene-subsample stability per axis | *ccng1* may carry an axis | Time |
| Library-specific clusters, integration check (`TASKS.md` confound-4) | Merged without batch correction | Not run |
| Expert review of the marker sets and the trap list (`holdout/traps.yaml`) | The expert knows which marker is the slide's trap | Held out by design |

## Update after C006 reading (16:10)

- Row 30: `myh7l`'s ventricular role is now sourced. Li 2025 Supplementary Data 2 lists it as a
  "Cardiomyocytes (V)" scRNA-seq marker (R-E19).
- Row 4: Li's Stereo-seq lists `myh7` (*vmhc*) as a border-zone-domain marker, which may explain
  its rarity in our CMs. This was found by reading, not tested.
- Expected proliferating tail: a 2026 resection study reports 1–2% of CMs cycling at 7 dpi (abstract
  only, R-E19). This is the size a proliferating axis should detect, and it is close to the 1%
  null tail used in C003–C005.
