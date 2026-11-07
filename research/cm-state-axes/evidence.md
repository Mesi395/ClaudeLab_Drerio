# Evidence log

Append-only. Each author uses its own IDs: SUP-E* for the supervisor, R-E* for the researcher.

## SUP-E1: the supplied .h5ad files embed the authors' labels (2026-11-07 14:10 CET)

- Source: `data/hu2022/Hu2022_CM_Ctrl_CMonly_processed_HVG.h5ad` and `..._allgenes.h5ad`, moved
  from ~/Downloads (MD5 3e0696df… and 2911317b…).
- Method: `anndata.read_h5ad(backed="r")`, printing only shape, obs column names, dtypes, the
  number of unique values per column, layer, obsm and uns keys, and whether X and the `counts`
  layer hold integers. No label values were printed.
- Result:
  - Both files are 32,044 cells; HVG has 2,000 genes and allgenes has 18,794.
  - obs columns: 7 technical or design columns (`orig.ident`, `nCount_RNA`, `nFeature_RNA`,
    `percent.mito`, `time`, `morphine`, `AV`) and 11 author annotation or clustering columns.
    The names of the 11 are withheld in this public repository.
  - `layers['counts']` is integer. X is normalised: scaled in HVG, log1p in allgenes.
  - obsm holds `X_pca`, `X_pca_harmony` and `X_umap`.
- Reading: the 11 annotation columns are close to the held-out annotation and must never reach an
  agent.
- Limits: the label meanings were inferred from column names only.

## SUP-E2: design columns are constant within each library (2026-11-07 14:11 CET)

- Method: for `time`, `AV`, `morphine`, one constant annotation column and the authors'
  clustering column, the maximum number of distinct values within one `orig.ident`.
- Result:
  - `time`, `AV` and `morphine` are constant per library, so they are design facts.
  - The clustering column varies per cell.
  - Per-library cell counts match the slide table in `docs/problem.md` exactly; the timepoint
    totals are Ctrl 7,048, 3dpi 6,185, 7dpi 14,094 and 30dpi 4,717.
- Consequence: the 564-cell filter question in `docs/data.md` is moot, because these files
  already hold the slide's 32,044 cells.

## SUP-E3: the zip duplicates the .h5ad files (2026-11-07 14:10 CET)

- Method: `unzip -Z1` with a name filter for label-like words. The filter hid 0 of 2 entries.
- Result: `20261107_ClaudeLabDay.zip` contains only the same two `.h5ad` file names. It was not
  extracted.

## SUP-E4: Li et al. 2025, partial read (2026-11-07 14:11 CET)

- Source: Li L. et al., "An organ-wide spatiotemporal transcriptomic and cellular atlas of the
  regenerating zebrafish heart", *Nat Commun* 16 (2025), doi:10.1038/s41467-025-59070-0,
  PMC12009352.
- Method: WebFetch summary of the PMC text, characters 0–100,000 only. nature.com redirected to
  login.
- Result, as the source claims it:
  - Injury: ventricular apex **resection**. Cryoinjury was used only for supporting experiments.
  - Stages: 6 hpa, 12 hpa and 1, 3, 7, 14 and 28 dpa.
  - Platforms: Stereo-seq spots and scRNA-seq.
  - Five stage-specific border-zone myocardial domains, with these named genes:
    - activated: *tnfrsf11b, hspb6, hspb11, nppa, nppb*
    - dedifferentiated: *nppa, nppb*
    - dedifferentiated-proliferating: *nppb, mustn1b, desma, pcna*
    - proliferating: *pcna, cdk1, ccng1, ccng2*
    - re-differentiating: *cox6a2, cox5b2, cox4i1l*
  - *tpm4a* is proposed in re-differentiation.
  - Full marker lists are in Supplementary Data 2 and 5 (not opened).
  - The paper reports correlation with Hu et al. clusters (Supplementary Fig. 3c, 14b).
- Limits: this is a summary by a fetch model, not a full-text read. The second half of the
  paper, the Data Availability section and the supplements were not read. Spot domains are
  multicellular, and the injury model differs from this dataset's cryoinjury.

## SUP-E5: supervisor checks of C001 (2026-11-07 15:28 CET)

- Command: `uv run python -I research/cm-state-axes/checks/c001_verify.py data/clean/hu2022_cm_allgenes_sanitized.h5ad`
- Output: `research/cm-state-axes/checks/c001_verify.out`. Commit 371e96d with a dirty tree.
- Results:
  - The sanitised file's obs holds only the 7 allowlisted columns; uns, obsm and obsp are empty.
  - Hr11 has the lowest median depth (659), but four libraries are within 120.
  - Median share of each cell's UMIs: *nppa* 0.00301 in Hr11 against 0.00075 in H7; *ttn.2*
    0.00140 against 0.
  - `myh7l` is detected in 0.4772 of CMs, `myh6` in 0.3168 and `myh7` in 0.0340.
  - The top gene by counts is `NC-002333.4` (0.0582 of UMIs); `hbba1.1` is in the top 15.
- Limits: the contig's identity is inferred from the RefSeq accession pattern, and the
  `myh7l` = ventricular myosin reading is inferred, not sourced.

## SUP-E6: supervisor checks of C002 (2026-11-07 15:35 CET)

- Command: `uv run python -I research/cm-state-axes/checks/c002_verify.py data/clean/hu2022_cm_allgenes_sanitized.h5ad`
  (10 s). Output: `research/cm-state-axes/checks/c002_verify.out`. Commit fea5089 with a dirty tree.
- Method: independent of the researcher's code. scanpy `score_genes` on log1p CP10k; own UCell
  (max rank 1500); own AUCell-style AUC (top 5% of genes, zero-count genes given random ranks);
  50 random sets matched on 25 mean-expression bins; nppa share binned on Hr11 depth quintiles.
- Results, border axis:
  - Median within-library |rho| with log10 counts: score_genes 0.287, UCell 0.468, AUCell-style
    0.427. Both rank scores exceed score_genes in 22/22 libraries.
  - P90 above the null's 95th percentile in 0/22 libraries.
- Results, *nppa* share in Hr11 against uninjured cells by depth quintile: ratios 2.53, 2.81, 3.10,
  3.13 and 1.53.
- Limits:
  - Border axis only for P1 and P4.
  - The AUCell tie handling is assumed, and the AUCell package was not run.
  - One seed.

## SUP-E7: supervisor check of C003, chamber confound (2026-11-07 15:46 CET)

- Command: `uv run python -I research/cm-state-axes/checks/c003_verify.py data/clean/hu2022_cm_allgenes_sanitized.h5ad results/c003_heart_scores.csv`.
  Output: `research/cm-state-axes/checks/c003_verify.out`.
- Results:
  - Border-zone score_genes median, atrium vs ventricle library of the same heart: H8 0.679 vs
    0.068; Hr6 0.359 vs 0.058; Hr7 0.345 vs −0.004.
  - Injury-responsive: 0.876 vs 0.404, 0.313 vs 0.189, 0.371 vs 0.167.
  - A crude atrial call (myh6 share > myh7l share) labels 0.661 of atrium-library cells and 0.021
    of ventricle-library cells atrial. The atrial fraction per heart ranges from 0.043 to 0.827.
  - Spearman of the per-heart score with atrial fraction: border 0.26, injury-responsive 0.42.
  - Border axis in ventricle-called cells: the top four hearts are injured (Hr9, Hr8, Hr11,
    Hr12).
- Limits:
  - The chamber call is my assumption and is not validated.
  - The ventricle-only contrast is exploratory and unthinned.
  - That myh7l is ventricular is inferred.

## SUP-E8: supervisor checks of C004 (2026-11-07 15:47 CET)

- Command: `uv run python -I research/cm-state-axes/checks/c004_verify.py data/clean/hu2022_cm_allgenes_sanitized.h5ad results/c004_agent_labels.csv`.
  Output: `research/cm-state-axes/checks/c004_verify.out`.
- Results:
  - The label file has 32,044 unique cells, equal to the sanitised obs.
  - Chamber call against AV in the six chamber-only libraries: atrium 1.000, ventricle 0.973
    among called cells; balanced 0.986; coverage 0.701.
  - `injury_responsive` calls: 885 atrial, 184 ventricular, 119 chamber-unassigned. Share by
    time: Ctrl 0.050, 3 dpi 0.027, 7 dpi 0.022, 30 dpi 0.074.
- Limits: the chamber validation uses the same six libraries the rule was chosen on.

## SUP-E9: supervisor checks of C005 (2026-11-07 15:59 CET)

- Command: `uv run python -I research/cm-state-axes/checks/c005_verify.py results/c005_agent_labels.csv results/c005_doublets_by_library.csv`.
  Output: `research/cm-state-axes/checks/c005_verify.out`.
- Results:
  - Every state call is in ventricle-called cells; the 8,816 atrial and 9,204
    chamber-unassigned cells are `not_scored`.
  - State calls per heart: Hr9 149, Hr11 69, Hr8 44.
  - Scrublet all-cell doublet rate 0.000–0.0004 at thresholds of 0.45–0.67 in libraries of up to
    16,587 cells, which is a threshold failure.
- Limits: the doublet scores themselves were not recomputed by the supervisor.

## SUP-E10: supervisor read of Li 2025 SD2 chamber rows (2026-11-07 16:06 CET)

- Command: `uv run python -I scripts/xlsx_peek.py research/cm-state-axes/sources/li2025_supp/MOESM5.xlsx '^(myh7l|myh6)$'`
- Results, sheet "scRNA-seq chamber":
  - `myh6` marks Cardiomyocytes (A): avg_log2FC 9.23, pct.1 0.969 against pct.2 0.054.
  - `myh7l` marks Cardiomyocytes (V): avg_log2FC 6.79, pct.1 0.996 against pct.2 0.290.
- Results, sheet "scRNA-seq regeneration": the same, and `myh7l` also marks Cardiomyocytes
  (dediff.).
- Limits: the MOESM5 = SD2 mapping is inferred from the file order and the sheet contents (R-E19).
  Li 2025 used resection and a newer annotation.

## SUP-E11: supervisor check of C007, cell-cycle signal per heart (2026-11-07 16:10 CET)

- Command: `uv run python -I research/cm-state-axes/checks/c007_verify.py data/clean/hu2022_cm_allgenes_sanitized.h5ad markers/li2025_sd2_sets.tsv results/c005_agent_labels.csv results/c001_library_qc.csv`.
  Output: `research/cm-state-axes/checks/c007_verify.out`.
- Results: share of ventricle cells with ≥2 genes of the Li SD2 proliferating set detected:
  Hr9 0.0296, Hr8 0.0104, the other 7 dpi hearts ≤0.0037, controls ≤0.0023. After thinning: Hr9
  0.0068, Hr8 0.0044.
- Limits: a detection count, not expression; one seed for thinning.

## SUP-E12: what the axes track within ventricle cells (2026-11-07 16:37 CET)

- Prompted by the user's question: what do the markers mark inside the ventricle?
- Command: `uv run python -I research/cm-state-axes/checks/c008_verify.py data/clean/hu2022_cm_allgenes_sanitized.h5ad results/c005_agent_labels.csv`.
  Output: `research/cm-state-axes/checks/c008_verify.out`.
- Scope: 14,024 ventricle-called cells; 19 libraries with at least 50 of them.
- Results, median within-library Spearman rho:
  - border, injury and dedifferentiated score_genes against the *nppa* share alone: 0.521, 0.577,
    0.524;
  - the same three against log10 depth: 0.215, 0.227, 0.232;
  - the same three against the mito share including NC: −0.130, −0.131, −0.150;
  - proliferating against depth, mito and *nppa*: 0.036, −0.037, −0.019.
- Results, border UCell against depth: 0.437, larger than score_genes in 19 of 19 libraries.
- Results, heart explains 0.261 of the cell-level variance of the border score.
- Reading:
  - Inside the ventricle, the three axes behave largely as one *nppa*-led axis, still partly tied
    to depth.
  - F1 holds within ventricle cells.
- Limits: unthinned scores; the crude per-cell heart grouping; no null comparison.

## SUP-E13: supervisor check of C008 (2026-11-07 16:42 CET)

- Command: `uv run python -I research/cm-state-axes/checks/c008b_verify.py data/clean/hu2022_cm_allgenes_sanitized.h5ad results/c005_agent_labels.csv`.
  Output: `research/cm-state-axes/checks/c008b_verify.out`.
- Results, ventricle cells:
  - Disjoint sets without *nppa* and *nppb*: pairwise median rho 0.085, 0.083, 0.121 over 19
    libraries.
  - Hr11 against control ventricle cells by Hr11 depth quintile: *nppa* excess of +1.42 to +2.02
    molecules per cell, at ratios of 4.37–7.62×; *myh6* 3.93–21.28×; *hbba1.1* 0.55–1.24×.
- Limits: the chamber call itself uses *myh6*; no empty droplets; one seed.

## R-E1: leakage guard and sanitised copy (2026-11-07 15:00 CET, researcher)

- Code: `src/zf/io.py::sanitize_h5ad`, test `tests/test_io_sanitize.py` (written first; failed on
  import, then passed: 3 passed). The function opens the source with h5py and reads only the
  allowlisted obs columns, the obs index, `var` and `layers/counts`; dropped columns are never read.
- Command: `uv run python -m zf.io data/hu2022/Hu2022_CM_Ctrl_CMonly_processed_allgenes.h5ad data/clean/hu2022_cm_allgenes_sanitized.h5ad`
- Input MD5 2911317bada824d1575d2476ed3322c7; output MD5 56c3c8b2ba8eb57a4ce2939d72174526;
  provenance sidecar `data/clean/hu2022_cm_allgenes_sanitized.h5ad.json` (names of dropped
  elements only). X in the output = raw integer counts.
- Check: the output's HDF5 keys are obs = the 7 allowlisted columns + index; uns, obsm empty;
  layers = counts only. 32,044 cells × 18,794 genes.
- No label value was printed or seen. The HVG file was not sanitised (not needed).
- Limits: `var` keeps the authors' HVG flags and dispersion stats (allowed as "var stats").
  Commit 371e96d with a dirty tree (nothing may be committed in auto mode).

## R-E2: Honkoop et al. 2019 eLife 8:e50163, doi:10.7554/eLife.50163 (accessed 2026-11-07)

- Method: WebFetch of elifesciences.org/articles/50163, characters 0–100,000 of 225,826 (the
  second half, methods and supplements, were not read).
- As the source claims: cryoinjury, 7 dpi, TgBAC(nppa:mCitrine) sorted CMs, SORT-seq, 352 cells.
  Border-zone cluster: *nppa*, *vmhc*, *mustn1b*; embryonic-like: *nkx2.5*, *tbx20*, *tnnt2a*;
  glycolysis up: *hk1*, *pkma*, *pdk2a* (ISH), *ldha* (qPCR), *slc2a1a/b*; mitochondrial-encoded and
  OXPHOS transcripts down; proliferation by PCNA/pH3 protein only.

## R-E3: Li et al. 2025, second read via Europe PMC full-text XML (2026-11-07 15:01 CET)

- Method: WebFetch of `europepmc/webservices/rest/PMC12009352/fullTextXML`, characters 0–100,000
  of 284,850. Data availability statement and the captions of Supplementary Data 2 and 5 are
  beyond what was read; the supplements themselves were not opened.
- Confirms SUP-E4's domain genes, all from **Stereo-seq spots**, resection. Adds stages: activated
  6 hpa; dedifferentiated 12 hpa/1 dpa; dediff-proliferating 3 dpa; proliferating 7 dpa;
  re-differentiating 14 dpa. Border-zone genes up from 12 hpa: *hand2, tbx20, mustn1b, desma,
  mef2cb, cacybp*. Only *tnfrsf11b* was checked in cryoinjury (6 h post-cryoinjury, Suppl. Fig. 7b).
- Consequence: the "activated" genes were defined at 6 hpa; our data's earliest timepoint is 3 dpi.

## R-E4: abstracts from the Zotero reading list (2026-11-07 14:58 CET)

- Source: `data/zotero/CM_zotero_papers.csv` (535 rows, 384 unique titles; the four CSVs are
  identical, SUP/Journal). Filter saved to `sources/zotero_zf_heart_regen_states.csv` (82 papers:
  zebrafish + heart + injury/regeneration + state terms).
- Abstract-level claims used: Boskovic 2026 (10.1016/j.cdev.2026.204082) *ankrd1a* reporter in
  border CMs from 15 hpi, cryoinjury; Meena 2026 (10.1016/j.ydbio.2026.10.004) *mustn1b* reporter
  in dedifferentiating CMs; Zuppo 2023 (10.1242/dev.201163) *foxm1* up in border-zone CMs, target
  *cenpf*; Constanty 2025 (10.1038/s41467-025-59169-4) *mmp14b* in several border-zone cell types;
  Ben-Yair 2019 (10.1242/dev.178632) sarcomere/structural genes silenced in proliferating CMs;
  Sallin 2015 (10.1016/j.ydbio.2014.12.002) embryonic ventricular myosin protein in the
  injury-abutting zone, cryoinjury; Forman-Rubinsky 2025 (10.1242/dev.204659) snRNA-seq CM
  populations incl. stress response, *cited4a* in mature CMs; Kikuchi 2010 (10.1038/nature08804)
  *gata4* — title only.
- Limits: abstracts only; injury model/timepoint not checked where the YAML says so.

## R-E5: Wu et al. 2016, doi:10.1016/j.devcel.2015.12.010 (search summary only)

- Method: WebSearch (ZFIN figure pages) + Zotero abstract. Claims: two border zones (BZmyo, BZinj);
  *nppb*, *mustn1* up in BZmyo; *crip2*, *ywhabl* in BZinj; BMP targets peak at 3 dpi; cryoinjury,
  tomo-seq tissue sections. Full text and GSE74652 not opened. Not decision-grade.

## R-E6: ttn.2 is named in the source paper

- *ttn.2* is named in the source paper (Hu et al. 2022), so any *ttn.2* result is not found unaided. It is listed in the YAML only as a circular comparison, never as an independent
  marker.

## R-E7: ZFIN: vmhc is a previous name of myh7 (accessed 2026-11-07 15:01 CET)

- https://zfin.org/ZDB-GENE-991123-5: current symbol *myh7*; previous names include *vmhc*.
  The data has `myh7`, `myh7l`, `myh7ba`, `myh7bb`; not `vmhc`.

## R-E8: per-library QC and marker detection (2026-11-07 15:01 CET)

- Command: `uv run python scripts/data_facts.py data/clean/hu2022_cm_allgenes_sanitized.h5ad markers results`
- Outputs: `results/c001_library_qc.csv`, `results/c001_marker_detection.csv`,
  `results/c001_data_facts.json` (commit, dirty-tree flag, command, input). A rerun gave a
  byte-identical detection table.
- Checks: counts are integers; 0 missing covariates; recomputed total counts equal `nCount_RNA`
  exactly in all libraries; 13 `mt-` genes.
- Results (quoted from the files):
  - Hr11 median total counts per CM = 659, the lowest of 22; but Hr7a (689), Hr2b (751), Hr2a (752)
    and H8v (776) are close. In CMs, Hr11's depth gap is smaller than the all-cell figure in
    `docs/data.md` suggests. Hr3 (2,317) is the outlier on the high side.
  - Atrium-only libraries (Hr7a, Hr6a, H8a… see file) have lower mito fraction except H8a.
  - Absent symbols: `slc2a1a`, `vmhc` (→ `myh7`, R-E7).
  - Detection ≤1% of all CMs: *ankrd1a, mmp14b, gata4, foxm1, cenpf, cdk1*; *pcna* 2.3%.
  - Observation, not a conclusion: in Hr11 (lowest depth), *nppa* (0.771 vs 0.561 in H7) and
    *ttn.2* (0.620 vs 0.288) are detected in more cells than in a control library, while most
    other genes are detected less. Either 3 dpi biology or a library-specific artefact (ambient);
    untested.
- Limits: detection fraction is not expression level; no within-library depth control yet.
  `myh7` detected in only 3.4% of CMs, which is implausible for the main ventricular myosin;
  the vmhc reads may map to another `myh7*` symbol — unresolved.

## R-E9: F1, the vmhc symbol (2026-11-07 15:24 CET, researcher)

- Sources: ZFIN ZDB-GENE-991123-5 (myh7; previous names include *vmhc*), R-E7; ZFIN ZDB-GENE-061027-393
  (myh7l, "myosin heavy chain 7-like"; previous names *vmhcl*, zgc:113832, zgc:152732, myh1), via a
  web search summary of the ZFIN pages (page not opened directly).
- Reading: Honkoop's *vmhc* is `myh7` by nomenclature; `myh7l` is a different paralog (*vmhcl*). The
  YAML keeps `myh7` and records the check. `myh7`'s low detection (supervisor SUP-E5) stays
  unexplained: an annotation or mapping issue in the authors' reference is possible but not
  shown. Not settled with Ensembl IDs: `var` holds symbols only. The teammate file
  `li2025_bzm.tsv` was not read (user: optional).

## R-E10: F2, NC-002333.* features (2026-11-07 15:25 CET)

- Source: NCBI E-utilities esummary, nuccore NC_002333: "Danio rerio mitochondrion, complete
  genome", 16,596 bp.
- Data: 19 features named `NC-002333.*` in the sanitised file, plus 13 `mt-` genes. The supplied
  `percent.mito` equals the `mt-` share exactly (median absolute difference 0), so it **excludes**
  the NC features. Their median share per cell is 0.0628 (max 0.65), with Spearman 0.56 against the
  `mt-` share. Command: inline script in this session; not saved as a result file. Recompute before
  use as a covariate: mito_all = mt- + NC-002333.*.

## R-E11: C002 primary Probe (2026-11-07 15:27 CET)

- Command: `uv run python scripts/c002_probe.py data/clean/hu2022_cm_allgenes_sanitized.h5ad markers results 100 100`
  (log `results/c002_probe.log`, provenance `results/c002_probe.json`, 60 s). Code: `src/zf/score.py`
  (UCell re-implemented, test-first: `tests/test_score.py`, matches a brute-force UCell).
  Reading rules fixed at 15:24, before computing: `cycles/C002/reading-rules.md`.
- Reproduced: rerun with the final script into /tmp gave byte-identical summary and per-library
  tables.
- Outputs: `results/c002_axis_summary.csv`, `c002_axis_by_library.csv`, `c002_scores.csv.gz`
  (per-cell scores), `c002_gene_depth.csv` (+ background), `c002_c9_depth_matched.csv`,
  `c002_heart_summary.csv` (18 hearts, from `scripts/c002_hearts_and_flags.py`).
- Results are quoted in the C002 report from these files.

## R-E12: exploratory null variant (2026-11-07 15:29 CET, exploratory)

- `C002_NULL_BINS=100 C002_NULL_EXCLUDE_TECH=1 uv run python scripts/c002_probe.py ... results/c002_exploratory_null 100 100`
- Purpose: test whether a too-wide top expression bin made the null too lenient or strict.
  Outputs in `results/c002_exploratory_null/`.

## R-E13: binomial thinning to 500 UMIs (2026-11-07 15:30 CET)

- `uv run python scripts/c002_thinning.py data/clean/hu2022_cm_allgenes_sanitized.h5ad markers results`
  (seed 0; 46 cells below target were kept unthinned). Outputs: `results/c002_thinning_summary.csv`,
  `c002_thinning_by_library.csv`, `c002_thinning.json`. Rule added at 15:30, before this run.

## R-E14: C003 Probe, positive control and ambient test (2026-11-07 15:37 CET)

- Rules fixed at 15:36, before computing: `cycles/C003/reading-rules.md`.
- Command: `uv run python scripts/c003_probe.py data/clean/hu2022_cm_allgenes_sanitized.h5ad markers results`
  (log `results/c003_probe.log`, provenance `results/c003_probe.json`).
- Outputs: `results/c003_heart_scores.csv` (per heart × axis: thinned median score, null 99th
  percentile, tail fraction), `results/c003_positive_control.csv` (contrasts, AUC, pass, with and
  without Hr11), `results/c003_ambient_by_bin.csv` (per group × depth bin × gene), and
  `results/c003_ambient_verdict.json` (rule components and verdict).
- Departure from the rules: `hbaa1.1` is absent from the data, so the non-CM panel has 12 genes.
  Hr24 and Hr27 have fewer than 20 cells in bin 0, so only bin 1 enters their shallow-bin
  statistics.
- Reproduced: a rerun into /tmp gave byte-identical outputs for all four result files (15:41).

## R-E15: C003 exploratory within-library cell correlations (2026-11-07 15:38 CET, exploratory)

- `uv run python scripts/c003_exploratory_cellcorr.py data/clean/hu2022_cm_allgenes_sanitized.h5ad results`
  gives `results/c003_exploratory_cellcorr.csv`. Added after seeing the C003 results.

## R-E16: C004 chamber call, ventricle positive control, G4, agent labels (2026-11-07 15:45 CET)

- Rules fixed at 15:43, before computing: `cycles/C004/reading-rules.md`.
- Chamber sources: Berdougo et al. 2003, doi:10.1242/dev.00838 (myh6/amhc atrial), and Yelon et al.
  1999, doi:10.1006/dbio.1999.9406 (vmhc/myh7 ventricular). Both are search-level (ZFIN and
  catalogue records), not full text. `myh7l` (formerly vmhcl) has no source for ventricular
  expression.
- Command: `uv run python scripts/c004_probe.py data/clean/hu2022_cm_allgenes_sanitized.h5ad markers results`
  (log `results/c004_probe.log`, provenance and label-file format `results/c004_probe.json`).
- Outputs: `results/c004_chamber_validation.csv`, `c004_agent_labels.csv` (32,044 rows),
  `c004_ventricle_heart_scores.csv`, `c004_ventricle_positive_control.csv`,
  `c004_g4_cooccurrence.csv`, `c004_g4.json`.
- Limits: the chamber rule was selected on the same six libraries it was validated on. The
  ventricle-only positive control uses the same thinning seed as C003, and the null calibration
  differs between hearts.
- Reproduced: a rerun into /tmp gave byte-identical outputs for all six result files.

## R-E17: C005 corrected label file, F10 (2026-11-07 15:49 CET)

- Rules fixed at 15:48, before computing: `cycles/C005/reading-rules.md`.
- Command: `uv run python scripts/c005_labels.py data/clean/hu2022_cm_allgenes_sanitized.h5ad markers results`
  (log `results/c005_labels.log`).
- Outputs: `results/c005_agent_labels.csv` (32,044 rows; format and provenance in
  `results/c005_agent_labels.json`), `results/c005_label_counts_vs_c004.csv`,
  `results/c005_labels_by_heart.csv`. `c004_agent_labels.csv` is unchanged.
- Note for the evaluate.py runner: `cycles/C005/note-for-evaluate.md`.

## R-E18: C005 doublets, R3.2 (2026-11-07 15:50 CET onward)

- Data: `bash scripts/get_data.sh` downloaded the 22 GEO filtered matrices to `data/raw/`
  (exit 0; log /tmp/c005_download.log, not kept in the repo).
- Attempt 1 failed in every library: scanpy's scrublet needs scikit-image for its automatic
  threshold. Kept as `results/c005_doublets_attempt1_failed.log` and `..._by_library.csv`.
  Fix: `uv add scikit-image` (0.26.0; pyproject.toml and uv.lock changed).
- Attempt 2: `uv run python -I scripts/c005_doublets.py data/raw results` (background job; log
  `results/c005_doublets.log`).
- Attempt 2 result (15:58 CET): exit 0, 22/22 libraries ok, mapping rate 1.0. Outputs:
  `results/c005_doublets_by_library.csv`, `c005_doublets_summary.json`, `c005_doublets_cells.csv`.
  The automatic threshold called no doublets (rate 0.000 in all libraries).
- Exploratory: `uv run python scripts/c005_doublets_exploratory.py results` gives
  `results/c005_doublets_exploratory_top10.csv` (+ .json).

## R-E19: C006 forward citations and Li 2025 supplements (2026-11-07 16:02–16:08 CET)

- OpenAlex (`api.openalex.org`, accessed 2026-11-07): Honkoop 2019 = W2995700709 (cited_by 292);
  Li 2025 = W4409591603 (cited_by 36). Citing lists (filtered as stated in map section 7) are in
  `sources/openalex/*.json`. Abstracts were read via OpenAlex for 10.3390/cells15161455 and
  10.7554/elife.110256.
- Li 2025 supplements: Springer static MOESM3–10 (MOESM7 not downloaded), saved to
  `sources/li2025_supp/`. The MOESM3 legend gives SD1–SD8. SD2 = "Markers of cell types and
  domains" (Seurat FindAllMarkers, Wilcoxon) and SD5 = "Markers of regeneration days". Mapping
  them to MOESM5 and MOESM8 is inferred from the file order and matches the sheet contents.
- Read with `uv run python -I scripts/xlsx_peek.py <file> [--clusters | regex]`:
  - SD2 sheets: scRNA-seq chamber (16 clusters), scRNA-seq regeneration (25 clusters, including
    "Cardiomyocytes (dediff.)" with 484 rows and "Proliferating cells" with 477), Stereo-seq
    regeneration (20 domains, including BZm activ., dediff., dediff.-prolif., prolif. and
    re-diff.), Stereo-seq 3D.
  - Chamber rows: `myh6` for Cardiomyocytes (A): pct.1 0.969 vs pct.2 0.054, log2FC 9.23.
    `myh7l` for Cardiomyocytes (V): 0.996 vs 0.290, log2FC 6.79. `myh7` for Cardiomyocytes (V):
    0.765 vs 0.012, log2FC 6.75.
  - Stereo-seq: `myh7` is a marker of BZm (dediff.-prolif.) (log2FC 2.65) and BZm (prolif.)
    (2.90).
- Consequence: the C004 chamber rule B's `myh7l` is now sourced (Li 2025 SD2, scRNA-seq,
  resection model). No new analysis was run.

## R-E20: C007 Probe of direction A (Li 2025 SD2 sets vs the C004 sets) (2026-11-07 16:06–16:10 CET)

- Rules fixed at 16:06, before building sets or computing: `cycles/C007/reading-rules.md`.
  Addendum at 16:08: an exploratory variant, decided after seeing the gene lists but before any
  scores.
- Sets: `uv run python -I scripts/c007_build_sets.py research/cm-state-axes/sources/li2025_supp/MOESM5.xlsx data/clean/hu2022_cm_allgenes_sanitized.h5ad markers/li2025_sd2_sets.tsv`
  gives `markers/li2025_sd2_sets.tsv` (+ .json): 25 genes per axis, disjoint; 3,304 genes dropped
  as non-CM or atrial markers.
- Probe: `uv run python scripts/c007_probe.py data/clean/hu2022_cm_allgenes_sanitized.h5ad markers markers/li2025_sd2_sets.tsv results c007`
  gives `results/c007_positive_control.csv`, `c007_heart_tails.csv`, `c007_probe.json` (log
  `results/c007_probe.log`, exit 0).
- Exploratory: `C007_DROP_CM_IDENTITY=1` build gives `markers/li2025_sd2_sets_exploratory_noCMid.tsv`
  (gene counts: proliferating 25, injury_responsive 8, dedifferentiated 3, border_zone_like 1).
  Probe output: `results/c007_exploratory_noCMid_positive_control.csv`, `..._heart_tails.csv`,
  `..._probe.json` (log in `results/c007_exploratory_noCMid/`, exit 0). No pass.
- Reproduced: after a style-only lint fix to `tails()`, a rerun into /tmp gave byte-identical
  `c007_positive_control.csv` and `c007_heart_tails.csv`.

## R-E21: C008 re-run of the checks within ventricle cells (2026-11-07 16:39 CET onward)

- Rules fixed at 16:38, before computing: `cycles/C008/reading-rules.md`. Marker sets read from
  `proposals/loop/` (moved there by the team's commit 4cc9e9d).
- Command: `uv run python scripts/c008_probe.py data/clean/hu2022_cm_allgenes_sanitized.h5ad proposals/loop results`
  (log `results/c008_probe.log`).
- Outputs: `results/c008_ties_by_library.csv`, `c008_null_by_library.csv`,
  `c008_program_by_library.csv`, `c008_ambient_ventricle.csv`, and `c008_summary.json`
  (provenance in its `meta` field).
- Result: exit 0 at 16:40. 14,024 ventricle cells, 19 libraries, 18 hearts. Reproduced: a rerun
  (after removing one unused variable) gave byte-identical per-library tables.
