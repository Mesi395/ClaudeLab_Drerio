# Research plan and first review (Hemang's lane)

Status: plan, not results.
Every claim below marked *to verify* is a lead, not a finding.

## Framing

This is an open problem with no answer key.
The paper's labels in the processed AnnData
(`pub_CMsubtypes`, `trashid`, `final.zoom`, `lineage.ident`, `work.ident3`, `plot.ident*`)
are one published opinion: Ventricle, Atrium, dediff./ttn.2, proliferating.
They are kept out of marker choice and score tuning,
and compared per state at the end only.

## Data facts (checked by opening the file)

- `Hu2022_CM_Ctrl_CMonly_processed_allgenes.h5ad`: 32,044 cells × 18,794 genes, 22 libraries (`orig.ident`).
- `X` is log-normalised; `layers['counts']` holds raw counts.
- `obs`: `time` (Ctrl, 3dpi, 7dpi, 30dpi), `AV` (Atrium, Ventricle, Wholeheart), `nCount_RNA`, `nFeature_RNA`, `percent.mito`, `leiden`.
- `obsm`: `X_pca`, `X_pca_harmony`, `X_umap`.
- Each library belongs to one timepoint, so library and time are confounded.

## Plan

1. **Markers with citations**, from a local paper corpus (the 535-paper Zotero export, ranked by abstract, top papers read in full) and from spatial data (Wu 2016 tomo-seq, Honkoop 2019 *nppa*+ cells, the 2025 spatial paper).
2. **Four scores per cell**, two methods (`score_genes` and a rank-based method).
3. **Expression-matched decoys**: about 1,000 random gene sets matched on mean expression and detection rate; report the percentile of the real score.
4. **Confound correction, not only detection**: depth residualisation within library (remove only the depth spline), binomial thinning to equal depth and rescoring, doublet and ambient-RNA flags. Never regress on library or time, which would erase biology.
5. **Per-library agreement with high-confidence cut-offs** instead of all-library agreement (see below).
6. **Co-occurrence** of the axes per timepoint and chamber, on corrected scores.
7. **Findings file in two tiers**: observations (unverified) and claims (re-run, survived every check, reviewed).

## First senior review (Claude Fable, to verify)

- Hr11 (3 dpi, about 750 UMIs per cell) holds two-thirds of the 3 dpi cells. "Rises at 3 dpi" and "shallow library" are hard to tell apart. Treat the library as the replication unit; Hr12, Hr24 and Hr27 replicate Hr11.
- Doublet and ambient checks on a cardiomyocyte-only matrix are blind to cardiomyocyte–fibroblast or cardiomyocyte–macrophage doublets. Run them on the full per-library matrices, then map back to cardiomyocyte barcodes.
- Immediate-early and heat-shock genes (fos, jun, egr1, atf3, hsp) are induced by tissue dissociation. Use a dissociation gene set as a negative-control axis before calling a cell injury-responsive.
- Marker sets that share genes, plus a shared dependence on depth, can produce co-occurrence on their own. Use one gene per axis and co-occurrence on corrected scores, against two nulls.
- "Loss of maturity" genes (*ttn.2*, sarcomere, *cox7a1*, *ckma*) lose detection with depth. Keep them secondary to "gain" genes for dedifferentiation, and gate them by the thinning test.
- Proposed high-confidence structure: per-cell empirical p against matched decoys (FDR within library); per-library effect; random-effects meta-regression across libraries (effect ~ timepoint + chemistry, report τ² and I²); Benjamini–Heller partial conjunction (for example 3 of 4 libraries at 3 dpi, 6 of 10 at 7 dpi); leave-one-library-out.

## Open questions for Emi

1. Is there a per-cell expert annotation besides the paper's labels? If so, keep it until the end; tell us its form only.
2. Which border-zone reference do you trust most: Wu 2016, Honkoop 2019, or the 2025 spatial data (amputation, not cryoinjury)?
3. Score within ventricle and atrium separately, or across all cardiomyocytes?
