# Research plan and first review

Status: plan, not results.
Every claim below marked *to verify* is a lead, not a finding.

## Framing

This is an open problem with no answer key.
The processed file carries published labels;
they are stripped before any agent sees the data and compared once, at the end.

## Units and confounds (from the processed file and GEO titles)

- 32,044 cardiomyocytes × 18,794 genes, 22 libraries; raw counts are available.
- **The unit of replication is the heart, not the library.**
  Hr2a/Hr2b, H8a/H8v, Hr6a/Hr6v and Hr7a/Hr7v are one heart each,
  which gives 4 uninjured hearts, 4 at 3 dpi, 7 at 7 dpi and 3 at 30 dpi.
- Each heart belongs to one timepoint, so heart and time are confounded.
- Chemistry (10x v3 for Hr24 and Hr27, v2 for the rest) varies only within 3 dpi,
  so it cannot be separated from 3 dpi in a model.
- The depth range inside 3 dpi is the natural test of the depth trap.
  Median UMIs per barcode, all cells: Hr11 713, Hr12 1,696, Hr27 3,424, Hr24 7,247.
- Chamber comes from the GEO titles, never from library suffixes
  (the `a` in Hr2a means "first sample", not atrium).
  Whole-heart libraries have no per-cell chamber.

## Plan

Builds on `config/axes.yaml` and `config/checks.yaml` (#2, #3),
which already hold matched decoys, the dissociation control axis, disjoint marker sets and partial-correlation confound checks.
New parts proposed here:

1. **Markers with citations** from a local paper corpus
   (the Zotero export, deduplicated to unique titles, ranked by abstract; abstracts are not committed)
   and from spatial data (Wu 2016 tomo-seq, Honkoop 2019 *nppa*+ cells, the 2025 spatial paper).
   Each gene in only one axis.
2. **Decoys matched on detection rate as well as mean expression**;
   report the percentile of the real score among about 1,000 decoys.
3. **A `depth_thinning` check**: binomially thin every cell to equal depth and rescore.
   A real state keeps its per-heart pattern; an artefact appears in every heart.
4. **Depth residualisation within heart**: remove only the depth spline, never heart or time.
5. **Doublet and ambient-RNA flags on the full per-library matrices**, mapped back to cardiomyocyte barcodes.
6. **Per-heart high-confidence rule** in the evaluation: per-cell empirical p against decoys (FDR within heart);
   per-heart effect; random-effects meta-analysis across hearts (τ² reported, not gated on, with 3–4 hearts per condition);
   Benjamini–Heller partial conjunction across hearts; leave-one-heart-out.
7. **Co-occurrence** of the axes on corrected scores.
8. **Findings in two tiers**: observations (unverified) and claims (re-run, survived every check, reviewed).

## First senior review (Claude Fable, to verify)

- Hr11 holds two-thirds of the 3 dpi cells and is the shallowest library:
  "rises at 3 dpi" and "shallow library" are hard to tell apart.
  Hr12, Hr24 and Hr27 replicate it at higher depth.
- Doublet and ambient checks on a cardiomyocyte-only matrix are blind to cardiomyocyte–fibroblast or cardiomyocyte–macrophage doublets.
- Immediate-early and heat-shock genes (fos, jun, egr1, atf3, hsp) are induced by tissue dissociation;
  use the dissociation control axis before calling a cell injury-responsive.
- Marker sets that share genes, plus a shared dependence on depth, can produce co-occurrence on their own.
- "Loss of maturity" genes (*ttn.2*, sarcomere, *cox7a1*, *ckma*) lose detection with depth;
  keep them secondary to "gain" genes for dedifferentiation, gated by the thinning test.

## Open questions for Eme

1. Is there a per-cell expert annotation besides the published labels? If so, keep it until the end; tell us its form only.
2. Which border-zone reference do you trust most: Wu 2016, Honkoop 2019, or the 2025 spatial data (amputation, not cryoinjury)?
3. For ventricle-only and atrium-only hearts, score per chamber or across all cardiomyocytes?
