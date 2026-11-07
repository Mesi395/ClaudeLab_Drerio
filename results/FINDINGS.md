# Findings

Two tiers. **Claims** survived every check and have an artifact. **Observations** are
unverified. Numbers live in the linked files; this page only points at them.

## Claims

None yet. No result has been re-run from a clean checkout and passed thinning, decoys,
leave-one-heart-out and the confound check together; everything below is an observation.

## Observations

### Latent-class comparison vs the published annotation (2026-11-07, hk/results)

Commit 87df55f. Command (ROG, 30 s): `python scripts/compare_to_annotation.py data/derived/labels_v1.tsv.gz
data/derived/splithalf_v1.tsv.gz <owner's original h5ad, outside the repo> --column <published column>
--dediff <substr> --prolif <substr>`. Inputs: Card B's per-cell calls (labels_v1, commit 0f67e67) and
`scripts/split_half_calls.py` output; both stay on the ROG. Output: `results/hk/compare_to_annotation.txt`.
22 libraries (18 hearts), 32,044 cells, 20 resamples over the 18 hearts, seen once. The first run resampled
libraries, not hearts (Hr2a/b, H8a/v, Hr6a/v, Hr7a/v are one heart each); it was re-run after the fix and only the
intervals moved (L1, L3, README table).

L1. **The published dedifferentiated cluster is 25 % of cells; the latent-class prevalence is 2.5 %.** This is
   relative to "the two protocol arms agree" (L5), not an estimate of the true state size. Arms a and b use a
   nominal per-cell p < 0.05 (448 and 579 positives); the protocol's own call, which also needs BH and both arms,
   makes 0. The cluster's estimated sensitivity is 0.43 [0.26, 0.47] and specificity 0.75 [0.69, 0.88].
L2. **That cluster sits on deeper libraries**: median depth 1.42× its own heart's median (other cells 0.91×);
   its sensitivity rises from 0.18–0.36 (low/mid depth tertile) to 0.70 (top tertile). Deeper cells can also be
   real (larger, transcriptionally active cells), so depth alone does not show the cluster is an artefact.
L3. **The published proliferating cluster is specific but small**: 279 cells, specificity 0.997, sensitivity 0.20 [0.05, 0.35].
L4. **κ vs the published annotation is 0 at every timepoint**, because protocol v1 makes no high-confidence
   positive call on either axis. Agreement is untestable on high-confidence calls; the protocol abstains.
L5. **The latent class is not identifiable from independent labellers yet**: the UCell and score_genes arms
   share genes (both ≈ 0.97–0.99 sensitive), and split-half sets fire on ≤ 5 cells. Treat 1–3 as relative to
   "the two protocol arms agree", not to truth.

### Labeling protocol v1, run 3 (2026-11-07, hk/scores)

Files: `results/hk/protocol_v1_summary.csv`, `results/hk/protocol_v1_gate.csv`,
`results/runs/hk-protocol-v1/manifest.yaml`. Protocol: `docs/labeling-protocol.md`.

1. **No state axis beats its detection-matched decoys in a majority of hearts at any
   timepoint.** Injury-responsive, border zone, dedifferentiated and proliferating are all
   undetectable per cell at this depth, in all cells and in the ventricle stratum (summary `verdict`
   column, which literally reads `not supported`, and `supported_times`); the per-cell protocol calls 0–8
   positives per axis out of 32,044 cells. This is a statement about power with 3–10 genes per cell at
   600 UMIs, not about the biology.
   Border zone is `depth_confounded` in the atrium stratum.
2. **The protocol can call a highly expressed program.** The dissociation / heat-shock control is
   `supported` at every timepoint in all cells and survives leave-one-heart-out at three of
   four (summary row `dissociation_control`). That set is highly expressed (in the 3 dpi ventricle pseudobulk,
   hsp70l has baseMean 1,086 against a median of 3.4 over 7,903 genes), so this shows the code runs; it does not
   show that a small state with few genes would be callable. No power-matched spike-in has been run.
   Its marker-bootstrap stability is low (`boot_positive_stable_80`, 5 resamples only).
3. **Thinning changes few calls**: Cohen's κ thinned vs unthinned is in the summary column
   `kappa_thinned_vs_full` for every axis × stratum.
4. **Morphine-study hearts (Hr8, Hr11) do not change any verdict** (`supported_without_morphine_study`).
5. **Dissociation is higher in atrium than in ventricle in every heart with both chambers**, and
   varies widely between hearts (`results/hk/chamber_contrast.csv`,
   `results/hk/dissociation_per_heart.csv`). Hr11, the shallowest library, has the highest
   median dissociation score but almost no dissociation-positive cells, because most of its
   cells are below 600 UMIs and get the rank call only.
6. **The injury-responsive score is higher in ventricle than atrium in 16 of 18 hearts,
   including all uninjured hearts** (`chamber_contrast.csv`), so the contrast is a chamber
   baseline, not a local injury response. Border zone is lower in ventricle in every heart.
7. **Injury does not depend on dissociation-high cells**: within-heart partial ρ of injury with
   dissociation given log UMI, and the injury verdict with and without the top 10%
   dissociation cells per heart, are in `results/hk/injury_vs_dissociation.csv` (verdict
   unchanged: undetectable per cell).

### Pseudobulk DE per heart x chamber (2026-11-07, hk/scores)

Files: `results/hk/pseudobulk_de_summary.csv`, `pseudobulk_thinning_robustness.csv`,
`pseudobulk_samples.csv`, `pseudobulk_dropped.csv`, `pseudobulk_loho_dpi3.csv`,
`pseudobulk_vs_literature.csv`, per-model tables `pseudobulk_de_{raw,thinned600}_*.csv.gz`. Script:
`scripts/pseudobulk_de.py` (pydeseq2 0.5.4, hearts as replicates).

8. **Ventricle 3 dpi vs uninjured (4 vs 4 hearts) gives a DEG list; 7 dpi vs uninjured (7 vs 4)
   gives almost none** (`pseudobulk_de_summary.csv`, after removing dissociation/heat-shock
   and ambient non-CM genes, padj < 0.05, |log2FC| > 1).
9. **A local (ventricle-specific) injury effect is not detectable with the myosin chamber call**: the
   interaction model `~heart + chamber + vent_injured` finds 0–1 DEGs at 3 and 7 dpi, raw or thinned. The model
   is underpowered (4 v 4 hearts, chamber call 57–87 % accurate, no chamber-sorted libraries at 3 dpi), so this is
   not a null. The ventricle 3 dpi DEGs may be heart-level (both chambers) effects, which the interaction would
   cancel, or ventricle effects it cannot resolve.
10. **Thinning to 600 UMIs keeps the log2FC ranking (Spearman in
    `pseudobulk_thinning_robustness.csv`) but the DEG count falls from 67 to 26** (DEG Jaccard there).
    Thinning cuts counts and so power for any effect; the ranking is kept.
11. **The data-derived 3 dpi up set shares no gene with the 4-gene literature injury set**
    (`pseudobulk_vs_literature.csv`); that set is too small to overlap anything. The up set holds genes
    named in border-zone and dedifferentiation work (for example pdk2a and myom1b) and also
    mitochondrial-genome transcripts (NC-002333.*, mt-co1) and unannotated loci, which is how a depth/quality
    difference looks. It has not yet been compared with the full border-zone and dedifferentiation marker sets
    or tested against a permutation null, so no injury claim is made either way.
12. Leave-one-heart-out top-30 up-set Jaccard against the full set is in `pseudobulk_loho_dpi3.csv`;
    dropping Hr11 gives the lowest overlap.

What 18 hearts cannot decide: a state carried by one or two hearts (7 dpi proliferation in
Hr8/Hr9) cannot pass a majority-of-hearts gate; whether ttn.2-high is "dedifferentiation" is
untestable because the gain genes are barely detected.

Not run in v1: partial conjunction (Simes), REML meta-analysis, dip test, taxonomy odds ratios,
multi-gene chamber score, 200-resample bootstrap.
