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
22 hearts, 32,044 cells, 20 heart-bootstrap resamples, seen once.

1. **The published dedifferentiated cluster is 25 % of cells; the latent class puts the state at 2.5 %.**
   Its estimated sensitivity is 0.43 [0.28, 0.51] and specificity 0.75 [0.68, 0.85].
2. **That cluster sits on deeper libraries**: median depth 1.42× its own heart's median (other cells 0.91×);
   its sensitivity rises from 0.18–0.36 (low/mid depth tertile) to 0.70 (top tertile).
3. **The published proliferating cluster is specific but small**: 279 cells, specificity 0.997, sensitivity 0.20 [0.03, 0.42].
4. **κ vs the published annotation is 0 at every timepoint**, because protocol v1 makes no high-confidence
   positive call on either axis. Agreement is untestable on high-confidence calls; the protocol abstains.
5. **The latent class is not identifiable from independent labellers yet**: the UCell and score_genes arms
   share genes (both ≈ 0.97–0.99 sensitive), and split-half sets fire on ≤ 5 cells. Treat 1–3 as relative to
   "the two protocol arms agree", not to truth.

### Labeling protocol v1, run 3 (2026-11-07, hk/scores)

Files: `results/hk/protocol_v1_summary.csv`, `results/hk/protocol_v1_gate.csv`,
`results/runs/hk-protocol-v1/manifest.yaml`. Protocol: `docs/labeling-protocol.md`.

1. **No state axis beats its detection-matched decoys in a majority of hearts at any
   timepoint.** Injury-responsive, border zone, dedifferentiated and proliferating are all
   `not supported` in all cells and in the ventricle stratum (summary `verdict`,
   `supported_times`); the per-cell protocol calls 0–8 positives per axis out of 32,044 cells.
   Border zone is `depth_confounded` in the atrium stratum.
2. **The protocol can call a real program.** The dissociation / heat-shock control is
   `supported` at every timepoint in all cells and survives leave-one-heart-out at three of
   four (summary row `dissociation_control`). So "not supported" for the four states is not a
   protocol that calls nothing. Its marker-bootstrap stability is low (`boot_positive_stable_80`,
   5 resamples only).
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
   unchanged: not supported).

What 18 hearts cannot decide: a state carried by one or two hearts (7 dpi proliferation in
Hr8/Hr9) cannot pass a majority-of-hearts gate; whether ttn.2-high is "dedifferentiation" is
untestable because the gain genes are barely detected.

Not run in v1: partial conjunction (Simes), REML meta-analysis, dip test, taxonomy odds ratios,
multi-gene chamber score, 200-resample bootstrap.
