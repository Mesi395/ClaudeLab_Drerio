# Findings

Two tiers. **Claims** survived every check and have an artifact. **Observations** are
unverified. Numbers live in the linked files; this page only points at them.

## Claims

None yet.

## Observations

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
