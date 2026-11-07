# Labeling protocol v1

A fixed decision procedure. For every cardiomyocyte and every state axis it gives a call
(`positive` / `negative` / `uncertain`), two scores, two p-values and the artefact flags that
vetoed it. Code: `src/zf/protocol.py`, entry point `label_axis(adata, genes)`. Runner:
`scripts/run_protocol.py`. Pre-registered in `JOURNAL.md` (2026-11-07 16:00) before any call
was looked at; amendments there are dated and say why.

The four states are hypotheses. An axis that does not beat its matched decoys in a majority of
hearts is reported as **not supported**, not as a weak score.

## Rules (per axis)

1. **Thin.** Binomially thin every cell to T = 600 UMIs. Cells below 600 keep their counts, are
   flagged `low_depth` and get only the rank call (rule 4).
2. **Score twice.** UCell on raw counts (rank-based, max rank 1,500) and the mean
   log-normalised expression on thinned counts.
3. **Null.** 500 decoy gene sets per library; each marker is replaced by a random gene from the
   same detection-rate bin (20 quantile bins) of the same library. p per method: the cell's z
   against its own decoys, ranked against the pooled decoy z's of all cells in its library
   (scDRS-style).
4. **Call.** Benjamini–Hochberg within each heart, per method. `positive`: q ≤ 0.05 in both
   methods (UCell only for `low_depth` cells). `negative`: p ≥ 0.10 in both. Else `uncertain`.
5. **Veto → uncertain.** `doublet_flag` (zf.confounds); ≥ 2 of the fibroblast control genes
   (col1a1a, col1a2, postnb, dcn) detected; hemoglobin fraction above the library's 95th
   percentile; pct_mito > 30.
6. **Gate.** Per heart, the positive fraction must beat the decoy false-positive rate (20 decoy
   sets re-called as pseudo-axes through rules 3–5; one-sided binomial, p < 0.05). A timepoint
   is **supported** when ≥ ⌊n/2⌋+1 of its n hearts pass (3/4 uninjured, 3/4 at 3 dpi, 4/7 at
   7 dpi, 2/3 at 30 dpi).
7. **High confidence.** A positive whose timepoint is supported and stays supported with any
   one heart left out.
8. **Eligibility.** If the within-heart partial Spearman |ρ| of the thinned score with log UMI,
   or with hemoglobin fraction, is ≥ 0.10, the axis is `depth_confounded` and emits no
   positives.

Axis verdict: `supported` (≥ 1 timepoint supported) · `not supported` · `depth_confounded`.

## Grading without ground truth

- Thinning agreement: Cohen's κ between calls with and without thinning.
- Marker bootstrap: the share of positives that stay positive in ≥ 80% of resamples of 80% of
  the markers.
- Leave-one-heart-out: rule 7.
- Decoy false-positive rate: rule 6.
- Independence from confounds: rule 8.

Sensitivity: the gate without Hr8 and Hr11 (the morphine study, GSE158919).

## Fixed, never tuned

T, decoy count, bins, q/p cuts, the u rule, vetoes, ρ cut, metric pass marks, disjoint marker
sets. The protocol never optimises the separation of uninjured from injured hearts, because
that rewards depth.

## Known limits

- With 18 hearts, a state carried by one or two hearts (7 dpi proliferation in Hr8/Hr9) cannot
  pass the gate; it is reported per heart, not as a replicated population claim.
- Whether ttn.2-high cells are "dedifferentiated" is untestable here: the gain genes (gata4,
  runx1, tbx20) are detected in ≤ 1% of cells.
- Chamber strata use Card A's two-gene myosin call (57–87% accurate on libraries with a known
  chamber) for whole-heart libraries.

## Dissociation control (domain input 2026-11-07 15:09)

Cells were dissociated enzymatically at 37 °C, which induces immediate-early and heat-shock
genes (van den Brink 2017, doi:10.1038/nmeth.4437; Neuschulz 2022, doi:10.15252/msb.202211147).
`dissociation_control` is scored through rules 1–8 like an axis, but it is never a state. Its
per-heart level is in `results/hk/dissociation_per_heart.csv`. The injury axis is not regressed
on it (immediate-early genes are part of the real border-zone program); instead the runner
reports the within-heart partial ρ of injury with dissociation and the injury verdict with and
without the top 10% dissociation cells of each heart (`results/hk/injury_vs_dissociation.csv`).
The ventricle − atrium contrast per heart (`results/hk/chamber_contrast.csv`) separates the two:
dissociation hits both chambers alike, a local injury response is ventricle-only.
