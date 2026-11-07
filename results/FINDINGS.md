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

### Pseudobulk DE per heart x chamber (2026-11-07, hk/scores)

Files: `results/hk/pseudobulk_de_summary.csv`, `pseudobulk_thinning_robustness.csv`,
`pseudobulk_samples.csv`, `pseudobulk_dropped.csv`, `pseudobulk_loho_dpi3.csv`,
`pseudobulk_vs_literature.csv`, per-model tables `pseudobulk_de_{raw,thinned600}_*.csv.gz`. Script:
`scripts/pseudobulk_de.py` (pydeseq2 0.5.4, hearts as replicates).

8. **Ventricle 3 dpi vs uninjured (4 vs 4 hearts) gives a DEG list; 7 dpi vs uninjured (7 vs 4)
   gives almost none** (`pseudobulk_de_summary.csv`, after removing dissociation/heat-shock
   and ambient non-CM genes, padj < 0.05, |log2FC| > 1).
9. **The local (ventricle-specific) injury effect is near zero**: the interaction model
   `~heart + chamber + vent_injured` finds 0–1 DEGs at 3 and 7 dpi, raw or thinned. So the
   ventricle 3 dpi DEGs are mostly heart-level (both chambers) effects: systemic, dissociation,
   depth or batch, which the interaction cancels.
10. **Thinning to 600 UMIs keeps the log2FC ranking (Spearman in
    `pseudobulk_thinning_robustness.csv`) but changes the DEG set** (DEG Jaccard there): about
    two-thirds of the raw 3 dpi DEGs need depth to reach significance.
11. **The data-derived 3 dpi up set shares no gene with the literature injury set**
    (`pseudobulk_vs_literature.csv`), and it includes mitochondrial-genome transcripts
    (NC-002333.*, mt-co1) and unannotated loci, which is how a depth/quality difference looks.
    Two independent routes do not agree, so there is no injury claim.
12. Leave-one-heart-out top-30 up-set Jaccard against the full set is in `pseudobulk_loho_dpi3.csv`;
    dropping Hr11 gives the lowest overlap.
13. **Set test, programme at sample level** (`pseudobulk_set_test.csv`; mean DESeq2 Wald stat
    vs 500 baseMean-matched decoy sets, mitochondrial genes dropped). Ventricle 3 dpi vs Ctrl,
    raw: dissociation control z 9.1 (p 0.002), memory border-zone core z 5.1 (p 0.002),
    dedifferentiated z 4.0 (p 0.002), glycolysis hk1/pkma/pdk2a z 3.4 (p 0.006), markers_v1
    border_zone z 2.2 (p 0.018); the core without AP-1 (nppb, mustn1b, myh7) z 1.7 (p 0.052);
    injury_responsive and proliferating do not beat their decoys. Interaction model, raw: memory
    core p 0.002, glycolysis 0.004, border-zone core and injury_responsive 0.016, dissociation
    control z 9.8; thinned, only memory core, glycolysis and dissociation control stay below
    0.05. The memory core's lead comes from its AP-1 genes, which are dissociation genes.
    The 3 dpi up set overlaps dedifferentiated (myom1b, pdk2a; hypergeometric p 0.001 on 18,794
    genes) and glycolysis (pdk2a), and no border-zone set (raw; thinned adds myh7).
14. **Permutation calibration** (`pseudobulk_permutation.csv`, all 70 label assignments of the
    4 vs 4 ventricles): 67 observed DEGs; only the observed labels and their exact mirror (69)
    reach that, the next is 27, so p = 2/70 = 0.029, the floor of this design. Memory core,
    dedifferentiated and glycolysis set stats are the highest of the 70 (p = 1/70, the floor);
    markers_v1 border_zone 0.043, core without AP-1 0.071, dissociation control 0.13.
15. **30 dpi vs Ctrl negative control** (3 vs 4 hearts): 5 DEGs, so no batch-sized DEG list.
    But the border-zone and AP-1 set stats at 30 dpi equal or exceed 3 dpi (memory core 2.80
    vs 2.65, core without AP-1 1.30 vs 1.38, markers_v1 border-zone core 1.62 vs 1.04,
    dissociation control 2.16 vs 2.29); only glycolysis (0.74 vs 2.68) and dedifferentiated
    (0.08 vs 1.28) fall back. The border-zone shift therefore does not follow injury time and
    looks like a control-group difference; the glycolysis/dedifferentiation programme at sample
    level is the one that behaves like a 3 dpi effect. No decoy test at 30 dpi; Hr19 has 55 cells.

What 18 hearts cannot decide: a state carried by one or two hearts (7 dpi proliferation in
Hr8/Hr9) cannot pass a majority-of-hearts gate; whether ttn.2-high is "dedifferentiation" is
untestable because the gain genes are barely detected.

Not run in v1: partial conjunction (Simes), REML meta-analysis, dip test, taxonomy odds ratios,
multi-gene chamber score, 200-resample bootstrap.
