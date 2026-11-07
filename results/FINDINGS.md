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

16. **Leave-one-heart-out, size-free** (`pseudobulk_loho_rank.csv`; 8 ventricle hearts, each left out
    in turn; null = 10 label shuffles, each against its own full fit). The Jaccard in item 12 was
    mostly a list-size artefact: median over left-out hearts, 93% of the full top-30 up genes stay
    in the left-out top 100, top-30 overlap coefficient 0.75, every full up DEG keeps its sign,
    Wald-stat Spearman 0.97 (0.89 over the DEGs). Shuffles give Spearman 0.94 (0.83 over DEGs),
    top-30 rank overlap 0.63 and top-100 containment 0.95. So the up genes are stable, but only a
    little above noise (Spearman, overlap) and not above it for top-100 containment: with 8 hearts,
    dropping one leaves 7/8 of the data, so this test cannot separate signal from noise; item 14
    (permutation) carries the evidence. 7 of 10 shuffled full fits have no DEG, so DEG-based
    shuffle medians rest on 3 shuffles (sign agreement and DEG-list overlap are not comparable).

17. **Per-heart set scores, exact heart relabelling** (`pseudobulk_set_score*.csv`; script
    pre-registered in f8ebb07 before its first run). One score per ventricle heart and set; hearts,
    not genes, are relabelled, so set genes that move together cannot inflate the result. One-sided
    exact p (70 or 35 relabellings, floors 0.014 / 0.029), BH over the six primary sets. Raw counts:
    - dedifferentiated: up at 3 dpi (p 0.014, q 0.043, Hedges' g 4.0), back to control level at
      30 dpi (p 0.46), 3 > 30 dpi p 0.029 (the floor; q 0.17). The only set that formally fades;
      nominal after BH.
    - glycolysis (hk1, pkma, pdk2a): up at 3 dpi (p 0.014, q 0.043, g 2.2) and at 7 dpi (p 0.006);
      30 dpi not above control (p 0.14), but 3 > 30 dpi only p 0.086, so its fade is not shown.
    - memory border-zone core: up at 3 dpi p 0.029 (q 0.057); as high at 30 dpi (means 0.41 vs
      0.45, 3 > 30 p 0.49, 30 > control p 0.057); without its four AP-1 genes p 0.071 at 3 dpi.
      The dissociation control has the same shape (3 dpi 0.40, 30 dpi 0.41, control -0.42). Item 15
      is now tested, not judged by eye: the border-zone shift separates uninjured from all injured
      hearts and does not follow injury time.
    - full border_zone, injury_responsive, proliferating: not shown at 3 dpi (p 0.057, 0.30, 0.086).
    - The dissociation control, z 9.1 in the competitive test (item 13), gives p 0.14 here: the
      competitive test overstated sets whose genes move together. Glycolysis and dedifferentiation
      hold under both tests.
    - Thinned to 600 UMIs: the same calls for dedifferentiated (fades) and glycolysis (3 dpi only).
      The 30 dpi arm is 3 hearts, one with 55 cells; the 3 dpi heart Hr27 has 47 cells.

18. **Statistical plan v2, pre-registered** (`docs/statistical-plan.md` and
    `scripts/pseudobulk_set_plan.py` committed in cadd8d6 before the run; `pseudobulk_set_plan*.csv`).
    Heart = unit; 7 dpi vs control is the primary family because Hr24 and Hr27 (3 dpi) are 10x v3
    per GEO and Hr11 comes from a second series. Welch t, exact relabelling, BH over six sets.
    Supersedes items 12, 13, 15 and 16 as evidence (plan section 4). Raw counts:
    - **Glycolysis (hk1, pkma, pdk2a) is the one injury programme shown.** 7 dpi +0.92 z [0.19, 1.66],
      p 0.006, q 0.036 (primary); 3 dpi +1.58, p 0.014, q 0.043 (chemistry/series confounded). It
      does not follow handling or ambient fibroblast RNA across hearts (within-condition r 0.03 and
      0.06), no single heart carries it, and it holds at 600 UMIs and without the two smallest
      hearts. It is ventricle-local at 7 dpi (ventricle minus atrium +0.67, p 0.006). Return by
      30 dpi not shown (3 vs 30 dpi p 0.086). Consistent with Honkoop et al. 2019, who confirmed
      these three genes by in situ hybridisation in border-zone cardiomyocytes at 7 dpi.
    - **Dedifferentiation: 3 dpi only, and not beyond glycolysis in a convincing way.** 3 dpi +0.58
      [0.36, 0.80], p 0.014, q 0.043 (confounded); 7 dpi p 0.082; 3 > 30 dpi p 0.029 (the floor).
      Without its nine glycolytic genes it still rises at 3 dpi (p 0.014), which passes the
      pre-registered rule, but without the four library-size genes as well (ttn.1, ttn.2, bves,
      synpo2lb) p 0.14: the non-glycolytic signal rests on genes that track library size.
    - **Border zone, memory core, injury-responsive, proliferating: not shown at heart level.** The
      memory core (r 0.91) and injury-responsive (0.77) follow handling across hearts and cannot be
      separated from it with 18 hearts. Exploratory, not a call: the ventricle-minus-atrium line,
      which cancels heart-level handling, rises at 7 dpi for the memory core (+0.58, p 0.009) and the
      border-zone set (+0.52, p 0.030; this set has no AP-1 genes).
    - **Ambient fibroblast RNA in the cardiomyocyte pseudobulk follows injury time** (col1a1a,
      col1a2, postnb, dcn: 7 dpi +1.17, p 0.006; 3 dpi +1.84, p 0.014; 3 > 30 dpi p 0.029). Any
      per-cell cardiomyocyte state call has to account for it.
    - Genome-wide shift diagnostic (200 random 10-gene sets, share with p ≤ 0.05): 7 dpi 0.05,
      3 dpi 0.095, 30 dpi 0.025. Hearts at 7 dpi do not differ overall; at 3 dpi they differ a
      little, as the chemistry mix predicts.
    - **Limits of item 18** (from the independent statistics review, research s3, 2026-11-07):
      the glycolysis rise at 7 dpi had already been seen before the plan was written (item 17,
      p 0.006, and an exploratory 33-test scan, q 0.26). The plan fixes the analysis, but this
      result is not a blind confirmation; that needs new hearts or an independent dataset. At gene
      level the 7 dpi rise is not significant (DESeq2 padj pkma 0.46, pdk2a 0.68); it shows only as
      a set. A different per-heart score (mean z of logCPM, simpler veto) gave the 3 dpi glycolysis
      rise p 0.071 instead of 0.014, so the 3 dpi p depends on the scoring. The ventricle-minus-atrium
      hint for the memory core is partly circular: whole-heart cells are assigned to a chamber by
      myh7/myh7l vs myh6 (`src/zf/io.py`), and myh7 is a memory-core gene. Only the border-zone
      set's hint (no myh7) remains, and it is exploratory.

19. **Every cell, or a subgroup? (plan section 2F, pre-registered in af1018a; NaN handling fixed in
    5ab8dba, same numbers)** (`cell_tail*.csv`). Per ventricle heart, cells thinned to exactly 600
    UMIs; tail excess (share above the heart's decoy 95th percentile, minus the decoys' share) and
    dispersion excess (variance-to-mean beyond detection-matched decoys). Both ignore a uniform shift
    (synthetic check in the script). Results:
    - **Glycolysis: every cell shifts** at 7 and 3 dpi: tail excess -0.014 and -0.016, dispersion
      excess -0.07 and -0.13 (p 0.84, 0.86), while the heart-level level rose (item 18). A small
      subgroup cannot be excluded: with 3 genes at 600 UMIs this test sees a subgroup only if it is
      large or strong (the synthetic 10%-of-cells x5 case gives tail +0.07).
    - **Dedifferentiated at 3 dpi: every cell shifts** (dispersion excess +0.01, p 0.49).
    - **Memory border-zone core and injury-responsive: hints of a subgroup, not shown.** Dispersion
      excess at 3 dpi +0.64 (p 0.029) and 7 dpi +0.49 (p 0.048) for the core; injury-responsive 7 dpi
      +0.17 (p 0.030); core tail excess at 3 dpi p 0.043. None reaches BH q 0.05 (min 0.12). A border
      zone is by definition a subgroup near the wound, so this is the expected direction; the core
      also shares genes with the handling set.
    - **Proliferating: not computable** (its genes are undetected in some hearts at 600 UMIs).
