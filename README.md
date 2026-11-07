# ClaudeLab_Drerio

Can an agent annotate cell states when there is no ground truth?

A research loop built at the [Claude Impact Lab Berlin 2026](https://github.com/xemantic/claude-impact-lab-2026)
on a problem from the [Junker lab](https://junkerlab.com/) (MDC Berlin).

## The problem

After cryoinjury, zebrafish cardiomyocytes take on several states:
injury-responsive, border-zone-like, dedifferentiated and proliferating.
Two questions are open.
Do these states co-occur in the same cells?
And what do published markers really mark?
For example, one supposed dedifferentiation marker mostly tracks library size, not biology.

The full problem statement is in [docs/problem.md](docs/problem.md).

## The data

The data is the public zebrafish heart single-cell data from Hu et al. (2022).
It has 32,044 cells × 18,794 genes from 22 libraries.
The libraries cover uninjured hearts and 3, 7 and 30 days post injury.
A held-out expert annotation is used only for a final agreement check.

## What the agent does

1. Finds marker genes for each state, from literature or spatial data.
2. Scores every cell on each state axis.
3. Tests the scores against technical confounds: library size, doublets and ambient RNA.
4. Tests whether the axes coincide or come apart.

## Answers to the two questions

### Short answer (statistical plan v2, 2026-11-07)

The heart is the unit of evidence: 18 hearts, not 32,044 cells. The plan, its decision rules and its
code were committed before the run ([docs/statistical-plan.md](docs/statistical-plan.md); results in
[results/FINDINGS.md](results/FINDINGS.md), items 18 and 19).

- **Shown:** glycolysis genes (hk1, pkma, pdk2a) rise in injured ventricles at 7 dpi (the primary
  test, q 0.036) and at 3 dpi. The rise is in every cell, not a subgroup. It is not explained by
  dissociation stress, fibroblast RNA, a single heart or sequencing depth. This is consistent with
  Honkoop et al. 2019. It is not a blind confirmation: the 7 dpi rise was seen before the plan was
  written, so confirming it needs new hearts.
- **Not shown:** border zone, injury-responsive and proliferating at heart level, and
  dedifferentiation beyond glycolysis (its non-glycolytic part rests on genes that track library
  size). Border-zone scores follow dissociation stress from heart to heart; hints of a border-zone
  subgroup do not survive multiple-testing correction.
- **Per cell:** no state can be called at about 600 UMIs per cell. That is a detection limit, not
  absence.
- **The design limits what can be decided:** two of the four 3 dpi hearts were sequenced with 10x v3
  and one comes from a second series (GEO records), so every 3 dpi number is confounded and 7 dpi
  carries the claims. 30 dpi has 3 hearts, so "fades" or "persists" can be estimated, not shown.

**Which annotation to trust:** one that (1) passes calibration controls, (2) replicates across
hearts, (3) follows injury time and can be separated from handling, ambient RNA and chemistry, and
(4) agrees with evidence that shares no genes with it. By these criteria the glycolytic shift is the
only label that earns trust in this dataset. The published dedifferentiated cluster (25 % of cells,
enriched in deep cells) does not.

**Checks added after the event:** a design table (chemistry, series) before any test; one
pre-registered heart-level test per gene set with exact relabelling of hearts and stated floors;
gene splits for overlapping sets; handling and ambient-RNA criteria; a mean-shift vs subgroup test.
An expert should still add new hearts or an independent dataset for a blind confirmation, in situ
staining to see whether glycolysis is wound-local, and a limma `fry` cross-check once R is available.

### How the statistics were chosen

Each choice below came from the design table first (which hearts, which chemistry, which series),
then from what the data can carry. The options are the ones we ran during the event (FINDINGS items
1-17) and the standard alternatives in the references of [docs/statistical-plan.md](docs/statistical-plan.md).
A senior design review (Fable) changed the first draft before anything ran.

| choice | options | chosen | why |
|---|---|---|---|
| unit of evidence | each cell (32,044) · each heart (18) | heart (pseudobulk) | Cells from one heart share its dissection, dissociation and sequencing run. Counting them as independent turns any heart-to-heart difference into a tiny p-value (Squair et al. 2021, Crowell et al. 2020). |
| what is tested | single genes · gene sets | gene sets for states; genes only as support (plan B) | The question is about states, i.e. programmes of many genes. With 4-7 hearts per group a single gene has little power, while a set pools many weak genes. |
| kind of set test | competitive: set genes vs other genes (item 13) · self-contained: injured hearts vs control hearts | self-contained | The competitive test treats genes as independent. Genes in a set move together, so the handling set reached z 9.1 there but only p 0.14 when hearts were relabelled (Goeman & Bühlmann 2007; Wu & Smyth 2012). |
| p-value | normal-theory t-test · every possible relabelling of hearts | Welch t, p from exact relabelling, one-sided | At 4-7 hearts a normal table is not trustworthy. Relabelling assumes no distribution, and its smallest possible p is known in advance (0.003 at 7 dpi, 0.014 at 3 dpi, 0.029 at 30 dpi). Welch, because control hearts vary less than injured ones. One-sided, because "injury raises the programme" was stated before the run. |
| primary contrast | 3 dpi · 7 dpi | 7 dpi vs control | Two of the four 3 dpi hearts are 10x v3 and one comes from a second series (GEO). Chemistry sits only inside 3 dpi, so it cannot be adjusted for. 7 dpi has 7 hearts and one chemistry. |
| nuisances (dissociation stress, fibroblast RNA) | regress them out · report them as criteria and split the gene sets | criteria (r ≥ 0.6 within condition) and gene splits | The handling genes overlap the tested sets, and handling also rises with injury itself, so regressing it out would remove real signal. The review removed the handling covariate and the chemistry term. |
| multiplicity | one correction over everything · an omnibus gate · a fixed sequence | Benjamini-Hochberg in a fixed sequence: 7 dpi family, then 3 dpi | This controls the share of false calls, and the confounded 3 dpi family cannot dilute the primary one. The review removed the omnibus gate. |
| every cell or a subgroup | cluster cells and count them · call each cell · shape test against decoy sets within each heart | shape test (plan 2F) | At about 600 UMIs, per-cell calls cannot see the states, and clustering puts every cell in some state anyway. Tail and dispersion excess ask directly whether some cells switch on, and they ignore a uniform shift. |
| order of work | test, then interpret (v1) · write the plan and code first | plan and script committed before the run, run once | Choosing among many tests after seeing results inflates false positives. The exception is stated: the 7 dpi glycolysis rise was seen before the plan. |

Retired from v1: leave-one-heart-out gene overlap (at 8 hearts it cannot separate signal from
noise) and latent-class agreement as a measure of trust (its labellers share marker genes, so it
breaks the independence it assumes). Planned, not run: DESeq2 per gene (B), detection power and
Milo neighbourhoods (C), and co-occurrence (D).

### How we decide which annotation to trust, with no ground truth

Two layers (Card C, during the event).

1. **Protocol stability** (Card B, [docs/labeling-protocol.md](docs/labeling-protocol.md)).
   A call counts only if it survives count thinning, detection-matched decoy gene sets and a per-heart gate.
2. **Latent-class model** (Dawid–Skene, `scripts/compare_to_annotation.py`).
   Five imperfect labellers per axis: (a) protocol UCell arm, (b) protocol score_genes arm on thinned counts,
   (c, d) two disjoint split-half marker sets, (e) the published annotation mapped to the axis it claims.
   EM estimates each labeller's sensitivity and specificity; 20 bootstrap resamples over hearts (not cells).
   Run once on the ROG, 30 s, 32,044 cells, 22 libraries (18 hearts). Output: [results/hk/compare_to_annotation.txt](results/hk/compare_to_annotation.txt).

What it says (observations, one run):

| Axis | Labeller | Sensitivity [90 % heart-bootstrap] | Specificity |
|---|---|---|---|
| dedifferentiated | a UCell | 0.99 [0.99, 1.00] | 0.998 |
| | b score_genes | 0.97 [0.92, 0.98] | 1.000 |
| | c split-half 1 | 0.00 | 1.000 |
| | e published annotation | 0.43 [0.26, 0.47] | 0.75 [0.69, 0.88] |
| proliferating | a UCell | 0.99 [0.98, 1.00] | 0.999 |
| | b score_genes | 0.78 [0.73, 0.86] | 0.998 |
| | c / d split-half | 0.03 / 0.01 | 1.000 |
| | e published annotation | 0.20 [0.05, 0.35] | 0.997 |

How to read it, honestly:
- **The latent class is mostly "a and b agree".** a and b share their genes, so conditional independence fails,
  and their near-perfect numbers are partly self-agreement. The split-half sets were meant to break this,
  but with half the genes they almost never fire (c on dedifferentiated: 0 positive cells; d had < 2 genes and was skipped).
- **The published dedifferentiated cluster disagrees with every marker-based labeller:** it holds 25 % of cells
  (8,137) against a latent prevalence of 2.5 %. Its cells sit at 1.42× their own heart's median depth
  (other cells 0.91×), and its estimated sensitivity rises from 0.18–0.36 in the low/mid depth tertiles to 0.70 in the top tertile.
  So the cluster is enriched in deep cells; depth and biology are not separable here. The published proliferating cluster is
  small (279 cells), specific (0.997) and insensitive (0.20).
- **Cohen's κ vs the published annotation is 0 at every timepoint**, because protocol v1 makes no
  high-confidence positive calls on either axis (dedifferentiated: 0 positive, 3,683 uncertain;
  proliferating: 8 positive, 0 high-confidence). Our protocol abstains rather than agrees or disagrees.
  Of the cells it calls negative, the published dedifferentiated cluster holds 68 % at 3 dpi, 8 % at 7 dpi, 32 % at 30 dpi and 13 % in controls.
- **Depth tertiles:** within tertiles UCell sensitivity drops to 0.26 for dedifferentiated in the deepest third,
  so arm (a) is depth-sensitive too. No labeller is depth-free.

### Which checks the agent ran itself, and which an expert should add

Ran by the agent: count thinning, detection-matched decoys, leave-one-heart-out gate, library-size / doublet /
ambient-RNA confounds, dissociation-stress control axis, within-heart partial correlation, ventricle-minus-atrium
contrast (Card B), and the latent-class comparison with a heart bootstrap and depth tertiles (this section).

An expert should add:
- **Evidence outside the transcriptome.** Dedifferentiation and proliferation are also marked by physical force
  on the cell, morphology and sarcomere disassembly (Beisaw et al. 2020, doi:10.1161/CIRCRESAHA.119.316167;
  Jopling et al. 2010, doi:10.1038/nature08899), protein markers (PCNA or EdU incorporation, phospho-histone H3)
  and position relative to the wound. None of these is in scRNA-seq.
- **Dissociation artefact.** The cells were dissociated with enzymes at 37 °C, per the data owner. Heat-shock and dissociation stress
  vary between samples and can compound the injury stress response. Our checks separate the dissociation axis
  itself (control axis), its within-heart correlation with each state (partial correlation), and chamber-specific
  signal (ventricle-minus-atrium; only the ventricle is injured). They cannot separate the **injury × dissociation
  interaction**: injured tissue may respond to dissociation differently. The fix is experimental: metabolic labelling
  of new transcripts separates dissociation-induced transcription (Neuschulz et al. 2022, *Mol Syst Biol*).
- **A truly independent labeller.** Split-half marker sets were too small to act as one; an expert annotation
  made blind to these markers, or a spatial/imaging label, would make the latent-class model identifiable.
  Emi can run `scripts/compare_to_annotation.py` on her own annotation without sharing it; it prints numbers only.

## Decisions and reasons

Every analysis choice, why it was made, and its source. Generated from [docs/decisions.md](docs/decisions.md).

| Decision | Reason | Source |
|---|---|---|
| No morphine covariate | constant ("Ctrl") for all 32,044 cells, and one heart per timepoint would make it the same as regressing on heart | GEO GSE158919 sample records (Hr8, Hr11 "morphine treatment: no"); data check of the `morphine` column; Lelek et al. 2020, doi:10.1101/2020.10.01.322560 |
| Sensitivity run without Hr8 and Hr11 | they come from a separate experiment; Hr11 is the shallowest library with the most ambient RNA | GEO GSE158919; `results/hk/library_flags.tsv` |
| Never regress on library, heart or timepoint; chemistry not a covariate | chemistry exists only at 3 dpi, so it is inseparable from biology | `config/dataset.yaml` |
| Thin to 600 UMIs before the score_genes-type score | cell-cycle and ttn.2 signal track depth within library | senior review 15:39 (observation tier); `results/hk/protocol_v1_summary.csv` kappa columns |
| Decoys matched on detection rate within each library (500 sets, 20 bins) | makes the null follow each library's depth; scDRS matches on mean and variance by default | Zhang et al. 2022 scDRS, doi:10.1038/s41588-022-01167-z |
| Pooled-z p-values instead of (1+#)/(1+B) | the plain empirical p is floored at 1/501, so BH within a heart of thousands of cells called 0 positives on every axis in run 1 | scDRS (same DOI); JOURNAL.md 2026-11-07 16:07 |
| Non-CM veto uses the 4 fibroblast control genes of `config/axes.yaml`, not the 106 soup genes | the soup-gene version vetoed 59% of cells in run 1 | `config/axes.yaml` controls; JOURNAL.md 2026-11-07 16:07 |
| Use continuous `ambient_frac` / hemoglobin fraction, not the ambient flag | at ~700 UMIs one count reads as 10–30% contamination | `src/zf/confounds.py` docstring |
| Chamber: Card A's myosin call (myh7, myh7l vs myh6; vmhc is myh7 in GRCz11, confirmed by Emi 16:32) for whole-heart cells; uncertain chamber cells are left out of both strata | the multi-gene thinned chamber score was not built in time | `results/hk/va_accuracy.tsv` |
| Emi's biology calls (16:32, 16:39): vmhc is myh7; nppa stays out of every score because it marks mainly the atrium (`ambiguous.yaml`); nppb is used in the border-zone score | domain expert decision | `results/hk/markers_v1/ambiguous.yaml`, `border_zone.yaml` |
| Partial conjunction by Simes, REML meta-analysis, dip test and taxonomy odds ratios: not run in v1 | time | JOURNAL.md 2026-11-07 16:00 |
| Dissociation / heat-shock control axis scored with the same protocol (Card A's cited dissociation set plus hsp70/hspa*, hsp90aa1*, dnajb1*, fos*, jun*, egr1, atf3 present in the data, minus every scored axis gene) | cells were dissociated enzymatically at 37 °C (per the data owner), which induces immediate-early and heat-shock genes that vary between samples and overlap the injury response | Neuschulz et al. 2022 MSB, doi:10.15252/msb.202211147; van den Brink et al. 2017 Nat Methods, doi:10.1038/nmeth.4437; domain input 2026-11-07 15:09 |
| Injury axis judged against dissociation, never regressed on it | IEGs are also part of the real border-zone program (AP-1), so removing the dissociation score would remove biology; instead report the within-heart partial ρ (given log UMI) and the injury verdict with and without the top 10% dissociation cells per heart | `results/hk/injury_vs_dissociation.csv` |
| Ventricle-minus-atrium contrast per heart for injury, border zone and dissociation | dissociation hits both chambers alike, a local injury response is ventricle-only | owner 15:52 + domain input 15:09; `results/hk/chamber_contrast.csv` |
| Marker bootstrap with 5 resamples in this run, not 200 | time (16:25 cut); the metric is reported with its resample count in `results/hk/protocol_v1_summary.csv` | JOURNAL.md 2026-11-07 16:12 |
| Latent-class (Dawid–Skene) model to judge annotations with no ground truth | standard way to estimate each imperfect labeller's sensitivity/specificity without truth; `uncertain` treated as missing | Dawid & Skene 1979, J R Stat Soc C 28:20; owner 15:43 |
| Bootstrap over hearts, 20 resamples | cells within a heart are not independent; 20 is the time cut (owner 16:20) | `results/hk/compare_to_annotation.txt` |
| Also fit within depth tertiles | depth drives every labeller, which breaks conditional independence | owner 15:43 |
| Published annotation read in memory from the owner's file outside the repo; category strings passed on the command line; only summary numbers written | repo goes public; Adish's review on PR #14 | `scripts/compare_to_annotation.py` docstring |
| Split-half labellers use alternate genes of the sorted marker list | disjoint genes make errors more independent; halves with < 2 genes are skipped | `scripts/split_half_calls.py` |
| Pseudobulk DE (sum of raw counts per heart × chamber, non-vetoed CMs, pydeseq2 0.5.4, hearts as replicates) instead of cell-level DE | cell-level tests treat cells as replicates and inflate false discoveries | Squair et al. 2021 Nat Commun, doi:10.1038/s41467-021-25960-2; Crowell et al. 2020 Nat Commun, doi:10.1038/s41467-020-19894-4; owner directive 15:11 |
| Pseudobulks with < 30 cells dropped and listed | DESeq2 size factors and dispersions are unreliable on a few cells | `results/hk/pseudobulk_dropped.csv` (Hr19 atrium only) |
| Whole-heart libraries split by the per-cell myosin chamber call | the only chamber annotation available; 57–87% accurate, so whole-heart chamber pseudobulks are mixed | `results/hk/va_accuracy.tsv` |
| Local injury effect = one explicit `vent_injured` column in `~heart + chamber + vent_injured` | formulaic's chamber:condition coding adds a column per chamber level that is collinear with heart (condition is nested in heart) and the fit was singular | `scripts/pseudobulk_de.py` |
| DEG cleaning: drop dissociation/heat-shock genes (control axis) and the 106 non-CM ambient genes from zf.confounds; padj < 0.05, \|log2FC\| > 1; top 30 up = data-derived injury set | owner directive 16:20 | `results/hk/pseudobulk_de_summary.csv` |
| Thinned re-run: every cell binomially thinned to 600 UMIs before summing | tests whether the DEGs need depth | `results/hk/pseudobulk_thinning_robustness.csv` |

## How to run

```shell
uv sync
scripts/get_data.sh        # downloads the public counts into data/
uv run pytest
```

## Repository layout

| Path | Contents |
|------|----------|
| `src/zf/` | Python package: loading, scoring, confound tests, evaluation |
| `scripts/` | Command-line entry points and the data download |
| `config/` | Dataset manifest, state axes with expected directions, pre-registered checks, loop settings |
| `markers/` | Marker gene sets, one YAML per state, each gene with a source |
| `literature/` | Reading list and notes on each paper |
| `docs/` | Problem statement and the original problem slide |
| `results/` | Tables and figures, each with the command and commit that made it |
| `tests/` | Tests, written before the code they test |
| `notebooks/` | Exploration only; nothing in `results/` comes from a notebook |

Team conventions are in [CONTRIBUTING.md](CONTRIBUTING.md), and the task list is in [TASKS.md](TASKS.md).

## License

Code is licensed under [Apache 2.0](LICENSE). Text is licensed under [CC BY 4.0](LICENSE-CC-BY-4.0).
