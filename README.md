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

### How we decide which annotation to trust, with no ground truth

Two layers.

1. **Protocol stability** (Card B, [docs/labeling-protocol.md](docs/labeling-protocol.md)).
   A call counts only if it survives count thinning, detection-matched decoy gene sets and a per-heart gate.
2. **Latent-class model** (Dawid–Skene, `scripts/compare_to_annotation.py`).
   Five imperfect labellers per axis: (a) protocol UCell arm, (b) protocol score_genes arm on thinned counts,
   (c, d) two disjoint split-half marker sets, (e) the published annotation mapped to the axis it claims.
   EM estimates each labeller's sensitivity and specificity; 20 bootstrap resamples over hearts (not cells).
   Run once on the ROG, 30 s, 32,044 cells, 22 hearts. Output: [results/hk/compare_to_annotation.txt](results/hk/compare_to_annotation.txt).

What it says (observations, one run):

| Axis | Labeller | Sensitivity [90 % heart-bootstrap] | Specificity |
|---|---|---|---|
| dedifferentiated | a UCell | 0.99 [0.99, 1.00] | 0.998 |
| | b score_genes | 0.97 [0.93, 0.98] | 1.000 |
| | c split-half 1 | 0.00 | 1.000 |
| | e published annotation | 0.43 [0.28, 0.51] | 0.75 [0.68, 0.85] |
| proliferating | a UCell | 0.99 [0.98, 1.00] | 0.999 |
| | b score_genes | 0.78 [0.74, 0.90] | 0.998 |
| | c / d split-half | 0.03 / 0.01 | 1.000 |
| | e published annotation | 0.20 [0.03, 0.42] | 0.997 |

How to read it, honestly:
- **The latent class is mostly "a and b agree".** a and b share their genes, so conditional independence fails,
  and their near-perfect numbers are partly self-agreement. The split-half sets were meant to break this,
  but with half the genes they almost never fire (c on dedifferentiated: 0 positive cells; d had < 2 genes and was skipped).
- **The published dedifferentiated cluster disagrees with every marker-based labeller:** it holds 25 % of cells
  (8,137) against a latent prevalence of 2.5 %. Its cells sit at 1.42× their own heart's median depth
  (other cells 0.91×), and its estimated sensitivity rises from 0.18–0.36 in the low/mid depth tertiles to 0.70 in the top tertile.
  So the cluster tracks library depth more than the marker programme. The published proliferating cluster is
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
  on the cell, morphology and sarcomere disassembly, protein markers (PCNA or EdU incorporation, phospho-histone H3)
  and position relative to the wound. None of these is in scRNA-seq.
- **Dissociation artefact.** The cells were dissociated with enzymes at 37 °C. Heat-shock and dissociation stress
  vary between samples and can compound the injury stress response. Our checks separate the dissociation axis
  itself (control axis), its within-heart correlation with each state (partial correlation), and chamber-specific
  signal (ventricle-minus-atrium; only the ventricle is injured). They cannot separate the **injury × dissociation
  interaction**: injured tissue may respond to dissociation differently. The fix is experimental: metabolic labelling
  of new transcripts separates dissociation-induced transcription (Neuschulz et al. 2022, *Mol Syst Biol*).
- **A truly independent labeller.** Split-half marker sets were too small to act as one; an expert annotation
  made blind to these markers, or a spatial/imaging label, would make the latent-class model identifiable.
  Emi can run `scripts/compare_to_annotation.py` on her own annotation without sharing it; it prints numbers only.

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
