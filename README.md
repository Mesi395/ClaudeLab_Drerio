# ClaudeLab_Drerio

An agent that checks itself when cell states have no ground truth.

Built at the [Claude Impact Lab Berlin 2026](https://github.com/xemantic/claude-impact-lab-2026)
on a problem from the [Junker lab](https://junkerlab.com/) (MDC Berlin).

## The problem

After a cryoinjury, zebrafish heart muscle cells take on several states:
injury-responsive, border-zone-like, dedifferentiated and proliferating.
Each paper marks these states with different genes, and no agreed annotation protocol exists.
So nobody knows which cells are in which state, or whether the states co-occur.
Nobody knows either whether published markers measure biology or technical artefacts such as sequencing depth.

The data is public: Hu et al. (2022), 32,044 cardiomyocytes from 22 libraries and 18 hearts.
Full statement: [docs/problem.md](docs/problem.md).

## Our approach

Without ground truth, an agent cannot be scored on its annotation, so we score its checks instead.
The loop proposes cited markers, scores the four states, and runs checks validated on planted defects.
Disagreements with human labels are reported as diagnoses, not errors.

| | What | Why it sits here |
|---|---|---|
| **Evaluation**: truth known | Planted defects: each check must catch its own defect and nothing else | We planted the defect |
| | Same verdict across seeds; cited sources exist and name the gene; reproducible runs; no heart in both train and test | Checkable by script |
| **Diagnosis**: no truth | Agreement with the held-out human labels; disagreement between agent labellers; confound flags on the real states; co-occurrence; time course | The human labels contain mistakes, so a disagreement is a case to examine |
| **Future** | Agreed state definitions and marker standards across labs; an annotation protocol measured by agreement between humans and between models; independent labels and spatial validation | Needs a community, not an afternoon |

## Results

- **Evaluation of the checks** on planted defects over 20 seeds:
  [results/detection/matrix.tsv](results/detection/matrix.tsv), [manifest](results/runs/detection-matrix/manifest.yaml) (PR #15).
- **Checks on the real data** (no labels read):
  [results/real/checks.tsv](results/real/checks.tsv), [manifest](results/runs/real-v0/manifest.yaml) (PR #15).
- **The research loop's study**, cycles C001–C007:
  [docs/research-report.md](docs/research-report.md) and [research/cm-state-axes/](research/cm-state-axes/) (PR #17).
- **Which annotation to trust**, a latent-class comparison with the published annotation, run once:
  "Answers to the two questions" in PR #18.

## The loop

- **Agents** (Claude Code): a Research Supervisor and a Researcher run cycles of assignment, probe, report and independent review.
  Method: [docs/research-methodology.md](docs/research-methodology.md); skills in `.claude/skills/` (PR #19).
- **Code decides** every number and every gate: `src/zf/score.py`, `src/zf/checks.py`. Agents propose and review.
- **State in files**: [PROMPT.md](PROMPT.md), [TASKS.md](TASKS.md), [JOURNAL.md](JOURNAL.md), one manifest per run in `results/runs/`.
- **Guardrails**: agents read only a label-stripped copy of the data ([scripts/make_agent_input.py](scripts/make_agent_input.py)).
  A write-guard hook, pre-commit and pre-push hooks, and `scripts/verify.py` block protected writes and check every run.
- Design: [docs/loop-architecture.md](docs/loop-architecture.md).

## How to run

```shell
git clone https://github.com/Mesi395/ClaudeLab_Drerio && cd ClaudeLab_Drerio
uv sync && uv run pre-commit install
scripts/get_data.sh                          # the 22 GEO libraries into data/raw/
uv run python scripts/detection_matrix.py    # evaluate the checks on planted defects
uv run python scripts/run_real.py            # run the checks on the 32,044 cells
uv run python scripts/verify.py              # exit 0 pass, 1 fail, 2 a check could not run
```

To continue the study unattended in Claude Code, start a session in the repository and run:

```
/goal continue the study in PROMPT.md until the latest manifest in results/runs/ has an outcome for every axis in config/axes.yaml and every check in config/checks.yaml, `uv run python scripts/verify.py` exits 0 or 2, and research/cm-state-axes/research-map.md names the next assignment
```

For a long computation, add a check on a timer:

```
/loop 10m if a new manifest exists in results/runs/, append its outcomes to JOURNAL.md
```

Today the loop ran as supervised cycles; `/goal` is the entry point for running it unattended.

## Next

Discovery loops are only as good as the protocols they check against.
With no agreed state definitions or annotation standards, an agent can diagnose disagreement but cannot be evaluated against truth.
Agreed protocols from scientists would turn these diagnoses into evaluations.

Next builds: a depth check between libraries, a null that can see rare states, and an independent labeller.

## Repository layout

| Path | Contents |
|------|----------|
| `src/zf/` | Python package: loading, scoring, checks, planted defects, provenance, guards |
| `scripts/` | Data download, agent input, detection matrix, real-data run, verify, hooks |
| `config/` | Dataset, cell list, state axes, pre-registered checks, loop and hillclimb rules |
| `markers/` | Marker gene sets, each gene with a source |
| `literature/` | Marker evidence and methods notes |
| `research/` | The loop's study record |
| `results/` | Tables and run manifests |
| `docs/` | Problem, data, loop architecture, methodology, report |

Team conventions: [CONTRIBUTING.md](CONTRIBUTING.md) and [MERGING.md](MERGING.md).

## License

Code: [Apache 2.0](LICENSE). Text: [CC BY 4.0](LICENSE-CC-BY-4.0).
