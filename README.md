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
| `markers/` | Marker gene sets, one YAML per state, each gene with a source |
| `literature/` | Reading list and notes on each paper |
| `docs/` | Problem statement and the original problem slide |
| `results/` | Tables and figures, each with the command and commit that made it |
| `tests/` | Tests, written before the code they test |
| `notebooks/` | Exploration only; nothing in `results/` comes from a notebook |

Team conventions are in [CONTRIBUTING.md](CONTRIBUTING.md), and the task list is in [TASKS.md](TASKS.md).

## License

Code is licensed under [Apache 2.0](LICENSE). Text is licensed under [CC BY 4.0](LICENSE-CC-BY-4.0).
