# Working in this repository

Four people, one afternoon, one `main`.
The rules below keep `main` runnable at every moment,
so that the loop can be demoed from it at 17:15.

## Hard deadlines

| Time  | Event |
|-------|-------|
| 15:30 | Integration: every lane merges what runs, even if partial |
| 16:15 | Feature freeze: only fixes, results and docs after this |
| 16:45 | Repository made **public**; submission issue drafted |
| 17:00 | Submission issue filed in `xemantic/claude-impact-lab-2026` |
| 17:15 | Demo, 3 minutes, starts with the team name |

## Branches

- `main` is protected by convention: never push to it directly, except the lane owner of `TASKS.md` updating task status.
- One branch per task, named `<initials>/<lane>-<short-topic>`, for example `ar/score-axes` or `ah/confound-libsize`.
- Lanes: `data`, `markers`, `score`, `confound`, `eval`, `loop`, `docs`.
- A branch lives at most 45 minutes. Merge something that runs, then branch again.
- Rebase on `main` before opening a pull request: `git pull --rebase origin main`.

## Pull requests

- Small: one task from `TASKS.md`, named in the title, e.g. `[score] AUCell axis scores`.
- Another team member approves; for marker or biology changes, Eme approves.
- Squash-merge. The squash message says what changed and which result file it produced.
- A pull request that adds a result also adds the code and the command that produced it.

## Who owns which files

One writer per file avoids merge conflicts.
If you need to change a file outside your lane, ask its owner in the team chat.

| Path | Owner lane |
|------|-----------|
| `data/` (git-ignored), `src/zf/io.py` | data |
| `markers/*.yaml`, `markers/SOURCES.md` | markers |
| `src/zf/score.py` | score |
| `src/zf/confounds.py` | confound |
| `src/zf/evaluate.py`, `holdout/` (git-ignored) | eval |
| `PROMPT.md`, `TASKS.md`, `JOURNAL.md`, `.claude/` | loop |
| `README.md`, `RESOURCES.md`, submission text | docs |

## Data rules

- Never commit data. `data/` and `holdout/` are in `.gitignore`.
- The held-out expert annotation lives only in `holdout/`.
  Only `src/zf/evaluate.py` reads it, and only after the agent has written its annotation.
  The agent must not see it, or the agreement check measures nothing.
- `scripts/get_data.sh` downloads the public data, so anyone can reproduce the run.

## Results

- Every figure or table in `results/` carries, in its file or a sidecar `.json`:
  the commit hash, the command, and the input file.
- No number is typed into `README.md` or the submission by hand.
  It is copied from a file in `results/`, and the file is linked.
- Negative results go in too. The judges score honesty.

## Commits

- Imperative mood, one line: `Add library-size partial correlation per axis`.
- Commit often on your branch; the squash keeps `main` clean.
