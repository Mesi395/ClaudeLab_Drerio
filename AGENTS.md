# AGENTS.md

Rules for any coding agent working in this repository.
Only what cannot be read from the code goes here.
Team conventions for humans are in [CONTRIBUTING.md](CONTRIBUTING.md).

## Never read the held-out annotation

`holdout/` holds the expert annotation and the list of known traps.
Do not open, list, grep or load anything under it.
Only `src/zf/evaluate.py` reads it, and only a human runs that script.
If you see its contents, the agreement check is void: say so in `JOURNAL.md`.

## Environment

- Python through uv only: `uv run python ...`, `uv run pytest`, `uv add <pkg>`. Never `pip install`.
- Lint and format with ruff before every commit: `uv run ruff check --fix . && uv run ruff format .`
- Data lives in `data/` and is never committed.

## How to work

- Write a failing test in `tests/` before the code that makes it pass.
- Run each atomic step in its own subagent; keep the main session for planning and review.
- Take the top open task from `TASKS.md`, do it, then update its status.
- Append every result and decision to `JOURNAL.md`, with the time and links to the output files.
- Keep `PROMPT.md` current, so a new session can continue the goal from it.
- Do not commit in auto mode. Leave changes in the working tree for a human to review.

## Results

- Every file in `results/` comes from a script in `scripts/` or `src/zf/`, never from a notebook.
- Each result file records the commit hash, the command and the input file, in the file or in a sidecar `.json`.
- Never type a number into prose. Quote it from a file in `results/` and link the file.
- Report a check that could not run as "not run", never as "passed".

## Markers

- Every marker gene in `markers/*.yaml` carries a source: a DOI, or a result file from this repository.
- Use zebrafish gene symbols as they appear in the data (`nppa`, not `NPPA`).
- A marker whose score tracks library size, doublet score or ambient RNA is flagged, not silently dropped.
