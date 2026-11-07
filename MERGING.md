# Merging — for people and their agents

Several people each run their own agents on this repository.
These rules keep `main` runnable and stop agents from undoing each other's work.
They extend [CONTRIBUTING.md](CONTRIBUTING.md); where the two disagree, CONTRIBUTING.md wins.

## Before you start a task: claim it

- Comment on the coordination issue (title starts with `[coord]`):
  `claiming <task id> on <branch>`.
  Check that nobody else claimed it first.
- One task, one branch, one pull request.
  Branch name: `<initials>/<lane>-<topic>`.
- Touch only files your lane owns (table in CONTRIBUTING.md).
  To change someone else's file, ask in the pull request or the coordination issue.

## Before you push

- `git fetch origin && git rebase origin/main` — never merge `main` into your branch.
- `uv run ruff check --fix . && uv run ruff format . && uv run pytest`: all green.
- `uv.lock` conflict: never edit it by hand.
  Take `main`'s version, then run `uv lock`.
- Push only your own branch.
  After a rebase use `git push --force-with-lease`, never `--force`, and never on someone else's branch.

## Reviews: A writes, B reviews

- Never approve or merge your own pull request, and never let your agent do it.
- Proposed review pairs (change them in the coordination issue): Hemang ↔ Adish, Ryan ↔ Annette.
  Any change to `markers/` or `literature/`, or any biological claim: Emi (`@Mesi395`) approves.
- An agent may write a review comment.
  An approval needs a person.

## Merging

- **Squash-merge only**, one pull request at a time, in the order of the merge queue in the coordination issue.
- Before you press merge: the branch is rebased on the current `main`, tests pass, and one approval from someone other than the author.
- After each merge, every other open branch rebases before its own merge.
- `main` must pass `uv run pytest` after every merge.
  If it breaks, the person who merged fixes it or reverts, straight away.

## Stacked pull requests

A stacked PR is based on another PR's branch, not on `main`, e.g. #1 ← #2 ← #3 ← #7.

- Merge the bottom one first.
- A squash-merge rewrites the parent's commits, so the child must drop them:
  `git rebase --onto origin/main <old-parent-branch> <child-branch>`, then push with `--force-with-lease`, then set the child's base to `main`.
- Merge the child only after that.

## Results and data

- Never commit data or anything under `holdout/`.
- A result file carries the commit hash, the command and the input file (see CONTRIBUTING.md).
- An agent writes what it found to `JOURNAL.md` through the loop owner, or in its own pull request description, never by editing another lane's file.
