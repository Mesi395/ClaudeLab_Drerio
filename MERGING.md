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

## Pull request body

Three headings, short (the same shape Omarchy uses):

- **What changes**: what was wrong or missing, and what this does about it.
- **Verification**: the commands you ran and what you saw. Say "not run" for anything you could not run.
- **Stacked on**: only if the base is not `main`: the parent PR number.

Close a pull request you no longer need with one line, `Closing in favour of #n`. Only the author closes their own.

## Reviews: A writes, B reviews

- Never approve or merge your own pull request, and never let your agent do it.
- After you push fixes for a "changes requested" review, re-request that review in the same hour.
  The block stays until the reviewer clears it.
- Proposed review pairs (change them in the coordination issue): Hemang ↔ Adish, Ryan ↔ Annette.
  Any change to `markers/` or `literature/`, or any biological claim: Emi (`@Mesi395`) approves.
- An agent may write a review comment.
  An approval needs a person.

## Merging

- **Squash-merge only**, one pull request at a time, in the order of the merge queue in the coordination issue.
  The reviewer who approved presses the button.
- Before you press merge: the branch is rebased on the current `main`, tests pass, and one approval from someone other than the author.
- After each merge, every other open branch rebases before its own merge.
- The merge queue in the coordination issue is the only list.
  Whoever merges or closes a pull request edits the queue in the same minute.
- `main` must pass `uv run pytest` after every merge.
  If it breaks, the person who merged fixes it or reverts, straight away.
  A revert is a normal pull request (`git revert <squash sha>`) and jumps the queue.
- Small beats perfect: when something fails, a small follow-up pull request that fixes it is the normal path, not a reason to hold the queue.

## Stacked pull requests

A stacked PR is based on another PR's branch, not on `main`, e.g. #1 ← #2 ← #3 ← #7.

- Merge the bottom one first.
- A squash-merge rewrites the parent's commits, so the child must drop them:
  `git rebase --onto origin/main <old-parent-branch> <child-branch>`, then push with `--force-with-lease`, then set the child's base to `main`.
- Do not delete the parent's branch until its child is rebased: the command above needs the old parent tip.
- Run `uv run pytest` on the rebased child before you ask for the merge.
- Merge the child only after that.
- Open new stacks at most three deep.
  A stack waits on every approval below it, and one blocked pull request (for example a biology review) blocks all the rest.

## Deadlines and the submission tag

The times are in CONTRIBUTING.md.

- 15:30: every approved pull request is merged or has a named blocker in the coordination issue.
- After 16:15 only fixes, results and docs merge.
- When the last merge is done, before 17:00, tag `main`: `git tag submission && git push origin submission`.
  Submit and demo from that tag, so a late merge cannot change what the judges see.

## Results and data

- Never commit data or anything under `holdout/`.
- A result file carries the commit hash, the command and the input file (see CONTRIBUTING.md).
- An agent writes what it found to `JOURNAL.md` through the loop owner, or in its own pull request description, never by editing another lane's file.
