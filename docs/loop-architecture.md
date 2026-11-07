# Research loop architecture

Status: spec for review. Nothing here is built yet except where a file is named as existing.

## Principle

Language models propose; code decides.
An LLM may propose a marker gene, judge a citation or explain a flag.
Only deterministic code computes a number, decides a gate, or accepts a change.
No agent reads `holdout/`.

## Parts

| Part | Kind | Input | Output | Decides |
|------|------|-------|--------|---------|
| Planner | main Claude Code session under `/goal` | `PROMPT.md`, `TASKS.md`, `JOURNAL.md`, latest run manifest | one step at a time, each sent to a subagent | the order of steps |
| Annotators (2–3) | LLM subagents, isolated contexts | literature (PubMed), spatial tables (`markers/sources/`) | `proposals/<annotator>/<axis>.yaml` in the markers schema, with DOIs | nothing |
| Citation check | code, then an LLM judge | each proposed gene and its DOI | `supported`, `unsupported` or `unresolvable`, with a quote | nothing; feeds the gate |
| Scorer and checks | code: `zf.score`, `zf.checks` | data, marker sets, `config/checks.yaml` | per axis: statistic and outcome per check | the outcomes |
| Gate | code | check outcomes | per axis: `trusted`, `flagged` or `not_run` | the gate |
| Critic | LLM subagent | run manifest, `JOURNAL.md`, train-heart results only | flags sorted by cause; one proposed change; escalation items | nothing |
| Hillclimb controller | code, rules in `config/hillclimb.yaml` | a proposed change, train and test results | keep or revert | accepting changes |

## One round

1. The planner reads the state files and picks the top task.
2. Each annotator proposes marker sets for every axis, without seeing the others.
3. The citation check runs. Unsupported genes are dropped from the proposal and logged.
4. Code merges the proposals into candidate sets. It computes disagreement per axis: Fleiss' kappa of cell calls and Jaccard of gene sets.
5. Code scores the cells and runs every check on the train hearts.
6. The gate assigns each axis an outcome. Checks run in order: signal, then confound, then robustness. If a signal check flags, the later checks on that axis are reported `not_run`. A confound flag on an axis with no signal cannot be interpreted.
7. The critic reads the train results and proposes one change, with a citation or a result file.
8. Code reruns the checks with the change on train and test hearts. The controller keeps or reverts it.
9. Code writes the run manifest. The planner appends the round to `JOURNAL.md` and updates `TASKS.md` and `PROMPT.md`.

## The gate

- An axis is `trusted` only if every gating check in `config/checks.yaml` passes on the test hearts.
- Outcomes are `pass`, `flag`, `within_noise`, `unstable` or `not_run`. Only `pass` counts as passed.
- The `/goal` condition must be checkable from files, for example:
  every axis in `config/axes.yaml` has an outcome for every check in the latest manifest, and `scripts/verify.py` exits 0 or 2.
- No LLM decides a gate.

## LLM evaluators

| Evaluator | Model rule | Rubric | Stability |
|-----------|-----------|--------|-----------|
| Annotators | at least two different models; different evidence per annotator (literature only, spatial only, recent papers only) | output must validate against the markers schema | not applicable |
| Citation judge | a different model from every annotator | checkable claims: the passage names the gene; it links the gene to the state; it states species and injury model | run twice; a changed verdict marks the item `unstable` |
| Critic | a different model from every annotator | each cause must cite a file and a value from the manifest | run twice; a changed cause is an escalation item |

Disagreement between annotators never overrides a failed check. It only ranks what goes to a human.

## Escalation

Code writes `results/runs/<id>/escalations.yaml` when:

- annotator agreement on an axis is below the threshold in `config/annotation.yaml` (to be written);
- an axis passes under one annotator's set and flags under another's;
- a gene is proposed for two different axes;
- the cells that change call between annotators come mostly from one library;
- a cited source does not support its gene.

The queue is capped and ranked by the number of cells affected. Today it is reported, not resolved.

## Self-improvement

**Inner loop: the marker sets.** The surface is `markers/*.yaml`, core genes only. One change per round. Keep it only if train and test hearts both improve. Revert if only train improves or either gets worse. When the objective stalls for two rounds, make no change: sort every remaining flag by check and write the table to `JOURNAL.md`.

**Outer loop: the checks.** `zf.synthetic` plants one known defect at a time, one for each gating check in `config/checks.yaml`. The detection matrix over many seeds measures each check's detection and false-alarm rate. A check that misfires is fixed and the matrix is rerun. Check thresholds change only by a human commit to `config/checks.yaml`.

## What agents check, and what they never do

Agents check:

- whether a paper supports a gene–state link, with a quote, and whether its species and injury model match ours;
- whether annotators disagree, and on which genes;
- what caused each flag, and which single change could fix it.

Agents never compute a number, decide a gate, edit `config/`, or read `holdout/`. The write guard in `.claude/settings.json` enforces the last two.

## State in files

| File | Writer | Content |
|------|--------|---------|
| `PROMPT.md` | planner | the goal and how to continue |
| `TASKS.md` | planner | prioritised tasks |
| `JOURNAL.md` | planner | every round: change, outcome, decision, links |
| `proposals/<annotator>/` | annotators | proposed marker sets |
| `results/runs/<id>/manifest.yaml` | code | inputs, versions, outcomes |
| `results/runs/<id>/escalations.yaml` | code | items for a human |

## Open questions for review

1. Should the orchestration be prompts and subagents in one Claude Code session, or a script that drives the Agent SDK?
2. Where do annotator model choices live: `config/annotation.yaml` or `config/loop.yaml`?
3. How do we keep one round under ten minutes on all cells? Options: subsample the train hearts, or cache scores per gene.
4. What is the smallest end-to-end demo by 16:15: one round with two annotators and no hillclimb?
