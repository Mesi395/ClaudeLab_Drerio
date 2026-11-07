# Research methodology: the supervisor–researcher loop

Status: **in use**. This describes how study `cm-state-axes` was actually run on 2026-11-07
(cycles C001–C007). The progress so far is in [research-report.md](research-report.md). The
team's planned architecture is in [loop-architecture.md](loop-architecture.md); the section
"Relation to the loop architecture" at the end maps one onto the other.

## Principle

The researcher proposes; the evidence decides; the supervisor checks the evidence a decision
rests on.

- Every number comes from code with a recorded command, and a claim is no broader than its
  evidence.
- Rules for reading a result are fixed before the result exists.
- No agent reads `holdout/`.

## Members

| Member | Kind | Responsibility | Input | Output | Decides |
|--------|------|----------------|-------|--------|---------|
| **Research Supervisor** | Claude agent (skill `research-supervision`) | Holds the user's goal and the direction. Turns the biggest gap into one assignment at a time. Rechecks the claims each decision rests on with its own code. Reviews the method. | The user's request, the research map, the researcher's report, raw data and code | `assignment.md`, `review.md`, `checks/*.py` and `.out`, decisions in the research map | Continue, complete, blocked, decision needed; which work comes next; whether a claim stands |
| **Researcher** | Claude agent (skill `evidence-research`) | Does the research: reads the literature, writes code test-first, runs the analyses, keeps the artefacts reproducible, reports claims with evidence for and against. | The assignment, the research map, the sanitised data, sources | `reading-rules.md`, `report.md`, code in `src/zf/` and `scripts/`, files in `results/`, evidence entries | Methods and setup within the assignment; it *proposes* conclusions but never completes the study |
| **User / team** | humans | Set the goal and scope. Run the agreement check against the held-out annotation (`src/zf/evaluate.py`). Review the biology. | Reports, the label file | Scope decisions, biology review | Scope; whether results are accepted into `main` |

Both agents follow one shared contract: how a study moves, what counts as evidence, and which
records each one owns.

## Kinds of work

| Kind | Purpose | Needs before it starts |
|------|---------|------------------------|
| **Understand** | Study existing work and the supplied material; build the research map | Nothing |
| **Probe** | A cheap, faithful test of whether a direction deserves deeper work | Nothing, so it can run at any time |
| **Deepen** | Serious building or validation of a chosen direction | A saved direction decision that compares the alternatives, and a method review by the supervisor |

A strong Probe signal earns Deepen. A weak or negative signal sends the study back to the map to
choose again, not into repairs.

## One cycle

1. **Assign.** The supervisor writes `cycles/CNNN/assignment.md`, containing:
   - the work kind and one question;
   - the goal part and the gap it addresses;
   - the decision it informs;
   - the findings to resolve first;
   - the evidence expected.
   It does not include the hoped-for answer or the method steps.
2. **Fix the rules.** Before computing, the researcher writes `reading-rules.md`, time-stamped:
   thresholds, what counts as "above null", "tracks depth" or "pass", and the predicted
   direction of each effect. Later changes are labelled exploratory.
3. **Do the work.** The researcher writes a failing test, then the code. It runs the analysis
   and writes results with a provenance sidecar (commit, command, input), then reruns to check
   reproducibility.
4. **Report.** The researcher writes `report.md`. It holds claims with evidence IDs for and
   against, requirement coverage, deviations, weaknesses and one status:
   - `Review requested`;
   - `Checkpoint`;
   - `Blocked`;
   - `Decision needed`.
5. **Check.** The supervisor reruns the decisive claims **with separate code** in `checks/`. It
   chooses checks that could falsify the conclusion: an independent implementation, another
   null, a confound the report did not test. Each claim is marked confirmed, corrected or not
   checked.
6. **Review and decide.** In `review.md`, the supervisor writes:
   - findings with IDs (gap, why it matters, what would close it);
   - a step back (what signal this gave, and whether the direction is still the best use of the
     work);
   - the decision and the next assignment.
   Findings that cannot change a goal decision go to `follow-ups.md`.
7. **Record.** Both agents update the research map. The supervisor appends to `JOURNAL.md` and
   `PROMPT.md`.

## Rules that guard the evidence

| Rule | How it is enforced in this study |
|------|----------------------------------|
| **No exposure to the answer key** | The supplied `.h5ad` carried the authors' labels. A tested allowlist copy was made first (`src/zf/io.py::sanitize_h5ad`), reading only the allowlisted columns. `holdout/` is denied in `.claude/settings.json`. *ttn.2* is named in the source paper (Hu et al. 2022), so any *ttn.2* result is not found unaided (FU-6). |
| **Rules before results** | `reading-rules.md` in every cycle, time-stamped before the primary run |
| **Right unit** | Statistics within library or within heart, with hearts (18) counted, never pooled cells |
| **Controls** | Matched random gene sets as the null; binomial thinning to equal depth; a positive control (the predicted rise after injury); chamber-only libraries as a metadata-labelled chamber set |
| **Independent rechecks** | `research/cm-state-axes/checks/c00*_verify.py` with `.out` outputs, and evidence entries SUP-E* |
| **Provenance** | A `results/*.json` sidecar per result; byte-identical reruns; the evidence log is append-only |
| **Honest status** | A check that did not run is reported "not run", never "passed"; negative results stay in the map |

## Records

| File | Writer | Content |
|------|--------|---------|
| `research/cm-state-axes/research-map.md` | both | The goal (as an abstract, with parts G1–G5 and requirements), how far we are, existing work, experiments, directions and decision, search trail, current state |
| `cycles/CNNN/assignment.md`, `review.md` | supervisor | One assignment and its review per cycle |
| `cycles/CNNN/reading-rules.md`, `report.md` | researcher | Rules fixed in advance; the claims |
| `evidence.md` | both | Append-only: source or run, command, result, limits (R-E* researcher, SUP-E* supervisor) |
| `checks/` | supervisor | Recheck scripts and their outputs |
| `follow-ups.md` | both | Out-of-scope work and restart conditions |
| `checks-inventory.md` | researcher | Every check, who initiated it, what it caught; what an expert would add |
| `JOURNAL.md`, `PROMPT.md` | supervisor | Repository log and the restart prompt |

## Relation to the loop architecture

| [loop-architecture.md](loop-architecture.md) | This study |
|------|------|
| Planner | Research Supervisor (also decides continuation and checks claims) |
| Annotators, citation check | Researcher (literature reading, sourced marker YAMLs); no second annotator yet |
| Scorer and checks (code) | `src/zf/score.py`, `scripts/c00*_*.py` |
| Gate (code) | Pre-registered rules in `reading-rules.md`, applied by code; the supervisor's independent recheck |
| Critic | The supervisor's review and step back |
| Hillclimb controller | Not built; a direction change is a supervisor decision, recorded and justified |

The main difference is that this loop puts a second agent between the result and the decision,
and that agent reruns the decisive numbers itself. That is how the chamber confound (F8) and the
atrial labels (F10) were found ([checks-inventory.md](../research/cm-state-axes/checks-inventory.md)).

## Limits

- The supervisor's method review is its own, not independent.
- Two agents from the same model family may share blind spots. Domain and statistical review by
  humans is still needed.
- The researcher's runtime had no subagent tool, so `AGENTS.md`'s "one subagent per step" was
  not followed (FU-4).
