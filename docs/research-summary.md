# Research summary: can an agent annotate cardiomyocyte states without ground truth?

Study `cm-state-axes`, 2026-11-07. Method: [research-methodology.md](research-methodology.md).
Full progress log: [research-report.md](research-report.md). Status: **not complete**.

## Question

After cryoinjury, zebrafish heart muscle cells are said to enter four states: injury-responsive,
border-zone-like, dedifferentiated and proliferating. No per-cell ground truth exists. We asked:

- Can published marker genes score these states in 32,044 cardiomyocytes from 18 hearts
  (Hu et al. 2022)?
- Do the scores survive technical confounds?
- Do the states co-occur in the same cells?

## How we worked

Two Claude agents worked in cycles: a **Research Supervisor**, which set each question and
re-ran the decisive numbers with its own code, and a **Researcher**, which read the literature,
wrote tested code and ran the analyses.

- The answer key was stripped before any analysis.
- Rules for reading each result were fixed before computing.
- Statistics count hearts, not cells. Every score is compared against random genes and against
  cells thinned to equal sequencing depth.
- Seven cycles: two of reading, five of quick tests.

## What we found

1. **Rank-based scoring is not safe from sequencing depth here; it is worse.**
   - UCell and an AUCell-style score follow depth more strongly than scanpy's `score_genes`, in
     22 of 22 libraries (border-zone axis |rho| 0.468 and 0.427 against 0.287,
     [check](../research/cm-state-axes/checks/c002_verify.out)).
   - Use `score_genes` with a depth-thinning control.
2. **The published "border-zone" and "injury" markers mainly mark the atrium.**
   - Atrium cells score 5–10× higher than ventricle cells of the same heart (H8 0.679 vs 0.068,
     [check](../research/cm-state-axes/checks/c003_verify.out)).
   - Hearts differ widely in their share of atrium cells, which hides any injury signal.
3. **Once depth and chamber are controlled, no marker set shows the predicted rise after injury**
   above random genes: 0 of 6 tests, for both the original and the improved Li 2025 marker sets
   ([results/c007_positive_control.csv](../results/c007_positive_control.csv)).
4. **What signal remains sits in two hearts, not in a timepoint.** Cell-cycle genes are found
   mainly in hearts Hr9 and Hr8 (7 dpi) ([check](../research/cm-state-axes/checks/c007_verify.out)).
   Why is still open.
5. **Hr11 (3 dpi, the most weakly measured sample) carries contamination and stress signals.**
   Its *nppa* excess (2.5–3.1× at equal depth) is unexplained.
6. **Whether the states co-occur cannot be decided with literature markers.** Shared genes make
   them overlap by construction, and disjoint sets look like random pairs.

## Answers to the problem's questions

| Question | Answer so far |
|---|---|
| What do published markers really mark? | Mainly chamber (atrium) and sequencing depth, not injury states |
| Do the states co-occur? | Cannot be answered with these markers |
| Which annotation can we trust without ground truth? | Only labels that survive depth, chamber and null controls. Our chamber call does: balanced accuracy 0.986 on six labelled samples ([check](../research/cm-state-axes/checks/c004_verify.out)). Our state calls do not; they are at chance. |
| Which checks did the agent run alone, and which would an expert add? | See [checks-inventory.md](../research/cm-state-axes/checks-inventory.md). The researcher caught the false comfort of method agreement and the Hr11 anomaly. The supervisor caught the embedded answer key, the chamber confound and an error in the labels. An expert would add ambient-RNA estimates from raw droplets, a second doublet method and spatial validation. |

## What is not established

- Whether the four states exist as separable programs in these cells. The tests lack power
  (4 control hearts), and their sensitivity is not validated.
- Doublets: one method's threshold failed.
- Ambient RNA: there is no raw-droplet profile.
- *ttn.2* is named in the source paper (Hu et al. 2022), so any *ttn.2* result is not found
  unaided.
- The agent labels ([results/c005_agent_labels.csv](../results/c005_agent_labels.csv)) still
  await the human-run agreement check.

## Next

- Let the data define programs in ventricle cells, with count splitting against circular
  reasoning and heart-level covariates.
- Validate the tests on a dataset with known border-zone cells (Honkoop 2019).
- Ask the expert whether chamber is the right unit for scoring.
