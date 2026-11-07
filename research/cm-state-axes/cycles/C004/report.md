# C004 report: Probe (researcher)

Written 2026-11-07 15:50 CET. Rules fixed at 15:43, before computing:
[reading-rules.md](reading-rules.md). Evidence R-E16. The numbers come from the files named.

## 1. Chamber call ([c004_chamber_validation.csv](../../../../results/c004_chamber_validation.csv))

| Rule | Genes | Balanced accuracy (6 libraries) | Mean coverage | Usable (≥0.9, ≥0.5)? |
|---|---|---|---|---|
| A | *myh6* vs *myh7* (both sourced) | 0.867 | 0.358 | No |
| B | *myh6* vs *myh7* + *myh7l* | 0.991 | 0.735 | **Yes, used** |

- Rule A fails because `myh7` is barely detected (F1). It assigns only 1.7–6.4% of the cells in
  the ventricle-only libraries.
- Rule B, per library: atrium libraries have coverage 0.57–0.80 and accuracy 0.999–1.000; ventricle
  libraries have coverage 0.75–0.87 and accuracy 0.960–1.000.
- Caveats:
  - `myh7l`'s ventricular role is unsourced. Only its former name, *vmhcl*, suggests it.
  - The rule was chosen on the same six libraries it was validated on, so 0.991 is optimistic.
  - About 29% of all cells are unassigned: 9,204 of 32,044 ([c004_probe.json](../../../../results/c004_probe.json)).
- Calls across all cells: ventricle 14,024, atrium 8,816, unassigned 9,204. F8 is confirmed as
  a real composition difference, and the chamber call now lets us control for it.

## 2. Positive control in ventricle-called cells ([c004_ventricle_positive_control.csv](../../../../results/c004_ventricle_positive_control.csv), per heart: [c004_ventricle_heart_scores.csv](../../../../results/c004_ventricle_heart_scores.csv))

All 18 hearts have at least 50 ventricle-called cells (Hr19 has 57, Hr27 60).

| Axis | Contrast | Tail AUC, difference of medians | Median-score AUC, difference | Rule verdict |
|---|---|---|---|---|
| border-zone | 7 dpi > Ctrl | 0.786, +0.002 | 0.643, +0.011 | fail (just under 0.8) |
| dedifferentiated | 3 dpi > Ctrl | 0.750, +0.001 | 0.500, −0.040 | fail |
| dedifferentiated | 7 dpi > Ctrl | 0.857, +0.001 | 0.571, −0.035 | **pass** |
| injury-responsive | 3 dpi > Ctrl | 0.812, +0.009 | 0.375, −0.158 | **pass** |
| proliferating | 7 dpi > Ctrl | 0.679, +0.006 | 0.786, +0.428 | fail |

**Claim C4-1.** In ventricular CMs, two axes pass the fixed rule, but both passes are weak. Do not
read them as validated states:

- The tail differences are 0.1–0.9 percentage points.
- Control hearts mostly sit **below** the null's 1% (border-zone tails in Ctrl are 0.000–0.001). An
  injured heart at about 1% therefore "beats" controls while being at null level. Only Hr9 is
  clearly above null on every axis (border-zone 0.034, dedifferentiated 0.022, injury-responsive
  0.041, proliferating 0.046).
- On the median score, both passing axes move the **wrong way** (−0.035, −0.158).

Compared with C003, removing chamber raised border-zone from 0.571 to 0.786 and dedifferentiated
at 7 dpi from 0.464 to 0.857, so F8 was hiding part of a time signal. The signal that remains is
close to the null level.

## 3. G4 co-occurrence ([c004_g4.json](../../../../results/c004_g4.json), [c004_g4_cooccurrence.csv](../../../../results/c004_g4_cooccurrence.csv))

Because two axes passed, G4 ran by the rule.

- Disjoint sets removed the shared genes *nppa*, *nppb*, *mustn1b* and *ankrd1a*.
  Border-zone is left with *myh7*, *hk1*, *pkma*, *pdk2a*, *ldha* and *mmp14b*; injury-responsive with
  *tnfrsf11b*, *hspb6* and *hspb11* (the full lists are in the json).
- **Every pair is "mixed":** at most 17% of hearts fall outside the null band in either direction.
  Median within-heart rho is −0.02 to 0.08.

**Claim C4-2.** In ventricular CMs, with disjoint sets, the four axes neither co-occur nor come
apart beyond random gene pairs.

This is **not** an answer to the biological question. The axes carry little signal of their own
(C002–C004), and the genes that made them overlap (*nppa*, *nppb*, *mustn1b*) were removed. The
honest G4 statement: **with these sourced sets, co-occurrence cannot be resolved.** Overlapping
sets would co-occur by construction (F7). Disjoint sets are indistinguishable from random.

## 4. Agent label file (R5.2): delivered

- File: [results/c004_agent_labels.csv](../../../../results/c004_agent_labels.csv), 32,044 rows.
- Format: `src/zf/evaluate.py` is a one-line stub that defines no format, so I defined one. It is
  recorded in `results/c004_probe.json` under `label_file_format`:
  - `cell` = obs_names of the sanitised file
  - `library`, `heart`
  - `agent_state` ∈ {border_zone_like, dedifferentiated, injury_responsive, proliferating,
    unassigned}
  - `agent_states_all` (semicolon-separated)
  - `agent_chamber` ∈ {atrium, ventricle, unassigned}
  - `score_<axis>`
- Rule: section 4 of the reading rules, fixed before computing. `holdout/` was never read.
- Counts: unassigned 30,064; injury-responsive 1,188; proliferating 407; dedifferentiated 199;
  border-zone 186. 198 cells are high on more than one axis.
- **Caveat, to state to the human who runs `evaluate.py`:** at a 1% null threshold per heart and
  axis, chance alone gives about 320 cells per axis (1% of 32,044). Border-zone, dedifferentiated
  and proliferating labels are therefore at or below chance level. Injury-responsive (about 3.7%)
  is above chance, but it concentrates in control H8 (208 cells) and in 30 dpi Hr3 (227), not at
  3 dpi. That points to chamber or library effects (these labels use all cells, not only
  ventricle-called ones), not injury. The `agent_chamber` column is the most defensible call in
  the file.

## Requirement coverage

- R5.2: done.
- G4: run; answer "cannot be resolved with these sets".
- F8: confirmed and controlled.
- Not done: doublets (R3.2); a validated positive-control set; leave-one-gene-out; chamber
  balancing in the label file.

## Weaknesses

- The chamber rule rests on an unsourced gene (`myh7l`) and on selection using the validation
  libraries.
- One thinning seed.
- The null's 1% is not calibrated per heart: control hearts sit below it. That calibration gap is
  what makes the tail AUC pass, so the AUC rule needs a calibrated threshold, for example the null's
  percentile measured in control hearts.
- No figure.

## Status

**Review requested.** Items 1–4 are done. A rerun reproduced all six result files byte-identically.
