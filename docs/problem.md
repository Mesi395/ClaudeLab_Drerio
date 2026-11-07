# Problem statement

Transcribed from the problem slide ([problem-slide.pptx](problem-slide.pptx)) by Emese Ratkai, Junker lab.

## Can an agent annotate cell states when there is no ground truth?

After cryoinjury, zebrafish cardiomyocytes take on several states:

- injury-responsive
- border-zone-like
- dedifferentiated
- proliferating

Open questions:

- Do these states co-occur in cells?
- What do published markers really mark?
  For example, one supposed "dedifferentiation" marker turns out to mostly track library size, not biology.

## The data

- Hu et al. (2022), public zebrafish heart single-cell data (Junker lab, MDC).
- 32,044 cells × 18,794 genes from 22 libraries.
- Counts and metadata.
- A held-out expert annotation for an agreement check.
- Known traps to test the agent against.

The 22 libraries, with the cell counts shown on the slide:

| Condition | Libraries (cells) | Cells |
|-----------|-------------------|-------|
| Uninjured | H5 (806), H6 (1,211), H7 (3,824), H8a (1,033), H8v (174) | 7,048 |
| 3 dpi | Hr11 (4,063), Hr12 (1,262), Hr24 (370), Hr27 (490) | 6,185 |
| 7 dpi | Hr1 (2,612), Hr13 (350), Hr2a (1,354), Hr2b (1,635), Hr6a (927), Hr6v (1,117), Hr7a (1,500), Hr7v (789), Hr8 (1,133), Hr9 (2,677) | 14,094 |
| 30 dpi | Hr19 (115), Hr3 (2,708), Hr4 (1,894) | 4,717 |
| Total | 22 libraries | 32,044 |

dpi: days post injury.

## Task of the agent

- Find marker genes for each state itself, from literature or spatial data.
- Score every cell on each axis.
- Test the scores against technical confounds: library size, doublets and ambient RNA.
- Test whether the axes coincide or come apart.

## Questions for the humans

- With no ground truth, how do we decide which annotation to trust?
- Which checks did the agent run on its own, and which would a domain expert add?

## Spatial reference on the slide

Wu et al. (2016) profiled cryoinjured hearts by tomo-seq and found three zones.
Example genes shown on the slide:

| Zone | Genes |
|------|-------|
| Injury area | *grn2, ctsd, lygl1, tmsb4x, grn1* |
| Wound border zone | *crip2, nppb, myl6, vwa1, psmd14* |
| Uninjured myocardium | *cox7a1, gapdh, pgam2, ckma, ckmt2b* |
