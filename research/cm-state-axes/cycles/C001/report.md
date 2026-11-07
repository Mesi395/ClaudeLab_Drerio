# C001 report: Understand (researcher)

Written 2026-11-07 15:20 CET. Status at the end.

## What was done

1. **Leakage guard (R5.3), first.** `src/zf/io.py::sanitize_h5ad` with its test written first
   (`tests/test_io_sanitize.py`, 3 passed). It reads only the allowlist via h5py; dropped columns
   are never loaded. Sanitised copy: `data/clean/hu2022_cm_allgenes_sanitized.h5ad` (R-E1).
   No label value was seen. `holdout/` was not opened.
2. **Literature**: Honkoop 2019 (first 100k of 226k chars), Li 2025 (first 100k of 285k chars of the
   Europe PMC XML, after SUP-E4), abstracts of 9 more papers from the Zotero list, Wu 2016 via a
   search summary only, ZFIN for a symbol rename (R-E2–R-E7).
3. **Data facts** from the sanitised file only: per-library depth, genes, mito, chemistry; presence
   and per-library detection of every candidate (R-E8; `results/c001_*.csv|json`).
4. Draft `markers/*.yaml` (4 files, status `candidate`, each gene with source, model, time,
   resolution, read depth and evidence ID).

## Claims

| # | Claim | For | Against / limits |
|---|---|---|---|
| C1 | The sanitised copy carries no author labels. | R-E1: HDF5 keys checked; test asserts a planted label string is absent from the output bytes. | HVG flags in `var` kept (allowed). HVG file not sanitised. |
| C2 | Only one cryoinjury **single-cell** marker source was read: Honkoop 2019 (7 dpi, 352 sorted cells). Everything else is spots (Li 2025, resection), tissue sections (Wu 2016) or reporters/abstracts. | R-E2–R-E5 | Search not exhaustive (see weaknesses). |
| C3 | Border-zone-like has the most independent support: *nppa*, *nppb*, *mustn1b* appear in ≥2 independent sources and are detected in 29–55% of CMs. | R-E2, R-E3, R-E4, R-E5; `results/c001_marker_detection.csv` | *nppa/nppb* are also natriuretic stress genes and appear in Li's activated and dedifferentiated domains too: the same genes sit on three axes. |
| C4 | The four states are **not separable by existing marker lists**: *nppa*, *nppb*, *mustn1b* are listed for injury-responsive, border-zone and dedifferentiated. | YAMLs; R-E3 | Disjoint sets would need a data-driven or negative (sarcomere-down) definition. |
| C5 | Injury-responsive has a timing mismatch: Li's "activated" genes were defined at 6 hpa; our earliest timepoint is 3 dpi. *tnfrsf11b* (the only gene checked in cryoinjury) is in 1.1% of CMs, *ankrd1a* 0.3%. | R-E3, R-E8 | Not tested whether 3 dpi cells retain the program. |
| C6 | Proliferating is the best-defined biologically but the least detectable: *cdk1* 0.3%, *foxm1* 0.02%, *cenpf* 0.06%, *pcna* 2.3% of CMs. *ccng1/ccng2* are in Li's list but carry conflicting meaning (p53 target / arrest). | R-E3, R-E8 | Generic cell-cycle sets (*mki67*, *top2a* present in data) not yet sourced; detection may be too sparse for score_genes. |
| C7 | Dedifferentiated is the thinnest: defined in the literature mostly by protein (embryonic myosin), sarcomere disassembly and silencing (down-genes), and reporters; the transcript markers overlap the border-zone set. *ttn.2* comes only from the authors' own label: circular. | R-E4, R-E6 | Ben-Yair 2019 / Sallin 2015 full texts not read. |
| C8 | Hr11 depth trap is smaller in CMs than suggested: Hr11 median 659 counts per CM is lowest, but four libraries are within ~120 counts. | `results/c001_library_qc.csv` | Medians only. |
| C9 | Lead for the Probe: in Hr11, *nppa* and *ttn.2* are detected in more cells than in control H7 despite lower depth. | `results/c001_marker_detection.csv` | Observation only; biology (3 dpi) vs ambient/library artefact not separated. Not used to select markers. |

## The user's AUCell claim (finding to resolve)

Sources read: `literature/methods-and-tools.md` cites Noureen 2022 (single-cell scores less
biased by library size than bulk scores, not unbiased). No source read claims full independence.
Reasoning (inference, untested): ranks ignore per-cell scaling, but with medians of 253–568 genes
detected per CM (`results/c001_library_qc.csv`), most genes are tied at zero, and the number of
ties depends on depth. The test stays R3.6. pyUCell implements UCell (Mann–Whitney U), not AUCell;
`decoupler` has AUCell.

## Requirement coverage

- R1.1 four candidate sets: **draft done**, each gene sourced. Two-independent-source bar met for
  border-zone (*nppa, nppb, mustn1b*) and, via shared genes, dedifferentiated; not met for
  injury-responsive (Li only, plus *ankrd1a* abstract) nor proliferating (Li + Zuppo abstract).
- R1.2 presence: done; missing `slc2a1a`; `vmhc` → `myh7` (but low detection, unresolved).
- R1.3 model/time/resolution: recorded per gene; several say "not checked".
- R5.3: sanitiser done and tested. Leakage log: JOURNAL entry.
- TASKS data-1 (load, 32,044 × 18,794 × 22): done from the sanitised file. data-2 covariates:
  done (0 missing).

## Deviations

- AGENTS.md asks for one subagent per atomic step; no subagent tool was available, so all steps
  ran in this session.
- Directly reading `.h5ad` key names (not values) to design the sanitiser.

## Weakest points of this method

- No full-text read of anything after character 100k; Li 2025 Supplementary Data 2 and 5 not
  opened; Wu 2016 only via a search summary; most others abstract-only.
- No forward-citation search (OpenAlex) from Honkoop 2019 or Li 2025; I relied on the user's Zotero
  list (82 matching papers), which may reflect the user's choices. Not searched: Beisaw 2020
  (AP-1, cryoinjury border-zone ATAC/RNA), Ma et al. single-cell regeneration atlases, Tirosh 2016
  cell-cycle sets mapped to zebrafish.
- Marker selection is by what sources name, so large-list sources (supplements) are missing.

## Next step proposed

Probe combining A and B on the border-zone axis first (strongest support), with the C9 lead and the
R3.6 depth test on Hr11 as the decisive check; in parallel, source a generic cell-cycle set
(Tirosh 2016 via ortholog mapping) for proliferating, and define dedifferentiation partly as a
sarcomere-down axis. Before scoring: open Li 2025 Supplementary Data 2/5 and Honkoop's second half.

## Status

**Checkpoint.** Missing: full reads (Li supplements, Honkoop second half, Wu 2016, Ben-Yair 2019,
Sallin 2015), a forward-citation search, a sourced generic cell-cycle set, and resolution of the
`myh7` detection anomaly. Map sections 3, 4, 5 and 7 are updated.
