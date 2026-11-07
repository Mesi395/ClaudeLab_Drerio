# hk corpus: papers behind markers v1

Titles and DOIs only. No abstracts in this repo. Every DOI below was resolved against the Crossref API on 2026-11-07 (title, first author and year matched); none failed.

## How the Zotero export dedupes

The team export `HR_1.csv` has 535 rows. Measured on 2026-11-07:

| count | value |
|---|---|
| rows | 535 |
| unique exact titles | 384 |
| unique titles after lower-casing and stripping punctuation | 379 (the team note says 380; the difference is one near-duplicate title) |
| unique DOIs | 329 |
| unique works after DOI-then-title merge (ranking pipeline) | 359 |

The rows repeat because the same paper sits in several Zotero collections. The ranking pipeline (BM25 over title + abstract + tags, five query sets: injury, border zone, dedifferentiation, proliferation, confounds, then a hand check of the top ~60) ran locally on the 359 works. Full texts were then read only for the papers marked "full text" below.

## The 8 must-read papers

| paper | year | DOI | in export | read for v1 | why |
|---|---|---|---|---|---|
| Wu et al., Dev Cell: Spatially resolved genome-wide transcriptional profiling identifies BMP signaling as essential regulator of zebrafish cardiomyocyte regeneration | 2016 | https://doi.org/10.1016/j.devcel.2015.12.010 | yes | no (supplement download blocked) | tomo-seq border-zone genes; v1 uses crip2, vwa1, psmd14 and nppb from the team draft, marked unverified |
| Honkoop et al., eLife: Single-cell analysis uncovers that metabolic reprogramming by ErbB2 signaling is essential for cardiomyocyte proliferation in the regenerating heart | 2019 | https://doi.org/10.7554/eLife.50163 | yes | full text + source data (cluster 7 vs cluster 2 table) | nppa+ sorted border-zone CMs; glycolysis ISH; most of the BZ and dedifferentiation claims |
| van Duijvenboden et al., Circulation: Conserved NPPB+ border zone switches from MEF2- to AP-1-driven gene program | 2019 | https://doi.org/10.1161/CIRCULATIONAHA.118.038944 | yes | no (no free copy found) | AP-1 genes are part of the real BZ program, so dissociation IEGs cannot simply be called artefact |
| Beisaw et al., Circ Res: AP-1 contributes to chromatin accessibility to promote sarcomere disassembly and cardiomyocyte protrusion during zebrafish heart regeneration | 2020 | https://doi.org/10.1161/CIRCRESAHA.119.316167 | yes | full text | gata4:EGFP as dedifferentiation reporter; runx1 named as dedifferentiation marker |
| Hu et al., Nat Genet: Origin and function of activated fibroblast states during zebrafish heart regeneration | 2022 | https://doi.org/10.1038/s41588-022-01129-5 | yes | full text (author manuscript) + supplementary tables | the dataset; ttn.1/ttn.2/bves/synpo2lb cluster genes; non-CM marker tables for the ambient list |
| Constanty et al., Nat Commun: Border-zone cardiomyocytes and macrophages regulate extracellular matrix remodeling to promote cardiomyocyte protrusion during heart regeneration | 2025 | https://doi.org/10.1038/s41467-025-59169-4 | yes | full text | names mustn1b, tagln, myl6 as BZ CM markers; pkma as glycolytic dedifferentiation marker |
| Bouwman et al., Nat Cardiovasc Res: Cross-species comparison reveals that Hmga1 reduces H3K27me3 levels to promote cardiomyocyte proliferation and cardiac regeneration | 2025 | https://doi.org/10.1038/s44161-024-00588-9 | yes | no | identified as "Bouwman 2025" by author match in the export (our identification); not used for v1 genes |
| van den Brink et al., Nat Methods: Single-cell sequencing reveals dissociation-induced gene expression in tissue subpopulations | 2017 | https://doi.org/10.1038/nmeth.4437 | yes | supplementary tables (via r2 controls) | dissociation IEG control axis |

## Other papers used for v1 genes

| paper | year | DOI | in export | read | used for |
|---|---|---|---|---|---|
| Li et al., Nat Commun: An organ-wide spatiotemporal transcriptomic and cellular atlas of the regenerating zebrafish heart | 2025 | https://doi.org/10.1038/s41467-025-59070-0 | yes | full text + Supplementary Data 2 | injury axis (activated BZ domain, Stereo-seq, no dissociation); BZ domains; pcna/cdk1 |
| Neuschulz et al., Mol Syst Biol: A single-cell RNA labeling strategy for measuring stress response upon tissue dissociation | 2022 | https://doi.org/10.15252/msb.202211147 | yes | datasets EV1/EV3/EV4 (via r2 controls) | dissociation control axis |
| Lekkos et al., Nat Cardiovasc Res: Oxidative phosphorylation is required for cardiomyocyte re-differentiation and long-term fish heart regeneration | 2025 | https://doi.org/10.1038/s44161-025-00718-x | yes | full text | mcm5 proliferation marker; hk1 early glycolysis |
| Kikuchi et al., Nature: Primary contribution to zebrafish heart regeneration by gata4+ cardiomyocytes | 2010 | https://doi.org/10.1038/nature08804 | yes | full text | gata4 |
| Jopling et al., Nature: Zebrafish heart regeneration occurs by cardiomyocyte dedifferentiation and proliferation | 2010 | https://doi.org/10.1038/nature08899 | yes | full text | plk1 |
| Fukuda et al., EMBO Rep: Stimulation of glycolysis promotes cardiomyocyte proliferation after injury in adult zebrafish | 2020 | https://doi.org/10.15252/embr.201949752 | no | full text | pkma |
| Ellman et al., Cells: Identifying transcription factors distinguishing regenerating cardiomyocytes in the zebrafish heart | 2026 | https://doi.org/10.3390/cells15161455 | no | full text | ccnb1, aurkb (mouse CMs) |
| Tirosh et al., Science: Dissecting the multicellular ecosystem of metastatic melanoma by single-cell RNA-seq | 2016 | https://doi.org/10.1126/science.aad0501 | no | no (list as redistributed in Seurat) | cell-cycle genes, ortholog mapping ours |
| Ogawa et al., Science: Krüppel-like factor 1 is a core cardiomyogenic trigger in zebrafish | 2021 | https://doi.org/10.1126/science.abe2762 | no | no | klf1 (team draft, unverified) |
| Tsedeke et al., Dev Biol: Cardiomyocyte heterogeneity during zebrafish development and regeneration | 2021 | https://doi.org/10.1016/j.ydbio.2021.03.014 | yes | no | tnnc2 (team draft, unverified) |
| Boskovic et al., Cells Dev: ankrd1a is consistently upregulated in cardiomyocytes bordering the injury or scar area | 2026 | https://doi.org/10.1016/j.cdev.2026.204082 | yes | title only | ankrd1a (unverified) |
| Fang et al., Front Cell Dev Biol: Tbx20 induction promotes zebrafish heart regeneration by inducing cardiomyocyte dedifferentiation | 2020 | https://doi.org/10.3389/fcell.2020.00738 | yes | title only | tbx20 (supporting) |
| de Sena-Tomás et al., Nat Commun: Activation of Nkx2.5 transcriptional program is required for adult myocardial repair | 2022 | https://doi.org/10.1038/s41467-022-30468-4 | yes | title only | nkx2.5 (supporting) |

## Gaps

- The export has no paper on doublet detection or ambient-RNA removal (title/abstract search: zero hits). Those method citations have to come from outside the corpus.
- Wu 2016 supplementary tables and van Duijvenboden 2019 were not read: every Wu-derived gene is marked unverified.
