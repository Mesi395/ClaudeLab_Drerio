# Marker genes per cardiomyocyte state

Each citation was checked against its PubMed or Europe PMC record on 2026-11-07.
**FT**: full text searched. **Abs**: abstract only. **Unverified**: no primary support found.
Gene symbols follow GRCz11 release 92, the annotation of the GEO count files.
The machine-readable sets are in [`../markers/`](../markers/).

## Cautions that apply to every state

- The four states share genes in the literature (*nppa*, *nppb*, *mustn1b*, glycolysis genes).
  Shared genes build co-occurrence into the scores, so the sets in `markers/` are disjoint.
- Hu et al. (2022) report about 1,600 transcripts per cardiomyocyte against about 4,500 for other cells.
  Very abundant cardiomyocyte genes such as *ttn.1* and *ttn.2* are the likeliest to follow library size.
  This is our hypothesis, not a published finding.
- Low *tnnt2a* is not a dedifferentiation marker: Honkoop et al. (2019, Fig. 2h) report it re-expressed in border-zone cells.
- Several papers use ventricular resection (Kikuchi 2010; Li 2025) or amputation (Jopling 2010), not cryoinjury.

## Injury-responsive

| Gene | Evidence | Source | Caveat |
|------|----------|--------|--------|
| *tnfrsf11b, hspb6, hspb11* | Stereo-seq, scRNA, FISH | Li et al. (2025) FT | Resection model |
| *jun, junba, junbb, fosab, fosl1a* | ATAC-seq, RT-qPCR, scRNA | Beisaw et al. (2020) FT; Constanty et al. (2025) FT | Dissociation induces these genes |
| *klf1* | Myocardial loss and gain of function | Ogawa et al. (2021) Abs | Likely low counts |
| *stat3, hmox1, atf3, klf2* | — | Unverified for zebrafish cardiomyocytes | *klf2* is endocardial (Rasouli et al., 2018) |

## Border-zone-like

| Gene | Evidence | Source | Caveat |
|------|----------|--------|--------|
| *nppa* | Reporter, FACS, SORT-seq, FISH | Honkoop et al. (2019) FT | Partly overlaps dedifferentiated cells at 7 dpi |
| *nppb, mustn1b* | ISH, scRNA | Honkoop et al. (2019) FT | *nppb* also tops a dedifferentiated domain (Li 2025) |
| *crip2, nppb, myl6, vwa1, psmd14* | Tomo-seq border-zone peak | Wu et al. (2016), figure on the problem slide | *myl6* is also high in non-cardiomyocytes |
| *tagln, myl6, ctgfa, cyr61* | scRNA differential expression | Constanty et al. (2025) FT | High in fibroblasts and smooth muscle |
| *vmhc* | ISH | Honkoop et al. (2019) FT | Absent from release 92; ask which symbol stands for it |
| Mammalian *Nppa, Nppb, Ankrd1, Xirp2* | Mouse and human spatial data | Calcagno et al. (2022) Abs; van Duijvenboden et al. (2019) Abs; Kuppe et al. (2022) FT | Zebrafish paralog choice unverified |

## Dedifferentiated

| Gene | Evidence | Source | Caveat |
|------|----------|--------|--------|
| *ttn.1, ttn.2, bves, synpo2lb* | 10x scRNA; *ttn.2* FISH | Hu et al. (2022) FT | Top library-size suspects |
| *myom1b* | ISH | Honkoop et al. (2019) FT | Seen at 7 dpi, absent at 3 dpi |
| *tnnc2* | Reporter line | Tsedeke et al. (2021) Abs | |
| *gata4* | Reporter and Cre lineage tracing | Kikuchi et al. (2010) Abs | Resection model; low mRNA |
| *hk1, pkma, pdk2a, ldha* | ISH, glucose uptake | Honkoop et al. (2019) FT; Fukuda et al. (2020) FT | Lekkos et al. (2025) report oxidative phosphorylation highest in border-zone cells; direction contested |
| Embryonic myosin | Immunostaining | Fukuda et al. (2020) FT; Lekkos et al. (2025) FT | Protein readout |

## Proliferating

| Gene | Evidence | Source | Caveat |
|------|----------|--------|--------|
| *pcna, cdk1* | Stereo-seq; PCNA immunostaining | Li et al. (2025) FT | Polyploid cells can score S phase without dividing |
| *mcm5* | scRNA reanalysis | Lekkos et al. (2025) FT | |
| *mki67, top2a, mcm2, ccnb1, cdc20, aurkb, ccna2* | Our ortholog mapping of general cell-cycle sets | Tirosh et al. (2016) | Mapping is ours; cycling non-cardiomyocytes can leak these |

## Technical confounds

- Hafemeister and Satija (2019): normalised expression should not correlate with sequencing depth, and one scaling factor fails for highly expressed genes.
- Tucker et al. (2020): in human heart single-nucleus data, they kept only cluster-specific genes to avoid genes present only in ambient RNA.
- We searched PubMed and Consensus for a zebrafish heart study that measures leakage of cardiomyocyte genes into other cells, and found none.

## References

- Beisaw, A., et al. (2020). AP-1 contributes to chromatin accessibility to promote sarcomere disassembly and cardiomyocyte protrusion during zebrafish heart regeneration. *Circulation Research, 126*(12), 1760–1778. https://doi.org/10.1161/CIRCRESAHA.119.316167
- Calcagno, D. M., et al. (2022). Single-cell and spatial transcriptomics of the infarcted heart define the dynamic onset of the border zone in response to mechanical destabilization. *Nature Cardiovascular Research, 1*(11), 1039–1055. https://doi.org/10.1038/s44161-022-00160-3
- Constanty, F., et al. (2025). Border-zone cardiomyocytes and macrophages regulate extracellular matrix remodeling to promote cardiomyocyte protrusion during cardiac regeneration. *Nature Communications.* https://doi.org/10.1038/s41467-025-59169-4
- Fukuda, R., et al. (2020). Stimulation of glycolysis promotes cardiomyocyte proliferation after injury in adult zebrafish. *EMBO Reports, 21*, e49752. https://doi.org/10.15252/embr.201949752
- Hafemeister, C., & Satija, R. (2019). Normalization and variance stabilization of single-cell RNA-seq data using regularized negative binomial regression. *Genome Biology, 20*, 296. https://doi.org/10.1186/s13059-019-1874-1
- Honkoop, H., et al. (2019). Single-cell analysis uncovers that metabolic reprogramming by ErbB2 signaling is essential for cardiomyocyte proliferation in the regenerating heart. *eLife, 8*, e50163. https://doi.org/10.7554/eLife.50163
- Hu, B., et al. (2022). Origin and function of activated fibroblast states during zebrafish heart regeneration. *Nature Genetics, 54*(8), 1227–1237. https://doi.org/10.1038/s41588-022-01129-5
- Jopling, C., et al. (2010). Zebrafish heart regeneration occurs by cardiomyocyte dedifferentiation and proliferation. *Nature, 464*, 606–609. https://doi.org/10.1038/nature08899
- Kikuchi, K., et al. (2010). Primary contribution to zebrafish heart regeneration by gata4+ cardiomyocytes. *Nature, 464*, 601–605. https://doi.org/10.1038/nature08804
- Kuppe, C., et al. (2022). Spatial multi-omic map of human myocardial infarction. *Nature, 608*, 766–777. https://doi.org/10.1038/s41586-022-05060-x
- Lekkos, K., et al. (2025). Oxidative phosphorylation is required for cardiomyocyte re-differentiation and long-term fish heart regeneration. *Nature Cardiovascular Research.* https://doi.org/10.1038/s44161-025-00718-x
- Li, L., et al. (2025). An organ-wide spatiotemporal transcriptomic and cellular atlas of the regenerating zebrafish heart. *Nature Communications.* https://doi.org/10.1038/s41467-025-59070-0
- Ogawa, M., et al. (2021). Krüppel-like factor 1 is a core cardiomyogenic trigger in zebrafish. *Science, 372*(6538), 201–205. https://doi.org/10.1126/science.abe2762
- Rasouli, S. J., et al. (2018). The flow responsive transcription factor Klf2 is required for myocardial wall integrity by modulating Fgf signaling. *eLife, 7*, e38889. https://doi.org/10.7554/eLife.38889
- Tirosh, I., et al. (2016). Dissecting the multicellular ecosystem of metastatic melanoma by single-cell RNA-seq. *Science, 352*, 189–196. https://doi.org/10.1126/science.aad0501
- Tsedeke, A. T., et al. (2021). Cardiomyocyte heterogeneity during zebrafish development and regeneration. *Developmental Biology, 476*, 259–271. https://doi.org/10.1016/j.ydbio.2021.03.014
- Tucker, N. R., et al. (2020). Transcriptional and cellular diversity of the human heart. *Circulation, 142*, 466–482. https://doi.org/10.1161/CIRCULATIONAHA.119.045401
- van Duijvenboden, K., et al. (2019). Conserved NPPB+ border zone switches from MEF2- to AP-1-driven gene program. *Circulation, 140*, 864–879. https://doi.org/10.1161/CIRCULATIONAHA.118.038944
- Wu, C.-C., et al. (2016). Spatially resolved genome-wide transcriptional profiling identifies BMP signaling as essential regulator of zebrafish cardiomyocyte regeneration. *Developmental Cell, 36*, 36–49. https://doi.org/10.1016/j.devcel.2015.12.010
