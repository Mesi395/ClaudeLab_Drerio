# Methods and tools sheet: scoring cell states with no ground truth

On 2026-11-07 I opened every link below, or matched it to a PubMed ID. "PMID" means PubMed returned that exact citation.

## 1. Gene-set scoring
- **scanpy `tl.score_genes`**: [docs](https://scanpy.readthedocs.io/en/stable/generated/scanpy.tl.score_genes.html). It returns the mean of the set minus the mean of control genes sampled from 25 expression bins (default `ctrl_size=50`). This gives the default axis score. Tirosh 2016, PMID [27124452](https://pubmed.ncbi.nlm.nih.gov/27124452/).
- **AUCell**: rank-based, so it is less sensitive to library size. In Python it is `decoupler.mt.aucell` or [pySCENIC](https://github.com/aertslab/pySCENIC). Aibar 2017, PMID [28991892](https://pubmed.ncbi.nlm.nih.gov/28991892/).
- **UCell / pyUCell**: a Mann-Whitney U rank score, installed with `pip install pyucell`. [pyucell](https://github.com/carmonalab/pyucell), [R](https://github.com/carmonalab/UCell). Andreatta 2021, PMID [34285779](https://pubmed.ncbi.nlm.nih.gov/34285779/).
- **VISION**: an R package. It tests signatures against random signatures matched for size and expression. Copy that null. [GitHub](https://github.com/YosefLab/VISION); DeTomaso 2019, PMID [31558714](https://pubmed.ncbi.nlm.nih.gov/31558714/).
- **decoupler**: one API for `ulm`, `mlm`, `aucell`, `zscore`, `gsva`, `ora` and `consensus`. `op.translate_net` maps a human gene set to orthologs, which matters for zebrafish. [docs](https://decoupler.scverse.org/en/latest/); Badia-i-Mompel 2022, [doi:10.1093/bioadv/vbac016](https://doi.org/10.1093/bioadv/vbac016).
- **Library-size caveat**: Noureen 2022 found that bulk-designed scores (ssGSEA, GSVA) are biased by dropout and library size, and single-cell methods less so. [eLife 71994](https://elifesciences.org/articles/71994), PMID [35212622](https://pubmed.ncbi.nlm.nih.gov/35212622/).

## 2. Confound tools
**Doublets**
- **scrublet**: a pure-Python port is built into scanpy as `pp.scrublet` ([docs](https://scanpy.readthedocs.io/en/stable/api/generated/scanpy.pp.scrublet.html), [original](https://github.com/swolock/scrublet)). Wolock 2019, PMID [30954476](https://pubmed.ncbi.nlm.nih.gov/30954476/).
- **scDblFinder**: an R package. Its paper reports that an independent benchmark found it outperformed the alternatives. [Bioconductor](https://bioconductor.org/packages/scDblFinder); Germain 2021, [doi:10.12688/f1000research.73600.2](https://doi.org/10.12688/f1000research.73600.2).
- **DoubletFinder**: R, Seurat-based. [GitHub](https://github.com/chris-mcginnis-ucsf/DoubletFinder); McGinnis 2019, PMID [30954475](https://pubmed.ncbi.nlm.nih.gov/30954475/).
- **solo**: `scvi.external.SOLO`. [docs](https://docs.scvi-tools.org/en/stable/api/reference/scvi.external.SOLO.html); Bernstein 2020, PMID [32592658](https://pubmed.ncbi.nlm.nih.gov/32592658/).
- **Benchmark**: Xi & Li 2021 found that DoubletFinder had the best detection accuracy. Cell Syst, PMID [33338399](https://pubmed.ncbi.nlm.nih.gov/33338399/).

**Ambient RNA**
- **SoupX**: R. [GitHub](https://github.com/constantAmateur/SoupX); Young 2020, GigaScience.
- **CellBender**: needs the raw droplet matrix, and a GPU is advisable. [docs](https://cellbender.readthedocs.io/en/latest/); Fleming 2023, PMID [37550580](https://pubmed.ncbi.nlm.nih.gov/37550580/).
- **DecontX**: R, in the [celda](https://bioconductor.org/packages/celda) package. Yang 2020, PMID [32138770](https://pubmed.ncbi.nlm.nih.gov/32138770/).
- **Benchmark**: Janssen 2023 measured background noise at 3-35% of UMIs per cell. Noise scaled with marker specificity. CellBender gave the most precise estimates, and clustering was fairly robust. Our inference: marker-based state scores carry more ambient risk than clusters do. [doi:10.1186/s13059-023-02978-x](https://doi.org/10.1186/s13059-023-02978-x).

**Normalisation**
- **sctransform**: R. [GitHub](https://github.com/satijalab/sctransform); Hafemeister 2019, PMID [31870423](https://pubmed.ncbi.nlm.nih.gov/31870423/).
- **scran pooling**: R. [Bioconductor](https://bioconductor.org/packages/scran); Lun 2016, PMID [27122128](https://pubmed.ncbi.nlm.nih.gov/27122128/).
- **Analytic Pearson residuals**: available in Python. [scanpy.experimental](https://scanpy.readthedocs.io/en/stable/generated/scanpy.experimental.pp.normalize_pearson_residuals.html); Lause 2021, PMID [34488842](https://pubmed.ncbi.nlm.nih.gov/34488842/).

**Batch and library effects**
- **Harmony**: [harmonypy](https://github.com/slowkow/harmonypy), run through [`sc.external.pp.harmony_integrate`](https://scanpy.readthedocs.io/en/stable/generated/scanpy.external.pp.harmony_integrate.html). Korsunsky 2019, PMID [31740819](https://pubmed.ncbi.nlm.nih.gov/31740819/).
- **scVI**: [docs](https://docs.scvi-tools.org/en/stable/); Lopez 2018, PMID [30504886](https://pubmed.ncbi.nlm.nih.gov/30504886/).
- **Integration benchmark**: Luecken 2022, Nat Methods, PMID [34949812](https://pubmed.ncbi.nlm.nih.gov/34949812/). For metrics, see [scib-metrics](https://github.com/YosefLab/scib-metrics).

## 3. Evaluating annotations without ground truth
- **Agreement**: use ARI, NMI and Cohen's kappa from scikit-learn ([kappa](https://scikit-learn.org/stable/modules/generated/sklearn.metrics.cohen_kappa_score.html)). Add a per-class confusion matrix.
- **Matched random gene-set null**: scDRS samples control sets matched on mean and variance, and it is the cleanest template. [GitHub](https://github.com/martinjzhang/scDRS); Zhang 2022, PMID [36050550](https://pubmed.ncbi.nlm.nih.gov/36050550/).
- **Circularity**: do not cluster and then test markers on the same counts. Count splitting gives valid inference. [countsplit](https://github.com/anna-neufeld/countsplit); Neufeld 2024, [doi:10.1093/biostatistics/kxac047](https://doi.org/10.1093/biostatistics/kxac047).
- **Is a cluster real**: sc-SHC tests cluster significance. [GitHub](https://github.com/igrabski/sc-SHC); Grabski 2023, PMID [37429993](https://pubmed.ncbi.nlm.nih.gov/37429993/).
- **Stability under resampling**: chooseR (Patterson-Cross 2021, PMID [33522897](https://pubmed.ncbi.nlm.nih.gov/33522897/)) and [clustree](https://github.com/lazappi/clustree) (Zappia 2018, PMID [30010766](https://pubmed.ncbi.nlm.nih.gov/30010766/)).
- **Tests to code yourself**: these need no package.
  - partial correlation of each score with log UMIs, mitochondrial fraction and doublet score, by library
  - leave-one-gene-out and 80% gene-subsample bootstraps of each set
  - label permutation within a library
  - a negative-control axis, such as a housekeeping set, which should not separate states
  - agreement across scoring methods (score_genes, AUCell, UCell, ULM)
- **Coverage gap**: I did not run a dedicated search for a general framework paper on judging annotation without ground truth. The nearest one I found is the LICT marker-credibility test (section 4).

## 4. LLM and agent annotation: prior art
- **GPTCelltype**: Hou & Ji 2024. The paper names four risks:
  - "undisclosed" training data
  - fine-tuning that hurts reproducibility
  - noisy differentially expressed genes
  - "AI hallucination"

  It recommends that experts validate the output. [doi:10.1038/s41592-024-02235-4](https://doi.org/10.1038/s41592-024-02235-4); [R pkg](https://github.com/Winnie09/GPTCelltype).
- **CASSIA**: a multi-agent annotator with a validator agent and optional RAG. Nat Commun 2025. [paper](https://www.nature.com/articles/s41467-025-67084-x), [GitHub](https://github.com/ElliotXie/CASSIA/).
- **LICT**: combines several models with a "talk-to-machine" loop. It scores an annotation as reliable when at least 4 suggested markers pass an expression threshold, which is a credibility test worth copying. Commun Biol 2025. [paper](https://www.nature.com/articles/s42003-025-08745-x), [GitHub](https://github.com/Glowworm-cell/LICT).
- **CellAgent**: a multi-agent pipeline for scRNA-seq analysis. [arXiv 2407.09811](https://arxiv.org/abs/2407.09811).
- **scAgent**: a universal LLM agent for cell annotation. [arXiv 2504.04698](https://arxiv.org/abs/2504.04698).
- **Biomni**: a general biomedical agent. [site](https://biomni.stanford.edu/), [GitHub](https://github.com/snap-stanford/Biomni), [bioRxiv](https://www.biorxiv.org/content/10.1101/2025.05.30.656746v1), and a Science version (doi 10.1126/science.adz4351, which I saw in search but did not open).
- **CellWhisperer**: a transcriptome-text embedding with chat, trained on human data. Nat Biotechnol 2025. [paper](https://www.nature.com/articles/s41587-025-02857-9).
- **Benchmark**: Xiao 2025 compared 7 LLMs with 3 classical tools on 34 human and mouse datasets. Markers ranked by log2 fold change worked best, and an ensemble beat single models. Brief Bioinform, [doi:10.1093/bib/bbaf677](https://doi.org/10.1093/bib/bbaf677). None of them is zebrafish.
- **Known failure modes**: invented markers, and cherry-picking real markers to fit a prior label. The second claim comes only from a vendor blog ([Nygen](https://www.nygen.io/resources/blog/cherry-picking-vs-hallucination)), not a paper. The guard against both is to check that every marker the LLM cites is in the data and enriched.

## 5. Claude-specific
- **Claude Science** ([overview](https://claude.com/docs/claude-science/overview)): a beta desktop app (macOS, Windows, Linux) on Pro, Max, Ultra, Team or Enterprise plans.
  - It runs Python, R or shell in a sandbox and saves versioned artifacts with provenance.
  - A background reviewer checks claims against the execution record. It does not re-run analyses.
  - It needs about 5 GB of disk.
- **`/goal`** ([docs](https://code.claude.com/docs/en/goal)): sets one completion condition per session. An evaluator, built on hooks, judges each turn. Use auto mode for unattended runs. Example: "/goal every axis has null, confound and stability results in results.json".
- **`/loop`** ([docs](https://code.claude.com/docs/en/scheduled-tasks)): repeats a prompt on an interval, or self-paces.
- **K-Dense scientific skills**: the repository is now [scientific-agent-skills](https://github.com/K-Dense-AI/scientific-agent-skills). It has 177 skills under the MIT licence, including Scanpy, AnnData and scvi-tools, but not decoupler. Install with `npx skills add K-Dense-AI/scientific-agent-skills`.

## Recommended minimal stack (laptop, pip only)
1. `pip install scanpy anndata scikit-learn` for QC, `pp.scrublet` doublets, Pearson residuals, `tl.score_genes`, and ARI, NMI and kappa.
2. `pip install decoupler pyucell` for ULM, AUCell and UCell scores, a consensus across them, and `translate_net` for zebrafish orthologs.
3. `pip install harmonypy` for a per-library correction check. Score on uncorrected data, then use Harmony only to test whether states survive the correction.
4. Code three nulls in about 100 lines of numpy. Use scDRS-style matched random sets, gene-subsample bootstraps, and per-library partial correlation with log UMIs.
5. Optional: `pip install scvi-tools` for SOLO doublets. Leave out CellBender and SoupX unless the team has raw droplet matrices and R. Use Janssen 2023 as the reason marker scores need an ambient check.
