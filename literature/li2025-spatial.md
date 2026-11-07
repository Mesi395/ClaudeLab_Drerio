# Border-zone markers from Li et al. (2025)

Li, L., Lu, M., Guo, L., et al. (2025). An organ-wide spatiotemporal transcriptomic and cellular atlas of the regenerating zebrafish heart.
*Nature Communications, 16*, 3716. https://doi.org/10.1038/s41467-025-59070-0

## What the paper provides

- Stereo-seq of amputated ventricles at eight stages, from uninjured to 28 days post amputation (dpa), with scRNA-seq at the same stages.
- Five border-zone myocardium (BZm) domains, in time order:
  - activated, at 6 hours
  - dedifferentiated, from 12 hours to 1 dpa
  - dedifferentiated-proliferating, at 3 dpa
  - proliferating, at 7 dpa
  - redifferentiated, at 14 dpa
- Supplementary Data 2: marker genes of every spatial domain and every scRNA-seq cell type.
- Processed data: [Zenodo 14991776](https://zenodo.org/records/14991776) and the [STOmics database](https://db.cngb.org/stomics/zebrafish_VRH). Raw data: CNGBdb CNP0005245 and NCBI PRJNA1233465.

## How we derive markers

`scripts/derive_li2025_markers.py` reads Supplementary Data 2 and writes [`markers/sources/li2025_bzm.tsv`](../markers/sources/li2025_bzm.tsv).
The thresholds, the input hash and the per-domain counts are in [`li2025_bzm.json`](../markers/sources/li2025_bzm.json).

A Stereo-seq spot holds about 3–8 cells, so a domain's markers include genes of neighbouring cells.
The script therefore drops any gene that the paper's own scRNA-seq calls a marker of a non-cardiomyocyte cell type.
Most domain markers fail this filter; the JSON gives the count.

## What it changes for our marker sets

| Finding | Consequence |
|---------|-------------|
| *nppb* is also a marker of activated endocardium in their scRNA-seq | Treat *nppb* as non-specific, not as a core border-zone marker |
| *pcna*, *ccng1*, *ccng2* also mark proliferating non-cardiomyocytes | The proliferating domain is mixed; the proliferating axis needs cardiomyocyte-specific genes |
| *nppa, csrp3, xirp1, flnca, rrad, aldoab, cyc1, ACTC1, CT027980.1* pass the filter in three or more BZm domains and are not remote-myocardium markers | Candidates for a border-zone core set. Several are abundant sarcomere genes, so they must pass the depth check |
| *tnfrsf11b, hspb6, hspb11* mark the activated domain at 6 hours | The first injured time point in Hu et al. (2022) is 3 days, so this state may be absent or faded in our data |
| *myh7, cox6a2, cox5b2* rise from the proliferating to the redifferentiated domain | A redifferentiation signal exists at 7–14 dpa; our four axes have no place for it yet |

## Limits

- The injury is ventricular amputation, not cryoinjury. The authors confirm *tnfrsf11b* and some later responses in cryoinjured hearts too.
- Li et al. annotate genes with a newer genome annotation. Genes missing from the Ensembl 92 annotation of Hu et al. (2022) are dropped; the JSON gives the count.
- In this model the BZm states follow each other in one place. That predicts that our border-zone axis co-occurs with each state axis, and that the state axes peak at different times.
