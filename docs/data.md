# The data

The problem dataset is the cardiomyocytes from 22 public libraries in Hu et al. (2022).
`scripts/get_data.sh` downloads those libraries from GEO.

## Source

Hu, B., Lelek, S., Spanjaard, B., et al., & Junker, J. P. (2022).
Origin and function of activated fibroblast states during zebrafish heart regeneration.
*Nature Genetics, 54*(8), 1227–1237. https://doi.org/10.1038/s41588-022-01129-5
(open access: [PMC7613248](https://pmc.ncbi.nlm.nih.gov/articles/PMC7613248/)).

| Resource | Where |
|----------|-------|
| Counts, 20 libraries | GEO [GSE159032](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE159032), one `*_filtered_matrix.h5` per library |
| Counts, Hr8 and Hr11 | GEO [GSE158919](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE158919), from the Lelek et al. morphine study |

The slide's library names map to GEO labels unchanged, except Hr2a = GEO `Hr2` and Hr2b = GEO `Hr2_morecells`.
The GEO files hold all cells of each library; the problem dataset is the cardiomyocyte subset.

## Technical traps already visible

| Trap | Evidence | Where checked |
|------|----------|---------------|
| Hr11 has low depth | Median 713 UMIs per barcode in Hr11; 1,069–7,247 in the other 21 libraries† | GEO filtered matrices, all barcodes |
| Chemistry differs | Hr24 and Hr27 (3 dpi) are 10x v3; the other 20 are v2 | GEO sample records |
| Two studies are mixed | Hr8 and Hr11 come from the morphine study (GSE158919) | GEO |
| Chamber-only libraries | H8a, Hr6a, Hr7a are atrium only; H8v, Hr6v, Hr7v are ventricle only | GEO sample titles |
| No ambient RNA correction | The paper reports none | paper, Methods |
| Cardiomyocytes have fewer transcripts | The authors report about 1,600 transcripts per cardiomyocyte against about 4,500 for other cells | paper, Results |

† Computed by us: the median of total counts per barcode in each GEO `*_filtered_matrix.h5`, over all cell types.

The low depth of Hr11, a 3 dpi library, means any score that tracks library size will also look like an injury response.

## Spatial reference data

| Dataset | Accession | Use |
|---------|-----------|-----|
| Wu et al. (2016), tomo-seq of 3 and 7 dpi hearts | [GSE74652](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE74652) | Border-zone genes from position along the injury axis |
| Honkoop et al. (2019), sorted *nppa*:mCitrine+ border-zone cardiomyocytes | [GSE139218](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE139218) | Direct border-zone reference |
