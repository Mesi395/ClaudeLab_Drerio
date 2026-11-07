# C003 reading rules (researcher), fixed 2026-11-07 15:36 CET before computing

Hearts: 18 (Hr2a/b, H8a/v, Hr6a/v, Hr7a/v pooled per heart). Ctrl = H5, H6, H7, H8 (4);
3 dpi = Hr11, Hr12, Hr24, Hr27 (4); 7 dpi = Hr1, Hr2, Hr6, Hr7, Hr8, Hr9, Hr13 (7);
30 dpi = Hr3, Hr4, Hr19 (3).

## Q1 Positive control (F5)

- Depth control: every cell binomially thinned to 500 UMIs (seed 0, as in C002), then
  score_genes (primary) on log1p CP10k of the thinned counts. 100 expression-matched null sets
  (25 bins, as C002) scored the same way.
- Per heart: (a) median thinned score; (b) **tail fraction** = share of the heart's cells whose
  score exceeds the 99th percentile of all null-set scores over that heart's cells (expected 0.01
  under the null).
- Predicted contrasts: border-zone and proliferating, 7 dpi > Ctrl; injury-responsive, 3 dpi > Ctrl;
  dedifferentiated, 3 dpi > Ctrl and 7 dpi > Ctrl (Li 2025 places it 12 hpa–3 dpa).
- Effect: difference of heart medians (injured − Ctrl) and AUC = P(injured heart > Ctrl heart)
  over all pairs. **Pass** if AUC ≥ 0.8 and the difference > 0, on the tail fraction (primary)
  and reported for the median score. Sensitivity, declared now: the same with Hr11 excluded.
- With 4 Ctrl hearts, AUC 0.8 is a descriptive bar, not significance.

## Q2 Hr11 excess: ambient RNA vs 3 dpi biology

- Depth bins on raw total counts, fixed: [0,600), [600,900), [900,1350), [1350,inf).
- Per library group (Hr11, Hr12, Hr24, Hr27, Ctrl pooled) and bin: mean share of UMIs per gene.
  Ratio = group / Ctrl in the same bin.
- Gene panels, fixed now (only genes present in the data are used):
  - target: nppa, ttn.2 (hinted, FU-6)
  - atrial CM genes (nppa is atrial-enriched, so atrial soup is one ambient source): myh6
  - abundant CM genes: the 10 genes with the highest mean share in Ctrl cells, excluding mt-, NC-,
    ribosomal (rp[ls]/RP[LS]) and haemoglobin genes, and the target genes
  - non-CM genes: hbba1.1, hbaa1.1, hbba1 (blood); col1a1a, col1a2, postnb, fn1b (fibroblast);
    kdrl, cdh5, fli1a (endothelial); mpeg1.1, lyz (immune); tcf21 (epicardial)
  - technical: mito share including NC-002333.* (F2)
- Predictions:
  - **Ambient** (Hr11 soup): Hr11-specific; at least half of the non-CM genes raised ≥ 2× in
    Hr11's shallow bins (bins 0–1); the excess ratio falls with depth while the excess molecules
    per cell (share difference × bin median depth) stay within a factor of 2 across bins.
  - **3 dpi biology**: the nppa excess (≥ 2× in bins 0–1) is shared by at least 2 of Hr12, Hr24,
    Hr27, and non-CM genes are not broadly raised in Hr11.
  - Anything else is **inconclusive**; say which part failed.
- The test can fail both ways: if non-CM genes are not raised and the other 3 dpi hearts share the
  excess, ambient loses.
