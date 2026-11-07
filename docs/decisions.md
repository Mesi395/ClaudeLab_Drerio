# Analysis decisions

One line each: decision · reason · source.

- No morphine covariate · constant ("Ctrl") for all 32,044 cells, and one heart per timepoint would make it the same as regressing on heart · GEO GSE158919 sample records (Hr8, Hr11 "morphine treatment: no"); data check of the `morphine` column; Lelek et al. 2020, doi:10.1101/2020.10.01.322560
- Sensitivity run without Hr8 and Hr11 · they come from a separate experiment; Hr11 is the shallowest library with the most ambient RNA · GEO GSE158919; `results/hk/library_flags.tsv`
- Never regress on library, heart or timepoint; chemistry not a covariate · chemistry exists only at 3 dpi, so it is inseparable from biology · `config/dataset.yaml`
- Thin to 600 UMIs before the score_genes-type score · cell-cycle and ttn.2 signal track depth within library · senior review 15:39 (observation tier); `results/hk/protocol_v1_summary.csv` kappa columns
- Decoys matched on detection rate within each library (500 sets, 20 bins) · makes the null follow each library's depth; scDRS matches on mean and variance by default · Zhang et al. 2022 scDRS, doi:10.1038/s41588-022-01167-z
- Pooled-z p-values instead of (1+#)/(1+B) · the plain empirical p is floored at 1/501, so BH within a heart of thousands of cells called 0 positives on every axis in run 1 · scDRS (same DOI); JOURNAL.md 2026-11-07 16:07
- Non-CM veto uses the 4 fibroblast control genes of `config/axes.yaml`, not the 106 soup genes · the soup-gene version vetoed 59% of cells in run 1 · `config/axes.yaml` controls; JOURNAL.md 2026-11-07 16:07
- Use continuous `ambient_frac` / hemoglobin fraction, not the ambient flag · at ~700 UMIs one count reads as 10–30% contamination · `src/zf/confounds.py` docstring
- Chamber: Card A's myosin call (myh7, myh7l vs myh6) for whole-heart cells; uncertain chamber cells are left out of both strata · the multi-gene thinned chamber score was not built in time · `results/hk/va_accuracy.tsv`
- Partial conjunction by Simes, REML meta-analysis, dip test and taxonomy odds ratios: not run in v1 · time · JOURNAL.md 2026-11-07 16:00
- Dissociation / heat-shock control axis scored with the same protocol (Card A's cited dissociation set plus hsp70/hspa*, hsp90aa1*, dnajb1*, fos*, jun*, egr1, atf3 present in the data, minus every scored axis gene) · cells were dissociated enzymatically at 37 °C, which induces immediate-early and heat-shock genes that vary between samples and overlap the injury response · Neuschulz et al. 2022 MSB, doi:10.15252/msb.202211147; van den Brink et al. 2017 Nat Methods, doi:10.1038/nmeth.4437; domain input 2026-11-07 15:09
- Injury axis judged against dissociation, never regressed on it · IEGs are also part of the real border-zone program (AP-1), so removing the dissociation score would remove biology; instead report the within-heart partial ρ (given log UMI) and the injury verdict with and without the top 10% dissociation cells per heart · `results/hk/injury_vs_dissociation.csv`
- Ventricle-minus-atrium contrast per heart for injury, border zone and dissociation · dissociation hits both chambers alike, a local injury response is ventricle-only · owner 15:52 + domain input 15:09; `results/hk/chamber_contrast.csv`
- Marker bootstrap with 5 resamples in this run, not 200 · time (16:25 cut); the metric is reported with its resample count in `results/hk/protocol_v1_summary.csv` · JOURNAL.md 2026-11-07 16:12
