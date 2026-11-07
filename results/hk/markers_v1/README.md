# hk markers v1

Cited marker sets for the four CM state axes, plus control and ambiguous sets. Every gene is in the dataset's gene list, carries an Ensembl ID from the dataset's own GRCz11/Ens92 mapping, and has a DOI source. Each claim behind a gene is one row in `markers/sources/hk_corpus_claims.tsv` (verified = yes means the passage or table row was read in the full text or supplementary file). Papers: `literature/hk-corpus.md`.

The authors' published clusters were not used to choose any gene. These are candidates (`status: candidate`); human-owned `markers/*.yaml` were not changed.

## Files

| file | axis | core | other roles |
|---|---|---|---|
| injury_responsive.yaml | injury_responsive | 7 | 2 non_specific, 1 dissociation_sensitive |
| border_zone.yaml | border_zone | 10 | 3 non_specific |
| dedifferentiated.yaml | dedifferentiated | 13 | 5 non_specific (maturity-loss, secondary), 4 library_size_suspect |
| proliferating.yaml | proliferating | 7 | 7 ortholog_mapping_ours |
| ambiguous.yaml | ambiguous | 0 | 5 non_specific (sensitivity analysis only) |
| controls.yaml | dissociation | 0 | 28 dissociation_sensitive |
| ambient.yaml | ambient | 0 | 155 non_specific (non-CM cluster markers) |

Rule used to split overlapping claims: a gene goes to the one state where its evidence is most specific. Early, activated-domain-restricted genes in spatial data go to injury. Genes found across all border-zone domains and in the sorted nppa+ cluster go to border zone. Embryonic TFs and glycolysis go to dedifferentiation. When two independent papers put the same gene in two states and nothing breaks the tie, it goes to ambiguous.yaml.

## Injury-responsive (IEG-free)

Core: tnfrsf11b, rln3a, hspb9, kcnh6a, sdc4, slc38a3b, klf1.

Why: Li 2025 (Stereo-seq, so no tissue dissociation) finds an activated border-zone myocardium domain at 6 h after injury. tnfrsf11b is its most enriched gene and is absent from uninjured hearts by FISH. The other five Li genes are activated-domain markers with low background (pct.2 0.02 to 0.17). klf1 comes from the team draft (Ogawa 2021) and is unverified.

Excluded:
- All immediate-early and heat-shock genes (fos*, jun*, egr*, atf3, ier*, dusp1, hsp70*, hspb1, socs3a, ubb, ...): every one is in the dissociation control set (van den Brink 2017; Neuschulz 2022). socs3a is up in Li's activated domain in situ, so it is kept in this file as `dissociation_sensitive` only.
- ifrd1 (immediate-early-like) and eif4a1a (translation, tracks depth): `non_specific`.
- hspb6, hspb11, rrad: also claimed for the border zone by Honkoop: moved to ambiguous.yaml.
- Wu 2016 injury-area genes (grn1/grn2/ctsd/tmsb4x/lygl1): immune/fibroblast genes, not CM states.

Caveats: Li used a resection model and the activated domain peaks at 6 h, but this dataset starts at 3 dpi. Li's activated genes also appear in its dedifferentiated domain (tnfrsf11b 1.44, rln3a 1.09), so some injury/dedifferentiation co-occurrence is built into the biology. Only 7 core genes: low-count genes will be noisy in shallow libraries.

## Border zone

Core: nppb, mustn1b, csrp3, flnca, xirp1, xirp2a, ankrd1a, crip2, vwa1, psmd14.

Why: nppb, mustn1b, csrp3, flnca and xirp1 each have two independent sources: Honkoop's nppa+ sorted cluster 7 (DE table) and Li's border-zone domains (all four or five domains). xirp2a has Honkoop only. ankrd1a has a 2026 title-level claim only. crip2, vwa1 and psmd14 come from Wu 2016 via the team draft; I could not reopen the Wu supplement (download blocked), so they are unverified.

Excluded:
- nppa: atrial identity gene and ventricular trabecular marker at baseline. A score weighted by nppa lights up atrial and whole-heart libraries and tracks depth. It is also claimed by dedifferentiation (Hu 2022) and injury (Li 2025), so it is in ambiguous.yaml.
- tagln: named a BZ CM marker (Constanty 2025, Honkoop), but it is a smooth-muscle/fibroblast gene (Li uses it for SMC domains): ambient risk, `non_specific`.
- myl6 (broad non-muscle myosin) and desma (also a remote-myocardium marker in Li): `non_specific`.
- vmhc (Honkoop BZ ISH) is not in this dataset's gene list.

Caveat: the Honkoop comparison is cluster 7 against cluster 2, and cluster 2 is defined by mitochondrial transcripts. Part of that DE table may be a quality/depth contrast. BZ genes are never zero in uninjured ventricle, so the useful contrast is the change from control ventricular CMs.

## Dedifferentiated (gain genes primary)

Core: gata4, hand2, nkx2.5, tbx20, runx1, pkma, hk1, ldha, pdk2a, pgk1, eno1a, tpi1b, myom1b.

Why: embryonic cardiac TFs re-expressed by regenerating CMs (gata4 reporter and lineage tracing, Kikuchi 2010; hand2 and nkx2.5 at the wound edge, Li 2025; tbx20 secondary; runx1 named as a dedifferentiation marker by Beisaw 2020, shown in rat CMs). Glycolysis induced in border-zone CMs: hk1, pkma, ldha, pdk2a by ISH, pgk1, eno1a, tpi1b in Honkoop's DE table (Honkoop, Fukuda 2020, Constanty 2025). myom1b is an embryonic sarcomere gene seen at 7 dpi but not 3 dpi.

Secondary (role `non_specific`, expected DOWN): cox7a1, ckma, pgam2. Papers disagree here. Honkoop's cluster 7 has these genes UP (2.5 to 2.7 log2FC), and Li puts them in the re-differentiating domain. Loss of detection also tracks depth. So maturity loss is not a reliable primary dedifferentiation readout; use it only after the thinning test.

library_size_suspect (out of core): ttn.1, ttn.2, bves, synpo2lb, the genes Hu 2022 used for the authors' small dedifferentiated cluster. In this annotation ttn.1 = ENSDARG00000000563 (ttnb) and ttn.2 = ENSDARG00000028213 (ttna). Both are also the top genes of Honkoop's cluster 7 (4.8 and 4.5 log2FC). Titin is the longest and most abundant CM transcript, so a high-ttn group is the likeliest depth artefact.

Other `non_specific`: gapdh (housekeeping), tnnc2 (team draft, unverified). aldoaa is in ambiguous.yaml: Hu 2022 uses it for the main ventricular cluster, while Honkoop and Li have it up in the border zone.

Caveat: Honkoop describes border-zone CMs as the dedifferentiating cells, so glycolysis genes are also border-zone genes by their reading. The two axes are expected to co-occur. Keeping the sets disjoint is what makes that co-occurrence measurable.

## Proliferating

Core: pcna, cdk1, plk1, mcm5, aurkb, ccnb1, mki67.
ortholog_mapping_ours: top2a, birc5a, cenpf, ube2c, nusap1, cdc20, mcm2 (Tirosh 2016 lists as redistributed in Seurat; human-to-zebrafish mapping is ours).

Why: pcna and cdk1 up in Li's proliferating border-zone domain; plk1 required for CM proliferation (Jopling 2010); mcm5 is the proliferation readout in Lekkos 2025; aurkb and ccnb1 mark division rather than endoreplication (Honkoop, Ellman 2026; both in mouse CMs).

Excluded or flagged: the Li 2025 passage that names mki67, smc2, smc4 and mad2l1 is about proliferating macrophages, not CMs. So mki67 is sourced only to Tirosh (unverified), and smc2/smc4/mad2l1 are not used. pcna also marks DNA repair. Long cell-cycle genes drop out in shallow cells. Proliferating non-CMs can form doublets with CMs. Cell-cycle mRNA is not division (polyploid CMs).

## Controls

- controls.yaml (axis dissociation): the 28 dissociation genes that appear in two or more independent dissociation sets, with at least one set measured directly in zebrafish (van den Brink 2017 mouse + zebrafish fin; Neuschulz 2022 zebrafish SLAM-seq). The injury axis must beat this control. It must not simply drop it, because van Duijvenboden 2019 shows AP-1 genes are part of the real BZ program.
- ambient.yaml (axis ambient): 155 non-CM top-50 cluster markers from Hu 2022 Supplementary Data 2 (endothelial/endocardial, fibroblast, epicardial, macrophage, neutrophil, T/B cell, perivascular). Use them to estimate ambient and doublet contamination per cell. Erythrocyte genes are left out: their only support was a name pattern, not a citation. Seurat ".1" duplicate symbols that have no Ensembl ID of their own are skipped.

## Disagreements with the team draft (`markers/*.yaml`)

- nppa: the draft has it as border-zone core. v1 moves it to ambiguous.
- hspb6, hspb11: draft injury core. v1 moves them to ambiguous, because Honkoop has both in the BZ cluster.
- myom1b, gata4: same axis as the draft. tnnc2: draft dedifferentiated core, v1 `non_specific` (unverified).
- pcna: draft proliferating core. v1 keeps it as core, with the DNA-repair caveat.
- jun/junba/junbb/fosab/fosl1a: draft injury `dissociation_sensitive`. v1 removes them from the injury file entirely and uses the dissociation axis in controls.yaml.
- Ensembl IDs follow the dataset mapping (GRCz11 Ens92 from the GEO matrix). All IDs in the team draft match it. Honkoop's source table uses older IDs for some genes (e.g. csrp3 ENSDARG00000069975 there, ENSDARG00000101706 here), so match by the dataset ID, not the paper's.
