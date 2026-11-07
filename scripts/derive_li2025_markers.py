"""Derive border-zone cardiomyocyte markers from Li et al. (2025), Supplementary Data 2.

Li, L., Lu, M., Guo, L., et al. (2025). An organ-wide spatiotemporal transcriptomic and cellular
atlas of the regenerating zebrafish heart. Nature Communications, 16, 3716.
https://doi.org/10.1038/s41467-025-59070-0

Stereo-seq spots hold about 3-8 cells, so a domain's markers include genes of neighbouring
non-cardiomyocytes. A gene is kept only if the paper's own scRNA-seq does not call it a marker
of any non-cardiomyocyte cell type. Every threshold is a constant below and is written to the
sidecar JSON with the input hash.

Run:
    uv run --with openpyxl python scripts/derive_li2025_markers.py
"""

import hashlib
import json
import subprocess
import urllib.request
from pathlib import Path

import h5py
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
URL = (
    "https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-025-59070-0/"
    "MediaObjects/41467_2025_59070_MOESM5_ESM.xlsx"
)
XLSX = ROOT / "data/external/li2025_supplementary_data_2.xlsx"
HU_FEATURES = ROOT / "data/raw/Hr11_filtered_matrix.h5"
OUT = ROOT / "markers/sources/li2025_bzm.tsv"

# Spatial domain markers (Seurat FindAllMarkers output in the sheet).
MAX_P_ADJ = 0.05
MIN_LOG2FC = 0.5
MIN_PCT_IN = 0.25
MIN_PCT_DIFF = 0.10
# A gene is a non-cardiomyocyte marker if it reaches this in any non-CM scRNA-seq cluster.
NON_CM_LOG2FC = 1.0
CM_CLUSTERS = {"Cardiomyocytes (A)", "Cardiomyocytes (V)", "Cardiomyocytes (dediff.)"}
EXCLUDED_SC_CLUSTERS = {"Others", "Red blood cells", "Proliferating cells"}

# Li et al.'s border-zone myocardium domains, in temporal order, and the axis each informs.
BZM = {
    "1:BZm (activ.)": "injury_responsive",
    "2:BZm (dediff.)": "dedifferentiated",
    "3:BZm (dediff.-prolif.)": "dedifferentiated",
    "4:BZm (prolif.)": "proliferating",
    "5:BZm (re-diff.)": "redifferentiated",
}
REMOTE = {"6:RZm1", "7:RZm2", "8:Vm (comp.)"}


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def hu_gene_ids():
    with h5py.File(HU_FEATURES) as f:
        names = [n.decode() for n in f["matrix/features/name"][:]]
        ids = [i.decode() for i in f["matrix/features/id"][:]]
    return dict(zip(names, ids, strict=True))


def main():
    if not XLSX.exists():
        XLSX.parent.mkdir(parents=True, exist_ok=True)
        urllib.request.urlretrieve(URL, XLSX)

    spatial = pd.read_excel(XLSX, sheet_name="Stereo-seq regeneration")
    sc = pd.read_excel(XLSX, sheet_name="scRNA-seq regeneration")

    passes = (
        (spatial.p_val_adj < MAX_P_ADJ)
        & (spatial.avg_log2FC >= MIN_LOG2FC)
        & (spatial["pct.1"] >= MIN_PCT_IN)
        & (spatial["pct.1"] - spatial["pct.2"] >= MIN_PCT_DIFF)
    )
    spatial = spatial[passes]

    sc_sig = sc[(sc.p_val_adj < MAX_P_ADJ) & (sc.avg_log2FC >= NON_CM_LOG2FC)]
    non_cm = sc_sig[~sc_sig.cluster.isin(CM_CLUSTERS | EXCLUDED_SC_CLUSTERS)]
    non_cm_of = non_cm.groupby("gene").cluster.apply(lambda c: ";".join(sorted(c)))
    cm_of = (
        sc_sig[sc_sig.cluster.isin(CM_CLUSTERS)]
        .groupby("gene")
        .cluster.apply(lambda c: ";".join(sorted(c)))
    )
    remote_genes = set(spatial[spatial.cluster.isin(REMOTE)].gene)

    bzm = spatial[spatial.cluster.isin(BZM)].copy()
    bzm["axis_hint"] = bzm.cluster.map(BZM)
    bzm["n_bzm_domains"] = bzm.groupby("gene").cluster.transform("nunique")
    bzm["also_remote_marker"] = bzm.gene.isin(remote_genes)
    bzm["non_cm_marker_of"] = bzm.gene.map(non_cm_of).fillna("")
    bzm["cm_scrna_marker_of"] = bzm.gene.map(cm_of).fillna("")
    ids = hu_gene_ids()
    bzm["ensembl_hu2022"] = bzm.gene.map(ids).fillna("")
    bzm["keep"] = (bzm.non_cm_marker_of == "") & (bzm.ensembl_hu2022 != "")

    cols = [
        "gene",
        "ensembl_hu2022",
        "cluster",
        "axis_hint",
        "avg_log2FC",
        "pct.1",
        "pct.2",
        "p_val_adj",
        "n_bzm_domains",
        "also_remote_marker",
        "cm_scrna_marker_of",
        "non_cm_marker_of",
        "keep",
    ]
    bzm = bzm[cols].sort_values(["cluster", "avg_log2FC"], ascending=[True, False])
    OUT.parent.mkdir(parents=True, exist_ok=True)
    bzm.to_csv(OUT, sep="\t", index=False, float_format="%.4g")

    commit = subprocess.run(
        ["git", "rev-parse", "HEAD"], cwd=ROOT, capture_output=True, text=True, check=False
    ).stdout.strip()
    summary = {
        "source": "https://doi.org/10.1038/s41467-025-59070-0",
        "input": {"url": URL, "sha256": sha256(XLSX)},
        "gene_ids_from": {
            "path": str(HU_FEATURES.relative_to(ROOT)),
            "sha256": sha256(HU_FEATURES),
        },
        "commit": commit,
        "command": "uv run --with openpyxl python scripts/derive_li2025_markers.py",
        "thresholds": {
            "max_p_adj": MAX_P_ADJ,
            "min_log2fc": MIN_LOG2FC,
            "min_pct_in": MIN_PCT_IN,
            "min_pct_diff": MIN_PCT_DIFF,
            "non_cm_log2fc": NON_CM_LOG2FC,
        },
        "rows": len(bzm),
        "kept_per_domain": bzm[bzm.keep].groupby("cluster").size().to_dict(),
        "dropped_non_cm": int((bzm.non_cm_marker_of != "").sum()),
        "dropped_not_in_hu2022": int(((bzm.non_cm_marker_of == "") & ~bzm.keep).sum()),
    }
    OUT.with_suffix(".json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
