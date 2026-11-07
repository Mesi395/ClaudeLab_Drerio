"""C007: build disjoint, CM-internal axis sets from Li 2025 Supplementary Data 2 (MOESM5.xlsx).

Rule: research/cm-state-axes/cycles/C007/reading-rules.md section 1 (fixed before building).
uv run python -I scripts/c007_build_sets.py research/cm-state-axes/sources/li2025_supp/MOESM5.xlsx \
    data/clean/hu2022_cm_allgenes_sanitized.h5ad markers/li2025_sd2_sets.tsv
"""

import json
import os
import re
import subprocess
import sys
import xml.etree.ElementTree as ET
import zipfile
from datetime import datetime
from pathlib import Path

import anndata as ad
import pandas as pd

NS = {"m": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}
RNS = "{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id"
TECH = re.compile(r"^(mt-|NC-|rp[ls]|RP[LS]|hb)")
N_GENES = 25
TIE_ORDER = ["proliferating", "dedifferentiated", "injury_responsive", "border_zone_like"]
BZM = [
    "1:BZm (activ.)",
    "2:BZm (dediff.)",
    "3:BZm (dediff.-prolif.)",
    "4:BZm (prolif.)",
    "5:BZm (re-diff.)",
]


def read_xlsx(path):
    z = zipfile.ZipFile(path)
    ss = [
        "".join(t.text or "" for t in si.iter(f"{{{NS['m']}}}t"))
        for si in ET.fromstring(z.read("xl/sharedStrings.xml")).findall("m:si", NS)
    ]
    wb = ET.fromstring(z.read("xl/workbook.xml"))
    rels = {
        r.get("Id"): r.get("Target") for r in ET.fromstring(z.read("xl/_rels/workbook.xml.rels"))
    }
    out = []
    for i, sh in enumerate(wb.find("m:sheets", NS)):
        t = rels[sh.get(RNS)].lstrip("/")
        t = t if t.startswith("xl/") else "xl/" + t
        rows = []
        for r in ET.fromstring(z.read(t)).find("m:sheetData", NS).findall("m:row", NS):
            vals = []
            for c in r.findall("m:c", NS):
                v = c.find("m:v", NS)
                x = v.text if v is not None else ""
                vals.append(ss[int(x)] if c.get("t") == "s" and x else x)
            rows.append(vals)
        out.append((i, sh.get("name").strip(), rows))
    return out


def to_df(rows):
    df = pd.DataFrame(rows[1:], columns=rows[0])
    for c in ("p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj"):
        df[c] = pd.to_numeric(df[c])
    return df


def main(xlsx, h5ad, out_tsv):
    sheets = read_xlsx(xlsx)
    # sheet 0 = scRNA-seq chamber, 1 = scRNA-seq regeneration, 2 = Stereo-seq regeneration
    names = [n for _, n, _ in sheets]
    assert names[:3] == ["scRNA-seq chamber", "scRNA-seq regeneration", "Stereo-seq regeneration"]
    cham, sreg, streg = (to_df(sheets[i][2]) for i in range(3))
    present = set(ad.read_h5ad(h5ad, backed="r").var_names)

    def is_marker(df, l2fc):
        return (df.p_val_adj < 0.05) & (df["pct.1"] >= 0.2) & (df.avg_log2FC >= l2fc)

    sig = (cham.p_val_adj < 0.05) & (cham.avg_log2FC >= 0.5)
    noncm = set(cham[sig & ~cham.cluster.isin(["Cardiomyocytes (A)", "Cardiomyocytes (V)"])].gene)
    atr = set(cham[sig & (cham.cluster == "Cardiomyocytes (A)")].gene)
    ven = set(cham[sig & (cham.cluster == "Cardiomyocytes (V)")].gene)
    drop = noncm | (atr - ven)
    if os.environ.get("C007_DROP_CM_IDENTITY") == "1":  # exploratory variant (rules addendum)
        drop = drop | atr | ven

    def keep(g):
        return g in present and not TECH.match(g) and g not in drop

    cand = {}
    for axis, df, cl, l2 in [
        ("dedifferentiated", sreg, "Cardiomyocytes (dediff.)", 0.5),
        ("proliferating", sreg, "Proliferating cells", 0.5),
        ("injury_responsive", streg, "1:BZm (activ.)", 0.25),
    ]:
        x = df[(df.cluster == cl) & is_marker(df, l2)].copy()
        x = x[x.gene.map(keep)].sort_values("avg_log2FC", ascending=False)
        x["sheet"] = "scRNA-seq regeneration" if df is sreg else "Stereo-seq regeneration"
        cand[axis] = x.drop_duplicates("gene").reset_index(drop=True)
    b = streg[streg.cluster.isin(BZM) & is_marker(streg, 0.25)]
    b = b[b.gene.map(keep)]
    g = b.groupby("gene").agg(
        n_domains=("cluster", "nunique"),
        avg_log2FC=("avg_log2FC", "mean"),
        pct1=("pct.1", "mean"),
        pct2=("pct.2", "mean"),
        p_val_adj=("p_val_adj", "max"),
        cluster=("cluster", lambda s: ";".join(sorted(s))),
    )
    g = g[g.n_domains >= 3].sort_values("avg_log2FC", ascending=False).reset_index()
    g = g.rename(columns={"pct1": "pct.1", "pct2": "pct.2"})
    g["sheet"] = "Stereo-seq regeneration"
    cand["border_zone_like"] = g
    # disjoint: best rank wins, ties by TIE_ORDER
    best = {}
    for axis in TIE_ORDER:
        for r, gene in enumerate(cand[axis].gene):
            if gene not in best or r < best[gene][0]:
                best[gene] = (r, axis)
    rows = []
    for axis in TIE_ORDER:
        x = cand[axis]
        x = x[x.gene.map(lambda q, a=axis: best[q][1] == a)].head(N_GENES).reset_index(drop=True)
        for i, r in x.iterrows():
            rows.append(
                {
                    "axis": axis,
                    "gene": r.gene,
                    "sheet": r.sheet,
                    "cluster": r.cluster,
                    "rank": i + 1,
                    "avg_log2FC": round(float(r.avg_log2FC), 4),
                    "pct.1": r["pct.1"],
                    "pct.2": r["pct.2"],
                    "p_val_adj": r.p_val_adj,
                    "n_candidates_before_disjoint": len(cand[axis]),
                }
            )
    out = pd.DataFrame(rows)
    out.to_csv(out_tsv, sep="\t", index=False)
    meta = {
        "commit": subprocess.check_output(["git", "rev-parse", "HEAD"], text=True).strip(),
        "dirty_tree": True,
        "command": "uv run python -I " + " ".join(sys.argv),
        "source": "Li et al. 2025 Nat Commun, doi:10.1038/s41467-025-59070-0, Supplementary "
        "Data 2 (MOESM5.xlsx, Springer static, downloaded 2026-11-07)",
        "n_genes_per_axis": out.groupby("axis").size().to_dict(),
        "n_dropped_noncm_or_atrial": len(drop),
        "rule": "research/cm-state-axes/cycles/C007/reading-rules.md section 1",
        "time": datetime.now().astimezone().isoformat(timespec="seconds"),
    }
    Path(out_tsv + ".json").write_text(json.dumps(meta, indent=2))
    print(json.dumps(meta, indent=1))
    for axis in TIE_ORDER:
        print(axis, list(out[out.axis == axis].gene))


if __name__ == "__main__":
    main(*sys.argv[1:4])
