"""Run labeling protocol v1 (src/zf/protocol.py) on the label-free cardiomyocytes.

Per stratum (all cells, ventricle, atrium) and axis: per-cell calls (kept on the compute host in
data/derived/), the per-heart gate, the axis verdict and three stability metrics: thinning
agreement (Cohen's kappa, thinned vs unthinned calls), marker bootstrap (80% of genes) and
leave-one-heart-out (the `robust` column of the gate). Sensitivity: gate without the two hearts
of the morphine study (Hr8, Hr11).

Usage: uv run python scripts/run_protocol.py [--boot N]
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
from sklearn.metrics import cohen_kappa_score

from zf import io, protocol, provenance

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "hk"
MORPHINE_STUDY = ["Hr8", "Hr11"]  # GEO GSE158919; morphine treatment: no


def axis_genes() -> dict:
    axes = yaml.safe_load((ROOT / "config" / "axes.yaml").read_text())["axes"]
    out = {}
    for ax in axes:
        m = yaml.safe_load((ROOT / ax["markers"]).read_text())
        out[ax["name"]] = [g["symbol"] for g in m["genes"] if g["role"] in ax["score_roles"]]
    return out


DISSOC_EXTRA = r"^(hsp70|hspa1[ab]?$|hsp90aa1|dnajb1[ab]?$|fos|jun|egr1$|atf3$)"  # heat-shock only


def dissoc_genes(var_names, exclude) -> list:
    """Dissociation / heat-shock control: Card A's cited set plus the named IEG and HSP
    families (domain input 15:09), disjoint from every scored axis."""
    m = yaml.safe_load((OUT / "markers_v1" / "controls.yaml").read_text())
    cited = {g["symbol"] for g in m["genes"]}
    extra = set(pd.Index(var_names)[pd.Index(var_names).str.match(DISSOC_EXTRA)])
    return sorted((cited | extra).intersection(var_names) - set(exclude))


def chamber_contrast(out, obs, cham, axis) -> pd.DataFrame:
    """Per heart with both chambers: ventricle minus atrium mean thinned score and positive
    fraction. Dissociation hits both chambers alike; a local injury response is ventricle-only."""
    df = pd.DataFrame(
        {
            "heart": obs["heart"].astype(str).to_numpy(),
            "time": obs["condition"].astype(str).to_numpy(),
            "chamber": cham,
            "score": out["score_mean_thinned"].to_numpy(),
            "pos": out["call"].eq("positive").to_numpy(),
        }
    )
    df = df[df["chamber"] != "uncertain"]
    t = df.pivot_table(index=["heart", "time"], columns="chamber", values=["score", "pos"])
    t = t.dropna()
    return pd.DataFrame(
        {
            "axis": axis,
            "score_v_minus_a": t[("score", "ventricle")] - t[("score", "atrium")],
            "pos_v_minus_a": t[("pos", "ventricle")] - t[("pos", "atrium")],
        }
    ).reset_index()


def strata(a) -> dict:
    """Ventricle = GEO ventricle libraries + whole-heart cells the myosin call puts there."""
    ch = a.obs["chamber"].astype(str)
    va = a.obs["va_call"].astype(str)
    whole = ch.eq("unknown")
    return {
        "all": np.ones(a.n_obs, dtype=bool),
        "ventricle": (ch.eq("ventricle") | (whole & va.eq("ventricle"))).to_numpy(),
        "atrium": (ch.eq("atrium") | (whole & va.eq("atrium"))).to_numpy(),
    }


def regate(g: pd.DataFrame, drop) -> dict:
    g = g[~g["heart"].isin(drop)]
    return {t: protocol.gate(d["pass"]) for t, d in g.groupby("time")}


def main(n_boot: int) -> None:
    a = io.load()
    # Non-CM genes for the veto: the fibroblast control of config/axes.yaml. The 106 soup genes
    # of zf.confounds are detected >= 2 in most shallow cells (59% vetoed in the first run).
    controls = yaml.safe_load((ROOT / "config" / "axes.yaml").read_text())["controls"]
    fib = next(c for c in controls if c["name"] == "ambient_fibroblast")
    ambient = [g["symbol"] for g in fib["genes"]]
    genes = axis_genes()
    scored = {g for gs in genes.values() for g in gs}
    genes["dissociation_control"] = dissoc_genes(a.var_names, scored)
    print("dissociation_control genes:", genes["dissociation_control"], flush=True)
    st = strata(a)
    cham = np.where(st["ventricle"], "ventricle", np.where(st["atrium"], "atrium", "uncertain"))
    contrasts = []
    all_out = {}
    rng = np.random.default_rng(protocol.SEED)
    summary, gates, calls = [], [], []
    for name, mask in st.items():
        sub = a[mask].copy()
        pre = protocol.prepare(sub)
        pre_full = {
            **pre,
            "thinned_log": protocol.lognorm(pre["counts"]),
            "low": pre["low"] * False,
        }
        for axis, gs in genes.items():
            kw = dict(ambient_genes=ambient)
            out = protocol.label_axis(sub, gs, precomputed=pre, **kw)
            ref = protocol.label_axis(sub, gs, precomputed=pre_full, **kw)
            pos = out["call"].eq("positive").to_numpy()
            agree = np.zeros(len(out))
            present = [g for g in gs if g in sub.var_names]
            k = max(2, int(round(0.8 * len(present))))
            for _ in range(n_boot):
                bs = list(rng.choice(present, size=k, replace=False))
                b = protocol.label_axis(sub, bs, precomputed=pre, **kw)
                agree += b["call"].eq("positive").to_numpy()
            agree /= max(n_boot, 1)
            g = pd.DataFrame(out.attrs["gate"]).assign(stratum=name, axis=axis)
            gates.append(g)
            sens = regate(g, MORPHINE_STUDY)
            full = g.groupby("time")["supported"].first().to_dict()
            summary.append(
                {
                    "stratum": name,
                    "axis": axis,
                    "n_genes": len(present),
                    "n_cells": len(out),
                    "verdict": out.attrs["verdict"],
                    "rho_log_umi": out.attrs["rho"]["log_umi"],
                    "rho_hb_frac": out.attrs["rho"]["hb_frac"],
                    "positive": int(pos.sum()),
                    "high_confidence": int(out["high_confidence"].sum()),
                    "uncertain": int(out["call"].eq("uncertain").sum()),
                    "low_depth": int(out["low_depth"].sum()),
                    "vetoed": int(out["vetoed"].sum()),
                    "kappa_thinned_vs_full": cohen_kappa_score(out["call"], ref["call"]),
                    "boot_n": n_boot,
                    "boot_positive_stable_80": float((agree[pos] >= 0.8).mean())
                    if pos.any()
                    else np.nan,
                    "decoy_fp_rate": float(g["decoy_rate"].mean()),
                    "supported_times": ";".join(t for t, s in full.items() if s) or "none",
                    "robust_loho_times": ";".join(
                        t for t, s in g.groupby("time")["robust"].first().items() if s
                    )
                    or "none",
                    "supported_without_morphine_study": ";".join(t for t, s in sens.items() if s)
                    or "none",
                }
            )
            if name == "all":
                all_out[axis] = out
                contrasts.append(chamber_contrast(out, sub.obs, cham, axis))
                calls.append(
                    out[
                        [
                            "call",
                            "score_ucell_raw",
                            "score_mean_thinned",
                            "p_ucell",
                            "p_score",
                            "high_confidence",
                            "low_depth",
                            "vetoed",
                        ]
                    ]
                    .assign(boot_agree=agree)
                    .add_prefix(f"{axis}.")
                )
            print(summary[-1], flush=True)
    # Dissociation vs injury (domain input 15:09).
    d = all_out["dissociation_control"]
    inj = all_out["injury_responsive"]
    heart = a.obs["heart"].astype(str).to_numpy()
    log_umi = np.log10(np.asarray(a.layers["counts"].sum(1)).ravel())
    rho = protocol._partial_spearman(
        inj["score_mean_thinned"].to_numpy(), d["score_mean_thinned"].to_numpy(), log_umi, heart
    )
    dhigh = (
        pd.Series(d["score_mean_thinned"].to_numpy())
        .groupby(heart)
        .transform(lambda v: v >= v.quantile(0.9))
        .to_numpy()
    )
    keep = a[~dhigh].copy()
    inj_wo = protocol.label_axis(keep, genes["injury_responsive"], ambient_genes=ambient)
    gate_wo = pd.DataFrame(inj_wo.attrs["gate"]).groupby("time")["supported"].first().to_dict()
    gate_w = pd.DataFrame(inj.attrs["gate"]).groupby("time")["supported"].first().to_dict()
    lvl = pd.DataFrame(
        {
            "heart": heart,
            "time": a.obs["condition"].astype(str).to_numpy(),
            "dissoc_score": d["score_mean_thinned"].to_numpy(),
            "dissoc_pos": d["call"].eq("positive").to_numpy(),
        }
    )
    dl = lvl.groupby(["heart", "time"]).agg(
        n_cells=("dissoc_score", "size"),
        dissoc_score_median=("dissoc_score", "median"),
        dissoc_pos_frac=("dissoc_pos", "mean"),
    )
    dl.reset_index().to_csv(OUT / "dissociation_per_heart.csv", index=False)
    pd.DataFrame(
        [
            {
                "partial_rho_injury_dissoc_within_heart_given_logumi": rho,
                "injury_verdict_all_cells": inj.attrs["verdict"],
                "injury_verdict_without_dissoc_top10": inj_wo.attrs["verdict"],
                "injury_supported_times_all": ";".join(t for t, v in gate_w.items() if v) or "none",
                "injury_supported_times_without_dissoc_top10": ";".join(
                    t for t, v in gate_wo.items() if v
                )
                or "none",
                "injury_positive_all": int(inj["call"].eq("positive").sum()),
                "injury_positive_without_dissoc_top10": int(inj_wo["call"].eq("positive").sum()),
                "cells_removed_dissoc_top10": int(dhigh.sum()),
            }
        ]
    ).to_csv(OUT / "injury_vs_dissociation.csv", index=False)
    pd.concat(contrasts).to_csv(OUT / "chamber_contrast.csv", index=False)
    print(dl.to_string(), flush=True)
    OUT.mkdir(parents=True, exist_ok=True)
    s = pd.DataFrame(summary)
    s.to_csv(OUT / "protocol_v1_summary.csv", index=False)
    pd.concat(gates).to_csv(OUT / "protocol_v1_gate.csv", index=False)
    # Per-cell tables carry barcodes: they stay on the compute host, never in git.
    pd.concat(calls, axis=1).to_csv(io.DATA / "derived" / "labels_v1.tsv.gz", sep="\t")
    outcomes = {
        r.axis: {
            "depth_eligibility": "flag" if r.verdict == "depth_confounded" else "pass",
            "decoy_gate": "pass" if r.verdict == "supported" else "flag",
            "loho": "pass" if r.robust_loho_times != "none" else "unstable",
        }
        for r in s[s["stratum"] == "all"].itertuples()
    }
    provenance.write_manifest(
        ROOT / "results" / "runs" / "hk-protocol-v1",
        command=" ".join(["uv run python scripts/run_protocol.py", *sys.argv[1:]]),
        inputs=[io.LABEL_FREE, io.FLAGS],
        params={
            "T": protocol.T,
            "B": protocol.B,
            "bins": protocol.N_BINS,
            "q_pos": protocol.Q_POS,
            "p_neg": protocol.P_NEG,
            "rho_max": protocol.RHO_MAX,
            "boot": n_boot,
            "mito_max": protocol.MITO_MAX,
            "seed": protocol.SEED,
        },
        outcomes=outcomes,
    )
    print(s.to_string())


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--boot", type=int, default=20)
    main(ap.parse_args().boot)
