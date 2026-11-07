"""C008: re-run the C002/C003/C005 checks within ventricle-called cells; test whether the three
nppa-led axes are one program (with and without nppa/nppb).

Rules: research/cm-state-axes/cycles/C008/reading-rules.md (fixed 16:38 CET, before computing).
uv run python scripts/c008_probe.py data/clean/hu2022_cm_allgenes_sanitized.h5ad proposals/loop \
    results
"""

import json
import re
import subprocess
import sys
from collections import Counter
from datetime import datetime
from itertools import combinations
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import yaml
from scipy.stats import spearmanr

from zf.io import ALLOWED_OBS
from zf.score import matched_random_sets, rank_matrix, ucell_from_ranks

sys.path.insert(0, str(Path(__file__).parent))
from c004_probe import heart  # noqa: E402

TARGET, SEED, N_NULL, MIN_CELLS, CIRCULAR = 500, 0, 100, 50, {"ttn.2"}
NPPA_AXES = ["border_zone_like", "injury_responsive", "dedifferentiated"]
NON_CM = [
    "hbba1.1",
    "hbba1",
    "col1a1a",
    "col1a2",
    "postnb",
    "fn1b",
    "kdrl",
    "cdh5",
    "fli1a",
    "mpeg1.1",
    "lyz",
    "tcf21",
]
BINS = [0, 600, 900, 1350, np.inf]


def rho(x, y):
    return float(spearmanr(x, y)[0])


def sg(b, genes, seed):
    sc.tl.score_genes(b, list(genes), random_state=seed, score_name="_s", use_raw=False)
    return b.obs.pop("_s").to_numpy()


def coherence(M):
    keep = M.std(0) > 0
    if keep.sum() < 2:
        return np.nan
    c = np.atleast_2d(spearmanr(M[:, keep])[0])
    return float(np.nanmean(c[np.triu_indices_from(c, 1)]))


def norm(X, names, obs_names):
    b = ad.AnnData(X=X.copy(), var=pd.DataFrame(index=names), obs=pd.DataFrame(index=obs_names))
    sc.pp.normalize_total(b, target_sum=1e4)
    sc.pp.log1p(b)
    return b


def main(src, marker_dir, out_dir):
    out = Path(out_dir)
    a = ad.read_h5ad(src)
    assert set(a.obs.columns) == set(ALLOWED_OBS)
    lab = pd.read_csv(out / "c004_agent_labels.csv")
    assert (lab.cell.to_numpy() == a.obs_names.to_numpy()).all()
    vm = lab.agent_chamber.to_numpy() == "ventricle"
    counts = a.layers["counts"].tocsr().astype(np.float32)
    tot = np.asarray(counts.sum(1)).ravel()
    vn = list(a.var_names)
    vidx = {g: i for i, g in enumerate(vn)}
    rng = np.random.default_rng(SEED)
    p = np.minimum(1.0, TARGET / tot)
    th = counts.copy()
    th.data = rng.binomial(th.data.astype(np.int64), np.repeat(p, np.diff(th.indptr))).astype(
        np.float32
    )
    th.eliminate_zeros()
    means = np.asarray(norm(counts, vn, a.obs_names).X.mean(0)).ravel()
    cv = counts[vm]
    bt = norm(th[vm], vn, a.obs_names[vm])
    bu = norm(cv, vn, a.obs_names[vm])
    R = rank_matrix(cv)
    lib = a.obs["orig.ident"].astype(str).to_numpy()[vm]
    t_v = tot[vm]
    logtot = np.log10(t_v)
    ngenes = np.asarray((cv > 0).sum(1)).ravel()
    nppa = cv[:, vidx["nppa"]].toarray().ravel() / t_v
    mito_i = [i for i, g in enumerate(vn) if re.match(r"^(mt-|NC-002333)", g)]
    mito = np.asarray(cv[:, mito_i].sum(1)).ravel() / t_v
    libs = [x for x, n in Counter(lib).items() if n >= MIN_CELLS]
    libs.sort()

    axes = {}
    for f in sorted(Path(marker_dir).glob("*.yaml")):
        spec = yaml.safe_load(f.read_text())
        axes[spec["state"]] = [
            m["symbol"]
            for m in spec["genes"]
            if m["symbol"] in vidx and m["symbol"] not in CIRCULAR
        ]
    excl = np.array(sorted({vidx[g] for gs in axes.values() for g in gs}))
    Xt = bt.X.tocsc()

    # ---- 1 & 2: ties and matched null ----
    tie_rows, null_rows, scores_t, nulls_t = [], [], {}, {}
    for k, (ax, genes) in enumerate(axes.items()):
        gi = np.array([vidx[g] for g in genes])
        s_th, s_un = sg(bt, genes, SEED), sg(bu, genes, SEED)
        s_uc = ucell_from_ranks(R, gi, 1500)
        scores_t[ax] = s_th
        sets = matched_random_sets(means, gi, N_NULL, 25, SEED + k, exclude=excl)
        nul = np.stack([sg(bt, [vn[i] for i in s], SEED + j) for j, s in enumerate(sets)], 1)
        nulls_t[ax] = nul
        for lb in libs:
            m = lib == lb
            r = {"axis": ax, "lib": lb, "heart": heart(lb), "n": int(m.sum())}
            for name, s in (("sg_unthinned", s_un), ("sg_thinned", s_th), ("ucell", s_uc)):
                r[f"{name}_rho_depth"] = rho(s[m], logtot[m])
                r[f"{name}_rho_ngenes"] = rho(s[m], ngenes[m])
                r[f"{name}_rho_nppa"] = rho(s[m], nppa[m])
                r[f"{name}_rho_mito"] = rho(s[m], mito[m])
            tie_rows.append(r)
            p90, np90 = np.percentile(s_th[m], 90), np.percentile(nul[m], 90, axis=0)
            G = Xt[:, gi][m].toarray()
            coh = coherence(G)
            ncoh = [coherence(Xt[:, s][m].toarray()) for s in sets]
            null_rows.append(
                {
                    "axis": ax,
                    "lib": lb,
                    "heart": heart(lb),
                    "p90": p90,
                    "p90_null95": np.percentile(np90, 95),
                    "p90_above": bool(p90 > np.percentile(np90, 95)),
                    "coherence": coh,
                    "coherence_null95": np.nanpercentile(ncoh, 95),
                    "coherence_above": bool(coh > np.nanpercentile(ncoh, 95)),
                }
            )
        print(ax, "ties+null done", flush=True)
    ties, nulls = pd.DataFrame(tie_rows), pd.DataFrame(null_rows)
    ties.to_csv(out / "c008_ties_by_library.csv", index=False)
    nulls.to_csv(out / "c008_null_by_library.csv", index=False)

    # ---- 3: one program? ----
    cnt = Counter(g for ax in NPPA_AXES for g in axes[ax] if g not in ("nppa", "nppb"))
    axes_b = {
        ax: [g for g in axes[ax] if g not in ("nppa", "nppb") and cnt[g] == 1] for ax in NPPA_AXES
    }
    sc_b, nul_b = {}, {}
    for k, ax in enumerate(NPPA_AXES):
        gi = np.array([vidx[g] for g in axes_b[ax]])
        sc_b[ax] = sg(bt, axes_b[ax], SEED)
        sets = matched_random_sets(means, gi, N_NULL, 25, SEED + 50 + k, exclude=excl)
        nul_b[ax] = np.stack([sg(bt, [vn[i] for i in s], SEED + j) for j, s in enumerate(sets)], 1)
    prog = []
    for variant, S, N in (("a_as_is", scores_t, nulls_t), ("b_no_nppa_nppb_disjoint", sc_b, nul_b)):
        for x1, x2 in combinations(NPPA_AXES, 2):
            for lb in libs:
                m = lib == lb
                rr = rho(S[x1][m], S[x2][m])
                nr = [rho(N[x1][m, j], N[x2][m, j]) for j in range(N_NULL)]
                prog.append(
                    {
                        "variant": variant,
                        "pair": f"{x1}|{x2}",
                        "lib": lb,
                        "heart": heart(lb),
                        "rho": rr,
                        "null_p95": np.nanpercentile(nr, 95),
                        "above": bool(rr > np.nanpercentile(nr, 95)),
                    }
                )
    pg = pd.DataFrame(prog)
    pg.to_csv(out / "c008_program_by_library.csv", index=False)
    b_nppa = {
        ax: float(np.median([rho(sc_b[ax][lib == lb], nppa[lib == lb]) for lb in libs]))
        for ax in NPPA_AXES
    }

    # ---- 4: ambient, ventricle cells ----
    tv = a.obs["time"].astype(str).to_numpy()[vm]
    grp = np.where(
        np.isin(lib, ["Hr11", "Hr12", "Hr24", "Hr27"]), lib, np.where(tv == "Ctrl", "Ctrl", "")
    )
    ctrl = grp == "Ctrl"
    tech = re.compile(r"^(mt-|NC-|rp[ls]|RP[LS]|hb)")
    share_c = np.asarray(cv[ctrl].sum(0)).ravel() / t_v[ctrl].sum()
    abundant = [
        vn[i]
        for i in np.argsort(-share_c)
        if not tech.match(vn[i]) and vn[i] not in ("nppa", "ttn.2", "myh6")
    ][:10]
    panels = {
        "target": ["nppa", "ttn.2"],
        "atrial": ["myh6"],
        "abundant_cm": abundant,
        "non_cm": [g for g in NON_CM if g in vidx],
    }
    genes = [g for gs in panels.values() for g in gs]
    Sh = cv[:, [vidx[g] for g in genes]].toarray() / t_v[:, None]
    bin_of = np.digitize(t_v, BINS[1:-1])
    amb = []
    for g_ in ["Hr11", "Hr12", "Hr24", "Hr27"]:
        for bb in range(4):
            gm, cm = (grp == g_) & (bin_of == bb), ctrl & (bin_of == bb)
            if gm.sum() < 20 or cm.sum() < 20:
                continue
            for pan, gs in panels.items():
                for g in gs:
                    j = genes.index(g)
                    x, y = Sh[gm, j].mean(), Sh[cm, j].mean()
                    amb.append(
                        {
                            "group": g_,
                            "depth_bin": bb,
                            "n_group": int(gm.sum()),
                            "n_ctrl": int(cm.sum()),
                            "panel": pan,
                            "gene": g,
                            "ratio": x / y if y > 0 else np.inf,
                            "excess_molecules": (x - y) * float(np.median(t_v[gm])),
                        }
                    )
    am = pd.DataFrame(amb)
    am.to_csv(out / "c008_ambient_ventricle.csv", index=False)
    sh = am[am.depth_bin.isin([0, 1])]
    ambv = {}
    for g_ in ["Hr11", "Hr12", "Hr24", "Hr27"]:
        x = sh[sh.group == g_]
        if not len(x):
            continue
        ambv[g_] = {
            "nppa_ratio_bins01": x[x.gene == "nppa"].ratio.round(3).tolist(),
            "ttn.2_ratio_bins01_hinted": x[x.gene == "ttn.2"].ratio.round(3).tolist(),
            "myh6_ratio_bins01": x[x.gene == "myh6"].ratio.round(3).tolist(),
            "non_cm_frac_ge2x": float(
                x[x.panel == "non_cm"].groupby("gene").ratio.min().ge(2).mean()
            ),
        }

    # ---- 5: doublets, ventricle cells ----
    d = pd.read_csv(out / "c005_doublets_cells.csv").set_index("cell")
    dv = d.reindex(a.obs_names[vm]).doublet_score.to_numpy()
    dbl = {
        ax: float(
            np.median(
                [
                    abs(
                        rho(
                            scores_t[ax][(lib == lb) & ~np.isnan(dv)],
                            dv[(lib == lb) & ~np.isnan(dv)],
                        )
                    )
                    for lb in libs
                ]
            )
        )
        for ax in axes
    }

    # ---- summary ----
    def hearts_all(df, col):
        return int(df.groupby("heart")[col].all().sum()), int(df.heart.nunique())

    summ = {
        "n_ventricle_cells": int(vm.sum()),
        "libraries": libs,
        "n_libraries": len(libs),
        "n_hearts": int(len(set(heart(x) for x in libs))),
        "ties_median": {},
        "null": {},
        "program": {},
        "variant_b_sets": axes_b,
        "variant_b_median_rho_nppa": b_nppa,
        "ambient_ventricle": ambv,
        "doublet_median_abs_rho_thinned_sg": dbl,
    }
    for ax in axes:
        t = ties[ties.axis == ax]
        summ["ties_median"][ax] = {
            c: round(float(t[c].median()), 3) for c in t.columns if "_rho_" in c
        }
        summ["ties_median"][ax]["ucell_abs_depth_gt_sg_unthinned_libs"] = int(
            (t.ucell_rho_depth.abs() > t.sg_unthinned_rho_depth.abs()).sum()
        )
        n = nulls[nulls.axis == ax]
        summ["null"][ax] = {
            "libs_p90_above": int(n.p90_above.sum()),
            "libs_coherence_above": int(n.coherence_above.sum()),
            "hearts_p90_above": hearts_all(n, "p90_above"),
            "hearts_coherence_above": hearts_all(n, "coherence_above"),
        }
    for (v, pr), x in pg.groupby(["variant", "pair"]):
        summ["program"][f"{v}:{pr}"] = {
            "libs_above": int(x.above.sum()),
            "n_libs": len(x),
            "hearts_above": hearts_all(x, "above"),
            "median_rho": round(float(x.rho.median()), 3),
            "median_null_p95": round(float(x.null_p95.median()), 3),
        }
    summ["meta"] = {
        "commit": subprocess.check_output(["git", "rev-parse", "HEAD"], text=True).strip(),
        "dirty_tree": True,
        "command": "uv run python " + " ".join(sys.argv),
        "input": src,
        "marker_dir": marker_dir,
        "time": datetime.now().astimezone().isoformat(timespec="seconds"),
    }
    (out / "c008_summary.json").write_text(json.dumps(summ, indent=1, default=str))
    print(json.dumps({k: v for k, v in summ.items() if k != "meta"}, indent=1, default=str))


if __name__ == "__main__":
    main(*sys.argv[1:4])
