"""Which annotation to trust, with no ground truth (Card C, steps 1-2).

Fits a two-class latent-class (Dawid-Skene) model per axis over five imperfect labellers:
  a  protocol v1, UCell arm      (p_ucell < 0.05 -> 1, > 0.2 -> 0, else missing)
  b  protocol v1, score_genes arm on thinned counts (p_score, same cut-offs)
  c  split-half marker set 1     (scripts/split_half_calls.py; uncertain -> missing)
  d  split-half marker set 2     (disjoint genes from c)
  e  the annotation under test   (cells in the matching cluster -> 1, all others -> 0)
and reports each labeller's sensitivity/specificity with a bootstrap over HEARTS, the same fit
within depth tertiles (depth breaks conditional independence), and Cohen's kappa + confusion
fractions of our high-confidence calls vs the annotation per timepoint.

It prints summary numbers ONLY. It never writes per-cell annotation anywhere. The annotation
category strings are passed on the command line, so none are stored in this repo.

Usage:
  python scripts/compare_to_annotation.py LABELS SPLITHALF ANNOTATION \
      [--column NAME] --dediff SUBSTR --prolif SUBSTR [--time-column C] [--depth-column C]
  LABELS     data/derived/labels_v1.tsv.gz (or .parquet)    SPLITHALF data/derived/splithalf_v1.tsv.gz
  ANNOTATION a .csv with columns cell,label[,time,depth] or an .h5ad read in memory with --column
"""

import argparse

import numpy as np
import pandas as pd

AXES = {"dedifferentiated": "dediff", "proliferating": "prolif"}
RNG = np.random.default_rng(0)


def read_table(path):
    return pd.read_parquet(path) if str(path).endswith(".parquet") else pd.read_csv(path, sep="\t", index_col=0)


def read_annotation(args) -> pd.DataFrame:
    if args.annotation.endswith(".h5ad"):
        import anndata as ad

        obs = ad.read_h5ad(args.annotation, backed="r").obs
        out = pd.DataFrame({"label": obs[args.column].astype(str)}, index=obs.index)
        for k, c in (("time", args.time_column), ("depth", args.depth_column)):
            if c in obs:
                out[k] = obs[c].to_numpy()
        return out
    return pd.read_csv(args.annotation, index_col=0)


def p_to_call(p):
    return np.where(p < 0.05, 1.0, np.where(p > 0.2, 0.0, np.nan))


def dawid_skene(X: np.ndarray, iters: int = 100):
    """X: cells x labellers, values 1/0/nan. Returns prevalence, sensitivity[j], specificity[j]."""
    obs = ~np.isnan(X)
    X0 = np.nan_to_num(X)
    t = np.nanmean(X, axis=1)  # majority-vote start
    t = np.nan_to_num(t, nan=0.0)
    for _ in range(iters):
        pi = t.mean()
        se = ((t[:, None] * X0 * obs).sum(0) + 1) / ((t[:, None] * obs).sum(0) + 2)
        sp = (((1 - t)[:, None] * (1 - X0) * obs).sum(0) + 1) / (((1 - t)[:, None] * obs).sum(0) + 2)
        l1 = np.log(pi + 1e-9) + (obs * (X0 * np.log(se) + (1 - X0) * np.log(1 - se))).sum(1)
        l0 = np.log(1 - pi + 1e-9) + (obs * ((1 - X0) * np.log(sp) + X0 * np.log(1 - sp))).sum(1)
        t = 1 / (1 + np.exp(l0 - l1))
    return pi, se, sp


def kappa(y1, y2) -> float:
    po = np.mean(y1 == y2)
    pe = np.mean(y1) * np.mean(y2) + np.mean(1 - y1) * np.mean(1 - y2)
    return float((po - pe) / (1 - pe)) if pe < 1 else float("nan")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("labels")
    ap.add_argument("splithalf")
    ap.add_argument("annotation")
    ap.add_argument("--column", default="label")
    ap.add_argument("--dediff", required=True, help="substring that marks the annotation's dediff cluster")
    ap.add_argument("--prolif", required=True, help="substring that marks the annotation's proliferating cluster")
    ap.add_argument("--time-column", default="time")
    ap.add_argument("--depth-column", default="nCount_RNA")
    ap.add_argument("--boot", type=int, default=20)
    args = ap.parse_args()

    lab, sh, ann = read_table(args.labels), read_table(args.splithalf), read_annotation(args)
    cells = lab.index.intersection(sh.index).intersection(ann.index)
    lab, sh, ann = lab.loc[cells], sh.loc[cells], ann.loc[cells]
    heart = pd.Index(cells).str.split("_").str[0].to_numpy()
    print(f"cells matched: {len(cells)}  hearts: {len(set(heart))}")
    names = ["a_ucell", "b_score", "c_half1", "d_half2", "e_annotation"]

    for axis, key in AXES.items():
        sub = getattr(args, key)
        e = ann["label"].str.contains(sub, regex=False).to_numpy().astype(float)
        half = lambda h: sh[f"{axis}.{h}.call"].map({"positive": 1.0, "negative": 0.0}).to_numpy() if f"{axis}.{h}.call" in sh else np.full(len(cells), np.nan)
        X = np.column_stack([p_to_call(lab[f"{axis}.p_ucell"].to_numpy()), p_to_call(lab[f"{axis}.p_score"].to_numpy()),
                             half("half1"), half("half2"), e])
        keep = [j for j in range(X.shape[1]) if np.nansum(~np.isnan(X[:, j])) > 0]
        X, nm = X[:, keep], [names[j] for j in keep]
        print(f"\n== {axis}: labellers {nm}; positives per labeller {np.nansum(X, 0).astype(int).tolist()}")
        pi, se, sp = dawid_skene(X)
        hs = np.unique(heart)
        boots = []
        for _ in range(args.boot):  # resample hearts, not cells
            idx = np.concatenate([np.flatnonzero(heart == h) for h in RNG.choice(hs, len(hs))])
            boots.append(np.r_[dawid_skene(X[idx])[1], dawid_skene(X[idx])[2]])
        lo, hi = np.percentile(boots, [5, 95], axis=0)
        k = len(nm)
        print(f"latent prevalence {pi:.4f}")
        for j, n in enumerate(nm):
            print(f"  {n:13s} sens {se[j]:.3f} [{lo[j]:.3f},{hi[j]:.3f}]  spec {sp[j]:.4f} [{lo[k+j]:.4f},{hi[k+j]:.4f}]")
        if "depth" in ann:
            tert = pd.qcut(ann["depth"].rank(method="first"), 3, labels=["low", "mid", "high"]).to_numpy()
            for t in ["low", "mid", "high"]:
                m = tert == t
                _, se_t, sp_t = dawid_skene(X[m])
                print(f"  depth {t:4s}: sens " + " ".join(f"{n[0]}={v:.2f}" for n, v in zip(nm, se_t))
                      + " | spec " + " ".join(f"{n[0]}={v:.3f}" for n, v in zip(nm, sp_t)))
            d = pd.DataFrame({"heart": heart, "depth": ann["depth"].to_numpy(), "e": e})
            rel = d["depth"] / d.groupby("heart")["depth"].transform("median")
            print(f"  annotation-positive cells: median depth / own-heart median = {rel[d['e'] == 1].median():.2f}"
                  f" (annotation-negative: {rel[d['e'] == 0].median():.2f})")
        # our high-confidence calls vs the annotation
        call, hc = lab[f"{axis}.call"].to_numpy(), lab[f"{axis}.high_confidence"].astype(bool).to_numpy()
        ours = np.where(call == "positive", 1.0, np.where(call == "negative", 0.0, np.nan))
        m = ~np.isnan(ours) & (hc | (call == "negative"))
        times = ann["time"].astype(str).to_numpy() if "time" in ann else np.full(len(cells), "all")
        for t in ["all"] + sorted(set(times) - {"all"}):
            mm = m & ((times == t) | (t == "all"))
            if mm.sum() == 0:
                continue
            y1, y2 = ours[mm], e[mm]
            cf = [np.mean((y1 == i) & (y2 == j)) for i in (1, 0) for j in (1, 0)]
            print(f"  kappa[{t:5s}] {kappa(y1, y2):+.3f}  n={mm.sum():6d}  confusion (ours+/ann+, ours+/ann-, ours-/ann+, ours-/ann-) "
                  + " ".join(f"{v:.4f}" for v in cf))


if __name__ == "__main__":
    main()
