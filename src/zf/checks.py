"""The gating checks of config/checks.yaml. Each returns a statistic and pass or flag."""

import numpy as np
from scipy import stats

from zf.panel import Panel
from zf.score import expression_bins, score_mean, score_rank

N_NULL = 200
CONFOUNDS = ("depth", "doublet", "ambient", "dissociation")


def _rank(x):
    return stats.rankdata(x)


def _residual(y, X):
    X = np.column_stack([np.ones(len(y)), X])
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    return y - X @ beta


def _within_condition_variance(panel, score):
    return np.mean([score[panel.condition == c].var() for c in np.unique(panel.condition)])


def null_matched_sets(panel: Panel, seed: int):
    rng = np.random.default_rng(seed)
    idx = panel.idx(panel.axis_genes)
    bins = expression_bins(panel)
    real = _within_condition_variance(panel, score_mean(panel, panel.axis_genes))
    null = []
    for _ in range(N_NULL):
        pick = [rng.choice(np.setdiff1d(np.flatnonzero(bins == bins[i]), idx)) for i in idx]
        null.append(
            _within_condition_variance(panel, score_mean(panel, [panel.genes[j] for j in pick]))
        )
    z = (real - np.mean(null)) / np.std(null)
    return {"z": float(z), "outcome": "flag" if z < 3 else "pass"}


def direction(panel: Panel, seed: int):
    if not panel.expect:
        return {"outcome": "not_run", "reason": "no expected direction for this axis"}
    score = score_mean(panel, panel.axis_genes)
    libs = np.unique(panel.library)
    lib_mean = {lib: score[panel.library == lib].mean() for lib in libs}
    lib_cond = {lib: panel.condition[panel.library == lib][0] for lib in libs}
    results = []
    for higher, than in panel.expect:
        a = [lib_mean[lib] for lib in libs if lib_cond[lib] == higher]
        b = [lib_mean[lib] for lib in libs if lib_cond[lib] == than]
        p = stats.mannwhitneyu(a, b, alternative="greater").pvalue
        results.append(float(p))
    worst = max(results) if results else 1.0
    return {"p": results, "outcome": "flag" if worst > 0.05 else "pass"}


def _partial(panel: Panel, target: str, threshold: float):
    if target not in panel.covariates:
        return {"outcome": "not_run", "reason": f"no {target} covariate"}
    score = score_mean(panel, panel.axis_genes)
    others = [c for c in CONFOUNDS if c != target and c in panel.covariates]
    rhos = []
    for lib in np.unique(panel.library):
        m = panel.library == lib
        controls = np.column_stack([_rank(panel.covariates[c][m]) for c in others])
        rs = _residual(_rank(score[m]), controls)
        rt = _residual(_rank(panel.covariates[target][m]), controls)
        rhos.append(float(np.corrcoef(rs, rt)[0, 1]))
    hits = sum(abs(r) > threshold for r in rhos)
    return {"rho": rhos, "outcome": "flag" if hits >= len(rhos) / 2 else "pass"}


def depth(panel, seed):
    return _partial(panel, "depth", 0.3)


def doublet(panel, seed):
    return _partial(panel, "doublet", 0.2)


def ambient(panel, seed):
    return _partial(panel, "ambient", 0.2)


def dissociation(panel, seed):
    return _partial(panel, "dissociation", 0.2)


def batch(panel: Panel, seed: int):
    score = score_mean(panel, panel.axis_genes)
    present = [c for c in CONFOUNDS if c in panel.covariates]
    resid = _residual(score, np.column_stack([panel.covariates[c] for c in present]))
    etas = {}
    for cond in np.unique(panel.condition):
        m = panel.condition == cond
        libs = np.unique(panel.library[m])
        if len(libs) < 2:
            continue
        r, grand = resid[m], resid[m].mean()
        between = sum(
            (panel.library[m] == lib).sum() * (r[panel.library[m] == lib].mean() - grand) ** 2
            for lib in libs
        )
        etas[str(cond)] = float(between / ((r - grand) ** 2).sum())
    return {"eta2": etas, "outcome": "flag" if any(e > 0.1 for e in etas.values()) else "pass"}


def leave_one_gene_out(panel: Panel, seed: int):
    full = score_mean(panel, panel.axis_genes)
    rhos = {
        g: float(
            stats.spearmanr(full, score_mean(panel, [h for h in panel.axis_genes if h != g]))[0]
        )
        for g in panel.axis_genes
    }
    return {"rho": rhos, "outcome": "flag" if min(rhos.values()) < 0.8 else "pass"}


def method_agreement(panel: Panel, seed: int):
    rho = float(
        stats.spearmanr(score_mean(panel, panel.axis_genes), score_rank(panel, panel.axis_genes))[0]
    )
    return {"rho": rho, "outcome": "flag" if rho < 0.6 else "pass"}


CHECKS = {
    "null_matched_sets": null_matched_sets,
    "direction": direction,
    "depth": depth,
    "doublet": doublet,
    "ambient": ambient,
    "dissociation": dissociation,
    "batch": batch,
    "leave_one_gene_out": leave_one_gene_out,
    "method_agreement": method_agreement,
}


def evaluate(check_id: str, panel: Panel, seed: int = 0) -> dict:
    return CHECKS[check_id](panel, seed)


def run(check_id: str, panel: Panel, seed: int = 0) -> str:
    return evaluate(check_id, panel, seed)["outcome"]


SIGNAL_GATED = ("depth", "doublet", "ambient", "dissociation", "batch")
SIGNAL_GATED += ("leave_one_gene_out", "method_agreement")


def run_all(panel: Panel, seed: int = 0) -> dict[str, str]:
    """Every check in order. If the null check flags, the axis has no signal of its own,
    so confound and robustness flags cannot be read: those checks report not_run."""
    out = {"null_matched_sets": run("null_matched_sets", panel, seed)}
    out["direction"] = run("direction", panel, seed)
    for check_id in SIGNAL_GATED:
        no_signal = out["null_matched_sets"] == "flag"
        out[check_id] = "not_run" if no_signal else run(check_id, panel, seed)
    return out
