"""Synthetic panels with one planted defect each: the ground truth the checks are evaluated on.

Clean data: a latent state, higher after injury, drives a co-expressed axis gene set.
Each defect adds exactly one artefact named in config/checks.yaml (`detects`).
"""

import numpy as np

from zf.panel import Panel
from zf.score import score_mean

N_GENES = 300
N_AXIS = 10
CELLS_PER_LIBRARY = 300
LIBRARIES = {"uninjured": 4, "dpi_7": 4}

AXIS = [f"axis{i}" for i in range(N_AXIS)]
AMBIENT = [f"fibro{i}" for i in range(10)]
DISSOCIATION = [f"ieg{i}" for i in range(10)]


def make(defect: str | None, seed: int = 0) -> Panel:
    rng = np.random.default_rng(seed)
    other = [f"g{i}" for i in range(N_GENES - len(AXIS) - len(AMBIENT) - len(DISSOCIATION))]
    genes = AXIS + AMBIENT + DISSOCIATION + other
    library, condition = [], []
    for cond, n_lib in LIBRARIES.items():
        for k in range(n_lib):
            library += [f"{cond}_{k}"] * CELLS_PER_LIBRARY
            condition += [cond] * CELLS_PER_LIBRARY
    library, condition = np.array(library), np.array(condition)
    n = len(library)
    injured = (condition != "uninjured").astype(float)

    base = rng.gamma(2.0, 1.0, len(genes)) * 0.5
    size = rng.lognormal(np.log(2000), 0.5, n)
    state = rng.normal(0, 1, n) + (0.0 if defect == "shuffled_conditions" else 1.0) * injured
    doublet = rng.normal(0, 1, n)
    ambient = rng.normal(0, 1, n)
    stress = rng.normal(0, 1, n)

    log_mu = np.tile(np.log(base), (n, 1))
    a = np.arange(N_AXIS)
    if defect == "random_gene_set":
        # each gene shifts with condition, but genes vary independently within a condition
        log_mu[:, a] += 0.8 * injured[:, None] + 0.15 * rng.normal(0, 1, (n, N_AXIS))
    elif defect == "single_gene_driver":
        log_mu[:, 0] += 2.0 * state
        log_mu[:, a[1:]] += 0.05 * state[:, None]
    else:
        log_mu[:, a] += 0.6 * state[:, None]

    log_mu[:, N_AXIS : N_AXIS + 10] += 0.8 * ambient[:, None]
    log_mu[:, N_AXIS + 10 : N_AXIS + 20] += 0.8 * stress[:, None]

    if defect == "depth_artefact":
        log_mu[:, a] += 1.0 * (np.log(size) - np.log(size).mean())[:, None]
    elif defect == "doublet_artefact":
        log_mu[:, a] += 0.6 * doublet[:, None]
    elif defect == "ambient_artefact":
        log_mu[:, a] += 0.6 * ambient[:, None]
    elif defect == "dissociation_artefact":
        log_mu[:, a] += 0.6 * stress[:, None]
    elif defect == "library_batch":
        log_mu[np.ix_(library == "dpi_7_0", a)] += 1.0
    elif defect == "method_artefact":
        # half the axis genes are far more abundant and vary independently of the state:
        # they dominate the mean-based score but always rank at the top, so the rank-based
        # score follows the other half
        noisy = a[: N_AXIS // 2]
        log_mu[:, noisy] = np.log(base.max() * 20) + 2.5 * rng.normal(0, 1, (n, len(noisy)))

    mu = np.exp(log_mu)
    mu = mu / mu.sum(axis=1, keepdims=True) * size[:, None]
    counts = rng.poisson(mu)

    panel = Panel(
        counts=counts,
        genes=genes,
        axis_genes=AXIS,
        library=library,
        condition=condition,
        covariates={},
        expect=[("dpi_7", "uninjured")],
    )
    panel.covariates = {
        "depth": np.log(counts.sum(axis=1)),
        "doublet": doublet + rng.normal(0, 0.3, n),
        "ambient": score_mean(panel, AMBIENT),
        "dissociation": score_mean(panel, DISSOCIATION),
    }
    return panel
