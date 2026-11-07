"""The input every check reads: counts, one axis gene set, and per-cell design and covariates."""

from dataclasses import dataclass, field

import numpy as np
import scipy.sparse as sp


@dataclass
class Panel:
    counts: np.ndarray | sp.spmatrix  # cells x genes, raw counts
    genes: list[str]
    axis_genes: list[str]
    library: np.ndarray
    condition: np.ndarray
    covariates: dict[str, np.ndarray]  # depth, doublet, ambient, dissociation
    expect: list[tuple[str, str]] = field(default_factory=list)  # (higher, than) conditions
    total: np.ndarray | None = None  # per-cell total over all genes, when counts is a gene subset

    def __post_init__(self):
        dense = self.counts.toarray() if sp.issparse(self.counts) else np.asarray(self.counts)
        if self.total is None:
            total = dense.sum(axis=1, keepdims=True)
        else:
            total = np.asarray(self.total, dtype=float).reshape(-1, 1)
        self.logn = np.log1p(dense / np.maximum(total, 1) * 1e4)
        self.index = {g: i for i, g in enumerate(self.genes)}
        self.bins = None  # expression bins, computed once by zf.score

    def idx(self, genes) -> np.ndarray:
        return np.array([self.index[g] for g in genes])
