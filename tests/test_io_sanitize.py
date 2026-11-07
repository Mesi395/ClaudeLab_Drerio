"""Leakage guard: the sanitised copy keeps only the allowlist (R5.3, study cm-state-axes C001)."""

import anndata as ad
import h5py
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

from zf.io import ALLOWED_OBS, sanitize_h5ad

LABEL_COLS = ["annot_a", "annot_b", "subtype_x", "clusters", "flag_x", "annot_c"]  # stand-ins


def _make_labelled(path):
    rng = np.random.default_rng(0)
    n, g = 12, 5
    counts = sp.csr_matrix(rng.poisson(1.0, size=(n, g)).astype(np.float32))
    obs = pd.DataFrame(
        {
            "orig.ident": pd.Categorical(["Hr1"] * 6 + ["Hr11"] * 6),
            "nCount_RNA": np.asarray(counts.sum(1)).ravel(),
            "nFeature_RNA": np.asarray((counts > 0).sum(1)).ravel(),
            "percent.mito": rng.random(n),
            "time": pd.Categorical(["Ctrl"] * 6 + ["3dpi"] * 6),
            "morphine": pd.Categorical(["no"] * n),
            "AV": pd.Categorical(["AV"] * n),
        },
        index=[f"c{i}" for i in range(n)],
    )
    for c in LABEL_COLS:
        obs[c] = pd.Categorical(["SECRET_LABEL"] * n)
    var = pd.DataFrame({"means": rng.random(g)}, index=[f"g{i}" for i in range(g)])
    a = ad.AnnData(X=np.log1p(counts), obs=obs, var=var)
    a.layers["counts"] = counts
    a.uns["clusters"] = {"params": {"resolution": 1.0}}
    a.obsm["X_umap"] = rng.random((n, 2))
    a.write_h5ad(path)
    return counts


def test_sanitize_keeps_only_allowlist(tmp_path):
    src, dst = tmp_path / "in.h5ad", tmp_path / "out.h5ad"
    counts = _make_labelled(src)
    sanitize_h5ad(src, dst)
    out = ad.read_h5ad(dst)
    assert list(out.obs.columns) == list(ALLOWED_OBS)
    assert out.uns == {} or "clusters" not in out.uns
    assert len(out.obsm) == 0
    # anndata 0.13 lists X as layers[None]; only named layer must be counts
    assert [k for k in out.layers if k is not None] == ["counts"]
    assert (out.layers["counts"] != counts).nnz == 0
    assert (counts != out.X).nnz == 0  # X is raw counts, not the authors' normalisation
    assert list(out.var_names) == [f"g{i}" for i in range(5)]
    assert "means" in out.var.columns
    with h5py.File(dst, "r") as f:
        # no label string anywhere in the file
        found = []
        f.visititems(lambda name, obj: found.append(name) if "SECRET" in name else None)
        assert not found
    assert b"SECRET_LABEL" not in dst.read_bytes()


def test_sanitize_refuses_missing_allowlisted_column(tmp_path):
    src, dst = tmp_path / "in.h5ad", tmp_path / "out.h5ad"
    _make_labelled(src)
    a = ad.read_h5ad(src)
    del a.obs["percent.mito"]
    a.write_h5ad(src)
    with pytest.raises(KeyError):
        sanitize_h5ad(src, dst)
