import importlib.util
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

spec = importlib.util.spec_from_file_location(
    "make_agent_input", Path(__file__).resolve().parents[1] / "scripts/make_agent_input.py"
)
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)


def toy():
    obs = pd.DataFrame(
        {
            "orig.ident": ["H5", "Hr11", "Hr11"],
            "nCount_RNA": [900.0, 700.0, 650.0],
            "nFeature_RNA": [400, 350, 300],
            "percent.mito": [5.0, 6.0, 7.0],
            "time": ["Ctrl", "3dpi", "3dpi"],
            "AV": ["Wholeheart"] * 3,
            "some_label": ["x", "y", "z"],
            "another_label": ["a", "b", "c"],
        },
        index=["c1", "c2", "c3"],
    )
    a = ad.AnnData(X=np.ones((3, 2)), obs=obs, var=pd.DataFrame(index=["nppa", "myl7"]))
    a.layers["counts"] = np.array([[1, 2], [3, 4], [5, 6]])
    a.obsm["X_umap"] = np.zeros((3, 2))
    a.uns["anything"] = {"k": 1}
    return a


def test_only_allowlisted_columns_survive():
    out = mod.strip(toy())
    assert list(out.obs.columns) == [c for c in mod.KEEP_OBS if c in toy().obs.columns]


def test_counts_kept_and_everything_derived_dropped():
    out = mod.strip(toy())
    np.testing.assert_array_equal(out.X, toy().layers["counts"])
    assert not out.obsm and not out.uns and not out.layers
    assert list(out.var.columns) == []


def test_cell_list_has_barcode_and_library_only():
    cells = mod.cell_list(toy())
    assert list(cells.columns) == ["cell", "library"]
    assert len(cells) == 3
