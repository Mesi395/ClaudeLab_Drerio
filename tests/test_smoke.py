import numpy as np
import scipy.sparse as sp

import zf
from zf import confounds, io


def test_package_imports():
    assert zf.__doc__


def test_hearts_per_condition():
    t = io.library_table()
    assert len(t) == 22
    hearts = t.groupby("condition")["heart"].nunique().to_dict()
    assert hearts == {"uninjured": 4, "dpi_3": 4, "dpi_7": 7, "dpi_30": 3}
    assert set(t.loc[t["chamber"] != "unknown"].index) == {
        "H8a", "H8v", "Hr6a", "Hr6v", "Hr7a", "Hr7v"
    }  # fmt: skip


def test_ambient_fraction_recovers_known_contamination():
    rng = np.random.default_rng(0)
    profile = rng.dirichlet(np.ones(50)) * 0.3  # 30% of soup counts fall on these genes
    rho = rng.uniform(0, 0.4, 2000)
    total = rng.integers(500, 5000, 2000)
    x = rng.poisson(np.outer(rho * total, profile))
    est = confounds.ambient_fraction(sp.csr_matrix(x), total, profile)
    assert np.corrcoef(est, rho)[0, 1] > 0.9
    assert abs(np.median(est - rho)) < 0.02
    clean = confounds.ambient_fraction(sp.csr_matrix(np.zeros((5, 50))), total[:5], profile)
    assert (clean == 0).all()
