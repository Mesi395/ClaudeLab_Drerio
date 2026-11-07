"""The contract for the checks: a diagonal detection matrix on synthetic data.

For each defect named in config/checks.yaml, zf.synthetic.make(defect, seed) builds a dataset
with that one defect injected. Every gating check runs on it.
The check that `detects` the defect must flag; every other check must not.
A clean dataset must pass every check.

If two checks keep firing on the same defect, they are not orthogonal: merge them or drop one,
then change config/checks.yaml by hand.
"""

from pathlib import Path

import pytest

from zf import config

checks_impl = pytest.importorskip("zf.checks", reason="zf.checks not implemented yet")
synthetic = pytest.importorskip("zf.synthetic", reason="zf.synthetic not implemented yet")

ROOT = Path(__file__).resolve().parents[1]
CHECKS = config.load_checks(ROOT / "config/checks.yaml")
GATING = [c for c in CHECKS.checks if c.family != "report"]


def flagged(data, seed) -> set[str]:
    return {c.id for c in GATING if checks_impl.run(c.id, data, seed=seed) == "flag"}


@pytest.mark.parametrize("seed", CHECKS.seeds)
def test_clean_data_passes_every_check(seed):
    assert flagged(synthetic.make(None, seed=seed), seed) == set()


@pytest.mark.parametrize("seed", CHECKS.seeds)
@pytest.mark.parametrize("check", GATING, ids=lambda c: c.detects)
def test_each_defect_raises_exactly_its_own_flag(check, seed):
    assert flagged(synthetic.make(check.detects, seed=seed), seed) == {check.id}
