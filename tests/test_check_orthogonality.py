"""The contract for the checks: a diagonal detection matrix on synthetic data.

For each defect named in config/checks.yaml, zf.synthetic.make(defect, seed) builds a dataset
with that one defect injected. zf.checks.run_all runs every gating check, signal checks first.
Over several seeds, the check that `detects` the defect must flag, and only it.
A clean dataset must never raise a flag.

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
SEEDS = range(5)
MIN_EXACT = 4  # of 5 seeds: the defect raises exactly its own flag


def flagged(defect, seed) -> set[str]:
    outcomes = checks_impl.run_all(synthetic.make(defect, seed=seed), seed=seed)
    return {k for k, v in outcomes.items() if v == "flag"}


def test_clean_data_never_raises_a_flag():
    for seed in SEEDS:
        assert flagged(None, seed) == set()


@pytest.mark.parametrize("check", GATING, ids=lambda c: c.detects)
def test_each_defect_raises_exactly_its_own_flag(check):
    exact = sum(flagged(check.detects, seed) == {check.id} for seed in SEEDS)
    assert exact >= MIN_EXACT
