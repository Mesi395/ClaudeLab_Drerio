import pytest

from zf import commit_rules


@pytest.mark.parametrize(
    "path",
    [
        "holdout/traps.yaml",
        "data/raw/H5_filtered_matrix.h5",
        "results/runs/r1/cells.h5ad",
        "notes/todo.md",
        "airdrop/HR_1.csv",
        "CLAUDE.local.md",
        ".claude/settings.local.json",
        ".env",
        ".env.local",
        "keys/deploy.pem",
    ],
)
def test_protected_paths_cannot_be_committed(path):
    assert commit_rules.reason(path)


@pytest.mark.parametrize(
    "path",
    ["src/zf/score.py", "markers/border_zone.yaml", "results/runs/r1/manifest.yaml", "README.md"],
)
def test_project_files_can_be_committed(path):
    assert commit_rules.reason(path) is None


def test_main_is_protected_from_direct_push():
    assert commit_rules.push_reason("refs/heads/main")
    assert commit_rules.push_reason("refs/heads/ad/score-axes") is None
