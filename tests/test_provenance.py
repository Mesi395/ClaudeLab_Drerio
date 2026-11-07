import pytest
import yaml

from zf import provenance


def test_manifest_records_inputs_outcomes_and_commit(tmp_path):
    data = tmp_path / "counts.h5"
    data.write_bytes(b"abc")
    path = provenance.write_manifest(
        tmp_path / "run1",
        command="uv run python scripts/run_checks.py",
        inputs=[data],
        params={"seed": 0},
        outcomes={"border_zone": {"library_size": "flag", "doublets": "not_run"}},
    )
    m = yaml.safe_load(path.read_text())
    assert m["schema"] == "zf.run/v1"
    assert m["inputs"][0]["sha256"].startswith("ba7816bf")
    assert m["outcomes"]["border_zone"]["doublets"] == "not_run"
    assert "commit" in m["git"]


def test_unknown_outcome_is_rejected(tmp_path):
    with pytest.raises(ValueError, match="passed"):
        provenance.write_manifest(
            tmp_path / "run2",
            command="x",
            inputs=[],
            params={},
            outcomes={"border_zone": {"library_size": "passed"}},
        )
