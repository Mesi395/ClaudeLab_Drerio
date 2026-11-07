from pathlib import Path

import pytest
import yaml

from zf import guard

ROOT = Path(__file__).resolve().parents[1]
MARKERS = ROOT / "markers/border_zone.yaml"


def write(path, content):
    return guard.decide("Write", {"file_path": str(path), "content": content}, ROOT, env={})


def test_unrelated_file_is_allowed():
    assert write(ROOT / "src/zf/score.py", "x = 1\n").allowed


def test_holdout_is_blocked():
    d = write(ROOT / "holdout/traps.yaml", "a: 1\n")
    assert not d.allowed and "holdout" in d.reason


def test_invalid_marker_file_is_blocked():
    d = write(MARKERS, "schema: zf.markers/v1\naxis: border_zone\n")
    assert not d.allowed and "schema" in d.reason


def test_marker_change_outside_flags_is_blocked():
    data = yaml.safe_load(MARKERS.read_text())
    data["genes"][0]["status"] = "accepted"
    d = write(MARKERS, yaml.safe_dump(data, sort_keys=False))
    assert not d.allowed and "flags" in d.reason


def test_adding_a_flag_is_allowed():
    data = yaml.safe_load(MARKERS.read_text())
    data["genes"][0]["flags"].append("library_size: results/runs/r1/manifest.yaml")
    assert write(MARKERS, yaml.safe_dump(data, sort_keys=False)).allowed


def test_human_override_allows_marker_edit():
    data = yaml.safe_load(MARKERS.read_text())
    data["genes"][0]["status"] = "accepted"
    d = guard.decide(
        "Write",
        {"file_path": str(MARKERS), "content": yaml.safe_dump(data, sort_keys=False)},
        ROOT,
        env={"ZF_HUMAN_EDIT": "1"},
    )
    assert d.allowed


def test_human_config_is_blocked_for_the_agent():
    d = write(ROOT / "config/checks.yaml", (ROOT / "config/checks.yaml").read_text())
    assert not d.allowed and "ZF_HUMAN_EDIT" in d.reason


def test_frozen_checks_are_blocked_even_for_humans():
    text = (
        (ROOT / "config/checks.yaml").read_text().replace("status: proposed", "status: frozen", 1)
    )
    current = guard.decide(
        "Write",
        {"file_path": str(ROOT / "config/checks.yaml"), "content": text},
        ROOT,
        env={"ZF_HUMAN_EDIT": "1"},
        current_text=text,
    )
    assert not current.allowed and "frozen" in current.reason


def test_manifest_must_come_from_write_manifest():
    d = write(ROOT / "results/runs/r1/manifest.yaml", "schema: zf.run/v1\n")
    assert not d.allowed and "write_manifest" in d.reason


def test_edit_is_validated_on_the_resulting_file():
    d = guard.decide(
        "Edit",
        {"file_path": str(MARKERS), "old_string": "role: core", "new_string": "role: maybe"},
        ROOT,
        env={"ZF_HUMAN_EDIT": "1"},
    )
    assert not d.allowed


@pytest.mark.parametrize("tool", ["Read", "Bash"])
def test_other_tools_pass_through(tool):
    assert guard.decide(tool, {}, ROOT, env={}).allowed
