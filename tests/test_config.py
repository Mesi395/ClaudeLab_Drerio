from pathlib import Path

import pytest
import yaml
from pydantic import ValidationError

from zf import config

ROOT = Path(__file__).resolve().parents[1]


def test_dataset_manifest_matches_expected_totals():
    ds = config.load_dataset(ROOT / "config/dataset.yaml")
    assert len(ds.libraries) == ds.expected.libraries
    assert sum(lib.expected_cells for lib in ds.libraries) == ds.expected.cells


def test_every_marker_file_validates_and_has_unique_ensembl_ids():
    for path in sorted((ROOT / "markers").glob("*.yaml")):
        markers = config.load_markers(path)
        ids = [g.ensembl for g in markers.genes]
        assert len(ids) == len(set(ids)), path.name


def test_core_marker_sets_are_disjoint_across_axes():
    seen: dict[str, str] = {}
    for path in sorted((ROOT / "markers").glob("*.yaml")):
        markers = config.load_markers(path)
        for gene in markers.genes:
            if gene.role == "core":
                assert gene.ensembl not in seen, f"{gene.symbol} in {seen.get(gene.ensembl)}"
                seen[gene.ensembl] = markers.axis


def test_axes_point_at_existing_marker_files_and_conditions():
    axes = config.load_axes(ROOT / "config/axes.yaml")
    conditions = {
        lib.condition for lib in config.load_dataset(ROOT / "config/dataset.yaml").libraries
    }
    for axis in axes.axes:
        assert (ROOT / axis.markers).exists()
        for exp in axis.expect:
            assert {exp.higher, exp.than} <= conditions


def test_checks_have_unique_ids():
    checks = config.load_checks(ROOT / "config/checks.yaml")
    ids = [c.id for c in checks.checks]
    assert len(ids) == len(set(ids))


def test_marker_without_source_is_rejected(tmp_path):
    bad = {
        "schema": "zf.markers/v1",
        "axis": "x",
        "annotation": "a",
        "genes": [{"symbol": "nppa", "ensembl": "ENSDARG00000052960", "role": "core"}],
    }
    path = tmp_path / "bad.yaml"
    path.write_text(yaml.safe_dump(bad))
    with pytest.raises(ValidationError):
        config.load_markers(path)
