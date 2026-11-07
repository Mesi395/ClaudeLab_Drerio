"""Card A gate: marker sets v1 are usable and the trust card exists before any score.

Every gene is in the data, core sets of the four axes are disjoint, every gene cites a DOI,
and JOURNAL.md carries the pre-registered trust card. Exit 1 on any failure.
Run where the label-free data exists: uv run python scripts/check_markers_v1.py
"""

import re
import sys
from pathlib import Path

import anndata as ad

from zf import config, io

ROOT = Path(__file__).resolve().parents[1]
SETS = ROOT / "results" / "hk" / "markers_v1"
AXES = ["injury_responsive", "border_zone", "dedifferentiated", "proliferating"]
TRUST_HEADER = "## Pre-registered trust card"


def main() -> int:
    genes = set(ad.read_h5ad(io.LABEL_FREE, backed="r").var_names)
    errors, core = [], {}
    for path in sorted(SETS.glob("*.yaml")):
        m = config.load_markers(path)
        for g in m.genes:
            if g.symbol not in genes:
                errors.append(f"{path.name}: {g.symbol} not in var_names")
            if not re.match(r"^https://doi\.org/10\.\S+$", g.source):
                errors.append(f"{path.name}: {g.symbol} has no DOI ({g.source!r})")
            if path.stem in AXES and g.role == "core":
                core.setdefault(g.symbol, []).append(path.stem)
    missing = [a for a in AXES if not (SETS / f"{a}.yaml").exists()]
    errors += [f"missing set {a}.yaml" for a in missing]
    errors += [f"{s} is core in {v}" for s, v in core.items() if len(v) > 1]
    if TRUST_HEADER not in (ROOT / "JOURNAL.md").read_text():
        errors.append("JOURNAL.md has no trust card")
    for e in errors:
        print("FAIL", e)
    print(f"{'FAIL' if errors else 'OK'}: {len(core)} core genes over {len(AXES)} axes")
    return 1 if errors else 0


if __name__ == "__main__":
    sys.exit(main())
