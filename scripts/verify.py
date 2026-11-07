"""Deterministic verification of the whole repository. Run before every commit and before a claim.

Exit codes: 0 every check ran and passed; 1 a check failed; 2 a check could not run.
A check that could not run is never reported as passed.
"""

import subprocess
import sys
from pathlib import Path

import yaml

from zf import config, provenance

ROOT = Path(__file__).resolve().parents[1]
# Files that name holdout/ only to block it; evaluate.py is its one reader.
NAMES_HOLDOUT_WITHOUT_READING = {"evaluate.py", "guard.py", "verify.py", "commit_rules.py"}
results: list[tuple[str, str, str]] = []  # (check, outcome, detail)


def record(name, outcome, detail=""):
    results.append((name, outcome, detail))


def run(name, cmd):
    try:
        r = subprocess.run(cmd, cwd=ROOT, capture_output=True, text=True, check=False)
    except FileNotFoundError as err:
        return record(name, "not_run", str(err))
    record(name, "pass" if r.returncode == 0 else "fail", (r.stdout + r.stderr).strip()[-400:])


def check_schemas():
    loaders = {
        "config/dataset.yaml": config.load_dataset,
        "config/axes.yaml": config.load_axes,
        "config/checks.yaml": config.load_checks,
        "config/hillclimb.yaml": config.load_hillclimb,
    }
    loaders |= {
        p.relative_to(ROOT).as_posix(): config.load_markers for p in ROOT.glob("markers/*.yaml")
    }
    for rel, load in sorted(loaders.items()):
        try:
            load(ROOT / rel)
            record(f"schema {rel}", "pass")
        except Exception as err:
            record(f"schema {rel}", "fail", str(err).splitlines()[0])


def check_manifests():
    manifests = sorted(ROOT.glob("results/runs/*/manifest.yaml"))
    if not manifests:
        return record("run manifests", "not_run", "no runs yet")
    for path in manifests:
        m = yaml.safe_load(path.read_text())
        rel = path.relative_to(ROOT).as_posix()
        stale = []
        for item in m["inputs"]:
            src = ROOT / item["path"]
            if not src.exists():
                record(f"inputs {rel}", "not_run", f"missing input {src}")
                break
            if provenance._sha256(src) != item["sha256"]:
                stale.append(str(src))
        else:
            record(f"inputs {rel}", "fail" if stale else "pass", ", ".join(stale))


def check_holdout_isolation():
    readers = [
        p.relative_to(ROOT).as_posix()
        for p in [*ROOT.glob("src/**/*.py"), *ROOT.glob("scripts/**/*.py")]
        if "holdout" in p.read_text() and p.name not in NAMES_HOLDOUT_WITHOUT_READING
    ]
    record("holdout read only by evaluate.py", "fail" if readers else "pass", ", ".join(readers))


def check_orthogonality_contract():
    import importlib.util

    missing = [m for m in ("zf.checks", "zf.synthetic") if importlib.util.find_spec(m) is None]
    if missing:
        return record("check orthogonality", "not_run", f"not implemented: {', '.join(missing)}")
    run(
        "check orthogonality",
        ["uv", "run", "--frozen", "pytest", "-q", "tests/test_check_orthogonality.py"],
    )


def main(allow_not_run=False):
    run("ruff check", ["uv", "run", "--frozen", "ruff", "check", "."])
    run("ruff format", ["uv", "run", "--frozen", "ruff", "format", "--check", "."])
    run("pytest", ["uv", "run", "--frozen", "pytest", "-q"])
    check_schemas()
    check_manifests()
    check_holdout_isolation()
    check_orthogonality_contract()

    for name, outcome, detail in results:
        print(
            f"{outcome.upper():8} {name}"
            + (f"  — {detail}" if detail and outcome != "pass" else "")
        )
    outcomes = {o for _, o, _ in results}
    if "fail" in outcomes:
        return 1
    if "not_run" in outcomes:
        if allow_not_run:
            print("Some checks could not run; they are not counted as passed.")
            return 0
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main(allow_not_run="--allow-not-run" in sys.argv[1:]))
