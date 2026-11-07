"""Run manifests: what produced a result, from which inputs, and what every check concluded."""

import hashlib
import subprocess
from datetime import datetime
from importlib.metadata import version
from pathlib import Path

import yaml

OUTCOMES = {"pass", "flag", "not_run"}


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _git() -> dict:
    def run(*args):
        r = subprocess.run(["git", *args], capture_output=True, text=True, check=False)
        return r.stdout.strip() if r.returncode == 0 else None

    return {"commit": run("rev-parse", "HEAD"), "dirty": bool(run("status", "--porcelain"))}


def write_manifest(run_dir, *, command, inputs, params, outcomes) -> Path:
    """Write results/runs/<id>/manifest.yaml. Outcomes map axis -> check id -> pass|flag|not_run."""
    for axis, checks in outcomes.items():
        bad = {o for o in checks.values()} - OUTCOMES
        if bad:
            raise ValueError(f"{axis}: unknown outcome {sorted(bad)}; use {sorted(OUTCOMES)}")
    run_dir = Path(run_dir)
    run_dir.mkdir(parents=True, exist_ok=True)
    manifest = {
        "schema": "zf.run/v1",
        "run_id": run_dir.name,
        "started": datetime.now().astimezone().isoformat(timespec="seconds"),
        "command": command,
        "git": _git(),
        "versions": {p: version(p) for p in ("zf", "scanpy", "anndata", "numpy")},
        "inputs": [{"path": str(p), "sha256": _sha256(Path(p))} for p in inputs],
        "params": params,
        "outcomes": outcomes,
    }
    path = run_dir / "manifest.yaml"
    path.write_text(yaml.safe_dump(manifest, sort_keys=False))
    return path
