"""Deterministic rules checked before Claude Code writes a file. Run as a PreToolUse hook.

The hook computes the file as it would be after the write, then applies the rules below.
A rule that blocks says why, so the agent can correct the write instead of guessing.
"""

import copy
from dataclasses import dataclass
from fnmatch import fnmatch
from pathlib import Path

import yaml
from pydantic import ValidationError

from zf import config

WRITE_TOOLS = {"Write", "Edit", "MultiEdit"}
HUMAN_CONFIG = {
    "config/dataset.yaml": config.Dataset,
    "config/axes.yaml": config.Axes,
    "config/checks.yaml": config.Checks,
}


@dataclass
class Decision:
    allowed: bool
    reason: str = ""


ALLOW = Decision(True)


def block(reason: str) -> Decision:
    return Decision(False, reason)


def _apply_edits(text: str, edits: list[dict]) -> str | None:
    for e in edits:
        old, new = e["old_string"], e["new_string"]
        if old not in text:
            return None
        text = text.replace(old, new) if e.get("replace_all") else text.replace(old, new, 1)
    return text


def _new_text(tool: str, tool_input: dict, current: str) -> str | None:
    if tool == "Write":
        return tool_input["content"]
    if tool == "Edit":
        return _apply_edits(current, [tool_input])
    return _apply_edits(current, tool_input["edits"])


def _parse(text: str, model) -> tuple[dict | None, str]:
    try:
        data = yaml.safe_load(text)
        model.model_validate(data)
        return data, ""
    except (yaml.YAMLError, ValidationError) as err:
        return None, str(err).splitlines()[0] + " ..."


def _without_flags(markers: dict) -> dict:
    stripped = copy.deepcopy(markers)
    for gene in stripped.get("genes", []):
        gene.pop("flags", None)
    return stripped


def _flags_append_only(old: dict, new: dict) -> bool:
    for before, after in zip(old["genes"], new["genes"], strict=True):
        if after["flags"][: len(before["flags"])] != before["flags"]:
            return False
    return True


def decide(tool, tool_input, root, env, current_text=None) -> Decision:
    if tool not in WRITE_TOOLS:
        return ALLOW
    root = Path(root).resolve()
    path = Path(tool_input["file_path"]).resolve()
    if not path.is_relative_to(root):
        return ALLOW
    rel = path.relative_to(root).as_posix()
    human = env.get("ZF_HUMAN_EDIT") == "1"

    if rel.startswith("holdout/"):
        return block("No agent may write to holdout/.")
    if fnmatch(rel, "results/runs/*/manifest.yaml"):
        return block("Run manifests are written only by zf.provenance.write_manifest.")

    is_config = rel in HUMAN_CONFIG
    is_markers = fnmatch(rel, "markers/*.yaml")
    if not (is_config or is_markers):
        return ALLOW

    if current_text is None:
        current_text = path.read_text() if path.exists() else ""
    new_text = _new_text(tool, tool_input, current_text)
    if new_text is None:
        return block("The edit does not apply to the current file; read it again.")

    if is_config:
        old = yaml.safe_load(current_text) or {}
        if old.get("status") == "frozen":
            return block(f"{rel} is frozen. Unfreeze it by hand, in a commit that says why.")
        if not human:
            return block(
                f"{rel} is written by humans. Propose the change in JOURNAL.md, "
                "or set ZF_HUMAN_EDIT=1 if a human asked for it."
            )
        _, err = _parse(new_text, HUMAN_CONFIG[rel])
        return block(f"{rel} would not match its schema: {err}") if err else ALLOW

    new, err = _parse(new_text, config.Markers)
    if err:
        return block(f"{rel} would not match the markers schema: {err}")
    if human or not current_text:
        return ALLOW
    old = yaml.safe_load(current_text)
    if _without_flags(old) != _without_flags(new):
        return block(f"The loop may only append to `flags` in {rel}; other fields are human-owned.")
    if not _flags_append_only(old, new):
        return block(f"Flags in {rel} are append-only; do not remove or reorder them.")
    return ALLOW
