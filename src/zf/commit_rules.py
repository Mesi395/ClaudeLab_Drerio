"""What may never be committed or pushed. Used by the pre-commit and pre-push hooks.

.gitignore stops accidental adds; these rules also stop `git add -f`.
"""

from fnmatch import fnmatch
from pathlib import PurePosixPath

PROTECTED_DIRS = ("holdout/", "data/", "notes/", "airdrop/")
PROTECTED_FILES = ("CLAUDE.local.md", ".claude/settings.local.json")
PROTECTED_PATTERNS = (
    "*.h5ad",
    "*.h5",
    "*.loom",
    "*.rds",
    "*.mtx",
    "*.mtx.gz",
    "*.pkl",
    "*.tar",
    "*.tar.gz",
    "*.bundle",
    ".env",
    ".env.*",
    "*.pem",
    "*.key",
)
PROTECTED_BRANCHES = ("refs/heads/main",)


def reason(path: str) -> str | None:
    """Why `path` may not be committed, or None if it may."""
    p = PurePosixPath(path).as_posix()
    if p.startswith(PROTECTED_DIRS):
        return f"{p}: {p.split('/')[0]}/ is never committed"
    if p in PROTECTED_FILES:
        return f"{p}: personal file, never committed"
    name = PurePosixPath(p).name
    for pattern in PROTECTED_PATTERNS:
        if fnmatch(name, pattern):
            return f"{p}: {pattern} files are never committed"
    return None


def push_reason(remote_ref: str) -> str | None:
    """Why a push to `remote_ref` is refused, or None if it is allowed."""
    if remote_ref in PROTECTED_BRANCHES:
        branch = remote_ref.removeprefix("refs/heads/")
        return f"Direct pushes to {branch} are not allowed; open a pull request."
    return None
