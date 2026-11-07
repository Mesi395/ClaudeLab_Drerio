"""PreToolUse hook: reads the tool call as JSON on stdin; exit 2 blocks it and shows the reason."""

import json
import os
import sys
from pathlib import Path

root = Path(os.environ.get("CLAUDE_PROJECT_DIR", Path(__file__).resolve().parents[2]))

try:
    from zf.guard import decide

    call = json.load(sys.stdin)
    decision = decide(call.get("tool_name", ""), call.get("tool_input", {}), root, os.environ)
except Exception as err:  # the guard fails closed: a write it cannot check is not made
    print(
        f"Write guard could not run ({type(err).__name__}: {err}); write blocked.", file=sys.stderr
    )
    sys.exit(2)

if not decision.allowed:
    print(decision.reason, file=sys.stderr)
    sys.exit(2)
