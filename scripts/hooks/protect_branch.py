"""pre-push hook: refuse direct pushes to protected branches."""

import os
import sys

from zf.commit_rules import push_reason

problem = push_reason(os.environ.get("PRE_COMMIT_REMOTE_BRANCH", ""))
if problem:
    print(problem, file=sys.stderr)
sys.exit(1 if problem else 0)
