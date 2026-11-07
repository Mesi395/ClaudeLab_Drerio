"""pre-commit hook: refuse staged files that must never be committed."""

import sys

from zf.commit_rules import reason

problems = [r for r in map(reason, sys.argv[1:]) if r]
for problem in problems:
    print(problem, file=sys.stderr)
sys.exit(1 if problems else 0)
