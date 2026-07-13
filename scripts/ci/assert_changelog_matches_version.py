#!/usr/bin/env python3
"""Fail unless CHANGELOG.md has a non-empty ## [version] matching src/config.h."""

from __future__ import annotations

import sys
from pathlib import Path

_REPO_ROOT = Path(__file__).resolve().parents[2]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from scripts.ci.changelog_release_section import section_body_for_version  # noqa: E402
from scripts.ci.gcta_version import read_gcta_version  # noqa: E402


def main() -> None:
    version = read_gcta_version()
    changelog = _REPO_ROOT / "CHANGELOG.md"
    body = section_body_for_version(version, changelog)
    if body is None:
        print(
            f"::error::CHANGELOG.md must have a non-empty section "
            f"## [{version}] matching GCTA_VERSION in src/config.h",
            file=sys.stderr,
        )
        sys.exit(1)
    print(f"OK: CHANGELOG.md has a section for GCTA_VERSION {version}.")


if __name__ == "__main__":
    main()
