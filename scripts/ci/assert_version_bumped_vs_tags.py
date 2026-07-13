#!/usr/bin/env python3
"""Fail unless src/config.h GCTA_VERSION is strictly greater than every X.Y.Z / vX.Y.Z tag."""

from __future__ import annotations

import re
import subprocess
import sys
from pathlib import Path

_REPO_ROOT = Path(__file__).resolve().parents[2]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from scripts.ci.gcta_version import read_gcta_version  # noqa: E402

_TAG_RE = re.compile(r"^v?(\d+)\.(\d+)\.(\d+)$")
_VERSION_RE = re.compile(r"^(\d+)\.(\d+)\.(\d+)$")


def _parse_triplet(v: str) -> tuple[int, int, int]:
    m = _VERSION_RE.match(v.strip())
    if not m:
        return (0, 0, 0)
    return (int(m.group(1)), int(m.group(2)), int(m.group(3)))


def _iter_semver_tags() -> list[tuple[int, int, int]]:
    raw = subprocess.run(
        ["git", "tag", "-l"],
        check=True,
        capture_output=True,
        text=True,
    ).stdout
    out: list[tuple[int, int, int]] = []
    for line in raw.splitlines():
        t = line.strip()
        m = _TAG_RE.match(t)
        if m:
            out.append((int(m.group(1)), int(m.group(2)), int(m.group(3))))
    return out


def main() -> None:
    about_raw = read_gcta_version()
    about = _parse_triplet(about_raw)
    if about == (0, 0, 0):
        print(
            f"GCTA_VERSION {about_raw!r} in src/config.h must match X.Y.Z (semver).",
            file=sys.stderr,
        )
        sys.exit(1)

    tags = _iter_semver_tags()
    latest = max(tags) if tags else (0, 0, 0)

    if about <= latest:
        print(
            (
                f"src/config.h GCTA_VERSION is {about_raw}, but the latest matching git tag is "
                f"{latest[0]}.{latest[1]}.{latest[2]}. Bump GCTA_VERSION to a value strictly "
                "greater than every existing X.Y.Z / vX.Y.Z tag before merging to main."
            ),
            file=sys.stderr,
        )
        sys.exit(1)

    if tags:
        print(
            f"OK: GCTA_VERSION {about_raw} is greater than latest tag "
            f"{latest[0]}.{latest[1]}.{latest[2]}."
        )
    else:
        print(f"OK: GCTA_VERSION {about_raw} (no X.Y.Z / vX.Y.Z tags in this repo yet).")


if __name__ == "__main__":
    main()
