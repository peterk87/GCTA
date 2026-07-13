"""Parse the Markdown body under ``## [<version>]`` in CHANGELOG-style files."""

from __future__ import annotations

import re
from pathlib import Path


def _section_pattern(version: str) -> re.Pattern[str]:
    return re.compile(
        rf"^## \[{re.escape(version)}\][^\n]*\n(.*?)(?=\n## \[|\Z)",
        re.MULTILINE | re.DOTALL,
    )


def section_body_for_version(version: str, changelog_path: Path) -> str | None:
    """Return body text under ``## [version]`` or ``None`` if missing or whitespace-only."""
    if not changelog_path.is_file():
        return None
    text = changelog_path.read_text(encoding="utf-8")
    match = _section_pattern(version).search(text)
    if not match:
        return None
    body = match.group(1).strip()
    return body if body else None


def heading_exists_for_version(version: str, text: str) -> bool:
    return bool(re.search(rf"^## \[{re.escape(version)}\]", text, re.MULTILINE))
