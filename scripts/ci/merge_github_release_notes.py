#!/usr/bin/env python3
"""Merge CHANGELOG.md section for a version with a GitHub generate-notes Markdown body."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

_REPO_ROOT = Path(__file__).resolve().parents[2]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from scripts.ci.changelog_release_section import (  # noqa: E402
    heading_exists_for_version,
    section_body_for_version,
)


def changelog_section(version: str, changelog_path: Path) -> str | None:
    if not changelog_path.is_file():
        print("::notice::CHANGELOG.md not present", file=sys.stderr)
        return None
    body = section_body_for_version(version, changelog_path)
    if body is None:
        text = changelog_path.read_text(encoding="utf-8")
        if heading_exists_for_version(version, text):
            print(
                f"::notice::CHANGELOG.md has ## [{version}] but no body under that heading",
                file=sys.stderr,
            )
        else:
            print(
                f"::notice::No CHANGELOG.md section found for version {version}",
                file=sys.stderr,
            )
        return None
    return "## Changelog\n\n" + body


def merge_notes(
    version: str,
    generated_path: Path,
    out_path: Path,
    changelog_path: Path,
    *,
    require_changelog_section: bool = False,
) -> None:
    if require_changelog_section:
        if not changelog_path.is_file():
            print(
                f"::error::CHANGELOG required but not found: {changelog_path}",
                file=sys.stderr,
            )
            sys.exit(1)
        body = section_body_for_version(version, changelog_path)
        if body is None:
            text = changelog_path.read_text(encoding="utf-8")
            if heading_exists_for_version(version, text):
                print(
                    f"::error::CHANGELOG must have a non-empty body under ## [{version}]",
                    file=sys.stderr,
                )
            else:
                err = (
                    f"::error::CHANGELOG must include heading ## [{version}] "
                    f"with a non-empty body for this release"
                )
                print(err, file=sys.stderr)
            sys.exit(1)
        block: str | None = "## Changelog\n\n" + body
    else:
        block = changelog_section(version, changelog_path)

    generated = generated_path.read_text(encoding="utf-8").strip()
    parts: list[str] = []
    if block is not None:
        parts.append(block)
    parts.append(generated)
    _ = out_path.write_text("\n\n---\n\n".join(parts) + "\n", encoding="utf-8")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=(__doc__ or "").strip())
    _ = parser.add_argument(
        "--release-version",
        required=True,
        metavar="SEMVER",
        help=(
            "Semantic version without leading v (e.g. 1.2.0); must match "
            "the CHANGELOG heading ## [<this value>]"
        ),
    )
    _ = parser.add_argument(
        "--generated",
        required=True,
        type=Path,
        metavar="FILE",
        help="GitHub generate-notes Markdown body written to this path",
    )
    _ = parser.add_argument(
        "--output",
        "-o",
        required=True,
        type=Path,
        metavar="FILE",
        help="Write merged release notes here",
    )
    _ = parser.add_argument(
        "--changelog",
        type=Path,
        default=_REPO_ROOT / "CHANGELOG.md",
        help="Changelog path (default: CHANGELOG.md at repository root)",
    )
    _ = parser.add_argument(
        "--require-changelog-section",
        action="store_true",
        help="Fail unless CHANGELOG has a non-empty ## [release-version] section",
    )
    return parser


def main() -> None:
    args = build_parser().parse_args()
    merge_notes(
        args.release_version,
        args.generated,
        args.output,
        args.changelog,
        require_changelog_section=args.require_changelog_section,
    )


if __name__ == "__main__":
    main()
