#!/usr/bin/env python3
"""Read GCTA_VERSION from src/config.h (prints X.Y.Z without a leading v)."""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

_REPO_ROOT = Path(__file__).resolve().parents[2]
_CONFIG_H = _REPO_ROOT / "src" / "config.h"
_VERSION_RE = re.compile(
    r'^\s*#\s*define\s+GCTA_VERSION\s+"v?(\d+\.\d+\.\d+)"\s*$',
    re.MULTILINE,
)


def read_gcta_version(config_path: Path = _CONFIG_H) -> str:
    if not config_path.is_file():
        raise FileNotFoundError(f"Missing GCTA config header: {config_path}")
    text = config_path.read_text(encoding="utf-8")
    match = _VERSION_RE.search(text)
    if not match:
        raise ValueError(
            f"Could not parse GCTA_VERSION from {config_path} "
            '(expected #define GCTA_VERSION "vX.Y.Z")'
        )
    return match.group(1)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=(__doc__ or "").strip())
    _ = parser.add_argument(
        "--config",
        type=Path,
        default=_CONFIG_H,
        help="Path to src/config.h (default: repo src/config.h)",
    )
    _ = parser.add_argument(
        "--tag",
        action="store_true",
        help="Print vX.Y.Z instead of X.Y.Z",
    )
    return parser


def main() -> None:
    args = build_parser().parse_args()
    try:
        version = read_gcta_version(args.config)
    except (OSError, ValueError) as exc:
        print(f"::error::{exc}", file=sys.stderr)
        sys.exit(1)
    print(f"v{version}" if args.tag else version)


if __name__ == "__main__":
    main()
