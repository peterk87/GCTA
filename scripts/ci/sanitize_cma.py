#!/usr/bin/env python3
"""Sanitize GCTA .cma.cojo conditional columns before golden compare.

Some COJO paths leave bC/bC_se/pC uninitialized for SNPs that were not
conditioned (should be NA). Those print as denormals / absurd magnitudes and
break byte-exact goldens across runs. Replace non-finite, denormal, or
absurdly large values in those three columns with NA. Other columns untouched.

Usage: sanitize_cma.py <infile >outfile
"""
import math
import sys

# Columns (0-based) after header: bC, bC_se, pC
COND_COLS = {10, 11, 12}
ABSURD = 1e100


def is_garbage(tok: str) -> bool:
    if tok == "NA":
        return False
    try:
        x = float(tok)
    except ValueError:
        return False
    if not math.isfinite(x):
        return True
    if abs(x) > ABSURD:
        return True
    # denormal junk (real printed p-values are normal floats or NA)
    if x != 0.0 and abs(x) < sys.float_info.min:
        return True
    return False


def main() -> int:
    lines = sys.stdin.readlines()
    if not lines:
        return 0
    sys.stdout.write(lines[0])
    for line in lines[1:]:
        parts = line.rstrip("\n").split("\t")
        for i in COND_COLS:
            if i < len(parts) and is_garbage(parts[i]):
                parts[i] = "NA"
        sys.stdout.write("\t".join(parts) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
