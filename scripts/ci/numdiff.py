#!/usr/bin/env python3
"""Token-wise numeric-tolerant file comparison for COJO .cma/.ldr golden tests.

Used by golden_test.sh in --tol mode: structure and non-numeric tokens
(SNP names, alleles, "NA", headers) must match exactly; numeric tokens may differ
within tolerance (BLAS/compiler reassociation perturbs low-order digits). Line and
per-line token counts must match exactly — a shape change is always a real diff.

Two tokens compare equal when both parse as floats and
    |a - b| <= atol + rtol * max(|a|, |b|)
otherwise they must be string-equal. NaN matches only NaN.

Exit 0 if equivalent, 1 otherwise. Prints a short summary + first mismatches.

Usage: numdiff.py A B [--tol ATOL] [--rtol RTOL] [--max-report N]
"""
import argparse
import math
import sys


def as_float(tok):
    try:
        return float(tok)
    except (ValueError, TypeError):
        return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("a")
    ap.add_argument("b")
    ap.add_argument("--tol", type=float, default=1e-6, help="absolute tolerance")
    ap.add_argument("--rtol", type=float, default=1e-6, help="relative tolerance")
    ap.add_argument("--max-report", type=int, default=10)
    args = ap.parse_args()

    with open(args.a) as fa, open(args.b) as fb:
        la, lb = fa.readlines(), fb.readlines()

    mismatches = []
    max_num_diff = 0.0

    if len(la) != len(lb):
        print(f"VERDICT: DIFFERENT (line count {len(la)} vs {len(lb)})")
        return 1

    for i, (ra, rb) in enumerate(zip(la, lb), 1):
        ta, tb = ra.split(), rb.split()
        if len(ta) != len(tb):
            mismatches.append(f"line {i}: token count {len(ta)} vs {len(tb)}")
            continue
        for j, (xa, xb) in enumerate(zip(ta, tb), 1):
            if xa == xb:
                continue
            fa_, fb_ = as_float(xa), as_float(xb)
            if fa_ is None or fb_ is None:
                mismatches.append(f"line {i} col {j}: {xa!r} != {xb!r}")
                continue
            if math.isnan(fa_) or math.isnan(fb_):
                if not (math.isnan(fa_) and math.isnan(fb_)):
                    mismatches.append(f"line {i} col {j}: {xa} != {xb} (NaN)")
                continue
            d = abs(fa_ - fb_)
            max_num_diff = max(max_num_diff, d)
            if d > args.tol + args.rtol * max(abs(fa_), abs(fb_)):
                mismatches.append(f"line {i} col {j}: {xa} != {xb} (|d|={d:.3g})")

    print(f"max numeric diff = {max_num_diff:.3g}; {len(mismatches)} mismatch(es) "
          f"(atol={args.tol}, rtol={args.rtol})")
    for m in mismatches[:args.max_report]:
        print("  " + m)
    if len(mismatches) > args.max_report:
        print(f"  ... and {len(mismatches) - args.max_report} more")

    if mismatches:
        print("VERDICT: DIFFERENT")
        return 1
    print("VERDICT: EQUIVALENT")
    return 0


if __name__ == "__main__":
    sys.exit(main())
