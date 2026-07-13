#!/usr/bin/env python3
"""Compare two GCTA-COJO .jma.cojo outputs for equivalence.

Exit 0 if SNP sets match and max |bJ| diff <= --bj-tol; else exit 1.
"""
import argparse
import csv
import math
import sys


def load(path):
    rows = {}
    with open(path) as f:
        r = csv.DictReader(f, delimiter="\t")
        for row in r:
            rows[row["SNP"]] = row
    return rows


def fnum(x):
    try:
        return float(x)
    except (ValueError, TypeError):
        return float("nan")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("a")
    ap.add_argument("b")
    ap.add_argument("--bj-tol", type=float, default=1e-3)
    args = ap.parse_args()

    A, B = load(args.a), load(args.b)
    sa, sb = set(A), set(B)
    shared = sa & sb

    print(f"A={args.a}: {len(sa)} signals")
    print(f"B={args.b}: {len(sb)} signals")
    print(f"shared: {len(shared)}   only-A: {len(sa - sb)}   only-B: {len(sb - sa)}")
    if sa - sb:
        print("  only in A:", sorted(sa - sb)[:20])
    if sb - sa:
        print("  only in B:", sorted(sb - sa)[:20])

    if sa != sb:
        print("VERDICT: DIFFERENT (SNP set mismatch)")
        return 1

    if not shared:
        print("VERDICT: EQUIVALENT (empty set)")
        return 0

    dbj = []
    for s in shared:
        ba, bb = fnum(A[s].get("bJ")), fnum(B[s].get("bJ"))
        dbj.append(abs(ba - bb))
    dbj.sort()
    maxd = dbj[-1]
    med = dbj[len(dbj) // 2]
    n_over = sum(1 for d in dbj if d > args.bj_tol or math.isnan(d))
    print(f"bJ abs-diff: median={med:.3g} max={maxd:.3g} "
          f"({n_over}/{len(shared)} exceed tol={args.bj_tol})")
    if n_over or math.isnan(maxd):
        print("VERDICT: DIFFERENT")
        return 1
    print("VERDICT: EQUIVALENT")
    return 0


if __name__ == "__main__":
    sys.exit(main())
