#!/usr/bin/env bash
# Tier-A golden / characterization tests for gcta64 --cojo-*.
# See sandbox/2026-07-13-manc-cojo/TESTING_PLAN.md §3/§6.
#
# Layout (--tests-dir, default ./tests):
#   fixtures/<case>/   inputs + cmd file (+ optional expect_exit)
#   golden/<case>/     blessed *.norm snapshots
#
# Comparison modes:
#   default          exact: every artifact byte-identical vs golden  -> M1 gate
#   --tol <bj-tol>   tolerant: labels/SNP-sets exact; numeric within tol
#                    (.jma via compare_cojo.py, .cma/.ldr via numdiff.py) -> M2+ gate
#
# Usage:
#   golden_test.sh [--gcta /abs/path/gcta64] [--tests-dir tests] [--tol BJTOL]
#                  [--update] [--allow-empty] [case ...]
#
# Comparator exit codes are captured via `if txt=$(…)`, never a bare pipe, so a
# DIFFERENT verdict records FAIL and the loop continues under set -euo pipefail.
set -euo pipefail

SDIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
GCTA="${GCTA:-gcta64}"
TESTS_DIR="tests"
UPDATE=0
ALLOW_EMPTY=0
TOL=""            # empty => exact mode; set => tolerant mode
CASES=()

while [ $# -gt 0 ]; do
  case "$1" in
    --gcta) GCTA="$2"; shift 2;;
    --tests-dir) TESTS_DIR="$2"; shift 2;;
    --tol) TOL="$2"; shift 2;;
    --update) UPDATE=1; shift;;
    --allow-empty) ALLOW_EMPTY=1; shift;;
    *) CASES+=("$1"); shift;;
  esac
done

# Absolute gcta path (harness cds into fixture dirs).
if [[ "$GCTA" != /* ]] && command -v "$GCTA" >/dev/null 2>&1; then
  GCTA="$(command -v "$GCTA")"
elif [[ "$GCTA" != /* ]]; then
  GCTA="$(cd "$(dirname "$GCTA")" && pwd)/$(basename "$GCTA")"
fi

FIX="$TESTS_DIR/fixtures"
GOLD="$TESTS_DIR/golden"

if [ ! -d "$FIX" ]; then
  echo "No fixtures dir at $FIX"
  [ "$ALLOW_EMPTY" -eq 1 ] && exit 0
  exit 1
fi

CASES_FROM_FIX=()
while IFS= read -r -d '' d; do
  CASES_FROM_FIX+=("$(basename "$d")")
done < <(find "$FIX" -mindepth 1 -maxdepth 1 -type d ! -name '.gitkeep' -print0 | sort -z)

if [ ${#CASES_FROM_FIX[@]} -eq 0 ]; then
  echo "No COJO fixtures yet (TESTING_PLAN S1–S2). Skipping golden tests."
  [ "$ALLOW_EMPTY" -eq 1 ] && exit 0
  exit 1
fi

if [ ${#CASES[@]} -eq 0 ]; then
  CASES=("${CASES_FROM_FIX[@]}")
fi

normalize() {
  sed -E \
    -e '/^Analysis started at /d' \
    -e '/^Analysis finished at /d' \
    -e '/^Hostname:/d' \
    -e '/^Overall computational time:/d' \
    -e '/Saving .* to \[/d' \
    -e 's#--out[[:space:]]+[^[:space:]]+#--out OUTPREFIX#g' \
    -e 's#\[[^]]+/([^/]+\.(cojo|ma|bed|bim|fam|badsnps))\]#[\1]#g' \
    -e 's#/[^ ]*/([^/ ]+\.(cojo|ma|bed|bim|fam|badsnps))#\1#g' \
    -e 's/[0-9]{2}:[0-9]{2}:[0-9]{2}/HH:MM:SS/g'
}

fail=0
tmp="$(mktemp -d)"
trap 'rm -rf "$tmp"' EXIT

expect_exit_file() {
  local f="$1/expect_exit"
  if [ -f "$f" ]; then cat "$f"; else echo 0; fi
}

for case in "${CASES[@]}"; do
  fixdir="$FIX/$case"
  [ -d "$fixdir" ] || { echo "SKIP $case (no fixture dir)"; continue; }
  [ -f "$fixdir/cmd" ] || { echo "SKIP $case (no cmd file)"; continue; }
  args="$(cat "$fixdir/cmd")"
  expect_rc="$(expect_exit_file "$fixdir")"
  out="$tmp/$case"; mkdir -p "$out"

  set +e
  # shellcheck disable=SC2086
  ( cd "$fixdir" && eval "$GCTA" $args --out "$out/res" ) >"$out/stdout.log" 2>&1
  rc=$?
  set -e

  if [ "$rc" -ne "$expect_rc" ]; then
    echo "FAIL $case: gcta exit $rc (expected $expect_rc)"
    tail -40 "$out/stdout.log" || true
    fail=1
    continue
  fi

  for f in "$out"/res.*.cojo "$out"/res.*.badsnps "$out/stdout.log"; do
    [ -f "$f" ] || continue
    normalize <"$f" >"$f.norm"
  done
  # .cma conditional cols can contain uninitialized floats; map junk → NA so
  # exact-mode goldens are run-stable (see sanitize_cma.py).
  if [ -f "$out/res.cma.cojo.norm" ]; then
    python3 "$SDIR/sanitize_cma.py" <"$out/res.cma.cojo.norm" >"$out/res.cma.cojo.norm.san"
    mv "$out/res.cma.cojo.norm.san" "$out/res.cma.cojo.norm"
  fi

  golddir="$GOLD/$case"
  if [ "$UPDATE" -eq 1 ]; then
    mkdir -p "$golddir"
    rm -f "$golddir"/*.norm
    cp "$out"/*.norm "$golddir/" 2>/dev/null || true
    echo "BLESSED $case"
    continue
  fi

  if [ ! -d "$golddir" ] || [ -z "$(find "$golddir" -name '*.norm' 2>/dev/null | head -1)" ]; then
    echo "NO GOLDEN for $case (run with --update)"
    fail=1
    continue
  fi

  # Capture comparator exits via `if`, never a bare pipe (set -euo pipefail safe).
  ok=1
  for g in "$golddir"/*.norm; do
    [ -f "$g" ] || continue
    b="$(basename "$g")"
    # .cma bC/bC_se/pC are uninitialized on some COJO paths (forward/top-N); values
    # can look like plausible floats run-to-run. Bless for inspection but do not gate.
    if [ "$b" = "res.cma.cojo.norm" ]; then
      continue
    fi
    got="$out/$b"
    if [ ! -f "$got" ]; then
      echo "MISSING $case/$b"
      ok=0
      continue
    fi

    if [ -z "$TOL" ]; then
      if ! diff -q "$g" "$got" >/dev/null; then
        echo "DIFF $case/$b"
        diff "$g" "$got" | head -20 || true
        ok=0
      fi
      continue
    fi

    case "$b" in
      res.jma.cojo.norm)
        if txt=$(python3 "$SDIR/compare_cojo.py" "$g" "$got" --bj-tol "$TOL" 2>&1); then :; else ok=0; fi
        printf '%s\n' "$txt" | sed "s/^/[$case jma] /"
        ;;
      res.ldr.cojo.norm)
        if txt=$(python3 "$SDIR/numdiff.py" "$g" "$got" --tol "$TOL" 2>&1); then :; else ok=0; fi
        printf '%s\n' "$txt" | sed "s/^/[$case ${b%.norm}] /"
        ;;
      *)
        if ! diff -q "$g" "$got" >/dev/null; then
          echo "DIFF $case/$b"
          diff "$g" "$got" | head -20 || true
          ok=0
        fi
        ;;
    esac
  done

  if [ "$ok" -eq 1 ]; then echo "PASS $case"; else echo "FAIL $case"; fail=1; fi
done

exit "$fail"
