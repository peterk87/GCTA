#!/usr/bin/env bash
# Package a mostly-static gcta64 for GitHub Releases.
#
# Produces:
#   gcta64-<ver>-linux-x86_64           stripped + .gnu_debuglink (default runtime)
#   gcta64-<ver>-linux-x86_64.debug     detached DWARF (use with gdb / eu-unstrip)
#   gcta64-<ver>-linux-x86_64.upx       UPX --lzma of the stripped binary (smallest)
#   SHA256SUMS
#   README.txt
#
# Debugging stays easy: keep the .debug file next to the stripped binary (or point
# gdb at it). Do not use the .upx asset for debugging — UPX compression breaks
# symbol loading until decompressed (`upx -d`).
set -euo pipefail

usage() {
  echo "Usage: $0 --bin PATH --outdir DIR --version X.Y.Z [--commit SHA]" >&2
  exit 2
}

BIN=""
OUTDIR=""
VERSION=""
COMMIT=""

while [ $# -gt 0 ]; do
  case "$1" in
    --bin) BIN="$2"; shift 2;;
    --outdir) OUTDIR="$2"; shift 2;;
    --version) VERSION="$2"; shift 2;;
    --commit) COMMIT="$2"; shift 2;;
    -h|--help) usage;;
    *) echo "Unknown arg: $1" >&2; usage;;
  esac
done

[ -n "$BIN" ] && [ -n "$OUTDIR" ] && [ -n "$VERSION" ] || usage
[ -x "$BIN" ] || { echo "Not an executable: $BIN" >&2; exit 1; }
command -v objcopy >/dev/null || { echo "objcopy required (binutils)" >&2; exit 1; }
command -v upx >/dev/null || { echo "upx required" >&2; exit 1; }

mkdir -p "$OUTDIR"
NAME="gcta64-${VERSION}-linux-x86_64"
STRIPPED="$OUTDIR/$NAME"
DEBUG="$OUTDIR/${NAME}.debug"
UPX_BIN="$OUTDIR/${NAME}.upx"

cp -a "$BIN" "$STRIPPED"
# Detach DWARF into a sibling .debug file, strip the runtime binary, then add
# .gnu_debuglink. The debuglink filename must resolve beside the binary, so we
# pass a basename and run objcopy with OUTDIR as cwd.
objcopy --only-keep-debug "$STRIPPED" "$DEBUG"
objcopy --strip-unneeded "$STRIPPED"
(
  cd "$OUTDIR"
  objcopy --add-gnu-debuglink="${NAME}.debug" "$NAME"
)
chmod +x "$STRIPPED"
chmod a-x "$DEBUG" || true

cp -a "$STRIPPED" "$UPX_BIN"
# Prefer lzma for size; --best is slower but fine for release packaging.
upx --best --lzma -q "$UPX_BIN"

{
  echo "# gcta64 ${VERSION} (Linux x86_64, mostly-static)"
  echo "#"
  echo "# Built with -DGCTA_STATIC_EXE=ON (static MKL + GSL + libgcc/libstdc++;"
  echo "# glibc and libgomp remain dynamic)."
  echo "#"
  echo "# Assets:"
  echo "#   ${NAME}         stripped runtime binary with .gnu_debuglink"
  echo "#   ${NAME}.debug   detached debug symbols (place beside the binary for gdb)"
  echo "#   ${NAME}.upx     UPX/LZMA compressed copy of the stripped binary"
  echo "#"
  echo "# Prefer ${NAME} (+ optional .debug) when debugging. Prefer ${NAME}.upx when"
  echo "# you only need a smaller download; decompress with: upx -d ${NAME}.upx"
  echo "#"
  if [ -n "$COMMIT" ]; then
    echo "# Commit: ${COMMIT}"
  fi
  echo "# Size (bytes):"
  echo "#   stripped: $(stat -c%s "$STRIPPED")"
  echo "#   debug:    $(stat -c%s "$DEBUG")"
  echo "#   upx:      $(stat -c%s "$UPX_BIN")"
} >"$OUTDIR/README.txt"

(
  cd "$OUTDIR"
  sha256sum "$NAME" "${NAME}.debug" "${NAME}.upx" >SHA256SUMS
)

echo "Packaged release binaries in $OUTDIR:"
ls -lh "$STRIPPED" "$DEBUG" "$UPX_BIN"
cat "$OUTDIR/SHA256SUMS"
