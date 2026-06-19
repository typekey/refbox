#!/usr/bin/env bash
# build.sh — one-shot driver around `refbox pull`:
#   download (if missing) → build → test → publish (flatten to <Assembly>.<name>)
#
# Usage:
#   ./build.sh                              # all enabled assemblies
#   ./build.sh Homo_sapiens                 # one species (all enabled assemblies)
#   ./build.sh Homo_sapiens GRCh38          # one species + assembly
#   ./build.sh -- --resource genome cytoband   # extra args after `--` go to refbox
#
# Output layout (published / flat — the default):
#   $REFBOX_OUT/<Species>/<Assembly>/<Assembly>.<name>   (no build/ or raw/)
#
# Environment:
#   REFBOX_OUT         output root (default: this directory's parent)
#   REFBOX_CONFIG      path to species.yaml (default: ./config/species.yaml)
#   FORCE=1            pass --force (re-download + rebuild)
#   NO_TEST=1          skip the validation step (--no-test)
#   NO_FLAT=1          keep the build/ + raw/ working tree (--no-flat)
#   INCLUDE_DISABLED=1 also process assemblies marked enabled: false

set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# Default output root = parent of the repo so reference data lives outside the
# git checkout. Override with REFBOX_OUT=/some/path or --out on the CLI.
export REFBOX_OUT="${REFBOX_OUT:-$(cd "$HERE/.." && pwd)}"

# ── parse args ───────────────────────────────────────────────────────────────
SPECIES=()
ASSEMBLY=()
EXTRA=()

# split "before --" / "after --"
seen_sep=0
for arg in "$@"; do
    if [[ "$arg" == "--" ]]; then seen_sep=1; continue; fi
    if (( seen_sep )); then
        EXTRA+=("$arg")
    elif (( ${#SPECIES[@]} == 0 )); then
        SPECIES+=("$arg")
    elif (( ${#ASSEMBLY[@]} == 0 )); then
        ASSEMBLY+=("$arg")
    else
        EXTRA+=("$arg")
    fi
done

FILTER=()
(( ${#SPECIES[@]}  )) && FILTER+=(--species  "${SPECIES[@]}")
(( ${#ASSEMBLY[@]} )) && FILTER+=(--assembly "${ASSEMBLY[@]}")

FLAGS=()
[[ "${FORCE:-}" == "1" ]]            && FLAGS+=(--force)
[[ "${NO_TEST:-}" == "1" ]]          && FLAGS+=(--no-test)
[[ "${NO_FLAT:-}" == "1" ]]          && FLAGS+=(--no-flat)
[[ "${INCLUDE_DISABLED:-}" == "1" ]] && FLAGS+=(--include-disabled)

# ── sanity: refbox installed? ────────────────────────────────────────────────
if ! command -v refbox >/dev/null; then
    echo ">>> refbox not on PATH; installing from $HERE in editable mode..."
    pip install -e "$HERE"
fi

# ── sanity: external tools ───────────────────────────────────────────────────
for tool in bgzip tabix samtools sort grep; do
    if ! command -v "$tool" >/dev/null; then
        echo "ERROR: required tool '$tool' not found in PATH" >&2
        exit 1
    fi
done
# cytoband bigBed conversion needs the UCSC kent tools; warn (don't fail) so
# runs that don't build cytoband still work.
for tool in bedToBigBed bigBedToBed; do
    command -v "$tool" >/dev/null || \
        echo "WARNING: '$tool' not found; the cytoband resource will be skipped" >&2
done

echo "============================================================"
echo " refbox driver  (download → build → test → publish)"
echo "   out      : $REFBOX_OUT"
echo "   species  : ${SPECIES[*]:-<all enabled>}"
echo "   assembly : ${ASSEMBLY[*]:-<all>}"
echo "   flags    : ${FLAGS[*]:-<none>}"
echo "   extra    : ${EXTRA[*]:-<none>}"
echo "============================================================"

refbox pull "${FILTER[@]}" "${FLAGS[@]}" "${EXTRA[@]}"

echo ">>> done."
