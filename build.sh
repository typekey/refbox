#!/usr/bin/env bash
# build.sh — one-shot driver: download → build → test reference files
#
# Usage:
#   ./build.sh                              # all enabled assemblies
#   ./build.sh Homo_sapiens                 # one species (all enabled assemblies)
#   ./build.sh Homo_sapiens GRCh38          # one species + assembly
#   ./build.sh -- --resource genome ccre    # extra args after `--` go to refbox
#
# Environment:
#   REFBOX_OUT       output root (default: this directory)
#   REFBOX_CONFIG    path to species.yaml (default: bundled / ./config/species.yaml)
#   STEPS            subset of steps to run, e.g. STEPS="download build" (default: "download build test")
#   FORCE=1          pass --force to download and build

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

FORCE_FLAG=()
[[ "${FORCE:-}" == "1" ]] && FORCE_FLAG=(--force)

STEPS="${STEPS:-download build test}"

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

echo "============================================================"
echo " refbox driver"
echo "   out      : $REFBOX_OUT"
echo "   species  : ${SPECIES[*]:-<all>}"
echo "   assembly : ${ASSEMBLY[*]:-<all>}"
echo "   steps    : $STEPS"
echo "   extra    : ${EXTRA[*]:-<none>}"
echo "============================================================"

for step in $STEPS; do
    case "$step" in
        download)
            echo ">>> [1/3] refbox download"
            refbox download "${FILTER[@]}" "${FORCE_FLAG[@]}" "${EXTRA[@]}"
            ;;
        build)
            echo ">>> [2/3] refbox build"
            refbox build    "${FILTER[@]}" "${FORCE_FLAG[@]}" "${EXTRA[@]}"
            ;;
        test)
            echo ">>> [3/3] refbox test"
            refbox test     "${FILTER[@]}" "${EXTRA[@]}"
            ;;
        *)
            echo "unknown step: $step" >&2; exit 2 ;;
    esac
done

echo ">>> done."
