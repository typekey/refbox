#!/usr/bin/env bash
# Build an RBrowser Index (.rba: SQLite + FTS5 search) for every species/assembly
# that has an annotation file, writing {Species}/{Assembly}/build/annotation.rba.
#
# Per assembly: prefer GFF3 over GTF; merge the local (already assembly-matched)
# RNAcentral GFF3 when present (normalized preferred); inject HGNC synonyms for
# Homo_sapiens only.
set -u
RESULTS=/home/leizheng/workspace/rbrowser/reference/results
HGNC=/home/leizheng/workspace/rbrowser/gene_info/rawdata/hgnc_complete_set.txt
BUILD=/home/leizheng/workspace/rbrowser/reference/refbox/script/build_rbrowser_sqlite_index.py
# Use an interpreter that has refbox importable (the src tree) and its deps.
# Override with e.g.  PY=python3 build_all_species_sqlite.sh
PY="${PY:-/home/leizheng/biotools/mamba/bin/python}"

ok=0; fail=0; skip=0
for raw in "$RESULTS"/*/*/raw; do
  asmdir=$(dirname "$raw")
  assembly=$(basename "$asmdir")
  species=$(basename "$(dirname "$asmdir")")

  # pick annotation: prefer gff3, else gtf
  annot=""
  [ -s "$raw/annotation_gff3.gff3" ] && annot="$raw/annotation_gff3.gff3"
  [ -z "$annot" ] && [ -s "$raw/annotation_gtf.gtf" ] && annot="$raw/annotation_gtf.gtf"
  if [ -z "$annot" ]; then continue; fi

  # pick rnacentral: prefer normalized
  rna=""
  [ -s "$raw/rnacentral.normalized.gff3" ] && rna="$raw/rnacentral.normalized.gff3"
  [ -z "$rna" ] && [ -s "$raw/rnacentral.gff3" ] && rna="$raw/rnacentral.gff3"

  out="$asmdir/build/annotation.rba"
  mkdir -p "$asmdir/build"

  args=(--input "$annot" --output "$out" --species "$species" --genome "$assembly"
        --source-name RBrowser --force --verbose)
  [ -n "$rna" ] && args+=(--rnacentral "$rna")
  [ "$species" = "Homo_sapiens" ] && [ -s "$HGNC" ] && args+=(--synonyms "$HGNC")

  echo "=========================================================="
  echo "[build] $species/$assembly  annot=$(basename "$annot")  rna=${rna:+$(basename "$rna")}"
  if "$PY" "$BUILD" "${args[@]}" > "$asmdir/build/build_sqlite.log" 2>&1; then
    sz=$(du -h "$out" | cut -f1)
    echo "  OK -> $out ($sz)"
    ok=$((ok+1))
  else
    echo "  FAILED (see $asmdir/build/build_sqlite.log)"
    tail -3 "$asmdir/build/build_sqlite.log"
    fail=$((fail+1))
  fi
done
echo "=========================================================="
echo "done: ok=$ok fail=$fail skip=$skip"
