#!/usr/bin/env bash
# Build an RBrowser Index (.rba: SQLite + FTS5 search) for every species/assembly
# that has an annotation file.
#
# Layout (post-reorg):
#   inputs : results/rawdata/<Species>/<Assembly>/<Assembly>.<file>
#   output : results/build/<Species>/<Assembly>/<Assembly>.annotation.rba
#
# Per assembly: prefer GFF3 over GTF; merge the local (already assembly-matched)
# RNAcentral GFF3 when present (normalized preferred); inject HGNC synonyms for
# Homo_sapiens only.
set -u
RESULTS=/home/leizheng/workspace/rbrowser/reference/results
RAWROOT="$RESULTS/rawdata"
BUILDROOT="$RESULTS/build"
HGNC=/home/leizheng/workspace/rbrowser/gene_info/rawdata/hgnc_complete_set.txt
BUILD=/home/leizheng/workspace/rbrowser/reference/refbox/script/build_rbrowser_sqlite_index.py
# Use an interpreter that has refbox importable (the src tree) and its deps.
# Override with e.g.  PY=python3 build_all_species_sqlite.sh
PY="${PY:-/home/leizheng/biotools/mamba/bin/python}"

ok=0; fail=0; skip=0
for rawdir in "$RAWROOT"/*/*; do
  [ -d "$rawdir" ] || continue
  assembly=$(basename "$rawdir")
  species=$(basename "$(dirname "$rawdir")")
  pre="$rawdir/$assembly"                       # <rawdir>/<Assembly>. prefix

  # pick annotation: prefer gff3, else gtf
  annot=""
  [ -s "$pre.annotation_gff3.gff3" ] && annot="$pre.annotation_gff3.gff3"
  [ -z "$annot" ] && [ -s "$pre.annotation_gtf.gtf" ] && annot="$pre.annotation_gtf.gtf"
  if [ -z "$annot" ]; then continue; fi

  # pick rnacentral: prefer normalized
  rna=""
  [ -s "$pre.rnacentral.normalized.gff3" ] && rna="$pre.rnacentral.normalized.gff3"
  [ -z "$rna" ] && [ -s "$pre.rnacentral.gff3" ] && rna="$pre.rnacentral.gff3"

  outdir="$BUILDROOT/$species/$assembly"
  out="$outdir/$assembly.annotation.rba"
  mkdir -p "$outdir"

  args=(--input "$annot" --output "$out" --species "$species" --genome "$assembly"
        --source-name RBrowser --force --verbose)
  [ -n "$rna" ] && args+=(--rnacentral "$rna")
  [ "$species" = "Homo_sapiens" ] && [ -s "$HGNC" ] && args+=(--synonyms "$HGNC")

  echo "=========================================================="
  echo "[build] $species/$assembly  annot=$(basename "$annot")  rna=${rna:+$(basename "$rna")}"
  if "$PY" "$BUILD" "${args[@]}" > "$outdir/$assembly.build_sqlite.log" 2>&1; then
    sz=$(du -h "$out" | cut -f1)
    echo "  OK -> $out ($sz)"
    ok=$((ok+1))
  else
    echo "  FAILED (see $outdir/$assembly.build_sqlite.log)"
    tail -3 "$outdir/$assembly.build_sqlite.log"
    fail=$((fail+1))
  fi
done
echo "=========================================================="
echo "done: ok=$ok fail=$fail skip=$skip"
