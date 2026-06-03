#!/usr/bin/env python3
"""Build a static RBrowser Index (.rbi) from a GTF/GFF/GFF3 file.

A ``.rbi`` ("RBrowser Index") is a self-contained SQLite + FTS5 search index,
served as a static file and queried in-browser via SQLite-WASM / HTTP Range VFS.

Standalone CLI wrapper around :mod:`refbox.sqlite_index`. It works whether or
not ``refbox`` is pip-installed: if the package import fails it falls back to
the ``src/`` tree next to this script.

Example
-------
    python build_rbrowser_sqlite_index.py \
        --input gencode.v45.annotation.gtf.gz \
        --output hg38.gencode.v45.transcript.rbi \
        --source-name GENCODE --species human --genome hg38 \
        --annotation-version v45 --force --verbose
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path


def _load_module():
    """Load ``refbox.sqlite_index`` — preferring an installed package, else
    loading the self-contained module file directly (so the script has no
    dependency on the rest of the package, which needs PyYAML/pysam)."""
    try:
        import refbox.sqlite_index as m  # type: ignore
        return m
    except Exception:
        import importlib.util
        path = Path(__file__).resolve().parent.parent / "src" / "refbox" / "sqlite_index.py"
        spec = importlib.util.spec_from_file_location("refbox_sqlite_index", path)
        m = importlib.util.module_from_spec(spec)
        sys.modules[spec.name] = m  # dataclasses needs the module registered
        spec.loader.exec_module(m)  # type: ignore[union-attr]
        return m


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(
        description="Build a read-only RBrowser Index (.rbi) from GTF/GFF3.")
    ap.add_argument("--input", required=True, help="GTF/GFF/GFF3 (.gz ok)")
    ap.add_argument("--output", default=None,
                    help="output .rbi path (default: <input>.rbi)")
    ap.add_argument("--source-name", default="", help="e.g. GENCODE, Ensembl")
    ap.add_argument("--species", default="", help="e.g. human, mouse")
    ap.add_argument("--genome", default="", help="e.g. hg38, GRCh38")
    ap.add_argument("--annotation-version", default="", help="e.g. v45, 112")
    ap.add_argument("--synonyms", default=None,
                    help="HGNC-style TSV (symbol/alias_symbol/prev_symbol/"
                         "ensembl_gene_id) → injected as gene_synonym aliases "
                         "(e.g. OCT4 → POU5F1)")
    ap.add_argument("--rnacentral", default=None,
                    help="RNAcentral genome-coordinates GFF3 (use the chrom-"
                         "normalized one) → merged in as ncRNA records")
    ap.add_argument("--fuzzy-scope", choices=["names", "all"], default="names",
                    help="trigram substring corpus: 'names' (default; gene/"
                         "transcript names + synonyms only, IDs excluded — much "
                         "smaller/faster) or 'all' (include IDs)")
    ap.add_argument("--force", action="store_true", help="overwrite existing output")
    ap.add_argument("--verbose", action="store_true", help="DEBUG logging")
    args = ap.parse_args(argv)

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s %(levelname)s %(name)s: %(message)s")

    build_sqlite_index = _load_module().build_sqlite_index
    out = build_sqlite_index(
        Path(args.input),
        Path(args.output) if args.output else None,
        source_name=args.source_name, species=args.species, genome=args.genome,
        annotation_version=args.annotation_version, synonyms=args.synonyms,
        rnacentral=args.rnacentral, fuzzy_scope=args.fuzzy_scope,
        force=args.force, verbose=args.verbose,
    )
    print(f"OK: {out} ({out.stat().st_size / 1e6:.1f} MB)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
