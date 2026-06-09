#!/usr/bin/env python3
"""Build a lightweight RBrowser Annotation Index (.rbai) from a GTF/GFF/GFF3.

A ``.rbai`` is a small, fast, static SQLite lookup index (no FTS5/trigram, no
exon/CDS structure) for gene/transcript search + position lookup. Stdlib only.

Example
-------
    python build_lite_rbrowser_annotation_index.py \
        --input annotation.gtf.gz --output annotation.rbai \
        --source-name GENCODE --species human --genome hg38 \
        --annotation-version v45 --enable-gram3 --force --verbose
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path


def _load_module():
    """Load ``refbox.lite_index`` (installed package, else the src/ file)."""
    try:
        import refbox.lite_index as m  # type: ignore
        return m
    except Exception:
        import importlib.util
        path = Path(__file__).resolve().parent.parent / "src" / "refbox" / "lite_index.py"
        spec = importlib.util.spec_from_file_location("refbox_lite_index", path)
        m = importlib.util.module_from_spec(spec)
        sys.modules[spec.name] = m
        spec.loader.exec_module(m)  # type: ignore[union-attr]
        return m


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(
        description="Build a lightweight RBrowser Annotation Index (.rbai).")
    ap.add_argument("--input", required=True, help="GTF/GFF/GFF3 (.gz ok)")
    ap.add_argument("--output", default=None,
                    help="output .rbai path (default: <input>.rbai)")
    ap.add_argument("--source-name", default="", help="e.g. GENCODE, Ensembl")
    ap.add_argument("--species", default="", help="e.g. human, mouse")
    ap.add_argument("--genome", default="", help="e.g. hg38, GRCh38")
    ap.add_argument("--annotation-version", default="", help="e.g. v45")
    ap.add_argument("--enable-gram3", action="store_true",
                    help="also build the 3-gram fuzzy-recall table")
    ap.add_argument("--force", action="store_true", help="overwrite output")
    ap.add_argument("--verbose", action="store_true")
    args = ap.parse_args(argv)

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s %(levelname)s %(name)s: %(message)s")

    m = _load_module()
    try:
        out = m.build_lite_index(
            args.input, args.output, source_name=args.source_name,
            species=args.species, genome=args.genome,
            annotation_version=args.annotation_version,
            enable_gram3=args.enable_gram3, force=args.force,
            verbose=args.verbose)
    except FileExistsError as e:
        print(f"error: {e}", file=sys.stderr)
        return 1
    print(f"OK: {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
