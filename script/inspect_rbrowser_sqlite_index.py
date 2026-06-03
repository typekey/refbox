#!/usr/bin/env python3
"""Inspect an RBrowser SQLite search index.

Prints file size, record counts, SQLite/FTS capabilities, the metadata table
and a few example transcript records.

Example
-------
    python inspect_rbrowser_sqlite_index.py --db hg38.gencode.v45.transcript.rba
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path


def _import_mod():
    """Load ``refbox.sqlite_index`` from the installed package, else directly
    from its self-contained source file (no PyYAML/pysam needed)."""
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


def _human(n: int) -> str:
    f = float(n)
    for unit in ("B", "KB", "MB", "GB", "TB"):
        if f < 1024:
            return f"{f:.1f} {unit}"
        f /= 1024
    return f"{f:.1f} PB"


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description="Inspect an RBrowser SQLite index.")
    ap.add_argument("--db", required=True, help="path to the .rba index")
    args = ap.parse_args(argv)

    m = _import_mod()
    info = m.inspect(Path(args.db))

    print(f"File:            {info['file']}")
    print(f"Size:            {_human(info['file_size_bytes'])} "
          f"({info['file_size_bytes']:,} bytes)")
    print(f"SQLite version:  {info['sqlite_version']}")
    print(f"FTS5 available:  {info['fts5_available']}")
    print(f"Has feature_fts: {info['has_fts']}")
    print(f"Has trigram:     {info['has_trigram']}")
    print()
    print(f"Genes:           {info['n_genes']:,}")
    print(f"Transcripts:     {info['n_transcripts']:,}")
    print(f"Features total:  {info['n_features']:,}")
    print(f"Aliases:         {info['n_aliases']:,}")
    print(f"FTS records:     {info['n_fts']:,}")
    print()
    print("Metadata:")
    for k, v in sorted(info["metadata"].items()):
        print(f"  {k:24s} {v}")
    print()
    print("Example transcript records:")
    for r in info["examples"]:
        print(f"  - {r['transcript_id']} ({r['transcript_name']}) "
              f"{r['gene_name']} {r['chrom']}:{r['start']}-{r['end']}"
              f"({r['strand']}) {r['biotype']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
