#!/usr/bin/env python3
"""Clean RNAcentral genome-coordinates GFF3 names.

The RNAcentral GFF3 carries only a long free-text ``description`` attribute
(no short symbol), e.g.::

    description=DEAD/H-box helicase 11 like 11 (pseudogene)%2C transcript variant 1 (DDX11L11)
    description=(human) non-protein coding lnc-OR4F29-11:7
    description=(human) mir-571 microRNA precursor family
    description=(human) hsa-miR-34a-5p
    description=(human) Homo_sapiens piRNA piR-hsa-4818588

This tool distills a short, recognizable ``gene_name`` per RNA type and writes
it back as a new ``gene_name=`` attribute (the original ``description`` is kept
untouched, so nothing is lost):

    DDX11L11 / lnc-OR4F29-11-7 / pre-mir-571 / miR-34a-5p / piR-hsa-4818588

It is stdlib-only and reuses the exact cleaning logic from
``refbox.sqlite_index.clean_rnacentral_name`` (so a standalone cleaned GFF and
the SQLite index agree byte-for-byte).

Examples
--------
    # rewrite a GFF3, injecting gene_name=
    python clean_rnacentral_gff.py --input rnacentral.gff3 --output rnacentral.clean.gff3

    # just preview the name mapping (URS, type, description -> short) as TSV
    python clean_rnacentral_gff.py --input rnacentral.gff3 --preview | head
"""

from __future__ import annotations

import argparse
import gzip
import re
import sys
from pathlib import Path


def _load_module():
    """Load ``refbox.sqlite_index`` — installed package preferred, else the
    self-contained module file next to this script (no PyYAML/pysam needed)."""
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


def _open_text(path: Path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "rt")


def _open_out(path: str | None):
    if path is None or path == "-":
        return sys.stdout
    if path.endswith(".gz"):
        return gzip.open(path, "wt")
    return open(path, "wt")


_ATTR = re.compile(r"(\w+)=([^;]*)")


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description="Clean RNAcentral GFF3 names.")
    ap.add_argument("--input", required=True, help="RNAcentral GFF3 (.gz ok)")
    ap.add_argument("--output", default="-",
                    help="cleaned GFF3 out (.gz ok; default: stdout)")
    ap.add_argument("--preview", action="store_true",
                    help="instead of a GFF, emit a TSV: urs<TAB>type<TAB>"
                         "short_name<TAB>description (one row per transcript)")
    args = ap.parse_args(argv)

    m = _load_module()
    clean = m.clean_rnacentral_name
    display_name = m.rnacentral_display_name

    inp = Path(args.input)
    n_tx = 0
    n_short = 0
    with _open_text(inp) as fh, _open_out(args.output) as out:
        for raw in fh:
            if not raw or raw[0] == "#":
                if not args.preview:
                    out.write(raw)
                continue
            cols = raw.rstrip("\n").split("\t")
            if len(cols) < 9:
                if not args.preview:
                    out.write(raw)
                continue
            attrs = dict(_ATTR.findall(cols[8]))
            rtype = attrs.get("type") or "ncRNA"
            desc = attrs.get("description") or ""
            urs = attrs.get("Name") or attrs.get("ID") or ""
            short, full = clean(desc, rtype)
            name = display_name(short, full, urs)   # what becomes gene_name

            if cols[2] == "transcript":
                n_tx += 1
                if name != urs:
                    n_short += 1
                if args.preview:
                    out.write(f"{urs}\t{rtype}\t{name}\t{full}\n")
                    continue

            if not args.preview:
                if name and "gene_name=" not in cols[8]:
                    # inject gene_name right after the description (or at end)
                    cols[8] = cols[8].rstrip(";") + f";gene_name={name}"
                out.write("\t".join(cols) + "\n")

    sys.stderr.write(
        f"transcripts: {n_tx}  named by symbol: {n_short} "
        f"(URS fallback: {n_tx - n_short}; "
        f"{100 * n_short // max(n_tx, 1)}% symbol)\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
