"""refbox CLI.

Subcommands
-----------
``refbox download`` — only fetch raw files configured in species.yaml.
``refbox pull``     — full pipeline for a configured assembly: download (if
                      missing) + build + test. This is what most users want.
``refbox test``     — re-run the validators against existing build/ outputs.
``refbox build``    — single-file / single-directory build for arbitrary
                      user-supplied inputs (no species.yaml entry needed):

    refbox build -fa  GENOME.fa [-o OUT.fa.gz]
    refbox build -gtf ANNOT.gtf [-o OUT.gtf.gz]
    refbox build -gff ANNOT.gff3 [-o OUT.gff3.gz]
    refbox build -bed FEATURES.bed [-o OUT.bed.gz] [--chrom-sizes FILE | --assembly NAME]
    refbox build -rmsk rmsk.txt.gz [-o OUT_DIR]
    refbox build -fa GENOME.fa -gtf ANNOT.gtf -o transcriptome.fa.gz # extract transcriptome
    refbox build -i DIR --assembly GRCh38 [--species NAME]            # directory ingest
    refbox build SOMEFILE        # auto-detect (-i for directories)
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

from .build import build_targets
from .config import find_species_by_assembly
from .download import download_targets
from .test import test_targets


def _resolve_species(
    species_args: list[str] | None, assembly_args: list[str] | None,
) -> list[str] | None:
    if species_args:
        return species_args
    if not assembly_args:
        return None
    found: list[str] = []
    for asm in assembly_args:
        sp = find_species_by_assembly(asm)
        if sp and sp not in found:
            found.append(sp)
    return found or None


def _common_filters(p: argparse.ArgumentParser) -> None:
    p.add_argument("--species", nargs="*",
                   help="species name (optional; inferred from --assembly)")
    p.add_argument("--assembly", nargs="*", default=None,
                   help="assembly identifier (e.g. GRCh38). Omit to run every "
                        "assembly that matches --species (or all assemblies if "
                        "neither is given).")
    p.add_argument("--resource", nargs="*", dest="resources",
                   help="subset of resources (genome, transcriptome, ...)")
    p.add_argument("--out", default=None,
                   help="output root (default: $REFBOX_OUT or current directory)")
    p.add_argument("--include-disabled", action="store_true",
                   help="include assemblies marked enabled: false in species.yaml")


def _parse_map(items: list[str] | None) -> dict[str, Path]:
    out: dict[str, Path] = {}
    for it in items or []:
        if ":" not in it:
            raise SystemExit(f"--map must be resource:path, got: {it}")
        k, v = it.split(":", 1)
        out[k.strip()] = Path(v).expanduser()
    return out


def _dispatch_build(args: argparse.Namespace) -> int:
    """Run the single-file / directory build (the ``refbox build`` command)."""
    from . import file_build as fb

    in_path: Path | None = Path(args.input) if args.input else None
    out_path: Path | None = Path(args.out) if args.out else None

    fa = Path(args.fa) if args.fa else None
    gtf = Path(args.gtf) if args.gtf else None
    gff = Path(args.gff) if args.gff else None
    bed = Path(args.bed) if args.bed else None
    rmsk = Path(args.rmsk) if args.rmsk else None
    ingest_dir = Path(args.ingest) if args.ingest else None

    if ingest_dir is None and in_path is not None and in_path.is_dir():
        ingest_dir = in_path
        in_path = None

    if ingest_dir is not None:
        if not args.assembly:
            raise SystemExit("-i/--ingest requires --assembly")
        from .ingest import ingest_directory
        sp = args.species or find_species_by_assembly(args.assembly) or args.assembly
        ingest_directory(
            src_dir=ingest_dir, assembly=args.assembly, species=sp,
            out=out_path, mapping=_parse_map(args.mapping),
            do_build=not args.no_build, force=args.force,
        )
        return 0

    sqlite_in = Path(args.sqlite) if args.sqlite else None

    annot = gtf or gff
    # Standalone SQLite index build: `refbox build -sqlite ANNOT.gtf[.gz]`
    if sqlite_in is not None:
        from .sqlite_index import build_sqlite_index
        build_sqlite_index(
            sqlite_in, out_path, source_name=args.source_name or "",
            species=args.species_name or "", genome=args.genome or "",
            annotation_version=args.annotation_version or "",
            synonyms=args.synonyms_file, force=args.force, verbose=args.verbose,
        )
        return 0
    if fa and annot:
        fb.build_transcriptome(fa, annot, out_path)
        return 0
    if fa:
        fb.build_fa(fa, out_path)
        return 0
    if annot:
        fb.build_gxf(
            annot, out_path, sqlite=args.with_sqlite,
            source_name=args.source_name or "", species=args.species_name or "",
            genome=args.genome or "", annotation_version=args.annotation_version or "",
            synonyms=args.synonyms_file, force=args.force,
        )
        return 0
    if bed:
        fb.build_bed(
            bed, out_path,
            chrom_sizes=Path(args.chrom_sizes) if args.chrom_sizes else None,
            assembly=args.assembly,
            make_bigbed=not args.no_bigbed,
        )
        return 0
    if rmsk:
        fb.build_rmsk(rmsk, out_path)
        return 0

    if in_path is None:
        raise SystemExit(
            "refbox build: no input. Provide -fa / -gtf / -gff / -bed / -rmsk, "
            "or pass a file/dir as the positional argument (auto-detect)."
        )

    # auto-detect mode on a single file
    kind = fb.auto_detect(in_path)
    if kind == "fa":
        fb.build_fa(in_path, out_path)
    elif kind == "gxf":
        fb.build_gxf(in_path, out_path)
    elif kind == "bed":
        fb.build_bed(
            in_path, out_path,
            chrom_sizes=Path(args.chrom_sizes) if args.chrom_sizes else None,
            assembly=args.assembly,
            make_bigbed=not args.no_bigbed,
        )
    elif kind == "rmsk":
        fb.build_rmsk(in_path, out_path)
    else:
        raise SystemExit(f"unsupported auto-detected kind: {kind}")
    return 0


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="refbox")
    parser.add_argument("-v", "--verbose", action="store_true")
    sub = parser.add_subparsers(dest="cmd", required=True)

    p_dl = sub.add_parser("download", help="copy/download raw files only")
    _common_filters(p_dl)
    p_dl.add_argument("--force", action="store_true")

    p_pull = sub.add_parser(
        "pull",
        help="full pipeline for configured assemblies: download + build + test",
    )
    _common_filters(p_pull)
    p_pull.add_argument("--force", action="store_true")
    p_pull.add_argument("--no-download", action="store_true",
                        help="skip auto-download of missing raw files")
    p_pull.add_argument("--no-test", action="store_true",
                        help="skip the post-build test step")

    p_tt = sub.add_parser("test", help="validate built outputs")
    _common_filters(p_tt)

    p_bd = sub.add_parser(
        "build",
        help="single-file / directory build for arbitrary inputs",
    )
    p_bd.add_argument("input", nargs="?", help="input file (auto-detect) or dir (-i)")
    p_bd.add_argument("-o", "--out", default=None, help="output path (file or dir)")
    p_bd.add_argument("-fa", "--fa", help="FASTA input (genome)")
    p_bd.add_argument("-gtf", "--gtf", help="GTF input")
    p_bd.add_argument("-gff", "--gff", help="GFF3 input")
    p_bd.add_argument("-bed", "--bed", help="BED input")
    p_bd.add_argument("-rmsk", "--rmsk", help="UCSC rmsk.txt[.gz] input")
    p_bd.add_argument("-sqlite", "--sqlite", default=None,
                      help="GTF/GFF3 input → build a standalone SQLite search index")
    p_bd.add_argument("-i", "--ingest", help="directory of user files to import")
    p_bd.add_argument("--assembly", default=None,
                      help="assembly identifier (folder name / chrom.sizes lookup)")
    p_bd.add_argument("--species", default=None, help="species name (-i only)")
    p_bd.add_argument("--with-sqlite", action="store_true",
                      help="for -gtf/-gff: also emit a SQLite search index "
                           "alongside the sorted/bgzip/tabix outputs")
    p_bd.add_argument("--source-name", default=None,
                      help="annotation source label stored in SQLite metadata "
                           "(e.g. GENCODE, Ensembl)")
    p_bd.add_argument("--species-name", default=None,
                      help="species label stored in SQLite metadata")
    p_bd.add_argument("--genome", default=None,
                      help="genome/assembly label stored in SQLite metadata (e.g. hg38)")
    p_bd.add_argument("--annotation-version", default=None,
                      help="annotation version stored in SQLite metadata (e.g. v45)")
    p_bd.add_argument("--synonyms-file", "--synonyms", default=None, dest="synonyms_file",
                      help="HGNC-style TSV (symbol/alias_symbol/prev_symbol/"
                           "ensembl_gene_id) to inject as gene_synonym aliases "
                           "(SQLite builds only)")
    p_bd.add_argument("--chrom-sizes", default=None,
                      help="chrom.sizes file for bigBed conversion")
    p_bd.add_argument("--no-bigbed", action="store_true",
                      help="skip bigBed generation for BED inputs")
    p_bd.add_argument("--map", nargs="*", dest="mapping",
                      help="-i only: resource:path overrides")
    p_bd.add_argument("--no-build", action="store_true",
                      help="-i only: copy raws but skip build pipeline")
    p_bd.add_argument("--force", action="store_true")

    args = parser.parse_args(argv)
    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
    )

    if args.cmd == "download":
        species = _resolve_species(args.species, args.assembly)
        download_targets(species=species, assembly=args.assembly,
                         resources=args.resources, out=args.out, force=args.force,
                         include_disabled=args.include_disabled)
        return 0
    if args.cmd == "pull":
        species = _resolve_species(args.species, args.assembly)
        build_targets(species=species, assembly=args.assembly,
                      resources=args.resources, out=args.out, force=args.force,
                      auto_download=not args.no_download,
                      include_disabled=args.include_disabled)
        if args.no_test:
            return 0
        failed = test_targets(species=species, assembly=args.assembly, out=args.out,
                              include_disabled=args.include_disabled)
        return 1 if failed else 0
    if args.cmd == "test":
        species = _resolve_species(args.species, args.assembly)
        failed = test_targets(species=species, assembly=args.assembly, out=args.out,
                              include_disabled=args.include_disabled)
        return 1 if failed else 0
    if args.cmd == "build":
        return _dispatch_build(args)
    return 0


if __name__ == "__main__":
    sys.exit(main())
