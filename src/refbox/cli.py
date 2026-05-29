"""refbox CLI: refbox {download|build|test|import} [filters]

``--assembly`` is required by every subcommand. ``--species`` is optional;
when omitted, the species is looked up from ``species.yaml`` based on the
assembly. If the assembly is not found in the config, the assembly identifier
itself is used as the species/folder name (useful for ``import``).
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
    """If ``--species`` was given, return it. Otherwise infer from --assembly."""
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


def _common_filters(p: argparse.ArgumentParser, *, assembly_required: bool = True) -> None:
    p.add_argument("--species", nargs="*",
                   help="species name (optional; inferred from --assembly via species.yaml)")
    p.add_argument("--assembly", nargs="*", required=assembly_required,
                   help="assembly identifier (e.g. GRCh38)")
    p.add_argument("--resource", nargs="*", dest="resources",
                   help="subset of resources (genome, transcriptome, ...)")
    p.add_argument("--out", default=None,
                   help="output root (default: $REFBOX_OUT or current directory)")


def _parse_map(items: list[str] | None) -> dict[str, Path]:
    out: dict[str, Path] = {}
    for it in items or []:
        if ":" not in it:
            raise SystemExit(f"--map must be resource:path, got: {it}")
        k, v = it.split(":", 1)
        out[k.strip()] = Path(v).expanduser()
    return out


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="refbox")
    parser.add_argument("-v", "--verbose", action="store_true")
    sub = parser.add_subparsers(dest="cmd", required=True)

    p_dl = sub.add_parser("download", help="copy/download raw files")
    _common_filters(p_dl)
    p_dl.add_argument("--force", action="store_true")

    p_bd = sub.add_parser(
        "build",
        help="build indexed outputs (auto-downloads missing raws, then runs test)",
    )
    _common_filters(p_bd)
    p_bd.add_argument("--force", action="store_true")
    p_bd.add_argument("--no-download", action="store_true",
                      help="skip auto-download of missing raw files")
    p_bd.add_argument("--no-test", action="store_true",
                      help="skip the post-build test step")

    p_tt = sub.add_parser("test", help="validate built outputs")
    _common_filters(p_tt)

    p_im = sub.add_parser(
        "import",
        help="ingest a directory of user-supplied files (auto-detected by extension)",
    )
    p_im.add_argument("src", help="directory containing reference files")
    p_im.add_argument("--assembly", required=True,
                      help="assembly identifier (output sub-folder)")
    p_im.add_argument("--species", default=None,
                      help="species name (defaults to species.yaml lookup or the assembly name)")
    p_im.add_argument("--out", default=None, help="output root")
    p_im.add_argument("--map", nargs="*", dest="mapping",
                      help="explicit overrides, e.g. --map genome:/path/to/my.fa")
    p_im.add_argument("--no-build", action="store_true",
                      help="copy raws only; skip the build step")
    p_im.add_argument("--force", action="store_true")

    args = parser.parse_args(argv)
    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
    )

    if args.cmd == "download":
        species = _resolve_species(args.species, args.assembly)
        download_targets(
            species=species, assembly=args.assembly,
            resources=args.resources, out=args.out, force=args.force,
        )
        return 0

    if args.cmd == "build":
        species = _resolve_species(args.species, args.assembly)
        build_targets(
            species=species, assembly=args.assembly,
            resources=args.resources, out=args.out, force=args.force,
            auto_download=not args.no_download,
        )
        if args.no_test:
            return 0
        failed = test_targets(species=species, assembly=args.assembly, out=args.out)
        return 1 if failed else 0

    if args.cmd == "test":
        species = _resolve_species(args.species, args.assembly)
        failed = test_targets(species=species, assembly=args.assembly, out=args.out)
        return 1 if failed else 0

    if args.cmd == "import":
        from .ingest import ingest_directory
        sp = args.species
        if sp is None:
            sp = find_species_by_assembly(args.assembly) or args.assembly
        ingest_directory(
            src_dir=Path(args.src),
            assembly=args.assembly,
            species=sp,
            out=Path(args.out) if args.out else None,
            mapping=_parse_map(args.mapping),
            do_build=not args.no_build,
            force=args.force,
        )
        return 0

    return 0


if __name__ == "__main__":
    sys.exit(main())
