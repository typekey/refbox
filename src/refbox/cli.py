"""refbox CLI: refbox {download|build|test} [filters]"""

from __future__ import annotations

import argparse
import logging
import sys

from .build import build_targets
from .download import download_targets
from .test import test_targets


def _common_filters(p: argparse.ArgumentParser) -> None:
    p.add_argument("--species", nargs="*", help="filter by species (e.g. Homo_sapiens)")
    p.add_argument("--assembly", nargs="*", help="filter by assembly (e.g. GRCh38)")
    p.add_argument(
        "--resource", nargs="*", dest="resources",
        help="subset of resources (genome, transcriptome, annotation_gtf, ...)",
    )
    p.add_argument(
        "--out", default=None,
        help="output root (default: $REFBOX_OUT or current directory)",
    )


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="refbox")
    parser.add_argument("-v", "--verbose", action="store_true")
    sub = parser.add_subparsers(dest="cmd", required=True)

    p_dl = sub.add_parser("download", help="copy/download raw files")
    _common_filters(p_dl)
    p_dl.add_argument("--force", action="store_true")

    p_bd = sub.add_parser("build", help="build indexed browser files from raw/")
    _common_filters(p_bd)
    p_bd.add_argument("--force", action="store_true")

    p_tt = sub.add_parser("test", help="validate built outputs")
    _common_filters(p_tt)

    args = parser.parse_args(argv)
    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
    )

    if args.cmd == "download":
        download_targets(
            species=args.species, assembly=args.assembly,
            resources=args.resources, out=args.out, force=args.force,
        )
    elif args.cmd == "build":
        build_targets(
            species=args.species, assembly=args.assembly,
            resources=args.resources, out=args.out, force=args.force,
        )
    elif args.cmd == "test":
        failed = test_targets(species=args.species, assembly=args.assembly, out=args.out)
        return 1 if failed else 0
    return 0


if __name__ == "__main__":
    sys.exit(main())
