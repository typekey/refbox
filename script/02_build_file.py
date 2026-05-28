#!/usr/bin/env python3
"""Step 02 — build bgzip/tabix/faidx outputs from raw/ into build/.

    refbox build --species Homo_sapiens --assembly GRCh38
"""
from __future__ import annotations

import logging

from refbox.build import build_targets

SPECIES: list[str] | None = None
ASSEMBLY: list[str] | None = None
RESOURCES: list[str] | None = None
FORCE = False

if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(name)s: %(message)s")
    build_targets(species=SPECIES, assembly=ASSEMBLY, resources=RESOURCES, force=FORCE)
