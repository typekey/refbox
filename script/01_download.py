#!/usr/bin/env python3
"""Step 01 — download/copy raw reference files into {Species}/{Assembly}/raw/.

Thin wrapper around `refbox.download.download_targets()`. Edit the SPECIES /
ASSEMBLY / RESOURCES lists below or call the package CLI directly:

    refbox download --species Homo_sapiens --assembly GRCh38
"""
from __future__ import annotations

import logging

from refbox.download import download_targets

SPECIES: list[str] | None = None      # None = all enabled
ASSEMBLY: list[str] | None = None
RESOURCES: list[str] | None = None
FORCE = False

if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(name)s: %(message)s")
    download_targets(species=SPECIES, assembly=ASSEMBLY, resources=RESOURCES, force=FORCE)
