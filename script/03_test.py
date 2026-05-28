#!/usr/bin/env python3
"""Step 03 — validate build/ outputs (file presence, indices, sample queries).

    refbox test --species Homo_sapiens --assembly GRCh38

Exit code: non-zero on any failure.
"""
from __future__ import annotations

import logging
import sys

from refbox.test import test_targets

SPECIES: list[str] | None = None
ASSEMBLY: list[str] | None = None

if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(name)s: %(message)s")
    failed = test_targets(species=SPECIES, assembly=ASSEMBLY)
    sys.exit(1 if failed else 0)
