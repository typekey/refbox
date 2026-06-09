#!/usr/bin/env python3
"""Search a lightweight RBrowser Annotation Index (.rbai).

Shows, per query, the resolution mode (exact / prefix / normalized / gram3), the
wall time in ms, the number of results, and the top hit(s). Stdlib only.

Example
-------
    python search_lite_rbrowser_annotation_index.py \
        --db annotation.rbai \
        --queries TP53 TP5 p53 ENST00000269305 ENST00000269305.10 BRCA1 protein_coding \
        --limit 10
"""

from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path


def _load_module():
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


_SHOW = ("feature_type", "gene_name", "gene_id", "transcript_name",
         "transcript_id", "transcript_biotype", "genome_position_str")


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description="Search a .rbai lite index.")
    ap.add_argument("--db", required=True, help="path to the .rbai index")
    ap.add_argument("--queries", nargs="+", required=True, help="query strings")
    ap.add_argument("--limit", type=int, default=10)
    ap.add_argument("--show", type=int, default=1,
                    help="how many top results to print per query")
    args = ap.parse_args(argv)

    m = _load_module()
    con = m.open_readonly(args.db)

    for q in args.queries:
        t0 = time.perf_counter()
        mode, results = m.search(con, q, limit=args.limit)
        dt_ms = (time.perf_counter() - t0) * 1000.0
        print(f"\nQuery: {q}")
        print(f"Mode: {mode}")
        print(f"Time: {dt_ms:.2f} ms")
        print(f"Results: {len(results)}")
        for i, r in enumerate(results[: args.show]):
            label = "Top hit:" if i == 0 else f"Hit {i + 1}:"
            print(label)
            for col in _SHOW:
                print(f"  {col}: {r.get(col) or ''}")
    con.close()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
