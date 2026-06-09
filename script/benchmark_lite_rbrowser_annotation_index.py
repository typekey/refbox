#!/usr/bin/env python3
"""Benchmark a lightweight RBrowser Annotation Index (.rbai).

Reports DB size, record/term/gram counts, and per-query latency statistics
(avg / median / p95 / max) plus an estimated queries-per-second. Stdlib only.

Example
-------
    python benchmark_lite_rbrowser_annotation_index.py \
        --db annotation.rbai \
        --queries TP53 TP5 p53 ENST00000269305 BRCA1 ACTB MALAT1 protein_coding \
        --repeat 100 --limit 10
"""

from __future__ import annotations

import argparse
import statistics
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


def _pct(values, p):
    """p-th percentile (nearest-rank) of a sorted-able list, in the same units."""
    if not values:
        return 0.0
    s = sorted(values)
    k = max(0, min(len(s) - 1, int(round(p / 100.0 * (len(s) - 1)))))
    return s[k]


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description="Benchmark a .rbai lite index.")
    ap.add_argument("--db", required=True)
    ap.add_argument("--queries", nargs="+", required=True)
    ap.add_argument("--repeat", type=int, default=100)
    ap.add_argument("--limit", type=int, default=10)
    args = ap.parse_args(argv)

    m = _load_module()
    info = m.inspect(args.db)
    con = m.open_readonly(args.db)

    # warm the cache / page in the indexes
    for q in args.queries:
        m.search(con, q, limit=args.limit)

    timings_ms: list[float] = []
    per_query: dict[str, list[float]] = {q: [] for q in args.queries}
    for _ in range(args.repeat):
        for q in args.queries:
            t0 = time.perf_counter()
            m.search(con, q, limit=args.limit)
            dt = (time.perf_counter() - t0) * 1000.0
            timings_ms.append(dt)
            per_query[q].append(dt)
    con.close()

    n = len(timings_ms)
    avg = statistics.mean(timings_ms)
    med = statistics.median(timings_ms)
    p95 = _pct(timings_ms, 95)
    mx = max(timings_ms)
    qps = 1000.0 / avg if avg > 0 else float("inf")

    print("=" * 60)
    print(f"DB: {args.db}")
    print(f"  file size:        {info['file_size'] / 1e6:.1f} MB "
          f"({info['file_size']:,} bytes)")
    print(f"  records:          {info['n_records']:,}")
    print(f"    genes:          {info['n_genes']:,}")
    print(f"    transcripts:    {info['n_transcripts']:,}")
    print(f"  terms:            {info['n_terms']:,}")
    print(f"  gram3 rows:       {info['n_gram3']:,}")
    print(f"  index_type:       {info['metadata'].get('index_type')}")
    print("-" * 60)
    print(f"queries:            {len(args.queries)}  × repeat {args.repeat} "
          f"= {n:,} lookups")
    print(f"  average:          {avg:.3f} ms")
    print(f"  median:           {med:.3f} ms")
    print(f"  p95:              {p95:.3f} ms")
    print(f"  max:              {mx:.3f} ms")
    print(f"  throughput:       {qps:,.0f} queries/sec (1 thread)")
    print("-" * 60)
    print("per-query average (ms):")
    for q in args.queries:
        qa = statistics.mean(per_query[q])
        mode, res = m.search(m.open_readonly(args.db), q, limit=args.limit)
        print(f"  {q:24} {qa:7.3f}   mode={mode:10} hits={len(res)}")
    print("=" * 60)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
