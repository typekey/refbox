#!/usr/bin/env python3
"""Benchmark + sanity-test RBrowser SQLite search.

Runs each query through :func:`refbox.sqlite_index.search` (the ranked search
the browser mirrors), repeats it to measure latency, and prints a per-query
summary plus an aggregate table.

Example
-------
    python test_rbrowser_sqlite_search.py \
        --db hg38.gencode.v45.rbrowser.sqlite \
        --queries TP53 ENST00000269305 p53 BRCA1 MALAT1 ACTB \
        --repeat 100 --limit 10
"""

from __future__ import annotations

import argparse
import statistics
import sys
import time
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


def _fmt_region(d: dict) -> str:
    return f"{d.get('chrom')}:{d.get('start')}-{d.get('end')}({d.get('strand')})"


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description="Benchmark RBrowser SQLite search.")
    ap.add_argument("--db", required=True, help="path to the .sqlite index")
    ap.add_argument("--queries", nargs="+", required=True, help="query strings")
    ap.add_argument("--repeat", type=int, default=100,
                    help="timed repetitions per query (default 100)")
    ap.add_argument("--limit", type=int, default=10, help="max results per query")
    args = ap.parse_args(argv)

    m = _import_mod()
    con = m.open_readonly(Path(args.db))

    rows: list[tuple] = []
    print(f"# DB: {args.db}")
    print(f"# repeat={args.repeat} limit={args.limit}\n")
    for q in args.queries:
        # cold-ish first run (also captures correctness)
        t0 = time.perf_counter()
        res = m.search(con, q, limit=args.limit)
        cold_ms = (time.perf_counter() - t0) * 1e3

        timings: list[float] = []
        for _ in range(max(1, args.repeat)):
            t = time.perf_counter()
            m.search(con, q, limit=args.limit)
            timings.append((time.perf_counter() - t) * 1e3)
        timings.sort()
        avg = statistics.mean(timings)
        med = statistics.median(timings)
        p95 = timings[min(len(timings) - 1, int(round(0.95 * len(timings))) - 1)]
        qps = 1000.0 / avg if avg else float("inf")

        print(f"Query: {q}")
        print(f"Results: {len(res)}")
        print(f"Cold: {cold_ms:.2f} ms")
        print(f"Average: {avg:.3f} ms")
        print(f"Median: {med:.3f} ms")
        print(f"P95: {p95:.3f} ms")
        print(f"QPS (est.): {qps:.0f}")
        if res:
            top = res[0]
            print("Top hit:")
            print(f"  gene_name: {top.get('gene_name')}")
            print(f"  gene_id: {top.get('gene_id')}")
            print(f"  transcript_id: {top.get('transcript_id')}")
            print(f"  transcript_name: {top.get('transcript_name')}")
            print(f"  region: {_fmt_region(top)}")
            print(f"  biotype: {top.get('biotype')}")
            print(f"  matched_field: {top.get('matched_field')}")
        else:
            print("  (no results)")
        print()
        rows.append((q, len(res), cold_ms, avg, med, p95, qps,
                     res[0].get("matched_field") if res else "-"))

    # aggregate table
    print("## Benchmark summary\n")
    print("| Query | Hits | Cold(ms) | Avg(ms) | Median(ms) | P95(ms) | QPS | matched_field |")
    print("| --- | --- | --- | --- | --- | --- | --- | --- |")
    for q, n, cold, avg, med, p95, qps, mf in rows:
        print(f"| {q} | {n} | {cold:.2f} | {avg:.3f} | {med:.3f} | {p95:.3f} "
              f"| {qps:.0f} | {mf} |")
    con.close()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
