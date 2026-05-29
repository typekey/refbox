"""Validate the build/ outputs are usable for genome browser loading."""

from __future__ import annotations

import logging
import subprocess
from dataclasses import dataclass
from pathlib import Path

from .config import Target, iter_targets

log = logging.getLogger(__name__)


@dataclass
class CheckResult:
    name: str
    ok: bool
    detail: str = ""


def _file_ok(p: Path) -> bool:
    return p.exists() and p.stat().st_size > 0


def _check_pair(out_gz: Path, idx_suffix: str, label: str) -> CheckResult:
    if not _file_ok(out_gz):
        return CheckResult(label, False, f"missing: {out_gz.name}")
    idx = Path(str(out_gz) + idx_suffix)
    if not _file_ok(idx):
        return CheckResult(label, False, f"missing index: {idx.name}")
    return CheckResult(label, True)


def _first_chrom(fai: Path) -> str | None:
    if not _file_ok(fai):
        return None
    with open(fai) as f:
        for line in f:
            return line.split("\t", 1)[0]
    return None


def _tabix_query(gz: Path, region: str) -> CheckResult:
    label = f"tabix:{gz.name}:{region}"
    try:
        r = subprocess.run(
            ["tabix", str(gz), region],
            check=True, capture_output=True, text=True, timeout=30,
        )
    except subprocess.CalledProcessError as e:
        return CheckResult(label, False, e.stderr.strip())
    n = sum(1 for _ in r.stdout.splitlines())
    return CheckResult(label, True, f"{n} records")


def _samtools_seq(gz: Path, region: str) -> CheckResult:
    label = f"faidx:{gz.name}:{region}"
    try:
        r = subprocess.run(
            ["samtools", "faidx", str(gz), region],
            check=True, capture_output=True, text=True, timeout=30,
        )
    except subprocess.CalledProcessError as e:
        return CheckResult(label, False, e.stderr.strip())
    lines = [ln for ln in r.stdout.splitlines() if not ln.startswith(">")]
    nbases = sum(len(ln) for ln in lines)
    return CheckResult(label, nbases > 0, f"{nbases} bp")


def check_target(target: Target) -> list[CheckResult]:
    b = target.build_dir
    results: list[CheckResult] = []

    # genome + chrom.sizes
    genome = b / "genome.fa.gz"
    fai = Path(f"{genome}.fai")
    if genome.exists():
        results.append(_check_pair(genome, ".fai", "genome+fai"))
        results.append(_check_pair(genome, ".gzi", "genome+gzi"))
        cs = b / "chrom.sizes"
        if _file_ok(cs):
            n_fai = sum(1 for _ in open(fai))
            n_cs = sum(1 for _ in open(cs))
            results.append(CheckResult(
                "chrom.sizes==fai", n_fai == n_cs, f"{n_cs} vs {n_fai}",
            ))
        else:
            results.append(CheckResult("chrom.sizes", False, "missing"))
        chrom = _first_chrom(fai)
        if chrom:
            results.append(_samtools_seq(genome, f"{chrom}:1-100"))

    # transcriptome (primary + optional derived)
    for name in ("transcriptome.fa.gz", "transcriptome.derived.fa.gz"):
        tx = b / name
        if tx.exists():
            results.append(_check_pair(tx, ".fai", f"{name}+fai"))

    # annotated indices + sample query
    chrom = _first_chrom(fai) if fai.exists() else None
    for name, suffix in [
        ("annotation.sorted.gtf.gz", ".tbi"),
        ("annotation.sorted.gff3.gz", ".tbi"),
        ("repeats.sorted.gtf.gz", ".tbi"),
        ("repeats.sorted.bed.gz", ".tbi"),
        ("rnacentral.sorted.gff3.gz", ".tbi"),
        ("ccre.sorted.bed.gz", ".tbi"),
    ]:
        f = b / name
        if not f.exists():
            continue
        results.append(_check_pair(f, suffix, name))
        if chrom:
            results.append(_tabix_query(f, f"{chrom}:1-1000000"))
    return results


def test_targets(
    species: list[str] | None = None,
    assembly: list[str] | None = None,
    *,
    out: str | None = None,
    include_disabled: bool = False,
) -> int:
    failures = 0
    for tgt in iter_targets(species=species, assembly=assembly, out_root=out,
                            include_disabled=include_disabled):
        print(f"\n=== {tgt.species} / {tgt.assembly} ===")
        for r in check_target(tgt):
            mark = "OK " if r.ok else "FAIL"
            print(f"  [{mark}] {r.name}  {r.detail}")
            if not r.ok:
                failures += 1
    return failures
