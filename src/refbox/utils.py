"""Thin wrappers around external CLI tools: bgzip, tabix, samtools, sort."""

from __future__ import annotations

import logging
import shutil
import subprocess
from pathlib import Path

log = logging.getLogger(__name__)


def run(cmd: list[str] | str, *, shell: bool = False, cwd: Path | None = None) -> None:
    """Run command and raise on failure. Logs the command."""
    if isinstance(cmd, list):
        printable = " ".join(map(str, cmd))
    else:
        printable = cmd
    log.info("$ %s", printable)
    subprocess.run(cmd, check=True, shell=shell, cwd=cwd)


def require_tool(name: str) -> None:
    if shutil.which(name) is None:
        raise RuntimeError(f"required external tool not found in PATH: {name}")


def bgzip_file(src: Path, dst: Path, *, force: bool = False, threads: int = 4) -> Path:
    """bgzip `src` into `dst` (.gz). `src` is preserved."""
    require_tool("bgzip")
    dst.parent.mkdir(parents=True, exist_ok=True)
    if dst.exists() and not force:
        return dst
    # bgzip -c writes to stdout
    with open(dst, "wb") as f:
        subprocess.run(
            ["bgzip", "-c", "-@", str(threads), str(src)],
            check=True, stdout=f,
        )
    return dst


def faidx(fa_gz: Path, *, force: bool = False) -> Path:
    """samtools faidx; produces .fai and (for bgzipped) .gzi."""
    require_tool("samtools")
    fai = Path(str(fa_gz) + ".fai")
    if fai.exists() and not force:
        return fai
    run(["samtools", "faidx", str(fa_gz)])
    return fai


def tabix(gz: Path, *, preset: str, force: bool = False) -> Path:
    """tabix index; preset = 'gff' | 'bed' | 'vcf'."""
    require_tool("tabix")
    tbi = Path(str(gz) + ".tbi")
    if tbi.exists() and not force:
        return tbi
    run(["tabix", "-f", "-p", preset, str(gz)])
    return tbi


def sort_gff(src: Path, dst: Path) -> Path:
    """Sort GFF/GTF by chrom, then start. Keeps comment header (## lines)."""
    dst.parent.mkdir(parents=True, exist_ok=True)
    # grep -v '^#' to drop headers from sort body, then prepend headers
    cmd = (
        f"(grep '^#' {src!s} || true; "
        f"grep -v '^#' {src!s} | sort -k1,1 -k4,4n) > {dst!s}"
    )
    run(cmd, shell=True)
    return dst


def sort_bed(src: Path, dst: Path) -> Path:
    """Sort BED by chrom then start. Drops lines beginning with '#' or 'track'."""
    dst.parent.mkdir(parents=True, exist_ok=True)
    cmd = (
        f"grep -v -E '^(#|track|browser)' {src!s} | "
        f"sort -k1,1 -k2,2n > {dst!s}"
    )
    run(cmd, shell=True)
    return dst


def write_chrom_sizes(fai: Path, dst: Path) -> Path:
    """Write chrom.sizes from a .fai file (first two columns)."""
    dst.parent.mkdir(parents=True, exist_ok=True)
    with open(fai) as fin, open(dst, "w") as fout:
        for line in fin:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 2:
                fout.write(f"{parts[0]}\t{parts[1]}\n")
    return dst
