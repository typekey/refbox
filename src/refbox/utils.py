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


def _is_gzip_magic(path: Path) -> bool:
    try:
        with open(path, "rb") as f:
            return f.read(2) == b"\x1f\x8b"
    except OSError:
        return False


def _is_bgzip(path: Path) -> bool:
    """A bgzip file is gzip with a BC subfield in its first member's extra field."""
    if not _is_gzip_magic(path):
        return False
    try:
        with open(path, "rb") as f:
            head = f.read(18)
        # gzip header: 0x1f 0x8b ID2(0x08) FLG MTIME(4) XFL OS  XLEN(2) SI1 SI2
        if len(head) < 18 or head[3] & 0x04 == 0:  # FEXTRA bit
            return False
        return head[12:14] == b"BC"
    except OSError:
        return False


def bgzip_file(src: Path, dst: Path, *, force: bool = False, threads: int = 4) -> Path:
    """bgzip ``src`` into ``dst`` (.gz). ``src`` is preserved.

    If ``src`` is already a bgzip file, it is copied directly.
    If ``src`` is a plain gzip (not bgzip), it is transparently re-bgzipped
    so downstream tabix/faidx work.
    """
    require_tool("bgzip")
    dst.parent.mkdir(parents=True, exist_ok=True)
    if dst.exists() and not force:
        return dst
    tmp = dst.with_suffix(dst.suffix + ".tmp")
    if tmp.exists():
        tmp.unlink()

    if _is_bgzip(src):
        # Already bgzip: just copy it (preserves byte-identical index compat).
        shutil.copyfile(src, tmp)
    elif _is_gzip_magic(src):
        # Plain gzip — pipe gunzip into bgzip.
        log.info("re-bgzip (input was plain gzip): %s", src.name)
        with open(tmp, "wb") as f:
            p1 = subprocess.Popen(["gunzip", "-c", str(src)], stdout=subprocess.PIPE)
            p2 = subprocess.Popen(
                ["bgzip", "-c", "-@", str(threads)],
                stdin=p1.stdout, stdout=f,
            )
            p1.stdout.close()  # type: ignore[union-attr]
            rc2 = p2.wait()
            rc1 = p1.wait()
            if rc1 != 0 or rc2 != 0:
                tmp.unlink(missing_ok=True)
                raise RuntimeError(f"re-bgzip failed (gunzip={rc1}, bgzip={rc2})")
    else:
        with open(tmp, "wb") as f:
            subprocess.run(
                ["bgzip", "-c", "-@", str(threads), str(src)],
                check=True, stdout=f,
            )
    tmp.replace(dst)
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
