"""Ingest a user-supplied directory of reference files into the canonical layout.

The user points us at a folder of pre-existing files (genome FASTA, GTF, BED,
etc.) for an arbitrary assembly. We classify each file by extension, copy it
into ``{Assembly}/raw/<resource>.<canonical-ext>`` under the standard name,
then run the same build pipeline used by configured assemblies.

This lets users provide their own reference (custom organism, in-house build,
unreleased assembly) without having to author a species.yaml entry.

Recognized file extensions (case-insensitive, after stripping ``.gz``/``.bgz``):

  ============= ===================================================
  Resource      Heuristic
  ============= ===================================================
  genome        ``*.fa`` / ``*.fasta`` / ``*.fna``  (largest one)
  transcriptome ``*.fa`` containing 'transcripts' / 'cdna' / 'cds'
  annotation_gtf ``*.gtf``
  annotation_gff3 ``*.gff3`` / ``*.gff``
  repeats_rmsk  ``*.rmsk.txt`` / ``*rmsk*.tsv``
  repeats_bed   ``*repeats*.bed`` / ``*rmsk*.bed``
  repeats_gtf   ``*repeats*.gtf`` / ``*rmsk*.gtf``
  repeats_fa    ``*.fa.out`` (RepeatMasker report)
  rnacentral    ``*rnacentral*.gff3``
  ccre          ``*ccre*.bed`` / ``*cCRE*.bed``
  ============= ===================================================

If the user passes ``--map name:path``, that wins over auto-detection.
"""

from __future__ import annotations

import gzip
import logging
import re
import shutil
from pathlib import Path
from typing import Iterable

from .build import build_targets
from .config import RESOURCE_EXT, RESOURCE_NAMES, Target
from .config import _DEFAULT_OUT as DEFAULT_OUT  # type: ignore[attr-defined]

log = logging.getLogger(__name__)


# ── classification ────────────────────────────────────────────────────────────

_FASTA_EXTS = {".fa", ".fasta", ".fna"}
_GTF_EXTS = {".gtf"}
_GFF3_EXTS = {".gff3", ".gff"}
_BED_EXTS = {".bed"}
_TSV_EXTS = {".tsv", ".txt"}


def _peel_gz(name: str) -> tuple[str, bool]:
    """Return (basename_without_gz, was_gzipped)."""
    lower = name.lower()
    for sfx in (".gz", ".bgz"):
        if lower.endswith(sfx):
            return name[: -len(sfx)], True
    return name, False


def _ext_of(name: str) -> str:
    """Extension after stripping ``.gz``/``.bgz`` (lowercased, with leading dot)."""
    base, _ = _peel_gz(name)
    # handle .fa.out (RepeatMasker)
    lower = base.lower()
    if lower.endswith(".fa.out"):
        return ".fa.out"
    p = Path(base)
    return p.suffix.lower()


def _matches(haystack: str, *needles: str) -> bool:
    h = haystack.lower()
    return any(n in h for n in needles)


def classify_file(path: Path) -> str | None:
    """Return the canonical resource name for ``path`` or None if unknown."""
    name = path.name
    ext = _ext_of(name)

    if ext == ".fa.out":
        return "repeats_fa"

    if ext in _FASTA_EXTS:
        if _matches(name, "transcript", "cdna", "cds", "rna"):
            return "transcriptome"
        return "genome"

    if ext in _GTF_EXTS:
        if _matches(name, "repeat", "rmsk"):
            return "repeats_gtf"
        return "annotation_gtf"

    if ext in _GFF3_EXTS:
        if _matches(name, "rnacentral"):
            return "rnacentral"
        if _matches(name, "repeat", "rmsk"):
            return "repeats_gtf"
        return "annotation_gff3"

    if ext in _BED_EXTS:
        if _matches(name, "ccre"):
            return "ccre"
        if _matches(name, "repeat", "rmsk"):
            return "repeats_bed"
        return "repeats_bed"  # default for unannotated BED

    if ext in _TSV_EXTS:
        if _matches(name, "rmsk"):
            return "repeats_rmsk"

    return None


def _pick_largest(files: Iterable[Path]) -> Path | None:
    files = [f for f in files if f.exists()]
    if not files:
        return None
    return max(files, key=lambda p: p.stat().st_size)


def scan_directory(src_dir: Path) -> dict[str, Path]:
    """Classify every regular file in ``src_dir`` (recursive) into resources.

    When multiple files map to the same resource, the largest one wins, except
    for ``genome`` vs ``transcriptome`` which are decided per-file by name.
    """
    candidates: dict[str, list[Path]] = {}
    for p in sorted(src_dir.rglob("*")):
        if not p.is_file():
            continue
        res = classify_file(p)
        if res:
            candidates.setdefault(res, []).append(p)

    chosen: dict[str, Path] = {}
    for res, paths in candidates.items():
        chosen[res] = _pick_largest(paths) or paths[0]
    return chosen


# ── ingest ────────────────────────────────────────────────────────────────────

def _decompress_to(src: Path, dst: Path) -> None:
    dst.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(src, "rb") as fin, open(dst, "wb") as fout:
        shutil.copyfileobj(fin, fout, length=1 << 20)


def _materialize(src: Path, dst: Path) -> None:
    """Copy (or decompress) ``src`` to the canonical raw ``dst`` path."""
    dst.parent.mkdir(parents=True, exist_ok=True)
    base, was_gz = _peel_gz(src.name)
    if was_gz:
        log.info("decompress %s -> %s", src.name, dst)
        _decompress_to(src, dst)
    else:
        log.info("copy %s -> %s", src.name, dst)
        shutil.copyfile(src, dst)


def ingest_directory(
    src_dir: Path,
    assembly: str,
    *,
    species: str = "Custom",
    out: Path | None = None,
    mapping: dict[str, Path] | None = None,
    do_build: bool = True,
    force: bool = False,
) -> Target:
    """Import a folder of user-supplied files for ``assembly`` and build it.

    Parameters
    ----------
    src_dir:
        Directory to scan for reference files.
    assembly:
        Assembly identifier (used as the output sub-folder).
    species:
        Optional species name; defaults to ``"Custom"``.
    out:
        Output root; falls back to ``$REFBOX_OUT`` or current directory.
    mapping:
        Optional explicit ``{resource: path}`` overrides that win over
        auto-detection.
    do_build:
        When True, run the standard build pipeline after copying raws.
    force:
        Forward ``--force`` to the build step (rebuild even if outputs exist).
    """
    src_dir = Path(src_dir).resolve()
    if not src_dir.is_dir():
        raise FileNotFoundError(f"not a directory: {src_dir}")
    out_root = Path(out).resolve() if out else DEFAULT_OUT

    detected = scan_directory(src_dir)
    if mapping:
        for k, v in mapping.items():
            if k not in RESOURCE_NAMES:
                raise ValueError(f"unknown resource: {k} (expected one of {RESOURCE_NAMES})")
            detected[k] = Path(v).resolve()

    if not detected:
        raise RuntimeError(f"no recognizable reference files in {src_dir}")

    log.info("=== ingest %s / %s (from %s) ===", species, assembly, src_dir)
    for r in RESOURCE_NAMES:
        src = detected.get(r)
        if src is None:
            log.info("  - %-15s  (none)", r)
        else:
            log.info("  + %-15s  %s", r, src)

    # Build a synthetic Target to use its directory layout.
    target = Target(
        species=species,
        assembly=assembly,
        enabled=True,
        resources={r: ({} if r in detected else None) for r in RESOURCE_NAMES},
        meta={},
        out_root=out_root,
    )
    target.raw_dir.mkdir(parents=True, exist_ok=True)

    # Copy/decompress each detected file into the canonical raw/ path.
    for r, src in detected.items():
        dst = target.raw_dir / f"{r}.{RESOURCE_EXT[r]}"
        if dst.exists() and not force:
            log.info("  skip %s (raw exists): %s", r, dst)
            continue
        _materialize(src, dst)

    if do_build:
        build_targets(
            species=[species], assembly=[assembly],
            out=str(out_root), force=force,
            # iter_targets won't find this synthetic Target via species.yaml,
            # so we pass it explicitly via the override hook.
            extra_targets=[target],
        )
    return target
