"""Build standardized browser-loadable files from raw/ into build/.

Each builder is idempotent: if the final indexed output already exists, it is
skipped unless `force=True`. Missing raw inputs are skipped silently.

Final output names (per assembly build/):
  genome.fa.gz + .gzi + .fai
  chrom.sizes
  transcripts.fa.gz + .fai
  annotation.sorted.gtf.gz + .tbi
  annotation.sorted.gff3.gz + .tbi
  repeats.sorted.gtf.gz + .tbi
  repeats.sorted.bed.gz + .tbi
  rnacentral.sorted.gff3.gz + .tbi
  ccre.sorted.bed.gz + .tbi
"""

from __future__ import annotations

import logging
from pathlib import Path

from .config import RESOURCE_NAMES, Target, iter_targets, raw_path
from .utils import (
    bgzip_file,
    faidx,
    sort_bed,
    sort_gff,
    tabix,
    write_chrom_sizes,
)

log = logging.getLogger(__name__)


def _exists(p: Path) -> bool:
    return p.exists() and p.stat().st_size > 0


# ── fasta builders ─────────────────────────────────────────────────────────────

def build_genome(target: Target, *, force: bool = False) -> None:
    src = raw_path(target, "genome")
    if not src.exists():
        return
    out = target.build_dir / "genome.fa.gz"
    if _exists(out) and _exists(Path(f"{out}.fai")) and not force:
        log.info("[%s/%s] genome up-to-date", target.species, target.assembly)
    else:
        bgzip_file(src, out, force=force)
        faidx(out, force=force)
    chrom_sizes = target.build_dir / "chrom.sizes"
    if not _exists(chrom_sizes) or force:
        write_chrom_sizes(Path(f"{out}.fai"), chrom_sizes)


def build_transcriptome(target: Target, *, force: bool = False) -> None:
    src = raw_path(target, "transcriptome")
    if not src.exists():
        return
    out = target.build_dir / "transcripts.fa.gz"
    if _exists(out) and _exists(Path(f"{out}.fai")) and not force:
        log.info("[%s/%s] transcripts up-to-date", target.species, target.assembly)
        return
    bgzip_file(src, out, force=force)
    faidx(out, force=force)


# ── annotation builders (generic GFF/GTF + BED) ────────────────────────────────

def _build_sorted_gff(src: Path, out_gz: Path, *, preset: str = "gff",
                      force: bool = False) -> None:
    if _exists(out_gz) and _exists(Path(f"{out_gz}.tbi")) and not force:
        return
    sorted_tmp = out_gz.with_suffix("")  # drop .gz
    sort_gff(src, sorted_tmp)
    bgzip_file(sorted_tmp, out_gz, force=True)
    sorted_tmp.unlink(missing_ok=True)
    tabix(out_gz, preset=preset, force=True)


def _build_sorted_bed(src: Path, out_gz: Path, *, force: bool = False) -> None:
    if _exists(out_gz) and _exists(Path(f"{out_gz}.tbi")) and not force:
        return
    sorted_tmp = out_gz.with_suffix("")
    sort_bed(src, sorted_tmp)
    bgzip_file(sorted_tmp, out_gz, force=True)
    sorted_tmp.unlink(missing_ok=True)
    tabix(out_gz, preset="bed", force=True)


def build_annotation_gtf(target: Target, *, force: bool = False) -> None:
    src = raw_path(target, "annotation_gtf")
    if not src.exists():
        return
    _build_sorted_gff(src, target.build_dir / "annotation.sorted.gtf.gz", force=force)


def build_annotation_gff3(target: Target, *, force: bool = False) -> None:
    src = raw_path(target, "annotation_gff3")
    if not src.exists():
        return
    _build_sorted_gff(src, target.build_dir / "annotation.sorted.gff3.gz", force=force)


def build_repeats_gtf(target: Target, *, force: bool = False) -> None:
    src = raw_path(target, "repeats_gtf")
    if not src.exists():
        return
    _build_sorted_gff(src, target.build_dir / "repeats.sorted.gtf.gz", force=force)


def build_repeats_bed(target: Target, *, force: bool = False) -> None:
    src = raw_path(target, "repeats_bed")
    if not src.exists():
        return
    _build_sorted_bed(src, target.build_dir / "repeats.sorted.bed.gz", force=force)


def build_rnacentral(target: Target, *, force: bool = False) -> None:
    src = raw_path(target, "rnacentral")
    if not src.exists():
        return
    # raw filename is rnacentral.gff3; some sources actually ship gtf-like.
    # We treat extension generically with tabix -p gff.
    _build_sorted_gff(src, target.build_dir / "rnacentral.sorted.gff3.gz", force=force)


def build_ccre(target: Target, *, force: bool = False) -> None:
    src = raw_path(target, "ccre")
    if not src.exists():
        return
    _build_sorted_bed(src, target.build_dir / "ccre.sorted.bed.gz", force=force)


BUILDERS = {
    "genome":          build_genome,
    "transcriptome":   build_transcriptome,
    "annotation_gtf":  build_annotation_gtf,
    "annotation_gff3": build_annotation_gff3,
    "repeats_rmsk":    None,   # raw TSV — converted to bed/gtf in repeats build step (TODO)
    "repeats_gtf":     build_repeats_gtf,
    "repeats_bed":     build_repeats_bed,
    "repeats_fa":      None,   # RepeatMasker .fa.out is a flat report — handled later
    "rnacentral":      build_rnacentral,
    "ccre":            build_ccre,
}


def build_targets(
    species: list[str] | None = None,
    assembly: list[str] | None = None,
    resources: list[str] | None = None,
    *,
    out: str | None = None,
    force: bool = False,
) -> None:
    resources = resources or RESOURCE_NAMES
    for tgt in iter_targets(species=species, assembly=assembly, out_root=out):
        log.info("=== build %s / %s ===", tgt.species, tgt.assembly)
        tgt.build_dir.mkdir(parents=True, exist_ok=True)
        for r in resources:
            fn = BUILDERS.get(r)
            if fn is None:
                continue
            try:
                fn(tgt, force=force)
            except Exception as e:
                log.error("[%s/%s] build %s FAILED: %s",
                          tgt.species, tgt.assembly, r, e)
