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
import shutil
import subprocess
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


# ── UCSC rmsk.txt.gz → BED6 / GTF converters ───────────────────────────────────
# rmsk columns (17): bin, swScore, milliDiv, milliDel, milliIns,
#   genoName, genoStart, genoEnd, genoLeft, strand,
#   repName, repClass, repFamily, repStart, repEnd, repLeft, id
# genoStart is already 0-based half-open (UCSC BED convention).

def _rmsk_iter(src: Path):
    with open(src) as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            c = line.rstrip("\n").split("\t")
            if len(c) < 13:
                continue
            yield c


def _rmsk_to_bed(src: Path, dst: Path) -> Path:
    """Convert rmsk.txt to BED6: chrom, start, end, repName, swScore, strand."""
    dst.parent.mkdir(parents=True, exist_ok=True)
    n = 0
    with open(dst, "w") as out:
        for c in _rmsk_iter(src):
            chrom, start, end, strand = c[5], c[6], c[7], c[9]
            name, score = c[10], c[1]
            out.write(f"{chrom}\t{start}\t{end}\t{name}\t{score}\t{strand}\n")
            n += 1
    log.info("wrote %d rows -> %s", n, dst.name)
    return dst


def _rmsk_to_gtf(src: Path, dst: Path) -> Path:
    """Convert rmsk.txt to GTF (one 'exon' feature per repeat)."""
    dst.parent.mkdir(parents=True, exist_ok=True)
    n = 0
    with open(dst, "w") as out:
        for c in _rmsk_iter(src):
            chrom = c[5]
            start = int(c[6]) + 1   # BED 0-based -> GTF 1-based
            end = c[7]
            strand = c[9] if c[9] in ("+", "-") else "."
            score = c[1]
            name, rclass, rfamily = c[10], c[11], c[12]
            attrs = (
                f'gene_id "{name}"; transcript_id "{name}"; '
                f'class "{rclass}"; family "{rfamily}";'
            )
            out.write(
                f"{chrom}\trmsk\texon\t{start}\t{end}\t{score}\t{strand}\t.\t{attrs}\n"
            )
            n += 1
    log.info("wrote %d rows -> %s", n, dst.name)
    return dst


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
    out = target.build_dir / "transcripts.fa.gz"
    if _exists(out) and _exists(Path(f"{out}.fai")) and not force:
        log.info("[%s/%s] transcripts up-to-date", target.species, target.assembly)
        return
    # If no raw transcripts FASTA, try to derive it from genome + GTF/GFF using gffread.
    if not src.exists():
        derived = _derive_transcripts(target)
        if derived is None:
            log.info("[%s/%s] no transcripts source (raw, url, or genome+gtf) — skip",
                     target.species, target.assembly)
            return
        src = derived
    bgzip_file(src, out, force=force)
    faidx(out, force=force)


def _derive_transcripts(target: Target) -> Path | None:
    """Build a transcripts FASTA from the genome + annotation when no direct
    transcriptome source is configured. Returns the produced raw .fa or None.
    """
    from .utils import extract_transcripts
    genome = raw_path(target, "genome")
    if not genome.exists():
        return None
    gtf = raw_path(target, "annotation_gtf")
    gff = raw_path(target, "annotation_gff3")
    annot = gtf if gtf.exists() else (gff if gff.exists() else None)
    if annot is None:
        return None
    dst = target.raw_dir / "transcriptome.fa"
    if dst.exists() and dst.stat().st_size > 0:
        return dst
    log.info("[%s/%s] deriving transcriptome via gffread from %s + %s",
             target.species, target.assembly, genome.name, annot.name)
    try:
        extract_transcripts(genome, annot, dst)
    except Exception as e:
        log.error("gffread failed: %s", e)
        return None
    return dst if dst.exists() and dst.stat().st_size > 0 else None


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
        # fall back to deriving from UCSC rmsk.txt.gz when available
        rmsk = raw_path(target, "repeats_rmsk")
        if rmsk.exists():
            log.info("[%s/%s] deriving repeats GTF from rmsk",
                     target.species, target.assembly)
            src = target.raw_dir / "repeats_from_rmsk.gtf"
            if not _exists(src) or force:
                _rmsk_to_gtf(rmsk, src)
        else:
            return
    _build_sorted_gff(src, target.build_dir / "repeats.sorted.gtf.gz", force=force)


def build_repeats_bed(target: Target, *, force: bool = False) -> None:
    src = raw_path(target, "repeats_bed")
    if not src.exists():
        rmsk = raw_path(target, "repeats_rmsk")
        if rmsk.exists():
            log.info("[%s/%s] deriving repeats BED from rmsk",
                     target.species, target.assembly)
            src = target.raw_dir / "repeats_from_rmsk.bed"
            if not _exists(src) or force:
                _rmsk_to_bed(rmsk, src)
        else:
            return
    _build_sorted_bed(src, target.build_dir / "repeats.sorted.bed.gz", force=force)


def _ensembl_to_ucsc_chrom(name: str) -> str:
    """Map a single Ensembl-style chrom to UCSC convention (chr1, chrX, chrM)."""
    if name.startswith("chr"):
        return name
    if name == "MT":
        return "chrM"
    return f"chr{name}"


def _normalize_gff3_to_ucsc(src: Path, dst: Path) -> int:
    """Rewrite a GFF3 file so chromosome names follow the UCSC ``chr*`` convention.

    Used as a pre-step before running UCSC ``liftOver`` (chain files use chr*).
    Comments are preserved; all data lines are rewritten. Returns row count.
    """
    dst.parent.mkdir(parents=True, exist_ok=True)
    n = 0
    with open(src) as fin, open(dst, "w") as fout:
        for line in fin:
            if not line.strip() or line.startswith("#"):
                fout.write(line)
                continue
            chrom, _, rest = line.partition("\t")
            fout.write(f"{_ensembl_to_ucsc_chrom(chrom)}\t{rest}")
            n += 1
    return n


def _liftover_gff(src: Path, chain: Path, dst: Path, unmapped: Path) -> None:
    """Run UCSC liftOver in -gff mode. Requires the ``liftOver`` binary on PATH."""
    if not shutil.which("liftOver"):
        raise RuntimeError(
            "liftOver binary not found on PATH; install UCSC tools or skip "
            "the liftover_from feature."
        )
    cmd = ["liftOver", "-gff", str(src), str(chain), str(dst), str(unmapped)]
    log.info("$ %s", " ".join(cmd))
    # liftOver writes progress + summary to stderr; let it stream.
    subprocess.run(cmd, check=True)


def _load_genome_chroms(target: Target) -> set[str]:
    """Return chromosome names present in the built genome.fa.gz.fai (or empty)."""
    fai = target.build_dir / "genome.fa.gz.fai"
    if not fai.exists():
        return set()
    chroms: set[str] = set()
    with open(fai) as fh:
        for line in fh:
            name = line.split("\t", 1)[0].strip()
            if name:
                chroms.add(name)
    return chroms


def _normalize_rnacentral_chroms(src: Path, dst: Path, genome_chroms: set[str]) -> int:
    """Rewrite rnacentral GFF3 so chromosome names match the target genome.

    RNAcentral uses Ensembl-style names (``1``, ``X``, ``MT``) while UCSC
    genomes use ``chr1``/``chrX``/``chrM``. Records whose chromosome cannot be
    mapped to a name present in the genome are dropped.

    Returns the number of records written (excluding comments/empty lines).
    """
    dst.parent.mkdir(parents=True, exist_ok=True)
    use_chr_prefix = any(c.startswith("chr") for c in genome_chroms)
    kept = 0
    dropped: dict[str, int] = {}
    with open(src) as fin, open(dst, "w") as fout:
        for line in fin:
            if not line.strip() or line.startswith("#"):
                fout.write(line)
                continue
            chrom, _, rest = line.partition("\t")
            mapped = chrom
            if mapped not in genome_chroms and use_chr_prefix:
                # MT -> chrM, otherwise prepend "chr"
                candidate = "chrM" if chrom == "MT" else f"chr{chrom}"
                if candidate in genome_chroms:
                    mapped = candidate
            if mapped not in genome_chroms:
                dropped[chrom] = dropped.get(chrom, 0) + 1
                continue
            fout.write(f"{mapped}\t{rest}")
            kept += 1
    if dropped:
        top = sorted(dropped.items(), key=lambda kv: -kv[1])[:5]
        log.warning(
            "rnacentral: dropped %d records on unmapped chromosomes (top: %s)",
            sum(dropped.values()),
            ", ".join(f"{c}={n}" for c, n in top),
        )
    return kept


def build_rnacentral(target: Target, *, force: bool = False) -> None:
    src = raw_path(target, "rnacentral")
    # If the canonical raw file is missing, attempt to derive it from a
    # source-assembly file via liftOver (configured by `rnacentral.liftover_from`).
    if not src.exists():
        src = _materialize_rnacentral_via_liftover(target, force=force)
        if src is None or not src.exists():
            return
    # raw filename is rnacentral.gff3; some sources actually ship gtf-like.
    # We treat extension generically with tabix -p gff.
    out_gz = target.build_dir / "rnacentral.sorted.gff3.gz"
    if _exists(out_gz) and _exists(Path(f"{out_gz}.tbi")) and not force:
        return
    # Normalize chromosome names to match the built genome.
    genome_chroms = _load_genome_chroms(target)
    if genome_chroms:
        normalized = target.raw_dir / "rnacentral.normalized.gff3"
        n = _normalize_rnacentral_chroms(src, normalized, genome_chroms)
        log.info("[%s/%s] rnacentral: %d records after chrom normalization",
                 target.species, target.assembly, n)
        src = normalized
    _build_sorted_gff(src, out_gz, force=force)


def _materialize_rnacentral_via_liftover(
    target: Target, *, force: bool = False
) -> Path | None:
    """Lift over rnacentral coords from a source assembly into target coords.

    Inputs (created by download step when `liftover_from` is configured):
      raw/rnacentral.source.gff3   -- GFF3 in the source assembly's coordinates
      raw/rnacentral.chain         -- UCSC liftOver chain (source -> target)

    Output:
      raw/rnacentral.gff3          -- GFF3 in this target's coordinates
      raw/rnacentral.unmapped.gff3 -- features that failed to lift (for inspection)
    """
    source_gff = target.raw_dir / "rnacentral.source.gff3"
    chain = target.raw_dir / "rnacentral.chain"
    if not (source_gff.exists() and chain.exists()):
        return None

    out = target.raw_dir / "rnacentral.gff3"
    unmapped = target.raw_dir / "rnacentral.unmapped.gff3"
    if out.exists() and not force:
        log.info("[%s/%s] rnacentral lifted file already exists: %s",
                 target.species, target.assembly, out)
        return out

    # Step 1: ensure chrom names in source match the chain (UCSC chr*).
    chrnorm = target.raw_dir / "rnacentral.source.chrnorm.gff3"
    n = _normalize_gff3_to_ucsc(source_gff, chrnorm)
    log.info("[%s/%s] rnacentral source -> UCSC chrom names: %d rows",
             target.species, target.assembly, n)

    # Step 2: liftOver in GFF mode.
    log.info("[%s/%s] lifting rnacentral via %s",
             target.species, target.assembly, chain.name)
    _liftover_gff(chrnorm, chain, out, unmapped)

    # Report mapping stats.
    try:
        kept = sum(1 for ln in open(out) if ln and not ln.startswith("#") and ln.strip())
        dropped = sum(1 for ln in open(unmapped) if ln and not ln.startswith("#") and ln.strip())
        log.info("[%s/%s] liftOver result: %d mapped, %d unmapped",
                 target.species, target.assembly, kept, dropped)
    except OSError:
        pass
    chrnorm.unlink(missing_ok=True)
    return out


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
    auto_download: bool = True,
    extra_targets: list[Target] | None = None,
) -> None:
    """Build indexed outputs for the requested targets.

    When ``auto_download`` is True (default), any selected resource whose raw
    input is missing will be fetched first by reusing the download pipeline.
    ``extra_targets`` injects ad-hoc Targets (e.g. from ``refbox import``)
    alongside whatever ``species.yaml`` yields.
    """
    resources = resources or RESOURCE_NAMES
    targets = list(iter_targets(species=species, assembly=assembly, out_root=out))
    if extra_targets:
        targets.extend(extra_targets)
    for tgt in targets:
        log.info("=== build %s / %s ===", tgt.species, tgt.assembly)
        tgt.build_dir.mkdir(parents=True, exist_ok=True)

        if auto_download:
            _ensure_raw_files(tgt, resources, force=False)

        for r in resources:
            fn = BUILDERS.get(r)
            if fn is None:
                continue
            try:
                fn(tgt, force=force)
            except Exception as e:
                log.error("[%s/%s] build %s FAILED: %s",
                          tgt.species, tgt.assembly, r, e)


def _ensure_raw_files(target: Target, resources: list[str], *, force: bool) -> None:
    """Download any missing raw files for ``target`` before the build step.

    Custom targets (those whose resources came from a directory ingest rather
    than from species.yaml) carry empty/None spec dicts; for those we just
    verify the raw file is on disk.
    """
    from .download import fetch_resource  # local import to avoid cycle
    for r in resources:
        if target.resource(r) is None:
            continue
        dst = raw_path(target, r)
        if dst.exists() and not force:
            continue
        try:
            fetch_resource(target, r, force=force)
        except Exception as e:
            log.error("[%s/%s] auto-download %s FAILED: %s",
                      target.species, target.assembly, r, e)
