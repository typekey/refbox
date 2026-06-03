"""Single-file build commands invoked by ``refbox build`` (no species config).

Each public ``build_*`` function takes user-supplied paths, validates them
(magic bytes for gzip detection, optional sort check), runs the canonical
transformation, and emits a stable set of output artifacts next to ``-o``.

These helpers exist so that users can bring arbitrary files (custom organism
data, in-house assemblies, pre-processed BEDs from a paper) and get the same
indexed outputs that the species.yaml-driven pipeline produces.
"""

from __future__ import annotations

import logging
import shutil
import subprocess
from pathlib import Path

from .utils import (
    bed_to_bigbed, bgzip_file, extract_transcripts, faidx,
    sort_bed, sort_gff, tabix, write_chrom_sizes,
)

log = logging.getLogger(__name__)


# ── shared helpers ────────────────────────────────────────────────────────────

def _open_text(path: Path):
    """Open path as text, transparently handling .gz."""
    if str(path).endswith(".gz"):
        import gzip
        return gzip.open(path, "rt")
    return open(path, "r")


def _is_sorted_bed(path: Path) -> bool:
    prev = None
    with _open_text(path) as f:
        for line in f:
            if not line or line.startswith(("#", "track", "browser")):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            try:
                key = (parts[0], int(parts[1]))
            except ValueError:
                return False
            if prev is not None and key < prev:
                return False
            prev = key
    return True


def _is_sorted_gff(path: Path) -> bool:
    prev = None
    with _open_text(path) as f:
        for line in f:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 5:
                continue
            try:
                key = (parts[0], int(parts[3]))
            except ValueError:
                return False
            if prev is not None and key < prev:
                return False
            prev = key
    return True


def _materialize_plain(src: Path, dst: Path) -> Path:
    """Copy ``src`` to ``dst`` decompressing if gzipped."""
    dst.parent.mkdir(parents=True, exist_ok=True)
    if str(src).endswith(".gz"):
        with _open_text(src) as fin, open(dst, "w") as fout:
            shutil.copyfileobj(fin, fout, length=1 << 20)
    elif src.resolve() != dst.resolve():
        shutil.copyfile(src, dst)
    return dst


# ── FASTA (genome) ────────────────────────────────────────────────────────────

def build_fa(input_fa: Path, output: Path | None = None) -> Path:
    """bgzip ``input_fa`` and build ``.fai``/``.gzi`` + chrom.sizes.

    ``output`` defaults to ``<input>.bgz.fa.gz`` next to the input. Whatever
    you pass is treated as the desired ``.fa.gz`` path; ``.fai``/``.gzi`` and
    ``<output>.chrom.sizes`` are written alongside.
    """
    if output is None:
        base = input_fa.name
        if base.endswith(".gz"):
            base = base[:-3]
        if not base.endswith((".fa", ".fasta", ".fna")):
            base += ".fa"
        output = input_fa.with_name(base + ".gz")
    output = Path(output)
    if not str(output).endswith(".gz"):
        output = output.with_suffix(output.suffix + ".gz")
    log.info("build_fa: %s -> %s", input_fa, output)
    bgzip_file(input_fa, output, force=True)
    fai = faidx(output, force=True)
    write_chrom_sizes(fai, output.parent / (output.stem + ".chrom.sizes"))
    return output


# ── transcriptome (genome + GTF/GFF -> transcripts.fa.gz) ─────────────────────

def build_transcriptome(
    genome_fa: Path, annotation: Path, output: Path | None = None,
) -> Path:
    """Extract a transcripts FASTA from ``genome_fa`` + ``annotation`` and
    produce the indexed ``transcriptome.fa.gz`` + ``.fai``/``.gzi``.
    """
    if output is None:
        output = annotation.with_name("transcriptome.fa.gz")
    output = Path(output)
    if not str(output).endswith(".gz"):
        output = output.with_suffix(output.suffix + ".gz")
    output.parent.mkdir(parents=True, exist_ok=True)
    log.info("build_transcriptome: %s + %s -> %s", genome_fa, annotation, output)
    # gffread needs plain (or bgzipped) fasta with .fai; ensure that.
    g = Path(genome_fa)
    if not str(g).endswith(".gz"):
        # leave as-is; gffread accepts plain fasta
        pass
    tmp_fa = output.with_suffix(".tmp.fa")
    extract_transcripts(g, Path(annotation), tmp_fa)
    bgzip_file(tmp_fa, output, force=True)
    faidx(output, force=True)
    tmp_fa.unlink(missing_ok=True)
    return output


# ── annotation (GTF / GFF3) ───────────────────────────────────────────────────

def build_gxf(
    input_gxf: Path,
    output: Path | None = None,
    *,
    sqlite: bool = False,
    sqlite_out: Path | None = None,
    source_name: str = "",
    species: str = "",
    genome: str = "",
    annotation_version: str = "",
    synonyms: str | Path | None = None,
    rnacentral: str | Path | None = None,
    fuzzy_scope: str = "names",
    force: bool = False,
) -> Path:
    """Sort + bgzip + tabix a GTF/GFF3 file. Returns the ``.gz`` path.

    When ``sqlite=True`` a read-only RBrowser Index (``.rba``) is also built next
    to the output (``<stem>.rba`` unless ``sqlite_out`` is given). The index
    powers the browser's transcript/gene search; the tabix file still serves
    positional range queries.
    """
    ext = ".gtf"
    lower = input_gxf.name.lower().rstrip(".gz")
    if lower.endswith((".gff3", ".gff")):
        ext = ".gff3"
    if output is None:
        stem = input_gxf.name
        if stem.endswith(".gz"):
            stem = stem[:-3]
        for s in (".gtf", ".gff3", ".gff"):
            if stem.lower().endswith(s):
                stem = stem[: -len(s)]
                break
        output = input_gxf.with_name(stem + ".sorted" + ext + ".gz")
    output = Path(output)
    if not str(output).endswith(".gz"):
        output = output.with_suffix(output.suffix + ".gz")
    output.parent.mkdir(parents=True, exist_ok=True)
    log.info("build_gxf: %s -> %s", input_gxf, output)

    # materialize plain (decompress if needed)
    plain = output.with_suffix(".tmp.plain")
    _materialize_plain(input_gxf, plain)

    if _is_sorted_gff(plain):
        sorted_path = plain
    else:
        log.info("input not sorted; sorting")
        sorted_path = output.with_suffix(".tmp.sorted")
        sort_gff(plain, sorted_path)

    bgzip_file(sorted_path, output, force=True)
    tabix(output, preset="gff", force=True)

    if sqlite:
        from .sqlite_index import build_sqlite_index
        if sqlite_out is None:
            base = output.name
            for s in (".gz", ".gtf", ".gff3", ".gff", ".sorted"):
                if base.endswith(s):
                    base = base[: -len(s)]
            sqlite_out = output.with_name(base + ".rba")
        # Index from the original input (full attribute set), not the sorted
        # plain temp which has already been removed below.
        build_sqlite_index(
            input_gxf, sqlite_out, source_name=source_name, species=species,
            genome=genome, annotation_version=annotation_version,
            synonyms=synonyms, rnacentral=rnacentral, fuzzy_scope=fuzzy_scope,
            force=force,
        )

    plain.unlink(missing_ok=True)
    if sorted_path != plain:
        sorted_path.unlink(missing_ok=True)
    return output


# ── BED ───────────────────────────────────────────────────────────────────────

def build_bed(
    input_bed: Path,
    output: Path | None = None,
    *,
    chrom_sizes: Path | None = None,
    assembly: str | None = None,
    make_bigbed: bool = True,
) -> Path:
    """Sort + bgzip + tabix a BED; also build bigBed when chrom_sizes are known.

    chrom_sizes resolution order:
      1. ``chrom_sizes`` argument (path to file)
      2. ``assembly`` argument -> uses ``zlbio.biofile.convert_bed_to_bigbed``
         which looks up the chrom.sizes from zlbio's config.
    """
    if output is None:
        stem = input_bed.name
        if stem.endswith(".gz"):
            stem = stem[:-3]
        if stem.lower().endswith(".bed"):
            stem = stem[:-4]
        output = input_bed.with_name(stem + ".sorted.bed.gz")
    output = Path(output)
    if not str(output).endswith(".gz"):
        output = output.with_suffix(output.suffix + ".gz")
    output.parent.mkdir(parents=True, exist_ok=True)
    log.info("build_bed: %s -> %s", input_bed, output)

    plain = output.with_suffix(".tmp.plain.bed")
    _materialize_plain(input_bed, plain)

    if _is_sorted_bed(plain):
        sorted_path = plain
    else:
        log.info("input not sorted; sorting")
        sorted_path = output.with_suffix(".tmp.sorted.bed")
        sort_bed(plain, sorted_path)

    bgzip_file(sorted_path, output, force=True)
    tabix(output, preset="bed", force=True)

    bb_path = output.with_suffix("").with_suffix(".bigBed")
    if make_bigbed:
        try:
            if chrom_sizes is not None:
                bed_to_bigbed(sorted_path, bb_path, Path(chrom_sizes))
            elif assembly is not None:
                # Delegate to zlbio which knows per-species chrom.sizes paths.
                from zlbio.biofile import convert_bed_to_bigbed  # type: ignore
                convert_bed_to_bigbed(
                    input_bed=str(sorted_path),
                    output_bigbed=str(bb_path),
                    species=assembly,
                )
            else:
                log.warning("no --chrom-sizes / --assembly given; skip bigBed")
        except Exception as e:
            log.error("bigBed build failed: %s", e)

    plain.unlink(missing_ok=True)
    if sorted_path != plain:
        sorted_path.unlink(missing_ok=True)
    return output


# ── repeats (rmsk.txt[.gz] -> bed + gtf, sorted/bgzipped/indexed) ─────────────

def build_rmsk(input_rmsk: Path, out_dir: Path | None = None) -> dict[str, Path]:
    """Convert a UCSC ``rmsk.txt[.gz]`` table to ``repeats.sorted.bed.gz`` and
    ``repeats.sorted.gtf.gz`` (both bgzipped + tabix-indexed).

    The UCSC rmsk table has these columns (positions used here):
      genoName(5) genoStart(6) genoEnd(7) strand(9) repName(10) repClass(11)
      repFamily(12) swScore(1)
    """
    if out_dir is None:
        out_dir = input_rmsk.parent
    out_dir = Path(out_dir); out_dir.mkdir(parents=True, exist_ok=True)
    log.info("build_rmsk: %s -> %s", input_rmsk, out_dir)

    bed_plain = out_dir / "repeats.tmp.bed"
    gtf_plain = out_dir / "repeats.tmp.gtf"
    with _open_text(input_rmsk) as fin, \
         open(bed_plain, "w") as bed_out, \
         open(gtf_plain, "w") as gtf_out:
        for line in fin:
            if not line.strip() or line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 13:
                continue
            try:
                chrom = f[5]
                start = int(f[6])
                end = int(f[7])
                strand = f[9]
                name = f[10]
                rclass = f[11]
                rfamily = f[12]
                score = int(float(f[1])) if f[1] else 0
            except ValueError:
                continue
            bed_out.write(f"{chrom}\t{start}\t{end}\t{name}|{rclass}|{rfamily}\t{score}\t{strand}\n")
            attrs = (f'gene_id "{name}"; transcript_id "{name}"; '
                     f'class "{rclass}"; family "{rfamily}";')
            gtf_out.write(f"{chrom}\trmsk\texon\t{start+1}\t{end}\t{score}\t{strand}\t.\t{attrs}\n")

    bed_sorted = out_dir / "repeats.sorted.bed"
    gtf_sorted = out_dir / "repeats.sorted.gtf"
    sort_bed(bed_plain, bed_sorted)
    sort_gff(gtf_plain, gtf_sorted)

    bed_gz = out_dir / "repeats.sorted.bed.gz"
    gtf_gz = out_dir / "repeats.sorted.gtf.gz"
    bgzip_file(bed_sorted, bed_gz, force=True)
    bgzip_file(gtf_sorted, gtf_gz, force=True)
    tabix(bed_gz, preset="bed", force=True)
    tabix(gtf_gz, preset="gff", force=True)

    for p in (bed_plain, gtf_plain, bed_sorted, gtf_sorted):
        p.unlink(missing_ok=True)
    return {"bed": bed_gz, "gtf": gtf_gz}


# ── auto-detect ───────────────────────────────────────────────────────────────

_GFF_EXTS = (".gtf", ".gff3", ".gff")
_FA_EXTS = (".fa", ".fasta", ".fna")


def auto_detect(input_path: Path) -> str:
    """Return the canonical builder kind for ``input_path`` based on extension."""
    name = input_path.name.lower()
    base = name[:-3] if name.endswith(".gz") else name
    if "rmsk" in base or base.endswith(".out"):
        if base.endswith(".out"):
            return "fa_out"
        return "rmsk"
    if base.endswith(_FA_EXTS):
        return "fa"
    if base.endswith(_GFF_EXTS):
        return "gxf"
    if base.endswith(".bed"):
        return "bed"
    raise ValueError(
        f"could not auto-detect file type for {input_path.name}; "
        "use -fa / -gtf / -gff / -bed / -rmsk explicitly"
    )
