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


def extract_transcripts(genome_fa: Path, annotation: Path, dst_fa: Path) -> Path:
    """Build a spliced-transcripts FASTA from genome + annotation with
    GENCODE-style headers.

    The output mirrors GENCODE ``*.transcripts.fa`` semantics: per transcript,
    its exons are taken from ``annotation`` (GTF or GFF3), splice-joined in
    transcribed (5'->3') order, reverse-complemented for minus-strand
    transcripts, and emitted under a pipe-delimited header::

        >{transcript_id}|{gene_id}|{strand}|{other_id}|{transcript_name}|{gene_name}|{length}|{transcript_type}|

    Missing optional fields are emitted as ``-``.

    ``genome_fa`` may be plain ``.fa`` or bgzip'd ``.fa.gz`` (.fai required;
    .gzi required for bgzip). ``annotation`` may be GTF/GFF3, plain or gzipped.
    """
    import gzip as _gzip
    import re as _re
    try:
        import pysam  # type: ignore
    except ImportError as e:  # pragma: no cover
        raise RuntimeError("pysam required for extract_transcripts") from e

    dst_fa.parent.mkdir(parents=True, exist_ok=True)

    # ---- 1) parse annotation: collect per-transcript exons + attributes ----
    is_gff3 = str(annotation).rstrip(".gz").endswith((".gff3", ".gff"))
    opener = _gzip.open if str(annotation).endswith(".gz") else open

    # GTF attr:  key "value"; key "value";
    # GFF3 attr: key=value;key=value
    _kv_gtf = _re.compile(r'(\w+)\s+"([^"]*)"')

    def _parse_attrs(s: str) -> dict:
        if is_gff3:
            out = {}
            for chunk in s.rstrip(";").split(";"):
                if "=" in chunk:
                    k, _, v = chunk.partition("=")
                    out[k.strip()] = v.strip()
            return out
        return {m.group(1): m.group(2) for m in _kv_gtf.finditer(s)}

    # transcripts[tid] = {"chrom","strand","exons":[(s,e)],"gene_id","gene_name",
    #                     "transcript_name","biotype","other_id"}
    transcripts: dict[str, dict] = {}
    order: list[str] = []
    with opener(annotation, "rt") as fh:  # type: ignore
        for line in fh:
            if not line or line[0] == "#":
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            ftype = f[2]
            attrs = _parse_attrs(f[8])
            if ftype == "transcript" or (is_gff3 and ftype in ("mRNA", "lnc_RNA",
                    "transcript", "ncRNA", "rRNA", "tRNA", "snoRNA", "snRNA",
                    "miRNA", "primary_transcript", "pseudogenic_transcript")):
                tid = attrs.get("transcript_id") or attrs.get("ID")
                if not tid:
                    continue
                rec = transcripts.setdefault(tid, {
                    "chrom": f[0], "strand": f[6], "exons": [],
                    "gene_id": "-", "gene_name": "-", "transcript_name": "-",
                    "biotype": "-", "other_id": "-",
                })
                if tid not in order:
                    order.append(tid)
                rec["chrom"] = f[0]
                rec["strand"] = f[6] or rec["strand"]
                rec["gene_id"] = (attrs.get("gene_id") or attrs.get("Parent")
                                  or rec["gene_id"])
                rec["gene_name"] = (attrs.get("gene_name") or attrs.get("Name")
                                    or rec["gene_name"])
                rec["transcript_name"] = (attrs.get("transcript_name")
                                          or rec["transcript_name"])
                rec["biotype"] = (attrs.get("transcript_biotype")
                                  or attrs.get("biotype")
                                  or attrs.get("transcript_type")
                                  or rec["biotype"])
                rec["other_id"] = (attrs.get("havana_transcript")
                                   or attrs.get("ccdsid")
                                   or rec["other_id"])
            elif ftype == "exon":
                tid = attrs.get("transcript_id") or attrs.get("Parent")
                if not tid:
                    continue
                # GFF3 Parent may list multiple ids
                for one in tid.split(","):
                    rec = transcripts.setdefault(one, {
                        "chrom": f[0], "strand": f[6], "exons": [],
                        "gene_id": "-", "gene_name": "-", "transcript_name": "-",
                        "biotype": "-", "other_id": "-",
                    })
                    if one not in order:
                        order.append(one)
                    rec["chrom"] = f[0]
                    rec["strand"] = f[6] or rec["strand"]
                    try:
                        rec["exons"].append((int(f[3]), int(f[4])))
                    except ValueError:
                        pass
                    # propagate gene_id if exon line carries it (GTF)
                    if rec["gene_id"] == "-" and attrs.get("gene_id"):
                        rec["gene_id"] = attrs["gene_id"]

    # ---- 2) splice + write fasta ----
    _COMP = str.maketrans("ACGTNacgtnRYKMSWBDHVryksmwbdhv",
                          "TGCANtgcanYRMKSWVHDByrmkswvhdb")

    def _revcomp(s: str) -> str:
        return s.translate(_COMP)[::-1]

    fa = pysam.FastaFile(str(genome_fa))
    fa_refs = set(fa.references)
    n_written = 0
    n_skipped = 0
    with open(dst_fa, "w") as fout:
        for tid in order:
            rec = transcripts.get(tid)
            if not rec or not rec["exons"]:
                continue
            chrom = rec["chrom"]
            if chrom not in fa_refs:
                # try common chr-prefix toggling
                alt = chrom[3:] if chrom.startswith("chr") else "chr" + chrom
                if alt in fa_refs:
                    chrom = alt
                else:
                    n_skipped += 1
                    continue
            exons = sorted(rec["exons"], key=lambda x: x[0])
            seq_parts = []
            try:
                for s, e in exons:
                    seq_parts.append(fa.fetch(chrom, s - 1, e))
            except Exception:
                n_skipped += 1
                continue
            seq = "".join(seq_parts).upper()
            if rec["strand"] == "-":
                seq = _revcomp(seq)
            header = "|".join([
                tid,
                rec["gene_id"],
                rec["strand"] or "-",
                rec["other_id"],
                rec["transcript_name"],
                rec["gene_name"],
                str(len(seq)),
                rec["biotype"],
                "",  # trailing pipe to match GENCODE style
            ])
            fout.write(f">{header}\n")
            for i in range(0, len(seq), 60):
                fout.write(seq[i:i + 60] + "\n")
            n_written += 1
    fa.close()
    log.info("extract_transcripts: wrote %d transcripts (skipped %d) -> %s",
             n_written, n_skipped, dst_fa)
    return dst_fa


def bed_to_bigbed(input_bed: Path, output_bb: Path, chrom_sizes: Path) -> Path:
    """Convert a sorted BED to bigBed using UCSC ``bedToBigBed``.

    Auto-detects column count to pick the right ``-type=bedN`` argument.
    """
    require_tool("bedToBigBed")
    output_bb.parent.mkdir(parents=True, exist_ok=True)
    n_cols = 3
    with open(input_bed) as f:
        for line in f:
            if line.startswith(("#", "track", "browser")) or not line.strip():
                continue
            n_cols = len(line.rstrip("\n").split("\t"))
            break
    n_cols = min(max(n_cols, 3), 12)
    run([
        "bedToBigBed", f"-type=bed{n_cols}",
        str(input_bed), str(chrom_sizes), str(output_bb),
    ])
    return output_bb
