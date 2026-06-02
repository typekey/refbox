"""Build a static, read-only SQLite search index from a GTF/GFF3 annotation.

The resulting ``*.sqlite`` file is a self-contained, backend-free search index
for a web genome / RNA browser. It is meant to be hosted as a static file and
queried directly from the browser via SQLite WASM + an HTTP Range VFS, so the
schema is optimised for *search* (exact / prefix / fuzzy / alias), not for the
positional range queries that ``tabix`` already serves.

Pipeline
--------
1. Stream-parse the annotation (GTF or GFF3, plain or ``.gz``); group exon /
   CDS / UTR features under their parent transcript and gene.
2. Emit one ``feature`` row per gene and per transcript, plus an ``alias`` row
   for every searchable synonym (names, IDs, versionless Ensembl IDs, RefSeq,
   HAVANA, Dbxref, gene_synonym, …).
3. Populate two FTS5 virtual tables — ``feature_fts`` (prefix / autocomplete)
   and ``feature_trigram`` (substring / fuzzy). If the installed SQLite lacks
   the trigram tokenizer the builder degrades gracefully (LIKE fallback).
4. Build secondary indexes, ``ANALYZE``, ``PRAGMA optimize`` and ``VACUUM`` for
   a compact, query-ready file.

Coordinate convention
----------------------
``start`` / ``end`` and all ``*_start`` / ``*_end`` columns are **1-based,
inclusive** (GTF/GFF convention) for user-facing display. ``chrom_start0`` /
``chrom_end0`` additionally expose the **0-based, half-open** span for browser
rendering. This is recorded in the ``metadata`` table under
``coord_convention``.

Pure standard library (``sqlite3``, ``gzip``, ``re``, ``json``) — no third
party dependency. ``tqdm`` is used only if importable.
"""

from __future__ import annotations

import gzip
import json
import logging
import re
import sqlite3
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable, Iterator

log = logging.getLogger(__name__)

# ── format / feature-type constants ───────────────────────────────────────────

# GFF3 feature types treated as "transcript-like" (each becomes one feature row).
_TRANSCRIPT_TYPES = {
    "mrna", "transcript", "lnc_rna", "ncrna", "rrna", "trna", "snrna",
    "snorna", "mirna", "pseudogenic_transcript", "primary_transcript",
    "scrna", " scarna", "scarna", "guide_rna", "rnase_p_rna", "rnase_mrp_rna",
    "srp_rna", "telomerase_rna", "vault_rna", "y_rna", "antisense_rna",
    "c_gene_segment", "d_gene_segment", "j_gene_segment", "v_gene_segment",
    "three_prime_overlapping_ncrna", "processed_transcript",
}
# Sub-features that define transcript structure.
_EXON_TYPES = {"exon"}
_CDS_TYPES = {"cds"}
_UTR5_TYPES = {"five_prime_utr", "5utr", "five_prime_utr_variant"}
_UTR3_TYPES = {"three_prime_utr", "3utr", "three_prime_utr_variant"}
_UTR_GENERIC_TYPES = {"utr"}
_GENE_TYPES = {
    "gene", "ncrna_gene", "pseudogene", "snrna_gene", "snorna_gene",
    "mirna_gene", "rrna_gene", "trna_gene", "lincrna_gene",
}

# GTF attribute regex:  key "value";   GFF3 is key=value;key=value
_GTF_ATTR = re.compile(r'(\w+)\s+"([^"]*)"')
_ENSEMBL_VERSION = re.compile(r"\.\d+$")
_SEPARATORS = re.compile(r"[\s_.\-]+")

# alias types that are *names* (worth substring/fuzzy matching), as opposed to
# IDs/accessions (gene_id, transcript_id, rnacentral_id/_db, havana, ccds,
# protein_id, hgnc, dbxref, refseq) which users search exactly or by prefix, not
# by fuzzy substring. The trigram (substring) index is built over names only
# when fuzzy_scope='names' — this both shrinks it dramatically and removes the
# pathological "numeric ID substring" query (e.g. 000003351 scanning all IDs).
_NAME_LIKE_ALIAS = {"gene_name", "transcript_name", "gene_synonym", "name", "alias"}


def _fuzzy_text(names: "list[str]", aliases: dict[str, str]) -> str:
    """Substring-search corpus for one feature: its name columns + name-like
    aliases (synonyms / Name / Alias), but NOT IDs."""
    parts = [x for x in names if x]
    parts += [a for a, ty in aliases.items() if ty in _NAME_LIKE_ALIAS]
    return " ".join(parts)


# ── normalization ─────────────────────────────────────────────────────────────

def strip_version(value: str) -> str:
    """``ENST00000335137.4`` → ``ENST00000335137`` (drops a trailing ``.<int>``)."""
    return _ENSEMBL_VERSION.sub("", value)


def normalize(value: str | None) -> str:
    """Normalize a token for fuzzy / separator-insensitive matching.

    lowercase → strip a trailing Ensembl ``.<version>`` → remove the common
    separators ``_ - . space``. Original values are always kept for display;
    this is only used to populate ``alias_norm`` and to canonicalize queries.

        ENSG00000141510.18 → ensg00000141510
        TP53-201           → tp53201
        p53                → p53
    """
    if not value:
        return ""
    s = value.strip().lower()
    s = _ENSEMBL_VERSION.sub("", s)
    s = _SEPARATORS.sub("", s)
    return s


# ── parsed records ────────────────────────────────────────────────────────────

@dataclass
class _Tx:
    """Accumulator for one transcript while streaming the annotation."""
    transcript_id: str
    chrom: str = ""
    strand: str = ""
    start: int = 0          # 1-based inclusive transcript span (from feature line)
    end: int = 0
    gene_id: str = ""
    gene_name: str = ""
    transcript_name: str = ""
    biotype: str = ""
    source: str = ""
    feature_type: str = "transcript"   # original SO term (mRNA, lnc_RNA, …)
    exons: list[tuple[int, int]] = field(default_factory=list)
    cds: list[tuple[int, int]] = field(default_factory=list)
    utr5: list[tuple[int, int]] = field(default_factory=list)
    utr3: list[tuple[int, int]] = field(default_factory=list)
    utr_generic: list[tuple[int, int]] = field(default_factory=list)
    aliases: dict[str, str] = field(default_factory=dict)  # alias -> alias_type


@dataclass
class _Gene:
    gene_id: str
    chrom: str = ""
    strand: str = ""
    start: int = 0
    end: int = 0
    gene_name: str = ""
    biotype: str = ""
    source: str = ""
    feature_type: str = "gene"
    aliases: dict[str, str] = field(default_factory=dict)
    has_gene_feature: bool = False     # True if an explicit gene line was seen


# ── attribute parsing ─────────────────────────────────────────────────────────

def _detect_gff3(path: Path, sample_attr: str | None) -> bool:
    """True if the file is GFF3 (``key=value``) rather than GTF (``key "v"``)."""
    name = path.name.lower().rstrip(".gz")
    if name.endswith((".gff3", ".gff")):
        return True
    if name.endswith(".gtf"):
        return False
    # Fall back to inspecting an attribute column.
    if sample_attr is not None:
        return "=" in sample_attr and '"' not in sample_attr.split(";")[0]
    return False


def _parse_attrs_gtf(col9: str) -> dict[str, str]:
    return {m.group(1): m.group(2) for m in _GTF_ATTR.finditer(col9)}


def _parse_attrs_gff3(col9: str) -> dict[str, str]:
    out: dict[str, str] = {}
    for chunk in col9.rstrip(";").split(";"):
        chunk = chunk.strip()
        if not chunk or "=" not in chunk:
            continue
        k, _, v = chunk.partition("=")
        out[k.strip()] = _gff3_unescape(v.strip())
    return out


def _gff3_unescape(value: str) -> str:
    """Decode the handful of percent-escapes GFF3 mandates in attribute values."""
    if "%" not in value:
        return value
    try:
        from urllib.parse import unquote
        return unquote(value)
    except Exception:
        return value


def _open_text(path: Path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return open(path, "r", encoding="utf-8", errors="replace")


def _maybe_tqdm(total_bytes: int):
    """Return a progress callback ``update(nbytes)`` and a ``close()``; no-op
    when ``tqdm`` is unavailable."""
    try:
        from tqdm import tqdm  # type: ignore
        bar = tqdm(total=total_bytes, unit="B", unit_scale=True,
                   desc="parse", disable=total_bytes <= 0)
        return bar.update, bar.close
    except Exception:
        return (lambda _n: None), (lambda: None)


# ── external synonym enrichment (HGNC etc.) ───────────────────────────────────

def load_synonyms(path: Path) -> tuple[dict[str, set[str]], dict[str, set[str]]]:
    """Load a gene-synonym table → ``(by_ensembl_id, by_symbol_upper)`` maps.

    Designed for the HGNC ``hgnc_complete_set.txt`` TSV (columns ``symbol``,
    ``alias_symbol``, ``prev_symbol``, ``ensembl_gene_id``; multi-value fields
    are ``|``-separated and may be double-quoted), but works with any TSV that
    has a ``symbol`` column plus any of those. Synonyms let common names that
    annotation files omit — ``OCT4`` → ``POU5F1``, ``p53`` → ``TP53`` — resolve
    as exact alias hits.
    """
    by_ensembl: dict[str, set[str]] = {}
    by_symbol: dict[str, set[str]] = {}
    with _open_text(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        idx = {name.strip(): i for i, name in enumerate(header)}
        c_sym = idx.get("symbol")
        c_alias = idx.get("alias_symbol")
        c_prev = idx.get("prev_symbol")
        c_ens = idx.get("ensembl_gene_id")
        if c_sym is None:
            raise ValueError(
                f"synonyms file {path} has no 'symbol' column (header={header[:5]}…)")

        def cell(fields: list[str], i: int | None) -> str:
            if i is None or i >= len(fields):
                return ""
            return fields[i].strip().strip('"')

        for line in fh:
            f = line.rstrip("\n").split("\t")
            sym = cell(f, c_sym)
            syns: set[str] = set()
            for col in (c_alias, c_prev):
                val = cell(f, col)
                if val:
                    for s in val.split("|"):
                        s = s.strip()
                        if s and s != sym:
                            syns.add(s)
            if not syns:
                continue
            ens = cell(f, c_ens)
            if ens:
                by_ensembl[strip_version(ens)] = syns
            if sym:
                by_symbol[sym.upper()] = syns
    return by_ensembl, by_symbol


def enrich_with_synonyms(
    genes: dict[str, "_Gene"], path: Path,
) -> int:
    """Inject external synonyms as ``gene_synonym`` aliases on matching genes.

    Matched first by (versionless) Ensembl gene ID, else by gene symbol. Returns
    the number of synonym aliases added.
    """
    by_ensembl, by_symbol = load_synonyms(path)
    added = 0
    for g in genes.values():
        syns = by_ensembl.get(strip_version(g.gene_id))
        if syns is None and g.gene_name:
            syns = by_symbol.get(g.gene_name.upper())
        if not syns:
            continue
        for s in syns:
            if s not in g.aliases:
                g.aliases[s] = "gene_synonym"
                added += 1
    log.info("synonyms: injected %d aliases from %s", added, path.name)
    return added


# ── RNAcentral merge ──────────────────────────────────────────────────────────

def parse_rnacentral(path: Path, *, verbose: bool = False) -> dict[str, "_Tx"]:
    """Parse an RNAcentral genome-coordinates GFF3 into transcript records.

    RNAcentral GFF3 uses ``transcript`` features (``ID=URS…_9606.N``,
    ``Name=URS…_9606``, ``type=lncRNA|piRNA|…``) with ``noncoding_exon`` children
    (``Parent=URS…_9606.N``). ``predicted_gene`` rows are ignored. Each record
    is searchable by its full ID, its URS accession (with and without the
    ``_taxid``), so e.g. ``URS000035F234`` resolves.

    Use the chromosome-normalized file (``chr``-prefixed) so coordinates display
    consistently with a GENCODE/UCSC main annotation.
    """
    txs: dict[str, _Tx] = {}
    n_lines = 0
    with _open_text(path) as fh:
        for raw in fh:
            n_lines += 1
            if not raw or raw[0] == "#":
                continue
            cols = raw.rstrip("\n").split("\t")
            if len(cols) < 9:
                continue
            ftype = cols[2]
            if ftype == "transcript":
                attrs = _parse_attrs_gff3(cols[8])
                tid = attrs.get("ID")
                if not tid:
                    continue
                try:
                    start, end = int(cols[3]), int(cols[4])
                except ValueError:
                    continue
                name = attrs.get("Name") or tid
                rtype = attrs.get("type") or "ncRNA"
                t = txs.get(tid)
                if t is None:
                    t = txs[tid] = _Tx(transcript_id=tid)
                t.chrom, t.strand = cols[0], cols[6]
                t.start, t.end = start, end
                t.source = "RNAcentral"
                t.feature_type = rtype
                t.biotype = rtype
                t.transcript_name = name
                _add_alias(t.aliases, tid, "rnacentral_id")
                _add_alias(t.aliases, name, "rnacentral_id")
                base = strip_version(tid)          # URS…_9606.1 -> URS…_9606
                if base != tid:
                    _add_alias(t.aliases, base, "rnacentral_id")
                urs = name.split("_")[0] if name else ""   # URS…_9606 -> URS…
                if urs and urs != name:
                    _add_alias(t.aliases, urs, "rnacentral_id_versionless")
                for db in (attrs.get("databases") or "").split(","):
                    _add_alias(t.aliases, db.strip(), "rnacentral_db")
            elif ftype == "noncoding_exon":
                attrs = _parse_attrs_gff3(cols[8])
                parent = attrs.get("Parent")
                if not parent:
                    continue
                try:
                    s, e = int(cols[3]), int(cols[4])
                except ValueError:
                    continue
                for one in parent.split(","):
                    t = txs.get(one)
                    if t is None:
                        t = txs[one] = _Tx(transcript_id=one)
                        t.chrom, t.strand = cols[0], cols[6]
                        t.source = "RNAcentral"
                        _add_alias(t.aliases, one, "rnacentral_id")
                    t.exons.append((s, e))
            # predicted_gene and others: ignored
    log.info("rnacentral: parsed %d lines -> %d transcripts", n_lines, len(txs))
    return txs


def _harmonize_chroms(rna: dict[str, "_Tx"], main_chroms: set[str]) -> None:
    """Rewrite RNAcentral chrom names in-place to match the main annotation's
    ``chr`` convention (add/strip the prefix; ``MT``↔``chrM``). No-op if the two
    sources already agree or the main set is empty."""
    if not main_chroms:
        return
    main_has_chr = sum(c.startswith("chr") for c in main_chroms) > len(main_chroms) / 2
    rna_chroms = {t.chrom for t in rna.values() if t.chrom}
    rna_has_chr = bool(rna_chroms) and \
        sum(c.startswith("chr") for c in rna_chroms) > len(rna_chroms) / 2
    if main_has_chr == rna_has_chr:
        return  # same convention already

    def remap(c: str) -> str:
        if main_has_chr and not c.startswith("chr"):
            return "chrM" if c in ("MT", "mt") else "chr" + c
        if not main_has_chr and c.startswith("chr"):
            return "MT" if c == "chrM" else c[3:]
        return c

    for t in rna.values():
        if t.chrom:
            t.chrom = remap(t.chrom)


# ── streaming parser ──────────────────────────────────────────────────────────

def parse_annotation(
    path: Path, *, verbose: bool = False,
) -> tuple[dict[str, _Gene], dict[str, _Tx], bool]:
    """Stream-parse a GTF/GFF3 file into gene/transcript accumulators.

    Returns ``(genes, transcripts, is_gff3)``. Memory scales with the number of
    genes/transcripts (and their exon counts), not file size — each data line is
    processed and discarded.
    """
    genes: dict[str, _Gene] = {}
    transcripts: dict[str, _Tx] = {}

    # peek one data line to decide the attribute dialect
    is_gff3 = None
    n_lines = 0
    try:
        total_bytes = path.stat().st_size
    except OSError:
        total_bytes = 0
    update, close = _maybe_tqdm(total_bytes if not str(path).endswith(".gz") else 0)

    with _open_text(path) as fh:
        for raw in fh:
            n_lines += 1
            if not str(path).endswith(".gz"):
                update(len(raw))
            if not raw or raw[0] == "#":
                continue
            line = raw.rstrip("\n")
            cols = line.split("\t")
            if len(cols) < 9:
                continue
            if is_gff3 is None:
                is_gff3 = _detect_gff3(path, cols[8])
                log.info("annotation dialect: %s", "GFF3" if is_gff3 else "GTF")
            chrom, source, ftype_raw = cols[0], cols[1], cols[2]
            ftype = ftype_raw.lower()
            try:
                start = int(cols[3])
                end = int(cols[4])
            except ValueError:
                continue
            strand = cols[6]
            attrs = (_parse_attrs_gff3 if is_gff3 else _parse_attrs_gtf)(cols[8])

            if ftype in _GENE_TYPES:
                _ingest_gene(genes, chrom, source, start, end, strand, attrs,
                             ftype_raw)
            elif ftype in _TRANSCRIPT_TYPES:
                _ingest_transcript(transcripts, genes, chrom, source, start, end,
                                   strand, attrs, ftype_raw, is_gff3)
            elif ftype in _EXON_TYPES:
                _ingest_subfeature(transcripts, attrs, start, end, "exon", is_gff3)
            elif ftype in _CDS_TYPES:
                _ingest_subfeature(transcripts, attrs, start, end, "cds", is_gff3)
            elif ftype in _UTR5_TYPES:
                _ingest_subfeature(transcripts, attrs, start, end, "utr5", is_gff3)
            elif ftype in _UTR3_TYPES:
                _ingest_subfeature(transcripts, attrs, start, end, "utr3", is_gff3)
            elif ftype in _UTR_GENERIC_TYPES:
                _ingest_subfeature(transcripts, attrs, start, end, "utr", is_gff3)
            # start_codon / stop_codon / Selenocysteine etc. are ignored: CDS
            # already bounds the coding region.

    close()
    log.info("parsed %d lines → %d genes, %d transcripts",
             n_lines, len(genes), len(transcripts))
    return genes, transcripts, bool(is_gff3)


def _add_alias(bucket: dict[str, str], value: str | None, atype: str) -> None:
    if value:
        v = value.strip()
        if v and v not in bucket:
            bucket[v] = atype


def _collect_common_aliases(bucket: dict[str, str], attrs: dict[str, str]) -> None:
    """Aliases shared by gene & transcript GFF3/GTF attribute conventions."""
    _add_alias(bucket, attrs.get("Name"), "name")
    for raw in (attrs.get("Alias") or "").split(","):
        _add_alias(bucket, raw, "alias")
    for raw in (attrs.get("Dbxref") or attrs.get("Xref") or "").split(","):
        raw = raw.strip()
        if not raw:
            continue
        # Dbxref values look like  RefSeq:NM_000546.6  or  HGNC:HGNC:11998
        low = raw.lower()
        if low.startswith(("refseq:", "ccds:", "ensembl:", "ucsc:")):
            _add_alias(bucket, raw.split(":", 1)[1], "refseq"
                       if low.startswith("refseq") else "dbxref")
        _add_alias(bucket, raw, "dbxref")
    for raw in (attrs.get("gene_synonym") or "").split(","):
        _add_alias(bucket, raw, "gene_synonym")


def _ingest_gene(genes, chrom, source, start, end, strand, attrs, ftype_raw):
    gid = attrs.get("gene_id") or attrs.get("ID")
    if not gid:
        return
    g = genes.get(gid)
    if g is None:
        g = genes[gid] = _Gene(gene_id=gid)
    g.chrom = chrom or g.chrom
    g.strand = strand or g.strand
    g.start = start
    g.end = end
    g.source = source or g.source
    g.feature_type = "gene"
    g.has_gene_feature = True
    g.gene_name = (attrs.get("gene_name") or attrs.get("Name")
                   or attrs.get("gene") or g.gene_name)
    g.biotype = (attrs.get("gene_type") or attrs.get("gene_biotype")
                 or attrs.get("biotype") or g.biotype)
    _add_alias(g.aliases, gid, "gene_id")
    # No separate gene_id_versionless row: normalize() already strips the
    # ``.<version>`` suffix, so the gene_id row's alias_norm IS the versionless
    # form — a versionless query matches it directly.
    _add_alias(g.aliases, g.gene_name, "gene_name")
    _add_alias(g.aliases, attrs.get("havana_gene"), "havana")
    _add_alias(g.aliases, attrs.get("hgnc_id"), "hgnc")
    _collect_common_aliases(g.aliases, attrs)


def _ingest_transcript(transcripts, genes, chrom, source, start, end, strand,
                       attrs, ftype_raw, is_gff3):
    tid = attrs.get("transcript_id") or attrs.get("ID")
    if not tid:
        return
    t = transcripts.get(tid)
    if t is None:
        t = transcripts[tid] = _Tx(transcript_id=tid)
    t.chrom = chrom or t.chrom
    t.strand = strand or t.strand
    t.start = start
    t.end = end
    t.source = source or t.source
    t.feature_type = ftype_raw
    t.gene_id = (attrs.get("gene_id") or attrs.get("Parent") or t.gene_id)
    t.gene_name = (attrs.get("gene_name") or attrs.get("gene")
                   or attrs.get("Name") or t.gene_name)
    t.transcript_name = (attrs.get("transcript_name") or t.transcript_name)
    t.biotype = (attrs.get("transcript_type") or attrs.get("transcript_biotype")
                 or attrs.get("biotype") or t.biotype)
    _add_alias(t.aliases, tid, "transcript_id")
    # (no transcript_id_versionless: alias_norm of transcript_id is already
    #  version-stripped, see _ingest_gene)
    _add_alias(t.aliases, t.transcript_name, "transcript_name")
    _add_alias(t.aliases, attrs.get("havana_transcript"), "havana")
    _add_alias(t.aliases, attrs.get("ccdsid") or attrs.get("ccds_id"), "ccds")
    _add_alias(t.aliases, attrs.get("protein_id"), "protein_id")
    _collect_common_aliases(t.aliases, attrs)
    # ensure the parent gene exists even if no explicit gene line precedes it
    if t.gene_id:
        for one in t.gene_id.split(","):
            g = genes.get(one)
            if g is None:
                g = genes[one] = _Gene(gene_id=one)
                _add_alias(g.aliases, one, "gene_id")
            if not g.gene_name and t.gene_name:
                g.gene_name = t.gene_name
                _add_alias(g.aliases, t.gene_name, "gene_name")
            if not g.chrom:
                g.chrom, g.strand, g.source = chrom, strand, source


def _ingest_subfeature(transcripts, attrs, start, end, kind, is_gff3):
    parent = (attrs.get("transcript_id") or attrs.get("Parent"))
    if not parent:
        return
    for tid in parent.split(","):
        t = transcripts.get(tid)
        if t is None:
            # sub-feature seen before its transcript line (rare ordering)
            t = transcripts[tid] = _Tx(transcript_id=tid)
            _add_alias(t.aliases, tid, "transcript_id")
            if not t.gene_id and attrs.get("gene_id"):
                t.gene_id = attrs["gene_id"]
        getattr(t, "exons" if kind == "exon" else
                "cds" if kind == "cds" else
                "utr5" if kind == "utr5" else
                "utr3" if kind == "utr3" else "utr_generic").append((start, end))


# ── transcript structure finalization ─────────────────────────────────────────

def _span(intervals: list[tuple[int, int]]) -> tuple[int | None, int | None]:
    if not intervals:
        return None, None
    return min(s for s, _ in intervals), max(e for _, e in intervals)


def _resolve_utrs(t: _Tx) -> tuple[tuple[int | None, int | None],
                                   tuple[int | None, int | None]]:
    """Return ((utr5_start, utr5_end), (utr3_start, utr3_end)), 1-based inclusive.

    Priority: explicit five/three_prime_UTR features → split generic ``UTR``
    features by CDS position → infer from transcript span vs CDS boundaries.
    Strand-aware. Returns ``(None, None)`` pairs when nothing can be inferred
    confidently (e.g. non-coding transcripts, or coding ones with no CDS).
    """
    if t.utr5 or t.utr3:
        return _span(t.utr5), _span(t.utr3)

    cds_min, cds_max = _span(t.cds)
    if cds_min is None:
        return (None, None), (None, None)   # no CDS → cannot place UTRs

    plus = t.strand != "-"
    if t.utr_generic:
        low = [iv for iv in t.utr_generic if iv[1] < cds_min]
        high = [iv for iv in t.utr_generic if iv[0] > cds_max]
        five, three = (low, high) if plus else (high, low)
        return _span(five), _span(three)

    # No UTR features at all: infer spans from transcript bounds vs CDS.
    tx_start = t.start or (min(s for s, _ in t.exons) if t.exons else cds_min)
    tx_end = t.end or (max(e for _, e in t.exons) if t.exons else cds_max)
    low = (tx_start, cds_min - 1) if tx_start <= cds_min - 1 else None
    high = (cds_max + 1, tx_end) if cds_max + 1 <= tx_end else None
    five, three = (low, high) if plus else (high, low)
    return (five or (None, None)), (three or (None, None))


# ── SQLite schema ─────────────────────────────────────────────────────────────

_SCHEMA = """
CREATE TABLE feature (
    id            INTEGER PRIMARY KEY,
    feature_type  TEXT,      -- 'gene' or a transcript SO term (mRNA, lnc_RNA…)
    gene_id       TEXT,
    gene_name     TEXT,
    transcript_id TEXT,
    transcript_name TEXT,
    chrom         TEXT,
    start         INTEGER,   -- 1-based inclusive
    end           INTEGER,   -- 1-based inclusive
    chrom_start0  INTEGER,   -- 0-based half-open
    chrom_end0    INTEGER,   -- 0-based half-open
    strand        TEXT,
    biotype       TEXT,
    source        TEXT,
    exon_count    INTEGER,
    exon_starts   TEXT,      -- comma-separated, 1-based, genomic-ascending
    exon_ends     TEXT,
    cds_start     INTEGER,
    cds_end       INTEGER,
    utr5_start    INTEGER,
    utr5_end      INTEGER,
    utr3_start    INTEGER,
    utr3_end      INTEGER,
    search_text   TEXT,      -- always NULL: the text lives in the FTS indexes
                             -- only (kept as a column for schema stability)
    payload_json  TEXT       -- lean: {"aliases":[...]} (other fields are columns)
);

CREATE TABLE alias (
    id          INTEGER PRIMARY KEY,
    feature_id  INTEGER NOT NULL,   -- the gene/transcript this alias points to
    alias_norm  TEXT,               -- normalized form, the only thing matched on
    alias_type  TEXT                -- gene_id | transcript_name | gene_synonym | rnacentral_id | …
    -- (original alias text & source dropped: the readable values live on the
    --  feature row / payload_json, so storing them again was pure overhead)
);

CREATE TABLE metadata (
    key   TEXT PRIMARY KEY,
    value TEXT
);
"""

# feature_fts only runs single-token prefix queries (``"tok"*``), so
# ``detail=none, columnsize=0`` is safe and shrinks it ~84% (measured 153→24 MB
# on GENCODE v44): we never use phrase / proximity / column-filtered / highlight
# queries, which is the only thing detail data is needed for.
_FTS_PREFIX = """
CREATE VIRTUAL TABLE feature_fts USING fts5(
    display_name, gene_id, gene_name, transcript_id, transcript_name,
    aliases, search_text,
    content='', tokenize='unicode61', prefix='2 3 4 5 6 7 8 9 10',
    detail=none, columnsize=0
);
"""

# feature_trigram MUST keep detail=full: a substring query longer than 3 chars
# becomes a *phrase* of consecutive trigram tokens, and FTS5 only supports phrase
# queries with detail=full. (detail=none would silently break e.g. "0000026930".)
# columnsize=0 is still safe (single column) and trims a little.
_FTS_TRIGRAM = """
CREATE VIRTUAL TABLE feature_trigram USING fts5(
    search_text, content='', tokenize='trigram', columnsize=0
);
"""

# The four exact-match columns are indexed COLLATE NOCASE so that
# case-insensitive lookups (``WHERE gene_name = ? COLLATE NOCASE``) use the
# index instead of scanning the whole table.
_INDEXES = [
    "CREATE INDEX idx_feature_gene_id        ON feature(gene_id COLLATE NOCASE)",
    "CREATE INDEX idx_feature_gene_name      ON feature(gene_name COLLATE NOCASE)",
    "CREATE INDEX idx_feature_transcript_id  ON feature(transcript_id COLLATE NOCASE)",
    "CREATE INDEX idx_feature_transcript_name ON feature(transcript_name COLLATE NOCASE)",
    "CREATE INDEX idx_feature_chrom_start_end ON feature(chrom, start, end)",
    # Covering index for the unified exact lookup:
    #   SELECT alias_type, feature_id FROM alias WHERE alias_norm = ?
    # Including alias_type + feature_id makes it index-only (no table row fetch),
    # so one indexed seek resolves any id/name/synonym match across all species
    # — over an HTTP Range VFS that is a handful of page reads instead of ~22.
    # (No idx_alias_feature: nothing queries alias by feature_id — a feature's
    #  own aliases are already inlined in feature.payload_json.)
    "CREATE INDEX idx_alias_norm ON alias(alias_norm, alias_type, feature_id)",
]


def _has_trigram(con: sqlite3.Connection) -> bool:
    try:
        con.execute("CREATE VIRTUAL TABLE _trgm_probe USING fts5(x, tokenize='trigram')")
        con.execute("DROP TABLE _trgm_probe")
        return True
    except sqlite3.OperationalError:
        return False


def _fts5_available(con: sqlite3.Connection) -> bool:
    try:
        con.execute("CREATE VIRTUAL TABLE _fts_probe USING fts5(x)")
        con.execute("DROP TABLE _fts_probe")
        return True
    except sqlite3.OperationalError:
        return False


# ── row materialization ───────────────────────────────────────────────────────

def _gene_row(g: _Gene) -> tuple:
    start = g.start or None
    end = g.end or None
    search_text = " ".join(filter(None, [
        g.gene_id, strip_version(g.gene_id), g.gene_name, g.biotype,
        *g.aliases.keys(),
    ]))
    # Lean payload: only what the typed columns don't already carry (the alias
    # list). Everything else the browser reconstructs from the SELECTed columns,
    # which avoids ~130 MB of column/JSON duplication.
    payload = {"aliases": sorted(g.aliases.keys())}
    return (
        "gene", g.gene_id, g.gene_name, None, None, g.chrom, start, end,
        (start - 1) if start else None, end, g.strand, g.biotype, g.source,
        # exon_count, exon_starts, exon_ends, cds_start/end, utr5/3_start/end (9)
        None, None, None, None, None, None, None, None, None,
        search_text, json.dumps(payload, ensure_ascii=False, sort_keys=True),
    )


def _tx_row(t: _Tx) -> tuple:
    exons = sorted(t.exons)
    exon_starts = ",".join(str(s) for s, _ in exons) if exons else None
    exon_ends = ",".join(str(e) for _, e in exons) if exons else None
    cds_start, cds_end = _span(t.cds)
    (utr5_start, utr5_end), (utr3_start, utr3_end) = _resolve_utrs(t)
    start = t.start or (exons[0][0] if exons else None)
    end = t.end or (exons[-1][1] if exons else None)
    gene_id = t.gene_id.split(",")[0] if t.gene_id else ""
    search_text = " ".join(filter(None, [
        t.transcript_id, strip_version(t.transcript_id), t.transcript_name,
        gene_id, strip_version(gene_id), t.gene_name, t.biotype,
        *t.aliases.keys(),
    ]))
    # Lean payload: only the alias list (not in any column); all structural
    # fields are reconstructed by the browser from the SELECTed columns.
    payload = {"aliases": sorted(t.aliases.keys())}
    return (
        t.feature_type, gene_id, t.gene_name, t.transcript_id, t.transcript_name,
        t.chrom, start, end, (start - 1) if start else None, end, t.strand,
        t.biotype, t.source, len(exons), exon_starts, exon_ends,
        cds_start, cds_end, utr5_start, utr5_end, utr3_start, utr3_end,
        search_text, json.dumps(payload, ensure_ascii=False, sort_keys=True),
    )


_FEATURE_COLS = (
    "feature_type, gene_id, gene_name, transcript_id, transcript_name, chrom, "
    "start, end, chrom_start0, chrom_end0, strand, biotype, source, exon_count, "
    "exon_starts, exon_ends, cds_start, cds_end, utr5_start, utr5_end, "
    "utr3_start, utr3_end, search_text, payload_json"
)


# ── builder ───────────────────────────────────────────────────────────────────

def build_sqlite_index(
    input_path: Path,
    output: Path | None = None,
    *,
    source_name: str = "",
    species: str = "",
    genome: str = "",
    annotation_version: str = "",
    synonyms: Path | str | None = None,
    rnacentral: Path | str | None = None,
    fuzzy_scope: str = "names",
    force: bool = False,
    verbose: bool = False,
) -> Path:
    """Build the read-only SQLite search index. Returns the output path.

    ``output`` defaults to ``<input-stem>.rbrowser.sqlite`` next to the input.
    ``synonyms`` optionally points at an HGNC-style TSV whose alias/prev symbols
    are injected as ``gene_synonym`` aliases (so ``OCT4`` resolves to POU5F1).
    ``rnacentral`` optionally points at an RNAcentral genome-coordinates GFF3
    (use the chromosome-normalized one) whose ncRNA transcripts are merged in as
    additional searchable records.
    ``fuzzy_scope`` controls the trigram (substring) corpus: ``'names'`` (default)
    indexes only gene/transcript names + synonyms — IDs are exact/prefix-only, so
    the trigram is ~6× smaller and the pathological numeric-ID-substring query
    disappears; ``'all'`` indexes the full search_text (IDs included) like before.
    """
    if fuzzy_scope not in ("names", "all"):
        raise ValueError("fuzzy_scope must be 'names' or 'all'")
    input_path = Path(input_path)
    if output is None:
        stem = input_path.name
        if stem.endswith(".gz"):
            stem = stem[:-3]
        for s in (".gtf", ".gff3", ".gff"):
            if stem.lower().endswith(s):
                stem = stem[: -len(s)]
                break
        output = input_path.with_name(stem + ".rbrowser.sqlite")
    output = Path(output)
    if output.exists() and not force:
        raise FileExistsError(
            f"{output} exists; pass force=True / --force to overwrite")
    output.parent.mkdir(parents=True, exist_ok=True)

    t0 = time.time()
    genes, transcripts, is_gff3 = parse_annotation(input_path, verbose=verbose)
    t_parse = time.time() - t0

    n_synonyms = 0
    if synonyms:
        n_synonyms = enrich_with_synonyms(genes, Path(synonyms))

    n_rnacentral = 0
    if rnacentral:
        rna = parse_rnacentral(Path(rnacentral), verbose=verbose)
        # Harmonize chromosome naming to the main annotation's style so ncRNA
        # regions display consistently (RNAcentral genome GFF3 is usually
        # Ensembl-style "1"; GENCODE/UCSC annotations are "chr1").
        main_chroms = {g.chrom for g in genes.values() if g.chrom}
        main_chroms |= {t.chrom for t in transcripts.values() if t.chrom}
        _harmonize_chroms(rna, main_chroms)
        for tid, t in rna.items():
            transcripts.setdefault(tid, t)   # URS ids won't collide with ENST
        n_rnacentral = len(rna)

    # Build into a temp file then atomically replace, so a crash can't leave a
    # half-written index in place.
    tmp = output.with_suffix(output.suffix + ".tmp")
    tmp.unlink(missing_ok=True)
    Path(str(tmp) + "-journal").unlink(missing_ok=True)

    con = sqlite3.connect(str(tmp))
    cur = con.cursor()
    # Build-time PRAGMAs (fast, unsafe — fine for a throwaway temp file).
    cur.execute("PRAGMA page_size=4096")
    cur.execute("PRAGMA journal_mode=OFF")
    cur.execute("PRAGMA synchronous=OFF")
    cur.execute("PRAGMA temp_store=MEMORY")
    cur.execute("PRAGMA cache_size=-200000")

    has_fts5 = _fts5_available(con)
    has_trigram = has_fts5 and _has_trigram(con)
    if not has_fts5:
        log.warning("FTS5 unavailable: building without FTS tables "
                    "(search falls back to LIKE on indexed columns)")
    elif not has_trigram:
        log.warning("trigram tokenizer unavailable: substring/fuzzy search "
                    "will use a LIKE fallback over search_text")

    cur.executescript(_SCHEMA)
    if has_fts5:
        cur.executescript(_FTS_PREFIX)
        if has_trigram:
            cur.executescript(_FTS_TRIGRAM)

    insert_feature = f"INSERT INTO feature(id, {_FEATURE_COLS}) VALUES (?, {','.join(['?'] * 24)})"
    insert_fts = ("INSERT INTO feature_fts(rowid, display_name, gene_id, "
                  "gene_name, transcript_id, transcript_name, aliases, "
                  "search_text) VALUES (?,?,?,?,?,?,?,?)")
    insert_trgm = "INSERT INTO feature_trigram(rowid, search_text) VALUES (?,?)"
    insert_alias = ("INSERT INTO alias(feature_id, alias_norm, alias_type) "
                    "VALUES (?,?,?)")

    def _alias_rows(fid, aliases):
        """Yield (fid, alias_norm, alias_type) deduped per feature — different
        original strings (OCT-4 / OCT4) can normalize to the same token, and we
        no longer store the original, so collapse exact (norm, type) duplicates."""
        seen = set()
        for alias, atype in aliases.items():
            nrm = normalize(alias)
            key = (nrm, atype)
            if not nrm or key in seen:
                continue
            seen.add(key)
            yield (fid, nrm, atype)

    # Deterministic order: genes first, then transcripts, each sorted by
    # (chrom, start, id). Feature ids are therefore stable across runs.
    gene_items = sorted(genes.values(),
                        key=lambda g: (g.chrom, g.start, g.gene_id))
    tx_items = sorted(transcripts.values(),
                      key=lambda t: (t.chrom, t.start, t.transcript_id))

    n_genes = n_tx = n_alias = 0
    cur.execute("BEGIN")
    fid = 0
    for g in gene_items:
        fid += 1
        row = _gene_row(g)
        # Null out feature.search_text in the table — it is fully redundant with
        # the FTS index copies and never read back at query time (saves ~52 MB).
        cur.execute(insert_feature, (fid, *row[:22], None, row[23]))
        if has_fts5:
            cur.execute(insert_fts, (fid, g.gene_name or g.gene_id, g.gene_id,
                                     g.gene_name, None, None,
                                     " ".join(g.aliases.keys()), row[22]))
            if has_trigram:
                trg = (row[22] if fuzzy_scope == "all"
                       else _fuzzy_text([g.gene_name], g.aliases))
                if trg:
                    cur.execute(insert_trgm, (fid, trg))
        for arow in _alias_rows(fid, g.aliases):
            cur.execute(insert_alias, arow)
            n_alias += 1
        n_genes += 1

    for t in tx_items:
        fid += 1
        row = _tx_row(t)
        cur.execute(insert_feature, (fid, *row[:22], None, row[23]))  # null search_text
        if has_fts5:
            display = t.transcript_name or t.transcript_id
            cur.execute(insert_fts, (fid, display, row[1], t.gene_name,
                                     t.transcript_id, t.transcript_name,
                                     " ".join(t.aliases.keys()), row[22]))
            if has_trigram:
                trg = (row[22] if fuzzy_scope == "all"
                       else _fuzzy_text([t.gene_name, t.transcript_name], t.aliases))
                if trg:
                    cur.execute(insert_trgm, (fid, trg))
        for arow in _alias_rows(fid, t.aliases):
            cur.execute(insert_alias, arow)
            n_alias += 1
        n_tx += 1
    con.commit()

    # metadata
    meta = {
        "schema_version": "1",
        "generator": "refbox.sqlite_index",
        "source_name": source_name,
        "species": species,
        "genome": genome,
        "annotation_version": annotation_version,
        "input_file": input_path.name,
        "input_format": "GFF3" if is_gff3 else "GTF",
        "synonyms_file": Path(synonyms).name if synonyms else "",
        "n_synonyms_injected": str(n_synonyms),
        "rnacentral_file": Path(rnacentral).name if rnacentral else "",
        "n_rnacentral_transcripts": str(n_rnacentral),
        "coord_convention": (
            "start/end and *_start/*_end are 1-based inclusive; "
            "chrom_start0/chrom_end0 are 0-based half-open"),
        "n_genes": str(n_genes),
        "n_transcripts": str(n_tx),
        "n_aliases": str(n_alias),
        "fts5": "1" if has_fts5 else "0",
        "trigram": "1" if has_trigram else "0",
        "fuzzy_scope": fuzzy_scope,
        "sqlite_version": sqlite3.sqlite_version,
        "build_seconds_parse": f"{t_parse:.2f}",
    }
    cur.executemany("INSERT OR REPLACE INTO metadata(key, value) VALUES (?,?)",
                    list(meta.items()))
    con.commit()

    # indexes after bulk load
    for stmt in _INDEXES:
        cur.execute(stmt)
    con.commit()

    # finalize: re-enable a journal so VACUUM is safe, then optimize/compact
    log.info("optimizing (ANALYZE / optimize / VACUUM)…")
    cur.execute("PRAGMA journal_mode=DELETE")
    cur.execute("ANALYZE")
    try:
        cur.execute("PRAGMA optimize")
    except sqlite3.OperationalError:
        pass
    con.commit()
    cur.execute("VACUUM")
    con.commit()
    con.close()

    tmp.replace(output)
    Path(str(tmp) + "-journal").unlink(missing_ok=True)
    log.info("wrote %s (%.1f MB) in %.1fs: %d genes, %d transcripts, %d aliases",
             output, output.stat().st_size / 1e6, time.time() - t0,
             n_genes, n_tx, n_alias)
    return output


# ── search (Python reference implementation for testing) ──────────────────────

# (rank score, matched_field label) per tier — lower score = better.
_RANK = {
    "transcript_id_exact": 1, "transcript_name_exact": 2,
    "gene_name_exact": 3, "gene_id_exact": 4, "alias_exact": 5,
    "prefix": 6, "trigram": 7, "like": 8,
}

# alias_type → (rank score, matched_field) for the unified exact lookup. Any
# type not listed (gene_synonym, hgnc, havana, ccds, protein_id, rnacentral_*,
# dbxref, refseq, name, alias, …) ranks as a generic alias_exact.
_ALIAS_FIELD = {
    "transcript_id":   (1, "transcript_id_exact"),
    "transcript_name": (2, "transcript_name_exact"),
    "gene_name":       (3, "gene_name_exact"),
    "gene_id":         (4, "gene_id_exact"),
}

_SELECT = (
    "id, feature_type, gene_name, gene_id, transcript_name, transcript_id, "
    "chrom, start, end, strand, biotype")


def _row_to_dict(row) -> dict:
    keys = ("id", "feature_type", "gene_name", "gene_id", "transcript_name",
            "transcript_id", "chrom", "start", "end", "strand", "biotype")
    return dict(zip(keys, row))


def _rerank(rows: list, norm_q: str) -> list:
    """Order fuzzy (prefix / trigram) candidates by relevance.

    SQLite returns ``MATCH`` candidates in rowid order, which buries the obvious
    hit (``p53`` would surface ``RNA5SP53`` before ``TP53``). Re-sort so that a
    name which *starts with* the query wins, then one that merely *contains* it,
    then the shortest name (tightest match) — e.g. ``p53`` → ``TP53`` first.
    """
    def key(r):
        name = r[4] or r[2] or r[5] or r[3] or ""   # tx_name|gene_name|tx_id|gene_id
        nn = normalize(name)
        starts = 0 if norm_q and nn.startswith(norm_q) else 1
        contains = 0 if norm_q and norm_q in nn else 1
        return (starts, contains, len(name), name)
    return sorted(rows, key=key)


def _fts_query_token(query: str) -> str:
    """Sanitize a user query into a single safe FTS5 token (quoted)."""
    cleaned = re.sub(r'["\']', " ", query).strip()
    cleaned = re.sub(r"\s+", " ", cleaned)
    return cleaned


def search(
    con: sqlite3.Connection, query: str, *, limit: int = 10,
) -> list[dict]:
    """Ranked search over the index. Returns up to ``limit`` result dicts with a
    ``rank_score`` and ``matched_field``. Implements the priority:

    exact transcript_id → transcript_name → gene_name → gene_id → alias exact →
    prefix (FTS) → trigram/substring → LIKE fallback.
    """
    query = query.strip()
    if not query:
        return []
    cur = con.cursor()
    norm_q = normalize(query)
    seen: set[int] = set()
    results: list[dict] = []

    def take(rows, matched_field):
        score = _RANK[matched_field]
        for row in rows:
            d = _row_to_dict(row)
            if d["id"] in seen:
                continue
            seen.add(d["id"])
            d["rank_score"] = score
            d["matched_field"] = matched_field
            results.append(d)
            if len(results) >= limit:
                return True
        return False

    # Tiers 1–5 collapsed into ONE index-only seek into
    # idx_alias_norm(alias_norm, alias_type, feature_id). A single lookup
    # resolves any exact match — transcript_id / transcript_name / gene_name /
    # gene_id / synonym / RefSeq / RNAcentral ID / … — across ALL species
    # (ENSMUSG, ENSDARG, FBgn, AT1G…), with alias_type giving both the rank and
    # the matched_field. Version-insensitive for free (alias_norm is already
    # version-stripped). ~8 page reads instead of ~22.
    cand = cur.execute(
        "SELECT alias_type, feature_id FROM alias WHERE alias_norm = ? LIMIT ?",
        (norm_q, max(limit * 8, 50))).fetchall()
    if cand:
        best: dict[int, tuple[int, str]] = {}   # feature_id -> (score, field)
        for atype, fid_ in cand:
            score, field = _ALIAS_FIELD.get(atype, (5, "alias_exact"))
            cur_best = best.get(fid_)
            if cur_best is None or score < cur_best[0]:
                best[fid_] = (score, field)
        # rank by (score, feature_id) — genes (lower ids) precede transcripts on ties
        ordered = sorted(best.items(), key=lambda kv: (kv[1][0], kv[0]))[:limit]
        ids = [fid_ for fid_, _ in ordered]
        ph = ",".join("?" * len(ids))
        rowmap = {r[0]: r for r in cur.execute(
            f"SELECT {_SELECT} FROM feature WHERE id IN ({ph})", ids).fetchall()}
        for fid_, (score, field) in ordered:
            r = rowmap.get(fid_)
            if r is None:
                continue
            d = _row_to_dict(r)
            d["rank_score"] = score
            d["matched_field"] = field
            seen.add(d["id"])
            results.append(d)
        # An exact hit must not be diluted (or slowed) by the fuzzy tiers below.
        return results

    # Tier 6: prefix / autocomplete via feature_fts (if present).
    has_fts = _table_exists(con, "feature_fts")
    if has_fts:
        token = _fts_query_token(query)
        if token:
            try:
                rows = cur.execute(
                    f"SELECT {_SELECT} FROM feature WHERE id IN "
                    f"(SELECT rowid FROM feature_fts WHERE feature_fts MATCH ?) "
                    f"LIMIT ?", (f'"{token}"*', limit * 4)).fetchall()
                if take(_rerank(rows, norm_q), "prefix"):
                    return results
            except sqlite3.OperationalError:
                pass

    # Tier 7: trigram substring search (if present). Skip it for queries with no
    # letter — the trigram holds only names/synonyms (never IDs), so a pure
    # digit/symbol fragment like "000003351" can never match a name; running it
    # would just scan postings for nothing (the old pathological case).
    if (_table_exists(con, "feature_trigram") and len(query) >= 3
            and any(c.isalpha() for c in query)):
        try:
            rows = cur.execute(
                f"SELECT {_SELECT} FROM feature WHERE id IN "
                f"(SELECT rowid FROM feature_trigram WHERE feature_trigram MATCH ?) "
                f"LIMIT ?", (f'"{query}"', limit * 4)).fetchall()
            if take(_rerank(rows, norm_q), "trigram"):
                return results
        except sqlite3.OperationalError:
            pass

    # Tier 8: LIKE fallback on the normalized alias table. This is a contains
    # scan, so only run it when no trigram index served the substring case
    # (otherwise it is redundant and slow). A prefix LIKE still uses idx_alias_norm.
    if not _table_exists(con, "feature_trigram"):
        rows = cur.execute(
            f"SELECT {_SELECT} FROM feature WHERE id IN "
            f"(SELECT feature_id FROM alias WHERE alias_norm LIKE ?) LIMIT ?",
            (f"{norm_q}%", limit)).fetchall()  # prefix → index-usable
        if not take(rows, "like"):
            rows = cur.execute(
                f"SELECT {_SELECT} FROM feature WHERE id IN "
                f"(SELECT feature_id FROM alias WHERE alias_norm LIKE ?) LIMIT ?",
                (f"%{norm_q}%", limit)).fetchall()
            take(rows, "like")
    return results


def _table_exists(con: sqlite3.Connection, name: str) -> bool:
    return con.execute(
        "SELECT 1 FROM sqlite_master WHERE type IN ('table','view') AND name=?",
        (name,)).fetchone() is not None


# ── inspection ────────────────────────────────────────────────────────────────

def inspect(db_path: Path) -> dict:
    """Return a summary dict describing the index (used by the inspect CLI)."""
    db_path = Path(db_path)
    con = sqlite3.connect(f"file:{db_path}?mode=ro", uri=True)
    cur = con.cursor()

    def scalar(sql, *args):
        r = cur.execute(sql, args).fetchone()
        return r[0] if r else None

    info: dict = {
        "file": str(db_path),
        "file_size_bytes": db_path.stat().st_size,
        "sqlite_version": sqlite3.sqlite_version,
        "fts5_available": _fts5_available_ro(con),
    }
    info["n_features"] = scalar("SELECT COUNT(*) FROM feature")
    info["n_genes"] = scalar("SELECT COUNT(*) FROM feature WHERE feature_type='gene'")
    info["n_transcripts"] = scalar(
        "SELECT COUNT(*) FROM feature WHERE feature_type<>'gene'")
    info["n_aliases"] = scalar("SELECT COUNT(*) FROM alias")
    info["has_fts"] = _table_exists(con, "feature_fts")
    info["has_trigram"] = _table_exists(con, "feature_trigram")
    # feature_fts is contentless + detail=none (does not support COUNT scans);
    # it is populated 1:1 with feature rows, so report the feature count.
    info["n_fts"] = info["n_features"] if info["has_fts"] else 0
    info["metadata"] = dict(cur.execute("SELECT key, value FROM metadata").fetchall())
    info["examples"] = [
        _row_to_dict(r) for r in cur.execute(
            f"SELECT {_SELECT} FROM feature WHERE transcript_id IS NOT NULL "
            f"ORDER BY id LIMIT 5").fetchall()
    ]
    con.close()
    return info


def _fts5_available_ro(con: sqlite3.Connection) -> bool:
    try:
        con.execute("SELECT fts5(?)", ("",))
    except sqlite3.OperationalError as e:
        return "no such function" not in str(e)
    except Exception:
        return True
    return True


def open_readonly(db_path: Path) -> sqlite3.Connection:
    """Open the index read-only (the way a static host would serve it)."""
    return sqlite3.connect(f"file:{Path(db_path)}?mode=ro", uri=True)
