"""Lightweight RBrowser Annotation Index (``.rbai``) — a small, fast, static
gene/transcript lookup index built from a GTF/GFF3 annotation.

Unlike the full RBrowser Index (``.rba``, SQLite + FTS5 + trigram, ~1 GB), this
is a *compact lookup* index (a few tens of MB): no FTS5, no trigram, no exon/CDS
structure, no payload JSON. Search resolves to plain B-tree seeks on a normalized
``term`` table, so exact / prefix / autocomplete queries are sub-millisecond.

Internally it is an ordinary SQLite database; the ``.rbai`` extension + the
``metadata`` rows (format_name = "RBrowser Annotation Index", RBAI) let RBrowser
verify the format.

Schema
------
``record``  one row per gene or transcript, carrying only display + position.
``term``    (term_norm, term_raw, field, record_id, priority) WITHOUT ROWID —
            the search index; PK leads with term_norm so exact (``= ?``) and
            prefix (``>= lo AND < hi``) are both B-tree range scans.
``gram3``   optional 3-gram → record map for fuzzy/partial recall.
``metadata``  key/value format + provenance.

Stdlib only (sqlite3, gzip, re). Streams the input; memory ∝ #genes+#transcripts.
"""

from __future__ import annotations

import gzip
import logging
import re
import sqlite3
import time
from pathlib import Path

log = logging.getLogger("refbox.lite_index")

SCHEMA_VERSION = 1
FORMAT_NAME = "RBrowser Annotation Index"
FORMAT_SHORT_NAME = "RBAI"

# field → priority (lower = ranked higher). Drives ORDER BY in search.
FIELD_PRIORITY = {
    "gene_name": 10,
    "gene_id": 20,
    "transcript_name": 30,
    "transcript_id": 40,
    "transcript_biotype": 80,
}
# fields that get 3-grams: names only. IDs and biotype are exact/prefix only —
# fuzzy-matching an accession or a biotype is pointless and bloats the table
# (e.g. "proteincoding" grams × every coding transcript).
GRAM3_FIELDS = ("gene_name", "transcript_name")

# GFF3 feature types treated as genes / transcripts.
_GFF_GENE_TYPES = {"gene", "ncrna_gene", "pseudogene"}
_GFF_TX_TYPES = {
    "mrna", "transcript", "lnc_rna", "lncrna", "ncrna", "rrna", "mt_rrna",
    "trna", "mt_trna", "snrna", "snorna", "scarna", "mirna", "misc_rna",
    "ribozyme", "srna", "scrna", "vault_rna", "y_rna", "antisense_rna",
    "pseudogenic_transcript", "processed_transcript", "primary_transcript",
    "rnase_mrp_rna", "rnase_p_rna", "telomerase_rna", "srp_rna", "tmrna",
    "guide_rna", "c_gene_segment", "d_gene_segment", "j_gene_segment",
    "v_gene_segment", "three_prime_overlapping_ncrna",
}
# child/part features used only to infer a parent span (never become records).
_PART_TYPES = {
    "exon", "cds", "noncoding_exon", "three_prime_utr", "five_prime_utr",
    "utr", "start_codon", "stop_codon", "intron", "selenocysteine",
    "stop_codon_redefined_as_selenocysteine",
}

# ── normalization ─────────────────────────────────────────────────────────────

_VERSION = re.compile(r"\.\d+$")          # trailing ".<int>" Ensembl version
_SEP = re.compile(r"[\s_.\-]+")           # separators removed for the sep-free form


def normalize(value: str) -> str:
    """Main normalized form: trim + lowercase (separators preserved).
    ``TP53`` → ``tp53``; ``TP53-201`` → ``tp53-201``; ``protein_coding`` →
    ``protein_coding``."""
    return value.strip().lower()


def strip_version(norm: str) -> str:
    """Drop a trailing Ensembl ``.<version>``: ``enst00000269305.10`` →
    ``enst00000269305``. Operates on an already-normalized string."""
    return _VERSION.sub("", norm)


def sepfree(norm: str) -> str:
    """Separator-free form: drop ``_ - . space``. ``tp53-201`` → ``tp53201``;
    ``protein_coding`` → ``proteincoding``."""
    return _SEP.sub("", norm)


def term_forms(raw: str) -> list[str]:
    """Distinct normalized lookup keys for one raw value: the normalized form,
    its version-stripped form, and its separator-free form (deduped, order
    preserved). ``term_raw`` keeps the original for display."""
    if not raw:
        return []
    n = normalize(raw)
    forms = [n]
    for cand in (strip_version(n), sepfree(n)):
        if cand and cand not in forms:
            forms.append(cand)
    return forms


def prefix_bounds(prefix_norm: str) -> "tuple[str, str]":
    """``[lo, hi)`` range for a B-tree prefix scan. ``tp`` → (``tp``, ``tq``);
    increments the last code point so ``term_norm >= lo AND term_norm < hi``
    selects exactly the rows starting with ``prefix_norm``."""
    lo = prefix_norm
    hi = lo[:-1] + chr(ord(lo[-1]) + 1) if lo else lo
    return lo, hi


# ── attribute parsing ─────────────────────────────────────────────────────────

_GTF_ATTR = re.compile(r'(\w+)\s+"([^"]*)"')
_GFF_ATTR = re.compile(r"([^=;]+)=([^;]*)")


def _parse_gtf_attrs(field: str) -> "dict[str, str]":
    return {k: v for k, v in _GTF_ATTR.findall(field)}


def _parse_gff_attrs(field: str) -> "dict[str, str]":
    out: dict[str, str] = {}
    for chunk in field.rstrip(";").split(";"):
        chunk = chunk.strip()
        if not chunk or "=" not in chunk:
            continue
        k, v = chunk.split("=", 1)
        out[k.strip()] = v.strip()
    return out


def _open_text(path: Path):
    """Open plain or gzip text transparently (sniff magic bytes, not just .gz)."""
    with open(path, "rb") as fh:
        magic = fh.read(2)
    if magic == b"\x1f\x8b":
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return open(path, "rt", encoding="utf-8", errors="replace")


def detect_dialect(path: Path) -> str:
    """Return 'gtf' or 'gff3' by sniffing the first data line's attribute style
    (``key "value";`` = GTF, ``key=value`` = GFF3); fall back to extension."""
    with _open_text(path) as fh:
        for line in fh:
            if not line or line[0] == "#":
                if line.startswith("##gff-version"):
                    return "gff3"
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9:
                continue
            attrs = cols[8]
            if "=" in attrs and '"' not in attrs.split("=", 1)[0]:
                return "gff3"
            if re.search(r'\w+\s+"', attrs):
                return "gtf"
    name = path.name.lower()
    return "gtf" if ".gtf" in name else "gff3"


def _strip_id_prefix(value: str) -> str:
    """Ensembl GFF3 ids look like ``gene:ENSG…`` / ``transcript:ENST…``; strip the
    ``<type>:`` prefix so they match the clean ``gene_id``/``transcript_id``."""
    if value and ":" in value:
        head, rest = value.split(":", 1)
        if head in ("gene", "transcript"):
            return rest
    return value


# ── parsed record containers (plain dicts for low overhead) ───────────────────

def _new_span(chrom: str, strand: str, start: int, end: int) -> dict:
    return {"chrom": chrom, "strand": strand, "start": start, "end": end}


def _grow(span: dict, start: int, end: int) -> None:
    if start < span["start"]:
        span["start"] = start
    if end > span["end"]:
        span["end"] = end


def parse_annotation(path: Path, dialect: str) -> "tuple[dict, dict]":
    """Stream the annotation into ``genes`` and ``txs`` dicts keyed by stable id.

    Returns ``(genes, txs)`` where each value carries name/id/biotype + a span
    (chrom, strand, start, end). Spans accumulate over all of a feature's lines
    (gene/transcript line + exon/CDS children), so a span is inferred even when
    the parent feature line itself is absent.
    """
    genes: dict[str, dict] = {}
    txs: dict[str, dict] = {}
    n_lines = 0

    if dialect == "gtf":
        with _open_text(path) as fh:
            for raw in fh:
                if not raw or raw[0] == "#":
                    continue
                cols = raw.rstrip("\n").split("\t")
                if len(cols) < 9:
                    continue
                n_lines += 1
                ftype = cols[2]
                chrom, strand = cols[0], cols[6]
                try:
                    start, end = int(cols[3]), int(cols[4])
                except ValueError:
                    continue
                a = _parse_gtf_attrs(cols[8])
                gid = a.get("gene_id")
                tid = a.get("transcript_id")
                gname = a.get("gene_name")
                gbio = a.get("gene_type") or a.get("gene_biotype")
                tbio = (a.get("transcript_type") or a.get("transcript_biotype")
                        or gbio)

                if gid:
                    g = genes.get(gid)
                    if g is None:
                        g = genes[gid] = {"gene_id": gid, "gene_name": gname,
                                          "biotype": gbio,
                                          **_new_span(chrom, strand, start, end)}
                    else:
                        _grow(g, start, end)
                    if ftype == "gene":
                        g["chrom"], g["strand"] = chrom, strand
                        if gname:
                            g["gene_name"] = gname
                        if gbio:
                            g["biotype"] = gbio
                    else:
                        if gname and not g["gene_name"]:
                            g["gene_name"] = gname
                        if gbio and not g["biotype"]:
                            g["biotype"] = gbio
                if tid:
                    t = txs.get(tid)
                    if t is None:
                        t = txs[tid] = {"transcript_id": tid,
                                        "transcript_name": a.get("transcript_name"),
                                        "gene_id": gid, "gene_name": gname,
                                        "biotype": tbio,
                                        **_new_span(chrom, strand, start, end)}
                    else:
                        _grow(t, start, end)
                    if ftype == "transcript":
                        t["chrom"], t["strand"] = chrom, strand
                        if a.get("transcript_name"):
                            t["transcript_name"] = a["transcript_name"]
                        if tbio:
                            t["biotype"] = tbio
                    if gid and not t["gene_id"]:
                        t["gene_id"] = gid
                    if gname and not t["gene_name"]:
                        t["gene_name"] = gname
                    if tbio and not t["biotype"]:
                        t["biotype"] = tbio
        _resolve_gtf_genes(genes, txs)
        log.info("lite: parsed %d GTF lines → %d genes, %d transcripts",
                 n_lines, len(genes), len(txs))
        return genes, txs

    # ── GFF3 ──
    gene_by_fid: dict[str, dict] = {}   # feature ID  -> gene dict
    tx_by_fid: dict[str, dict] = {}     # feature ID  -> tx dict (+ parent)
    with _open_text(path) as fh:
        for raw in fh:
            if not raw or raw[0] == "#":
                continue
            cols = raw.rstrip("\n").split("\t")
            if len(cols) < 9:
                continue
            n_lines += 1
            ftype = cols[2]
            low = ftype.lower()
            chrom, strand = cols[0], cols[6]
            try:
                start, end = int(cols[3]), int(cols[4])
            except ValueError:
                continue
            a = _parse_gff_attrs(cols[8])
            fid = a.get("ID")
            parent = a.get("Parent")

            if low in _GFF_GENE_TYPES:
                gid = a.get("gene_id") or _strip_id_prefix(fid or "")
                if not gid:
                    continue
                gname = a.get("gene_name") or a.get("Name") or a.get("gene")
                gbio = (a.get("gene_type") or a.get("gene_biotype")
                        or a.get("biotype"))
                rec = {"gene_id": gid, "gene_name": gname, "biotype": gbio,
                       **_new_span(chrom, strand, start, end)}
                if fid:
                    gene_by_fid[fid] = rec
                genes.setdefault(gid, rec)
            elif low in _GFF_TX_TYPES or (
                    a.get("transcript_id") and parent and low not in _PART_TYPES):
                tid = a.get("transcript_id") or _strip_id_prefix(fid or "")
                if not tid:
                    continue
                tbio = (a.get("transcript_type") or a.get("transcript_biotype")
                        or a.get("biotype"))
                rec = {"transcript_id": tid,
                       "transcript_name": a.get("transcript_name") or a.get("Name"),
                       "gene_id": a.get("gene_id"),
                       "gene_name": a.get("gene_name"),
                       "biotype": tbio, "parent": parent,
                       **_new_span(chrom, strand, start, end)}
                if fid:
                    tx_by_fid[fid] = rec
                txs[tid] = rec
            elif low in _PART_TYPES and parent:
                # extend the parent transcript span (exon/CDS/UTR)
                for p in parent.split(","):
                    t = tx_by_fid.get(p)
                    if t:
                        _grow(t, start, end)

    _resolve_gff_genes(genes, gene_by_fid, tx_by_fid)
    log.info("lite: parsed %d GFF3 lines → %d genes, %d transcripts",
             n_lines, len(genes), len(txs))
    return genes, txs


def _resolve_gtf_genes(genes: dict, txs: dict) -> None:
    """Fill any missing transcript→gene name from the gene record; synthesize a
    gene span from transcripts for gene_ids that never had a gene feature."""
    for t in txs.values():
        gid = t.get("gene_id")
        if gid:
            g = genes.get(gid)
            if g is None:
                g = genes[gid] = {"gene_id": gid, "gene_name": t.get("gene_name"),
                                  "biotype": t.get("biotype"),
                                  **_new_span(t["chrom"], t["strand"],
                                              t["start"], t["end"])}
            else:
                _grow(g, t["start"], t["end"])
            if not t.get("gene_name") and g.get("gene_name"):
                t["gene_name"] = g["gene_name"]


def _resolve_gff_genes(genes: dict, gene_by_fid: dict, tx_by_fid: dict) -> None:
    """Resolve each transcript's gene via its Parent feature ID; synthesize gene
    records (span from transcripts) for parents that aren't gene features."""
    for t in tx_by_fid.values():
        g = None
        parent = t.get("parent")
        if parent:
            for p in parent.split(","):
                g = gene_by_fid.get(p)
                if g:
                    break
        if g:
            if not t.get("gene_id"):
                t["gene_id"] = g["gene_id"]
            if not t.get("gene_name"):
                t["gene_name"] = g.get("gene_name")
            _grow(g, t["start"], t["end"])
        gid = t.get("gene_id")
        if gid and gid not in genes:
            genes[gid] = {"gene_id": gid, "gene_name": t.get("gene_name"),
                          "biotype": t.get("biotype"),
                          **_new_span(t["chrom"], t["strand"],
                                      t["start"], t["end"])}
        elif gid:
            _grow(genes[gid], t["start"], t["end"])


# ── position string ────────────────────────────────────────────────────────────

def position_str(chrom: str, start: int, end: int, strand: str) -> str:
    """``chrom:start-end`` (+ ``(strand)`` when strand is + or -)."""
    base = f"{chrom}:{start}-{end}"
    return f"{base}({strand})" if strand in ("+", "-") else base


# ── schema ─────────────────────────────────────────────────────────────────────

_SCHEMA = """
CREATE TABLE record (
    id INTEGER PRIMARY KEY,
    feature_type TEXT NOT NULL,
    gene_name TEXT,
    gene_id TEXT,
    transcript_name TEXT,
    transcript_id TEXT,
    transcript_biotype TEXT,
    chrom TEXT NOT NULL,
    start INTEGER NOT NULL,
    end INTEGER NOT NULL,
    strand TEXT,
    genome_position_str TEXT NOT NULL
);
CREATE TABLE term (
    term_norm TEXT NOT NULL,
    term_raw TEXT NOT NULL,
    field TEXT NOT NULL,
    record_id INTEGER NOT NULL,
    priority INTEGER NOT NULL,
    PRIMARY KEY (term_norm, priority, record_id, field)
) WITHOUT ROWID;
CREATE TABLE gram3 (
    gram TEXT NOT NULL,
    record_id INTEGER NOT NULL,
    term_norm TEXT NOT NULL,
    field TEXT NOT NULL,
    PRIMARY KEY (gram, record_id, term_norm, field)
) WITHOUT ROWID;
CREATE TABLE metadata (
    key TEXT PRIMARY KEY,
    value TEXT
);
"""

_INDEXES = """
CREATE INDEX idx_term_record_id ON term(record_id);
CREATE INDEX idx_record_gene_name ON record(gene_name);
CREATE INDEX idx_record_gene_id ON record(gene_id);
CREATE INDEX idx_record_transcript_name ON record(transcript_name);
CREATE INDEX idx_record_transcript_id ON record(transcript_id);
CREATE INDEX idx_record_biotype ON record(transcript_biotype);
"""


def _grams(norm: str) -> "set[str]":
    """3-grams of the separator-free form (so grams ignore _ - . spaces)."""
    s = sepfree(norm)
    return {s[i:i + 3] for i in range(len(s) - 2)} if len(s) >= 3 else set()


# ── build ──────────────────────────────────────────────────────────────────────

def build_lite_index(input_path, output=None, *, source_name="", species="",
                     genome="", annotation_version="", enable_gram3=False,
                     force=False, verbose=False) -> Path:
    """Build a ``.rbai`` lite index from a GTF/GFF/GFF3 (.gz ok) annotation."""
    input_path = Path(input_path)
    if output is None:
        stem = input_path.name
        for ext in (".gz", ".gtf", ".gff3", ".gff"):
            if stem.lower().endswith(ext):
                stem = stem[: -len(ext)]
        output = input_path.with_name(stem + ".rbai")
    output = Path(output)
    if output.exists() and not force:
        raise FileExistsError(f"{output} exists; pass force=True / --force")
    if output.exists():
        output.unlink()

    t0 = time.perf_counter()
    dialect = detect_dialect(input_path)
    log.info("lite: annotation dialect: %s", dialect.upper())
    genes, txs = parse_annotation(input_path, dialect)

    con = sqlite3.connect(output)
    cur = con.cursor()
    # fast bulk-load pragmas
    cur.execute("PRAGMA journal_mode=OFF")
    cur.execute("PRAGMA synchronous=OFF")
    cur.execute("PRAGMA temp_store=MEMORY")
    cur.execute("PRAGMA cache_size=-200000")
    cur.execute("PRAGMA page_size=4096")
    cur.executescript(_SCHEMA)

    # Deterministic record ids: genes (sorted by gene_id) then transcripts
    # (sorted by transcript_id). Build records + terms + grams in batches.
    rec_rows: list[tuple] = []
    term_rows: list[tuple] = []
    gram_rows: list[tuple] = []
    rid = 0

    def add_terms(record_id: int, field: str, raw: str) -> None:
        if not raw:
            return
        pr = FIELD_PRIORITY[field]
        for form in term_forms(raw):
            term_rows.append((form, raw, field, record_id, pr))
        if enable_gram3 and field in GRAM3_FIELDS:
            n = normalize(raw)
            for g in _grams(n):
                gram_rows.append((g, record_id, n, field))

    def emit_gene(gid: int, g: dict) -> None:
        nonlocal rid
        rid += 1
        pos = position_str(g["chrom"], g["start"], g["end"], g["strand"])
        rec_rows.append((rid, "gene", g.get("gene_name"), gid, None, None, None,
                         g["chrom"], g["start"], g["end"], g["strand"], pos))
        add_terms(rid, "gene_name", g.get("gene_name") or "")
        add_terms(rid, "gene_id", gid)

    def emit_tx(tid: int, t: dict) -> None:
        nonlocal rid
        rid += 1
        pos = position_str(t["chrom"], t["start"], t["end"], t["strand"])
        rec_rows.append((rid, "transcript", t.get("gene_name"), t.get("gene_id"),
                         t.get("transcript_name"), tid, t.get("biotype"),
                         t["chrom"], t["start"], t["end"], t["strand"], pos))
        add_terms(rid, "gene_name", t.get("gene_name") or "")
        add_terms(rid, "gene_id", t.get("gene_id") or "")
        add_terms(rid, "transcript_name", t.get("transcript_name") or "")
        add_terms(rid, "transcript_id", tid)
        add_terms(rid, "transcript_biotype", t.get("biotype") or "")

    # Cluster each gene immediately followed by ITS transcripts (both sorted),
    # so a gene-name match lands on contiguous record ids → contiguous pages →
    # far fewer HTTP Range fetches than genes-block-then-transcripts-block.
    tx_by_gene: dict[str, list] = {}
    orphans: list[str] = []
    for tid in sorted(txs):
        gid = txs[tid].get("gene_id")
        if gid in genes:
            tx_by_gene.setdefault(gid, []).append(tid)
        else:
            orphans.append(tid)
    for gid in sorted(genes):
        emit_gene(gid, genes[gid])
        for tid in tx_by_gene.get(gid, ()):
            emit_tx(tid, txs[tid])
    for tid in orphans:                       # transcripts with no resolvable gene
        emit_tx(tid, txs[tid])

    with con:
        cur.executemany(
            "INSERT INTO record VALUES (?,?,?,?,?,?,?,?,?,?,?,?)", rec_rows)
        cur.executemany(
            "INSERT OR IGNORE INTO term VALUES (?,?,?,?,?)", term_rows)
        if enable_gram3 and gram_rows:
            cur.executemany(
                "INSERT OR IGNORE INTO gram3 VALUES (?,?,?,?)", gram_rows)

    n_terms = cur.execute("SELECT count(*) FROM term").fetchone()[0]
    n_grams = cur.execute("SELECT count(*) FROM gram3").fetchone()[0]
    n_genes = sum(1 for r in rec_rows if r[1] == "gene")
    n_tx = len(rec_rows) - n_genes

    cur.executescript(_INDEXES)

    meta = {
        "format_name": FORMAT_NAME,
        "format_short_name": FORMAT_SHORT_NAME,
        "schema_version": str(SCHEMA_VERSION),
        "storage_engine": "SQLite",
        "index_type": "lite_btree_lookup",
        "description": "Lightweight gene/transcript search index for RBrowser",
        "source_name": source_name,
        "species": species,
        "genome": genome,
        "annotation_version": annotation_version,
        "input_file": input_path.name,
        "input_format": dialect.upper(),
        "gram3_enabled": "1" if enable_gram3 else "0",
        "n_records": str(len(rec_rows)),
        "n_genes": str(n_genes),
        "n_transcripts": str(n_tx),
        "n_terms": str(n_terms),
        "n_gram3": str(n_grams),
        "sqlite_version": sqlite3.sqlite_version,
        "generator": "refbox.lite_index",
    }
    with con:
        cur.executemany("INSERT OR REPLACE INTO metadata VALUES (?,?)",
                        sorted(meta.items()))

    log.info("lite: optimizing (ANALYZE / optimize / VACUUM)…")
    cur.execute("ANALYZE")
    cur.execute("PRAGMA optimize")
    con.commit()
    cur.execute("VACUUM")
    con.close()

    size_mb = output.stat().st_size / 1e6
    log.info("lite: wrote %s (%.1f MB) in %.1fs: %d genes, %d transcripts, "
             "%d terms, %d grams", output, size_mb, time.perf_counter() - t0,
             n_genes, n_tx, n_terms, n_grams)
    return output


# ── read side ──────────────────────────────────────────────────────────────────

_RECORD_COLS = ("id", "feature_type", "gene_name", "gene_id", "transcript_name",
                "transcript_id", "transcript_biotype", "chrom", "start", "end",
                "strand", "genome_position_str")


def open_readonly(path) -> sqlite3.Connection:
    con = sqlite3.connect(f"file:{Path(path)}?mode=ro&immutable=1", uri=True)
    con.row_factory = sqlite3.Row
    return con


# NOTE on performance: a low-cardinality term like "protein_coding" matches
# ~130k records. Sorting those by joined record columns would be O(N log N) per
# query. Instead we LIMIT inside the term index *first* — for a fixed term_norm
# the PK (term_norm, priority, record_id, …) is already ordered by priority then
# record_id, so the inner SELECT is an index-only range scan with no sort — then
# join only the top-k records and re-rank just those k rows in the outer query.
_EXACT_SQL = (
    "SELECT r.*, t.field AS matched_field, t.term_raw, t.priority "
    "FROM (SELECT record_id, field, term_raw, priority FROM term "
    "      WHERE term_norm = ? ORDER BY priority, record_id LIMIT ?) t "
    "JOIN record r ON r.id = t.record_id "
    "ORDER BY t.priority, r.feature_type, r.gene_name, r.transcript_name"
)
_PREFIX_SQL = (
    "SELECT r.*, t.field AS matched_field, t.term_raw, t.priority, t.term_norm "
    "FROM (SELECT record_id, field, term_raw, priority, term_norm FROM term "
    "      WHERE term_norm >= ? AND term_norm < ? "
    "      ORDER BY term_norm, priority LIMIT ?) t "
    "JOIN record r ON r.id = t.record_id "
    "ORDER BY t.priority, t.term_norm, r.feature_type"
)


def _dedup(rows, limit):
    """Keep the best (first, since pre-ordered) row per record id."""
    seen = set()
    out = []
    for r in rows:
        rid = r["id"]
        if rid in seen:
            continue
        seen.add(rid)
        out.append(dict(r))
        if len(out) >= limit:
            break
    return out


def search(con: sqlite3.Connection, query: str, limit: int = 10):
    """Ranked lookup. Returns ``(mode, results)`` where mode is one of
    ``exact`` / ``prefix`` / ``normalized`` / ``gram3`` / ``none``.

    Tier order (stops at the first tier that yields rows):
      1. exact term_norm
      2. prefix B-tree range
      3. separator-free exact
      4. version-stripped exact
      5. gram3 fuzzy candidate recall (if present), ranked in Python
    No full-table ``LIKE '%q%'`` ever runs.
    """
    q = (query or "").strip()
    if not q:
        return "none", []
    nq = normalize(q)

    # Inner LIMIT is kept close to `limit` (just a small buffer for the rare
    # case of one record matching a term via several normalized forms). Over-
    # fetching here directly costs scattered record-page reads over HTTP Range,
    # so limit*2 (not *4) roughly halves page fetches for gene-name matches that
    # fan out to many transcripts.
    rows = con.execute(_EXACT_SQL, (nq, limit * 2)).fetchall()
    if rows:
        return "exact", _dedup(rows, limit)

    lo, hi = prefix_bounds(nq)
    if hi:
        rows = con.execute(_PREFIX_SQL, (lo, hi, limit * 3)).fetchall()
        if rows:
            return "prefix", _dedup(rows, limit)

    sf = sepfree(nq)
    if sf and sf != nq:
        rows = con.execute(_EXACT_SQL, (sf, limit * 2)).fetchall()
        if rows:
            return "normalized", _dedup(rows, limit)

    sv = strip_version(nq)
    if sv != nq:
        rows = con.execute(_EXACT_SQL, (sv, limit * 2)).fetchall()
        if rows:
            return "normalized", _dedup(rows, limit)

    res = _gram3_search(con, nq, limit)
    if res:
        return "gram3", res
    return "none", []


def _gram3_search(con, nq, limit):
    """Fuzzy recall: pull candidate records sharing 3-grams with the query, then
    rank in Python by exact/startswith/contains/gram-overlap/priority. Only runs
    when the query is ≥3 chars and the gram3 table is populated."""
    grams = sorted(_grams(nq))
    if not grams:
        return []
    ph = ",".join("?" * len(grams))
    cand = con.execute(
        f"SELECT record_id, COUNT(*) AS hits FROM gram3 WHERE gram IN ({ph}) "
        "GROUP BY record_id ORDER BY hits DESC LIMIT 200", grams).fetchall()
    if not cand:
        return []
    hit_by_id = {r["record_id"]: r["hits"] for r in cand}
    ids = list(hit_by_id)
    rec_ph = ",".join("?" * len(ids))
    recs = con.execute(
        f"SELECT * FROM record WHERE id IN ({rec_ph})", ids).fetchall()
    sf = sepfree(nq)

    def rank(r):
        name = normalize(r["transcript_name"] or r["gene_name"]
                         or r["transcript_id"] or r["gene_id"] or "")
        nsf = sepfree(name)
        exact = 0 if name == nq or nsf == sf else 1
        starts = 0 if nsf.startswith(sf) else 1
        contains = 0 if sf in nsf else 1
        pr = 10 if r["feature_type"] == "gene" else 30
        return (exact, starts, contains, -hit_by_id[r["id"]], pr, len(name))

    ranked = sorted(recs, key=rank)[:limit]
    out = []
    for r in ranked:
        d = dict(r)
        d["matched_field"] = "gram3"
        out.append(d)
    return out


def inspect(path) -> dict:
    """Summary stats for benchmarking/verification (no FTS scan needed)."""
    con = open_readonly(path)
    meta = {k: v for k, v in con.execute("SELECT key, value FROM metadata")}
    info = {
        "file_size": Path(path).stat().st_size,
        "n_records": con.execute("SELECT count(*) FROM record").fetchone()[0],
        "n_genes": con.execute(
            "SELECT count(*) FROM record WHERE feature_type='gene'").fetchone()[0],
        "n_transcripts": con.execute(
            "SELECT count(*) FROM record WHERE feature_type='transcript'"
        ).fetchone()[0],
        "n_terms": con.execute("SELECT count(*) FROM term").fetchone()[0],
        "n_gram3": con.execute("SELECT count(*) FROM gram3").fetchone()[0],
        "metadata": meta,
    }
    con.close()
    return info
