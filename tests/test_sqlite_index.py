"""Tests for the SQLite search-index builder (``refbox.sqlite_index``).

These run without the external bioinformatics tools (pure ``sqlite3``), so they
exercise parsing, schema, coordinate conventions, alias extraction, and the
ranked search over a tiny synthetic GTF and GFF3.
"""

from __future__ import annotations

import gzip
import sqlite3
from pathlib import Path

import pytest

from refbox import sqlite_index as si


# A two-gene fixture covering: protein-coding gene with CDS+UTR, a lncRNA, a
# minus-strand transcript, and a HAVANA/CCDS/HGNC alias spread.
GTF = """\
##provider: TEST
chr17\tHAVANA\tgene\t100\t500\t.\t-\t.\tgene_id "ENSG00000141510.18"; gene_type "protein_coding"; gene_name "TP53"; havana_gene "OTTHUMG1"; hgnc_id "HGNC:11998";
chr17\tHAVANA\ttranscript\t100\t500\t.\t-\t.\tgene_id "ENSG00000141510.18"; transcript_id "ENST00000269305.9"; gene_name "TP53"; transcript_type "protein_coding"; transcript_name "TP53-201"; havana_transcript "OTTHUMT1"; ccdsid "CCDS11118.1";
chr17\tHAVANA\texon\t400\t500\t.\t-\t.\tgene_id "ENSG00000141510.18"; transcript_id "ENST00000269305.9";
chr17\tHAVANA\texon\t100\t200\t.\t-\t.\tgene_id "ENSG00000141510.18"; transcript_id "ENST00000269305.9";
chr17\tHAVANA\tCDS\t150\t450\t.\t-\t0\tgene_id "ENSG00000141510.18"; transcript_id "ENST00000269305.9";
chr17\tHAVANA\tUTR\t451\t500\t.\t-\t.\tgene_id "ENSG00000141510.18"; transcript_id "ENST00000269305.9";
chr17\tHAVANA\tUTR\t100\t149\t.\t-\t.\tgene_id "ENSG00000141510.18"; transcript_id "ENST00000269305.9";
chr11\tENSEMBL\tgene\t1000\t2000\t.\t+\t.\tgene_id "ENSG00000251562.8"; gene_type "lncRNA"; gene_name "MALAT1";
chr11\tENSEMBL\ttranscript\t1000\t2000\t.\t+\t.\tgene_id "ENSG00000251562.8"; transcript_id "ENST00000534336.1"; gene_name "MALAT1"; transcript_type "lncRNA"; transcript_name "MALAT1-201";
chr11\tENSEMBL\texon\t1000\t2000\t.\t+\t.\tgene_id "ENSG00000251562.8"; transcript_id "ENST00000534336.1";
"""

GFF3 = """\
##gff-version 3
chr17\tHAVANA\tgene\t100\t500\t.\t-\t.\tID=ENSG00000141510.18;gene_id=ENSG00000141510.18;gene_type=protein_coding;gene_name=TP53;Alias=p53,LFS1
chr17\tHAVANA\tmRNA\t100\t500\t.\t-\t.\tID=ENST00000269305.9;Parent=ENSG00000141510.18;transcript_id=ENST00000269305.9;gene_name=TP53;transcript_type=protein_coding;transcript_name=TP53-201;Dbxref=RefSeq:NM_000546.6
chr17\tHAVANA\texon\t400\t500\t.\t-\t.\tID=e1;Parent=ENST00000269305.9
chr17\tHAVANA\texon\t100\t200\t.\t-\t.\tID=e2;Parent=ENST00000269305.9
chr17\tHAVANA\tCDS\t150\t450\t.\t-\t0\tID=c1;Parent=ENST00000269305.9
chr17\tHAVANA\tfive_prime_UTR\t451\t500\t.\t-\t.\tID=u5;Parent=ENST00000269305.9
chr17\tHAVANA\tthree_prime_UTR\t100\t149\t.\t-\t.\tID=u3;Parent=ENST00000269305.9
"""


@pytest.fixture
def gtf_db(tmp_path: Path) -> Path:
    src = tmp_path / "annot.gtf"
    src.write_text(GTF)
    out = tmp_path / "idx.sqlite"
    si.build_sqlite_index(src, out, source_name="GENCODE", species="human",
                          genome="hg38", annotation_version="vTEST", force=True)
    return out


# ── normalization ─────────────────────────────────────────────────────────────

@pytest.mark.parametrize("raw,expected", [
    ("ENST00000335137.4", "enst00000335137"),
    ("ENSG00000141510.18", "ensg00000141510"),
    ("TP53-201", "tp53201"),
    ("p53", "p53"),
    ("NM_000546.6", "nm000546"),       # trailing version stripped
])
def test_normalize(raw: str, expected: str):
    assert si.normalize(raw) == expected


def test_strip_version():
    assert si.strip_version("ENST00000269305.9") == "ENST00000269305"
    assert si.strip_version("TP53") == "TP53"


# ── schema / counts ───────────────────────────────────────────────────────────

def test_counts_and_metadata(gtf_db: Path):
    info = si.inspect(gtf_db)
    assert info["n_genes"] == 2
    assert info["n_transcripts"] == 2
    assert info["n_aliases"] > 0
    assert info["has_fts"] and info["has_trigram"]
    assert info["metadata"]["source_name"] == "GENCODE"
    assert info["metadata"]["input_format"] == "GTF"


def test_coordinate_conventions(gtf_db: Path):
    con = sqlite3.connect(gtf_db)
    row = con.execute(
        "SELECT start, end, chrom_start0, chrom_end0, cds_start, cds_end, "
        "utr5_start, utr5_end, utr3_start, utr3_end, exon_count, exon_starts, "
        "exon_ends FROM feature WHERE transcript_id='ENST00000269305.9'"
    ).fetchone()
    start, end, s0, e0, cds_s, cds_e, u5s, u5e, u3s, u3e, ec, es, ee = row
    assert (start, end) == (100, 500)          # 1-based inclusive
    assert (s0, e0) == (99, 500)               # 0-based half-open
    assert (cds_s, cds_e) == (150, 450)
    # minus strand: 5'UTR is on the high-coordinate side
    assert (u5s, u5e) == (451, 500)
    assert (u3s, u3e) == (100, 149)
    assert ec == 2
    assert es == "100,400" and ee == "200,500"  # genomic-ascending
    con.close()


# ── no full table scans (critical for HTTP Range VFS hosting) ─────────────────

# These mirror the SQL templates inside ``search()``. Over an HTTP Range VFS a
# real table SCAN downloads the whole table, so every served query MUST resolve
# to an index seek. A virtual-table (FTS5) "SCAN ... VIRTUAL TABLE INDEX" is the
# normal MATCH plan (index-driven, not a table scan) and is allowed.
_TIER_QUERIES = [
    ("transcript_id",
     "SELECT id FROM feature WHERE transcript_id = ? COLLATE NOCASE LIMIT 10",
     ("ENST00000269305.9",)),
    ("versionless",
     "SELECT id FROM feature WHERE id IN (SELECT feature_id FROM alias "
     "WHERE alias_norm=? AND alias_type=?) LIMIT 10",
     ("enst00000269305", "transcript_id_versionless")),
    ("transcript_name",
     "SELECT id FROM feature WHERE transcript_name = ? COLLATE NOCASE LIMIT 10",
     ("TP53-201",)),
    ("gene_name",
     "SELECT id FROM feature WHERE gene_name = ? COLLATE NOCASE LIMIT 10",
     ("TP53",)),
    ("gene_id",
     "SELECT id FROM feature WHERE gene_id = ? COLLATE NOCASE LIMIT 10",
     ("ENSG00000141510.18",)),
    ("alias_exact",
     "SELECT id FROM feature WHERE id IN (SELECT feature_id FROM alias "
     "WHERE alias_norm=?) LIMIT 10",
     ("tp53",)),
    ("prefix_fts",
     "SELECT id FROM feature WHERE id IN (SELECT rowid FROM feature_fts "
     "WHERE feature_fts MATCH ?) LIMIT 80",
     ('"tp"*',)),
    ("trigram",
     "SELECT id FROM feature WHERE id IN (SELECT rowid FROM feature_trigram "
     "WHERE feature_trigram MATCH ?) LIMIT 80",
     ('"p53"',)),
]


def _is_full_scan(detail: str) -> bool:
    """True for any full scan that downloads a whole table *or* a whole index
    over an HTTP Range VFS. Only an FTS5 virtual-table MATCH ("SCAN ... VIRTUAL
    TABLE INDEX") is index-driven and safe; a "SCAN ... USING COVERING INDEX"
    still reads the entire index and is a hazard.
    """
    d = detail.upper()
    return d.startswith("SCAN") and "VIRTUAL TABLE" not in d


@pytest.fixture
def big_db(tmp_path: Path) -> Path:
    """A few hundred synthetic features so the planner makes production-like
    choices (on a 4-row table it scans everything regardless of indexes)."""
    lines = ["##synthetic"]
    for i in range(400):
        gid = f"ENSG{i:011d}.1"
        tid = f"ENST{i:011d}.2"
        gn = f"GENE{i}"
        lines.append(f'chr1\tT\tgene\t{i*100+1}\t{i*100+90}\t.\t+\t.\t'
                     f'gene_id "{gid}"; gene_type "protein_coding"; gene_name "{gn}";')
        lines.append(f'chr1\tT\ttranscript\t{i*100+1}\t{i*100+90}\t.\t+\t.\t'
                     f'gene_id "{gid}"; transcript_id "{tid}"; gene_name "{gn}"; '
                     f'transcript_name "{gn}-201";')
        lines.append(f'chr1\tT\texon\t{i*100+1}\t{i*100+90}\t.\t+\t.\t'
                     f'gene_id "{gid}"; transcript_id "{tid}";')
    src = tmp_path / "big.gtf"
    src.write_text("\n".join(lines) + "\n")
    out = tmp_path / "big.sqlite"
    si.build_sqlite_index(src, out, force=True)
    return out


def test_no_full_table_scans(big_db: Path):
    con = sqlite3.connect(big_db)
    for label, sql, params in _TIER_QUERIES:
        plan = con.execute("EXPLAIN QUERY PLAN " + sql, params).fetchall()
        offenders = [r[3] for r in plan if _is_full_scan(r[3])]
        assert not offenders, f"tier {label!r} would full-scan: {offenders}"
    con.close()


def test_anti_patterns_do_scan(big_db: Path):
    """Documents the three query shapes the browser must NEVER use — each one
    degrades to a full scan (i.e. a full-DB/index download over HTTP Range)."""
    con = sqlite3.connect(big_db)
    bad = [
        ("SELECT id FROM feature WHERE lower(gene_name) = ? LIMIT 10", ("gene1",)),
        ("SELECT id FROM feature WHERE gene_name LIKE ? LIMIT 10", ("%GENE1%",)),
        ("SELECT id FROM feature WHERE gene_name = ? LIMIT 10", ("gene1",)),  # no COLLATE
    ]
    for sql, a in bad:
        plan = con.execute("EXPLAIN QUERY PLAN " + sql, a).fetchall()
        assert any(_is_full_scan(r[3]) for r in plan), f"expected a scan for: {sql}"
    con.close()


# ── ranked search ─────────────────────────────────────────────────────────────

@pytest.mark.parametrize("query,expect_field,expect_tid", [
    ("ENST00000269305.9", "transcript_id_exact", "ENST00000269305.9"),
    ("ENST00000269305", "transcript_id_exact", "ENST00000269305.9"),  # versionless
    ("TP53-201", "transcript_name_exact", "ENST00000269305.9"),
    ("TP53", "gene_name_exact", None),
    ("ENSG00000141510", "gene_id_exact", None),
    ("CCDS11118.1", "alias_exact", "ENST00000269305.9"),
    ("TP5", "prefix", None),               # autocomplete prefix
    ("0000026930", "trigram", "ENST00000269305.9"),  # substring of the ID digits
])
def test_search_tiers(gtf_db: Path, query, expect_field, expect_tid):
    con = si.open_readonly(gtf_db)
    res = si.search(con, query, limit=10)
    con.close()
    assert res, f"no results for {query!r}"
    assert res[0]["matched_field"] == expect_field
    if expect_tid is not None:
        assert any(r["transcript_id"] == expect_tid for r in res)


def test_search_empty(gtf_db: Path):
    con = si.open_readonly(gtf_db)
    assert si.search(con, "   ", limit=5) == []
    assert si.search(con, "ZZZZNOTAGENE", limit=5) == []
    con.close()


# ── GFF3 parsing parity ───────────────────────────────────────────────────────

def test_gff3_aliases_and_utr(tmp_path: Path):
    src = tmp_path / "annot.gff3"
    src.write_text(GFF3)
    out = tmp_path / "idx.sqlite"
    si.build_sqlite_index(src, out, source_name="Ensembl", force=True)
    info = si.inspect(out)
    assert info["metadata"]["input_format"] == "GFF3"

    con = si.open_readonly(out)
    # Alias from Alias= and Dbxref=RefSeq:
    assert si.search(con, "p53", limit=5)
    refseq = si.search(con, "NM_000546.6", limit=5)
    assert refseq and refseq[0]["matched_field"] == "alias_exact"
    # explicit five/three_prime_UTR on minus strand
    row = con.execute(
        "SELECT utr5_start, utr5_end, utr3_start, utr3_end FROM feature "
        "WHERE transcript_id='ENST00000269305.9'").fetchone()
    assert row == (451, 500, 100, 149)
    con.close()


# ── synonym enrichment (HGNC) ─────────────────────────────────────────────────

POU5F1_GTF = (
    'chr6\tHAVANA\tgene\t31164337\t31180731\t.\t-\t.\t'
    'gene_id "ENSG00000204531.21"; gene_type "protein_coding"; gene_name "POU5F1";\n'
    'chr6\tHAVANA\ttranscript\t31164337\t31170682\t.\t-\t.\t'
    'gene_id "ENSG00000204531.21"; transcript_id "ENST00000259915.13"; '
    'gene_name "POU5F1"; transcript_name "POU5F1-201";\n'
    'chr6\tHAVANA\texon\t31164337\t31170682\t.\t-\t.\t'
    'gene_id "ENSG00000204531.21"; transcript_id "ENST00000259915.13";\n'
)
HGNC_TSV = (
    "symbol\talias_symbol\tprev_symbol\tensembl_gene_id\n"
    'POU5F1\t"OCT3|Oct4|OCT-4"\tOTF3\tENSG00000204531\n'
)


@pytest.mark.parametrize("query", ["oct4", "OCT4", "OCT-4", "OTF3", "oct3"])
def test_synonym_injection(tmp_path: Path, query: str):
    src = tmp_path / "annot.gtf"
    src.write_text(POU5F1_GTF)
    syn = tmp_path / "hgnc.tsv"
    syn.write_text(HGNC_TSV)
    out = tmp_path / "idx.sqlite"
    si.build_sqlite_index(src, out, synonyms=syn, force=True)
    # metadata records the injection
    assert int(si.inspect(out)["metadata"]["n_synonyms_injected"]) == 4
    con = si.open_readonly(out)
    res = si.search(con, query, limit=5)
    con.close()
    assert res, f"{query!r} did not resolve"
    assert res[0]["gene_name"] == "POU5F1"
    assert res[0]["matched_field"] == "alias_exact"


def test_no_synonyms_no_oct4(tmp_path: Path):
    """Without the synonym feed, OCT4 is absent (it is not in the annotation)."""
    src = tmp_path / "annot.gtf"
    src.write_text(POU5F1_GTF)
    out = tmp_path / "idx.sqlite"
    si.build_sqlite_index(src, out, force=True)
    con = si.open_readonly(out)
    assert si.search(con, "OCT4", limit=5) == []
    assert si.search(con, "POU5F1", limit=5)        # the real symbol still works
    con.close()


# ── gzip input ────────────────────────────────────────────────────────────────

def test_gzip_input(tmp_path: Path):
    src = tmp_path / "annot.gtf.gz"
    with gzip.open(src, "wt") as fh:
        fh.write(GTF)
    out = tmp_path / "idx.sqlite"
    si.build_sqlite_index(src, out, force=True)
    info = si.inspect(out)
    assert info["n_transcripts"] == 2


# ── determinism ───────────────────────────────────────────────────────────────

def test_deterministic_ids(tmp_path: Path):
    src = tmp_path / "annot.gtf"
    src.write_text(GTF)
    a = tmp_path / "a.sqlite"
    b = tmp_path / "b.sqlite"
    si.build_sqlite_index(src, a, force=True)
    si.build_sqlite_index(src, b, force=True)
    con_a = sqlite3.connect(a)
    con_b = sqlite3.connect(b)
    rows_a = con_a.execute(
        "SELECT id, feature_type, gene_id, transcript_id FROM feature ORDER BY id"
    ).fetchall()
    rows_b = con_b.execute(
        "SELECT id, feature_type, gene_id, transcript_id FROM feature ORDER BY id"
    ).fetchall()
    assert rows_a == rows_b
    con_a.close()
    con_b.close()
