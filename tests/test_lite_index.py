"""Tests for the lightweight RBrowser Index builder (``refbox.lite_index``).

The ``.rbi`` lite index is a compact SQLite B-tree lookup (no FTS5/trigram, no
exon/CDS structure). These tests run on pure ``sqlite3`` — no external tools —
and exercise parsing, the ``term`` table, the ranked tiered search, the optional
3-gram fuzzy recall, and the ``.rbi`` extension / format metadata.

They also assert the split between the two index families:
    refbox build -rba  → full RBrowser Annotation index (.rba, FTS5 + structure)
    refbox build -rbi  → lightweight RBrowser Index    (.rbi, B-tree lookup)
"""

from __future__ import annotations

import gzip
import sqlite3
from pathlib import Path

import pytest

from refbox import lite_index as li


# Same two-gene shape as the .rba tests: a minus-strand protein-coding gene and a
# plus-strand lncRNA, with versioned Ensembl IDs and hyphenated transcript names.
GTF = """\
##provider: TEST
chr17\tHAVANA\tgene\t100\t500\t.\t-\t.\tgene_id "ENSG00000141510.18"; gene_type "protein_coding"; gene_name "TP53";
chr17\tHAVANA\ttranscript\t100\t500\t.\t-\t.\tgene_id "ENSG00000141510.18"; transcript_id "ENST00000269305.9"; gene_name "TP53"; transcript_type "protein_coding"; transcript_name "TP53-201";
chr17\tHAVANA\texon\t400\t500\t.\t-\t.\tgene_id "ENSG00000141510.18"; transcript_id "ENST00000269305.9";
chr17\tHAVANA\texon\t100\t200\t.\t-\t.\tgene_id "ENSG00000141510.18"; transcript_id "ENST00000269305.9";
chr11\tENSEMBL\tgene\t1000\t2000\t.\t+\t.\tgene_id "ENSG00000251562.8"; gene_type "lncRNA"; gene_name "MALAT1";
chr11\tENSEMBL\ttranscript\t1000\t2000\t.\t+\t.\tgene_id "ENSG00000251562.8"; transcript_id "ENST00000534336.1"; gene_name "MALAT1"; transcript_type "lncRNA"; transcript_name "MALAT1-201";
chr11\tENSEMBL\texon\t1000\t2000\t.\t+\t.\tgene_id "ENSG00000251562.8"; transcript_id "ENST00000534336.1";
"""

GFF3 = """\
##gff-version 3
chr17\tHAVANA\tgene\t100\t500\t.\t-\t.\tID=ENSG00000141510.18;gene_id=ENSG00000141510.18;gene_type=protein_coding;gene_name=TP53
chr17\tHAVANA\tmRNA\t100\t500\t.\t-\t.\tID=ENST00000269305.9;Parent=ENSG00000141510.18;transcript_id=ENST00000269305.9;gene_name=TP53;transcript_name=TP53-201
chr17\tHAVANA\texon\t400\t500\t.\t-\t.\tID=e1;Parent=ENST00000269305.9
chr17\tHAVANA\texon\t100\t200\t.\t-\t.\tID=e2;Parent=ENST00000269305.9
"""


@pytest.fixture
def rbi(tmp_path: Path) -> Path:
    src = tmp_path / "annot.gtf"
    src.write_text(GTF)
    out = tmp_path / "idx.rbi"
    li.build_lite_index(src, out, source_name="GENCODE", species="human",
                        genome="hg38", annotation_version="vTEST",
                        enable_gram3=True, force=True)
    return out


# ── normalization helpers ───────────────────────────────────────────────────────

@pytest.mark.parametrize("raw,expected", [
    ("TP53", "tp53"),
    ("TP53-201", "tp53-201"),                 # separators preserved in the main form
    ("  MALAT1 ", "malat1"),
])
def test_normalize(raw: str, expected: str):
    assert li.normalize(raw) == expected


def test_strip_version_and_sepfree():
    assert li.strip_version("enst00000269305.9") == "enst00000269305"
    assert li.strip_version("tp53") == "tp53"
    assert li.sepfree("tp53-201") == "tp53201"
    assert li.sepfree("protein_coding") == "proteincoding"


def test_term_forms_dedup():
    forms = li.term_forms("ENST00000269305.9")
    # normalized, version-stripped, separator-free — distinct, order preserved
    assert forms[0] == "enst00000269305.9"
    assert "enst00000269305" in forms


def test_prefix_bounds():
    lo, hi = li.prefix_bounds("tp")
    assert lo == "tp" and hi == "tq"


# ── build / schema / metadata ───────────────────────────────────────────────────

def test_extension_and_default_name(tmp_path: Path):
    """Default output is ``<stem>.rbi`` (NOT .rbai)."""
    src = tmp_path / "gencode.v45.annotation.gtf"
    src.write_text(GTF)
    out = li.build_lite_index(src, force=True)
    assert out.name == "gencode.v45.annotation.rbi"
    assert out.suffix == ".rbi"


def test_counts_and_format_metadata(rbi: Path):
    info = li.inspect(rbi)
    assert info["n_genes"] == 2
    assert info["n_transcripts"] == 2
    assert info["n_terms"] > 0
    meta = info["metadata"]
    # the format identity is RBI / "RBrowser Index" (lightweight), not RBAI/.rbai
    assert meta["format_short_name"] == "RBI"
    assert meta["format_name"] == "RBrowser Index"
    assert "rbai" not in meta["format_short_name"].lower()
    assert meta["source_name"] == "GENCODE"
    assert meta["input_format"] == "GTF"
    assert meta["index_type"] == "lite_btree_lookup"


def test_no_fts_or_trigram_tables(rbi: Path):
    """The lite index must NOT carry the heavy FTS5/trigram tables — that is the
    whole point of ``.rbi`` vs the full ``.rba``."""
    con = sqlite3.connect(rbi)
    tables = {r[0] for r in con.execute(
        "SELECT name FROM sqlite_master WHERE type='table'")}
    con.close()
    assert "record" in tables and "term" in tables
    assert "feature_fts" not in tables
    assert "feature_trigram" not in tables


def test_record_holds_position(rbi: Path):
    con = sqlite3.connect(rbi)
    row = con.execute(
        "SELECT chrom, start, end, strand FROM record "
        "WHERE transcript_id='ENST00000269305.9'").fetchone()
    con.close()
    assert row == ("chr17", 100, 500, "-")


# ── ranked tiered search ─────────────────────────────────────────────────────────

@pytest.mark.parametrize("query,mode,expect_name", [
    ("TP53", "exact", "TP53"),                       # gene name exact
    ("ENST00000269305.9", "exact", "TP53"),          # transcript_id exact
    ("ENST00000269305", "exact", "TP53"),            # versionless → version-stripped term
    ("TP53-201", "exact", "TP53"),                   # transcript_name exact
    ("MALAT1", "exact", "MALAT1"),
])
def test_search_exact_tiers(rbi: Path, query, mode, expect_name):
    con = li.open_readonly(rbi)
    got_mode, res = li.search(con, query, limit=10)
    con.close()
    assert res, f"no results for {query!r}"
    assert got_mode == mode
    assert any((r["gene_name"] == expect_name) or
               (r["transcript_name"] == expect_name) or
               (r["gene_name"] == expect_name) for r in res)


def test_search_prefix(rbi: Path):
    con = li.open_readonly(rbi)
    mode, res = li.search(con, "TP5", limit=10)        # autocomplete
    con.close()
    assert mode == "prefix"
    assert any(r["gene_name"] == "TP53" for r in res)


def test_search_gram3_fuzzy(rbi: Path):
    """With gram3 enabled, an interior name substring recalls via the 3-gram tier
    (no exact/prefix term matches 'ALAT')."""
    con = li.open_readonly(rbi)
    mode, res = li.search(con, "ALAT", limit=10)
    con.close()
    assert mode == "gram3"
    assert any(r["gene_name"] == "MALAT1" for r in res)


def test_search_empty_and_miss(rbi: Path):
    con = li.open_readonly(rbi)
    assert li.search(con, "   ", limit=5) == ("none", [])
    assert li.search(con, "ZZZNOTAGENE", limit=5) == ("none", [])
    con.close()


# ── gram3 toggle ─────────────────────────────────────────────────────────────────

def test_no_gram3_smaller_and_no_fuzzy(tmp_path: Path):
    src = tmp_path / "a.gtf"
    src.write_text(GTF)
    out = tmp_path / "nogram.rbi"
    li.build_lite_index(src, out, enable_gram3=False, force=True)
    info = li.inspect(out)
    assert info["n_gram3"] == 0
    assert info["metadata"]["gram3_enabled"] == "0"
    con = li.open_readonly(out)
    # exact/prefix still work; interior-substring fuzzy no longer recalls
    assert li.search(con, "MALAT1", limit=5)[1]
    assert li.search(con, "ALAT", limit=5) == ("none", [])
    con.close()


# ── GFF3 parity + gzip input ─────────────────────────────────────────────────────

def test_gff3_input(tmp_path: Path):
    src = tmp_path / "annot.gff3"
    src.write_text(GFF3)
    out = tmp_path / "idx.rbi"
    li.build_lite_index(src, out, force=True)
    info = li.inspect(out)
    assert info["metadata"]["input_format"] == "GFF3"
    con = li.open_readonly(out)
    assert li.search(con, "TP53", limit=5)[1]
    con.close()


def test_gzip_input(tmp_path: Path):
    src = tmp_path / "annot.gtf.gz"
    with gzip.open(src, "wt") as fh:
        fh.write(GTF)
    out = tmp_path / "idx.rbi"
    li.build_lite_index(src, out, force=True)
    assert li.inspect(out)["n_transcripts"] == 2


def test_force_required_to_overwrite(tmp_path: Path):
    src = tmp_path / "a.gtf"
    src.write_text(GTF)
    out = tmp_path / "idx.rbi"
    li.build_lite_index(src, out, force=True)
    with pytest.raises(FileExistsError):
        li.build_lite_index(src, out, force=False)


# ── CLI dispatch: -rba (full) vs -rbi (lite) are distinct builders ───────────────

def test_cli_build_rbi(tmp_path: Path):
    from refbox.cli import main
    src = tmp_path / "annot.gtf"
    src.write_text(GTF)
    out = tmp_path / "out.rbi"
    rc = main(["build", "-rbi", str(src), "-o", str(out),
               "--source-name", "TEST", "--genome", "hg38", "--force"])
    assert rc == 0
    assert out.exists()
    assert li.inspect(out)["metadata"]["format_short_name"] == "RBI"


def test_cli_build_rba_is_full_index(tmp_path: Path):
    """``-rba`` routes to the full SQLite+FTS5 builder (has FTS tables); ``-rbi``
    routes to the lite builder (no FTS tables). Proves they are not aliases."""
    from refbox.cli import main
    from refbox import sqlite_index as si
    src = tmp_path / "annot.gtf"
    src.write_text(GTF)
    rba = tmp_path / "out.rba"
    rc = main(["build", "-rba", str(src), "-o", str(rba),
               "--source-name", "TEST", "--genome", "hg38", "--force"])
    assert rc == 0
    info = si.inspect(rba)
    assert info["has_fts"]                       # full index carries FTS5
    con = si.open_readonly(rba)
    assert si.search(con, "TP53", limit=5)[0]["matched_field"] == "gene_name_exact"
    con.close()


def test_cli_with_rbi_alongside_tabix(tmp_path: Path):
    """``refbox build -gtf … --with-rbi`` emits a sorted/bgzip/tabix GTF AND a
    lightweight .rbi next to it."""
    from refbox.cli import main
    src = tmp_path / "annot.gtf"
    src.write_text(GTF)
    rc = main(["build", "-gtf", str(src), "--with-rbi",
               "--source-name", "TEST", "--genome", "hg38", "--force"])
    assert rc == 0
    rbi_out = tmp_path / "annot.rbi"
    assert rbi_out.exists()
    assert li.inspect(rbi_out)["metadata"]["format_short_name"] == "RBI"
