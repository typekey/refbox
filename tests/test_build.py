"""Tests for ``refbox build`` single-file modes.

Each test starts from a tiny synthetic genome (3 chromosomes, 200 bp each)
and a matching GTF / BED / rmsk table, then verifies the produced indexed
artifacts can be read back with samtools/tabix.
"""

from __future__ import annotations

import gzip
import shutil
import subprocess
from pathlib import Path

import pytest

from refbox import file_build as fb
from refbox.utils import _is_bgzip


GENOME = (
    ">chr1\n" + ("A" * 60 + "\n") * 3 + "A" * 20 + "\n"
    ">chr2\n" + ("C" * 60 + "\n") * 3 + "C" * 20 + "\n"
    ">chr3\n" + ("G" * 60 + "\n") * 3 + "G" * 20 + "\n"
)

# Two transcripts on chr1, one on chr2.  Coordinates are 1-based, inclusive.
GTF = (
    "##sample\n"
    'chr1\ttest\ttranscript\t10\t100\t.\t+\t.\tgene_id "g1"; transcript_id "t1";\n'
    'chr1\ttest\texon\t10\t60\t.\t+\t.\tgene_id "g1"; transcript_id "t1";\n'
    'chr1\ttest\texon\t80\t100\t.\t+\t.\tgene_id "g1"; transcript_id "t1";\n'
    'chr2\ttest\ttranscript\t1\t150\t.\t-\t.\tgene_id "g2"; transcript_id "t2";\n'
    'chr2\ttest\texon\t1\t150\t.\t-\t.\tgene_id "g2"; transcript_id "t2";\n'
)

# Unsorted BED on purpose to exercise the sort path.
BED_UNSORTED = (
    "chr2\t10\t50\tfeatB\t100\t+\n"
    "chr1\t100\t200\tfeatA1\t50\t-\n"
    "chr1\t5\t15\tfeatA0\t10\t+\n"
    "chr3\t0\t90\tfeatC\t999\t-\n"
)

# UCSC rmsk.txt schema (16 cols).  We only fill what build_rmsk reads.
RMSK_ROWS = [
    # bin swScore milliDiv milliDel milliIns genoName genoStart genoEnd
    # genoLeft strand repName repClass repFamily repStart repEnd repLeft id
    ("0", "1234", "0", "0", "0", "chr1", "20", "40", "0", "+",
     "L1HS", "LINE", "L1", "0", "20", "0", "1"),
    ("0", "999",  "0", "0", "0", "chr1", "60", "80", "0", "-",
     "AluY", "SINE", "Alu", "0", "20", "0", "2"),
    ("0", "500",  "0", "0", "0", "chr2", "0",  "30", "0", "+",
     "MIR",  "SINE", "MIR", "0", "30", "0", "3"),
]


@pytest.fixture
def workspace(tmp_path: Path) -> Path:
    g = tmp_path / "genome.fa"
    g.write_text(GENOME)
    a = tmp_path / "annot.gtf"
    a.write_text(GTF)
    b = tmp_path / "feats.bed"
    b.write_text(BED_UNSORTED)
    r = tmp_path / "rmsk.txt"
    r.write_text("\n".join("\t".join(row) for row in RMSK_ROWS) + "\n")
    # also gzipped rmsk to exercise the .gz path
    rgz = tmp_path / "rmsk.txt.gz"
    with gzip.open(rgz, "wt") as fh:
        fh.write("\n".join("\t".join(row) for row in RMSK_ROWS) + "\n")
    return tmp_path


# ── build_fa ──────────────────────────────────────────────────────────────────

def test_build_fa(workspace: Path):
    out = fb.build_fa(workspace / "genome.fa")
    assert out.exists() and out.stat().st_size > 0
    assert _is_bgzip(out), "output FASTA must be bgzip-compressed"
    assert Path(str(out) + ".fai").exists()
    assert Path(str(out) + ".gzi").exists()
    cs = out.parent / (out.stem + ".chrom.sizes")
    assert cs.exists()
    chroms = [l.split("\t")[0] for l in cs.read_text().splitlines()]
    assert chroms == ["chr1", "chr2", "chr3"]

    # samtools faidx must be able to slice a region
    r = subprocess.run(
        ["samtools", "faidx", str(out), "chr1:1-10"],
        capture_output=True, text=True, check=True,
    )
    assert "AAAAAAAAAA" in r.stdout


# ── build_gxf (GTF) ───────────────────────────────────────────────────────────

def test_build_gxf_unsorted(workspace: Path):
    # write an unsorted GTF to verify the sort branch runs
    src = workspace / "annot_unsorted.gtf"
    src.write_text(
        'chr2\ttest\texon\t1\t150\t.\t-\t.\tgene_id "g2"; transcript_id "t2";\n'
        'chr1\ttest\texon\t10\t60\t.\t+\t.\tgene_id "g1"; transcript_id "t1";\n'
    )
    out = fb.build_gxf(src)
    assert out.exists() and _is_bgzip(out)
    assert Path(str(out) + ".tbi").exists()
    r = subprocess.run(["tabix", str(out), "chr1:1-1000"],
                       capture_output=True, text=True, check=True)
    assert "g1" in r.stdout


# ── build_bed ─────────────────────────────────────────────────────────────────

def test_build_bed(workspace: Path):
    cs_path = workspace / "chrom.sizes"
    cs_path.write_text("chr1\t200\nchr2\t200\nchr3\t200\n")
    out = fb.build_bed(workspace / "feats.bed", chrom_sizes=cs_path)
    assert out.exists() and _is_bgzip(out)
    assert Path(str(out) + ".tbi").exists()
    bb = out.with_suffix("").with_suffix(".bigBed")
    assert bb.exists() and bb.stat().st_size > 0, "bigBed should be created"
    # tabix query
    r = subprocess.run(["tabix", str(out), "chr1:1-300"],
                       capture_output=True, text=True, check=True)
    assert "featA0" in r.stdout and "featA1" in r.stdout


# ── build_rmsk (plain + .gz) ──────────────────────────────────────────────────

@pytest.mark.parametrize("name", ["rmsk.txt", "rmsk.txt.gz"])
def test_build_rmsk(workspace: Path, name: str):
    out = fb.build_rmsk(workspace / name, out_dir=workspace / f"rmsk_out_{name}")
    bed_gz = out["bed"]
    gtf_gz = out["gtf"]
    for f in (bed_gz, gtf_gz):
        assert f.exists() and _is_bgzip(f)
        assert Path(str(f) + ".tbi").exists()
    r = subprocess.run(["tabix", str(bed_gz), "chr1:1-100"],
                       capture_output=True, text=True, check=True)
    assert "L1HS" in r.stdout
    assert "AluY" in r.stdout


# ── build_transcriptome (genome + gtf -> transcripts) ─────────────────────────

def test_build_transcriptome(workspace: Path):
    # gffread needs a .fai for the genome — produce one first.
    g = workspace / "genome.fa"
    subprocess.run(["samtools", "faidx", str(g)], check=True)
    out = fb.build_transcriptome(g, workspace / "annot.gtf",
                                 workspace / "transcriptome.fa.gz")
    assert out.exists() and _is_bgzip(out)
    assert Path(str(out) + ".fai").exists()
    fai_lines = Path(str(out) + ".fai").read_text().strip().splitlines()
    names = [l.split("\t")[0] for l in fai_lines]
    # two transcripts defined in the GTF; headers are GENCODE-style pipe-delimited
    tids = {n.split("|")[0] for n in names}
    assert "t1" in tids and "t2" in tids


# ── auto-detect ───────────────────────────────────────────────────────────────

@pytest.mark.parametrize("fname,expected", [
    ("genome.fa", "fa"),
    ("annot.gtf", "gxf"),
    ("annot.gff3", "gxf"),
    ("feats.bed", "bed"),
    ("rmsk.txt", "rmsk"),
    ("rmsk.txt.gz", "rmsk"),
])
def test_auto_detect(workspace: Path, fname: str, expected: str):
    assert fb.auto_detect(workspace / fname) == expected


def test_auto_detect_unknown(tmp_path: Path):
    p = tmp_path / "mystery.dat"
    p.write_text("x")
    with pytest.raises(ValueError):
        fb.auto_detect(p)


# ── CLI integration smoke test ────────────────────────────────────────────────

def test_cli_build_fa(workspace: Path, monkeypatch):
    from refbox.cli import main
    out = workspace / "out.fa.gz"
    rc = main(["build", "-fa", str(workspace / "genome.fa"), "-o", str(out)])
    assert rc == 0
    assert out.exists() and _is_bgzip(out)
