# Installation

## Install the Python package

```bash
pip install refbox
```

Requires **Python ≥ 3.9**. The package pulls in `pyyaml`, `requests`, `tqdm`,
and `pysam`.

### From source (development)

```bash
git clone https://github.com/typekey/refbox.git
cd refbox
pip install -e .
pytest -q          # unit tests for `refbox build` single-file modes
refbox --help
```

## External command-line tools

`refbox` shells out to standard bioinformatics tools. Install them once with
conda / mamba:

```bash
mamba install -c bioconda htslib samtools gffread ucsc-bedtobigbed
```

| Tool | Used for | Required for |
|---|---|---|
| `bgzip` | block-gzip compression | every build |
| `tabix` | positional indexing of GFF/BED | every annotation/BED build |
| `samtools` | `faidx` on FASTA | genome / transcriptome |
| `gffread` | extract transcripts from genome + GTF | derived transcriptome |
| `liftOver` | RNAcentral cross-assembly liftover | assemblies using `liftover_from` |
| `bedToBigBed` / `bigBedToBed` | bigBed conversion / expansion | `refbox build -bed`, cytoband |
| GNU `sort`, `grep` | sorting / filtering | every build |
