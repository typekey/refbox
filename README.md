# refbox

[![PyPI version](https://img.shields.io/pypi/v/refbox.svg)](https://pypi.org/project/refbox/)
[![Python versions](https://img.shields.io/pypi/pyversions/refbox.svg)](https://pypi.org/project/refbox/)
[![Tests](https://github.com/typekey/refbox/actions/workflows/tests.yml/badge.svg)](https://github.com/typekey/refbox/actions/workflows/tests.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

**Build standardized, indexed reference files for genome browsers — in one command.**

`refbox` turns a YAML registry of species/assemblies into ready-to-load browser
inputs:

- Genome FASTA → `bgzip` + `samtools faidx` (`.fa.gz` + `.fai` + `.gzi`) + `chrom.sizes`
- Transcriptome FASTA → `bgzip` + `samtools faidx` (auto-derived from genome + GTF via `gffread` when no upstream URL exists)
- GTF / GFF3 annotations → sorted + `bgzip` + `tabix`
- Repeats (UCSC RepeatMasker `rmsk.txt.gz` → BED + GTF; `.fa.out.gz` report)
- RNAcentral non-coding RNA annotations (direct download _or_ liftover from another assembly)
- ENCODE SCREEN cCREs

It ships a registry of **26 species / 42 assemblies** (human, mouse, rat, dog,
cow, pig, chimp, gorilla, zebrafish, fly, worm, sea urchin, yeast, plants,
bacteria, viruses…) covering GENCODE, Ensembl, Ensembl Genomes, UCSC golden
path, NCBI, RNAcentral, and ENCODE SCREEN.

---

## Install

```bash
pip install refbox
```

External tools (install once via conda / mamba):

```bash
mamba install -c bioconda htslib samtools gffread ucsc-bedtobigbed
```

Required: `bgzip`, `tabix`, `samtools`, `gffread`, `liftOver` (for RNAcentral
cross-assembly), `bedToBigBed` (for `refbox build -bed`), GNU `sort`/`grep`.

---

## CLI overview

```
refbox download   # only fetch raw files configured in species.yaml
refbox pull       # full pipeline: download (if missing) + build + test
refbox test       # validate build/ outputs
refbox build      # single-file / directory build for arbitrary inputs
```

### `refbox pull` — the registry-driven pipeline

```bash
# Fetch + build + validate Human GRCh38 (uses bundled species.yaml).
refbox pull --species Homo_sapiens --assembly GRCh38

# --species is optional; it is inferred from --assembly via the registry.
refbox pull --assembly GRCm38

# Loop every assembly in the registry, including ones marked enabled: false.
refbox pull --include-disabled --resource genome transcriptome \
            annotation_gtf annotation_gff3 repeats_rmsk rnacentral

# Pull only specific resources.
refbox pull --assembly GRCh38 --resource ccre
```

| Flag | Meaning |
|---|---|
| `--species` | optional filter; inferred from `--assembly` when omitted |
| `--assembly` | optional filter; omit to run every assembly that matches `--species` |
| `--resource` | subset of `genome transcriptome annotation_gtf annotation_gff3 repeats_rmsk repeats_bed repeats_gtf repeats_fa rnacentral ccre` |
| `--include-disabled` | also process assemblies with `enabled: false` |
| `--out DIR` | output root (default: `$REFBOX_OUT` or `$PWD`) |
| `--force` | rebuild even when outputs exist |
| `--no-download` | skip the auto-download phase |
| `--no-test` | skip the post-build validation |
| `-v` | verbose / DEBUG logging |

### `refbox build` — single-file / directory build

For custom data that does not have a `species.yaml` entry. Each call detects
the input by extension (or by explicit flag), verifies bgzip / sort order,
runs the canonical transformation, and emits indexed outputs.

```bash
# Genome FASTA → bgzipped + faidx + chrom.sizes
refbox build -fa  GENOME.fa [-o OUT.fa.gz]

# GTF or GFF3 → sorted + bgzip + tabix
refbox build -gtf ANNOT.gtf [-o OUT.gtf.gz]
refbox build -gff ANNOT.gff3 [-o OUT.gff3.gz]

# Genome + annotation → transcripts FASTA (via gffread) + faidx
refbox build -fa GENOME.fa -gtf ANNOT.gtf -o transcripts.fa.gz

# BED → sorted + bgzip + tabix + bigBed
#   chrom.sizes resolved from --chrom-sizes FILE or --assembly NAME
#   (the latter delegates to zlbio's per-species lookup).
refbox build -bed FEATURES.bed [--chrom-sizes FILE | --assembly NAME]

# UCSC rmsk.txt[.gz] → repeats.sorted.bed.gz + repeats.sorted.gtf.gz (+ .tbi)
refbox build -rmsk rmsk.txt.gz [-o OUT_DIR]

# Directory of user files (auto-classified by extension) → full layout under {Species}/{Assembly}/
refbox build -i DIR --assembly NAME [--species NAME]

# Auto-detect a single file by extension
refbox build SOMEFILE
```

Inputs may be plain or `.gz`. The builder transparently re-bgzips a plain
gzip, copies an already-bgzipped file as-is, and re-sorts an unsorted GFF/BED
before indexing.

### `refbox download` / `refbox test`

```bash
refbox download --assembly GRCh38                  # raw files only
refbox test     --assembly GRCh38                  # validate existing build/
refbox test     --include-disabled                 # everything in the registry
```

---

## Output layout

```
{REFBOX_OUT}/
  {Species}/
    {Assembly}/
      raw/                          # original downloads / copies
        genome.fa
        transcriptome.fa             # may be derived from genome + GTF
        annotation_gtf.gtf
        annotation_gff3.gff3
        repeats_rmsk.tsv
        repeats_fa.fa
        rnacentral.gff3              # may be lifted from another assembly
        ccre.bed
      build/                        # browser-loadable
        genome.fa.gz                 + .fai + .gzi
        chrom.sizes
        transcripts.fa.gz            + .fai
        annotation.sorted.gtf.gz     + .tbi
        annotation.sorted.gff3.gz    + .tbi
        repeats.sorted.bed.gz        + .tbi
        repeats.sorted.gtf.gz        + .tbi
        rnacentral.sorted.gff3.gz    + .tbi
        ccre.sorted.bed.gz           + .tbi
```

---

## Programmatic API

```python
from refbox.config   import load_config, iter_targets
from refbox.download import download_targets
from refbox.build    import build_targets
from refbox.test     import test_targets
from refbox          import file_build as fb     # single-file builders
from refbox.ingest   import ingest_directory     # directory ingest
from refbox.report   import build_report         # Markdown status report
```

---

## Tutorial — adding a new assembly

Edit [`config/species.yaml`](config/species.yaml). Each entry is **three
levels deep**: species → assembly → resource. All 10 resource keys must
appear; use `null` for ones with no upstream source.

```yaml
species:
  Homo_sapiens:
    GRCh38:
      enabled: true
      gencode_version: 44
      ucsc_db: hg38
      genome:
        url: https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/GRCh38.primary_assembly.genome.fa.gz
      transcriptome:
        url: https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.transcripts.fa.gz
      annotation_gtf:
        url: https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.annotation.gtf.gz
      annotation_gff3:
        url: https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.annotation.gff3.gz
      repeats_rmsk:
        url: https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/rmsk.txt.gz
      repeats_bed: null
      repeats_gtf: null
      repeats_fa:
        url: https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.out.gz
      rnacentral:
        url: https://ftp.ebi.ac.uk/pub/databases/RNAcentral/current_release/genome_coordinates/gff3/homo_sapiens.GRCh38.gff3.gz
      ccre:
        url: https://downloads.wenglab.org/V3/GRCh38-cCREs.bed
```

### Fallbacks

```yaml
my_resource:
  local_path: /path/on/disk/file.fa.gz     # copy if present
  url:        https://.../file.fa.gz       # else download
```

```yaml
# RNAcentral cross-assembly liftover (no direct URL upstream):
rnacentral:
  liftover_from:
    source_assembly: GRCh38
    url:       https://.../homo_sapiens.GRCh38.gff3.gz
    chain_url: https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz
```

```yaml
# Transcriptome auto-derivation: leave it null and refbox will build
# transcripts.fa.gz from genome + GTF/GFF via gffread.
transcriptome: null
```

### Canonical resource names

| Name | raw/ | build/ |
|---|---|---|
| `genome` | `genome.fa` | `genome.fa.gz` + `.fai` + `.gzi`, `chrom.sizes` |
| `transcriptome` | `transcriptome.fa` | `transcripts.fa.gz` + `.fai` |
| `annotation_gtf` | `annotation_gtf.gtf` | `annotation.sorted.gtf.gz` + `.tbi` |
| `annotation_gff3` | `annotation_gff3.gff3` | `annotation.sorted.gff3.gz` + `.tbi` |
| `repeats_rmsk` | `repeats_rmsk.tsv` | `repeats.sorted.bed.gz` + `repeats.sorted.gtf.gz` |
| `repeats_bed` | `repeats_bed.bed` | `repeats.sorted.bed.gz` + `.tbi` |
| `repeats_gtf` | `repeats_gtf.gtf` | `repeats.sorted.gtf.gz` + `.tbi` |
| `repeats_fa` | `repeats_fa.fa` | _(RepeatMasker `.fa.out` report)_ |
| `rnacentral` | `rnacentral.gff3` | `rnacentral.sorted.gff3.gz` + `.tbi` |
| `ccre` | `ccre.bed` | `ccre.sorted.bed.gz` + `.tbi` |

---

## Environment

| Variable | Meaning |
|---|---|
| `REFBOX_OUT`    | default output root for `{Species}/{Assembly}/{raw,build}/` |
| `REFBOX_CONFIG` | path to a custom `species.yaml` (overrides the bundled registry) |

---

## Status reporter

```bash
python -m refbox.report --out /path/to/reference > report.md
```

Walks the output tree and emits a Chinese-language Markdown report listing,
per assembly, the status (✓ done / ⚠ missing index / ✗ missing) and size of
each artifact.

---

## Development

```bash
git clone https://github.com/typekey/refbox.git
cd refbox
pip install -e .
pytest -q                 # unit tests for `refbox build` single-file modes
refbox --help
```

### Release

Tags matching `v*` automatically build and publish to PyPI via GitHub Actions
([`.github/workflows/workflow.yml`](.github/workflows/workflow.yml)) using
[PyPI trusted publishing](https://docs.pypi.org/trusted-publishers/).

```bash
git tag v0.3.0
git push origin v0.3.0
```

---

## Changelog

- **v0.3.0** — CLI refactor: `build`→`pull`; new `refbox build` for arbitrary
  single-file/directory inputs (with auto bgzip / sort / bigBed); transcriptome
  auto-derivation via `gffread`; unit-tested + CI; `--include-disabled` flag;
  Chinese `refbox.report` generator.
- **v0.2.0** — `refbox import` subcommand; optional `--species`; `refbox build`
  auto-downloads missing raws and chains into `refbox test`; full canonical
  resource fields for all 42 assemblies; RNAcentral liftover entries for
  GRCh37 / mm9 / rn6.
- **v0.1.5** — RNAcentral `liftover_from` (download chain + source GFF,
  liftOver + chrom-name normalization) for assemblies without upstream files.
- **v0.1.4** — RNAcentral chrom-name normalization (Ensembl → UCSC).
- **v0.1.3** — atomic bgzip writes (`.tmp` + rename) to prevent truncated outputs.

## License

MIT — see [LICENSE](LICENSE).
