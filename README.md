# refbox

[![PyPI version](https://img.shields.io/pypi/v/refbox.svg)](https://pypi.org/project/refbox/)
[![Python versions](https://img.shields.io/pypi/pyversions/refbox.svg)](https://pypi.org/project/refbox/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

**Build standardized, indexed reference files for genome browsers — in one command.**

`refbox` turns a YAML registry of species/assemblies into ready-to-load browser
inputs:

- Genome FASTA → `bgzip` + `samtools faidx` (`.gz` + `.fai` + `.gzi`)
- Transcriptome FASTA → `bgzip` + `samtools faidx`
- GTF / GFF3 annotations → sorted + `bgzip` + `tabix`
- Repeats (UCSC RepeatMasker `rmsk.txt.gz` + `.fa.out.gz`)
- RNAcentral non-coding RNA annotations
- ENCODE SCREEN cCREs
- `chrom.sizes` derived from the genome `.fai`

It ships a registry of **26 species / 42 assemblies** (human, mouse, rat, dog,
cow, pig, chimp, gorilla, zebrafish, fly, worm, sea urchin, yeast, plants,
bacteria, viruses…) covering GENCODE, Ensembl, Ensembl Genomes, UCSC golden
path, NCBI, RNAcentral, and ENCODE SCREEN.

---

## Why

Genome browsers (e.g. [WashU Epigenome Browser](https://epigenomegateway.wustl.edu/),
[rbrowser](https://github.com/typekey/rbrowser)) expect a strict set of file
formats and indices. Manually downloading, sorting, bgzipping, tabixing, and
checking each one is tedious and error-prone. `refbox` makes the process:

- **declarative** — one YAML, every URL pinned
- **idempotent** — re-running skips finished files; `--force` to rebuild
- **filterable** — `--species`, `--assembly`, `--resource` to scope work
- **verifiable** — a `test` subcommand runs real `samtools faidx` / `tabix` queries

---

## Install

```bash
pip install refbox
```

`refbox` shells out to the standard htslib tooling. Install once via
conda/mamba (recommended):

```bash
mamba install -c bioconda htslib samtools
```

Required CLI tools: `bgzip`, `tabix`, `samtools`, GNU `sort`, `grep`.

---

## Quick start

```bash
# Download → build → validate Human GRCh38 (uses bundled species.yaml)
refbox download --species Homo_sapiens --assembly GRCh38
refbox build    --species Homo_sapiens --assembly GRCh38
refbox test     --species Homo_sapiens --assembly GRCh38
```

Or use the one-shot driver:

```bash
git clone https://github.com/typekey/refbox.git
cd refbox
./build.sh Homo_sapiens GRCh38
```

---

## CLI reference

```
refbox download [--species ...] [--assembly ...] [--resource ...] [--out DIR] [--force]
refbox build    [--species ...] [--assembly ...] [--resource ...] [--out DIR] [--force]
refbox test     [--species ...] [--assembly ...] [--out DIR]
```

| Flag | Default | Meaning |
|---|---|---|
| `--species` | all enabled | filter to one or more species (e.g. `Homo_sapiens`) |
| `--assembly` | all enabled | filter to one or more assemblies (e.g. `GRCh38`) |
| `--resource` | all 10 | subset: `genome transcriptome annotation_gtf annotation_gff3 repeats_rmsk repeats_bed repeats_gtf repeats_fa rnacentral ccre` |
| `--out` | `$REFBOX_OUT` or `$PWD` | output root |
| `--force` | off | rebuild even when outputs exist |
| `-v / --verbose` | | enable DEBUG logging |

Environment variables:

| Name | Meaning |
|---|---|
| `REFBOX_OUT` | default output root for `{Species}/{Assembly}/{raw,build}/` |
| `REFBOX_CONFIG` | path to a custom `species.yaml` (overrides the bundled registry) |

---

## Output layout

```
{REFBOX_OUT}/
  {Species}/
    {Assembly}/
      raw/                          # original downloads / copies
        genome.fa
        transcriptome.fa
        annotation_gtf.gtf
        annotation_gff3.gff3
        repeats_rmsk.tsv
        repeats_fa.fa
        rnacentral.gff3
        ccre.bed
      build/                        # browser-loadable
        genome.fa.gz                + .fai + .gzi
        chrom.sizes
        transcripts.fa.gz           + .fai
        annotation.sorted.gtf.gz    + .tbi
        annotation.sorted.gff3.gz   + .tbi
        repeats.sorted.bed.gz       + .tbi
        repeats.sorted.gtf.gz       + .tbi
        rnacentral.sorted.gff3.gz   + .tbi
        ccre.sorted.bed.gz          + .tbi
```

---

## Examples

### 1. Build only cCREs for Human GRCh38

```bash
refbox download --species Homo_sapiens --assembly GRCh38 --resource ccre
refbox build    --species Homo_sapiens --assembly GRCh38 --resource ccre
refbox test     --species Homo_sapiens --assembly GRCh38
```

### 2. Build a full reference into a specific directory

```bash
refbox download --species Mus_musculus --assembly GRCm38 --out /data/refs
refbox build    --species Mus_musculus --assembly GRCm38 --out /data/refs
refbox test     --species Mus_musculus --assembly GRCm38 --out /data/refs
```

### 3. Use a private YAML registry

```bash
export REFBOX_CONFIG=/path/to/my_species.yaml
refbox build
```

### 4. Drive everything from `build.sh`

```bash
./build.sh                              # all enabled assemblies
./build.sh Homo_sapiens                 # one species, all enabled assemblies
./build.sh Homo_sapiens GRCh38          # one species + assembly
./build.sh Homo_sapiens GRCh38 -- --resource genome ccre
FORCE=1 ./build.sh Mus_musculus GRCm38  # rebuild even when outputs exist
STEPS="test" ./build.sh                 # only run validation
```

---

## Tutorial — adding a new assembly

The registry lives in [`config/species.yaml`](config/species.yaml). Each entry
is **three levels deep**: species → assembly → resource.

### Step 1. Add an assembly block

```yaml
species:
  Homo_sapiens:
    GRCh38:
      enabled: true                  # set false to keep it idle
      gencode_version: 44
      ucsc_db: hg38

      genome:
        url: https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/GRCh38.primary_assembly.genome.fa.gz

      annotation_gtf:
        url: https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.annotation.gtf.gz

      repeats_rmsk:
        url: https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/rmsk.txt.gz

      ccre:
        url: https://downloads.wenglab.org/V3/GRCh38-cCREs.bed
```

### Step 2. Each resource entry follows a fallback rule

```yaml
my_resource:
  local_path: /path/on/disk/file.fa.gz   # used if it exists
  url:        https://.../file.fa.gz     # else downloaded
# omit / null = skipped silently
```

- `local_path` exists → **copy** into `raw/` (auto-gunzips `.gz`)
- otherwise `url` → **download** into `raw/`
- otherwise → **skip** (no error)

### Step 3. Canonical resource names

| Name | Output (raw) | Output (build) |
|---|---|---|
| `genome` | `genome.fa` | `genome.fa.gz` + `.fai` + `.gzi`, `chrom.sizes` |
| `transcriptome` | `transcriptome.fa` | `transcripts.fa.gz` + `.fai` |
| `annotation_gtf` | `annotation_gtf.gtf` | `annotation.sorted.gtf.gz` + `.tbi` |
| `annotation_gff3` | `annotation_gff3.gff3` | `annotation.sorted.gff3.gz` + `.tbi` |
| `repeats_rmsk` | `repeats_rmsk.tsv` | _(raw input for repeats_bed/gtf — derivation TODO)_ |
| `repeats_bed` | `repeats_bed.bed` | `repeats.sorted.bed.gz` + `.tbi` |
| `repeats_gtf` | `repeats_gtf.gtf` | `repeats.sorted.gtf.gz` + `.tbi` |
| `repeats_fa` | `repeats_fa.fa` | _(RepeatMasker `.fa.out` report)_ |
| `rnacentral` | `rnacentral.gff3` | `rnacentral.sorted.gff3.gz` + `.tbi` |
| `ccre` | `ccre.bed` | `ccre.sorted.bed.gz` + `.tbi` |

### Step 4. Run

```bash
refbox download --species Homo_sapiens --assembly GRCh38
refbox build    --species Homo_sapiens --assembly GRCh38
refbox test     --species Homo_sapiens --assembly GRCh38
```

---

## Programmatic API

```python
from refbox.config import load_config, iter_targets
from refbox.download import download_targets
from refbox.build import build_targets
from refbox.test import test_targets

cfg = load_config()
for t in iter_targets(cfg, species=["Homo_sapiens"]):
    print(t.species, t.assembly, list(t.resources))

download_targets(species=["Homo_sapiens"], assembly=["GRCh38"], out="/data/refs")
build_targets(   species=["Homo_sapiens"], assembly=["GRCh38"], out="/data/refs")
test_targets(    species=["Homo_sapiens"], assembly=["GRCh38"], out="/data/refs")
```

---

## Development

```bash
git clone https://github.com/typekey/refbox.git
cd refbox
pip install -e .
refbox --help
```

### Release

Tags matching `v*` automatically build and publish to PyPI via GitHub Actions
([`.github/workflows/workflow.yml`](.github/workflows/workflow.yml)) using
[PyPI trusted publishing](https://docs.pypi.org/trusted-publishers/) (no API
token required).

```bash
git tag v0.1.0
git push origin v0.1.0
```

Manual build (no upload):

```bash
./release.sh build       # writes dist/
```

---

## License

MIT — see [LICENSE](LICENSE).
