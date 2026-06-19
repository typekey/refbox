# Quick Start

## 1. Build a configured assembly end-to-end

The registry-driven pipeline (`pull`) does **download → build → test → publish**
in one shot, using the bundled `species.yaml`:

```bash
# Fetch + build + validate + publish Human GRCh38
refbox pull --species Homo_sapiens --assembly GRCh38

# --species is optional; it is inferred from --assembly
refbox pull --assembly GRCm38
```

Outputs land under `$REFBOX_OUT/<Species>/<Assembly>/` as flattened
`<Assembly>.<name>` files (e.g. `GRCh38.genome.fa.gz`). See
[Output Layout](config/output-layout.md).

## 2. Build a single arbitrary file

No `species.yaml` entry needed — bring your own data:

```bash
# Genome FASTA → bgzip + faidx + chrom.sizes
refbox build -fa my_genome.fa

# GTF/GFF3 → sorted + bgzip + tabix
refbox build -gtf my_annotation.gtf

# BED → sorted + bgzip + tabix + bigBed
refbox build -bed peaks.bed --assembly GRCh38
```

See the full [Build Command](commands/build.md) reference.

## 3. Build a search index (the question everyone asks)

=== "Full index (.rba)"

    ```bash
    # Standalone RBrowser Index (SQLite + FTS5 search)
    refbox build -rba gencode.v45.annotation.gtf.gz -o hg38.gencode.v45.rba \
        --source-name GENCODE --species human --genome hg38 \
        --annotation-version v45 --force
    ```

    Or emit it **alongside** the tabix output in one command:

    ```bash
    refbox build -gtf gencode.v45.annotation.gtf.gz --with-rba \
        --source-name GENCODE --genome hg38 --annotation-version v45
    ```

    Details → [RBrowser Index `.rba`](indexes/rba.md).

=== "Lite index (.rbi)"

    ```bash
    refbox build -rbi gencode.v45.annotation.gtf.gz -o hg38.gencode.v45.rbi \
        --source-name GENCODE --species human --genome hg38 \
        --annotation-version v45 --force
    ```

    Details → [Lite Index `.rbi`](indexes/rbi.md).

## 4. The one-shot driver script

`build.sh` wraps `refbox pull` with sensible defaults (output root = repo parent):

```bash
./build.sh                            # all enabled assemblies
./build.sh Homo_sapiens               # one species (all enabled assemblies)
./build.sh Homo_sapiens GRCh38        # one species + assembly
./build.sh -- --resource genome cytoband   # extra args after `--` go to refbox

# environment toggles
FORCE=1 ./build.sh Homo_sapiens GRCh38     # re-download + rebuild
NO_TEST=1 ./build.sh                       # skip validation
INCLUDE_DISABLED=1 ./build.sh              # also process enabled: false entries
```
