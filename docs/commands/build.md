# Build Command (`refbox build`)

`refbox build` builds **arbitrary user-supplied files** that have no
`species.yaml` entry — custom organisms, in-house assemblies, pre-processed BEDs
from a paper. Each call detects the input by extension (or by explicit flag),
verifies bgzip / sort order, runs the canonical transformation, and emits indexed
outputs.

Inputs may be plain or `.gz`. The builder transparently re-bgzips a plain gzip,
copies an already-bgzipped file as-is, and re-sorts an unsorted GFF/BED before
indexing.

## Synopsis

```text
refbox build [INPUT] [-o OUT] [INPUT-TYPE FLAG] [OPTIONS]
```

## Input selectors

Pick the input with one of these flags, or pass a single positional `INPUT` file
to auto-detect its type by extension.

| Flag | Input | Output |
|---|---|---|
| `-fa`, `--fa FILE` | genome FASTA | `.fa.gz` + `.fai` + `.gzi` + `chrom.sizes` |
| `-gtf`, `--gtf FILE` | GTF annotation | sorted `.gtf.gz` + `.tbi` |
| `-gff`, `--gff FILE` | GFF3 annotation | sorted `.gff3.gz` + `.tbi` |
| `-bed`, `--bed FILE` | BED features | sorted `.bed.gz` + `.tbi` (+ bigBed) |
| `-rmsk`, `--rmsk FILE` | UCSC `rmsk.txt[.gz]` | `repeats.sorted.bed.gz` + `repeats.sorted.gtf.gz` (+ `.tbi`) |
| `-rba`, `--rba FILE` | GTF/GFF3 → standalone **full index** (`.rba`) | `.rba` (SQLite + FTS5 + annotation) |
| `-rbi`, `--rbi FILE` | GTF/GFF3 → standalone **lite index** (`.rbi`) | `.rbi` (compact B-tree lookup) |
| `-i`, `--ingest DIR` | directory of user files | full `{Species}/{Assembly}/` layout |
| `INPUT` (positional) | single file, auto-detected | as per detected type |
| `-fa … -gtf …` together | genome + annotation | transcriptome FASTA (via gffread) |

!!! note "`.rba` and `.rbi` are different builders"
    `-rba` builds the full FTS5 + annotation index; `-rbi` builds the compact
    lite index. They are **not** aliases. `-sqlite` / `--sqlite` and
    `--with-sqlite` remain accepted aliases for `-rba` / `--with-rba`.
    See [Full Index `.rba`](../indexes/rba.md) and [Lite Index `.rbi`](../indexes/rbi.md).

## Options

| Flag | Applies to | Meaning |
|---|---|---|
| `-o`, `--out PATH` | all | output file or directory |
| `--with-rba` | `-gtf` / `-gff` | also emit a full `.rba` index alongside the tabix output |
| `--with-rbi` | `-gtf` / `-gff` | also emit a lite `.rbi` index alongside the tabix output |
| `--source-name NAME` | `.rba` / `.rbi` | annotation source label in metadata (e.g. `GENCODE`, `Ensembl`) |
| `--species-name NAME` | `.rba` / `.rbi` | species label stored in metadata |
| `--genome NAME` | `.rba` / `.rbi` | genome/assembly label in metadata (e.g. `hg38`) |
| `--annotation-version V` | `.rba` / `.rbi` | annotation version in metadata (e.g. `v45`) |
| `--no-gram3` | `.rbi` builds | skip the 3-gram table (smallest index; exact + prefix only) |
| `--synonyms`, `--synonyms-file FILE` | `.rba` builds | HGNC-style TSV injected as `gene_synonym` aliases |
| `--rnacentral FILE` | `.rba` builds | RNAcentral genome-coordinates GFF3 merged in as ncRNAs |
| `--fuzzy-scope {names,all}` | `.rba` builds | trigram corpus: `names` (default) or `all` (include IDs) |
| `--chrom-sizes FILE` | `-bed` | chrom.sizes for bigBed conversion |
| `--assembly NAME` | `-bed`, `-i` | assembly id (chrom.sizes lookup / folder name) |
| `--no-bigbed` | `-bed` | skip bigBed generation |
| `--species NAME` | `-i` | species name (ingest only) |
| `--map RES:PATH ...` | `-i` | resource → path overrides (ingest only) |
| `--no-build` | `-i` | copy raws but skip the build pipeline (ingest only) |
| `--force` | all | overwrite existing outputs |

---

## Examples by input type

### Genome FASTA

```bash
# → bgzip + faidx + chrom.sizes  (output defaults to <input>.fa.gz)
refbox build -fa GENOME.fa
refbox build -fa GENOME.fa -o my_genome.fa.gz
```

### GTF / GFF3 annotation

```bash
refbox build -gtf ANNOT.gtf                 # → ANNOT.sorted.gtf.gz + .tbi
refbox build -gff ANNOT.gff3 -o out.gff3.gz # → sorted + bgzip + tabix

# also emit a search index alongside the tabix output
refbox build -gtf ANNOT.gtf --with-rba \
             --source-name GENCODE --genome hg38 --annotation-version v45
```

### Transcriptome (genome + annotation)

```bash
# extract transcript sequences via gffread, then bgzip + faidx
refbox build -fa GENOME.fa -gtf ANNOT.gtf -o transcriptome.fa.gz
```

### BED features

```bash
# chrom.sizes resolved from an explicit file…
refbox build -bed FEATURES.bed --chrom-sizes hg38.chrom.sizes

# …or from a known assembly (delegates to zlbio's per-species lookup)
refbox build -bed FEATURES.bed --assembly GRCh38

# skip bigBed (tabix only)
refbox build -bed FEATURES.bed --no-bigbed
```

### UCSC RepeatMasker table

```bash
# rmsk.txt[.gz] → repeats.sorted.bed.gz + repeats.sorted.gtf.gz (+ .tbi)
refbox build -rmsk rmsk.txt.gz -o OUT_DIR
```

### Standalone search index

```bash
# GTF/GFF3 → full RBrowser Annotation index (.rba: FTS5 + annotation structure)
refbox build -rba ANNOT.gtf.gz -o OUT.rba \
             --source-name GENCODE --species human --genome hg38 \
             --annotation-version v45 --force

# GTF/GFF3 → lightweight RBrowser Index (.rbi: compact name->position lookup)
refbox build -rbi ANNOT.gtf.gz -o OUT.rbi \
             --source-name GENCODE --species human --genome hg38 \
             --annotation-version v45 --force
```

See [Full Index `.rba`](../indexes/rba.md) and [Lite Index `.rbi`](../indexes/rbi.md).

### Directory ingest

```bash
# auto-classify a directory of files by extension into the full layout
refbox build -i DIR --assembly GRCh38 --species Homo_sapiens

# override which file maps to which resource
refbox build -i DIR --assembly GRCh38 \
             --map genome:genome.fa --map annotation_gtf:genes.gtf

# copy raws but skip the build pipeline
refbox build -i DIR --assembly GRCh38 --no-build
```

### Auto-detect

```bash
# infer type from the extension (rmsk/.out/.fa/.gtf/.gff3/.bed)
refbox build SOMEFILE
```
