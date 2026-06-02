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

# …and also emit a static SQLite search index alongside the tabix output:
refbox build -gtf ANNOT.gtf --with-sqlite \
             --source-name GENCODE --genome hg38 --annotation-version v45

# Standalone: GTF/GFF3 → read-only SQLite search index (no tabix/bgzip)
refbox build -sqlite ANNOT.gtf.gz -o OUT.rbrowser.sqlite \
             --source-name GENCODE --species human --genome hg38 \
             --annotation-version v45 --force

# Genome + annotation → transcriptome FASTA (via gffread) + faidx
refbox build -fa GENOME.fa -gtf ANNOT.gtf -o transcriptome.fa.gz

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

## SQLite search index (for in-browser transcript/gene search)

`tabix` answers *positional* range queries; it cannot answer "what is TP53?".
For that, `refbox` can build a **standalone, read-only SQLite file** from a
GTF/GFF3 that powers exact / prefix-autocomplete / fuzzy-substring / alias
search. It is meant to be hosted as a static file and queried directly from the
browser via **SQLite WASM + an HTTP Range VFS** — no backend service required.

```bash
# build it next to the tabix output (one command)
refbox build -gtf gencode.v45.annotation.gtf.gz --with-sqlite \
             --source-name GENCODE --species human --genome hg38 \
             --annotation-version v45

# or standalone, anywhere
refbox build -sqlite gencode.v45.annotation.gtf.gz -o hg38.gencode.v45.rbrowser.sqlite \
             --source-name GENCODE --genome hg38 --annotation-version v45 --force

# with HGNC synonyms so common names resolve (OCT4 -> POU5F1, OTF3 -> POU5F1):
refbox build -sqlite gencode.v45.annotation.gtf.gz -o OUT.rbrowser.sqlite \
             --synonyms hgnc_complete_set.txt --force
```

> GENCODE/Ensembl annotations do **not** carry common gene synonyms (`OCT4`,
> `p53`, `Nanog`…). Pass `--synonyms` with an HGNC `hgnc_complete_set.txt` to
> inject `alias_symbol` / `prev_symbol` as searchable `gene_synonym` aliases.

```bash
# merge GENCODE genes/transcripts + RNAcentral ncRNAs + HGNC synonyms into one index
refbox build -sqlite gencode.v45.annotation.gff3.gz -o human.rbrowser.sqlite \
             --rnacentral rnacentral.normalized.gff3 \
             --synonyms hgnc_complete_set.txt --force
```

> `--rnacentral` merges an RNAcentral genome-coordinates GFF3 (use the
> chromosome-**normalized** file so coordinates display consistently). ncRNAs
> become searchable by full ID (`URS…_9606.N`), URS accession (`URS…_9606`) and
> versionless URS (`URS000035F234`).

The three helper scripts in [`script/`](script/) are dependency-free (Python
standard library only — they load the self-contained `refbox.sqlite_index`
module directly, so they run even without the rest of refbox installed):

```bash
python script/build_rbrowser_sqlite_index.py  --input ANNOT.gtf.gz --output idx.sqlite \
        --source-name GENCODE --species human --genome hg38 --annotation-version v45 --force
python script/inspect_rbrowser_sqlite_index.py --db idx.sqlite
python script/test_rbrowser_sqlite_search.py   --db idx.sqlite \
        --queries TP53 ENST00000269305 p53 BRCA1 MALAT1 ACTB --repeat 100 --limit 10
```

### What it indexes

One `feature` row per **gene** and per **transcript**, with the full structure
(exon starts/ends, CDS span, 5′/3′ UTR spans, biotype, source), a `search_text`
blob, and a `payload_json` for the browser to render without extra joins. Every
searchable synonym becomes an `alias` row: gene/transcript names & IDs, the
**versionless** Ensembl IDs (`ENST00000269305.9` → `ENST00000269305`), HAVANA /
CCDS / HGNC / protein IDs, GFF3 `Alias` / `Dbxref` (RefSeq) / `gene_synonym`.

| Table | Purpose |
|---|---|
| `feature` | gene + transcript records (1-based inclusive coords; `chrom_start0/end0` give the 0-based half-open span) |
| `alias` | `(feature_id, alias, alias_norm, alias_type, source)` — `alias_norm` is lowercased, version-stripped, separator-free |
| `metadata` | `key/value` provenance (source, species, genome, coord convention, counts, SQLite/FTS capabilities) |
| `feature_fts` | FTS5 prefix/autocomplete (`prefix='2 3 … 10'`) |
| `feature_trigram` | FTS5 `trigram` tokenizer for substring/fuzzy search (graceful LIKE fallback if unavailable) |

### Search ranking (mirrored by `test_rbrowser_sqlite_search.py`)

exact `transcript_id` → `transcript_name` → `gene_name` → `gene_id` → alias
exact → FTS prefix → trigram substring → LIKE fallback. The browser should run
the same tiers (the Python `refbox.sqlite_index.search()` is the reference
implementation). Exact-match columns are indexed `COLLATE NOCASE` so lookups are
index seeks, not table scans.

### Coordinate convention

`start` / `end` and all `*_start` / `*_end` columns are **1-based inclusive**
(GTF/GFF convention) for display; `chrom_start0` / `chrom_end0` add the
**0-based half-open** span for rendering. This is recorded in `metadata`.

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
        transcriptome.fa.gz          + .fai
        transcriptome.derived.fa.gz  + .fai   # gffread-extracted, when GTF available
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
from refbox.sqlite_index import (                # static SQLite search index
    build_sqlite_index, search, inspect, open_readonly, normalize)
```

```python
# Build a search index and query it the way the browser will.
db = build_sqlite_index("gencode.v45.annotation.gtf.gz", "idx.sqlite",
                        source_name="GENCODE", genome="hg38", force=True)
con = open_readonly(db)
for hit in search(con, "TP53", limit=10):
    print(hit["matched_field"], hit["gene_name"], hit["transcript_id"])
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
# Concatenate multiple upstream files into one raw file (e.g. Ensembl cdna+ncrna):
transcriptome:
  url:        https://ftp.ensembl.org/pub/release-111/fasta/danio_rerio/cdna/Danio_rerio.GRCz11.cdna.all.fa.gz
  extra_urls:
    - https://ftp.ensembl.org/pub/release-111/fasta/danio_rerio/ncrna/Danio_rerio.GRCz11.ncrna.fa.gz
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
# transcriptome.fa.gz from genome + GTF/GFF via gffread.
transcriptome: null
```

### Canonical resource names

| Name | raw/ | build/ |
|---|---|---|
| `genome` | `genome.fa` | `genome.fa.gz` + `.fai` + `.gzi`, `chrom.sizes` |
| `transcriptome` | `transcriptome.fa` | `transcriptome.fa.gz` + `.fai`<br>`transcriptome.derived.fa.gz` + `.fai` |
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

- **v0.5.6** — RNAcentral records are now searchable by their **description**:
  the merge reads the GFF3 `description` attribute (e.g. `(human) tRNA-Ala`),
  strips the leading `(species)` tag, and stores it as the ncRNA's `gene_name`
  (so it displays) and as a name-like alias (so it resolves by exact / prefix /
  substring — `tRNA-Ala`, `tRNA`, `piR-hsa-…`). Previously only the URS
  accession was searchable.
- **v0.5.5** — `--fuzzy-scope names` (new default): the trigram (substring) index
  is built over gene/transcript **names + synonyms only**, not IDs. IDs
  (gene_id / transcript_id / RNAcentral URS) are searched exactly or by prefix —
  never by fuzzy substring — so the trigram shrinks ~6× (measured 371 → 60 MB on
  GENCODE v44 + RNAcentral) and the pathological numeric-ID-substring query (e.g.
  `000003351` scanning every ID) disappears. Name/synonym substring search
  (`p53`, `malat`) is unchanged. Pass `--fuzzy-scope all` for the old
  IDs-included behavior.
- **v0.5.4** — SQLite search: unified exact lookup + leaner `alias` table.
  Exact matching (transcript_id / transcript_name / gene_name / gene_id /
  synonym / RNAcentral ID / …) is now a **single index-only seek** into
  `idx_alias_norm(alias_norm, alias_type, feature_id)` ranked by `alias_type`
  — species-agnostic (ENSMUSG / ENSDARG / FBgn all just work) and ~9 page reads
  instead of ~22 over an HTTP Range VFS. The `alias` table drops the original
  `alias` text and `source` columns (the readable values live on the feature row
  / `payload_json`) and no longer stores redundant `gene_id_versionless` /
  `transcript_id_versionless` rows (normalization already strips the version;
  `rnacentral_id_versionless` is kept — it strips the taxon suffix, so it is
  *not* redundant). Dropped the unused `idx_alias_feature`. Net on GENCODE v44 +
  RNAcentral + HGNC: 1638 → 1385 MB (−15%), aliases 6.5M → 5.4M. Behavior note:
  a gene-name/gene-id query now returns the **gene** record (use
  `WHERE gene_id = ? COLLATE NOCASE` to list its transcripts).
- **v0.5.3** — Optional `--rnacentral` feed for the SQLite index: an RNAcentral
  genome-coordinates GFF3 (use the chromosome-normalized one) is merged in as
  additional ncRNA transcript records (`transcript` + `noncoding_exon`;
  `predicted_gene` ignored), searchable by full ID, URS accession and
  versionless URS (`URS000035F234`). `refbox build -sqlite … --rnacentral
  rnacentral.normalized.gff3`.
- **v0.5.2** — Optional `--synonyms` feed for the SQLite index: an HGNC-style
  TSV (`symbol` / `alias_symbol` / `prev_symbol` / `ensembl_gene_id`) is injected
  as `gene_synonym` aliases (matched by Ensembl ID, fallback by symbol), so common
  names that GENCODE/Ensembl omit resolve as exact alias hits — e.g. `OCT4` →
  POU5F1, `OTF3` → POU5F1. `refbox build -sqlite … --synonyms hgnc_complete_set.txt`.
- **v0.5.1** — SQLite index size optimizations (no capability loss): `feature_fts`
  built with `detail=none, columnsize=0` (−~84%); `feature.search_text` no longer
  stored (lives only in the FTS index); `payload_json` slimmed to the alias list
  (other fields reconstructed from columns); covering `idx_alias_norm(alias_norm,
  feature_id)` for index-only alias lookups. `feature_trigram` stays `detail=full`
  (required for multi-trigram substring *phrase* queries). Net: GENCODE v44 index
  ~745 MB → ~430 MB.
- **v0.5.0** — Static **SQLite search index** for in-browser transcript/gene
  search. New `refbox.sqlite_index` module + `refbox build -sqlite` /
  `refbox build -gtf … --with-sqlite`; three stdlib-only helper scripts
  (`build_/inspect_/test_rbrowser_sqlite_index.py`). Emits a read-only
  `*.rbrowser.sqlite` (gene + transcript records, alias table, FTS5 prefix +
  trigram tables, metadata) suitable for SQLite-WASM + HTTP Range hosting.
  Ranked exact / prefix / fuzzy / alias search with version-insensitive Ensembl
  IDs; exact-match columns indexed `COLLATE NOCASE` for index-seek lookups.
- **v0.3.2** — Robust download backend fallback. `_download` now tries `axel
  → aria2c → wget → wget --no-check-certificate → requests → requests verify=False`
  in order, so a single broken TLS host (e.g. `ftp.ensemblgenomes.org` whose
  cert is not valid for its own hostname) no longer aborts a resource — it
  silently retries with the next backend. Also silenced axel/aria2c/wget
  progress noise.
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
