# Full Index (`.rba`)

`tabix` answers *positional* range queries; it cannot answer "**what is TP53?**".
For that, refbox builds a standalone, read-only **full RBrowser Annotation index**
— a `.rba` file — from a GTF/GFF3. Besides the search index it also carries the
**annotation content** (per gene/transcript exon/CDS/UTR structure), powering
exact / prefix-autocomplete / fuzzy-substring / alias search **and** in-browser
rendering without extra requests.

| | |
|---|---|
| **Name** | RBrowser Indexed Annotation Database |
| **File extension** | `.rba` |
| **Storage engine** | SQLite |
| **Search engine** | SQLite FTS5 (prefix + trigram) |
| **Access** | static file + SQLite WASM / HTTP Range VFS (no backend) |
| **Purpose** | gene / transcript / alias search **plus** annotation structure for rendering |

For a much smaller index that only does name→position lookup (no FTS5, no
annotation structure), see the [Lite Index `.rbi`](rbi.md).

## How to build it

=== "Standalone"

    ```bash
    refbox build -rba gencode.v45.annotation.gtf.gz -o hg38.gencode.v45.transcript.rba \
        --source-name GENCODE --species human --genome hg38 \
        --annotation-version v45 --force
    ```

=== "Alongside tabix"

    ```bash
    refbox build -gtf gencode.v45.annotation.gtf.gz --with-rba \
        --source-name GENCODE --species human --genome hg38 \
        --annotation-version v45
    ```

!!! tip "Aliases & default name"
    `-sqlite` / `--with-sqlite` are accepted aliases for `-rba` / `--with-rba`.
    The default output name is `<input-stem>.rba`.

    Note: `-rbi` is **not** an alias for `-rba` — it builds the separate
    [Lite Index `.rbi`](rbi.md).

## Metadata flags

All optional; they are stored verbatim in the index `metadata` table so the
browser can display provenance.

| Flag | Example | Meaning |
|---|---|---|
| `--source-name` | `GENCODE` | annotation source label |
| `--species-name` | `human` | species label |
| `--genome` | `hg38` | genome / assembly label |
| `--annotation-version` | `v45` | annotation version |

## Synonyms (HGNC)

GENCODE/Ensembl annotations do **not** carry common gene synonyms (`OCT4`, `p53`,
`Nanog`…). Pass `--synonyms` with an HGNC `hgnc_complete_set.txt` to inject
`alias_symbol` / `prev_symbol` as searchable `gene_synonym` aliases (matched by
Ensembl ID, fallback by symbol). Then `OCT4` → POU5F1, `OTF3` → POU5F1.

```bash
refbox build -rba gencode.v45.annotation.gtf.gz -o OUT.rba \
             --synonyms hgnc_complete_set.txt --force
```

## RNAcentral ncRNAs

`--rnacentral` merges an RNAcentral genome-coordinates GFF3 (use the
chromosome-**normalized** file so coordinates display consistently) as additional
ncRNA records. They become searchable by full ID (`URS…_9606.N`), URS accession
(`URS…_9606`) and versionless URS (`URS000035F234`).

```bash
# GENCODE genes/transcripts + RNAcentral ncRNAs + HGNC synonyms → one index
refbox build -rba gencode.v45.annotation.gff3.gz -o human.rba \
             --rnacentral rnacentral.normalized.gff3 \
             --synonyms hgnc_complete_set.txt --force
```

RNAcentral's only label is a long free-text `description`, so the display name is
distilled to a short, recognizable symbol — `DDX11L11`, `miR-34a-5p`,
`pre-mir-571`, `piR-hsa-4818588`, `5S-rRNA`, `Metazoan-SRP-RNA`. The full original
description is kept as a searchable alias, so recall is unchanged.

## Fuzzy scope

`--fuzzy-scope` controls what the trigram (substring) index covers:

| Value | Behavior |
|---|---|
| `names` *(default)* | trigram over gene/transcript **names + synonyms only**; IDs are exact/prefix-only. Trigram is ~6× smaller; the pathological numeric-ID-substring scan disappears. |
| `all` | include IDs in the trigram corpus (legacy behavior). |

## What it indexes

One `feature` row per **gene** and per **transcript**, with full structure (exon
starts/ends, CDS span, 5′/3′ UTR spans, biotype, source), a `search_text` blob,
and a `payload_json` for the browser to render without extra joins. Every
searchable synonym becomes an `alias` row: names & IDs, **versionless** Ensembl
IDs (`ENST00000269305.9` → `ENST00000269305`), HAVANA / CCDS / HGNC / protein IDs,
GFF3 `Alias` / `Dbxref` (RefSeq) / `gene_synonym`.

| Table | Purpose |
|---|---|
| `feature` | gene + transcript records (1-based inclusive coords; `chrom_start0/end0` give the 0-based half-open span) |
| `alias` | `(feature_id, alias, alias_norm, alias_type, source)` — `alias_norm` is lowercased, version-stripped, separator-free |
| `metadata` | `key/value` provenance (source, species, genome, coord convention, counts, capabilities) |
| `feature_fts` | FTS5 prefix/autocomplete (`prefix='2 3 … 10'`) |
| `feature_trigram` | FTS5 `trigram` tokenizer for substring/fuzzy search (graceful `LIKE` fallback if unavailable) |

## Search ranking

```text
exact transcript_id → transcript_name → gene_name → gene_id → alias exact
→ FTS prefix → trigram substring → LIKE fallback
```

The browser should run the same tiers; `refbox.sqlite_index.search()` is the
reference implementation. Exact-match columns are indexed `COLLATE NOCASE` so
lookups are index seeks, not table scans.

## Coordinate convention

`start` / `end` and all `*_start` / `*_end` columns are **1-based inclusive**
(GTF/GFF convention) for display; `chrom_start0` / `chrom_end0` add the **0-based
half-open** span for rendering. This is recorded in `metadata` under
`coord_convention`.

## Inspect & query

```python
from refbox.sqlite_index import build_sqlite_index, open_readonly, search, inspect

print(inspect("idx.rba")["metadata"])          # provenance + capabilities
con = open_readonly("idx.rba")
for hit in search(con, "TP53", limit=10):
    print(hit["matched_field"], hit["gene_name"], hit["transcript_id"])
```
