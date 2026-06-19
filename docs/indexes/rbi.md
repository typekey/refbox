# Lite Index (`.rbi`)

The full [`.rba` index](rba.md) (SQLite + FTS5 + trigram + annotation structure)
is ~1 GB for GENCODE. When you only need **gene/transcript name → position**
autocomplete — not substring/fuzzy alias search or the full transcript structure —
the **lightweight RBrowser Index** (`.rbi`) is a fraction of the size (a few tens
of MB) and resolves queries with plain B-tree seeks.

| | |
|---|---|
| **Name** | RBrowser Index |
| **Short name** | RBI |
| **File extension** | `.rbi` |
| **Storage engine** | SQLite (no FTS5, no trigram by default) |
| **Search** | B-tree range scans on a normalized `term` table (exact / prefix / autocomplete); optional 3-gram table for fuzzy recall |
| **Holds** | display name + genomic position per gene/transcript (no exon/CDS structure, no payload JSON) |

## How to build it

```bash
# Standalone lightweight index
refbox build -rbi gencode.v45.annotation.gtf.gz -o hg38.gencode.v45.rbi \
    --source-name GENCODE --species human --genome hg38 \
    --annotation-version v45 --force
```

Or emit it **alongside** a tabix annotation in one command:

```bash
refbox build -gtf gencode.v45.annotation.gtf.gz --with-rbi \
    --source-name GENCODE --genome hg38 --annotation-version v45
```

### Options

| Flag | Meaning |
|---|---|
| `-rbi FILE` | GTF/GFF3 input → build a standalone `.rbi` (default output `<stem>.rbi`) |
| `--with-rbi` | with `-gtf`/`-gff`: also emit a `.rbi` next to the sorted/bgzip/tabix output |
| `-o PATH` | output path |
| `--source-name` / `--species-name` / `--genome` / `--annotation-version` | metadata labels stored in the index |
| `--no-gram3` | skip the 3-gram table (smallest index; exact + prefix lookup only) |
| `--force` | overwrite an existing output |

!!! note "3-gram fuzzy recall is on by default"
    By default the `.rbi` includes a 3-gram table over gene/transcript **names**,
    enabling substring/fuzzy recall (e.g. `alat` → `MALAT1`). Pass `--no-gram3`
    for the smallest possible file when you only need exact + prefix lookup.

## Schema

| Table | Purpose |
|---|---|
| `record` | one row per gene/transcript — display name + genomic position only |
| `term` | `(term_norm, term_raw, field, record_id, priority)` `WITHOUT ROWID`; PK leads with `term_norm` so exact (`= ?`) and prefix (`>= lo AND < hi`) are B-tree range scans |
| `gram3` | optional 3-gram → record map for fuzzy/partial recall (default on; `--no-gram3` to skip) |
| `metadata` | format (`RBI` / "RBrowser Index") + provenance (source / species / genome / version / counts) |

### Term normalization

For each raw value the builder stores up to three lookup keys in `term`:

1. **normalized** — trimmed + lowercased (separators preserved): `TP53` → `tp53`
2. **version-stripped** — drop a trailing Ensembl `.<version>`: `ENST…305.10` → `enst…305`
3. **separator-free** — drop `_ - . space`: `TP53-201` → `tp53201`

Field priority drives ranking (lower = higher): `gene_name` (10) → `gene_id` (20)
→ `transcript_name` (30) → `transcript_id` (40) → `transcript_biotype` (80).
Only names get 3-grams — fuzzy-matching an accession or biotype is pointless and
bloats the table.

### Search tiers

`search()` stops at the first tier that yields rows:

```text
1. exact term_norm
2. prefix B-tree range
3. separator-free exact
4. version-stripped exact
5. gram3 fuzzy recall (if the table is present)
```

No full-table `LIKE '%q%'` ever runs.

## Programmatic use

```python
from refbox.lite_index import build_lite_index, open_readonly, search

db = build_lite_index("gencode.v45.annotation.gtf.gz", "hg38.gencode.v45.rbi",
                      source_name="GENCODE", genome="hg38", force=True)
con = open_readonly(db)
mode, results = search(con, "TP53", limit=10)   # mode: exact / prefix / normalized / gram3 / none
for hit in results:
    print(hit["gene_name"], hit["chrom"], hit["start"], hit["end"])
```

## `.rba` vs `.rbi` — which to build?

| Need | Use |
|---|---|
| substring/fuzzy alias search, synonyms, RNAcentral ncRNAs, full transcript structure (exons/CDS/UTR) to render without extra requests | **`.rba`** |
| fast name→position autocomplete, smallest static file | **`.rbi`** |
