# Python API

Everything the CLI does is also available programmatically.

## Imports

```python
from refbox.config      import load_config, iter_targets, RESOURCE_NAMES
from refbox.download    import download_targets
from refbox.build       import build_targets, publish_targets
from refbox.test        import test_targets
from refbox             import file_build as fb        # single-file builders
from refbox.ingest      import ingest_directory        # directory ingest
from refbox.report      import build_report            # Markdown status report
from refbox.sqlite_index import (                      # full RBrowser Index (.rba)
    build_sqlite_index, search, inspect, open_readonly, normalize)
from refbox.lite_index  import (                       # lite index (.rbi)
    build_lite_index, search as lite_search, open_readonly as lite_open, inspect as lite_inspect)
```

## Registry pipeline

```python
# Equivalent to: refbox pull --assembly GRCh38
build_targets(species=["Homo_sapiens"], assembly=["GRCh38"], auto_download=True)
test_targets(species=["Homo_sapiens"], assembly=["GRCh38"])
publish_targets(species=["Homo_sapiens"], assembly=["GRCh38"])
```

## Single-file builders

```python
from pathlib import Path
from refbox import file_build as fb

fb.build_fa(Path("genome.fa"))                                   # → genome.fa.gz + .fai + .gzi
fb.build_gxf(Path("annot.gtf"))                                  # → annot.sorted.gtf.gz + .tbi
fb.build_bed(Path("peaks.bed"), assembly="GRCh38")              # → sorted + tabix + bigBed
fb.build_rmsk(Path("rmsk.txt.gz"))                              # → repeats.sorted.{bed,gtf}.gz
fb.build_transcriptome(Path("genome.fa"), Path("annot.gtf"))   # → transcriptome.fa.gz
```

## Build & query a full index (`.rba`)

```python
from refbox.sqlite_index import build_sqlite_index, open_readonly, search

db = build_sqlite_index("gencode.v45.annotation.gtf.gz", "idx.rba",
                        source_name="GENCODE", genome="hg38", force=True)
con = open_readonly(db)
for hit in search(con, "TP53", limit=10):
    print(hit["matched_field"], hit["gene_name"], hit["transcript_id"])
```

`build_sqlite_index` signature:

```python
build_sqlite_index(
    input_path, output=None, *,
    source_name="", species="", genome="", annotation_version="",
    synonyms=None,          # HGNC-style TSV path
    rnacentral=None,        # RNAcentral genome-coordinates GFF3 path
    fuzzy_scope="names",    # or "all"
    force=False, verbose=False,
) -> Path
```

## Build & query a lite index (`.rbi`)

```python
from refbox.lite_index import build_lite_index, open_readonly, search

db = build_lite_index("gencode.v45.annotation.gtf.gz", "idx.rbi",
                      source_name="GENCODE", genome="hg38",
                      enable_gram3=True, force=True)
con = open_readonly(db)
mode, results = search(con, "TP53", limit=10)   # NOTE: returns (mode, results)
for hit in results:
    print(hit["gene_name"], hit["chrom"], hit["start"], hit["end"])
```

`build_lite_index` signature:

```python
build_lite_index(
    input_path, output=None, *,
    source_name="", species="", genome="", annotation_version="",
    enable_gram3=False, force=False, verbose=False,
) -> Path
```

!!! note "Reference search implementation"
    The browser-side search should mirror the tiered ranking in
    `refbox.sqlite_index.search()` (for `.rba`, returns a list of hits) and
    `refbox.lite_index.search()` (for `.rbi`, returns a `(mode, results)` tuple) —
    these Python functions are the reference implementations.
