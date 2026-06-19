# refbox

**Build standardized, indexed reference files for genome browsers — in one command.**

`refbox` turns a YAML registry of species/assemblies into ready-to-load browser
inputs, and can also build static, backend-free **search indexes** for in-browser
gene/transcript lookup.

## What it produces

| Resource | Output |
|---|---|
| Genome FASTA | `bgzip` + `samtools faidx` (`.fa.gz` + `.fai` + `.gzi`) + `chrom.sizes` |
| Transcriptome FASTA | `bgzip` + `faidx` (auto-derived from genome + GTF when no upstream URL exists) |
| GTF / GFF3 annotations | sorted + `bgzip` + `tabix` |
| Repeats | UCSC RepeatMasker `rmsk.txt.gz` → BED + GTF; `.fa.out.gz` report |
| RNAcentral ncRNAs | direct download *or* liftover from another assembly |
| ENCODE SCREEN cCREs | sorted + `bgzip` + `tabix` |
| UCSC cytogenetic bands | tabix BED + bigBed ideogram |
| **Full Index** (`.rba`) | SQLite + FTS5 + trigram + annotation structure — exact / prefix / fuzzy / alias search and rendering |
| **Lite Index** (`.rbi`) | compact SQLite B-tree lookup — fast name→position autocomplete |

It ships a registry of **26 species / 42 assemblies** (human, mouse, rat, dog,
cow, pig, chimp, gorilla, zebrafish, fly, worm, sea urchin, yeast, plants,
bacteria, viruses…) covering GENCODE, Ensembl, Ensembl Genomes, UCSC golden
path, NCBI, RNAcentral, and ENCODE SCREEN.

## Command surface

```text
refbox download   # only fetch raw files configured in species.yaml
refbox pull       # full pipeline: download (if missing) + build + test + publish
refbox publish    # flatten an existing build/ tree to <Assembly>.<name>
refbox test       # validate build/ outputs
refbox build      # single-file / directory build for arbitrary inputs
```

## Where to go next

<div class="grid cards" markdown>

- :material-download: **[Installation](installation.md)** — pip + external bioinformatics tools.
- :material-rocket-launch: **[Quick Start](quickstart.md)** — the most common workflows.
- :material-console: **[Pipeline Commands](commands/pipeline.md)** — `download` / `pull` / `publish` / `test`.
- :material-file-cog: **[Build Command](commands/build.md)** — build any file or directory by hand.
- :material-database-search: **[Full Index `.rba`](indexes/rba.md)** — FTS5 search + annotation structure.
- :material-database: **[Lite Index `.rbi`](indexes/rbi.md)** — compact name→position lookup.

</div>

!!! tip "Looking for how to build a search index?"
    The short answer: `refbox build -rba ANNOT.gtf.gz -o OUT.rba` (full index) or
    `refbox build -rbi ANNOT.gtf.gz -o OUT.rbi` (lite index). See
    **[Full Index](indexes/rba.md)** and **[Lite Index](indexes/rbi.md)**.
