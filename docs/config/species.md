# `species.yaml` Registry

The registry that drives `download` / `pull` / `publish` / `test` lives in
[`config/species.yaml`](https://github.com/typekey/refbox/blob/main/config/species.yaml).
Override it with the `REFBOX_CONFIG` environment variable.

Each entry is **three levels deep**: species → assembly → resource. All resource
keys should appear; use `null` for resources with no upstream source.

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
      cytoband:
        url: https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/cytoBand.txt.gz
      hgnc:
        url: https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt
```

The bundled registry covers **26 species / 42 assemblies** (only a few are
`enabled: true` by default; pass `--include-disabled` to process the rest).

## Canonical resource names

| Name | `raw/` | `build/` |
|---|---|---|
| `genome` | `genome.fa` | `genome.fa.gz` + `.fai` + `.gzi`, `chrom.sizes` |
| `transcriptome` | `transcriptome.fa` | `transcriptome.fa.gz` + `.fai`; `transcriptome.derived.fa.gz` + `.fai` |
| `annotation_gtf` | `annotation_gtf.gtf` | `annotation.sorted.gtf.gz` + `.tbi` |
| `annotation_gff3` | `annotation_gff3.gff3` | `annotation.sorted.gff3.gz` + `.tbi` |
| `repeats_rmsk` | `repeats_rmsk.tsv` | `repeats.sorted.bed.gz` + `repeats.sorted.gtf.gz` |
| `repeats_bed` | `repeats_bed.bed` | `repeats.sorted.bed.gz` + `.tbi` |
| `repeats_gtf` | `repeats_gtf.gtf` | `repeats.sorted.gtf.gz` + `.tbi` |
| `repeats_fa` | `repeats_fa.fa` | *(RepeatMasker `.fa.out` report)* |
| `rnacentral` | `rnacentral.gff3` | `rnacentral.sorted.gff3.gz` + `.tbi` |
| `ccre` | `ccre.bed` | `ccre.sorted.bed.gz` + `.tbi` |
| `cytoband` | `cytoband.tsv` | `cytoband.sorted.bed.gz` + `.tbi`, `cytoband.bb` |

`cytoband` pulls UCSC `cytoBand.txt.gz` (real Giemsa bands) where available,
falling back to `cytoBandIdeo.txt.gz`. T2T / `hs1` has no text table — its
`cytoBandMapped.bb` bigBed is downloaded and expanded to the same TSV before
indexing. Both a tabix `cytoband.sorted.bed.gz` (region queries) and a
`cytoband.bb` bigBed (whole-chromosome ideograms) are produced.

## Fallbacks & special sources

### Local file before download

```yaml
my_resource:
  local_path: /path/on/disk/file.fa.gz     # copied if present
  url:        https://.../file.fa.gz        # else downloaded
```

### Concatenate multiple upstream files

```yaml
# e.g. Ensembl cdna + ncrna into one transcriptome raw
transcriptome:
  url: https://ftp.ensembl.org/pub/release-111/fasta/danio_rerio/cdna/Danio_rerio.GRCz11.cdna.all.fa.gz
  extra_urls:
    - https://ftp.ensembl.org/pub/release-111/fasta/danio_rerio/ncrna/Danio_rerio.GRCz11.ncrna.fa.gz
```

### RNAcentral cross-assembly liftover

```yaml
# no direct URL upstream → lift from another assembly
rnacentral:
  liftover_from:
    source_assembly: GRCh38
    url:       https://.../homo_sapiens.GRCh38.gff3.gz
    chain_url: https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz
```

### Transcriptome auto-derivation

```yaml
# leave it null and refbox builds transcriptome.fa.gz from genome + GTF via gffread
transcriptome: null
```

## Adding a new assembly

1. Add a `species → assembly → resource` block to `species.yaml` (or your own
   file pointed to by `REFBOX_CONFIG`).
2. Set `enabled: true`.
3. Fill in each resource URL (or `null` / `liftover_from` / `local_path`).
4. Run `refbox pull --assembly <NAME>`.
