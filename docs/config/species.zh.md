# `species.yaml` 注册表

驱动 `download` / `pull` / `publish` / `test` 的注册表位于
[`config/species.yaml`](https://github.com/typekey/refbox/blob/main/config/species.yaml)。
可用环境变量 `REFBOX_CONFIG` 覆盖它。

每个条目**三层深**：物种 → 基因组版本 → 资源。所有资源键都应出现；无上游来源的资源用
`null`。

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

内置注册表覆盖 **26 个物种 / 42 个基因组版本**（默认仅少数 `enabled: true`；传
`--include-disabled` 处理其余）。

## 标准资源名称

| 名称 | `raw/` | `build/` |
|---|---|---|
| `genome` | `genome.fa` | `genome.fa.gz` + `.fai` + `.gzi`，`chrom.sizes` |
| `transcriptome` | `transcriptome.fa` | `transcriptome.fa.gz` + `.fai`；`transcriptome.derived.fa.gz` + `.fai` |
| `annotation_gtf` | `annotation_gtf.gtf` | `annotation.sorted.gtf.gz` + `.tbi` |
| `annotation_gff3` | `annotation_gff3.gff3` | `annotation.sorted.gff3.gz` + `.tbi` |
| `repeats_rmsk` | `repeats_rmsk.tsv` | `repeats.sorted.bed.gz` + `repeats.sorted.gtf.gz` |
| `repeats_bed` | `repeats_bed.bed` | `repeats.sorted.bed.gz` + `.tbi` |
| `repeats_gtf` | `repeats_gtf.gtf` | `repeats.sorted.gtf.gz` + `.tbi` |
| `repeats_fa` | `repeats_fa.fa` | *（RepeatMasker `.fa.out` 报告）* |
| `rnacentral` | `rnacentral.gff3` | `rnacentral.sorted.gff3.gz` + `.tbi` |
| `ccre` | `ccre.bed` | `ccre.sorted.bed.gz` + `.tbi` |
| `cytoband` | `cytoband.tsv` | `cytoband.sorted.bed.gz` + `.tbi`，`cytoband.bb` |

`cytoband` 在可用时拉取 UCSC `cytoBand.txt.gz`（真实 Giemsa 带型），否则回退到
`cytoBandIdeo.txt.gz`。T2T / `hs1` 没有文本表 —— 会下载其 `cytoBandMapped.bb` bigBed 并在
建索引前展开为相同的 TSV。同时生成 tabix `cytoband.sorted.bed.gz`（区域查询）和
`cytoband.bb` bigBed（整条染色体示意图）。

## 回退与特殊来源

### 优先使用本地文件

```yaml
my_resource:
  local_path: /path/on/disk/file.fa.gz     # 存在则复制
  url:        https://.../file.fa.gz        # 否则下载
```

### 拼接多个上游文件

```yaml
# 例如把 Ensembl cdna + ncrna 拼成一个转录组原始文件
transcriptome:
  url: https://ftp.ensembl.org/pub/release-111/fasta/danio_rerio/cdna/Danio_rerio.GRCz11.cdna.all.fa.gz
  extra_urls:
    - https://ftp.ensembl.org/pub/release-111/fasta/danio_rerio/ncrna/Danio_rerio.GRCz11.ncrna.fa.gz
```

### RNAcentral 跨版本 liftover

```yaml
# 上游无直接 URL → 从其他版本 lift
rnacentral:
  liftover_from:
    source_assembly: GRCh38
    url:       https://.../homo_sapiens.GRCh38.gff3.gz
    chain_url: https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz
```

### 转录组自动推导

```yaml
# 留空为 null，refbox 会通过 gffread 从基因组 + GTF 构建 transcriptome.fa.gz
transcriptome: null
```

## 添加一个新基因组版本

1. 在 `species.yaml`（或由 `REFBOX_CONFIG` 指向的自定义文件）中添加一个
   `物种 → 版本 → 资源` 块。
2. 设置 `enabled: true`。
3. 填入每个资源的 URL（或 `null` / `liftover_from` / `local_path`）。
4. 运行 `refbox pull --assembly <NAME>`。
