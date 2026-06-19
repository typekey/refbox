# build 命令（`refbox build`）

`refbox build` 用于构建**任意用户提供的文件**，这些文件在 `species.yaml` 中没有条目 ——
自定义物种、内部组装、论文中预处理过的 BED 等。每次调用会按扩展名（或显式参数）识别输入，
校验 bgzip / 排序顺序，执行标准转换，并输出带索引的结果。

输入可以是明文或 `.gz`。构建器会透明地把明文 gzip 重新 bgzip、把已经 bgzip 的文件原样复制，
并在建索引前对未排序的 GFF/BED 重新排序。

## 用法

```text
refbox build [INPUT] [-o OUT] [输入类型参数] [选项]
```

## 输入选择器

用下列参数之一指定输入，或传入单个位置参数 `INPUT` 文件以按扩展名自动识别类型。

| 参数 | 输入 | 输出 |
|---|---|---|
| `-fa`、`--fa FILE` | 基因组 FASTA | `.fa.gz` + `.fai` + `.gzi` + `chrom.sizes` |
| `-gtf`、`--gtf FILE` | GTF 注释 | 排序后的 `.gtf.gz` + `.tbi` |
| `-gff`、`--gff FILE` | GFF3 注释 | 排序后的 `.gff3.gz` + `.tbi` |
| `-bed`、`--bed FILE` | BED 特征 | 排序后的 `.bed.gz` + `.tbi`（+ bigBed） |
| `-rmsk`、`--rmsk FILE` | UCSC `rmsk.txt[.gz]` | `repeats.sorted.bed.gz` + `repeats.sorted.gtf.gz`（+ `.tbi`） |
| `-rba`、`--rba FILE` | GTF/GFF3 → 独立 **完整索引**（`.rba`） | `.rba`（SQLite + FTS5 + 注释） |
| `-rbi`、`--rbi FILE` | GTF/GFF3 → 独立 **轻量索引**（`.rbi`） | `.rbi`（紧凑 B-tree 查找） |
| `-i`、`--ingest DIR` | 用户文件目录 | 完整的 `{Species}/{Assembly}/` 布局 |
| `INPUT`（位置参数） | 单个文件，自动识别 | 取决于识别出的类型 |
| `-fa … -gtf …` 同时给 | 基因组 + 注释 | 转录组 FASTA（经 gffread） |

!!! note "`.rba` 与 `.rbi` 是不同的构建器"
    `-rba` 构建完整的 FTS5 + 注释索引；`-rbi` 构建紧凑的轻量索引。二者**不是**别名。
    `-sqlite` / `--sqlite` 和 `--with-sqlite` 仍是 `-rba` / `--with-rba` 的别名。
    参见 [完整索引 `.rba`](../indexes/rba.md) 与 [轻量索引 `.rbi`](../indexes/rbi.md)。

## 选项

| 参数 | 适用于 | 含义 |
|---|---|---|
| `-o`、`--out PATH` | 全部 | 输出文件或目录 |
| `--with-rba` | `-gtf` / `-gff` | 在 tabix 输出旁额外生成完整 `.rba` 索引 |
| `--with-rbi` | `-gtf` / `-gff` | 在 tabix 输出旁额外生成轻量 `.rbi` 索引 |
| `--source-name NAME` | `.rba` / `.rbi` | 注释来源标签（如 `GENCODE`、`Ensembl`） |
| `--species-name NAME` | `.rba` / `.rbi` | 物种标签（写入元数据） |
| `--genome NAME` | `.rba` / `.rbi` | 基因组 / 版本标签（如 `hg38`） |
| `--annotation-version V` | `.rba` / `.rbi` | 注释版本（如 `v45`） |
| `--no-gram3` | `.rbi` 构建 | 跳过 3-gram 表（最小索引；仅精确 + 前缀） |
| `--synonyms`、`--synonyms-file FILE` | `.rba` 构建 | HGNC 风格 TSV，注入为 `gene_synonym` 别名 |
| `--rnacentral FILE` | `.rba` 构建 | RNAcentral 基因组坐标 GFF3，作为 ncRNA 合并进来 |
| `--fuzzy-scope {names,all}` | `.rba` 构建 | trigram 语料：`names`（默认）或 `all`（含 ID） |
| `--chrom-sizes FILE` | `-bed` | 用于 bigBed 转换的 chrom.sizes |
| `--assembly NAME` | `-bed`、`-i` | 版本标识（chrom.sizes 查找 / 目录名） |
| `--no-bigbed` | `-bed` | 跳过 bigBed 生成 |
| `--species NAME` | `-i` | 物种名（仅 ingest） |
| `--map RES:PATH ...` | `-i` | 资源 → 路径覆盖（仅 ingest） |
| `--no-build` | `-i` | 复制原始文件但跳过构建流水线（仅 ingest） |
| `--force` | 全部 | 覆盖已有输出 |

---

## 按输入类型分的示例

### 基因组 FASTA

```bash
# → bgzip + faidx + chrom.sizes（输出默认为 <input>.fa.gz）
refbox build -fa GENOME.fa
refbox build -fa GENOME.fa -o my_genome.fa.gz
```

### GTF / GFF3 注释

```bash
refbox build -gtf ANNOT.gtf                 # → ANNOT.sorted.gtf.gz + .tbi
refbox build -gff ANNOT.gff3 -o out.gff3.gz # → 排序 + bgzip + tabix

# 在 tabix 输出旁同时生成搜索索引
refbox build -gtf ANNOT.gtf --with-rba \
             --source-name GENCODE --genome hg38 --annotation-version v45
```

### 转录组（基因组 + 注释）

```bash
# 通过 gffread 提取转录本序列，再 bgzip + faidx
refbox build -fa GENOME.fa -gtf ANNOT.gtf -o transcriptome.fa.gz
```

### BED 特征

```bash
# chrom.sizes 由显式文件解析……
refbox build -bed FEATURES.bed --chrom-sizes hg38.chrom.sizes

# ……或由已知版本解析（委托给 zlbio 的物种级查找）
refbox build -bed FEATURES.bed --assembly GRCh38

# 跳过 bigBed（仅 tabix）
refbox build -bed FEATURES.bed --no-bigbed
```

### UCSC RepeatMasker 表

```bash
# rmsk.txt[.gz] → repeats.sorted.bed.gz + repeats.sorted.gtf.gz（+ .tbi）
refbox build -rmsk rmsk.txt.gz -o OUT_DIR
```

### 独立搜索索引

```bash
# GTF/GFF3 → 完整 RBrowser 注释索引（.rba：FTS5 + 注释结构）
refbox build -rba ANNOT.gtf.gz -o OUT.rba \
             --source-name GENCODE --species human --genome hg38 \
             --annotation-version v45 --force

# GTF/GFF3 → 轻量 RBrowser 索引（.rbi：紧凑的名称→位置查找）
refbox build -rbi ANNOT.gtf.gz -o OUT.rbi \
             --source-name GENCODE --species human --genome hg38 \
             --annotation-version v45 --force
```

详见 [完整索引 `.rba`](../indexes/rba.md) 与 [轻量索引 `.rbi`](../indexes/rbi.md)。

### 目录导入（ingest）

```bash
# 按扩展名自动归类一个目录的文件，构建为完整布局
refbox build -i DIR --assembly GRCh38 --species Homo_sapiens

# 覆盖某个文件对应哪个资源
refbox build -i DIR --assembly GRCh38 \
             --map genome:genome.fa --map annotation_gtf:genes.gtf

# 复制原始文件但跳过构建流水线
refbox build -i DIR --assembly GRCh38 --no-build
```

### 自动识别

```bash
# 按扩展名推断类型（rmsk/.out/.fa/.gtf/.gff3/.bed）
refbox build SOMEFILE
```
