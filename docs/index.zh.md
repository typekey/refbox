# refbox

**一条命令，为基因组浏览器构建标准化、带索引的参考文件。**

`refbox` 把一份描述物种 / 基因组版本的 YAML 注册表，转化为可直接加载到浏览器的输入文件，
并且能够构建静态的、无需后端的**搜索索引**，用于在浏览器内进行基因 / 转录本查询。

## 它能生成什么

| 资源 | 输出 |
|---|---|
| 基因组 FASTA | `bgzip` + `samtools faidx`（`.fa.gz` + `.fai` + `.gzi`）+ `chrom.sizes` |
| 转录组 FASTA | `bgzip` + `faidx`（当上游没有 URL 时，自动从基因组 + GTF 推导） |
| GTF / GFF3 注释 | 排序 + `bgzip` + `tabix` |
| 重复序列 | UCSC RepeatMasker `rmsk.txt.gz` → BED + GTF；`.fa.out.gz` 报告 |
| RNAcentral 非编码 RNA | 直接下载 *或* 从其他版本 liftover |
| ENCODE SCREEN cCRE | 排序 + `bgzip` + `tabix` |
| UCSC 染色体带型 | tabix BED + bigBed 染色体示意图 |
| **完整索引**（`.rba`） | SQLite + FTS5 + trigram + 注释结构 —— 精确 / 前缀 / 模糊 / 别名搜索与渲染 |
| **轻量索引**（`.rbi`） | 紧凑的 SQLite B-tree 查询 —— 快速的「名称 → 位置」自动补全 |

内置一份覆盖 **26 个物种 / 42 个基因组版本** 的注册表（人、小鼠、大鼠、狗、牛、猪、
黑猩猩、大猩猩、斑马鱼、果蝇、线虫、海胆、酵母、植物、细菌、病毒……），数据来源涵盖
GENCODE、Ensembl、Ensembl Genomes、UCSC golden path、NCBI、RNAcentral 和 ENCODE SCREEN。

## 命令一览

```text
refbox download   # 仅下载 species.yaml 中配置的原始文件
refbox pull       # 完整流水线：下载（缺失时）+ 构建 + 校验 + 发布
refbox publish    # 把已有的 build/ 目录扁平化为 <Assembly>.<name>
refbox test       # 校验 build/ 输出
refbox build      # 针对任意输入的单文件 / 单目录构建
```

## 接下来读什么

<div class="grid cards" markdown>

- :material-download: **[安装](installation.md)** —— pip + 外部生信工具。
- :material-rocket-launch: **[快速上手](quickstart.md)** —— 最常用的工作流。
- :material-console: **[流水线命令](commands/pipeline.md)** —— `download` / `pull` / `publish` / `test`。
- :material-file-cog: **[build 命令](commands/build.md)** —— 手动构建任意文件或目录。
- :material-database-search: **[完整索引 `.rba`](indexes/rba.md)** —— FTS5 搜索 + 注释结构。
- :material-database: **[轻量索引 `.rbi`](indexes/rbi.md)** —— 紧凑的名称→位置查询。

</div>

!!! tip "想知道怎么构建搜索索引？"
    简短回答：`refbox build -rba ANNOT.gtf.gz -o OUT.rba`（完整索引），或
    `refbox build -rbi ANNOT.gtf.gz -o OUT.rbi`（轻量索引）。详见
    **[完整索引](indexes/rba.md)** 和 **[轻量索引](indexes/rbi.md)**。
