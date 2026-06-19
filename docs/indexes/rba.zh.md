# 完整索引（`.rba`）

`tabix` 回答的是*位置*范围查询；它无法回答「**TP53 是什么？**」。为此，refbox 可以从
GTF/GFF3 构建一个独立的、只读的 **完整 RBrowser 注释索引** —— 一个 `.rba` 文件。除了搜索
索引，它还包含**注释内容**（每个基因/转录本的外显子/CDS/UTR 结构），既支持精确 / 前缀自动
补全 / 模糊子串 / 别名搜索，也支持在浏览器内无额外请求地渲染。

| | |
|---|---|
| **名称** | RBrowser Indexed Annotation Database |
| **文件扩展名** | `.rba` |
| **存储引擎** | SQLite |
| **搜索引擎** | SQLite FTS5（prefix + trigram） |
| **访问方式** | 静态文件 + SQLite WASM / HTTP Range VFS（无后端） |
| **用途** | 基因 / 转录本 / 别名搜索 **加** 用于渲染的注释结构 |

如果只需要「名称→位置」查找的更小索引（无 FTS5、无注释结构），见
[轻量索引 `.rbi`](rbi.md)。

## 如何构建

=== "独立构建"

    ```bash
    refbox build -rba gencode.v45.annotation.gtf.gz -o hg38.gencode.v45.transcript.rba \
        --source-name GENCODE --species human --genome hg38 \
        --annotation-version v45 --force
    ```

=== "与 tabix 一起"

    ```bash
    refbox build -gtf gencode.v45.annotation.gtf.gz --with-rba \
        --source-name GENCODE --species human --genome hg38 \
        --annotation-version v45
    ```

!!! tip "别名与默认名"
    `-sqlite` / `--with-sqlite` 是 `-rba` / `--with-rba` 的别名。默认输出名为
    `<input-stem>.rba`。

    注意：`-rbi` **不是** `-rba` 的别名 —— 它构建独立的
    [轻量索引 `.rbi`](rbi.md)。

## 元数据参数

均为可选；它们会原样存入索引的 `metadata` 表，供浏览器展示来源信息。

| 参数 | 示例 | 含义 |
|---|---|---|
| `--source-name` | `GENCODE` | 注释来源标签 |
| `--species-name` | `human` | 物种标签 |
| `--genome` | `hg38` | 基因组 / 版本标签 |
| `--annotation-version` | `v45` | 注释版本 |

## 同义词（HGNC）

GENCODE/Ensembl 注释**不**包含常用基因同义词（`OCT4`、`p53`、`Nanog`……）。用 `--synonyms`
传入 HGNC 的 `hgnc_complete_set.txt`，即可把 `alias_symbol` / `prev_symbol` 注入为可搜索的
`gene_synonym` 别名（按 Ensembl ID 匹配，回退到 symbol）。于是 `OCT4` → POU5F1，
`OTF3` → POU5F1。

```bash
refbox build -rba gencode.v45.annotation.gtf.gz -o OUT.rba \
             --synonyms hgnc_complete_set.txt --force
```

## RNAcentral 非编码 RNA

`--rnacentral` 会把一个 RNAcentral 基因组坐标 GFF3（请用染色体名**已归一化**的文件，使坐标
显示一致）作为额外的 ncRNA 记录合并进来。它们可按完整 ID（`URS…_9606.N`）、URS 登录号
（`URS…_9606`）和去版本 URS（`URS000035F234`）搜索。

```bash
# GENCODE 基因/转录本 + RNAcentral ncRNA + HGNC 同义词 → 一个索引
refbox build -rba gencode.v45.annotation.gff3.gz -o human.rba \
             --rnacentral rnacentral.normalized.gff3 \
             --synonyms hgnc_complete_set.txt --force
```

RNAcentral 唯一的标签是一段冗长的自由文本 `description`，因此显示名会被提炼成简短、易识别的
符号 —— `DDX11L11`、`miR-34a-5p`、`pre-mir-571`、`piR-hsa-4818588`、`5S-rRNA`、
`Metazoan-SRP-RNA`。完整原始描述会保留为可搜索别名，所以召回率不变。

## 模糊范围（fuzzy scope）

`--fuzzy-scope` 控制 trigram（子串）索引覆盖的内容：

| 取值 | 行为 |
|---|---|
| `names`（默认） | trigram 只覆盖基因/转录本**名称 + 同义词**；ID 仅精确/前缀。trigram 体积约缩小 6×，消除病态的数字 ID 子串全表扫描。 |
| `all` | trigram 语料中包含 ID（旧行为）。 |

## 它索引了什么

每个**基因**和**转录本**各一行 `feature`，包含完整结构（外显子起止、CDS 跨度、5′/3′ UTR
跨度、biotype、来源）、一个 `search_text` 文本块，以及一个 `payload_json` 供浏览器无需额外
join 即可渲染。每个可搜索同义词成为一行 `alias`：名称与 ID、**去版本** Ensembl ID
（`ENST00000269305.9` → `ENST00000269305`）、HAVANA / CCDS / HGNC / 蛋白 ID、GFF3 `Alias` /
`Dbxref`（RefSeq）/ `gene_synonym`。

| 表 | 用途 |
|---|---|
| `feature` | 基因 + 转录本记录（1-based 闭区间坐标；`chrom_start0/end0` 给出 0-based 半开区间） |
| `alias` | `(feature_id, alias, alias_norm, alias_type, source)` —— `alias_norm` 为小写、去版本、去分隔符形式 |
| `metadata` | `key/value` 来源信息（来源、物种、基因组、坐标约定、计数、能力） |
| `feature_fts` | FTS5 前缀/自动补全（`prefix='2 3 … 10'`） |
| `feature_trigram` | FTS5 `trigram` 分词器，用于子串/模糊搜索（不可用时优雅退化为 `LIKE`） |

## 搜索排序

```text
精确 transcript_id → transcript_name → gene_name → gene_id → alias 精确
→ FTS 前缀 → trigram 子串 → LIKE 回退
```

浏览器应运行同样的层级；`refbox.sqlite_index.search()` 是参考实现。精确匹配列以
`COLLATE NOCASE` 建索引，因此查询是索引查找而非全表扫描。

## 坐标约定

`start` / `end` 以及所有 `*_start` / `*_end` 列为**1-based 闭区间**（GTF/GFF 约定），用于
显示；`chrom_start0` / `chrom_end0` 额外给出**0-based 半开区间**用于渲染。该约定记录在
`metadata` 的 `coord_convention` 中。

## 检查与查询

```python
from refbox.sqlite_index import build_sqlite_index, open_readonly, search, inspect

print(inspect("idx.rba")["metadata"])          # 来源信息 + 能力
con = open_readonly("idx.rba")
for hit in search(con, "TP53", limit=10):
    print(hit["matched_field"], hit["gene_name"], hit["transcript_id"])
```
