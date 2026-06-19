# 轻量索引（`.rbi`）

完整的 [`.rba` 索引](rba.md)（SQLite + FTS5 + trigram + 注释结构）对 GENCODE 而言约 1 GB。
当你只需要**基因/转录本名称 → 位置**的自动补全 —— 而不需要子串/模糊别名搜索或完整的转录本
结构 —— **轻量 RBrowser 索引**（`.rbi`）只占其中很小一部分（几十 MB），并通过纯 B-tree
查找解析查询。

| | |
|---|---|
| **名称** | RBrowser Index |
| **简称** | RBI |
| **文件扩展名** | `.rbi` |
| **存储引擎** | SQLite（默认无 FTS5、无 trigram） |
| **搜索** | 在归一化 `term` 表上做 B-tree 范围扫描（精确 / 前缀 / 自动补全）；可选 3-gram 表用于模糊召回 |
| **存储内容** | 每个基因/转录本的显示名 + 基因组位置（无外显子/CDS 结构，无 payload JSON） |

## 如何构建

```bash
# 独立的轻量索引
refbox build -rbi gencode.v45.annotation.gtf.gz -o hg38.gencode.v45.rbi \
    --source-name GENCODE --species human --genome hg38 \
    --annotation-version v45 --force
```

或在一条命令里与 tabix 注释**一起**生成：

```bash
refbox build -gtf gencode.v45.annotation.gtf.gz --with-rbi \
    --source-name GENCODE --genome hg38 --annotation-version v45
```

### 参数

| 参数 | 含义 |
|---|---|
| `-rbi FILE` | GTF/GFF3 输入 → 构建独立 `.rbi`（默认输出 `<stem>.rbi`） |
| `--with-rbi` | 配合 `-gtf`/`-gff`：在 sorted/bgzip/tabix 输出旁额外生成 `.rbi` |
| `-o PATH` | 输出路径 |
| `--source-name` / `--species-name` / `--genome` / `--annotation-version` | 写入索引的元数据标签 |
| `--no-gram3` | 跳过 3-gram 表（最小索引；仅精确 + 前缀查找） |
| `--force` | 覆盖已有输出 |

!!! note "3-gram 模糊召回默认开启"
    默认情况下 `.rbi` 会包含基于基因/转录本**名称**的 3-gram 表，支持子串/模糊召回
    （例如 `alat` → `MALAT1`）。当你只需要精确 + 前缀查找、希望文件尽可能小时，传
    `--no-gram3`。

## 表结构

| 表 | 用途 |
|---|---|
| `record` | 每个基因/转录本一行 —— 仅显示名 + 基因组位置 |
| `term` | `(term_norm, term_raw, field, record_id, priority)` `WITHOUT ROWID`；主键以 `term_norm` 开头，使精确（`= ?`）和前缀（`>= lo AND < hi`）都是 B-tree 范围扫描 |
| `gram3` | 可选的 3-gram → 记录映射，用于模糊/部分召回（默认开启；`--no-gram3` 跳过） |
| `metadata` | 格式（`RBI` / "RBrowser Index"）+ 来源信息（来源 / 物种 / 基因组 / 版本 / 计数） |

### term 归一化

对于每个原始值，构建器在 `term` 中最多存储三个查找键：

1. **归一化** —— 去空白 + 小写（保留分隔符）：`TP53` → `tp53`
2. **去版本** —— 去掉末尾的 Ensembl `.<version>`：`ENST…305.10` → `enst…305`
3. **去分隔符** —— 去掉 `_ - . 空格`：`TP53-201` → `tp53201`

字段优先级决定排序（越小越靠前）：`gene_name`（10）→ `gene_id`（20）→ `transcript_name`（30）
→ `transcript_id`（40）→ `transcript_biotype`（80）。只有名称会生成 3-gram —— 对登录号或
biotype 做模糊匹配毫无意义，还会让表膨胀。

### 搜索层级

`search()` 在第一个产生结果的层级停止：

```text
1. 精确 term_norm
2. 前缀 B-tree 范围
3. 去分隔符精确
4. 去版本精确
5. gram3 模糊召回（若该表存在）
```

绝不会运行全表 `LIKE '%q%'`。

## 编程使用

```python
from refbox.lite_index import build_lite_index, open_readonly, search

db = build_lite_index("gencode.v45.annotation.gtf.gz", "hg38.gencode.v45.rbi",
                      source_name="GENCODE", genome="hg38", force=True)
con = open_readonly(db)
mode, results = search(con, "TP53", limit=10)   # mode: exact / prefix / normalized / gram3 / none
for hit in results:
    print(hit["gene_name"], hit["chrom"], hit["start"], hit["end"])
```

## `.rba` 还是 `.rbi`？

| 需求 | 选择 |
|---|---|
| 子串/模糊别名搜索、同义词、RNAcentral ncRNA、需要完整转录本结构（外显子/CDS/UTR）以无额外请求渲染 | **`.rba`** |
| 快速的名称→位置自动补全、最小的静态文件 | **`.rbi`** |
