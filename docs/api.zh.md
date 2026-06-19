# Python API

CLI 能做的一切，也都可以通过编程方式完成。

## 导入

```python
from refbox.config      import load_config, iter_targets, RESOURCE_NAMES
from refbox.download    import download_targets
from refbox.build       import build_targets, publish_targets
from refbox.test        import test_targets
from refbox             import file_build as fb        # 单文件构建器
from refbox.ingest      import ingest_directory        # 目录导入
from refbox.report      import build_report            # Markdown 状态报告
from refbox.sqlite_index import (                      # 完整 RBrowser 索引 (.rba)
    build_sqlite_index, search, inspect, open_readonly, normalize)
from refbox.lite_index  import (                       # 轻量索引 (.rbi)
    build_lite_index, search as lite_search, open_readonly as lite_open, inspect as lite_inspect)
```

## 注册表流水线

```python
# 等价于：refbox pull --assembly GRCh38
build_targets(species=["Homo_sapiens"], assembly=["GRCh38"], auto_download=True)
test_targets(species=["Homo_sapiens"], assembly=["GRCh38"])
publish_targets(species=["Homo_sapiens"], assembly=["GRCh38"])
```

## 单文件构建器

```python
from pathlib import Path
from refbox import file_build as fb

fb.build_fa(Path("genome.fa"))                                   # → genome.fa.gz + .fai + .gzi
fb.build_gxf(Path("annot.gtf"))                                  # → annot.sorted.gtf.gz + .tbi
fb.build_bed(Path("peaks.bed"), assembly="GRCh38")              # → 排序 + tabix + bigBed
fb.build_rmsk(Path("rmsk.txt.gz"))                              # → repeats.sorted.{bed,gtf}.gz
fb.build_transcriptome(Path("genome.fa"), Path("annot.gtf"))   # → transcriptome.fa.gz
```

## 构建并查询完整索引（`.rba`）

```python
from refbox.sqlite_index import build_sqlite_index, open_readonly, search

db = build_sqlite_index("gencode.v45.annotation.gtf.gz", "idx.rba",
                        source_name="GENCODE", genome="hg38", force=True)
con = open_readonly(db)
for hit in search(con, "TP53", limit=10):
    print(hit["matched_field"], hit["gene_name"], hit["transcript_id"])
```

`build_sqlite_index` 签名：

```python
build_sqlite_index(
    input_path, output=None, *,
    source_name="", species="", genome="", annotation_version="",
    synonyms=None,          # HGNC 风格 TSV 路径
    rnacentral=None,        # RNAcentral 基因组坐标 GFF3 路径
    fuzzy_scope="names",    # 或 "all"
    force=False, verbose=False,
) -> Path
```

## 构建并查询轻量索引（`.rbi`）

```python
from refbox.lite_index import build_lite_index, open_readonly, search

db = build_lite_index("gencode.v45.annotation.gtf.gz", "idx.rbi",
                      source_name="GENCODE", genome="hg38",
                      enable_gram3=True, force=True)
con = open_readonly(db)
mode, results = search(con, "TP53", limit=10)   # 注意：返回 (mode, results)
for hit in results:
    print(hit["gene_name"], hit["chrom"], hit["start"], hit["end"])
```

`build_lite_index` 签名：

```python
build_lite_index(
    input_path, output=None, *,
    source_name="", species="", genome="", annotation_version="",
    enable_gram3=False, force=False, verbose=False,
) -> Path
```

!!! note "参考搜索实现"
    浏览器端的搜索应当镜像 `refbox.sqlite_index.search()`（针对 `.rba`，返回命中列表）和
    `refbox.lite_index.search()`（针对 `.rbi`，返回 `(mode, results)` 元组）中的分层排序
    逻辑 —— 这两个 Python 函数就是参考实现。
