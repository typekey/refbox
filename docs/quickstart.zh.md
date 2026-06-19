# 快速上手

## 1. 端到端构建一个已配置的基因组版本

注册表驱动的流水线（`pull`）一条命令完成 **下载 → 构建 → 校验 → 发布**，使用内置的
`species.yaml`：

```bash
# 下载 + 构建 + 校验 + 发布 人类 GRCh38
refbox pull --species Homo_sapiens --assembly GRCh38

# --species 可省略；会从 --assembly 推断
refbox pull --assembly GRCm38
```

输出位于 `$REFBOX_OUT/<Species>/<Assembly>/`，以扁平化的 `<Assembly>.<name>` 文件呈现
（例如 `GRCh38.genome.fa.gz`）。参见 [输出布局](config/output-layout.md)。

## 2. 构建单个任意文件

无需 `species.yaml` 条目 —— 直接处理自己的数据：

```bash
# 基因组 FASTA → bgzip + faidx + chrom.sizes
refbox build -fa my_genome.fa

# GTF/GFF3 → 排序 + bgzip + tabix
refbox build -gtf my_annotation.gtf

# BED → 排序 + bgzip + tabix + bigBed
refbox build -bed peaks.bed --assembly GRCh38
```

完整说明见 [build 命令](commands/build.md)。

## 3. 构建搜索索引（最常被问到的问题）

=== "完整索引 (.rba)"

    ```bash
    # 独立的 RBrowser 索引（SQLite + FTS5 搜索）
    refbox build -rba gencode.v45.annotation.gtf.gz -o hg38.gencode.v45.rba \
        --source-name GENCODE --species human --genome hg38 \
        --annotation-version v45 --force
    ```

    或者在一条命令里与 tabix 输出**一起**生成：

    ```bash
    refbox build -gtf gencode.v45.annotation.gtf.gz --with-rba \
        --source-name GENCODE --genome hg38 --annotation-version v45
    ```

    详见 [RBrowser 索引 `.rba`](indexes/rba.md)。

=== "轻量索引 (.rbi)"

    ```bash
    refbox build -rbi gencode.v45.annotation.gtf.gz -o hg38.gencode.v45.rbi \
        --source-name GENCODE --species human --genome hg38 \
        --annotation-version v45 --force
    ```

    详见 [轻量索引 `.rbi`](indexes/rbi.md)。

## 4. 一键驱动脚本

`build.sh` 封装了 `refbox pull` 并提供合理的默认值（输出根目录 = 仓库的上级目录）：

```bash
./build.sh                            # 所有 enabled 的版本
./build.sh Homo_sapiens               # 一个物种（所有 enabled 的版本）
./build.sh Homo_sapiens GRCh38        # 一个物种 + 版本
./build.sh -- --resource genome cytoband   # `--` 之后的参数会传给 refbox

# 环境变量开关
FORCE=1 ./build.sh Homo_sapiens GRCh38     # 重新下载 + 重新构建
NO_TEST=1 ./build.sh                       # 跳过校验
INCLUDE_DISABLED=1 ./build.sh              # 同时处理 enabled: false 的条目
```
