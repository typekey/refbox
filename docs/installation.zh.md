# 安装

## 安装 Python 包

```bash
pip install refbox
```

需要 **Python ≥ 3.9**。该包会自动安装 `pyyaml`、`requests`、`tqdm` 和 `pysam`。

### 从源码安装（开发）

```bash
git clone https://github.com/typekey/refbox.git
cd refbox
pip install -e .
pytest -q          # refbox build 单文件模式的单元测试
refbox --help
```

## 外部命令行工具

`refbox` 会调用标准生信工具。用 conda / mamba 一次性安装：

```bash
mamba install -c bioconda htslib samtools gffread ucsc-bedtobigbed
```

| 工具 | 用途 | 何时需要 |
|---|---|---|
| `bgzip` | 块状 gzip 压缩 | 每次构建 |
| `tabix` | 对 GFF/BED 建立位置索引 | 每次注释 / BED 构建 |
| `samtools` | 对 FASTA 执行 `faidx` | 基因组 / 转录组 |
| `gffread` | 从基因组 + GTF 提取转录本 | 推导转录组 |
| `liftOver` | RNAcentral 跨版本 liftover | 使用 `liftover_from` 的版本 |
| `bedToBigBed` / `bigBedToBed` | bigBed 转换 / 展开 | `refbox build -bed`、cytoband |
| GNU `sort`、`grep` | 排序 / 过滤 | 每次构建 |
