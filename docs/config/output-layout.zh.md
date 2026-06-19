# 输出布局

## 工作目录（构建过程中）

构建过程中，每个版本都有 `raw/`（下载）和 `build/`（可加载到浏览器）两个工作目录。使用
`refbox pull --no-flat` 或 `refbox download` 会得到这种布局。

```text
{REFBOX_OUT}/
  {Species}/
    {Assembly}/
      raw/                          # 原始下载 / 复制
        genome.fa
        transcriptome.fa             # 可能由基因组 + GTF 推导
        annotation_gtf.gtf
        annotation_gff3.gff3
        repeats_rmsk.tsv
        repeats_fa.fa
        rnacentral.gff3              # 可能从其他版本 lift
        ccre.bed
        cytoband.tsv                 # UCSC cytoBand[Ideo] 文本，或 bigBed（hs1）
      build/                        # 可加载到浏览器
        genome.fa.gz                 + .fai + .gzi
        chrom.sizes
        transcriptome.fa.gz          + .fai
        transcriptome.derived.fa.gz  + .fai   # gffread 提取，当有 GTF 时
        annotation.sorted.gtf.gz     + .tbi
        annotation.sorted.gff3.gz    + .tbi
        repeats.sorted.bed.gz        + .tbi
        repeats.sorted.gtf.gz        + .tbi
        rnacentral.sorted.gff3.gz    + .tbi
        ccre.sorted.bed.gz           + .tbi
        cytoband.sorted.bed.gz       + .tbi
        cytoband.bb                          # bigBed（bed4+1，含 gieStain）
```

## 发布布局（默认）

`pull` 以 **publish** 步骤收尾（用 `--no-flat` 跳过），它把每个 `build/<name>` 扁平化为
版本目录下的 `<Assembly>.<name>`，并删除 `build/` / `raw/`：

```text
{REFBOX_OUT}/
  {Species}/
    {Assembly}/
      {Assembly}.genome.fa.gz            + .fai + .gzi
      {Assembly}.chrom.sizes
      {Assembly}.annotation.sorted.gtf.gz   + .tbi
      {Assembly}.annotation.sorted.gff3.gz  + .tbi
      {Assembly}.repeats.sorted.bed.gz      + .tbi
      {Assembly}.repeats.sorted.gtf.gz      + .tbi
      {Assembly}.rnacentral.sorted.gff3.gz  + .tbi
      {Assembly}.ccre.sorted.bed.gz         + .tbi
      {Assembly}.cytoband.sorted.bed.gz     + .tbi
      {Assembly}.cytoband.bb
```

## 环境变量

| 变量 | 含义 |
|---|---|
| `REFBOX_OUT` | `{Species}/{Assembly}/{raw,build}/` 的默认输出根目录（否则为当前目录） |
| `REFBOX_CONFIG` | 自定义 `species.yaml` 的路径（覆盖内置注册表） |

## 状态报告器

```bash
python -m refbox.report --out /path/to/reference > report.md
```

遍历输出目录树，生成一份 Markdown 报告，按版本列出每个产物的状态
（✓ 完成 / ⚠ 缺索引 / ✗ 缺失）及大小。
