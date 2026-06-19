# 流水线命令

这四个子命令操作 **`species.yaml` 注册表** —— 它们知道物种、基因组版本以及 11 个标准资源。
若要构建任意用户文件，请参见 [build 命令](build.md)。

```text
refbox download   # 仅下载 species.yaml 中配置的原始文件
refbox pull       # 完整流水线：下载（缺失时）+ 构建 + 校验 + 发布
refbox publish    # 把已有的 build/ 目录扁平化为 <Assembly>.<name>
refbox test       # 校验 build/ 输出
```

## 通用过滤参数

每个流水线命令都接受相同的选择性过滤参数：

| 参数 | 含义 |
|---|---|
| `--species [NAME ...]` | 物种过滤（可选；可从 `--assembly` 推断） |
| `--assembly [NAME ...]` | 版本过滤（如 `GRCh38`）；省略则运行匹配 `--species` 的**全部**版本（都不给则运行所有） |
| `--resource [NAME ...]` | 资源子集（见下表） |
| `--out DIR` | 输出根目录（默认：`$REFBOX_OUT` 或当前目录） |
| `--include-disabled` | 同时处理 `species.yaml` 中标记为 `enabled: false` 的版本 |
| `-v`、`--verbose` | DEBUG 级别日志（全局参数，放在子命令之前） |

**资源名称**（`--resource`）：`genome` · `transcriptome` · `annotation_gtf` ·
`annotation_gff3` · `repeats_rmsk` · `repeats_bed` · `repeats_gtf` · `repeats_fa` ·
`rnacentral` · `ccre` · `cytoband`。

---

## `refbox pull`

注册表驱动的流水线：**下载 → 构建 → 校验 → 发布**。这是大多数用户想要的。

默认情况下，每个版本的 `build/` 输出会被**扁平化**为发布布局
`<out>/<Species>/<Assembly>/<Assembly>.<name>`，并删除 `build/` + `raw/` 工作目录。
传入 `--no-flat` 可保留它们。

### 额外参数

| 参数 | 含义 |
|---|---|
| `--force` | 即使输出已存在也重新构建 |
| `--no-download` | 跳过自动下载阶段（使用已有的 `raw/`） |
| `--no-test` | 跳过构建后的校验步骤 |
| `--no-flat` | 保留 `build/` + `raw/` 布局（跳过发布 / 扁平化） |

### 示例

```bash
# 下载 + 构建 + 校验 + 发布 人类 GRCh38（内置 species.yaml）
refbox pull --species Homo_sapiens --assembly GRCh38

# --species 从 --assembly 经注册表推断
refbox pull --assembly GRCm38

# 注册表中所有版本（含 enabled: false），指定资源
refbox pull --include-disabled \
            --resource genome transcriptome annotation_gtf annotation_gff3 \
                       repeats_rmsk rnacentral cytoband

# 只构建一种资源
refbox pull --assembly GRCh38 --resource ccre

# 保留 build/ + raw/ 工作目录，而非扁平布局
refbox pull --assembly GRCh38 --no-flat

# 从头重新构建
refbox pull --assembly GRCh38 --force
```

如果校验发现问题（且未给 `--no-test`），退出码为 `1`，否则为 `0`。

---

## `refbox download`

仅下载 `species.yaml` 中配置的原始文件并**停止** —— 不构建。

| 参数 | 含义 |
|---|---|
| 通用过滤参数 | 见上文 |
| `--force` | 即使原始文件已存在也重新下载 |

下载后端按顺序尝试 `axel → aria2c → wget → wget --no-check-certificate → requests →
requests (verify=False)`，因此单个 TLS 损坏的主机不会中断某个资源。

```bash
refbox download --assembly GRCh38                      # 仅原始文件
refbox download --assembly GRCh38 --resource genome    # 一种资源
refbox download --include-disabled                     # 全部
```

---

## `refbox publish`

把已有的 `build/` 目录扁平化为 `<Assembly>.<name>` 并删除 `raw/`。除非给了 `--no-flat`，
否则 `pull` 会自动执行此步骤 —— 单独运行它可在 `--no-flat` 构建后再做扁平化。

| 参数 | 含义 |
|---|---|
| 通用过滤参数 | 见上文 |
| `--keep-build` | 不删除（已清空的）`build/` 目录 |
| `--keep-raw` | 不删除 `raw/` 目录 |

```bash
refbox publish --assembly GRCh38                  # 一个版本
refbox publish --include-disabled                 # 注册表中全部
refbox publish --assembly GRCh38 --keep-raw       # 保留 raw/ 下载
```

---

## `refbox test`

针对已有的 `build/` 输出重新运行校验器。检查每个预期的 `.gz` 及其索引
（`.fai`/`.gzi`/`.tbi`）是否存在、示例 tabix 区域查询是否返回行、`samtools faidx`
是否返回序列。

```bash
refbox test --assembly GRCh38           # 校验一个版本
refbox test --include-disabled          # 校验注册表中全部
```

任一检查失败退出码为 `1`，否则为 `0`。
