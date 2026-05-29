"""Scan a refbox output tree and produce a Chinese Markdown status report.

Usage:
    python -m refbox.report --out /path/to/reference > report.md

For each species/assembly under ``--out``, reports the presence of the five
build artifacts requested by the user (genome / transcriptome / annotation
GTF+GFF3 / repeats / rnacentral) and the file sizes.
"""

from __future__ import annotations

import argparse
from pathlib import Path

from .config import load_config


# (display_name, build_dir filename used as the "primary" artifact, [index siblings])
ARTIFACTS = [
    ("基因组 (Genome)",        "genome.fa.gz",              [".fai", ".gzi"]),
    ("转录组 (Transcriptome)", "transcriptome.fa.gz",         [".fai"]),
    ("转录组 derived (gffread)", "transcriptome.derived.fa.gz", [".fai"]),
    ("基因注释 GTF",           "annotation.sorted.gtf.gz",  [".tbi"]),
    ("基因注释 GFF3",          "annotation.sorted.gff3.gz", [".tbi"]),
    ("重复元件 BED",           "repeats.sorted.bed.gz",     [".tbi"]),
    ("重复元件 GTF",           "repeats.sorted.gtf.gz",     [".tbi"]),
    ("RNAcentral",            "rnacentral.sorted.gff3.gz", [".tbi"]),
]


def _human(n: int) -> str:
    for unit in ("B", "KB", "MB", "GB", "TB"):
        if n < 1024:
            return f"{n:.1f} {unit}"
        n /= 1024
    return f"{n:.1f} PB"


def _status(build_dir: Path, fname: str, idx_suffixes: list[str]) -> tuple[str, str]:
    """Return (status emoji+label, size string) for one artifact."""
    main = build_dir / fname
    if not main.exists() or main.stat().st_size == 0:
        return ("✗ 缺失", "—")
    for sfx in idx_suffixes:
        idx = Path(str(main) + sfx)
        if not idx.exists():
            return ("⚠ 缺索引", _human(main.stat().st_size))
    return ("✓ 完成", _human(main.stat().st_size))


def build_report(out_root: Path) -> str:
    cfg = load_config()
    lines: list[str] = []
    lines.append("# refbox 参考数据构建报告")
    lines.append("")
    lines.append(f"- 输出根目录：`{out_root}`")
    lines.append(f"- refbox 配置：26 物种 / 42 装配版本")
    lines.append("- 资源类别：基因组、转录组、基因注释 (GTF/GFF3)、重复元件、RNAcentral")
    lines.append("")
    lines.append("## 统计")
    lines.append("")

    # ── per-assembly detail tables ────────────────────────────────────────────
    detail_blocks: list[str] = []
    overall_done = 0
    overall_total = 0
    for sp, asm_map in cfg.get("species", {}).items():
        for asm_name in (asm_map or {}):
            build_dir = out_root / sp / asm_name / "build"
            raw_dir = out_root / sp / asm_name / "raw"

            rows = ["| 资源类别 | 状态 | 大小 | 文件 |", "| --- | --- | --- | --- |"]
            n_done = 0
            for label, fname, idx in ARTIFACTS:
                st, sz = _status(build_dir, fname, idx)
                if st.startswith("✓"):
                    n_done += 1
                file_path = (build_dir / fname).relative_to(out_root)
                rows.append(f"| {label} | {st} | {sz} | `{file_path}` |")
            overall_done += n_done
            overall_total += len(ARTIFACTS)

            header = f"### {sp} / {asm_name}  ({n_done}/{len(ARTIFACTS)})"
            tree = []
            if build_dir.exists():
                tree.append(f"- `build/` ：{sum(1 for _ in build_dir.iterdir())} 个文件")
            else:
                tree.append("- `build/` ：未创建")
            if raw_dir.exists():
                tree.append(f"- `raw/`   ：{sum(1 for _ in raw_dir.iterdir())} 个文件")
            else:
                tree.append("- `raw/`   ：未创建")

            detail_blocks.append("\n".join([header, "", *tree, "", *rows, ""]))

    pct = (overall_done / overall_total * 100) if overall_total else 0.0
    lines.append(f"- 总条目：{overall_total}")
    lines.append(f"- 已完成：{overall_done}（{pct:.1f}%）")
    lines.append("")
    lines.append("## 目录结构")
    lines.append("")
    lines.append("```")
    lines.append(f"{out_root.name}/")
    lines.append("├── {Species}/")
    lines.append("│   └── {Assembly}/")
    lines.append("│       ├── raw/        # 原始下载文件 (genome.fa, annotation_gtf.gtf, ...)")
    lines.append("│       └── build/      # 已索引的最终产物 (genome.fa.gz + .fai + .gzi, ...)")
    lines.append("├── _logs/              # 每个 assembly 一份 refbox pull 日志")
    lines.append("└── report.md           # 本报告")
    lines.append("```")
    lines.append("")
    lines.append("## 逐 Assembly 明细")
    lines.append("")
    lines.extend(detail_blocks)
    return "\n".join(lines)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(prog="refbox.report")
    ap.add_argument("--out", required=True, help="refbox output root")
    args = ap.parse_args(argv)
    print(build_report(Path(args.out).resolve()))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
