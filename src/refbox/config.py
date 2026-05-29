"""Load species.yaml and provide iteration helpers."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterator

import yaml

# Canonical resource names — order matters for download/build progress display.
RESOURCE_NAMES = [
    "genome",
    "transcriptome",
    "annotation_gtf",
    "annotation_gff3",
    "repeats_rmsk",
    "repeats_bed",
    "repeats_gtf",
    "repeats_fa",
    "rnacentral",
    "ccre",
]

# File extension each resource is normalized to inside raw/
RESOURCE_EXT = {
    "genome":          "fa",
    "transcriptome":   "fa",
    "annotation_gtf":  "gtf",
    "annotation_gff3": "gff3",
    "repeats_rmsk":    "tsv",
    "repeats_bed":     "bed",
    "repeats_gtf":     "gtf",
    "repeats_fa":      "fa",
    "rnacentral":      "gff3",
    "ccre":            "bed",
}

import os

# Output root: where {Species}/{Assembly}/raw|build/ are written.
# Override with environment variable REFBOX_OUT or pass --out on the CLI.
_DEFAULT_OUT = Path(os.environ.get("REFBOX_OUT", Path.cwd())).resolve()

# Config search order:
#   1. $REFBOX_CONFIG env var
#   2. ./config/species.yaml (when running from a checkout)
#   3. bundled species.yaml inside the installed package
_PKG_DIR = Path(__file__).resolve().parent
_BUNDLED_CONFIG = _PKG_DIR / "config" / "species.yaml"
_REPO_CONFIG = _PKG_DIR.parents[1] / "config" / "species.yaml"  # only exists in source checkout


def _default_config_path() -> Path:
    env = os.environ.get("REFBOX_CONFIG")
    if env:
        return Path(env).expanduser().resolve()
    if _REPO_CONFIG.exists():
        return _REPO_CONFIG
    return _BUNDLED_CONFIG


DEFAULT_CONFIG = _default_config_path()


@dataclass
class Target:
    species: str
    assembly: str
    enabled: bool
    resources: dict[str, dict[str, Any] | None]   # resource_name -> spec or None
    meta: dict[str, Any]                          # extras like gencode_version
    out_root: Path = _DEFAULT_OUT

    @property
    def raw_dir(self) -> Path:
        return self.out_root / self.species / self.assembly / "raw"

    @property
    def build_dir(self) -> Path:
        return self.out_root / self.species / self.assembly / "build"

    def resource(self, name: str) -> dict[str, Any] | None:
        spec = self.resources.get(name)
        if spec is None:
            return None
        return spec


def load_config(path: str | Path = DEFAULT_CONFIG) -> dict[str, Any]:
    with open(path) as f:
        return yaml.safe_load(f)


def iter_targets(
    config: dict[str, Any] | None = None,
    species: list[str] | None = None,
    assembly: list[str] | None = None,
    include_disabled: bool = False,
    out_root: Path | str | None = None,
) -> Iterator[Target]:
    """Yield (species, assembly) Target objects honoring filters and enabled flag."""
    cfg = config if config is not None else load_config()
    root = Path(out_root).resolve() if out_root else _DEFAULT_OUT
    for sp_name, asm_map in cfg.get("species", {}).items():
        if species and sp_name not in species:
            continue
        for asm_name, asm_cfg in (asm_map or {}).items():
            if assembly and asm_name not in assembly:
                continue
            asm_cfg = asm_cfg or {}
            enabled = bool(asm_cfg.get("enabled", True))
            if not enabled and not include_disabled:
                continue
            resources: dict[str, dict[str, Any] | None] = {}
            meta: dict[str, Any] = {}
            for k, v in asm_cfg.items():
                if k == "enabled":
                    continue
                if k in RESOURCE_NAMES:
                    resources[k] = v
                else:
                    meta[k] = v
            yield Target(
                species=sp_name,
                assembly=asm_name,
                enabled=enabled,
                resources=resources,
                meta=meta,
                out_root=root,
            )


def raw_path(target: Target, resource: str) -> Path:
    """Canonical raw file path (uncompressed extension)."""
    ext = RESOURCE_EXT[resource]
    return target.raw_dir / f"{resource}.{ext}"


def find_species_by_assembly(
    assembly: str, config: dict[str, Any] | None = None
) -> str | None:
    """Return the species name that owns ``assembly``, or None if unknown."""
    cfg = config if config is not None else load_config()
    for sp_name, asm_map in cfg.get("species", {}).items():
        if assembly in (asm_map or {}):
            return sp_name
    return None
