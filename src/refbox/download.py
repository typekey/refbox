"""Download or copy raw resource files into {Species}/{Assembly}/raw/.

Policy per resource:
  - local_path exists -> copy/decompress into raw/<resource>.<ext>
  - else url is set   -> download (and decompress if .gz) into raw/<resource>.<ext>
  - else null         -> skip silently
"""

from __future__ import annotations

import gzip
import logging
import shutil
from pathlib import Path

import requests
from tqdm import tqdm

from .config import RESOURCE_NAMES, Target, iter_targets, raw_path

log = logging.getLogger(__name__)

CHUNK = 1 << 20  # 1 MiB


def _copy_or_decompress(src: Path, dst: Path) -> None:
    dst.parent.mkdir(parents=True, exist_ok=True)
    if str(src).endswith(".gz"):
        log.info("decompressing %s -> %s", src, dst)
        with gzip.open(src, "rb") as fin, open(dst, "wb") as fout:
            shutil.copyfileobj(fin, fout, length=CHUNK)
    else:
        log.info("copying %s -> %s", src, dst)
        shutil.copyfile(src, dst)


def _download(url: str, dst: Path) -> None:
    dst.parent.mkdir(parents=True, exist_ok=True)
    tmp = dst.with_suffix(dst.suffix + ".part")
    log.info("downloading %s -> %s", url, dst)
    with requests.get(url, stream=True, timeout=60) as r:
        r.raise_for_status()
        total = int(r.headers.get("Content-Length", 0))
        with open(tmp, "wb") as f, tqdm(
            total=total or None, unit="B", unit_scale=True, desc=dst.name
        ) as bar:
            for chunk in r.iter_content(chunk_size=CHUNK):
                f.write(chunk)
                bar.update(len(chunk))
    if url.endswith(".gz") and not str(dst).endswith(".gz"):
        log.info("decompressing %s -> %s", tmp, dst)
        with gzip.open(tmp, "rb") as fin, open(dst, "wb") as fout:
            shutil.copyfileobj(fin, fout, length=CHUNK)
        tmp.unlink()
    else:
        tmp.rename(dst)


def fetch_resource(target: Target, resource: str, *, force: bool = False) -> Path | None:
    spec = target.resource(resource)
    if spec is None:
        log.info("[%s/%s] skip %s (null)", target.species, target.assembly, resource)
        return None
    dst = raw_path(target, resource)
    if dst.exists() and not force:
        log.info("[%s/%s] %s already exists: %s",
                 target.species, target.assembly, resource, dst)
        return dst

    local = spec.get("local_path")
    url = spec.get("url")
    if local and Path(local).exists():
        _copy_or_decompress(Path(local), dst)
        return dst
    if url:
        _download(url, dst)
        return dst
    log.warning("[%s/%s] %s has neither readable local_path nor url",
                target.species, target.assembly, resource)
    return None


def download_targets(
    species: list[str] | None = None,
    assembly: list[str] | None = None,
    resources: list[str] | None = None,
    *,
    out: str | None = None,
    force: bool = False,
) -> None:
    resources = resources or RESOURCE_NAMES
    for tgt in iter_targets(species=species, assembly=assembly, out_root=out):
        log.info("=== download %s / %s ===", tgt.species, tgt.assembly)
        for r in resources:
            try:
                fetch_resource(tgt, r, force=force)
            except Exception as e:
                log.error("[%s/%s] %s FAILED: %s",
                          tgt.species, tgt.assembly, r, e)
