"""Download or copy raw resource files into {Species}/{Assembly}/raw/.

Policy per resource:
  - local_path exists -> copy/decompress into raw/<resource>.<ext>
  - else url is set   -> download (and decompress if .gz) into raw/<resource>.<ext>
  - else null         -> skip silently

Download backend priority (first available wins):
  1. axel     -- fastest, 16 parallel connections (-n REFBOX_CONNECTIONS)
  2. aria2c   -- multi-connection, resumable
  3. wget     -- single-connection, reliable
  4. requests -- built-in pure-Python fallback
"""

from __future__ import annotations

import gzip
import logging
import os
import shutil
import subprocess
from pathlib import Path

import requests
from tqdm import tqdm

from .config import RESOURCE_NAMES, Target, iter_targets, raw_path

log = logging.getLogger(__name__)

CHUNK = 1 << 20  # 1 MiB
# Parallel connections for axel / aria2c; override with REFBOX_CONNECTIONS=N
CONNECTIONS = int(os.environ.get("REFBOX_CONNECTIONS", "16"))


# ── backend detection ─────────────────────────────────────────────────────────

def _has(tool: str) -> bool:
    return shutil.which(tool) is not None


def _download_backend() -> str:
    """Return the best available CLI download backend name."""
    for tool in ("axel", "aria2c", "wget"):
        if _has(tool):
            return tool
    return "requests"


# ── local copy ────────────────────────────────────────────────────────────────

def _copy_or_decompress(src: Path, dst: Path) -> None:
    dst.parent.mkdir(parents=True, exist_ok=True)
    if str(src).endswith(".gz"):
        log.info("decompressing %s -> %s", src, dst)
        with gzip.open(src, "rb") as fin, open(dst, "wb") as fout:
            shutil.copyfileobj(fin, fout, length=CHUNK)
    else:
        log.info("copying %s -> %s", src, dst)
        shutil.copyfile(src, dst)


# ── CLI download wrappers ─────────────────────────────────────────────────────

def _run(cmd: list[str]) -> None:
    log.info("$ %s", " ".join(cmd))
    subprocess.run(cmd, check=True)


def _dl_axel(url: str, tmp: Path) -> None:
    _run(["axel", "-n", str(CONNECTIONS), "-a", "-o", str(tmp), url])


def _dl_aria2c(url: str, tmp: Path) -> None:
    _run([
        "aria2c",
        "-x", str(CONNECTIONS), "-s", str(CONNECTIONS),
        "-k", "10M",
        "--file-allocation=none",
        "-o", tmp.name, "-d", str(tmp.parent),
        url,
    ])


def _dl_wget(url: str, tmp: Path) -> None:
    _run(["wget", "--no-verbose", "--show-progress", "-O", str(tmp), url])


def _dl_requests(url: str, tmp: Path) -> None:
    log.info("downloading (requests) %s -> %s", url, tmp)
    with requests.get(url, stream=True, timeout=60) as r:
        r.raise_for_status()
        total = int(r.headers.get("Content-Length", 0))
        with open(tmp, "wb") as f, tqdm(
            total=total or None, unit="B", unit_scale=True, desc=tmp.name
        ) as bar:
            for chunk in r.iter_content(chunk_size=CHUNK):
                f.write(chunk)
                bar.update(len(chunk))


_BACKENDS = {
    "axel":     _dl_axel,
    "aria2c":   _dl_aria2c,
    "wget":     _dl_wget,
    "requests": _dl_requests,
}


def _download(url: str, dst: Path) -> None:
    dst.parent.mkdir(parents=True, exist_ok=True)
    tmp = dst.with_suffix(dst.suffix + ".part")

    backend = _download_backend()
    log.info("[%s] %s -> %s", backend, url, dst)
    _BACKENDS[backend](url, tmp)

    # decompress on-the-fly if the source URL is .gz but the target ext is not
    if url.endswith(".gz") and not str(dst).endswith(".gz"):
        log.info("decompressing %s -> %s", tmp, dst)
        with gzip.open(tmp, "rb") as fin, open(dst, "wb") as fout:
            shutil.copyfileobj(fin, fout, length=CHUNK)
        tmp.unlink()
    else:
        tmp.rename(dst)


# ── public API ────────────────────────────────────────────────────────────────

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
    log.info("download backend: %s  connections: %s", _download_backend(), CONNECTIONS)
    for tgt in iter_targets(species=species, assembly=assembly, out_root=out):
        log.info("=== download %s / %s ===", tgt.species, tgt.assembly)
        for r in resources:
            try:
                fetch_resource(tgt, r, force=force)
            except Exception as e:
                log.error("[%s/%s] %s FAILED: %s",
                          tgt.species, tgt.assembly, r, e)
