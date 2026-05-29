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


def _download_backends() -> list[str]:
    """Return all available download backends in priority order."""
    backends = [t for t in ("axel", "aria2c", "wget") if _has(t)]
    if _has("wget"):
        backends.append("wget-insecure")
    backends.append("requests")
    backends.append("requests-insecure")
    return backends


def _download_backend() -> str:
    """Return the highest-priority backend (kept for log line / API)."""
    return _download_backends()[0]


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
    _run(["axel", "-q", "-n", str(CONNECTIONS), "-a", "-o", str(tmp), url])


def _dl_aria2c(url: str, tmp: Path) -> None:
    _run([
        "aria2c",
        "-q",
        "-x", str(CONNECTIONS), "-s", str(CONNECTIONS),
        "-k", "10M",
        "--file-allocation=none",
        "-o", tmp.name, "-d", str(tmp.parent),
        url,
    ])


def _dl_wget(url: str, tmp: Path) -> None:
    _run(["wget", "--quiet", "--tries=3", "--timeout=60", "-O", str(tmp), url])


def _dl_wget_insecure(url: str, tmp: Path) -> None:
    _run(["wget", "--quiet", "--tries=3", "--timeout=60",
          "--no-check-certificate", "-O", str(tmp), url])


def _dl_requests(url: str, tmp: Path) -> None:
    _requests_get(url, tmp, verify=True)


def _dl_requests_insecure(url: str, tmp: Path) -> None:
    _requests_get(url, tmp, verify=False)


def _requests_get(url: str, tmp: Path, *, verify: bool) -> None:
    log.info("downloading (requests, verify=%s) %s -> %s", verify, url, tmp)
    with requests.get(url, stream=True, timeout=60, verify=verify) as r:
        r.raise_for_status()
        total = int(r.headers.get("Content-Length", 0))
        with open(tmp, "wb") as f, tqdm(
            total=total or None, unit="B", unit_scale=True,
            desc=tmp.name, disable=not log.isEnabledFor(logging.DEBUG),
        ) as bar:
            for chunk in r.iter_content(chunk_size=CHUNK):
                f.write(chunk)
                bar.update(len(chunk))


_BACKENDS = {
    "axel":               _dl_axel,
    "aria2c":             _dl_aria2c,
    "wget":               _dl_wget,
    "wget-insecure":      _dl_wget_insecure,
    "requests":           _dl_requests,
    "requests-insecure":  _dl_requests_insecure,
}


def _download(url: str, dst: Path) -> None:
    dst.parent.mkdir(parents=True, exist_ok=True)
    tmp = dst.with_suffix(dst.suffix + ".part")

    backends = _download_backends()
    last_err: Exception | None = None
    for backend in backends:
        # Clean up any partial leftover from a previous backend.
        if tmp.exists():
            try:
                tmp.unlink()
            except OSError:
                pass
        # axel also writes a sidecar .st state file; remove it too.
        st = Path(str(tmp) + ".st")
        if st.exists():
            try:
                st.unlink()
            except OSError:
                pass
        log.info("[%s] %s -> %s", backend, url, dst)
        try:
            _BACKENDS[backend](url, tmp)
            break
        except (subprocess.CalledProcessError, requests.RequestException) as e:
            last_err = e
            log.warning("[%s] failed for %s (%s); trying next backend",
                        backend, url, e)
    else:
        raise RuntimeError(f"all download backends failed for {url}: {last_err}")

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

    # Special path: resource is to be derived by lifting over coordinates from
    # another assembly's upstream file. We fetch the source GFF + chain instead
    # of the canonical raw file; build_<resource>() runs liftOver later.
    lift = spec.get("liftover_from") if isinstance(spec, dict) else None
    if lift:
        return _fetch_liftover_inputs(target, resource, lift, force=force)

    dst = raw_path(target, resource)
    if dst.exists() and not force:
        log.info("[%s/%s] %s already exists: %s",
                 target.species, target.assembly, resource, dst)
        return dst

    local = spec.get("local_path")
    url = spec.get("url")
    extras = spec.get("extra_urls") or []
    if local and Path(local).exists():
        _copy_or_decompress(Path(local), dst)
        return dst
    if url:
        _download(url, dst)
        # Append any extra URLs (e.g. Ensembl ncrna alongside cdna). We download
        # each to a sibling .part file, decompress if needed, and concatenate.
        for extra_url in extras:
            extra_tmp = dst.with_suffix(dst.suffix + ".extra.part")
            try:
                _download(extra_url, extra_tmp)
                # _download decompresses .gz when dst extension is not .gz, so
                # extra_tmp matches dst's compression state.
                with open(dst, "ab") as out, open(extra_tmp, "rb") as fin:
                    shutil.copyfileobj(fin, out, length=CHUNK)
                log.info("[%s/%s] %s: appended %s",
                         target.species, target.assembly, resource, extra_url)
            finally:
                if extra_tmp.exists():
                    extra_tmp.unlink()
        return dst
    log.warning("[%s/%s] %s has neither readable local_path nor url",
                target.species, target.assembly, resource)
    return None


def _fetch_liftover_inputs(
    target: Target,
    resource: str,
    lift: dict,
    *,
    force: bool = False,
) -> Path | None:
    """Download the source-assembly file and chain file for a lifted resource.

    Files are written to:
      raw/<resource>.source.<ext>   (e.g. rnacentral.source.gff3)
      raw/<resource>.chain          (UCSC liftOver chain, decompressed)
    """
    from .config import RESOURCE_EXT  # local import to avoid cycle in tests
    ext = RESOURCE_EXT[resource]
    src_dst = target.raw_dir / f"{resource}.source.{ext}"
    chain_dst = target.raw_dir / f"{resource}.chain"

    url = lift.get("url")
    chain_url = lift.get("chain_url")
    if not url or not chain_url:
        log.warning("[%s/%s] %s liftover_from requires both url and chain_url",
                    target.species, target.assembly, resource)
        return None

    if not src_dst.exists() or force:
        _download(url, src_dst)
    else:
        log.info("[%s/%s] %s source already exists: %s",
                 target.species, target.assembly, resource, src_dst)

    if not chain_dst.exists() or force:
        _download(chain_url, chain_dst)
    else:
        log.info("[%s/%s] %s chain already exists: %s",
                 target.species, target.assembly, resource, chain_dst)
    return src_dst


def download_targets(
    species: list[str] | None = None,
    assembly: list[str] | None = None,
    resources: list[str] | None = None,
    *,
    out: str | None = None,
    force: bool = False,
    include_disabled: bool = False,
) -> None:
    resources = resources or RESOURCE_NAMES
    log.info("download backend: %s  connections: %s", _download_backend(), CONNECTIONS)
    for tgt in iter_targets(species=species, assembly=assembly, out_root=out,
                            include_disabled=include_disabled):
        log.info("=== download %s / %s ===", tgt.species, tgt.assembly)
        for r in resources:
            try:
                fetch_resource(tgt, r, force=force)
            except Exception as e:
                log.error("[%s/%s] %s FAILED: %s",
                          tgt.species, tgt.assembly, r, e)
