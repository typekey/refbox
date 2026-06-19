# Pipeline Commands

These four subcommands operate on the **`species.yaml` registry** — they know
about species, assemblies and the 11 canonical resources. For building arbitrary
user files, see the [Build Command](build.md).

```text
refbox download   # only fetch raw files configured in species.yaml
refbox pull       # full pipeline: download (if missing) + build + test + publish
refbox publish    # flatten an existing build/ tree to <Assembly>.<name>
refbox test       # validate build/ outputs
```

## Common filters

Every pipeline command accepts the same selection filters:

| Flag | Meaning |
|---|---|
| `--species [NAME ...]` | species filter (optional; inferred from `--assembly`) |
| `--assembly [NAME ...]` | assembly filter (e.g. `GRCh38`); omit to run **every** assembly matching `--species` (or all if neither is given) |
| `--resource [NAME ...]` | subset of resources (see list below) |
| `--out DIR` | output root (default: `$REFBOX_OUT` or the current directory) |
| `--include-disabled` | also process assemblies marked `enabled: false` in `species.yaml` |
| `-v`, `--verbose` | DEBUG logging (global flag, place before the subcommand) |

**Resource names** (`--resource`): `genome` · `transcriptome` ·
`annotation_gtf` · `annotation_gff3` · `repeats_rmsk` · `repeats_bed` ·
`repeats_gtf` · `repeats_fa` · `rnacentral` · `ccre` · `cytoband`.

---

## `refbox pull`

The registry-driven pipeline: **download → build → test → publish**. This is
what most users want.

By default the per-assembly `build/` outputs are **flattened** to the published
layout `<out>/<Species>/<Assembly>/<Assembly>.<name>` and the `build/` + `raw/`
working directories are removed. Pass `--no-flat` to keep them.

### Extra flags

| Flag | Meaning |
|---|---|
| `--force` | rebuild even when outputs already exist |
| `--no-download` | skip the auto-download phase (use existing `raw/`) |
| `--no-test` | skip the post-build validation step |
| `--no-flat` | keep the `build/` + `raw/` layout (skip publish/flatten) |

### Examples

```bash
# Fetch + build + validate + publish Human GRCh38 (bundled species.yaml)
refbox pull --species Homo_sapiens --assembly GRCh38

# --species inferred from --assembly via the registry
refbox pull --assembly GRCm38

# Every assembly in the registry, including enabled: false, selected resources
refbox pull --include-disabled \
            --resource genome transcriptome annotation_gtf annotation_gff3 \
                       repeats_rmsk rnacentral cytoband

# Only one resource
refbox pull --assembly GRCh38 --resource ccre

# Keep the build/ + raw/ working tree instead of the flat layout
refbox pull --assembly GRCh38 --no-flat

# Rebuild from scratch
refbox pull --assembly GRCh38 --force
```

The exit code is `1` if validation found problems (and `--no-test` was not given),
otherwise `0`.

---

## `refbox download`

Fetch the raw files configured in `species.yaml` and **stop** — no build.

| Flag | Meaning |
|---|---|
| common filters | see above |
| `--force` | re-download even if the raw file already exists |

The download backend tries, in order, `axel → aria2c → wget → wget --no-check-certificate
→ requests → requests (verify=False)`, so a single broken TLS host does not abort
a resource.

```bash
refbox download --assembly GRCh38                      # raw files only
refbox download --assembly GRCh38 --resource genome    # one resource
refbox download --include-disabled                     # everything
```

---

## `refbox publish`

Flatten an existing `build/` tree to `<Assembly>.<name>` and drop `raw/`. `pull`
runs this automatically unless `--no-flat` was given — run it standalone to
flatten after a `--no-flat` build.

| Flag | Meaning |
|---|---|
| common filters | see above |
| `--keep-build` | do not remove the (now-empty) `build/` directory |
| `--keep-raw` | do not remove the `raw/` directory |

```bash
refbox publish --assembly GRCh38                  # one assembly
refbox publish --include-disabled                 # everything in the registry
refbox publish --assembly GRCh38 --keep-raw       # keep raw/ downloads
```

---

## `refbox test`

Re-run the validators against existing `build/` outputs. Checks that each
expected `.gz` exists with its index (`.fai`/`.gzi`/`.tbi`), that a sample tabix
region query returns rows, and that `samtools faidx` returns sequence.

```bash
refbox test --assembly GRCh38           # validate one assembly
refbox test --include-disabled          # validate everything in the registry
```

Exit code is `1` if any check failed, otherwise `0`.
