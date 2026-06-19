# Output Layout

## Working tree (during a build)

During a run, each assembly has a `raw/` (downloads) and `build/`
(browser-loadable) working tree. This is what you get with `refbox pull --no-flat`
or `refbox download`.

```text
{REFBOX_OUT}/
  {Species}/
    {Assembly}/
      raw/                          # original downloads / copies
        genome.fa
        transcriptome.fa             # may be derived from genome + GTF
        annotation_gtf.gtf
        annotation_gff3.gff3
        repeats_rmsk.tsv
        repeats_fa.fa
        rnacentral.gff3              # may be lifted from another assembly
        ccre.bed
        cytoband.tsv                 # UCSC cytoBand[Ideo] text, or a bigBed (hs1)
      build/                        # browser-loadable
        genome.fa.gz                 + .fai + .gzi
        chrom.sizes
        transcriptome.fa.gz          + .fai
        transcriptome.derived.fa.gz  + .fai   # gffread-extracted, when GTF available
        annotation.sorted.gtf.gz     + .tbi
        annotation.sorted.gff3.gz    + .tbi
        repeats.sorted.bed.gz        + .tbi
        repeats.sorted.gtf.gz        + .tbi
        rnacentral.sorted.gff3.gz    + .tbi
        ccre.sorted.bed.gz           + .tbi
        cytoband.sorted.bed.gz       + .tbi
        cytoband.bb                          # bigBed (bed4+1 with gieStain)
```

## Published layout (default)

`pull` finishes with the **publish** step (skip it with `--no-flat`), which
flattens each `build/<name>` to `<Assembly>.<name>` directly under the assembly
directory and removes `build/` / `raw/`:

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

## Environment variables

| Variable | Meaning |
|---|---|
| `REFBOX_OUT` | default output root for `{Species}/{Assembly}/{raw,build}/` (else the current directory) |
| `REFBOX_CONFIG` | path to a custom `species.yaml` (overrides the bundled registry) |

## Status reporter

```bash
python -m refbox.report --out /path/to/reference > report.md
```

Walks the output tree and emits a Markdown report listing, per assembly, the
status (✓ done / ⚠ missing index / ✗ missing) and size of each artifact.
