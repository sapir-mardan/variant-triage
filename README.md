# variant-triage

A command-line tool for validating sequencing data and comparing variant callsets.

Variant callers disagree, and standard tools count those disagreements without
explaining them. `variant-triage` is being built to explain them — which calls to
trust, and why. This repository is active work; the sections below describe what
currently runs.

## What it does now

**FASTQ validation and QC** — an installable CLI that reads a FASTQ file and
reports record count, read lengths, GC content, and generates distribution plots.
Malformed and truncated records raise errors rather than passing silently.

```bash
uv run variant-triage data/01_raw/SRR2589044_1.fastq
```
```
{'num_reads': 1107090, 'min_len': 150, 'max_len': 150, 'mean_len': 150.0, 'gc_content': 50.54}
```

**A documented calling pipeline** — alignment through variant calling under two
different callers, with every command and decision recorded in
[`data/PROVENANCE.md`](data/PROVENANCE.md).

## Current result

*E. coli* (LTEE clone SRR2589044) against the REL606 ancestral reference, 72×
coverage, paired-end, duplicates marked, haploid ploidy set explicitly:

| | bcftools | GATK |
|---|---|---|
| Calls | 13 | 15 |

| Shared between callers | 7 |
|---|---|
| **Distinct sites** | **21** |

The similar totals conceal substantial disagreement: only **7 of 21 sites** are
called by both. Those 14 discordant sites are what the project is built to explain.

## Install

Python dependencies via [uv](https://docs.astral.sh/uv/):

```bash
git clone https://github.com/sapir-mardan/variant-triage.git
cd variant-triage
uv sync
```

Bioinformatics tools are compiled binaries and live in a separate conda
environment:

```bash
conda create -n variant-triage-bio -c bioconda -c conda-forge \
  bwa samtools bcftools gatk4 sra-tools ncbi-datasets-cli
conda activate variant-triage-bio
```

## Reproducing the analysis

Sequencing data is not committed — it is public and downloadable. Every command,
input accession and decision is in [`data/PROVENANCE.md`](data/PROVENANCE.md),
which runs from raw reads to the caller comparison.

Callsets are committed under `data/03_results/` so results can be inspected
without re-running the pipeline.

## Tests

```bash
uv run pytest
```

13 tests covering record counting, FASTQ format validation (malformed headers,
mismatched sequence and quality lengths, truncated records) and GC content edge
cases. Ruff and pytest also run on every push via GitHub Actions.

## Known limitations

- Quality lines are not parsed, so mean Phred score, per-position quality and
  adapter contamination cannot be reported — the tool cannot yet answer "should
  these reads be trimmed?"
- No truth set, so caller concordance is measurable but precision and recall are
  not.
- Calls are unfiltered; no QUAL, depth or strand-bias thresholds applied.

## Roadmap

- **Stage C** — the pipeline as a portable Nextflow workflow with containers and
  full provenance
- **Stage D** — biological annotation, calls stored queryably
- **Stage E** — per-call evidence extraction, discordance explained by feature,
  confidence scoring, triaged review list

## Layout

```
src/variant_triage/   package and CLI
tests/                pytest suite
analysis/             notebooks
data/                 01_raw → 02_processed → 03_results, plus reference/
plots/                generated figures
```