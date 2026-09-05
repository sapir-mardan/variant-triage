# Data Provenance

*All commands run from the project root.*

**Contents:** [1. Inputs](#1-inputs) · [2. Pipeline](#2-pipeline) · [3. Decisions](#3-decisions) · [4. Limitations](#4-limitations) · [5. Environment](#5-environment) · [6. Previous run](#6-previous-run-superseded)

---

## 1. Inputs

### 1.1 Reads

- **Organism / study:** *E. coli*, Lenski Long-Term Evolution Experiment (LTEE)
- **Accession:** SRR2589044 (ENA/SRA)
- **ENA record:** https://www.ebi.ac.uk/ena/browser/view/SRR2589044
- **SRA record:** https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&acc=SRR2589044&display=metadata

```bash
# requires: conda install -c bioconda -c conda-forge sra-tools
fasterq-dump SRR2589044 --split-files -O data/01_raw/
```
```
spots read    : 1,107,090
reads read    : 2,214,180
reads written : 2,214,180
```

### 1.2 Reference genome

- **Organism/strain:** *E. coli* B str. REL606
- **Assembly:** GCF_000017985.1 (ASM1798v1)
- **Source:** https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000017985.1/
- **Why REL606:** the LTEE ancestral strain, so differences from it are mutations
  that arose during the experiment. A generic K-12 reference would add tens of
  thousands of pre-existing strain differences.

```bash
# requires: conda install -c conda-forge ncbi-datasets-cli
datasets download genome accession GCF_000017985.1 --include genome
unzip ncbi_dataset.zip
mv ncbi_dataset/data/GCF_000017985.1/GCF_000017985.1_ASM1798v1_genomic.fna data/reference/REL606.fasta
rm -rf ncbi_dataset ncbi_dataset.zip README.md md5sum.txt   # remove unzipped extras
```

### 1.3 Reference indexes

```bash
bwa index data/reference/REL606.fasta        # alignment index (.amb .ann .bwt .pac .sa)
samtools faidx data/reference/REL606.fasta   # position index (.fai)
gatk CreateSequenceDictionary -R data/reference/REL606.fasta   # GATK needs .dict
```

Genome size: **4,629,812 bp** — column 2 of `REL606.fasta.fai`.

### 1.4 Coverage

`coverage = reads × read length ÷ genome size`

```bash
head -2 data/01_raw/SRR2589044_1.fastq | tail -1 | wc -c   # → 151 (150 bp + newline)
python3 -c "print(2214180 * 150 / 4629812)"                # → 71.7
```

**Result: ~72×** (first Stage B run: ~0.3×)

---

## 2. Pipeline

### 2.1 QC

```bash
uv run variant-triage data/01_raw/SRR2589044_1.fastq
uv run variant-triage data/01_raw/SRR2589044_2.fastq
```
```
_1: {'num_reads': 1107090, 'min_len': 150, 'max_len': 150, 'mean_len': 150.0, 'gc_content': 50.54}
_2: {'num_reads': 1107090, 'min_len': 150, 'max_len': 150, 'mean_len': 150.0, 'gc_content': 50.646}
```

Read counts match `fasterq-dump` and each other — pairing intact, nothing truncated.
GC ~50.5% is consistent with *E. coli* (~50.8%). Uniform 150 bp confirms reads are
raw and untrimmed.

### 2.2 Alignment (paired-end)

```bash
bwa mem -R '@RG\tID:SRR2589044\tSM:ecoli_LTEE' \
  data/reference/REL606.fasta \
  data/01_raw/SRR2589044_1.fastq \
  data/01_raw/SRR2589044_2.fastq \
  > data/02_processed/aligned_pe.sam
```

`-R` sets the read group; `SM` (sample name) is required by GATK, which fails
without it. Runtime 62 s.

### 2.3 Sort and index

```bash
samtools sort -o data/02_processed/aligned_pe.sorted.bam data/02_processed/aligned_pe.sam
samtools index data/02_processed/aligned_pe.sorted.bam
samtools flagstat data/02_processed/aligned_pe.sorted.bam
```
```
2241568 + 0 in total
2214180 + 0 primary
  27388 + 0 supplementary
      0 + 0 duplicates
2201840 + 0 primary mapped (99.44%)
2110968 + 0 properly paired (95.34%)
  12022 + 0 singletons (0.54%)
```

99.44% mapped confirms the reference choice. 95.34% properly paired indicates a
normal library. Duplicates show as 0 because nothing has marked them yet.

This sorted BAM is a checkpoint for `flagstat` only — duplicate marking below
starts again from the SAM.

### 2.4 Mark PCR duplicates

PCR amplification during library prep produces multiple reads from the same original
molecule. They are not independent evidence: an error introduced in an early
amplification cycle appears in every copy and can mimic a real variant.

Two sorts are needed because the tools have different requirements: `fixmate` needs
the two mates adjacent (name order), `markdup` needs reads at the same locus grouped
(position order).

```bash
samtools sort -n -o data/02_processed/aligned_pe.namesort.bam data/02_processed/aligned_pe.sam
samtools fixmate -m data/02_processed/aligned_pe.namesort.bam data/02_processed/aligned_pe.fixmate.bam
samtools sort -o data/02_processed/aligned_pe.possort.bam data/02_processed/aligned_pe.fixmate.bam
samtools markdup data/02_processed/aligned_pe.possort.bam data/02_processed/aligned_pe.markdup.bam
samtools index data/02_processed/aligned_pe.markdup.bam
samtools flagstat data/02_processed/aligned_pe.markdup.bam
```
```
29597 + 0 duplicates
2201840 + 0 primary mapped (99.44%)
2110968 + 0 properly paired (95.34%)
```

**1.34% duplicates** — low, indicating a library that was not over-amplified (under
5% is healthy; 20%+ suggests too few starting molecules). Mapping and pairing rates
are unchanged: marking sets a flag rather than removing reads, so callers skip them
but the data stays.

`-m` on `fixmate` adds mate score tags, used by `markdup` to choose which read in a
duplicate set to keep.

`gatk MarkDuplicates` does the same job in a single command and is the more common
choice in GATK-based pipelines. The samtools route was used here.

Intermediates removed once `markdup.bam` was verified:

```bash
rm data/02_processed/aligned_pe.namesort.bam \
   data/02_processed/aligned_pe.fixmate.bam \
   data/02_processed/aligned_pe.possort.bam
```

`data/02_processed/aligned_pe.markdup.bam` is the input to both callers.

### 2.5 Variant calling

Both callers default to diploid. *E. coli* is haploid, so ploidy is set explicitly:
a diploid model would allow heterozygous genotypes that cannot exist in this
organism. The first Stage B run left the default in place.

**bcftools** — pileup-based, Bayesian:

```bash
bcftools mpileup -f data/reference/REL606.fasta data/02_processed/aligned_pe.markdup.bam \
  | bcftools call -mv --ploidy 1 -o data/03_results/calls_pe.vcf
```
```
[mpileup] maximum number of reads per input file set to -d 250
```

`-m` multiallelic model, `-v` variants only. The `-d 250` cap is bcftools' default
maximum reads per position — not reached at 72×, but relevant on high-depth
targeted data.

**GATK** — local haplotype assembly:

```bash
gatk HaplotypeCaller \
  -R data/reference/REL606.fasta \
  -I data/02_processed/aligned_pe.markdup.bam \
  -O data/03_results/calls_gatk_pe.vcf \
  --sample-ploidy 1
```

Slower than bcftools: it reassembles reads into candidate haplotypes rather than
tallying bases position by position.

```bash
grep -vc "^#" data/03_results/calls_pe.vcf        # → 13
grep -vc "^#" data/03_results/calls_gatk_pe.vcf   # → 15
```

| | First run (~0.3×) | This run (~72×) |
|---|---|---|
| bcftools | 178 | **13** |
| GATK | 2 | **15** |

The 178 → 13 drop reflects removal of coverage noise: at ~0.3× a single read
carrying a sequencing error was often the only evidence at a position, with nothing
to contradict it.

### 2.6 Callset comparison

Similar counts do not mean the same calls. `bcftools isec` splits two callsets into
private and shared records.

```bash
bgzip -k data/03_results/calls_pe.vcf && bcftools index data/03_results/calls_pe.vcf.gz
bgzip -k data/03_results/calls_gatk_pe.vcf && bcftools index data/03_results/calls_gatk_pe.vcf.gz

bcftools isec -p data/03_results/isec_pe \
  data/03_results/calls_pe.vcf.gz \
  data/03_results/calls_gatk_pe.vcf.gz

for f in data/03_results/isec_pe/000*.vcf; do echo -n "$f: "; grep -vc "^#" $f; done
```
```
0000.vcf: 6     bcftools only
0001.vcf: 8     GATK only
0002.vcf: 7     shared (file 1 representation)
0003.vcf: 7     shared (file 2 representation)
```

`bgzip` is block compression, which unlike plain gzip allows random access;
indexing enables position lookups rather than linear scanning. Both are
prerequisites for `isec`.

| Category | Count |
|---|---|
| bcftools only | 6 |
| GATK only | 8 |
| Shared | 7 |
| **Distinct sites** | **21** |

Check: 6 + 7 = 13 bcftools, 8 + 7 = 15 GATK.

**Only 7 of 21 sites are called by both.** The similar totals (13 vs 15) conceal
substantial disagreement — two-thirds of calls are made by one caller and not the
other. The 14 discordant sites are the material this project exists to explain.

---

## 3. Decisions

| Decision | Reason |
|---|---|
| **Paired-end** (first run was single-end, R1 only) | Mates constrain each other's placement, improving mapping in repetitive regions; also yields insert-size and proper-pair metrics unavailable single-end. |
| **Duplicates marked** | At 72×, PCR duplicates inflate apparent support for a call. At ~0.3× this was irrelevant — almost no position had more than one read. |
| **`samtools markdup` over `gatk MarkDuplicates`** | Equivalent outcome; samtools is lighter. GATK's is the more common choice in GATK-based pipelines and produces a metrics file. |
| **`--ploidy 1` on both callers** | *E. coli* is haploid. Both callers default to diploid; the first run left that default unset and unrecorded. |
| **No trimming** | Reads are uniform 150 bp and 99.44% mapped, so trimming was not indicated. Note this rests on length and mapping rate, not on quality scores — see Limitations. |
| **Two callers retained** | bcftools (pileup, Bayesian) and GATK (local haplotype assembly) use different algorithms. Comparing them is the purpose of the project. |
| **Reference REL606, not K-12** | REL606 is the LTEE ancestral strain. A K-12 reference would add tens of thousands of pre-existing strain differences, burying the experimental mutations. |

---

## 4. Limitations

- **Three variables changed at once** — depth, paired-end, and ploidy. Depth is
  almost certainly dominant, but their effects cannot be separated from a single
  run. Isolating depth would require downsampling with the other two held constant.
- **No truth set.** Concordance between callers is measurable; precision and recall
  are not. That requires a benchmarked sample such as GIAB.
- **No quality-based QC.** The Stage A tool reads sequence lines but not FASTQ
  quality lines, so mean Phred, per-position quality and adapter contamination are
  unassessed. The "no trimming" decision rests on read-length uniformity and mapping
  rate rather than on quality scores.
- **Calls are unfiltered.** No QUAL, depth or strand-bias thresholds have been
  applied to either callset.

---

## 5. Environment

Bio tools are installed in an isolated conda env (`variant-triage-bio`), separate
from the project's uv environment, which holds only Python dependencies
(`pyproject.toml` / `uv.lock`).

```
bwa                 0.7.19-r1273
samtools            1.24
bcftools            1.24
gatk                4.6.2.0
fasterq-dump        3.4.1
ncbi-datasets-cli
```

```bash
conda create -n variant-triage-bio -c bioconda -c conda-forge \
  bwa samtools bcftools gatk4 sra-tools ncbi-datasets-cli
```

---

## 6. Previous run (superseded)

An initial Stage B run used a 10,000-read single-end subset (~0.3× coverage). It
produced 178 bcftools calls and 2 GATK calls — a coverage artefact rather than an
algorithmic difference: at that depth most positions had one read or none, so
sequencing errors and real variants were indistinguishable. Ploidy was left at the
diploid default. Those callsets have been removed; this run supersedes them.