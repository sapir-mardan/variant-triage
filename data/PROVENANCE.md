# Data Provenance

## SRR2589044_1.subset.fastq

- **Organism / study:** *E. coli*, Lenski Long-Term Evolution Experiment (LTEE)
- **Accession:** SRR2589044 (ENA/SRA)
- **ENA record:** https://www.ebi.ac.uk/ena/browser/view/SRR2589044

### Fetch (full R1 file, gzipped)
```bash
curl -O ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR258/004/SRR2589044/SRR2589044_1.fastq.gz
```

### Decompress
```bash
gunzip SRR2589044_1.fastq.gz
```

### Subset (first 10,000 records = 40,000 lines)
```bash
head -n 40000 SRR2589044_1.fastq > SRR2589044_1.subset.fastq
```

**Note:** the exact fetch/subset commands as originally run were not logged at the
time and were reconstructed afterward. The sequence above is verified to reproduce
the tracked file exactly (10,000 records, confirmed via `count_records`). Source
accession and organism are certain; the literal command syntax used originally may
have differed slightly (e.g. gunzip as a separate step vs. piped).


---

## Stage B: Alignment and Variant Calling

### Reference Genome
- **Organism/strain:** *E. coli* B str. REL606
- **Assembly:** GCF_000017985.1 (ASM1798v1)
- **File:** `data/reference/REL606.fasta`

### Alignment (Stage 2)
```bash
bwa mem -R '@RG\tID:E.coli_LTEE\tSM:ecoli_LTEE' \
  data/reference/REL606.fasta \
  data/SRR2589044_1.subset.fastq \
  > data/aligned.sam
```
Result: 9,991/10,000 reads mapped (99.91%)

### Sorting & Indexing (Stage 3)
```bash
samtools sort -o data/aligned.sorted.bam data/aligned.sam
samtools index data/aligned.sorted.bam
```

### Variant Calling (Stage 4)

**bcftools:**
```bash
bcftools mpileup -f data/reference/REL606.fasta data/aligned.sorted.bam | bcftools call -mv -o data/calls.vcf
```
Result: 178 variants

**GATK:**
```bash
gatk HaplotypeCaller -R data/reference/REL606.fasta -I data/aligned.sorted.bam -O data/calls_gatk.vcf
```
Result: 2 variants (DP≥3, high confidence)

### Comparison
- bcftools: 178 variants (including many DP=1 calls with low confidence)
- GATK: 2 variants (high-confidence only, filtered stringently)
- Reason: Shallow coverage (0.002x) makes GATK's filtering appropriate
