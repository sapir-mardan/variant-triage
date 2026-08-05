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
