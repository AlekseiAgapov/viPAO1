# Mutation Calling Pipeline

This workflow calls mutations in bacterial whole genome sequencing data using breseq. It is intended for trimmed paired-end reads plus orphan reads produced during trimming.

## Input Layout

For each sample, provide four gzipped FASTQ files:

``` text
reads/SAMPLE_1_trimmed.fastq.gz
reads/SAMPLE_2_trimmed.fastq.gz
reads/SAMPLE_U1_trimmed.fastq.gz
reads/SAMPLE_U2_trimmed.fastq.gz
```

References may be FASTA or GenBank. GenBank is preferable when gene-level mutation annotation is required in breseq reports.

## Environment

``` bash
mamba env create -f env/environment.yml
mamba activate breseq-wgs
```

## Run One Sample

``` bash
bash scripts/run_sample_breseq.sh \
  --sample SAMPLE \
  --ref path/to/reference.gbk \
  --outdir breseq_results \
  --threads 6 \
  --reads-dir reads
```

## Main Outputs

``` text
breseq_results/SAMPLE/output/index.html
breseq_results/SAMPLE/output/output.gd
breseq_results/SAMPLE/output/output.vcf
breseq_results/SAMPLE/SAMPLE.breseq.pipeline.log
```

The HTML report should be inspected manually, especially for structural variants, junction calls, and mixed evidence.
