# Pooled Environmental and Induced Phage Analysis Pipeline

This workflow processes pooled long-read libraries and calls phage presence against a curated reference panel. It applies to both environmental phage pools and induced prophage pools.

## Environment

```bash
mamba env create -f env/environment.yml
mamba activate phage-pools
```

Medaka can be installed in the same environment or, if that does not work, it can be run from a separate environment by passing `-e ENV_NAME` to `pooled_phage_ont_pipeline.sh`.

## 1. Read QC, Filtering, Host Depletion, Assembly, and Polishing

Use the pooled ONT pipeline:

```bash
bash scripts/pooled_phage_ont_pipeline.sh \
  -i reads.fastq.gz \
  -H host_reference.fasta \
  -s sample_name \
  -t 8
```

The pipeline performs:

```text
NanoPlot raw-read QC
chopper -q 10 -l 1000
minimap2 -ax map-ont against the host reference
samtools extraction of unmapped reads
assembly of all host-depleted reads with Flye --meta --nano-hq by default
Medaka polishing
candidate contig export
```

## 2. Environmental Contig Curation

Polished environmental assemblies were analysed with geNomAD v1.12.0 and CheckV v1.0.3 using the `end-to-end` workflow in both cases. Contigs were manually curated to retain medium-to-high confidence phage-derived contigs, defined as contigs with a geNomAD viral score greater than 0.9 and CheckV-estimated phage genome completeness of at least 50%.

Typical commands:

```bash
genomad end-to-end --cleanup --splits 8 candidate_contigs.fasta genomad_out genomad_db

checkv end_to_end candidate_contigs.fasta checkv_out -d checkv-db-v1.5 -t 8
```

Retained contigs were compared all-versus-all with minimap2 v2.30, and near-duplicate contigs were clustered when pairwise identity was at least 90% and alignment coverage of the shorter contig was at least 50%.

Typical all-versus-all command:

```bash
minimap2 -x asm20 -c --eqx -N 50 retained_contigs.fasta retained_contigs.fasta \
  > retained_contigs.all_vs_all.asm20.paf
```

## 3. Induced Prophage Candidate Discovery

For induced libraries, host-depleted reads were mapped back to the closed bacterial genomes corresponding to the induction host background. Covered regions with elevated read depth were treated as candidate induced prophage loci. These candidates should be inspected manually using coverage plots and genomic context, then representative prophage sequences should be extracted from the bacterial reference genome.

Example mapping command:

```bash
minimap2 -a -x map-ont -t 8 host_genome.fasta sample.no_host.fastq \
  | samtools sort -@ 8 -o sample.to_host_genome.bam

samtools index sample.to_host_genome.bam
samtools depth -Q 20 -G SUPPLEMENTARY sample.to_host_genome.bam \
  > sample.to_host_genome.mapq20.depth.tsv
```

## 4. Read Recruitment to a Curated Panel

Prepare a sample table:

```text
sample	read_fastq
library_1	path/to/library_1.no_host.fastq
library_2	path/to/library_2.no_host.fastq
```

Map host-depleted reads to the curated panel:

```bash
bash scripts/map_no_host_reads_to_panel.sh \
  --ref curated_panel.fasta \
  --samples sample_inputs.tsv \
  --outdir panel_mapping \
  --threads 8
```

The script suppresses secondary alignments so each read contributes one primary placement. It creates MAPQ >= 20 primary alignments and calls presence/absence from covered length and read support.

The output folder stores the reference panel and sample table in a standard layout:

```text
panel_mapping/ref/panel.fasta
panel_mapping/sample_inputs.tsv
panel_mapping/mapping/
panel_mapping/results/
```

## 5. Presence/Absence Calling

A phage genome is considered present if:

```text
primary MAPQ >=20 alignments cover at least max(50% of reference length, 10000 bp)
and at least 5 MAPQ >=20 primary reads map to the reference
```

Final output:

```text
panel_mapping/results/present_max50pct_10kb_mq20_primary_matrix.tsv
```

## 6. Presence/Absence Heatmap

The binary presence/absence matrix can be plotted as a PDF heatmap:

```bash
cd panel_mapping/results

python /path/to/scripts/plot_phage_presence_absence.py
```
