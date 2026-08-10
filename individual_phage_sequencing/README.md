# Individual Phage Sequencing Pipeline

This workflow assembles individual phage genomes from Oxford Nanopore reads generated from phage lysate DNA.

For phages recovered during screening that formed plaques but did not form stable lysogens under the conditions tested, phage DNA was extracted directly from lysates using the Norgen Biotek phage DNA extraction kit according to the manufacturer's instructions. DNA from these lysates was sequenced using long-read sequencing by Plasmidsaurus.

## Workflow

Oxford Nanopore reads are quality assessed with NanoPlot, filtered with Chopper, depleted of residual host-derived reads by mapping to the host genome with minimap2, assembled with Flye, and polished with Medaka.

## Environment

``` bash
mamba env create -f env/environment.yml
mamba activate individual-phage-ont
```

Medaka can be installed in the same environment or run from a separate environment by passing `-e ENV_NAME` to `individual_phage_ont_pipeline.sh`.

Run:

``` bash
bash scripts/individual_phage_ont_pipeline.sh \
  -i reads.fastq.gz \
  -H host_reference.fasta \
  -s sample_name \
  -t 6
```

Optional reference comparison:

``` bash
bash scripts/individual_phage_ont_pipeline.sh \
  -i reads.fastq.gz \
  -H host_reference.fasta \
  -s sample_name \
  -r expected_phage_reference.fasta \
  -t 6
```

## Main Steps

``` text
NanoPlot raw-read QC
chopper -q 10 -l 1000
minimap2 -ax map-ont against the host reference
samtools extraction of unmapped reads
assembly of all host-depleted reads by default
optional Filtlong subsetting if a target size is explicitly supplied
Flye --nano-hq assembly
Medaka polishing
SeqKit assembly statistics
read mapping back to the polished assembly
optional comparison to an expected reference
```

## Main Outputs

``` text
sample_name_control_pipeline/02_filter_host/sample_name.filtered.fastq
sample_name_control_pipeline/02_filter_host/sample_name.no_host.fastq
sample_name_control_pipeline/03_flye/assembly.fasta
sample_name_control_pipeline/04_medaka/consensus.fasta
sample_name_control_pipeline/05_reports/
```

The polished `consensus.fasta` is the final assembled phage sequence.
