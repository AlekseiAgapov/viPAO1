# Pseudomonas Phage and Mutant Analysis Pipelines

This repository contains analysis pipelines used in the paper **"Host-associated barriers shape experimentally recoverable phage diversity in Pseudomonas aeruginosa"**. It consists of three parts:

1.  `mutation_calling/`: bacterial mutant whole genome sequencing and mutation calling with breseq.
2.  `pooled_phage_analysis/`: pooled environmental and induced phage library processing, reference panel read recruitment, and phage presence/absence calling.
3.  `individual_phage_sequencing/`: long-read sequencing, host-read depletion, assembly, and polishing of individual phage lysates.

## Requirements

The workflows use standard command-line bioinformatics tools, including minimap2, samtools, NanoPlot, chopper, Flye, Medaka, SeqKit, breseq, and Python 3. Tool-specific environment files and requirements are listed in the subdirectory README files.

## Directory Layout

``` text
mutation_calling/
pooled_phage_analysis/
individual_phage_sequencing/
```

Each subdirectory contains a short README and scripts.
