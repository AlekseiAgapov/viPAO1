#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'USAGE'
Usage:
  individual_phage_ont_pipeline.sh -i reads.fastq.gz -H host.fasta -s sample_name [options]

Required:
  -i FASTQ      Raw ONT reads, optionally gzipped.
  -H FASTA      Host genome FASTA used to remove host reads.
  -s NAME       Sample name. Output goes to NAME_phage_pipeline/.

Optional:
  -r FASTA      Expected phage reference FASTA for control comparison.
  -t INT        Threads. Default: 6.
  -q INT        Minimum read quality for chopper. Default: 10.
  -l INT        Minimum read length for chopper. Default: 1000.
  -b INT        Optional target bases for filtlong subsetting. Default: no subsetting.
  -m MODEL      Medaka model. Default: medaka default auto/model choice.
  -e ENV        Mamba/conda env containing medaka. Default: medaka.
  -h            Show this help.

Assumes core tools are available in the active environment:
  NanoPlot, chopper, minimap2, samtools, flye, seqkit

If -b is supplied, filtlong must also be available.

Medaka is run via:
  mamba run -n <ENV> medaka_consensus
USAGE
}

READS=""
HOST=""
SAMPLE=""
REF=""
THREADS=6
MIN_Q=10
MIN_LEN=1000
TARGET_BASES=""
MEDAKA_ENV="medaka"
MEDAKA_MODEL=""

while getopts ":i:H:s:r:t:q:l:b:m:e:h" opt; do
  case "$opt" in
    i) READS="$OPTARG" ;;
    H) HOST="$OPTARG" ;;
    s) SAMPLE="$OPTARG" ;;
    r) REF="$OPTARG" ;;
    t) THREADS="$OPTARG" ;;
    q) MIN_Q="$OPTARG" ;;
    l) MIN_LEN="$OPTARG" ;;
    b) TARGET_BASES="$OPTARG" ;;
    m) MEDAKA_MODEL="$OPTARG" ;;
    e) MEDAKA_ENV="$OPTARG" ;;
    h) usage; exit 0 ;;
    :) echo "Missing argument for -$OPTARG" >&2; usage; exit 2 ;;
    \?) echo "Unknown option: -$OPTARG" >&2; usage; exit 2 ;;
  esac
done

if [[ -z "$READS" || -z "$HOST" || -z "$SAMPLE" ]]; then
  usage
  exit 2
fi

for f in "$READS" "$HOST"; do
  if [[ ! -s "$f" ]]; then
    echo "Input file not found or empty: $f" >&2
    exit 1
  fi
done

if [[ -n "$REF" && ! -s "$REF" ]]; then
  echo "Reference file not found or empty: $REF" >&2
  exit 1
fi

need_tool() {
  if ! command -v "$1" >/dev/null 2>&1; then
    echo "Required tool not found in PATH: $1" >&2
    exit 1
  fi
}

for tool in NanoPlot chopper minimap2 samtools flye seqkit mamba; do
  need_tool "$tool"
done
if [[ -n "$TARGET_BASES" ]]; then
  need_tool filtlong
fi

OUT="${SAMPLE}_phage_pipeline"
QC_DIR="$OUT/01_qc_raw"
FILTER_DIR="$OUT/02_filter_host"
ASM_DIR="$OUT/03_flye"
POLISH_DIR="$OUT/04_medaka"
REPORT_DIR="$OUT/05_reports"

mkdir -p "$QC_DIR" "$FILTER_DIR" "$ASM_DIR" "$POLISH_DIR" "$REPORT_DIR"

FILTERED="$FILTER_DIR/${SAMPLE}.filtered.fastq"
HOST_BAM="$FILTER_DIR/${SAMPLE}.host_alignments.bam"
NO_HOST="$FILTER_DIR/${SAMPLE}.no_host.fastq"
SUBSET=""
if [[ -n "$TARGET_BASES" ]]; then
  SUBSET="$FILTER_DIR/${SAMPLE}.no_host.${TARGET_BASES}bp.fastq"
fi
ASSEMBLY="$ASM_DIR/assembly.fasta"
POLISHED="$POLISH_DIR/consensus.fasta"

log() {
  printf '[%s] %s\n' "$(date '+%Y-%m-%d %H:%M:%S')" "$*"
}

log "Writing outputs to $OUT"

log "Raw read QC with NanoPlot"
NanoPlot --fastq "$READS" -o "$QC_DIR"

log "Filtering reads with chopper: min_q=$MIN_Q min_len=$MIN_LEN"
chopper -q "$MIN_Q" -l "$MIN_LEN" -i "$READS" > "$FILTERED"

log "Removing host reads with minimap2/samtools"
minimap2 -ax map-ont -t "$THREADS" "$HOST" "$FILTERED" \
  | samtools view -@ "$THREADS" -b -f 4 \
  | samtools fastq -@ "$THREADS" - > "$NO_HOST"

log "Read stats after filtering and host depletion"
seqkit stats "$READS" "$FILTERED" "$NO_HOST" > "$REPORT_DIR/read_stats.tsv"
cat "$REPORT_DIR/read_stats.tsv"

ASSEMBLY_READS="$NO_HOST"
if [[ -n "$TARGET_BASES" ]]; then
  log "Subsetting no-host reads to target_bases=$TARGET_BASES with filtlong"
  filtlong --target_bases "$TARGET_BASES" "$NO_HOST" > "$SUBSET"
  ASSEMBLY_READS="$SUBSET"
  seqkit stats "$NO_HOST" "$SUBSET" > "$REPORT_DIR/subset_stats.tsv"
  cat "$REPORT_DIR/subset_stats.tsv"
else
  log "No filtlong target specified; assembling all no-host reads"
fi

if [[ -n "$REF" ]]; then
  log "Mapping assembly reads to expected phage reference"
  minimap2 -ax map-ont -t "$THREADS" "$REF" "$ASSEMBLY_READS" \
    | samtools sort -@ "$THREADS" -o "$REPORT_DIR/reads_to_reference.bam"
  samtools index "$REPORT_DIR/reads_to_reference.bam"
  samtools coverage "$REPORT_DIR/reads_to_reference.bam" > "$REPORT_DIR/reads_to_reference.coverage.tsv"
  cat "$REPORT_DIR/reads_to_reference.coverage.tsv"
fi

log "Assembling with Flye"
flye --nano-hq "$ASSEMBLY_READS" -o "$ASM_DIR" --threads "$THREADS"

log "Assembly stats"
seqkit stats "$ASSEMBLY" > "$REPORT_DIR/assembly_stats.tsv"
cat "$REPORT_DIR/assembly_stats.tsv"
cp "$ASM_DIR/assembly_info.txt" "$REPORT_DIR/assembly_info.txt"
cat "$REPORT_DIR/assembly_info.txt"

log "Polishing with Medaka from env '$MEDAKA_ENV'"
MEDAKA_ARGS=(-i "$ASSEMBLY_READS" -d "$ASSEMBLY" -o "$POLISH_DIR" -t "$THREADS" -f)
if [[ -n "$MEDAKA_MODEL" ]]; then
  MEDAKA_ARGS+=(-m "$MEDAKA_MODEL")
fi
mamba run -n "$MEDAKA_ENV" medaka_consensus "${MEDAKA_ARGS[@]}"

log "Polished assembly stats"
seqkit stats "$POLISHED" > "$REPORT_DIR/polished_stats.tsv"
cat "$REPORT_DIR/polished_stats.tsv"
seqkit fx2tab -n -i -l -g "$POLISHED" > "$REPORT_DIR/polished_lengths_gc.tsv"
cat "$REPORT_DIR/polished_lengths_gc.tsv"

log "Mapping assembly reads back to polished assembly"
minimap2 -ax map-ont -t "$THREADS" "$POLISHED" "$ASSEMBLY_READS" \
  | samtools sort -@ "$THREADS" -o "$REPORT_DIR/reads_to_polished.bam"
samtools index "$REPORT_DIR/reads_to_polished.bam"
samtools coverage "$REPORT_DIR/reads_to_polished.bam" > "$REPORT_DIR/reads_to_polished.coverage.tsv"
cat "$REPORT_DIR/reads_to_polished.coverage.tsv"

if [[ -n "$REF" ]]; then
  log "Comparing polished assembly to expected reference"
  minimap2 -cx asm5 --cs "$REF" "$POLISHED" > "$REPORT_DIR/polished_vs_reference.paf"
  awk 'BEGIN { OFS="\t"; print "query","query_len","query_start","query_end","strand","reference","reference_len","reference_start","reference_end","matches","aln_len","identity","mapq" }
       { print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$10/$11,$12 }' \
    "$REPORT_DIR/polished_vs_reference.paf" > "$REPORT_DIR/polished_vs_reference.summary.tsv"
  cat "$REPORT_DIR/polished_vs_reference.summary.tsv"

  if command -v dnadiff >/dev/null 2>&1; then
    log "Running dnadiff against expected reference"
    dnadiff -p "$REPORT_DIR/reference_vs_polished" "$REF" "$POLISHED"
  else
    log "dnadiff not found; skipping MUMmer comparison"
  fi
fi

N_CONTIGS=$(seqkit seq -n "$POLISHED" | wc -l | tr -d ' ')
log "Finished. Polished contigs: $N_CONTIGS"

if [[ "$N_CONTIGS" -ne 1 ]]; then
  echo
  echo "NOTE: Individual phage assemblies often should produce one main contig."
  echo "This run produced $N_CONTIGS polished contigs. Inspect:"
  echo "  $REPORT_DIR/assembly_info.txt"
  echo "  $REPORT_DIR/reads_to_polished.coverage.tsv"
  echo "  $REPORT_DIR/polished_vs_reference.summary.tsv, if a reference was supplied"
fi
