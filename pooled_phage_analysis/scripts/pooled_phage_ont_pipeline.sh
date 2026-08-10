#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'USAGE'
Usage:
  pooled_phage_ont_pipeline.sh -i reads.fastq.gz -H host.fasta -s sample_name [options]

Required:
  -i FASTQ      Raw ONT reads, optionally gzipped.
  -H FASTA      Host genome FASTA used to remove host reads.
  -s NAME       Sample name. Output goes to NAME_pooled_pipeline/.

Optional:
  -t INT        Threads. Default: 6.
  -q INT        Minimum read quality for chopper. Default: 10.
  -l INT        Minimum read length for chopper. Default: 1000.
  -c INT        Minimum contig length for candidate FASTA. Default: 10000.
  -F MODE       Flye read mode: nano-hq or nano-raw. Default: nano-hq.
  -m MODEL      Medaka model. Default: medaka default auto/model choice.
  -e ENV        Mamba/conda env containing medaka. Default: medaka.
  -P            Skip Medaka polishing.
  -h            Show this help.

Assumes core tools are available in the active environment:
  NanoPlot, chopper, minimap2, samtools, flye, seqkit

Medaka is run via:
  mamba run -n <ENV> medaka_consensus

USAGE
}

READS=""
HOST=""
SAMPLE=""
THREADS=6
MIN_Q=10
MIN_LEN=1000
MIN_CONTIG_LEN=10000
FLYE_MODE="nano-hq"
MEDAKA_ENV="medaka"
MEDAKA_MODEL=""
SKIP_POLISH=0

while getopts ":i:H:s:t:q:l:c:F:m:e:Ph" opt; do
  case "$opt" in
    i) READS="$OPTARG" ;;
    H) HOST="$OPTARG" ;;
    s) SAMPLE="$OPTARG" ;;
    t) THREADS="$OPTARG" ;;
    q) MIN_Q="$OPTARG" ;;
    l) MIN_LEN="$OPTARG" ;;
    c) MIN_CONTIG_LEN="$OPTARG" ;;
    F) FLYE_MODE="$OPTARG" ;;
    m) MEDAKA_MODEL="$OPTARG" ;;
    e) MEDAKA_ENV="$OPTARG" ;;
    P) SKIP_POLISH=1 ;;
    h) usage; exit 0 ;;
    :) echo "Missing argument for -$OPTARG" >&2; usage; exit 2 ;;
    \?) echo "Unknown option: -$OPTARG" >&2; usage; exit 2 ;;
  esac
done

if [[ -z "$READS" || -z "$HOST" || -z "$SAMPLE" ]]; then
  usage
  exit 2
fi

case "$FLYE_MODE" in
  nano-hq|nano-raw) ;;
  *)
    echo "Invalid Flye mode: $FLYE_MODE. Use nano-hq or nano-raw." >&2
    exit 2
    ;;
esac

for f in "$READS" "$HOST"; do
  if [[ ! -s "$f" ]]; then
    echo "Input file not found or empty: $f" >&2
    exit 1
  fi
done

need_tool() {
  if ! command -v "$1" >/dev/null 2>&1; then
    echo "Required tool not found in PATH: $1" >&2
    exit 1
  fi
}

for tool in NanoPlot chopper minimap2 samtools flye seqkit; do
  need_tool "$tool"
done

if [[ "$SKIP_POLISH" -eq 0 ]]; then
  need_tool mamba
fi

OUT="${SAMPLE}_pooled_pipeline"
QC_DIR="$OUT/01_qc_raw"
FILTER_DIR="$OUT/02_filter_host"
ASM_DIR="$OUT/03_flye"
POLISH_DIR="$OUT/04_medaka"
REPORT_DIR="$OUT/05_reports"
FINAL_DIR="$OUT/06_final"

mkdir -p "$QC_DIR" "$FILTER_DIR" "$ASM_DIR" "$POLISH_DIR" "$REPORT_DIR" "$FINAL_DIR"

FILTERED="$FILTER_DIR/${SAMPLE}.filtered.fastq"
NO_HOST="$FILTER_DIR/${SAMPLE}.no_host.fastq"
ASSEMBLY="$ASM_DIR/assembly.fasta"
POLISHED="$POLISH_DIR/consensus.fasta"
FINAL_FASTA="$FINAL_DIR/${SAMPLE}.final.fasta"
CANDIDATES="$FINAL_DIR/${SAMPLE}.candidate_phage_contigs.min${MIN_CONTIG_LEN}.fasta"

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
log "Assembling all no-host reads directly"

log "Assembling pooled phage reads with Flye --meta --$FLYE_MODE"
if [[ "$FLYE_MODE" == "nano-hq" ]]; then
  flye --nano-hq "$ASSEMBLY_READS" --meta -o "$ASM_DIR" --threads "$THREADS"
else
  flye --nano-raw "$ASSEMBLY_READS" --meta -o "$ASM_DIR" --threads "$THREADS"
fi

log "Assembly stats"
seqkit stats "$ASSEMBLY" > "$REPORT_DIR/assembly_stats.tsv"
cat "$REPORT_DIR/assembly_stats.tsv"
cp "$ASM_DIR/assembly_info.txt" "$REPORT_DIR/assembly_info.txt"
cat "$REPORT_DIR/assembly_info.txt"

if [[ "$SKIP_POLISH" -eq 0 ]]; then
  log "Polishing with Medaka from env '$MEDAKA_ENV'"
  MEDAKA_ARGS=(-i "$ASSEMBLY_READS" -d "$ASSEMBLY" -o "$POLISH_DIR" -t "$THREADS" -f)
  if [[ -n "$MEDAKA_MODEL" ]]; then
    MEDAKA_ARGS+=(-m "$MEDAKA_MODEL")
  fi
  mamba run -n "$MEDAKA_ENV" medaka_consensus "${MEDAKA_ARGS[@]}"
  cp "$POLISHED" "$FINAL_FASTA"
else
  log "Skipping Medaka polishing"
  cp "$ASSEMBLY" "$FINAL_FASTA"
fi

log "Final assembly stats"
seqkit stats "$FINAL_FASTA" > "$REPORT_DIR/final_stats.tsv"
cat "$REPORT_DIR/final_stats.tsv"
seqkit fx2tab -n -i -l -g "$FINAL_FASTA" > "$REPORT_DIR/final_lengths_gc.tsv"
cat "$REPORT_DIR/final_lengths_gc.tsv"

log "Extracting candidate phage-sized contigs with min length $MIN_CONTIG_LEN"
seqkit seq -m "$MIN_CONTIG_LEN" "$FINAL_FASTA" > "$CANDIDATES"
seqkit stats "$CANDIDATES" > "$REPORT_DIR/candidate_contig_stats.tsv"
cat "$REPORT_DIR/candidate_contig_stats.tsv"

log "Mapping assembly reads back to final assembly"
minimap2 -ax map-ont -t "$THREADS" "$FINAL_FASTA" "$ASSEMBLY_READS" \
  | samtools sort -@ "$THREADS" -o "$REPORT_DIR/reads_to_final.bam"
samtools index "$REPORT_DIR/reads_to_final.bam"
samtools coverage "$REPORT_DIR/reads_to_final.bam" > "$REPORT_DIR/reads_to_final.coverage.tsv"
cat "$REPORT_DIR/reads_to_final.coverage.tsv"

log "Creating all-vs-all contig similarity table"
minimap2 -cx asm5 --secondary=yes -N 50 "$FINAL_FASTA" "$FINAL_FASTA" > "$REPORT_DIR/final_self.paf"
awk 'BEGIN { OFS="\t"; print "query","query_len","query_start","query_end","strand","target","target_len","target_start","target_end","matches","aln_len","identity","mapq" }
     $1 != $6 { print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$10/$11,$12 }' \
  "$REPORT_DIR/final_self.paf" > "$REPORT_DIR/final_self.nonself.summary.tsv"
cat "$REPORT_DIR/final_self.nonself.summary.tsv"

N_FINAL=$(seqkit seq -n "$FINAL_FASTA" | wc -l | tr -d ' ')
N_CANDIDATES=$(seqkit seq -n "$CANDIDATES" | wc -l | tr -d ' ')

log "Finished. Final contigs: $N_FINAL. Candidate contigs >= ${MIN_CONTIG_LEN} bp: $N_CANDIDATES"
echo
echo "Inspect these first:"
echo "  $REPORT_DIR/assembly_info.txt"
echo "  $REPORT_DIR/reads_to_final.coverage.tsv"
echo "  $REPORT_DIR/final_lengths_gc.tsv"
echo "  $REPORT_DIR/final_self.nonself.summary.tsv"
echo "  $CANDIDATES"
