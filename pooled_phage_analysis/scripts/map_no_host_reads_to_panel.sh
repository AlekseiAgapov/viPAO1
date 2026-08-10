#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'USAGE'
Usage:
  map_no_host_reads_to_panel.sh \
    --ref panel.fasta \
    --samples sample_inputs.tsv \
    --outdir mapping_analysis \
    [--threads 6]

Required:
  --ref       Multifasta reference panel used for read recruitment.
  --samples   Tab-separated file with columns: sample read_fastq
  --outdir    Output directory for mapping and result tables.

Optional:
  --threads   Number of minimap2 threads. Default: 6.

The sample table should contain one row per library:
  sample    read_fastq
  library1  path/to/library1.no_host.fastq

Presence is called downstream as:
  MAPQ >=20 primary alignments cover at least max(50% of reference length, 10000 bp)
  and at least 5 MAPQ >=20 primary reads map to the reference.

Final output:
  OUTDIR/results/present_max50pct_10kb_mq20_primary_matrix.tsv
USAGE
}

die() {
  echo "ERROR: $*" >&2
  exit 1
}

REF=""
SAMPLES=""
OUTDIR=""
THREADS=6
MINIMAP2="${MINIMAP2:-minimap2}"
SAMTOOLS="${SAMTOOLS:-samtools}"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --ref) REF="${2:-}"; shift 2 ;;
    --samples) SAMPLES="${2:-}"; shift 2 ;;
    --outdir) OUTDIR="${2:-}"; shift 2 ;;
    --threads) THREADS="${2:-}"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) die "Unknown option: $1" ;;
  esac
done

[[ -n "$REF" ]] || die "Missing --ref"
[[ -n "$SAMPLES" ]] || die "Missing --samples"
[[ -n "$OUTDIR" ]] || die "Missing --outdir"
[[ -s "$REF" ]] || die "Reference FASTA not found or empty: $REF"
[[ -s "$SAMPLES" ]] || die "Sample table not found or empty: $SAMPLES"
[[ "$THREADS" =~ ^[1-9][0-9]*$ ]] || die "--threads must be a positive integer"

command -v "$MINIMAP2" >/dev/null 2>&1 || die "minimap2 not found: $MINIMAP2"
command -v "$SAMTOOLS" >/dev/null 2>&1 || die "samtools not found: $SAMTOOLS"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
mkdir -p "$OUTDIR/ref" "$OUTDIR/mapping" "$OUTDIR/results"
cp "$REF" "$OUTDIR/ref/panel.fasta"
cp "$SAMPLES" "$OUTDIR/sample_inputs.tsv"

tail -n +2 "$SAMPLES" | while IFS=$'\t' read -r sample reads; do
  [[ -n "${sample:-}" ]] || continue
  [[ -s "$reads" ]] || die "Read file not found or empty for $sample: $reads"

  full_bam="$OUTDIR/mapping/${sample}.panel.primary.bam"
  mq20_bam="$OUTDIR/mapping/${sample}.panel.mq20_primary.bam"
  depth="$OUTDIR/mapping/${sample}.panel.mq20_primary.depth.tsv"
  read_counts="$OUTDIR/mapping/${sample}.panel.mq20_primary_read_counts.tsv"

  "$MINIMAP2" -a -x map-ont --secondary=no -t "$THREADS" "$REF" "$reads" \
    | "$SAMTOOLS" sort -@ 3 -o "$full_bam" -

  "$SAMTOOLS" index "$full_bam"
  "$SAMTOOLS" view -b -q 20 -F 2308 "$full_bam" -o "$mq20_bam"
  "$SAMTOOLS" index "$mq20_bam"
  "$SAMTOOLS" depth -Q 20 -G SUPPLEMENTARY "$mq20_bam" > "$depth"
  "$SAMTOOLS" view "$mq20_bam" \
    | awk 'BEGIN{OFS="\t"} {reads[$3]++; bases[$3]+=length($10)} END{for (r in reads) print r, reads[r], bases[r]}' \
    | sort -k1,1 > "$read_counts"
done

(
  cd "$OUTDIR"
  python3 "$SCRIPT_DIR/summarise_panel_mapping.py"
)
