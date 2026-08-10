#!/usr/bin/env bash
set -euo pipefail

export LC_ALL=C
export LANG=C

usage() {
  cat <<'EOF'
Usage:
  run_sample_breseq.sh -s SAMPLE -r REF -o OUTDIR [-t THREADS] [--reads-dir READS_DIR]

Required:
  -s, --sample     Sample/read stem, e.g. 279991_pao1plps4
  -r, --ref        Reference FASTA or GenBank file
  -o, --outdir     Output directory

Optional:
  -t, --threads    Number of threads [default: 6]
      --reads-dir  Directory containing trimmed FASTQ files [default: reads]
  -h, --help       Show this help

Expected read files:
  READS_DIR/SAMPLE_1_trimmed.fastq.gz
  READS_DIR/SAMPLE_2_trimmed.fastq.gz
  READS_DIR/SAMPLE_U1_trimmed.fastq.gz
  READS_DIR/SAMPLE_U2_trimmed.fastq.gz
EOF
}

die() {
  echo "ERROR: $*" >&2
  exit 1
}

sample=""
threads="6"
outdir=""
ref=""
reads_dir="reads"

while [[ $# -gt 0 ]]; do
  case "$1" in
    -s|--sample)
      sample="${2:-}"
      shift 2
      ;;
    -t|--threads)
      threads="${2:-}"
      shift 2
      ;;
    -o|--outdir)
      outdir="${2:-}"
      shift 2
      ;;
    -r|--ref)
      ref="${2:-}"
      shift 2
      ;;
    --reads-dir)
      reads_dir="${2:-}"
      shift 2
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      die "Unknown option: $1"
      ;;
  esac
done

[[ -n "$sample" ]] || die "Missing -s/--sample"
[[ -n "$ref" ]] || die "Missing -r/--ref"
[[ -n "$outdir" ]] || die "Missing -o/--outdir"
[[ "$threads" =~ ^[1-9][0-9]*$ ]] || die "Threads must be a positive integer"

# Accept either a bare sample stem or a path-like stem such as reads/sample.
sample="${sample%/}"
sample="${sample##*/}"
reads_dir="${reads_dir%/}"

command -v breseq >/dev/null 2>&1 || die "'breseq' is not on PATH"

r1="${reads_dir}/${sample}_1_trimmed.fastq.gz"
r2="${reads_dir}/${sample}_2_trimmed.fastq.gz"
u1="${reads_dir}/${sample}_U1_trimmed.fastq.gz"
u2="${reads_dir}/${sample}_U2_trimmed.fastq.gz"

[[ -s "$ref" ]] || die "Reference not found or empty: $ref"
[[ -s "$r1" ]] || die "R1 FASTQ not found or empty: $r1"
[[ -s "$r2" ]] || die "R2 FASTQ not found or empty: $r2"
[[ -s "$u1" ]] || die "U1 FASTQ not found or empty: $u1"
[[ -s "$u2" ]] || die "U2 FASTQ not found or empty: $u2"

sample_dir="${outdir%/}/${sample}"
log="$sample_dir/${sample}.breseq.pipeline.log"

mkdir -p "$sample_dir"

{
  echo "Sample:    $sample"
  echo "Reference: $ref"
  echo "R1:        $r1"
  echo "R2:        $r2"
  echo "U1:        $u1"
  echo "U2:        $u2"
  echo "Threads:   $threads"
  echo "Output:    $sample_dir"
  echo
} | tee "$log"

echo "Running breseq with paired and unpaired reads..." | tee -a "$log"
breseq \
  -j "$threads" \
  -o "$sample_dir" \
  -r "$ref" \
  "$r1" "$r2" "$u1" "$u2" \
  2>&1 | tee -a "$log"

echo
echo "Done."
echo "breseq output: $sample_dir"
