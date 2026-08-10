#!/usr/bin/env python3
import csv
import math
from pathlib import Path


PANEL_FASTA = Path("ref/panel.fasta")
SAMPLE_TABLE = Path("sample_inputs.tsv")
MAPPING_DIR = Path("mapping")
RESULTS_DIR = Path("results")
OUTPUT = RESULTS_DIR / "present_max50pct_10kb_mq20_primary_matrix.tsv"


def read_fasta_lengths(path):
    names = []
    lengths = {}
    name = None
    length = 0

    with open(path) as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    names.append(name)
                    lengths[name] = length
                name = line[1:].split()[0]
                length = 0
            else:
                length += len(line)

    if name is not None:
        names.append(name)
        lengths[name] = length

    return names, lengths


def read_samples(path):
    with open(path) as handle:
        return [row["sample"] for row in csv.DictReader(handle, delimiter="\t")]


def read_mapped_read_counts(path):
    counts = {}
    if not path.exists() or path.stat().st_size == 0:
        return counts

    with open(path) as handle:
        for row in csv.reader(handle, delimiter="\t"):
            if len(row) >= 2:
                counts[row[0]] = int(row[1])

    return counts


def read_covered_bases(path, refs):
    covered = {ref: 0 for ref in refs}
    if not path.exists() or path.stat().st_size == 0:
        return covered

    with open(path) as handle:
        for row in csv.reader(handle, delimiter="\t"):
            if len(row) < 3:
                continue
            ref = row[0]
            depth = int(row[2])
            if ref in covered and depth >= 1:
                covered[ref] += 1

    return covered


refs, lengths = read_fasta_lengths(PANEL_FASTA)
samples = read_samples(SAMPLE_TABLE)
presence = {}

for sample in samples:
    read_counts = read_mapped_read_counts(
        MAPPING_DIR / f"{sample}.panel.mq20_primary_read_counts.tsv"
    )
    covered_bases = read_covered_bases(
        MAPPING_DIR / f"{sample}.panel.mq20_primary.depth.tsv",
        refs,
    )

    for ref in refs:
        required_covered_bases = max(math.ceil(0.50 * lengths[ref]), 10_000)
        presence[(sample, ref)] = int(
            covered_bases[ref] >= required_covered_bases
            and read_counts.get(ref, 0) >= 5
        )

RESULTS_DIR.mkdir(exist_ok=True)

with open(OUTPUT, "w", newline="") as handle:
    writer = csv.writer(handle, delimiter="\t")
    writer.writerow(["sample", *refs])
    for sample in samples:
        writer.writerow([sample, *[presence[(sample, ref)] for ref in refs]])
