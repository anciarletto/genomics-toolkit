#!/bin/bash
set -euo pipefail  # Pipeline stops if HISAT2 fails

# =========================================================
# run_hisat2.sh
#
# Purpose:
#   Aligns paired-end RNA-seq reads to the GRCh38 genome
#   using HISAT2 and outputs sorted BAM files + indexes.
#
# Outputs:
#   - alignments/*.sorted.bam
#   - alignments/*.sorted.bam.bai
#   - alignments/*_summary.txt
#
# Key Steps:
#   1. Download HISAT2 index if missing
#   2. Align reads using HISAT2
#   3. Sort BAM with samtools
#   4. Index BAM for downstream QC/tools
# =========================================================

# --- CONFIGURATION ---
SCRIPT_DIR=$(dirname "$(realpath "$0")")
REPO_ROOT="$SCRIPT_DIR/.."

FASTQ_DIR="$REPO_ROOT/examples/bulk_rnaseq/fastq"
ALIGN_DIR="$REPO_ROOT/examples/bulk_rnaseq/alignments"
REFERENCE_DIR="$REPO_ROOT/examples/bulk_rnaseq/reference"
INDEX_DIR="$REFERENCE_DIR/grch38"
INDEX_PREFIX="$INDEX_DIR/genome"

GENOME_URL="https://genome-idx.s3.amazonaws.com/hisat/grch38_genome.tar.gz"

mkdir -p "$ALIGN_DIR"
mkdir -p "$INDEX_DIR"

# STEP 1: Download HISAT2 index if missing
if [ ! -f "${INDEX_PREFIX}.1.ht2" ]; then
    echo "HISAT2 index not found. Downloading..."
    wget -O "$INDEX_DIR/grch38_genome.tar.gz" "$GENOME_URL"
    tar -xzf "$INDEX_DIR/grch38_genome.tar.gz" -C "$INDEX_DIR"
    echo "HISAT2 index ready."
else
    echo "HISAT2 index already exists. Skipping download."
fi

# STEP 2: Align + sort + index
THREADS=6

for fq1 in "$FASTQ_DIR"/*_1.fastq*; do
    [ -e "$fq1" ] || { echo "No FASTQ files found in $FASTQ_DIR"; exit 1; }

    base=$(basename "$fq1" "_1.fastq.gz")
    fq2="$FASTQ_DIR/${base}_2.fastq.gz"

    out_bam="$ALIGN_DIR/${base}.sorted.bam"
    summary_file="$ALIGN_DIR/${base}_summary.txt"

    # index file path
    bam_index="$out_bam.bai"

    if [ -f "$out_bam" ]; then
        echo "$out_bam already exists. Skipping alignment."
    else
        echo "=== Aligning $base ==="

        hisat2 -p "$THREADS" \
            -x "$INDEX_PREFIX" \
            -1 "$fq1" -2 "$fq2" \
            --summary-file "$summary_file" \
        | samtools sort -@ "$THREADS" -o "$out_bam"
    fi

    # BAM indexing (critical for samtools idxstats, IGV, QC)
    if [ ! -f "$bam_index" ]; then
        echo "Indexing $out_bam ..."
        samtools index "$out_bam"
    else
        echo "Index already exists for $base"
    fi

done

echo "=== Alignment + indexing complete. BAM files ready in $ALIGN_DIR ==="