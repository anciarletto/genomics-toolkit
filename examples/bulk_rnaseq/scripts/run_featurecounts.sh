#!/bin/bash
set -euo pipefail

# =========================================================
# run_featurecounts.sh
#
# Purpose:
#   Generates gene-level raw read counts from aligned BAM
#   files using Subread featureCounts.
#
# Inputs:
#   - BAM files: alignments/*.sorted.bam
#   - Gene annotation: Ensembl GTF (GRCh38)
#
# Outputs:
#   - counts/gene_counts.txt (raw count matrix)
#
# Key Steps:
#   1. Downloads Ensembl GTF annotation if missing
#   2. Uncompresses and prepares annotation file
#   3. Counts reads overlapping genes using featureCounts
#
# Key Parameters:
#   - Paired-end mode (-p)
#   - Requires both mates mapped (-B)
#   - Removes chimeric fragments (-C)
#   - Unstranded protocol (-s 0)
#
# Dependencies:
#   - featureCounts (Subread)
#   - wget
#   - gunzip
#
# Reference:
#   Ensembl GRCh38 release 109
#
# Notes:
#   - Assumes unstranded RNA-seq library prep
#   - Uses gene-level summarization (not exon-level)
#   - Compatible with DESeq2 downstream analysis
# =========================================================

SCRIPT_DIR=$(dirname "$(realpath "$0")")
REPO_ROOT="$SCRIPT_DIR/.."

BAM_DIR="$REPO_ROOT/examples/bulk_rnaseq/alignments"
OUTPUT_DIR="$REPO_ROOT/examples/bulk_rnaseq/counts"
REFERENCE_DIR="$REPO_ROOT/examples/bulk_rnaseq/reference"

GTF_GZ="$REFERENCE_DIR/Homo_sapiens.GRCh38.109.gtf.gz"
GTF="$REFERENCE_DIR/Homo_sapiens.GRCh38.109.gtf"

THREADS=4

mkdir -p "$OUTPUT_DIR"
mkdir -p "$REFERENCE_DIR"

# Check dependencies
for tool in featureCounts wget; do
    if ! command -v $tool &> /dev/null; then
        echo "$tool not found"
        exit 1
    fi
done

# Download GTF if needed
if [ ! -f "$GTF" ]; then
    echo "Downloading GTF..."

    if [ ! -f "$GTF_GZ" ]; then
        wget -O "$GTF_GZ" \
        ftp://ftp.ensembl.org/pub/release-109/gtf/homo_sapiens/Homo_sapiens.GRCh38.109.gtf.gz
    fi

    gunzip -f "$GTF_GZ"
fi

# Check BAMs
if ! ls "$BAM_DIR"/*.sorted.bam 1> /dev/null 2>&1; then
    echo "No BAM files found"
    exit 1
fi

echo "Running featureCounts..."

featureCounts \
    -T "$THREADS" \
    -a "$GTF" \
    -o "$OUTPUT_DIR/gene_counts.txt" \
    -p \
    -B \
    -C \
    -s 0 \
    "$BAM_DIR"/*.sorted.bam

echo "=== Done ==="