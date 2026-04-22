#!/bin/bash
set -euo pipefail

# =========================================================
# download_fastq.sh
#
# Purpose:
#   Downloads SRA sequencing data using a list of SRR
#   accessions, converts them to FASTQ, compresses output,
#   and performs basic cleanup.
#
# Inputs:
#   - raw_data/SRR_Acc_List.txt (list of SRR IDs)
#
# Outputs:
#   - raw_data/fastq/*.fastq.gz
#   - raw_data/download_log.txt
#
# Requirements:
#   - SRA Toolkit (prefetch, fasterq-dump)
#   - pigz (optional, for compression)
#
# Notes:
#   - Uses --split-files for paired-end reads
#   - Removes SRA cache after processing each sample
#   - Designed for reproducible bulk RNA-seq workflows
# =========================================================

# --- CONFIGURATION ---
RAW_DIR="raw_data"
SRR_LIST="$RAW_DIR/SRR_Acc_List.txt"

FASTQ_DIR="$RAW_DIR/fastq"
LOG_FILE="$RAW_DIR/download_log.txt"

mkdir -p "$FASTQ_DIR"
touch "$LOG_FILE"

# --- PROCESS ---
while IFS= read -r srr; do
    [[ -z "$srr" || "$srr" == \#* ]] && continue

    echo "=== Processing $srr ==="

    # download
    if ! prefetch "$srr"; then
        echo "prefetch FAILED: $srr" >> "$LOG_FILE"
        continue
    fi

    # convert
    if ! fasterq-dump "$srr" --split-files -O "$FASTQ_DIR"; then
        echo "fasterq-dump FAILED: $srr" >> "$LOG_FILE"
        continue
    fi

    # verify
    if ! ls "$FASTQ_DIR/${srr}"*.fastq* >/dev/null 2>&1; then
        echo "WARNING missing FASTQ: $srr" >> "$LOG_FILE"
        continue
    fi

    # compress
    if command -v pigz >/dev/null 2>&1; then
        pigz "$FASTQ_DIR/${srr}"*.fastq
    fi

    # cleanup SRA cache
    rm -rf "$HOME/ncbi/public/sra/$srr"

    echo "$srr completed" >> "$LOG_FILE"

done < "$SRR_LIST"

echo "All downloads complete"