#!/bin/bash
# run_fastqc.sh
# Runs FastQC on all FASTQ files and stores outputs in qc/fastqc

# --- CONFIGURATION ---
FASTQ_DIR="raw_data/fastq"
QC_DIR="qc"
FASTQC_DIR="$QC_DIR/fastqc"

# --- CREATE DIRECTORIES IF THEY DON'T EXIST ---
mkdir -p "$FASTQC_DIR"

# --- CHECK INPUT DIRECTORY ---
if [ ! -d "$FASTQ_DIR" ]; then
    echo "ERROR: FASTQ directory '$FASTQ_DIR' not found"
    exit 1
fi

# --- FIND FASTQ FILES ---
shopt -s nullglob
FASTQ_FILES=("$FASTQ_DIR"/*.fastq.gz)

if [ ${#FASTQ_FILES[@]} -eq 0 ]; then
    echo "ERROR: No FASTQ files found in $FASTQ_DIR"
    exit 1
fi

# --- RUN FASTQC ---
for fq in "${FASTQ_FILES[@]}"; do
    echo "Running FastQC on $fq ..."
    fastqc -o "$FASTQC_DIR" "$fq" || {
        echo "FastQC failed for $fq"
        continue
    }
done

echo "=== FastQC completed ==="
echo "Results in: $FASTQC_DIR"