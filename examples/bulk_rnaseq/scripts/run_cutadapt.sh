#!/bin/bash
# run_cutadapt.sh
# Trims adapters and low-quality bases from FASTQ files using Cutadapt

# --- CONFIGURATION ---
INPUT_DIR="fastq"             # Folder with original FASTQ files
OUTPUT_DIR="trimmed_fastq"    # Folder for trimmed FASTQ files
ADAPTER_SEQ="AGATCGGAAGAGC"   # Common Illumina adapter sequence (adjust if needed based on multiQC report)
QUALITY=30                    # Minimum Phred quality score to keep
MIN_LEN=20                     # Minimum read length to keep

mkdir -p "$OUTPUT_DIR"

# Check if Cutadapt is installed
if ! command -v cutadapt &> /dev/null; then
    echo "Cutadapt could not be found. Please install it and ensure it's in your PATH."
    exit 1
fi

# --- PROCESS FASTQ FILES ---
for fq in "$INPUT_DIR"/*.fastq*; do
    [ -e "$fq" ] || { echo "No FASTQ files found in $INPUT_DIR"; exit 1; }

    base=$(basename "$fq")
    output_file="$OUTPUT_DIR/$base"

    echo "Trimming $fq ..."
    cutadapt -a "$ADAPTER_SEQ" -q "$QUALITY" -m "$MIN_LEN" -o "$output_file" "$fq" \
        || { echo "Cutadapt failed for $fq"; continue; }
done

echo "=== Cutadapt trimming complete. Trimmed files in $OUTPUT_DIR ==="