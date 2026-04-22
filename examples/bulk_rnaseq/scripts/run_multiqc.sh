#!/bin/bash
# run_multiqc.sh
# Summarizes FastQC reports using MultiQC

# --- CONFIGURATION ---
QC_DIR="qc"
FASTQC_DIR="$QC_DIR/fastqc"
MULTIQC_DIR="$QC_DIR/multiqc"

# --- CREATE OUTPUT DIRECTORY ---
mkdir -p "$MULTIQC_DIR"

# --- CHECK INPUT DIRECTORY ---
if [ ! -d "$FASTQC_DIR" ]; then
    echo "ERROR: FastQC directory not found: $FASTQC_DIR"
    exit 1
fi

# --- CHECK MULTIQC INSTALLATION ---
if ! command -v multiqc &> /dev/null; then
    echo "ERROR: MultiQC not installed or not in PATH"
    exit 1
fi

echo "Running MultiQC on $FASTQC_DIR ..."

multiqc "$FASTQC_DIR" -o "$MULTIQC_DIR" || {
    echo "MultiQC failed"
    exit 1
}

echo "=== MultiQC complete ==="
echo "Results in: $MULTIQC_DIR"