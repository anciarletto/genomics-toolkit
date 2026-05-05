#!/bin/bash
# scripts/summarize_featurecounts.sh
#
# Purpose:
#   Extracts gene assignment statistics from featureCounts summary files
#   and consolidates them into a single QC table.
#
# Function:
#   Parses "Assigned" read counts from featureCounts output summaries
#   across all samples to evaluate library complexity and gene detection efficiency.
#
# Output:
#   - qc/assigned_reads_summary.txt
#
# Key QC metric:
#   - Number of reads successfully assigned to annotated genes
#
# Dependencies:
#   - bash
#   - grep
#
# Notes:
#   - Requires featureCounts summary files in counts/
#   - Used for post-quantification library quality assessment
#   - Can be extended to include Unassigned categories (e.g., Multi-mapping, No feature)

grep "Assigned" counts/*summary > qc/assigned_reads_summary.txt