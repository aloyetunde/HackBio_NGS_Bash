#!/bin/bash
set -euo pipefail

echo "=== QC on trimmed reads started at $(date) ==="

mkdir -p qc/trimmed_qc report

# Run FastQC on all trimmed files
fastqc trimmed/*.fastq.gz -o qc/trimmed_qc

# Run MultiQC to aggregate results
multiqc qc/trimmed_qc -o report

echo "✅ QC on trimmed reads completed at $(date)"
