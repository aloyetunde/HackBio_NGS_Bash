#!/bin/bash
# 01_qc.sh - Quality control of raw FASTQ files

# Exit if any command fails
set -e

# Timestamp for log files
timestamp=$(date +"%Y%m%d_%H%M%S")

# Define directories
raw_dir="raw"
qc_dir="qc"
log_dir="logs"

# Create output directories if missing
mkdir -p $qc_dir $log_dir

echo "[$(date)] Starting FastQC..." | tee $log_dir/qc_$timestamp.log

# Run FastQC on all raw FASTQ files
fastqc $raw_dir/*.fastq.gz -o $qc_dir 2>&1 | tee -a $log_dir/qc_$timestamp.log

echo "[$(date)] FastQC completed." | tee -a $log_dir/qc_$timestamp.log

# Run MultiQC to aggregate results
echo "[$(date)] Running MultiQC..." | tee -a $log_dir/qc_$timestamp.log
multiqc $qc_dir -o $qc_dir 2>&1 | tee -a $log_dir/qc_$timestamp.log

echo "[$(date)] MultiQC completed." | tee -a $log_dir/qc_$timestamp.log
echo "[$(date)] QC pipeline finished successfully." | tee -a $log_dir/qc_$timestamp.log
