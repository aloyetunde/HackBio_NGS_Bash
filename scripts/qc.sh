#!/bin/bash
# Quality Control Script for NGS data


mkdir -p qc

#  FastQC on all raw FASTQ files
fastqc raw/*.fastq.gz -o qc/

# MultiQC to summarize FastQC reports
multiqc qc/ -o report/
