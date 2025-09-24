#!/bin/bash
# trimming files
mkdir -p trim

for filename in clean_files/*.fq; do
    base=$(basename "$filename" .fq)
    fastp -i "$filename" -o "trim/${base}_trim.fq" -h "trim/${base}_report.html"
done

# Run MultiQC on the trim folder to summarize fastp reports
multiqc trim -o trim
