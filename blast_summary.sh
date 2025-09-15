#!/bin/bash
# blast_summary.sh
# Summarize BLAST results into a single CSV table

set -euo pipefail

echo "Sample,Top_Hit,%_Identity,Alignment_Length,Query_Length,Subject_Length,Evalue,Bitscore" > report/blast_summary.csv

for file in report/*.blast.out
do
    base=$(basename "$file" .blast.out)
    # Take the best hit (first line, highest bitscore)
    top_hit=$(head -n 1 "$file")
    echo "$base,$top_hit" >> report/blast_summary.csv
done

echo "✅ BLAST summary saved to report/blast_summary.csv"
