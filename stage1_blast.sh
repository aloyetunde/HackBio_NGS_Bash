#!/bin/bash
# stage1_blast.sh
# BLAST organism identification for SA Polony Project

set -euo pipefail

# Directories
mkdir -p reference fasta_reads report

# Step 1: Create BLAST database (only once)
makeblastdb -in NZ_CP007492.1.fasta -dbtype nucl -out reference/polony_db

# Step 2: Loop through all FASTA reads and run BLAST
for file in fasta_reads/*.fasta
do
    base=$(basename "$file" .fasta)
    echo "Running BLAST for $base ..."
    blastn -query "$file" \
           -db reference/polony_db \
           -out report/${base}.blast.out \
           -evalue 1e-10 -outfmt "6 qseqid sseqid pident length qlen slen evalue bitscore" \
           -num_threads 4
done

echo "✅ BLAST finished for all samples."


