#!/bin/bash

# =========================================================
# This script is used to perform quality control for our RNA-seq data for C. elegans
# =========================================================

# Step 1: Create output directory
mkdir -p qc

# Step 2: Run FastQC on each .fq file
for filename in clean_files/*.fq; do
    echo "🔍 Running FastQC on $filename ..."
    fastqc -o qc/ "$filename"
done

echo "✅ Quality control complete. Results saved in qc/"

# Step 3: Git add, commit, and push results
git add qc
git commit -m "Added FastQC reports for RNA-Seq data: $(date '+%Y-%m-%d %H:%M:%S')"
git push origin main
