#!/bin/bash
# Script to count reads per gene using featureCounts

# create output directory if not exists
mkdir -p counts

# run featureCounts
featureCounts -O -t gene -g ID \
  -a genome/c_elegans.gff3 \
  -o counts/counts.txt IGV/*.bam

echo "✅ Counting complete. Results saved in counts/counts.txt"
