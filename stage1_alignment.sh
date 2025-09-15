#!/bin/bash
# stage1_alignment.sh
# Align trimmed reads to reference genome and generate sorted, indexed BAM files

# Exit if any command fails
set -euo pipefail

# Paths
REF="reference/reference.fasta"
TRIMMED_DIR="trimmed"
MAPPED_DIR="mapped"

# Loop over all forward reads
for fq1 in ${TRIMMED_DIR}/*_1.trimmed.fastq.gz; do
    sample=$(basename "$fq1" _1.trimmed.fastq.gz)
    fq2="${TRIMMED_DIR}/${sample}_2.trimmed.fastq.gz"

    echo ">>> Aligning sample: $sample"

    bwa mem "$REF" "$fq1" "$fq2" | \
        samtools view -bS - | samtools sort -o ${MAPPED_DIR}/${sample}.sorted.bam

    samtools index ${MAPPED_DIR}/${sample}.sorted.bam

    echo ">>> Finished: ${sample}"
    echo
done

echo "✅ Alignment completed for all samples."
