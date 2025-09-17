#!/bin/bash

# 1. Define directories
# -----------------------------
TRIM_DIR="trimmed"
ASSEMBLY_DIR="assembly"
REPORT_DIR="report/abricate"

mkdir -p $ASSEMBLY_DIR $REPORT_DIR

# -----------------------------
# 2. Assemble genomes with SPAdes
# -----------------------------
echo "🔧 Running SPAdes assembly for all samples..."

for fq1 in ${TRIM_DIR}/*_1.trimmed.fastq.gz; do
    fq2=${fq1/_1.trimmed.fastq.gz/_2.trimmed.fastq.gz}
    sample=$(basename $fq1 _1.trimmed.fastq.gz)

    echo "Assembling sample: $sample"
    spades.py -1 $fq1 -2 $fq2 -o ${ASSEMBLY_DIR}/${sample} --threads 8 --memory 32
done

# -----------------------------
# 3. Run Abricate (CARD DB)
# -----------------------------
echo "🧬 Running abricate (CARD DB) for AMR gene detection..."

for asm in ${ASSEMBLY_DIR}/*/contigs.fasta; do
    sample=$(basename $(dirname $asm))
    abricate --db card $asm > ${REPORT_DIR}/${sample}_card.tab
done

# Summarize across all samples
abricate --summary ${REPORT_DIR}/*_card.tab > ${REPORT_DIR}/AMR_summary_CARD.txt

# -----------------------------
# 4. Run Abricate (VFDB DB - Toxin genes)
# -----------------------------
echo "☣️ Running abricate (VFDB DB) for toxin gene detection..."

for asm in ${ASSEMBLY_DIR}/*/contigs.fasta; do
    sample=$(basename $(dirname $asm))
    abricate --db vfdb $asm > ${REPORT_DIR}/${sample}_vfdb.tab
done

# Summarize toxin results
abricate --summary ${REPORT_DIR}/*_vfdb.tab > ${REPORT_DIR}/Toxin_summary_VFDB.txt

# -----------------------------
# 5. Final message
# -----------------------------
echo "✅ Stage 2 completed!"
echo "AMR summary: ${REPORT_DIR}/AMR_summary_CARD.txt"
echo "Toxin summary: ${REPORT_DIR}/Toxin_summary_VFDB.txt"
