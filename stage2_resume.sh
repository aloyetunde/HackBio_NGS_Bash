#!/bin/bash

# -----------------------------
# Stage 2: Assembly + AMR/Toxin detection
# -----------------------------

TRIM_DIR="trimmed"
ASSEMBLY_DIR="assembly"
REPORT_DIR="report/abricate"

mkdir -p $ASSEMBLY_DIR $REPORT_DIR

# -----------------------------
# 1. Assemble genomes with SPAdes (resume mode)
# -----------------------------
echo "🔧 Checking SPAdes assemblies..."

for fq1 in ${TRIM_DIR}/*_1.trimmed.fastq.gz; do
    fq2=${fq1/_1.trimmed.fastq.gz/_2.trimmed.fastq.gz}
    sample=$(basename $fq1 _1.trimmed.fastq.gz)
    outdir=${ASSEMBLY_DIR}/${sample}

    if [ -f "${outdir}/contigs.fasta" ]; then
        echo "✅ Assembly already exists for $sample, skipping..."
    else
        echo "🛠️ Assembling $sample..."
        spades.py -1 $fq1 -2 $fq2 -o $outdir --threads 8 --memory 32
    fi
done

# -----------------------------
# 2. Run Abricate (CARD DB for AMR)
# -----------------------------
echo "🧬 Running abricate (CARD DB)..."

for asm in ${ASSEMBLY_DIR}/*/contigs.fasta; do
    sample=$(basename $(dirname $asm))
    outfile=${REPORT_DIR}/${sample}_card.tab

    if [ -f "$outfile" ]; then
        echo "✅ AMR (CARD) already done for $sample, skipping..."
    else
        echo "🔎 Detecting AMR genes in $sample..."
        abricate --db card $asm > $outfile
    fi
done

abricate --summary ${REPORT_DIR}/*_card.tab > ${REPORT_DIR}/AMR_summary_CARD.txt

# -----------------------------
# 3. Run Abricate (VFDB DB for Toxins)
# -----------------------------
echo "☣️ Running abricate (VFDB DB)..."

for asm in ${ASSEMBLY_DIR}/*/contigs.fasta; do
    sample=$(basename $(dirname $asm))
    outfile=${REPORT_DIR}/${sample}_vfdb.tab

    if [ -f "$outfile" ]; then
        echo "✅ Toxin (VFDB) already done for $sample, skipping..."
    else
        echo "☢️ Detecting toxins in $sample..."
        abricate --db vfdb $asm > $outfile
    fi
done

abricate --summary ${REPORT_DIR}/*_vfdb.tab > ${REPORT_DIR}/Toxin_summary_VFDB.txt

# -----------------------------
# 4. Done
# -----------------------------
echo "🎉 Stage 2 (resume) complete!"
echo "AMR summary: ${REPORT_DIR}/AMR_summary_CARD.txt"
echo "Toxin summary: ${REPORT_DIR}/Toxin_summary_VFDB.txt"
