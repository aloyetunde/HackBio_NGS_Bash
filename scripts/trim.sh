#!/bin/bash
set -euo pipefail

# Directories
RAW="raw"
TRIM="trimmed"
QC="qc"
LOG="logs"

mkdir -p "$TRIM" "$QC" "$LOG"

# Log file with date & time
LOGFILE="$LOG/trim_$(date +'%Y%m%d_%H%M%S').log"

echo "=== Trimming started at $(date) ===" | tee -a "$LOGFILE"

for R1 in ${RAW}/*_1.fastq.gz; do
    SAMPLE=$(basename "$R1" | sed 's/_.*//')
    R2="${RAW}/${SAMPLE}_Genome_Sequencing_of_Listeria_monocytogenes_SA_outbreak_2017_2.fastq.gz"

    echo "Processing sample: $SAMPLE" | tee -a "$LOGFILE"

    # Check if R2 exists
    if [[ ! -f "$R2" ]]; then
        echo "⚠️  Missing file: $R2" | tee -a "$LOGFILE"
        echo "Re-running download script..." | tee -a "$LOGFILE"
        bash SA_Polony_100_download.sh
    fi

    # Check again after re-download
    if [[ ! -f "$R2" ]]; then
        echo "❌ Still missing $R2. Skipping $SAMPLE." | tee -a "$LOGFILE"
        continue
    fi

    # Run fastp
    fastp \
      -i "$R1" \
      -I "$R2" \
      -o "${TRIM}/${SAMPLE}_1.trimmed.fastq.gz" \
      -O "${TRIM}/${SAMPLE}_2.trimmed.fastq.gz" \
      -h "${QC}/${SAMPLE}_fastp.html" \
      -j "${QC}/${SAMPLE}_fastp.json" \
      2>&1 | tee -a "$LOGFILE"

    echo "✅ Finished trimming $SAMPLE" | tee -a "$LOGFILE"
done

echo "=== Trimming completed at $(date) ===" | tee -a "$LOGFILE"
