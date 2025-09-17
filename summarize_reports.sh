#!/bin/bash
# summarize_reports.sh
# Script to clean abricate AMR & toxin summaries and compute prevalence

set -euo pipefail

cd report/abricate || { echo "abricate report directory not found"; exit 1; }

### --- Function to process a summary file ---
process_summary() {
    local infile=$1
    local prefix=$2

    echo "Processing $infile ..."

    # Clean spacing and make CSV
    awk 'NR==1{gsub(/\t/," "); gsub(/  +/," "); print; next} {gsub(/\t/," "); gsub(/  +/," "); print}' "$infile" > "${prefix}.cleaned.txt"
    sed 's/ \+/,/g' "${prefix}.cleaned.txt" > "${prefix}.csv"

    # Run python prevalence calculation
    python3 - <<PY > ${prefix}_prevalence.csv
import csv
from collections import Counter

rows=[]
with open('${prefix}.csv') as f:
    rdr=csv.reader(f)
    header=next(rdr)
    genes = header[2:]
    for r in rdr:
        if not r: continue
        r += ["0"]*(len(header)-len(r))
        rows.append(r)

counts = Counter()
for r in rows:
    for i,g in enumerate(genes, start=2):
        try:
            if float(r[i]) > 0:
                counts[g] += 1
        except Exception:
            pass

print("Gene,Count,Prevalence_percent")
n = len(rows)
for g in genes:
    c = counts[g]
    print(f"{g},{c},{(c/n*100):.1f}")
PY
}

### --- Run for AMR (CARD) ---
process_summary AMR_summary_CARD.txt AMR_summary_CARD

### --- Run for Toxins (VFDB) ---
process_summary Toxin_summary_VFDB.txt Toxin_summary_VFDB
