# S_aureus_PJI_RNAseq — Full Pipeline & Report

**Project:** Transcriptomic profiling of *Staphylococcus aureus* during acute vs. chronic PJI  

**Repository path (working directory for commands):**  
`~/Mary/HackBio_NGS_Bash/S_aureus_PJI_RNAseq/`

**Metadata file (full path):**  
`~/Mary/HackBio_NGS_Bash/S_aureus_PJI_RNAseq/metadata.tsv`

# Table of Contents

1. **Project Summary**  
2. **Files & Folders Created (What You Have Now)**  
3. **Exact Commands & Scripts**  
   - Download → QC → Trimming → Mapping → Counts → DE  
   - `scripts/qc_script.sh`  
   - `scripts/trim_script.sh`  
   - `scripts/mapping_script.sh`  
   - `featureCounts` command used to generate `counts_USA300_clean.txt`  
   - `scripts/deseq2_USA300.R` (full DESeq2 R script without EnhancedVolcano)  
4. **Troubleshooting**  
   - Sample-name mismatch: commands we used  
5. **What to Push to GitHub & `.gitignore` (Exact Content)**  
6. **Short Biological / Results Notes & How to Create the Clinical Report**  
7. **Useful References & Final Checklist**

# 1 — Project Summary (Short)

We profiled **8 RNA-seq samples** from **PRJNA867318** (4 chronic, 4 acute).  

### Steps performed:
1. **Download SRA FASTQs** → stored under `raw_sra/`  
2. **Pre-QC with FastQC** → aggregated with **MultiQC**  
3. **Trim adapters / low-quality bases with fastp** → output in `trim/`  
4. **Post-trim QC** → FastQC + MultiQC  
5. **Index + align trimmed reads** to USA300 (FPR3757) reference using `bwa mem` → `mapped_reads/`  
6. **Index BAMs for IGV** → stored in `IGV/`  
7. **Count reads per gene** with `featureCounts` → `counts/counts_USA300_clean.txt`  
8. **Differential expression with DESeq2** (`scripts/deseq2_USA300.R`) → results: PCA, volcano, heatmap in `results/`  

# 2 — Files / Folders (What We Created)

Typical content of project root:

S_aureus_PJI_RNAseq/
├─ raw_sra/                      # downloaded SRA fastq.gz files (do NOT push)
├─ clean_fastq/                  # intermediate cleaned fastq (optional)
├─ trim/                         # output from fastp (trimmed fastq.gz)
├─ mapped_reads/                 # sorted .bam (do NOT push binary BAMs)
├─ IGV/                          # BAMs + index for IGV viewing (do NOT push big files)
├─ genome/                       # reference FASTA & annotation (do NOT push big files)
├─ genome_index/                 # bwa/bowtie/STAR indexes (do NOT push)
├─ counts/
│   ├─ counts_USA300_clean.txt
│   └─ counts_USA300_clean.txt.summary
├─ qc/
│   ├─ fastqc/                    # fastqc outputs (safe to push)
│   └─ multiqc_report.html
├─ scripts/
│   ├─ qc_script.sh
│   ├─ trim_script.sh
│   ├─ mapping_script.sh
│   └─ deseq2_USA300.R
├─ results/
│   ├─ PCA_USA300.pdf
│   ├─ Volcano_USA300.pdf
│   ├─ Heatmap_USA300.pdf
│   └─ DESeq2_results_USA300.csv
├─ metadata.tsv
└─ S_aureus_PJI_RNAseq_REPORT.md  # this file
