# S_aureus_PJI_RNAseq — Full Pipeline & Report

**Project:** Transcriptomic profiling of *Staphylococcus aureus* during acute vs. chronic PJI  
## Background and Rationale

Periprosthetic joint infections (PJIs) are severe complications following orthopedic implant surgery. *Staphylococcus aureus* is a leading pathogen in these infections, capable of switching between **acute (toxin-producing, planktonic)** and **chronic (biofilm-adapted, stress-tolerant)** phases. Understanding transcriptional changes underlying this shift is critical for diagnostics and therapeutic targeting.

RNA-seq was employed to compare gene expression profiles of *S. aureus* isolates from acute vs chronic PJIs.

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

# Project Root Structure

```text
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
```

# 3 — Exact Commands & Scripts

Below are the exact scripts & commands we used (copy/paste).  
The scripts are intended to be placed in `scripts/` and made executable:

```bash
chmod +x scripts/*.sh
```

## 3.0 — Quick: Clone Repo / Change to Project Folder

```bash
# on my server 
cd ~/Mary/HackBio_NGS_Bash

git clone https://github.com/aloyetunde/HackBio_NGS_Bash.git

cd S_aureus_PJI_RNAseq
```
## 3.1 — Download SRA in `raw_sra/`

```bash
# create folder
mkdir -p raw_sra

# Used SRA-Explorer generated direct URLs with wget
wget -P raw_sra <link-to-SRR20959676_1.fastq.gz>
wget -P raw_sra <link-to-SRR20959676_2.fastq.gz>
```
### `raw_sra/` Files

SRR20959676_1.fastq.gz SRR20959678_2.fastq.gz SRR20959681_1.fastq.gz
SRR20959676_2.fastq.gz SRR20959679_1.fastq.gz SRR20959681_2.fastq.gz
SRR20959677_1.fastq.gz SRR20959679_2.fastq.gz SRR20959682_1.fastq.gz
SRR20959677_2.fastq.gz SRR20959680_1.fastq.gz SRR20959682_2.fastq.gz
SRR20959678_1.fastq.gz SRR20959680_2.fastq.gz

## 3.2 — QC (FastQC + MultiQC) — `scripts/qc_script.sh`

Create `scripts/qc_script.sh`:

```bash
#!/bin/bash
set -euo pipefail

mkdir -p qc/fastqc

# Run FastQC for all raw fastq.gz files
fastqc -o qc/fastqc raw_sra/*.fastq.gz

# Aggregate with MultiQC
multiqc qc/fastqc -o qc/
```
### Run the script:
```bash
chmod +x scripts/qc_script.sh
./scripts/qc_script.sh
```
### QC Outputs:

`qc/fastqc/*html`

`qc/multiqc_report.html`

## 3.3 — Trimming with fastp — `scripts/trim_script.sh`

Create `scripts/trim_script.sh` 
```bash
#!/bin/bash
set -euo pipefail

mkdir -p trim

for r1 in raw_sra/*_1.fastq.gz; do
  r2=${r1/_1.fastq.gz/_2.fastq.gz}
  sample=$(basename "$r1" _1.fastq.gz)
  echo "Trimming $sample ..."
  fastp -i "$r1" -I "$r2" \
        -o "trim/${sample}_R1_trim.fastq.gz" -O "trim/${sample}_R2_trim.fastq.gz" \
        -h "trim/${sample}_fastp.html" -j "trim/${sample}_fastp.json" \
        --detect_adapter_for_pe -w 4
done
```
### Run Trimming Script

```bash
chmod +x scripts/trim_script.sh
./scripts/trim_script.sh
```
### Outputs:

`trim/*_R1_trim.fastq.gz`

`trim/*_R2_trim.fastq.gz`

`trim/*_fastp.html`

## 3.4 — Post-trim QC (FastQC + MultiQC)

```bash
# Create folder for post-trim QC
mkdir -p qc/post_trim

# Run FastQC on trimmed reads
fastqc -o qc/post_trim trim/*.fastq.gz

# Aggregate with MultiQC
multiqc qc/post_trim -o qc/post_trim/
```

## 3.5 — Download USA300 Reference & Annotation (if needed)

```bash
# Create genome folder
mkdir -p genome

# Replace URL with the RefSeq you chose (example below uses USA300 / GCF_000013465.1)
wget -O genome/s_aureus_USA300.fa.gz \
  ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/013/465/GCF_000013465.1_ASM1346v1/GCF_000013465.1_ASM1346v1_genomic.fna.gz

wget -O genome/USA300_FPR3757_genomic.gff3.gz \
  ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/013/465/GCF_000013465.1_ASM1346v1/GCF_000013465.1_ASM1346v1_genomic.gff.gz

# Unzip downloaded files
gunzip -f genome/*.gz
```
```text
# Final files:
# genome/s_aureus_USA300.fa
# genome/USA300_FPR3757_genomic.gff3
```
### 3.6 — Build BWA Index for Mapping

```bash
cd genome
bwa index s_aureus_USA300.fa
cd ..
```

## 3.7 — Mapping Trimmed Reads (Paired-End) — `scripts/mapping_script.sh`

Create `scripts/mapping_script.sh`:

```bash
#!/bin/bash
set -euo pipefail

mkdir -p mapped_reads

# Loop over samples based on trimmed filenames
for r1 in trim/*_R1_trim.fastq.gz; do
  sample=$(basename "$r1" _R1_trim.fastq.gz)
  r2="trim/${sample}_R2_trim.fastq.gz"
  echo "Mapping ${sample} ..."
  bwa mem -t 8 genome/s_aureus_USA300.fa "$r1" "$r2" \
    | samtools view -bS - \
    | samtools sort -o mapped_reads/${sample}.bam -
  samtools index mapped_reads/${sample}.bam
done
```
Excute
```bash
chmod +x scripts/mapping_script.sh
./scripts/mapping_script.sh
```
### Mapping Output

- `mapped_reads/{sample}.bam`  
- `mapped_reads/{sample}.bam.bai`
---

### Quick Mapping Check

```bash
samtools flagstat mapped_reads/fm_1.bam > mapped_reads/fm_1.flagstat
cat mapped_reads/fm_1.flagstat
```
## 3.8 — Move BAMs to IGV Folder

```bash
# Create IGV folder
mkdir -p IGV

# Move BAMs and indexes
mv mapped_reads/*.bam IGV/
mv mapped_reads/*.bai IGV/
```
`#  kept the indexes for visualization because of the size`

`Index BAMs `

```bash
for bam in IGV/*.bam; do
  samtools index "$bam"
done
```
## 3.9 — featureCounts: Generate Gene Counts

Ran `featureCounts` with paired-end flags `-p -B -C` and counted at gene level using the `locus_tag` attribute in the GFF.

```bash
mkdir -p counts

featureCounts -p -B -C -O -t gene -g locus_tag \
  -a genome/USA300_FPR3757_genomic.gff3 \
  -o counts/counts_USA300_clean.txt \
  mapped_reads/fm_1.bam mapped_reads/fm_2.bam mapped_reads/fm_3.bam mapped_reads/fm_4.bam \
  mapped_reads/m_1.bam mapped_reads/m_2.bam mapped_reads/m_3.bam mapped_reads/m_4.bam
```

### Post-featureCounts Check

Examine the summary and first few lines of the counts table:

```bash
cat counts/counts_USA300_clean.txt.summary
head -n 5 counts/counts_USA300_clean.txt
```
### Counts Files Overview

- `counts/counts_USA300_clean.txt.summary`  
  Contains the **assigned/unassigned counts** per sample.

- `counts/counts_USA300_clean.txt`  
  Wide table with **Geneid**, coordinates, and **per-sample count columns**.

## 3.10 — DESeq2 Analysis R Script (No EnhancedVolcano)

**File:** `scripts/deseq2_USA300.R`
# scripts/deseq2_USA300.R

```r
# Run: Rscript scripts/deseq2_USA300.R
library(DESeq2)
library(ggplot2)
library(pheatmap)

# --- 1. Load data
counts <- read.table("counts/counts_USA300_clean.txt",
                     header = TRUE, row.names = 1, sep = "\t", comment.char = "#")

metadata <- read.table("metadata.tsv",
                       header = TRUE, row.names = 1, sep = "\t", stringsAsFactors = TRUE)

# Make sure 'condition' is factor in the order we want (chronic then acute)
metadata$condition <- factor(metadata$condition, levels = c("chronic", "acute"))

# --- 2. Prepare count matrix (drop annotation columns)
# featureCounts output: first columns are Geneid, Chr, Start, End, Strand, Length (or similar)
# keep only numeric count columns (starting from 6)
counts_mat <- counts[, -(1:5)]

# Clean sample column names to match metadata rownames
colnames(counts_mat) <- gsub("^X", "", colnames(counts_mat))               # remove leading X if R added it
colnames(counts_mat) <- gsub("mapped_reads\\.", "", colnames(counts_mat))  # remove prefix
colnames(counts_mat) <- gsub("\\.bam$", "", colnames(counts_mat))          # remove .bam extension
colnames(counts_mat) <- gsub("\\..*$", "", colnames(counts_mat))          # remove potential extra suffix

# Reorder columns to match metadata
if (!all(colnames(counts_mat) %in% rownames(metadata))) {
  stop("Some sample columns in counts are not in metadata. Please check names.")
}
counts_mat <- counts_mat[, rownames(metadata)]

# Final sanity check
stopifnot(all(colnames(counts_mat) == rownames(metadata)))

# --- 3. Create DESeq2 object
dds <- DESeqDataSetFromMatrix(countData = counts_mat,
                              colData   = metadata,
                              design    = ~ condition)

# Filter low-count genes
dds <- dds[rowSums(counts(dds)) > 10, ]

# Run DESeq2
dds <- DESeq(dds)

# --- 4. Results
res <- results(dds, contrast = c("condition", "acute", "chronic"))
res <- res[order(res$padj), ]

# Save full results table
write.csv(as.data.frame(res), file = "results/DESeq2_results_USA300.csv", row.names = TRUE)

# rlog for visualization
rld <- rlog(dds, blind = FALSE)

# --- 5. PCA
png("results/PCA_USA300.png", width = 6, height = 5, units = "in", res = 300)
plotPCA(rld, intgroup = "condition") + ggtitle("PCA: Acute vs Chronic")
dev.off()

# --- 6. Sample distance heatmap
sampleDists <- dist(t(assay(rld)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- colnames(sampleDistMatrix) <- rownames(metadata)
png("results/SampleDistanceHeatmap_USA300.png", width = 6, height = 6, units = "in", res = 300)
pheatmap(sampleDistMatrix, annotation_col = metadata)
dev.off()

# --- 7. Volcano plot (ggplot2)
volcano_data <- as.data.frame(res)
volcano_data$gene <- rownames(volcano_data)

volcano_data$threshold <- "Not Sig"
volcano_data$threshold[volcano_data$padj < 0.05 & volcano_data$log2FoldChange > 1] <- "Up"
volcano_data$threshold[volcano_data$padj < 0.05 & volcano_data$log2FoldChange < -1] <- "Down"

png("results/Volcano_USA300.png", width = 6, height = 5, units = "in", res = 300)
ggplot(volcano_data, aes(x = log2FoldChange, y = -log10(padj), color = threshold)) +
  geom_point(alpha = 0.6, size = 1.5, na.rm = TRUE) +
  scale_color_manual(values = c("Up" = "red", "Down" = "blue", "Not Sig" = "grey")) +
  theme_minimal() +
  xlab("log2 Fold Change (Acute vs Chronic)") +
  ylab("-log10 Adjusted p-value") +
  ggtitle("Volcano Plot: Acute vs Chronic PJI")
dev.off()

# --- 8. Heatmap of top 30 DEGs
topGenes <- head(order(res$padj), 30)
png("results/Heatmap_USA300_top30.png", width = 7, height = 8, units = "in", res = 300)
pheatmap(assay(rld)[topGenes, , drop = FALSE],
         annotation_col = metadata,
         show_rownames = TRUE, cluster_cols = TRUE)
dev.off()

# --- 9. Export normalized counts
norm_counts <- counts(dds, normalized = TRUE)
write.csv(as.data.frame(norm_counts), "results/NormalizedCounts_USA300.csv", row.names = TRUE)
```

## DESeq2 Analysis Outputs

All outputs from `scripts/deseq2_USA300.R` are saved in `results/`:

- **PCA plot:** `results/PCA_USA300.png`  
- **Sample distance heatmap:** `results/SampleDistanceHeatmap_USA300.png`  
- **Volcano plot:** `results/Volcano_USA300.png`  
- **Top 30 DEGs heatmap:** `results/Heatmap_USA300_top30.png`  
- **Full DESeq2 results table:** `results/DESeq2_results_USA300.csv`  
- **Normalized counts table:** `results/NormalizedCounts_USA300.csv`

# 4 — Troubleshooting — Sample Name Mismatches

 encountered the error:

```r
stopifnot(all(colnames(counts_mat) == rownames(metadata)))
```
run these checks in R:

```r
# Show names to inspect mismatch
colnames(counts)[1:10]       # raw header fields
colnames(counts_mat)         # sample columns after dropping first 5 columns
rownames(metadata)[1:20]     # metadata rownames

# Differences
setdiff(colnames(counts_mat), rownames(metadata))
setdiff(rownames(metadata), colnames(counts_mat))
```

when columns have prefixes/suffixes (e.g., mapped_reads.fm_1.bam), cleaned them using:

```r
colnames(counts_mat) <- gsub("mapped_reads\\.", "", colnames(counts_mat))
colnames(counts_mat) <- gsub("\\.bam$", "", colnames(counts_mat))
colnames(counts_mat) <- gsub("^X", "", colnames(counts_mat))
colnames(counts_mat) <- gsub("\\..*$", "", colnames(counts_mat))
```

Then reorder columns to match metadata:

```r
counts_mat <- counts_mat[, rownames(metadata)]
stopifnot(all(colnames(counts_mat) == rownames(metadata)))
```

Used these exact steps interactively; once cleaned, stopifnot passed successfully.

# 5 — GitHub: What to Push & `.gitignore`

**Do not push large binary files:**  

- `raw_sra/`  
- `mapped_reads/*.bam`  
- `genome_index/`  
- `genome/` FASTA if >50MB  
- `*.bam` or index files >100MB  

> when accidentally committed large files, remove them with `git rm --cached` and force-push.

---

## Recommended `.gitignore` Content

 edit `.gitignore` at the repo root:

```gitignore
# Raw and big data
raw_sra/
mapped_reads/
IGV/
clean_fastq/

# Trimmed fastq (optional: you may want to keep)
trim/*_R1_trim.fastq.gz
trim/*_R2_trim.fastq.gz

# Genomes & indexes (do NOT track)
genome/
genome_index/
RNA_Seq/genome/genomeIndex/

# Binary & intermediate
*.bam
*.bai
*.sam
*.sra

# R / RStudio
.Rhistory
.RData
.Rproj.user/

# Large counts / annotation if too big
# (But counts_clean.txt is allowed if under size limits)
```

### Add / Commit / Push Safe Files to GitHub

```bash
git add scripts/ docs/ counts/counts_USA300_clean.txt counts/counts_USA300_clean.txt.summary \
    metadata.tsv results/ qc/* scripts/deseq2_USA300.R

git commit -m "Add S. aureus PJI RNA-seq pipeline scripts, cleaned counts, results, metadata"

git pull --rebase origin main   # sync with remote

git push origin main
```

### Removing Large Files Already Committed

If accidentally committed large files and received remote rejected errors, follow these steps carefully:

```bash
# Remove files from git index but keep locally
git rm --cached -r genome/genomeIndex
git rm --cached genome/c_elegans.fa
git rm --cached counts/c_elegans.gff3

git commit -m "Remove large genome/index files from tracking and add .gitignore rules"

# Optional cleanup (rewriting history, be cautious)
rm -rf .git/refs/original/
git reflog expire --expire=now --all
git gc --prune=now --aggressive

# Force push if necessary
git push origin main --force
