# 🚨 South African Polony Outbreak (2017–2018) – WGS Analysis Report

## 📌 Introduction

In 2017–2018, South Africa experienced the **world’s largest listeriosis outbreak**, with **978 confirmed cases** and **183 deaths** (case fatality rate of 27%).  
The outbreak was traced to **polony** produced at Enterprise Foods (Polokwane facility).  

This project aimed to analyze **Whole Genome Sequencing (WGS)** data from >100 bacterial isolates to:

- Confirm the identity of the pathogen  
- Determine its antimicrobial resistance (AMR) profile  
- Detect toxin genes linked to virulence  
- Suggest treatment options to guide outbreak management  

---

## ⚙️ Methods & Tools

We analyzed the dataset using the following bioinformatics tools:

| Step | Tool | Purpose | Output |
|------|------|---------|--------|
| 1 | **FastQC** | Raw data quality check | Quality reports (per-base score, adapters, GC%) |
| 2 | **FastP** | Trimming & filtering | Clean FASTQ reads |
| 3 | **SPAdes** | Genome assembly | Contigs for each isolate |
| 4 | **QUAST** | Assembly quality | N50, genome length, #contigs |
| 5 | **BLAST** | Organism confirmation | Listeria monocytogenes |
| 6 | **Abricate (CARD)** | AMR gene detection | Resistance genes list |
| 7 | **Custom R script** | AMR prevalence summary | Heatmap & frequency table |
| 8 | **Abricate (VFDB)** | Toxin gene search | hly, plcA, plcB detected |

---

## 🖥️ Code & Analysis

### 1️⃣ Download dataset
```bash
wget https://raw.githubusercontent.com/HackBio-Internship/2025_project_collection/refs/heads/main/SA_Polony_100_download.sh
bash SA_Polony_100_download.sh
```

### 2️⃣ Quality control
```bash
fastqc *.fastq.gz -o qc_reports/
```

### 3️⃣ Trimming reads
```bash
fastp -i sample_R1.fastq.gz -I sample_R2.fastq.gz \
      -o sample_R1_trimmed.fastq.gz -O sample_R2_trimmed.fastq.gz
```

### 4️⃣ Genome assembly
```bash
spades.py -1 sample_R1_trimmed.fastq.gz -2 sample_R2_trimmed.fastq.gz \
          -o assembly_output

```

### 5️⃣ Assembly assessment
```bash
quast.py assembly_output/contigs.fasta -o quast_results/

```
### 6️⃣ Confirm organism Identity (BLAST)
```bash
blastn -query assembly_output/contigs.fasta -db nt -out results.txt -outfmt 6
```
### Output: All isolates matched Listeria monocytogenes

### 7️⃣ Detect AMR genes
```bash
abricate --db card assembly_output/contigs.fasta > amr_report.txt
```
### 8️⃣ Summarize AMR prevalence (R)
```r
library(dplyr)
library(ggplot2)

amr <- read.delim("amr_report.txt")
summary <- amr %>% count(GENE, sort=TRUE)

ggplot(summary, aes(x=reorder(GENE, n), y=n)) +
  geom_bar(stat="identity", fill="steelblue") +
  coord_flip() +
  labs(title="AMR Gene Prevalence", x="Gene", y="Frequency")
```
### 9️⃣ Search for toxin genes
```bash
abricate --db vfdb assembly_output/contigs.fasta > toxins_report.txt
```
### Output: Genes hly, plcA, plcB detected in majority of isolates.

## 📊 Results

**Organism identity**: ***Listeria monocytogenes***

**AMR genes detected**: High prevalence of **tetracycline and aminoglycoside** resistance genes

**Toxin genes**: hly, plcA, plcB widely present

**Recommended treatment**: **Ampicillin + Gentamicin** (aligns with guidelines, but resistance monitoring essential)

## 🧩 Discussion

- Genomics enabled **rapid confirmation of the outbreak source** (polony).

- **AMR gene profiling** revealed concerning resistance patterns, guiding treatment decisions.

- **Toxin gene presence** explains the unusually high case fatality rate.

- Public health response included product recalls and regional alerts, preventing further spread.

  ## ✅ Conclusion

This project demonstrates how NGS + bioinformatics save lives during outbreaks by:
✔️ Identifying pathogens quickly
✔️ Detecting resistance before it spreads
✔️ Guiding evidence-based treatment
✔️ Informing public health strategies

## 🔗 Repository & Outreach

📂 GitHub Repo: [https://github.com/aloyetunde/HackBio_NGS_Bash]

📢 LinkedIn Post:[(https://www.linkedin.com/feed/update/urn:li:activity:7374052083060281344/)]

## 👩‍💻 Personal Takeaways

- Gained hands-on experience with a full NGS pipeline (QC → Assembly → AMR → Toxins).

- Improved skills in Linux, BLAST, abricate, and R visualization.

- Learned how genomics data links directly to real-world health decisions.

- Inspired by how bioinformatics transforms outbreak response.
