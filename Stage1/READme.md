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

## 1. Confirmed Identity of the Pathogen
Using BLAST alignment, all isolates were confirmed to be **_Listeria monocytogenes_**, the causative agent of the South African polony outbreak (2017–2018).

---

## 2. Antimicrobial Resistance (AMR) Profiles

The CARD database (`abricate --db card`) was used to identify AMR genes.  
Below is a summary of the most prevalent resistance genes across isolates:

### Antimicrobial Resistance (AMR) Gene Prevalence

The following table summarizes resistance genes detected across the outbreak isolates using **abricate (CARD database)**.

| AMR Gene | Function / Resistance Mechanism | Count (n=10) | Prevalence (%) |
|----------|---------------------------------|--------------|----------------|
| **FosX** | Fosfomycin resistance            | 10           | 100%           |
| **mprF** | Membrane modification; antimicrobial peptide resistance (Listeria survival factor) | 10 | 100% |
| **lin**  | Lincosamide resistance (e.g., clindamycin) | 10 | 100% |
| **norB** | Efflux pump; fluoroquinolone resistance | 10 | 100% |

---

**Interpretation:**  
All isolates carried resistance determinants to **fosfomycin, lincosamides, and fluoroquinolones**, along with `mprF`, which contributes to resistance against host antimicrobial peptides.  
This indicates that these drugs are not viable treatment options for the outbreak strains.  


## 3. Toxin Gene Profiles (VFDB)

The VFDB database (`abricate --db vfdb`) was used to detect virulence and toxin genes.
### Toxin Gene Prevalence

The following table summarizes virulence/toxin genes detected across the outbreak isolates using **abricate (VFDB database)**.

| Toxin Gene | Function / Role in Pathogenesis | Count (n=10) | Prevalence (%) |
|------------|---------------------------------|--------------|----------------|
| **actA**   | Actin polymerization, intracellular motility | 10 | 100% |
| **bsh**    | Bile salt hydrolase, survival in intestine  | 10 | 100% |
| **clpC**   | Stress response chaperone protease         | 10 | 100% |
| **clpE**   | Stress tolerance & virulence regulator     | 10 | 100% |
| **clpP**   | Protease, virulence factor                 | 10 | 100% |
| **fbpA**   | Fibronectin-binding protein (adhesion)     | 10 | 100% |
| **gtcA**   | Teichoic acid glycosylation (cell wall modification) | 9  | 90%  |
| **hly**    | Listeriolysin O (major pore-forming toxin) | 10 | 100% |
| **hpt**    | Hexose phosphate transporter (intracellular survival) | 10 | 100% |

---

**Interpretation:**  
- Core virulence factors such as **hly (Listeriolysin O)**, **actA**, and **bsh** were found in all isolates, explaining the **high pathogenicity** of the outbreak strain.  
- The presence of **clp genes** indicates stress tolerance and survival advantages in hostile environments (e.g., food processing, host immune system).  
- **gtcA** was slightly less prevalent (90%), suggesting some variability in cell wall glycosylation among isolates.  



---

## 4. Recommended Treatment Options

Based on the AMR profile:  
- **Ampicillin + Gentamicin** remains the gold standard for *Listeria* infections.  
- Alternatives: **Trimethoprim-sulfamethoxazole** may be effective in severe cases or in patients allergic to beta-lactams.  
- Avoid tetracyclines, macrolides, and fosfomycin due to resistance markers.

---

## 5. Public Health Implications

- WGS confirmed the outbreak strain as *Listeria monocytogenes*.  
- The AMR and toxin profile highlights why this outbreak had a high mortality rate, especially among neonates, pregnant women, and immunocompromised patients.  
- Results underscore the importance of **food safety monitoring** and **genomic surveillance** in preventing future outbreaks.  






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
