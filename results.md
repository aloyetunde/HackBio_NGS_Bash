# Results Summary
# Results Summary

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









