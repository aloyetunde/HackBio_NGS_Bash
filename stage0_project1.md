# BASh Basic Project 1 Solution - Alo Yetunde

This document details the step-by-step bash commands used to complete the Project 1: BASh Basic, as per the project requirements. The script was created and tested in Google Cloud Shell.

---

## Tasks and Bash Commands

### 1. Print your name
echo "Alo Yetunde"

### 2. Create a folder titled your name
mkdir -p "Alo Yetunde"


### 3. Create another new directory titled `biocomputing` and change to that directory
mkdir -p biocomputing && cd biocomputing

### 4. Download these 3 files
wget https://raw.githubusercontent.com/josoga2/dataset-repos/main/wildtype.fna
wget https://raw.githubusercontent.com/josoga2/dataset-repos/main/wildtype.gbk
wget https://raw.githubusercontent.com/josoga2/dataset-repos/main/wildtype.gbk

text

### 5. Move the `.fna` file to the folder titled your name
mv wildtype.fna ../"Alo Yetunde"

text

### 6. Delete the duplicate `.gbk` file
Since the duplicate `.gbk` file is overwritten by the second download, no additional deletion command is necessary.

### 7. Confirm if the `.fna` file is mutant or wild type (tatatata vs tata)
cd ../"Alo Yetunde"
if grep -q 'tatatata' wildtype.fna; then
echo "Mutant"
else
echo "Wild type"
fi

text

### 8. If mutant, print all matching lines into a new file
grep 'tatatata' wildtype.fna > mutant_matches.txt

text

### 9. Count number of lines (excluding header) in the `.gbk` file
cd ../biocomputing
grep -v '^LOCUS' wildtype.gbk | wc -l

text

### 10. Print the sequence length of the `.gbk` file (Use the LOCUS tag in the first line)
head -1 wildtype.gbk | awk '{print $3}'

text

### 11. Print the source organism of the `.gbk` file (Use the SOURCE tag)
grep '^SOURCE' wildtype.gbk | head -1 | awk '{$1=""; print $0}'

text

### 12. List all the gene names of the `.gbk` file
grep "/gene=" wildtype.gbk

text

### 13. Clear your terminal space and print all commands used today
clear
history

text

### 14. List the files in the two folders
cd ../
echo "Files in 'Alo Yetunde':"
ls -l "Alo Yetunde"
echo "Files in 'biocomputing':"
ls -l biocomputing


