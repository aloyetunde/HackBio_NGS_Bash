#!/bin/bash

# 1. Print your name
echo "Alo Yetunde"

# 2. Create a folder titled with your name
mkdir -p "Alo Yetunde"

# 3. Create another new directory titled 'biocomputing' and change to that directory
mkdir -p biocomputing && cd biocomputing

# 4. Download these 3 files
wget https://raw.githubusercontent.com/josoga2/dataset-repos/main/wildtype.fna
wget https://raw.githubusercontent.com/josoga2/dataset-repos/main/wildtype.gbk
wget https://raw.githubusercontent.com/josoga2/dataset-repos/main/wildtype.gbk -O wildtype_dup.gbk

# 5. Move the `.fna` file to the folder titled your name
mv wildtype.fna ../"Alo Yetunde"

# 6. Delete the duplicate `.gbk` file
  rm wildtype_dup.gbk

# 7. Confirm if the `.fna` file is mutant or wild type (tatatata vs tata)
cd ../"Alo Yetunde"
if grep -q "tatatata" wildtype.fna; then     echo "Mutant"; else     echo "Wild type"; fi
   
# 8. If mutant, print all matching lines into a new file
   grep 'tatatata' wildtype.fna > mutant_matches.txt

# 9. Count number of lines (excluding header) in the .gbk file
cd ../biocomputing
# Header line starts with LOCUS or so, count non-header lines (exclude lines starting with LOCUS)
grep -v '^LOCUS' wildtype.gbk | wc -l

# 10. Print the sequence length of the .gbk file. 
head -1 wildtype.gbk | awk '{print $3}'

# 11. Print the source organism of the .gbk file. (Use the SOURCE tag)
grep '^SOURCE' wildtype.gbk | head -1 | awk '{$1=""; print $0}'

# 12. List all the gene names of the .gbk file (grep '/gene=')
grep "/gene=" wildtype.gbk

# 13. Clear your terminal space and print all commands used today
clear
history

# 14. List the files in the two folders
cd ../
echo "Files in '${PWD}/Alo Yetunde':"
ls -l "Alo Yetunde"
echo "Files in '${PWD}/biocomputing':"
ls -l biocomputing

#!/bin/bash
# Stage 0 - Project 1: BASh Basic

# 1. Print your name
echo "Alo Adejumoke"

# 2. Create a folder titled your name
mkdir -p "Alo_Adejumoke"

# 3. Create another new directory titled biocomputing and change to that directory with one line of command
mkdir -p biocomputing && cd biocomputing

# 4. Download the 3 files
wget https://raw.githubusercontent.com/josoga2/dataset-repos/main/wildtype.fna
wget https://raw.githubusercontent.com/josoga2/dataset-repos/main/wildtype.gbk
wget https://raw.githubusercontent.com/josoga2/dataset-repos/main/wildtype.gbk -O wildtype_dup.gbk

# 5. Move the .fna file to the folder titled your name
mv wildtype.fna ../Alo_Adejumoke/

# 6. Delete the duplicate gbk file
rm wildtype_dup.gbk

# 7. Confirm if the .fna file is mutant or wild type (tatatata vs tata)
echo "Checking if mutant (tatatata) or wild type (tata)..."
grep "tatatata" ../Alo_Adejumoke/wildtype.fna && echo "Mutant found" || echo "Wild type"

# 8. If mutant, print all matching lines into a new file
grep "tatatata" ../Alo_Adejumoke/wildtype.fna > mutant_lines.txt

# 9. Count number of lines (excluding header) in the .gbk file
echo "Number of lines (excluding header) in .gbk:"
grep -v "^LOCUS" wildtype.gbk | wc -l

# 10. Print the sequence length of the .gbk file. (Use LOCUS tag)
echo "Sequence length:"
grep "LOCUS" wildtype.gbk | awk '{print $3}'

# 11. Print the source organism of the .gbk file. (Use SOURCE tag)
echo "Source organism:"
grep "SOURCE" wildtype.gbk | head -n 1

# 12. List all the gene names of the .gbk file.
echo "List of gene names:"
grep "/gene=" wildtype.gbk

# 13. Clear your terminal space and print all commands used today
history > commands_used.txt
clear
cat commands_used.txt

# 14. List the files in the two folders
echo "Files in Alo_Adejumoke folder:"
ls ../Alo_Adejumoke
echo "Files in biocomputing folder:"
ls .
