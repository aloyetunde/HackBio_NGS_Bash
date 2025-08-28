#!/bin/bash
# Stage 0 - Project 1: BASh Basic

# 1. Print your name
echo "Alo Yetunde"

# 2. Create a folder titled your name
mkdir -p "Alo Yetunde"

# 3. Create another new directory titled biocomputing and change to that directory
mkdir -p biocomputing && cd biocomputing

# 4. Download the 3 files
wget https://raw.githubusercontent.com/josoga2/dataset-repos/main/wildtype.fna
wget https://raw.githubusercontent.com/josoga2/dataset-repos/main/wildtype.gbk
wget https://raw.githubusercontent.com/josoga2/dataset-repos/main/wildtype.gbk -O wildtype_dup.gbk

# 5. Move the .fna file to the folder titled your name
mv wildtype.fna ../"Alo Yetunde"/

# 6. Delete the duplicate gbk file
rm wildtype_dup.gbk

# 7. Confirm if the .fna file is mutant or wild type (tatatata vs tata)
cd ../"Alo Yetunde"
if grep -q "tatatata" wildtype.fna; then
    echo "Mutant"
    # 8. If mutant, print all matching lines into a new file
    grep "tatatata" wildtype.fna > mutant_matches.txt
else
    echo "Wild type"
fi

# 9. Count number of lines (excluding header) in the .gbk file
cd ../biocomputing
echo "Number of lines (excluding header) in .gbk:"
grep -v "^LOCUS" wildtype.gbk | wc -l

# 10. Print the sequence length of the .gbk file. (Use LOCUS tag)
echo "Sequence length:"
grep "^LOCUS" wildtype.gbk | awk '{print $3}'

# 11. Print the source organism of the .gbk file. (Use SOURCE tag)
echo "Source organism:"
grep "^SOURCE" wildtype.gbk | head -1 | awk '{$1=""; print $0}'

# 12. List all the gene names of the .gbk file
echo "Gene names:"
grep "/gene=" wildtype.gbk

# 13. Clear your terminal space and print all commands used today
clear
history

# 14. List the files in the two folders
cd ..
echo "Files in Alo Yetunde folder:"
ls -l "Alo Yetunde"
echo "Files in biocomputing folder:"
ls -l biocomputing
