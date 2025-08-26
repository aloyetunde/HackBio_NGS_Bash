#!/bin/bash

# 1. Print your name
echo "Alo Yetunde"

# 2. Create a folder titled with your name
mkdir -p "Alo Yetunde"

# 3. Create another new directory titled 'biocomputing' and change to that directory
mkdir -p biocomputing && cd biocomputing

# 4. Download the 3 files
wget https://raw.githubusercontent.com/josoga2/dataset-repos/main/wildtype.fna
wget https://raw.githubusercontent.com/josoga2/dataset-repos/main/wildtype.gbk
wget https://raw.githubusercontent.com/josoga2/dataset-repos/main/wildtype.gbk  # Duplicate intentionally as per project, will handle next 

# 5. Move the .fna file to the folder titled your name
mv wildtype.fna ../"Alo Yetunde"

# 6. Delete the duplicate .gbk file
# Assuming the duplicate created by the second wget above, remove the second file if exists (files have the same name)
# Since it's downloaded twice with the same name, the second one overwrites the first, so no duplicate exists physically,
# but to ensure, just keep one .gbk file
# We'll keep one .gbk as is, no action needed in bash unless duplicates exist

# 7. Confirm if the .fna file is mutant or wild type (tatatata vs tata)
# We'll grep for 'tatatata' in the .fna file
cd ../"Alo Yetunde"
if grep -q 'tatatata' wildtype.fna; then
  echo "Mutant"
  # 8. If mutant, print all matching lines into a new file
  grep 'tatatata' wildtype.fna > mutant_matches.txt
else
  echo "Wild type"
fi

# 9. Count number of lines (excluding header) in the .gbk file
cd ../biocomputing
# Header line starts with LOCUS or so, count non-header lines (exclude lines starting with LOCUS)
grep -v '^LOCUS' wildtype.gbk | wc -l

# 10. Print the sequence length of the .gbk file. (Use the LOCUS tag in the first line)
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

