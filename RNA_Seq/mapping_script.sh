# create directories for output
# script for mapping reads


mkdir -p mapped
#! /bin/bash
for infile in trim/*.fq ; do
	outfile=$(basename "$infile".fq)
	STAR --genomeDir genome/genomeIndex --readFilesIn $infile --outFileNamePrefix mapped/$outfile --outSAMtype BAM SortedByCoordinate --outSAMattributes All
done

