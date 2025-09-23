# create directories for output
# script for mapping reads


mkdir -p mapped_reads
#! /bin/bash
for infile in trim/*.fq ; do
	outfile=$(basename "$infile".fq)
	STAR --genomeDir genome/genomeIndex --readFilesIn $infile --outFileNamePrefix mapped_reads/$outfile --outSAMtype BAM SortedByCoordinate --outSAMattributes All
done

