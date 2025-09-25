#!/bin/bash
mkdir -p mapped_reads

for sample in $(ls trim | sed -n 's/_R1_trim.fastq.gz$//p' | sort -u); do
  if [ -f trim/${sample}_R1_trim.fastq.gz ] && [ -f trim/${sample}_R2_trim.fastq.gz ]; then
    bwa mem -t 8 genome/s_aureus_USA300.fa trim/${sample}_R1_trim.fastq.gz trim/${sample}_R2_trim.fastq.gz \
      | samtools view -bS - \
      | samtools sort -o mapped_reads/${sample}.bam
  else
    bwa mem -t 8 genome/s_aureus_USA300.fa trim/${sample}_trim.fastq.gz \
      | samtools view -bS - \
      | samtools sort -o mapped_reads/${sample}.bam
  fi
  samtools index mapped_reads/${sample}.bam
done
