mkdir -p trim
for base in $(ls clean_fastq/*_R1.fastq.gz | sed -E 's/_R1.fastq.gz$//' ); do
  sample=$(basename "$base")
  fastp -i clean_fastq/${sample}_R1.fastq.gz -I clean_fastq/${sample}_R2.fastq.gz \
        -o trim/${sample}_R1_trim.fastq.gz -O trim/${sample}_R2_trim.fastq.gz \
        -h trim/${sample}_fastp.html -j trim/${sample}_fastp.json --detect_adapter_for_pe
done
