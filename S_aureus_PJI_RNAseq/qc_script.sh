fastqc -o qc/ clean_fastq/*.gz
multiqc qc/ -o qc/multiqc_report
