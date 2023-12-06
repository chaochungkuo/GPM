################## test script ##################
# PATH_NEXTFLOW run nf-core/rnaseq -profile test,docker


#default
PATH_NEXTFLOW run nf-core/bacass -profile docker  --input samplesheet.csv --kraken2db "https://genome-idx.s3.amazonaws.com/kraken/k2_standard_8gb_20210517.tar.gz" -bg

#wo kraken
PATH_NEXTFLOW run nf-core/bacass -profile docker  --input samplesheet.csv  -bg --skip_kraken2

#
PATH_NEXTFLOW run nf-core/bacass -profile docker  --input samplesheet.csv --kraken2db "https://genome-idx.s3.amazonaws.com/kraken/k2_standard_8gb_20210517.tar.gz" --multiqc_title "Genomics_Facility_IZKF_Aachen_Sequencing_Quality_Report" -bg

# Options for --genome:
# gencode_hg38, gencode_mm10, hg38, mm10
# Useful options: --removeRiboRNA
