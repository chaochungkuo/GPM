################## test script ##################
# PATH_NEXTFLOW run nf-core/bacass -profile test,docker
# https://nf-co.re/bacass/2.0.0/usage

PATH_NEXTFLOW run nf-core/bacass -r 2.0.0 -profile docker \
     --input samplesheet.csv

# samplesheet.csv
# ID    R1    R2    LongFastQ    Fast5    GenomeSize
# shortreads    ./data/S1_R1.fastq.gz    ./data/S1_R2.fastq.gz    NA    NA    NA
# longreads    NA    NA    ./data/S1_long_fastq.gz    ./data/FAST5    2.8m
# shortNlong    ./data/S1_R1.fastq.gz    ./data/S1_R2.fastq.gz    ./data/S1_long_fastq.gz    ./data/FAST5    2.8m