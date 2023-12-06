################## test script ##################
# PATH_NEXTFLOW run nf-core/rnaseq -profile test,docker

PATH_NEXTFLOW run nf-core/ampliseq -profile docker \
    --input samplesheet.csv \
    --FW_primer GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAGcgtgacgtagaaagtaataa \
    --RV_primer TCGTCGGCAGCGTCAGATGTGTATAAGAGACAGgactagccttattttaacttgct \
    --metadata metadata.tsv
    --outdir ./results
