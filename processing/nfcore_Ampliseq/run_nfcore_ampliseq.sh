################## test script ##################
# nextflow run nf-core/rnaseq -profile test,docker

nextflow run nf-core/ampliseq -profile docker \
    --input samplesheet.csv \
    --FW_primer GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAGcgtgacgtagaaagtaataa \
    --RV_primer TCGTCGGCAGCGTCAGATGTGTATAAGAGACAGgactagccttattttaacttgct \
    --metadata metadata.tsv
    --outdir ./results
