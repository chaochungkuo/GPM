################## test script ##################
# nextflow run nf-core/ampliseq -profile test,docker --outdir ./results

nextflow run nf-core/ampliseq  -profile docker \
    --input samplesheet.csv \
    --metadata metadata.tsv \
    --FW_primer "YYYRYGRDDBVCWSCA" \
    --RV_primer "GAHTACNVRRGTNTCTAAKYY" \
#    --FW_primer "CCTACGGGDGGCWGCAG", "CCTAYGGGGYGCWGCAG" \
#    --RV_primer "GACTACNVGGGTMTCTAATCC" \
#   --FW_primer "GTGYCAGCMGCCGCGGTAA" \ 
#   --RV_primer "GGACTACNVGGGTWTCTAAT" \
    --outdir ./results
