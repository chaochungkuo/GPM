################## test script ##################
# PATH_NEXTFLOW run nf-core/ampliseq -profile test,docker --outdir ./results
# switch to 16s enviroment before run this pipe

PATH_NEXTFLOW run nf-core/ampliseq  -profile docker \
    --input samplesheet.csv \
    --metadata metadata.tsv \
    --FW_primer "YYYRYGRDDBVCWSCA" \
    --RV_primer "GAHTACNVRRGTNTCTAAKYY" \
    --min_frequency 10 \
    --outdir ./results_16s 

#    --FW_primer "CCTACGGGDGGCWGCAG", "CCTAYGGGGYGCWGCAG" \
#    --RV_primer "GACTACNVGGGTMTCTAATCC" \
#   --FW_primer "GTGYCAGCMGCCGCGGTAA" \ 
#   --RV_primer "GGACTACNVGGGTWTCTAAT" \
