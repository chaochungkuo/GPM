################## test script ##################
# PATH_NEXTFLOW run nf-core/atacseq -profile test,docker

PATH_NEXTFLOW run nf-core/atacseq -profile docker \
     --input design.csv \
     --single_end \
     --genome gencode_hg38 \ # Please define the genome ID: hg38, mm10
     -name TITLE_NAME

# Other optional parameters:
# --narrow_peak
# --save_macs_pileup

# Narrow Peak Calling
# If your target protein is a transcription factor, you should probably choose narrow peak calling. You can also try the narrow peak calling workflows for the following histone marks:
# H3K4me3
# H3K4me2
# H3K9-14ac
# H3K27ac
# H2A.Z

# Broad Peak Calling
# You should try the broad peak calling workflows for the following histone marks:
# H3K36me3
# H3K79me2
# H3K27me3
# H3K9me3
# H3K9me1

# Special Cases
# In some scenarios, H3K4me1, H3K9me2 and H3K9me3 might behave between narrow and broad shape, you might need to look into each peak region and consult experts.