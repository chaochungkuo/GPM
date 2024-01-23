################## test script ##################
# PATH_NEXTFLOW run nf-core/chipseq -profile test,docker

PATH_NEXTFLOW run nf-core/chipseq -r 2.0.0 -profile docker \
     --input samplesheet.csv \
     --single_end --fragment_size 300 \
     --blacklist /data/genomes/chipseq_blacklists/v3.0/hg38-blacklist.v3.bed \
     --genome GENCODE_GRCh38_v44 --macs_gsize 2.9e9

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

# Genome Effective size
# GRCh37
# 2864785220
# GRCh38
# 2913022398
# T2T/CHM13CAT_v2
# 3117292070
# GRCm37
# 2620345972
# GRCm38
# 2652783500
# dm3
# 162367812
# dm6
# 142573017
# GRCz10
# 1369631918
# GRCz11
# 1368780147
# WBcel235
# 100286401
# TAIR10
# 119482012

