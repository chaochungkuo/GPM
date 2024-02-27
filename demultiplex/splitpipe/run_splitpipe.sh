############################################################
# Demultiplex the sublibrary FASTQ into samples
# REF_GENOME=/data/genomes/spipe/GRCh38
# REF_GENOME=/data/genomes/spipe/GRCm39
# REF_GENOME=/data/genomes/spipe/GRCm38

# Sublibrary 1
PATH_SPLITPIPE \
    --nthreads N_CORES \
    --mode all \
    --chemistry v2 \
    --genome_dir /newvolume/genomes/hg38_GRCm39/ \
    --fq1 ./sublibraries_FASTQ/s1_S1_R1_001.fastq.gz \
    --output_dir ./results_s1/ \
    --sample sample1 A1-A3 \
    --sample sample2 A4-A6 \
    --sample sample3 A7-A9 \
    --sample sample4 A10-A12

# Sublibrary 2
PATH_SPLITPIPE \
    --nthreads N_CORES \
    --mode all \
    --chemistry v2 \
    --genome_dir /newvolume/genomes/hg38_GRCm39/ \
    --fq1 ./sublibraries_FASTQ/s2_S2_R1_001.fastq.gz \
    --output_dir ./results_s2/ \
    --sample sample1 A1-A3 \
    --sample sample2 A4-A6 \
    --sample sample3 A7-A9 \
    --sample sample4 A10-A12

# --sample <sample name> <wells> specifies the sample name and corresponding wells in which each sample was loaded into the first round of barcoding. Each of these lines consists of the sample command line flag followed by a name and a specification of which wells to include. To see formatting rules for specifying well subsets, call split-pipe --explain. Only valid wells for the chosen kit are allowed. By default, a combined analysis of all samples (named 'all-sample') is also performed by the pipeline.

# After running each sublibrary through the pipeline, the outputs of each sublibrary can be combined into a single dataset using split-pipe --mode combine.
PATH_SPLITPIPE\
    --mode comb \
    --sublibraries ./results_s1/ \
                   ./results_s2/ \
    --output_dir ./results_combined/
