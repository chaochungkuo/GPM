# screen -S fastq
# bash run_splitpipe.sh

############################################################
# Demultiplex the reads for sublibraries
PATH_BCL2FASTQ \
  --no-lane-splitting \
  --runfolder-dir PROJECT_BCL_PATH \
  --output-dir ./sublibraries_FASTQ \
  --interop-dir ./InterOp/ \
  --sample-sheet ./samplesheet_bcl2fastq.csv \
  --processing-threads N_CORES

###### Running FASTQC ######################################
mkdir -p ./fastqc
find * -maxdepth 1 -name "*.fastq.gz" | parallel -j N_CORES "fastqc {} -o ./fastqc"

###### Running fastq_screen ###############################
mkdir -p ./fastq_screen
fastq_screen --outdir ./fastq_screen --threads N_CORES */*.fastq.gz *.fastq.gz

###### Running MultiQC #####################################
mkdir -p multiqc
multiqc -f ./fastq_screen/ . -o ./multiqc

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
    --sample NAME1 A1-A6 \
    --sample NAME2 C1-C6 \
    --sample NAME3 D1-D6

# Sublibrary 2 and so on...
# PATH_SPLITPIPE \
#     --nthreads N_CORES \
#     --mode all \
#     --chemistry v2 \
#     --genome_dir /newvolume/genomes/hg38_GRCm39/ \
#     --fq1 ./sublibraries_FASTQ/s2_S1_R1_001.fastq.gz \
#     --output_dir ./results_s1/ \
#     --sample NAME1 A1-A6 \
#     --sample NAME2 C1-C6 \
#     --sample NAME3 D1-D6

# --sample <sample name> <wells> specifies the sample name and corresponding wells in which each sample was loaded into the first round of barcoding. Each of these lines consists of the sample command line flag followed by a name and a specification of which wells to include. To see formatting rules for specifying well subsets, call split-pipe --explain. Only valid wells for the chosen kit are allowed. By default, a combined analysis of all samples (named 'all-sample') is also performed by the pipeline.

# After running each sublibrary through the pipeline, the outputs of each sublibrary can be combined into a single dataset using split-pipe --mode combine.
PATH_SPLITPIPE\
    --mode comb \
    --sublibraries ./results_s1/ \
                   ./results_s2/ \
    --output_dir ./results_combined/

