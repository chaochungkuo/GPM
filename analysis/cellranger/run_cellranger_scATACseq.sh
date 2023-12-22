#!/bin/bash

# this script will run the scATACseq analysis
# make sure cellranger is available in your environment
# otherwise scource shared enviroment from /data/shared_env/shared_paths.sh
# this script can be run directly from shell

mkdir -p ./scATACseq_output

Local_Cores=30
refGenome="/data/shared_env/10xGenomics/refdata-cellranger-arc-mm10-2020-A-2.0.0"

cd ./scATACseq_output

sampleid=""
fastqDir=""
cellranger-atac count --sample=$sampleid \
					--fastqs=$fastqDir \
					--id=$sampleid \
					--reference=$refGenome \
					--localcores=${Local_Cores} \
					--localmem=120 \

# Soft link multiQC
# Because FASTQ path for scVDJseq targets to the directory of FASTQ files, we have to go upward to the parent directory
PATH_FASTQ=$(cut -d "/" -f 1,2,3,4 <<< PROJECT_FASTQ_DIR)
ln -s ${PATH_FASTQ}/multiqc/ ../multiqc