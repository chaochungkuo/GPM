#!/bin/bash

# this script will run the scATACseq analysis
# make sure cellranger is available in your environment
# otherwise scource shared enviroment from /data/shared_env/shared_paths.sh
# this script can be run directly from shell

Local_Cores=N_CORES
refGenome="REFDATA_cellranger_arc_mm10_2020_A"

mkdir -p ./results
cd ./results

sampleid="sample"
fastqDir="FASTQ_PATH"
PATH_CELLRANGER_ATAC count --sample=$sampleid \
					--fastqs=$fastqDir \
					--id=$sampleid \
					--reference=$refGenome \
					--localcores=${Local_Cores} \
					--localmem=120 \
