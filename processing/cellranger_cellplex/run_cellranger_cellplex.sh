#!/bin/bash

# this script will run the cellranger multi  analysis for cells multiplexed with CellPlex
# make sure cellranger is available in your environment
# otherwise scource shared enviroment from /data/shared_env/shared_paths.sh
# config file has to be prepared before running the analysis


# first generate template for config file and modify it 

cellranger multi-template --output=/path/to/FILE.csv

# run analysis 

Local_Cores=N_CORES

mkdir -p ./results
cd ./results

sampleid="sample"


CELLRANGER_BASE multi --id=$sampleid \
                    --csv="PATH_TO_CONFIG_FILE" \
                    --localcores=${Local_Cores} \
                    --localmem=120 