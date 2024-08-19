#!/bin/bash

# this script will run Fixed RNA Profiling with Cell Ranger multi
# https://www.10xgenomics.com/support/software/cell-ranger/latest/analysis/running-pipelines/cr-flex-multi-frp
# make sure cellranger is available in your environment
# otherwise scource shared enviroment from /data/shared_env/shared_paths.sh
# this script can be run directly from shell

Local_Cores=30

cellranger multi --id="FRP" --csv="multi_config.csv" --localcores=$Local_Cores

# If multiple libraries are included in the multi_config.csv file, you don't have to run cellranger aggr.