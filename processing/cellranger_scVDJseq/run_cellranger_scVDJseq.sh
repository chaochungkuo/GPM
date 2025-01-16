#!/bin/bash

# this script will run the scVDJ analysis
# make sure cellranger is available in your environment
# otherwise scource shared enviroment from /data/shared_env/shared_paths.sh
# this script can be run directly from shell

# Each sample needs to have one multi_config.csv and one execution command.
# For multiple samples, please duplicate the command below:
cellranger multi --id=<ID-matching-mkfastq> \
--csv multi_config.csv --localcores 20

# Soft link multiQC
# Because FASTQ path for scVDJseq targets to the directory of FASTQ files, we have to go upward to the parent directory
PATH_FASTQ=$(cut -d "/" -f 1,2,3,4 <<< PROJECT_FASTQ_PATH)
ln -s ${PATH_FASTQ}/multiqc/ ../multiqc