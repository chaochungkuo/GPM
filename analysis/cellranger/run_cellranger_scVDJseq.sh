#!/bin/bash

# this script will run the scVDJ analysis
# make sure cellranger is available in your environment
# otherwise scource shared enviroment from /data/shared_env/shared_paths.sh
# this script can be run directly from shell

mkdir -p ./multi_output

Local_Cores=30

cd ./multi_output

Sample_ID=`ls ../${multi_config}*.csv`
for sample in ${Sample_ID}; do
  sample_name=${sample#*config_}
  sample_name=${sample_name%.csv*}
  Config_File="../${sample}.csv"
  cellranger multi --id=$sample_name --csv=$Config_File --localcores=$Local_Cores
done
cd ../

# Soft link multiQC
# Because FASTQ path for scVDJseq targets to the directory of FASTQ files, we have to go upward to the parent directory
PATH_FASTQ=$(cut -d "/" -f 1,2,3,4 <<< FASTQ_DIR)
ln -s ${PATH_FASTQ}/multiqc/ ../multiqc