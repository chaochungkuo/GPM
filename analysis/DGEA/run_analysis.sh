#!/bin/bash

################## anaylsis #####################
# Before Executing this script, please make sure to correctly: 
# (For extended setting up instractions please review the README.md)
# 1. Adjust config.yml with the appropriate parameters for your project
# 2. In the samplesheet.csv file rename the generated columns accordingly to your experimental structure.
# 3. Add your comparisons to the contrasts.csv file.
#
#
# Following that, please run"
# screen -S analysis
# bash run_analysis.sh
#
#
################## Generate the Full analysis report #####################

python rmd_report_extender.py

################## render the report #####################
# Store the previously activated conda environment
previous_env=$(conda info --envs | grep "*" | awk '{print $1}')

source activate /opt/miniconda3/envs/rstudio
Rscript -e 'rmarkdown::render("Analysis_Report_RNAseq.Rmd", output_format = "html_document", output_file = "Analysis_Report_RNAseq.html")'
# Activate the previous conda environment
source activate $previous_env
