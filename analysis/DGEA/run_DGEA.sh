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


