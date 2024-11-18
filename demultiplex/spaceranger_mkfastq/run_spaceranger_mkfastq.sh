# Please make sure spaceranger-2.1.0 is available in your environment
# You should review the paramenters and run:
# screen -S spaceranger
# bash run_spaceranger.sh

# Please execute this command in the directory OUTPUT_DIR
spaceranger mkfastq --id=mkfastq --localcores=40 \
                   --run=PROJECT_BCL_PATH \
                   --csv=./samplesheet_spaceranger.csv

bash run_merge_lanes.sh mkfastq multiqc