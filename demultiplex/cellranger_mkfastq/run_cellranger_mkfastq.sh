# Please make sure cellranger-7.1.0 is available in your environment
# You should review the paramenters and run:
# screen -S cellranger
# bash run_cellranger.sh

# Please execute this command in the directory OUTPUT_DIR
PATH_CELLRANGER mkfastq --id=mkfastq --localcores=N_CORES \
                   --run=BCL_PATH \
                   --csv=./samplesheet_cellranger.csv

bash run_merge_lanes.sh